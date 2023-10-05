# Wrapper for the Klujn et al. 2015 flux footprint model

import os
import utm_zone
import numpy as np
import pandas as pd
import configparser
import geopandas as gpd
from functools import partial
import matplotlib.pyplot as plt
from multiprocessing import Pool
from Klujn_2015_Model import FFP
from shapely.geometry import Polygon

import rasterio
from rasterio import features
from rasterio.transform import from_origin


class RunClimatology():

    def __init__(self,Site_code,Site_Config,Subsets=None,Basemap=None,Basemap_Class=None):
        self.Subsets = Subsets
        self.Site_code = Site_code
        self.ini = configparser.ConfigParser()
        self.ini.read('../MicrometPy.ini')
        self.ini.read('configuration.ini')
        self.ini.read(Site_Config)
        if Basemap is not None and Basemap_Class is not None:
            self.ini[self.Site_code]['basemap'] = Basemap
            self.ini[self.Site_code]['basemap_class'] = Basemap_Class
        else:
            self.ini[self.Site_code]['basemap'] = 'N/A'
            self.ini[self.Site_code]['basemap_class'] = 'N/A'
        
        # Dump Site to a Dataframe
        Site = pd.DataFrame(data=dict(self.ini[self.Site_code]),index=[0])
        for c in Site.columns:
            try:
                Site[c]=Site[c].astype('float64')
            except:
                pass
        
        self.lon_lat = list(Site[['longitude','latitude']].values[0])
        self.EPSG = utm_zone.epsg(self.lon_lat)
        
        self.Site_WGS = gpd.GeoDataFrame(Site, geometry=gpd.points_from_xy(Site.longitude, Site.latitude), crs="EPSG:4326")
        self.Site_UTM = self.Site_WGS.to_crs(self.EPSG)
        self.Site_UTM = self.Site_UTM.reset_index(drop=True)

        self.z = float(self.ini[self.Site_code]['zm'])

        # =====================================================================================
        # Define grid parameters for model
        # Domain is the upwind_fetch (m) in all directions
        # Will create a grid centered on [0 0 zm]
        self.domain = int(self.ini['FFP_Parameters']['upwind_fetch'])
        # dx Cell size of domain [m], Small dx results in higher spatial resolution and higher computing time
        self.dx = int(self.ini['FFP_Parameters']['resolution'])
        # Percentage of source area for which to provide contours, must be between 10% and 90%        
        self.rs = [float(rs) for rs in self.ini['FFP_Parameters']['rs'].split(',')]

        self.nx = int(self.domain*2 / self.dx)

        x = np.linspace(-self.domain, self.domain, self.nx)# + 1)
        self.x_2d, self.y_2d = np.meshgrid(x, x)

        # Polar coordinates
        # Set theta such that North is pointing upwards and angles increase clockwise
        self.rho = np.sqrt(self.x_2d**2 + self.y_2d**2)
        self.theta = np.arctan2(self.x_2d, self.y_2d)

        # Apply a symmetric mask to restrict summations to a radius of upwind_fetch around [0 0 zm]
        symetric_Mask = self.rho.copy()
        symetric_Mask[self.rho>self.domain] = np.nan
        self.symetric_Mask = symetric_Mask*0 + 1

        # initialize rasters for footprint climatology
        self.fclim_2d_empty = np.zeros(self.x_2d.shape)


        # basemap is an optional input, requires a 'path to vector layer' pluss a 'classification' key
        self.rasterizeBasemap(self.ini[self.Site_code]['basemap'],self.ini[self.Site_code]['basemap_class'])

        # ==================================
        # Define input keys

        self.read_Met()
        
    def read_Met(self):

        self.vars = self.ini['Input_Variable_Names']

        self.vars_metadata = self.ini['Input_Variable_Definitions']
        
        if self.ini['FFP_Parameters']['verbose'] == "True":
            print('Requires the following inputs, expecting them to be named as specified:\n')
            for key,value in self.vars.items():
                print(self.vars_metadata[key])
                print(f'Labelled as "{value}" in input self.dataset\n')

        self.Data = pd.read_csv(self.ini['Input_Files']['dpath'],
                 parse_dates=[self.ini['Input_Files']['timestamp']],
                 )

        # # Exclude timesteps with *No* Data
        # self.Data.dropna(how='all')

        if self.vars['canopy_height'] not in self.Data:
            self.Data['canopy_height']=self.Site_UTM['canopy_height'][0]
        else:
            self.Data[self.vars['canopy_height']] = self.Data[self.vars['canopy_height']].fillna(self.Site_UTM['canopy_height'][0])
        if self.vars['z0'] != '':
            self.Data['z0'] = self.Data[self.vars['z0']].copy()
        if self.ini['Assumptions']['roughness_length'] != 'None':
            self.Data['z0'] = self.Data[self.vars['canopy_height']]*float(self.ini['Assumptions']['roughness_length'])
        elif self.ini['Assumptions']['roughness_length'] == 'None':
            self.Data['z0'] = None
        self.Data['zm-d'] = self.z-(self.Data[self.vars['canopy_height']]*float(self.ini['Assumptions']['displacement_height']))
        self.Data['zm/ol'] = self.Data['zm-d']/self.Data[self.vars['ol']]

        if self.Subsets is not None:
            self.Data['Subset']=self.Data[self.Subsets]
            self.Data['Subset'] = self.Data['Subset'].fillna('N/A')
        else:
           self.Data['Subset'] = 'Climatology'
            
        self.Data[self.Fc_Names] = np.nan
        self.Subset_Climatology = {}
        for self.sub_name in self.Data['Subset'].unique():
            if self.sub_name != 'N/A':
                print(self.sub_name)
                self.run(sub_data=self.Data.loc[self.Data['Subset']==self.sub_name])
                for c in self.Fc_Names:
                    self.Data[c] = self.Data[c].fillna(self.sub_data[c])

                # Reverse the ordering of the array so that north is up
                self.Subset_Climatology[self.sub_name] = self.fclim_2d[::-1]
        
        self.summarizeClimatology()

    def rasterizeBasemap(self,basemap,basemap_class):
        x,y = self.Site_UTM.geometry.x[0],self.Site_UTM.geometry.y[0]
        west = x-(self.nx*self.dx)/2
        north = y+(self.nx*self.dx)/2
        self.Transform = from_origin(west,north,self.dx,self.dx)
        
        if os.path.isfile(basemap):
            print('Rasterizing bassemap')
            # Read basemap layer and reproject if not already in the proper WGS 1984 Zone
            self.baseVector = gpd.read_file(basemap).to_crs(self.EPSG)
            self.baseVector = gpd.clip(self.baseVector,self.Site_UTM.buffer(self.domain))
            if basemap_class != 'None' and basemap_class != '':
                self.baseVector = self.baseVector.dissolve(by=basemap_class).reset_index()
            else:
                self.baseVector = self.baseVector.dissolve().reset_index(drop=True)
                basemap_class = 'aoi'
                self.baseVector[basemap_class] = basemap_class
            self.baseVector.index+=1
            self.baseRasterKey = self.baseVector[basemap_class].to_dict()
            self.Fc_Names = [self.baseRasterKey[i]+'_Fc' for i in self.baseVector.index.values]

            shapes = ((geom,value) for geom,value in zip(self.baseVector['geometry'],self.baseVector.index))

            with rasterio.open(f"{self.ini['Output']['raster_output']}/Footprint_Basemap_{self.dx}m.tif",'w+',driver='GTiff',width = self.nx, height = self.nx,#+1,
                            count = 1,dtype=np.float32,transform = self.Transform,crs = ({'init': f'EPSG:{self.EPSG}'})) as out:
                out_arr = out.read(1)
                self.baseRaster = features.rasterize(shapes=shapes,fill = 100,out = out_arr,transform = self.Transform,default_value=-1)
                self.baseRaster = self.baseRaster * self.symetric_Mask
                out.write(self.baseRaster,1)
        else: 
            print('Basemap not provided, creating default')
            self.baseRaster = self.symetric_Mask
            self.Fc_Names = []
            self.baseRasterKey = {f'Contribution within {self.domain} m':''}
        
    def run(self,sub_data):         
        self.sub_data = sub_data.copy()
        self.fclim_2d = self.fclim_2d_empty.copy()

        self.Filter()
        print(f"Processing: {self.sub_data.loc[self.sub_data['process']==1].shape[0]} out of {self.sub_data.shape[0]} input records")
        self.sub_data = self.sub_data.loc[self.sub_data['process'] == 1]

        if (__name__ == 'FFP_Asssment' or __name__ == '__main__') and int(self.ini['Multi_Processing']['processes'])>1:
            
            batchsize=int(np.ceil(self.sub_data.shape[0]))
            if batchsize > int(self.ini['Multi_Processing']['batchSize']):
                batchsize = int(self.ini['Multi_Processing']['batchSize'])

            ix = 0
            while self.sub_data[ix:ix+batchsize].shape[0]>0:
                batch = self.sub_data[ix:ix+batchsize].copy()
                print(f'Processing Batch {ix}:{ix+batchsize}')
                umean = batch[self.vars['umean']]
                ustar = batch[self.vars['ustar']]
                sigmav = batch[self.vars['sigmav']]
                h = batch[self.vars['h']]
                ol = batch[self.vars['ol']]
                wind_dir = batch[self.vars['wind_dir']]
                z0 = batch['z0']
                zm = batch['zm-d']
                index = batch.index
                ix += batchsize
                
                with Pool(processes=int(self.ini['Multi_Processing']['processes'])) as pool:
                    for out in pool.starmap(partial(FFP,theta=self.theta,rho=self.rho,x_2d=self.x_2d,basemap=self.baseRaster),
                                        zip(index,umean,ustar,sigmav,h,ol,wind_dir,z0,zm)):
                        self.processOutputs(out)
                    pool.close()

        else:
            for i,row in self.sub_data.iterrows():
                print((i,row[self.vars['ustar']],row[self.vars['sigmav']],row[self.vars['h']],
                    row[self.vars['ol']],row[self.vars['wind_dir']],row['z0'],row['zm-d']))
                out = FFP(i,row[self.vars['umean']],row[self.vars['ustar']],row[self.vars['sigmav']],row[self.vars['h']],
                    row[self.vars['ol']],row[self.vars['wind_dir']],row['z0'],row['zm-d'],
                    self.theta,self.rho,self.x_2d,basemap=self.baseRaster)
                self.processOutputs(out)

    def Filter(self):
        d = int(self.ini['FFP_Parameters']['exclude_wake'])
        b = float(self.ini[self.Site_code]['bearing'])
        Exclude = {
            'under':{
                'zm/ol':-15.5,
                'zm-d':self.sub_data['z0']*12.5,
                self.vars['ustar']:.1,
                self.vars['sigmav']:0,
                self.vars['h']:10,
                self.vars['h']:self.sub_data['zm-d'],
                self.vars['wind_dir']:0,
            },
            'over':{
                'zm-d':self.sub_data[self.vars['h']],
                self.vars['wind_dir']:360
            },
            'between':{
                self.vars['wind_dir']:[b-180-d,b-180+d,b+180-d,b+180+d]
            }
        }
        self.sub_data['process'] = 1

        for key,val in self.vars.items():
            if val != '':
                self.sub_data.loc[self.sub_data[val].isna(),'process']=-1

        NA_flag = (self.sub_data.loc[self.sub_data['process']==-1].shape[0])
        print(f'{NA_flag} records flagged for missing data')

        for key,value in Exclude['under'].items():
            flagged = self.sub_data.loc[self.sub_data[key]<value].shape[0]
            if flagged > 0:
                print(f'{flagged} records flagged for low {key}')  
            self.sub_data.loc[self.sub_data[key]<value,'process']=0
            
        for key,value in Exclude['over'].items():
            flagged = self.sub_data.loc[self.sub_data[key]>value].shape[0]
            if flagged > 0:
                print(f'{flagged} records flagged for high {key}')
            self.sub_data.loc[self.sub_data[key]>value,'process']=0

        for key,value in Exclude['between'].items():
            flagged = self.sub_data.loc[(((self.sub_data[key]>value[0]) & (self.sub_data[key]<value[1]))|
                                        ((self.sub_data[key]>value[2]) & (self.sub_data[key]<value[3])))].shape[0]
            if flagged > 0:
                print(f'{flagged} records flagged for unacceptable {key}')
            self.sub_data.loc[(((self.sub_data[key]>value[0]) & (self.sub_data[key]<value[1]))|
                                        ((self.sub_data[key]>value[2]) & (self.sub_data[key]<value[3]))),'process']=0

    def processOutputs(self,out):
        self.fclim_2d = self.fclim_2d + out[1] * self.symetric_Mask
        if len(out) >2 and len(self.Fc_Names) > 0:
            self.Data.loc[self.Data.index==out[0],self.Fc_Names]=out[2]
            self.Data.loc[self.Data.index==out[0],f'Contribution within {self.domain} m']=out[1].sum()
            self.Data.loc[self.Data.index==out[0],'process']=1
            
        else:
            self.Data.loc[self.Data.index==out[0],f'Contribution within {self.domain} m']=out[1].sum()
            self.Data.loc[self.Data.index==out[0],'process']=1

    def summarizeClimatology(self):
        self.fclim_2d_Full = self.fclim_2d_empty.copy()
        N = 0

        
        self.feature = {"type": "Feature", "properties": {}, "geometry": {}}
        self.FeatureCollection = {
            "type": "FeatureCollection",
            "features": []
        }


        with rasterio.open(f"{self.ini['Output']['raster_output']}{self.Site_code}_FFP_Clim_{self.dx}m.tif",'w+',driver='GTiff',
                           width = self.nx, height = self.nx,count = len(self.sub_name)+1,dtype=np.float32,
                           transform = self.Transform,crs = ({'init': f'EPSG:{self.EPSG}'})) as raster_out:
            for i,(self.sub_name,self.fclim_2d) in enumerate(self.Subset_Climatology.items()):
                self.fclim_2d_Full += self.fclim_2d.copy()
                n = self.Data.loc[((self.Data['Subset']==self.sub_name)&(self.Data['process']==1))].shape[0]
                N += n
                self.fclim_2d = self.fclim_2d/n
                self.contours_from_rio(i+1)
                raster_out.write(self.fclim_2d,i+1)
                # self.countours()
            # print('Adj')
            # print(N)
            # print(np.nansum(self.fclim_2d_Full))
            self.fclim_2d = self.fclim_2d_Full/N
            raster_out.write(self.fclim_2d,i+2)
            # print(np.nansum(self.fclim_2d))
            self.sub_name = 'Climatology'
            # self.countours()
            self.contours_from_rio(i+2)
            self.contour_levels = gpd.GeoDataFrame.from_features(self.FeatureCollection["features"],crs=self.EPSG)
            
            self.contour_levels.to_file(f"{self.ini['Output']['raster_output']}{self.Site_code}_FFP_Clim_{self.dx}m.shp")

            self.webMap()
            # self.gdf.geometry = gdf.geometry.simplify(FFP.dx/S_F).buffer(FFP.dx/S_F, join_style=1).buffer(FFP.dx/S_F, join_style=1)
            # gdf.geometry = gdf.geometry.buffer(FFP.dx, cap_style=3).buffer(-FFP.dx, cap_style=3)
            # self.gdf[self.gdf['r']>.5].plot(color='None',edgecolor='k')

    def contours_from_rio(self,i):
        fclim_2d_r = self.fclim_2d*0
        sf = np.sort(self.fclim_2d, axis=None)[::-1]
        msf = np.ma.masked_array(sf, mask=(np.isnan(sf) | np.isinf(sf))) 
        csf = msf.cumsum().filled(np.nan)

        self.rs.sort(reverse=True)
        for r in self.rs:
            sfv = sf[np.where(csf>=r)].max()
            fclim_2d_r[np.where(self.fclim_2d>=sfv)]=r

        fclim_2d_r[np.isnan(fclim_2d_r)]=0
        Mask = np.array(fclim_2d_r,dtype=bool)
        fclim_2d_r = np.float32(fclim_2d_r)

        shapes = features.shapes(fclim_2d_r,mask=Mask,transform=self.Transform)
        for s in shapes:
            self.feature["properties"] = {'r':s[1],'Interval':self.sub_name,'BandID':i}
            self.feature["geometry"] = s[0]
            self.FeatureCollection['features'].append(self.feature.copy())
    
    def webMap(self):
        self.WGS = self.contour_levels.to_crs('WGS1984')
        self.WGS['info'] = (self.WGS['r']*100).astype(int).astype(str)+ ' % Flux Source Area Contour<br>For: '+self.WGS['Interval'].astype(str)       
        self.WGS = self.WGS.sort_values(by='r',ascending=False)

        MapTemplate = open(self.ini['Input_Files']['MapTemplate'],'r')
        MapFmt = MapTemplate.read().replace('Tower_Coords',str(self.lon_lat))

        rep_VarList = 'var FFP_Contour_Levels = FFP_json'
        VarList = ''
        rep_SourceList = MapFmt.split('// Source_Template_Start')[-1].split('// Source_Template_End')[0]
        SourceList = ''
        rep_StyleList = MapFmt.split('// Style_Template_Start')[-1].split('// Style_Template_End')[0]
        StyleList = ''

        Interval_list = "["
        for interval in self.WGS['Interval'].unique():
            Interval_list += f"'{interval}',"
            VarList += rep_VarList.replace('FFP_Contour_Levels',interval+'_Contour_Levels').replace('FFP_json',self.WGS.loc[self.WGS['Interval']==interval].to_json())+'\n\n    '
            SourceList += rep_SourceList.replace('FFP_Contour_Levels',interval+'_Contour_Levels')+'\n\n    '
            Temp_Style = rep_StyleList.replace('FFP_Contour_Levels',interval+'_Contour_Levels')+'\n\n    '
            if interval == 'Climatology':
                Temp_Style = Temp_Style.replace("'visibility':'none'","'visibility':'visible'")
            StyleList += Temp_Style
        Interval_list += "]"

        MapFmt = MapFmt.replace("const FFP_Contour_Intervals = [];",f"const FFP_Contour_Intervals = {Interval_list};")
        MapFmt = MapFmt.replace(rep_VarList,VarList)
        MapFmt = MapFmt.replace(rep_SourceList,SourceList)
        MapFmt = MapFmt.replace(rep_StyleList,StyleList)
        MapFmt = MapFmt.replace('SITE_json',self.Site_WGS.to_json())
        with open(f"{self.ini['Output']['webmapoutput']}{self.Site_code}_FFP_Clim_{self.dx}m.html",'w+') as out:
            out.write(MapFmt)

        self.WGS.to_file(f"{self.ini['Output']['webmapoutput']}{self.Site_code}_FFP_Clim_{self.dx}m.geojson",driver='GeoJSON')

    # def countours(self):
    #     pclevs = np.empty(len(self.rs))
    #     pclevs[:] = np.nan
    #     ars = np.empty(len(self.rs))
    #     ars[:] = np.nan
        
    #     sf = np.sort(self.fclim_2d, axis=None)[::-1]
    #     msf = np.ma.masked_array(sf, mask=(np.isnan(sf) | np.isinf(sf))) 
        
    #     csf = msf.cumsum().filled(np.nan)

    #     for ix, r in enumerate(self.rs):
    #         dcsf = np.abs(csf - r)
    #         pclevs[ix] = sf[np.nanargmin(dcsf)]
    #         ars[ix] = csf[np.nanargmin(dcsf)]
    #     self.contour_levels = {'r':[],'r_true':[],'geometry':[]}

    #     for r, r_thresh, lev in zip(self.rs, ars, pclevs):
    #         geom = self.getGeom(lev)
    #         if geom is not None:
    #             self.contour_levels['r'].append(r)
    #             self.contour_levels['r_true'].append(r_thresh)
    #             self.contour_levels['geometry'].append(Polygon((geom)))
    
    #     self.contour_levels = gpd.GeoDataFrame(data = {'r':self.contour_levels['r'],'r_true':self.contour_levels['r_true']},geometry=self.contour_levels['geometry'],crs=self.EPSG)
        
    #     if self.ini['Output']['shapefile_output']!='None':
    #         self.contour_levels.to_file(f"{self.ini['Output']['raster_output']}{self.Site_code}_FFP_Clim_{self.dx}m_{self.sub_name.replace(':','')}.shp")
        
    #     if os.path.isdir(self.ini['Output']['webmapoutput']):
    #         self.WGS = self.contour_levels.to_crs('WGS1984')
    #         self.WGS['info'] = (self.WGS['r']*100).astype(int).astype(str)+ ' % Flux Source Area Contour'       
    #         self.WGS = self.WGS.sort_values(by='r',ascending=False)

    #         MapTemplate = open(self.ini['Input_Files']['MapTemplate'],'r')
    #         MapFmt = MapTemplate.read().replace('Tower_Coords',str(self.lon_lat))

    #         MapFmt = MapFmt.replace('FFP_json',self.WGS.to_json())
    #         MapFmt = MapFmt.replace('SITE_json',self.Site_WGS.to_json())
    #         with open(f"{self.ini['Output']['webmapoutput']}{self.Site_code}_FFP_Clim_{self.dx}m_{self.sub_name.replace(':','')}.html",'w+') as out:
    #             out.write(MapFmt)

    #         self.WGS.to_file(f"{self.ini['Output']['webmapoutput']}{self.Site_code}_FFP_Clim_{self.dx}m_{self.sub_name.replace(':','')}.geojson",driver='GeoJSON')

    # def getGeom(self,lev):
    #     cs = plt.contour(self.x_2d, self.y_2d, self.fclim_2d, [lev])
    #     plt.close()
    #     segs = cs.allsegs[0]#[0]

    #     print('Update implementation')
    #     print(cs.levels,len(cs.allsegs))
    #     print(self.sub_name)
    #     print()
    #     xr = [vert[0] for vert in segs]
    #     yr = [vert[1] for vert in segs]
    #     #Set contour to None if it's found to reach the physical domain
    #     if self.x_2d.min() >= min(segs[:, 0]) or max(segs[:, 0]) >= self.x_2d.max() or \
    #     self.y_2d.min() >= min(segs[:, 1]) or max(segs[:, 1]) >= self.y_2d.max():
    #         return None
    #     else:
    #         return([[x+self.Site_UTM.geometry.x[0], y+self.Site_UTM.geometry.y[0]] for x,y in zip(xr,yr)])

        