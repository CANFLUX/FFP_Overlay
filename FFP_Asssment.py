# Wrapper for the Klujn et al. 2015 flux footprint model
import re
import os
import utm_zone
import numpy as np
import pandas as pd
import configparser
import geopandas as gpd
from functools import partial
from multiprocessing import Pool
from Klujn_2015_Model import FFP
from collections import defaultdict

import rasterio
from rasterio import features
from rasterio.transform import from_origin

from HelperFunctions import progressbar

import random

import matplotlib.pyplot as plt

class RunClimatology():

    def __init__(self,Site_code,Site_Config,Subsets=None,Basemap=None,Basemap_Class=None):
        self.Subsets = Subsets
            
        if isinstance(Site_code, str):
            Site_code=[Site_code]
        self.Site_code = Site_code
        self.ini = configparser.ConfigParser()
        self.ini.read('../MicrometPy.ini')
        self.ini.read('configuration.ini')
        self.ini.read(Site_Config)
        for sc in self.Site_code:
            if Basemap is not None and Basemap_Class is not None:
                self.ini[sc]['basemap'] = Basemap
                self.ini[sc]['basemap_class'] = Basemap_Class
            else:
                self.ini[sc]['basemap'] = 'N/A'
                self.ini[sc]['basemap_class'] = 'N/A'
        
        # Dump Site to a Dataframe
        # Site = pd.DataFrame(data=dict(self.ini[self.Site_code]),index=[0])
        Site={sc:dict(self.ini[sc]) for i,sc in enumerate(self.Site_code)}
        Site = pd.DataFrame(Site).T
        Site.index=[sc for sc in self.Site_code]

        for c in Site.columns:
            try:
                Site[c]=Site[c].astype('float64')
            except:
                pass
        # print
        # Default centerpoint to first site provided
        self.lon_lat = list(Site[['longitude','latitude']].values[0])
        self.EPSG = utm_zone.epsg(self.lon_lat)
        
        self.Site_WGS = gpd.GeoDataFrame(Site, geometry=gpd.points_from_xy(Site.longitude, Site.latitude), crs="EPSG:4326")
        self.Site_UTM = self.Site_WGS.to_crs(self.EPSG)
        # self.Site_UTM = self.Site_UTM.reset_index(drop=True)

        # self.z = float(self.ini[self.Site_code]['zm'])

        # =====================================================================================
        # Define grid parameters for model
        # Domain is the upwind_fetch (m) in all directions
        # Will create a grid centered on [0 0 zm]
        self.domain = int(self.ini['FFP_Parameters']['domain'])
        # dx Cell size of domain [m], Small dx results in higher spatial resolution and higher computing time
        self.dx = int(self.ini['FFP_Parameters']['dx'])
        # Percentage of source area for which to provide contours, must be between 10% and 90%        
        self.rs = [float(rs) for rs in self.ini['FFP_Parameters']['rs'].split(',')]

        self.nx = int(self.domain*2 / self.dx)

        x = np.linspace(-self.domain, self.domain, self.nx, dtype=np.float32)# + 1)
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
        self.fclim_2d_empty = np.zeros(self.x_2d.shape,dtype=np.float32)

        # basemap is an optional input, requires a 'path to vector layer' pluss a 'classification' key
        # update in future to allow for multiple basemaps
        self.rasterizeBasemap(self.ini[self.Site_code[0]]['basemap'],self.ini[self.Site_code[0]]['basemap_class'])

        self.FFP_Climatology = {}
        self.read_Met()
        if self.FFP_Climatology != {}:
            self.summarizeClimatology()
        self.Data.to_csv(f"{self.ini['Output']['dpath']}/FFP_Source_Area_{self.dx}m.csv")

    def rasterizeBasemap(self,basemap,basemap_class):
        x,y = self.Site_UTM.geometry.x.iloc[0],self.Site_UTM.geometry.y.iloc[0]
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
                basemap_class = 'Unit'
                self.baseVector[basemap_class]=self.baseVector.index
                self.baseVector = self.baseVector.dissolve().reset_index(drop=True)
            self.baseVector.index+=1
            self.baseRasterKey = self.baseVector[basemap_class].to_dict()
            self.Fc_Names = [self.baseRasterKey[i]+'_Fc' for i in self.baseVector.index.values]

            shapes = ((geom,value) for geom,value in zip(self.baseVector['geometry'],self.baseVector.index))

            with rasterio.open(f"{self.ini['Output']['dpath']}/Footprint_Basemap_{self.dx}m.tif",'w+',driver='GTiff',width = self.nx, height = self.nx,#+1,
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
        
        
    def read_Met(self):

        self.vars = self.ini['Input_Variable_Names']

        self.Data = pd.read_csv(self.ini['Input_Data']['dpath'],dtype=np.object_)

        if self.Subsets is not None:
            self.Data['Subset']=self.Data[self.Subsets]
            self.Data['Subset'] = self.Data['Subset'].fillna('N/A')
        else:
           self.Data['Subset'] = 'Climatology'

        self.Data = self.Data.merge(pd.Series(self.Site_code,name='Site_Code'),how='cross')
        self.Data.loc[self.Data['Subset']!='N/A','Subset']+=' '+self.Data.loc[self.Data['Subset']!='N/A','Site_Code']
        
        # Set FFP inputs to floats
        for key,val in self.vars.items():
            if val in self.Data.columns:
                self.Data[val]=self.Data[val].astype(np.float_)

        for static_variables in ['canopy_height','zm','bearing']:
            if self.vars[static_variables] not in self.Data:
                for sc in self.Site_code:
                    self.Data.loc[self.Data['Site_Code']==sc,self.vars[static_variables]]=self.Site_UTM.loc[sc,self.vars[static_variables]]
            else:
                for sc in self.Site_code:
                    self.Data.loc[self.Data['Site_Code']==sc,self.vars[static_variables]] = self.Data.loc[self.Data['Site_Code']==sc,self.vars[static_variables]].fillna(self.Site_UTM.loc[sc,static_variables])

        if self.vars['z0'] != '':
            self.Data['z0'] = self.Data[self.vars['z0']].copy()
        if self.ini['Assumptions']['roughness_length'] != 'None':
            self.Data['z0'] = self.Data[self.vars['canopy_height']]*float(self.ini['Assumptions']['roughness_length'])
        elif self.ini['Assumptions']['roughness_length'] == 'None':
            self.Data['z0'] = None
        self.Data['zm-d'] = self.Data['zm']-(self.Data[self.vars['canopy_height']]*float(self.ini['Assumptions']['displacement_height']))
        self.Data['zm/ol'] = self.Data['zm-d']/self.Data[self.vars['ol']]

            
        self.Data[self.Fc_Names] = np.nan
        for self.n_sub,self.sub_name in enumerate(self.Data['Subset'].unique()):
            if self.sub_name != 'N/A':
                print('\nProcessing: ',self.sub_name)
                self.sub_data = self.Data.loc[self.Data['Subset']==self.sub_name].copy()
                C = self.Filter()
                if C == True:
                    self.run()
                    for c in self.Fc_Names:
                        self.Data[c] = self.Data[c].fillna(self.sub_data[c])
                    sc = self.sub_data['Site_Code'].values[0]
                    self.FFP_Climatology[sc]={}
                    # Reverse the ordering of the array so that north is up
                    self.FFP_Climatology[sc][self.sub_name] = self.fclim_2d[::-1]
                else:
                    print(f'No valid input records for {self.sub_name}')

    def Filter(self):
        d = int(self.ini['FFP_Parameters']['exclude_wake'])
        # b = float(self.ini[self.Site_code]['bearing'])
        Exclude = {
            'under':{
                'zm/ol':-15.5,
                'zm-d':self.sub_data['z0']*12.5,
                self.vars['ustar']:.1,
                self.vars['sigmav']:0,
                self.vars['h']:10,
                # self.vars['h']:self.sub_data['zm-d'],
                self.vars['wind_dir']:0,
            },
            'over':{
                'zm-d':self.sub_data[self.vars['h']],
                self.vars['wind_dir']:360
            },
            # 'between':{
            #     self.vars['wind_dir']:[b-180-d,b-180+d,b+180-d,b+180+d]
            # }
        }
        self.sub_data['process'] = 1

        for key,val in self.vars.items():
            if val != '':
                self.sub_data.loc[self.sub_data[val].isna(),'process']=-1

        NA_flag = (self.sub_data.loc[self.sub_data['process']==-1].shape[0])
        if self.sub_data.loc[self.sub_data['process']==-1].shape[0]>0:
            print(f'{NA_flag} records skipped: missing data')

        if self.sub_data.loc[self.sub_data['process']==1].shape[0]>0:
            for key,value in Exclude['under'].items():
                flagged = self.sub_data.loc[self.sub_data[key]<value].shape[0]
                if flagged > 0:
                    print(f'{flagged} records skipped: low {key}')  
                self.sub_data.loc[self.sub_data[key]<value,'process']=0
                
            for key,value in Exclude['over'].items():
                flagged = self.sub_data.loc[self.sub_data[key]>value].shape[0]
                if flagged > 0:
                    print(f'{flagged} records skipped: high {key}')
                self.sub_data.loc[self.sub_data[key]>value,'process']=0

            # for key,value in Exclude['between'].items():
            
                
                # self.vars['wind_dir']:[b-180-d,b-180+d,b+180-d,b+180+d]
            flagged = self.sub_data.loc[(((self.sub_data['wind_dir']>self.sub_data['bearing']-180-d) & (self.sub_data['wind_dir']<self.sub_data['bearing']-180+d))|
                                        ((self.sub_data['wind_dir']>self.sub_data['bearing']+180-d) & (self.sub_data['wind_dir']<self.sub_data['bearing']+180+d))),'process'].shape[0]
            if flagged > 0:
                print(f'{flagged} records skipped: unacceptable {key}')
            self.sub_data.loc[(((self.sub_data['wind_dir']>self.sub_data['bearing']-180-d) & (self.sub_data['wind_dir']<self.sub_data['bearing']-180+d))|
                                            ((self.sub_data['wind_dir']>self.sub_data['bearing']+180-d) & (self.sub_data['wind_dir']<self.sub_data['bearing']+180+d))),'process']=0
        return(self.sub_data.loc[self.sub_data['process']==1].shape[0]>0)
        # else:
        #     return(False)
                
    def run(self):         
        print(f"Processing: {self.sub_data.loc[self.sub_data['process']==1].shape[0]} out of {self.sub_data.shape[0]} input records for {self.sub_name}")
        self.fclim_2d = self.fclim_2d_empty.copy()
        self.sub_data = self.sub_data.loc[self.sub_data['process'] == 1].copy()
        
        # Restricts multi-processing for small sub-sets
        n_processes = min(int(self.ini['Multi_Processing']['processes']),self.sub_data.shape[0])

        if (__name__ == 'FFP_Asssment' or __name__ == '__main__') and int(n_processes)>1:
            
            batchsize=min(self.sub_data.shape[0],int(self.ini['Multi_Processing']['batchSize']))

            ix = 0
            pb = progressbar(self.sub_data.shape[0],'Processing FFP')
            while self.sub_data[ix:ix+batchsize].shape[0]>0:
                batch = self.sub_data[ix:ix+batchsize].copy()
                # print(f' Batch {ix}:{ix+batchsize}')
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
                
                with Pool(processes=int(n_processes)) as pool:
                    for out in pool.starmap(partial(FFP,theta=self.theta,rho=self.rho,x_2d=self.x_2d,basemap=self.baseRaster),
                                        zip(index,umean,ustar,sigmav,h,ol,wind_dir,z0,zm)):
                        self.processOutputs(out)
                    pool.close()
                pb.step(batchsize)
            pb.close()

        else:
            for i,row in self.sub_data.iterrows():
                out = FFP(i,row[self.vars['umean']],row[self.vars['ustar']],row[self.vars['sigmav']],row[self.vars['h']],
                    row[self.vars['ol']],row[self.vars['wind_dir']],row['z0'],row['zm-d'],
                    self.theta,self.rho,self.x_2d,basemap=self.baseRaster)
                self.processOutputs(out)

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
        self.All_Contours = gpd.GeoDataFrame()
        for sc in self.Site_code:
            self.FeatureCollection = {
                "type": "FeatureCollection",
                "features": []
            }
            with rasterio.open(f"{self.ini['Output']['dpath']}{sc}_FFP_{self.dx}m.tif",'w+',driver='GTiff',
                            width = self.nx, height = self.nx,count = self.n_sub+2,dtype=np.float32,
                            transform = self.Transform,crs = ({'init': f'EPSG:{self.EPSG}'})) as raster_out:
                for i,(self.sub_name,self.fclim_2d) in enumerate(self.FFP_Climatology[sc].items()):
                    self.fclim_2d_Full += self.fclim_2d.copy()
                    n = self.Data.loc[((self.Data['Subset']==self.sub_name)&(self.Data['process']==1))].shape[0]
                    N += n
                    if len(self.FFP_Climatology[sc].items()) >0:
                        self.fclim_2d = self.fclim_2d/n
                        self.contours_from_rio(i+1)
                        raster_out.write(self.fclim_2d,i+1)
                if i > 0:
                    self.fclim_2d = self.fclim_2d_Full/N
                    self.contours_from_rio(i+2)
                    raster_out.write(self.fclim_2d,i+2)
                    self.sub_name = 'Climatology '+sc
                self.contour_levels = gpd.GeoDataFrame.from_features(self.FeatureCollection["features"],crs=self.EPSG)
                # Dissolve to get small "corner cells" merged into main shape
                self.contour_levels = self.contour_levels.dissolve(by=self.GDF_columns).reset_index()

                if self.ini['Output']['smoothing_factor']!='':
                    S_F = float(self.ini['Output']['smoothing_factor'])
                    self.contour_levels.geometry = self.contour_levels.geometry.simplify(self.dx*S_F).buffer(self.dx*S_F, join_style=1).buffer(self.dx*S_F, join_style=1)
                
                self.contour_levels = self.contour_levels.dissolve(by=self.GDF_columns).reset_index()

                for i,row in self.contour_levels.iterrows():
                    Dissolved = self.contour_levels.loc[((self.contour_levels['BandID'] == row['BandID'])&
                                                    (self.contour_levels['r'] <= row['r']))].dissolve().geometry
                    self.contour_levels.loc[self.contour_levels.index == i,'geometry'] = [Dissolved[0]]
            
                self.contour_levels['Area'] = self.contour_levels.area
                self.contour_levels['Site_Code']=sc
                self.All_Contours = gpd.GeoDataFrame(pd.concat([self.All_Contours,self.contour_levels], ignore_index=True), crs=self.contour_levels.crs)
            self.All_Contours.to_file(f"{self.ini['Output']['dpath']}{self.Site_code[0]}_FFP_{self.dx}m.shp")

            self.webMap()

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
        int_name = self.sub_name.replace(':','').replace('-','')
        n_Obs = self.Data.loc[self.Data['Subset']==self.sub_name].shape[0]
        int_name = re.sub('[^0-9a-zA-Z]+', '_', int_name)
        if int_name == 'Climatology':
            n_Obs = self.Data.loc[self.Data['Subset']!='N/A'].shape[0]
        if int_name[0].isdigit():
            int_name = 'I'+int_name
        for s in shapes:
            self.feature["properties"] = {'r':s[1],'Interval':int_name,'BandID':i,'n_Obs':n_Obs}
            self.feature["geometry"] = s[0]
            # print(self.feature)
            self.FeatureCollection['features'].append(self.feature.copy())
        self.GDF_columns=list(self.feature["properties"].keys())
    
    def webMap(self):
        self.WGS = self.All_Contours.to_crs('WGS1984')
        self.WGS['info'] = '<h3>'+self.WGS['Interval'].astype(str)+'</h3><br>'+np.round(self.WGS['r']*100).astype(int).astype(str)+ ' % FFP Contour<br>'+ self.WGS['Area'].apply(lambda n: f'{n:.2e}')+' m<sup>2</sup>'     
        self.WGS = self.WGS.sort_values(by='r',ascending=False)

        MapTemplate = open('MapTemplate.html','r')
        MapFmt = MapTemplate.read().replace('Tower_Coords',str(self.lon_lat))

        rep_VarList = 'var FFP_Contour_Levels = FFP_json'
        VarList = ''
        rep_SourceList = MapFmt.split('// Source_Template_Start')[-1].split('// Source_Template_End')[0]
        SourceList = ''
        rep_StyleList = MapFmt.split('// Style_Template_Start')[-1].split('// Style_Template_End')[0]
        StyleList = ''
        Interval_list = "["
        Sorted = self.WGS.groupby('Interval').first().sort_values(by='BandID')
        for interval in Sorted.index:
            Interval_list += f"'{interval}',"
            VarList += rep_VarList.replace('FFP_Contour_Levels',interval).replace('FFP_json',self.WGS.loc[self.WGS['Interval']==interval].to_json())+'\n\n    '
            SourceList += rep_SourceList.replace('FFP_Contour_Levels',interval)+'\n\n    '
            Temp_Style = rep_StyleList.replace('FFP_Contour_Levels',interval)+'\n\n    '
            if interval == 'Climatology':
                Temp_Style = Temp_Style.replace("'visibility':'none'","'visibility':'visible'")
                c = '#f76605'
            else:
                c = "#%06x" % random.randint(0, 0xFFFFFF)
            Temp_Style = Temp_Style.replace('Fill_Color',c)
            StyleList += Temp_Style
        Interval_list += "]"

        MapFmt = MapFmt.replace("const FFP_Contour_Intervals = [];",f"const FFP_Contour_Intervals = {Interval_list};")
        MapFmt = MapFmt.replace(rep_VarList,VarList)
        MapFmt = MapFmt.replace(rep_SourceList,SourceList)
        MapFmt = MapFmt.replace(rep_StyleList,StyleList)
        MapFmt = MapFmt.replace('SITE_json',self.Site_WGS.to_json())
        with open(f"{self.ini['Output']['dpath']}{self.Site_code[0]}_FFP_{self.dx}m.html",'w+') as out:
            out.write(MapFmt)

        self.WGS.to_file(f"{self.ini['Output']['dpath']}{self.Site_code[0]}_FFP_{self.dx}m.geojson",driver='GeoJSON')
