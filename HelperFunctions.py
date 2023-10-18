import sys
import time

class progressbar():

    def __init__(self,items,prefix,size=60,out=sys.stdout):
        self.items = items
        self.out = out
        self.i = 0
        self.prefix=prefix
        self.size=size
        self.T1 = time.time()
        self.show(0)

    def show(self,j):
        x = int(self.size*j/self.items)
        e = time.time()-self.T1
        print(f"{self.prefix}[{u'█'*x}{('.'*(self.size-x))}] {min(j,self.items)}/{self.items} {round(e,1)} Seconds elapsed", end='\r', file=self.out, flush=True)

    def step(self,step_size=1):
        self.i+=step_size
        self.show(self.i)

    def close(self):
        print()
