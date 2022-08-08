import numpy as np
import math

resolution=0.25
tmax=40
numpoints=int(math.floor(tmax/resolution))
a=np.linspace(resolution,tmax,num=numpoints,endpoint=True)
fp=open("timesteps.txt","w")
for i in range(len(a)):
    fp.write(str(a[i])+"\n")
fp.close()
