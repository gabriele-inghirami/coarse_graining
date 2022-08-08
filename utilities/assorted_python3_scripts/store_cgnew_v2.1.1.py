# store_cg_new.py - version 2.1.0 - 22/05/2020
# it reads the energy and number densities and returns the temperatures, with/without the anisotropic correction as in PhysRevC.83.034907
# we assume that, in addition to tensor_densities_* files, there are corresponding tensor_Tmunu_* files
# in this version, we save also the number of hadrons in each cell and the total baryon + antibaryon density

import fileinput
import math
import numpy as np
import sys
import os
import pickle
from scipy import interpolate

#number threshold of particles to accept a cell
num_thresh=20

#we get the name of input and output files
N_input_files=len(sys.argv)-1

if(N_input_files<2):
   print ('Syntax: ./store_cg.py <coarse file data 1> [coarse file data 2] ... <outputfile>')
   print ("coarse file data 1,2,3...N are the density files produced by the coarse graining code")
   print ("outputfile is obviously the name of the output file with the results of the postprocessing")
   sys.exit(1)

coarsefile=sys.argv[1:N_input_files]
outputfile=sys.argv[N_input_files]

Tmunufile=[]
for i in range(len(coarsefile)):
    Tmunufile.append(coarsefile[i].replace("densities","Tmunu"))

#these are the unit values of energy density and pressure
e0=0.14651751415742
#these are the unit values of the net baryon density and entropy density
n0=0.15891

datadir="EOS_HG_UrQMD/"

fstd=datadir+"hadgas_eos.dat"
Ne_std=2001
Nn_std=401
en_std_max=1000.
rho_std_max=40.
enarr_std=np.linspace(0.,en_std_max,num=Ne_std)
rhoarr_std=np.linspace(0.,rho_std_max,num=Nn_std)

fmed=datadir+"hg_eos_small.dat"
Ne_med=201
Nn_med=201
en_med_max=10.
rho_med_max=2.
enarr_med=np.linspace(0.,en_med_max,num=Ne_med)
rhoarr_med=np.linspace(0.,rho_med_max,num=Nn_med)

fmin=datadir+"hg_eos_mini.dat"
Ne_min=201
Nn_min=201
en_min_max=0.1
rho_min_max=0.02
enarr_min=np.linspace(0.,en_min_max,num=Ne_min)
rhoarr_min=np.linspace(0.,rho_min_max,num=Nn_min)

temparr_std=np.zeros((Ne_std,Nn_std),dtype=np.float64)
muarr_std=np.zeros((Ne_std,Nn_std),dtype=np.float64)
sarr_std=np.zeros((Ne_std,Nn_std),dtype=np.float64)
parr_std=np.zeros((Ne_std,Nn_std),dtype=np.float64)

temparr_med=np.zeros((Ne_med,Nn_med),dtype=np.float64)
muarr_med=np.zeros((Ne_med,Nn_med),dtype=np.float64)
sarr_med=np.zeros((Ne_med,Nn_med),dtype=np.float64)
parr_med=np.zeros((Ne_med,Nn_med),dtype=np.float64)

temparr_min=np.zeros((Ne_min,Nn_min),dtype=np.float64)
muarr_min=np.zeros((Ne_min,Nn_min),dtype=np.float64)
sarr_min=np.zeros((Ne_min,Nn_min),dtype=np.float64)
parr_min=np.zeros((Ne_min,Nn_min),dtype=np.float64)

def readeos(ff,tarr,marr,parr,sarr,nx,ny):
    for j in range(ny):
        for i in range(nx):
            stuff=ff.readline().split()
            tarr[i,j],marr[i,j],parr[i,j],sarr[i,j]=np.float64(stuff[0]),np.float64(stuff[1]),np.float64(stuff[3]),np.float64(stuff[5])

print("Reading the tabulated EoS from the files")

with open(fstd,"r") as infile:
     readeos(infile,temparr_std,muarr_std,parr_std,sarr_std,Ne_std,Nn_std)
     temp_interp_std=interpolate.interp2d(enarr_std, rhoarr_std, temparr_std.transpose(), kind='linear')
     muB_interp_std=interpolate.interp2d(enarr_std, rhoarr_std, muarr_std.transpose(), kind='linear')
     p_interp_std=interpolate.interp2d(enarr_std, rhoarr_std, parr_std.transpose(), kind='linear')
     s_interp_std=interpolate.interp2d(enarr_std, rhoarr_std, sarr_std.transpose(), kind='linear')

with open(fmed,"r") as infile:
     readeos(infile,temparr_med,muarr_med,parr_med,sarr_med,Ne_med,Nn_med)
     temp_interp_med=interpolate.interp2d(enarr_med, rhoarr_med, temparr_med.transpose(), kind='linear')
     muB_interp_med=interpolate.interp2d(enarr_med, rhoarr_med, muarr_med.transpose(), kind='linear')
     p_interp_med=interpolate.interp2d(enarr_med, rhoarr_med, parr_med.transpose(), kind='linear')
     s_interp_med=interpolate.interp2d(enarr_med, rhoarr_med, sarr_med.transpose(), kind='linear')

with open(fmin,"r") as infile:
     readeos(infile,temparr_min,muarr_min,parr_min,sarr_min,Ne_min,Nn_min)
     temp_interp_min=interpolate.interp2d(enarr_min, rhoarr_min, temparr_min.transpose(), kind='linear')
     muB_interp_min=interpolate.interp2d(enarr_min, rhoarr_min, muarr_min.transpose(), kind='linear')
     p_interp_min=interpolate.interp2d(enarr_min, rhoarr_min, parr_min.transpose(), kind='linear')
     s_interp_min=interpolate.interp2d(enarr_min, rhoarr_min, sarr_min.transpose(), kind='linear')

print("Done.\n")

def get_mub_T(rhoB_input_w_sign,edens):
     #before callin this function we already checked that both arguments are > 0
     compute=True
     rhoB=np.abs(rhoB_input_w_sign)
     edens=edens/e0
     rhoB=rhoB/n0
     if(edens<=en_std_max):
       if((edens<en_min_max ) and (rhoB<rho_min_max)):
         ftemp=temp_interp_min
         fmuB=muB_interp_min
         fpr=p_interp_min
         fs=s_interp_min
       if(edens<en_med_max ) and (rhoB<rho_med_max) and ((edens>=en_min_max ) or (rhoB>=rho_min_max)):     
         ftemp=temp_interp_med
         fmuB=muB_interp_med
         fpr=p_interp_med
         fs=s_interp_med
       if((edens>=en_med_max ) or (rhoB>=rho_med_max)):     
         if(rhoB>rho_std_max):
           print("Net baryon density exceeding the maximum of the table. Changed from "+str(rhoB)+" to "+str(rho_std_max*0.999999))
           rhoB=rho_std_max*0.999999
         ftemp=temp_interp_std
         fmuB=muB_interp_std
         fpr=p_interp_std
         fs=s_interp_std
     elif(edens>en_std_max):    
         compute=False
#         temperature=350./1000.
#         muB=3./1000.
         temperature=0.
         muB=0.
         pressure=0.
         entr_dens=0.
     else:        
         compute=False
         temperature=0.
         muB=0.
         pressure=0.
         entr_dens=0.
     
     if(compute): #all is expressed in GeV
       temperature=ftemp(edens,rhoB)[0]/1000.
       muB=3*fmuB(edens,rhoB)[0]/1000. 
       pressure=fpr(edens,rhoB)[0]*e0
       entr_dens=fs(edens,rhoB)[0]*n0

     return muB, temperature, pressure, entr_dens


#the number of time intervals should be equal to the number of coarse graining datafiles
nt=len(coarsefile)
tt=np.zeros(nt,dtype=np.float64)

#auxiliary funcion kronecker's delta
def kron(ii,jj):
    if(ii==jj):
      return 1
    return 0

print("Reading coarse data files... ")
#we open the coarse graining data files
first_time=True
for ff in range(nt):
    if(first_time):
      first_time=False
      cdata = open(coarsefile[ff],"rb")
      print("Extracting data from "+coarsefile[ff])
      nevents=np.fromfile(cdata,dtype=np.int64,count=1)[0]
      print("Number of events: "+str(nevents))
      tt[ff]=np.fromfile(cdata,dtype=np.float64,count=1)[0] 
      print("Time array index: "+str(ff)+", with time value: "+str(tt[ff]))
      num_part=np.fromfile(cdata,dtype=np.int32,count=1)[0]
      print("Number of hadrons in the coarse graining data: "+str(num_part))
      nx,ny,nz=np.fromfile(cdata,dtype=np.int32,count=3)
      print("Grid size: "+str(nx)+", "+str(ny)+", "+str(nz))
      dx,dy,dz=np.fromfile(cdata,dtype=np.float64,count=3) 
      print("Grid resolution: "+str(dx)+", "+str(dy)+", "+str(dz))
      datas=np.fromfile(cdata,dtype=np.float64,count=nx*ny*nz*(15+3*num_part)).reshape([nx,ny,nz,15+3*num_part])
      cdata.close()

      #we calculate the extrema of the grid, assuming that it symmetric
      lx=nx*dx
      xstart=-lx/2.+dx/2.
      xend=-xstart
      ly=ny*dy
      ystart=-ly/2.+dy/2.
      yend=-ystart
      lz=nz*dz
      zstart=-lz/2.+dz/2.
      zend=-zstart

      #we get the extrema of the grid and the position of the cells
      xx=np.linspace(xstart,xend,num=nx)
      yy=np.linspace(ystart,yend,num=ny)
      zz=np.linspace(zstart,zend,num=nz)

      #arrays to store coarse graining data
      print("Allocating and initializing arrays... ")
      vx=np.zeros((nt,nx,ny,nz),dtype=np.float64)
      vy=np.zeros((nt,nx,ny,nz),dtype=np.float64)
      vz=np.zeros((nt,nx,ny,nz),dtype=np.float64)
      en=np.zeros((nt,nx,ny,nz),dtype=np.float64)
      rho=np.zeros((nt,nx,ny,nz),dtype=np.float64)
      temp=np.zeros((nt,nx,ny,nz),dtype=np.float64)
      mub=np.zeros((nt,nx,ny,nz),dtype=np.float64)
      tempA=np.zeros((nt,nx,ny,nz),dtype=np.float64)
      mubA=np.zeros((nt,nx,ny,nz),dtype=np.float64)
      ptra=np.zeros((nt,nx,ny,nz),dtype=np.float64)
      ppar=np.zeros((nt,nx,ny,nz),dtype=np.float64)
      s_entr_dens=np.zeros((nt,nx,ny,nz),dtype=np.float64)
      press=np.zeros((nt,nx,ny,nz),dtype=np.float64)
      s_entr_densA=np.zeros((nt,nx,ny,nz),dtype=np.float64)
      pressA=np.zeros((nt,nx,ny,nz),dtype=np.float64)
      Tmunu=np.zeros((4,4),dtype=np.float64)
      Lambda=np.zeros((4,4),dtype=np.float64)
      LambdaF=np.zeros((4,4),dtype=np.float64)
      uvel=np.zeros(4,dtype=np.float64)
      dens_hadrons=np.zeros((nt,nx,ny,nz),dtype=np.float64)
      rho_bab=np.zeros((nt,nx,ny,nz),dtype=np.float64) #total density of baryons and anti-baryons

      print("Done.")
    else:
      cdata = open(coarsefile[ff],"rb")
      print("Extracting data from "+coarsefile[ff])
      nevents_new=np.fromfile(cdata,dtype=np.int64,count=1)[0]
      if(nevents_new != nevents):
        print("WARNING: different number of events!!! Expected: "+str(nevents)+", read now: "+str(nevents_new))
      tt[ff]=np.fromfile(cdata,dtype=np.float64,count=1)[0]
      print("Time array index: "+str(ff)+", with time value: "+str(tt[ff]))
      num_part_new=np.fromfile(cdata,dtype=np.int32,count=1)[0]
      if(num_part_new != num_part):
        print("FATAL ERROR: different number of hadrons!!! Expected: "+str(num_part)+", read now: "+str(num_part_new))
        sys.exit(2)
      nx_new,ny_new,nz_new=np.fromfile(cdata,dtype=np.int32,count=3)
      if((nx_new != nx) or (ny_new != ny) or (nz_new != nz)):
        print("FATAL ERROR: different grid structure!!! Expected: "+str(nx)+", "+str(ny)+", "+str(nz)+", read now: "+str(nx_new)+", "+str(ny_new)+", "+str(nz_new))
        sys.exit(2)
      dx_new,dy_new,dz_new=np.fromfile(cdata,dtype=np.float64,count=3) 
      if((dx_new != dx) or (dy_new != dy) or (dz_new != dz)):
        print("FATAL ERROR: different grid resolution!!! Expected: "+str(dx)+", "+str(dy)+", "+str(dz)+", read now: "+str(dx_new)+", "+str(dy_new)+", "+str(dz_new))
        sys.exit(2)
      datas=np.fromfile(cdata,dtype=np.float64,count=nx*ny*nz*(15+3*num_part)).reshape([nx,ny,nz,15+3*num_part])
      cdata.close()

    #now we read Tmunu
    Tdata = open(Tmunufile[ff],"rb")
    print("Extracting data from "+Tmunufile[ff])
    nevents_T=np.fromfile(Tdata,dtype=np.int64,count=1)[0]
    if(nevents != nevents_T):
      print("WARNING: different number of events in "+Tmunufile[ff]+" ("+str(nevents_T)+") and "+coarsefile[ff]+" ("+str(nevents)+")")
    time_Tmunu=np.fromfile(Tdata,dtype=np.float64,count=1)[0]
    if(time_Tmunu != tt[ff]): 
      print("FATAL ERROR: different times in "+Tmunufile[ff]+" ("+str(time_Tmunu)+") and "+coarsefile[ff]+" ("+str(tt[ff])+")")
      sys.exit(2) 
    num_part_Tmunu=np.fromfile(Tdata,dtype=np.int32,count=1)[0]
    if(num_part_Tmunu != num_part):
       print("FATAL ERROR: different number of hadrons!!! "+Tmunufile[ff]+" ("+str(num_part_Tmunu)+") and "+coarsefile[ff]+" ("+str(num_part)+")")
       sys.exit(2)
    nx_T,ny_T,nz_T=np.fromfile(Tdata,dtype=np.int32,count=3)
    if((nx_T != nx) or (ny_T != ny) or (nz_T != nz)):
      print("FATAL ERROR: different grid structure between energy-momentum and densities tensors!!! Expected: "+str(nx)+", "+str(ny)+", "+str(nz)+", read now: "+str(nx_T)+", "+str(ny_T)+", "+str(nz_T))
      sys.exit(2)
    dx_T,dy_T,dz_T=np.fromfile(Tdata,dtype=np.float64,count=3) 
    if((dx_T != dx) or (dy_T != dy) or (dz_T != dz)):
      print("FATAL ERROR: different grid resolution between energy-momentum and densities tensors!!! Expected: "+str(dx)+", "+str(dy)+", "+str(dz)+", read now: "+str(dx_T)+", "+str(dy_T)+", "+str(dz_T))
      sys.exit(2)
    Tp=np.fromfile(Tdata,dtype=np.float64,count=nx*ny*nz*num_part*10).reshape([nx,ny,nz,num_part,10])
    Jb=np.fromfile(Tdata,dtype=np.float64,count=nx*ny*nz*4,offset=nx*ny*nz*num_part*4*8).reshape([nx,ny,nz,4])
     #regarding the line above: offset option available only with numpy>=1.17, expressed in bytes (so, with float64, we must multiply by 8)
    Tdata.close()    

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                glf_arg=Jb[i,j,k,0]*Jb[i,j,k,0]-Jb[i,j,k,1]*Jb[i,j,k,1]-Jb[i,j,k,2]*Jb[i,j,k,2]-Jb[i,j,k,3]*Jb[i,j,k,3]
                if(glf_arg<=0):
                  continue
                glf=np.sqrt(glf_arg)
                uvel[:]=Jb[i,j,k,:]/glf
                Lambda[0,0]=uvel[0]
                Lambda[0,1:]=-uvel[1:]
                Lambda[1:,0]=-uvel[1:]
                for aa in range(1,4):
                    for bb in range(1,4):
                        Lambda[aa,bb]=kron(aa,bb)+uvel[aa]*uvel[bb]/(1+uvel[0])
                LambdaF[0,0]=uvel[0]
                LambdaF[0,1:]=uvel[1:]
                LambdaF[1:,0]=uvel[1:]
                for aa in range(1,4):
                    for bb in range(1,4):
                        LambdaF[aa,bb]=kron(aa,bb)+uvel[aa]*uvel[bb]/(1+uvel[0])
                Tmunu[0,0]=np.sum(Tp[i,j,k,:,0])
                for aa in range(1,4):
                    Tmunu[0,aa]=np.sum(Tp[i,j,k,:,aa])
                    Tmunu[aa,0]=Tmunu[0,aa]
                Tmunu[1,1]=np.sum(Tp[i,j,k,:,4])
                Tmunu[1,2]=np.sum(Tp[i,j,k,:,5])
                Tmunu[2,1]=Tmunu[1,2]
                Tmunu[1,3]=np.sum(Tp[i,j,k,:,6])
                Tmunu[3,1]=Tmunu[1,3]
                Tmunu[2,2]=np.sum(Tp[i,j,k,:,7])
                Tmunu[2,3]=np.sum(Tp[i,j,k,:,8])
                Tmunu[3,2]=Tmunu[2,3]
                Tmunu[3,3]=np.sum(Tp[i,j,k,:,9])
                Tmunu_lrf=np.matmul(Lambda.transpose(),np.matmul(Tmunu,Lambda)) 
                #Tmunu_F=np.matmul(LambdaF.transpose(),np.matmul(Tmunu_lrf,LambdaF))
                #print("\n"+str(i)+"  "+str(j)+"  "+str(k))
                #print("boosted")
                #print(str(Tmunu_lrf[:,:]))
                #print("original")
                #print(str(Tmunu[:,:]))
                #print("reboosted")
                #print(str(Tmunu_F[:,:]))
                pressure_transverse=(Tmunu_lrf[1,1]+Tmunu_lrf[2,2])/2.
                #if(pressure_transverse<=0):
                #  print("LOOK HERE")
                #continue
                pressure_parallel=Tmunu_lrf[3,3]
                if((pressure_transverse<=0) or (pressure_parallel<=0)):
                #  print("Warning, negative or null pressure_transverse in "+str(i)+",  "+str(j)+",  "+str(k)+" : "+str(pressure_transverse))
                #  print(str(Tmunu[:,:]))
                #  print(str(Tmunu_lrf[:,:]))
                #  print(str(uvel[:]))
                  continue
                aniso_ratio=pressure_transverse/pressure_parallel
                #print("aniso_ratio: "+str(aniso_ratio))


                rho[ff,i,j,k]=datas[i,j,k,0]
                vx[ff,i,j,k]=datas[i,j,k,1]
                vy[ff,i,j,k]=datas[i,j,k,2]
                vz[ff,i,j,k]=datas[i,j,k,3]
                rho_bab[ff,i,j,k]=datas[i,j,k,4]
                edens=0
                num_hadrons=0
                for p in range(num_part):
                    num_hadrons=datas[i,j,k,15+3*p]+num_hadrons
                    dens_hadrons[ff,i,j,k]=datas[i,j,k,16+3*p]+dens_hadrons[ff,i,j,k]
                    edens=datas[i,j,k,17+3*p]+edens

                x_aniso=aniso_ratio**(4./3.)
                if(x_aniso<1):
                   r_aniso=1/2.*x_aniso**(-1./3.)*(1.+(x_aniso*np.arctanh(np.sqrt(1.-x_aniso)))/(np.sqrt(1.-x_aniso)))
                elif(x_aniso>1):
                   r_aniso=1/2.*x_aniso**(-1./3.)*(1.+(x_aniso*np.arctan(np.sqrt(x_aniso-1.)))/(np.sqrt(x_aniso-1.)))
                else: #x_aniso==1
                   r_aniso=1.
                edens_aniso=edens/r_aniso
                  
                en[ff,i,j,k]=edens

                if((en[ff,i,j,k]>0) and (rho[ff,i,j,k]>0) and (num_hadrons>num_thresh)):
                  mub[ff,i,j,k], temp[ff,i,j,k] , press[ff,i,j,k], s_entr_dens[ff,i,j,k] = get_mub_T(rho[ff,i,j,k],en[ff,i,j,k])
                  mubA[ff,i,j,k], tempA[ff,i,j,k] , pressA[ff,i,j,k], s_entr_densA[ff,i,j,k] = get_mub_T(rho[ff,i,j,k],edens_aniso)
                  #print("Eos results: "+str(rho[ff,i,j,k])+", "+str(en[ff,i,j,k])+", "+str(mub[ff,i,j,k])+", "+str(temp[ff,i,j,k]))

                ptra[ff,i,j,k]=pressure_transverse                   
                ppar[ff,i,j,k]=pressure_parallel
  
print("Done.")

print("Pickling tt,xx,yy,zz,vx,vy,vz,rho,en,mub,temp,muA,tempA,ptra,ppar,dens_hadrons,press,s_entr_dens,pressA,s_entr_densA,rho_bab")
with open(outputfile,"wb") as po:
     pickle.dump((tt[0:nt],xx,yy,zz,vx[0:nt,:,:,:],vy[0:nt,:,:,:],vz[0:nt,:,:,:],rho[0:nt,:,:,:],en[0:nt,:,:,:],mub[0:nt,:,:,:],temp[0:nt,:,:,:],mubA[0:nt,:,:,:],tempA[0:nt,:,:,:],ptra[0:nt,:,:,:],ppar[0:nt,:,:,:],dens_hadrons[0:nt,:,:,:],press[0:nt,:,:,:],s_entr_dens[0:nt,:,:,:],pressA[0:nt,:,:,:],s_entr_densA[0:nt,:,:,:],rho_bab[0:nt,:,:,:]),po)
print("All done.")



