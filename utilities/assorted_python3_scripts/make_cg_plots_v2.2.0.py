# version 2.2.0 - 06/08/2020

# changelog: 2.2.0
# this version works with store_cg_2.1 and it plots also the entropy density
# changelog: 2.1.0
# this version works with store_cg_2.1 and it assumes that the entropy and the total baryon densities are also included in the input files
# changelog: 1.8.0
# this version works with store_cg_1.8 and it assumes that anisotropic corrections data are also included in the input files
# changelog: 1.6.1
# log plots for T and mu for animations
# changelog: 1.6.0
# this version works with store_cg_v1.6 and it assumes that density and energy density are also stored in the input files

import fileinput
import math
import numpy as np
import sys
import os
import pickle
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt

use_aniso=True #if True it uses the anisotropic correction

z_of_interest=[0.5]

time_min=0
time_max=30

#we get the name of input and output files
N_input_files=len(sys.argv)-1

if(N_input_files!=2):
   print ('Syntax: ./make_cg_plots.py <pickle binary data file> <output directory>')
   sys.exit(1)

inputfile=sys.argv[1]
od=sys.argv[2]
if(not os.path.exists(od)):
  os.mkdir(od)

with open(inputfile,"rb") as pi:
     data=pickle.load(pi)

tt,xx,yy,zz,vx,vy,vz,rho,en,muS,tempS,muA,tempA,ptra,ppar,dens_had,pressS,s_entr_densS,pressA,s_entr_densA,rho_bab=data[:]

hc3=(0.197326**3)

zz_selected=[] #it contains the indexes of zz correponding to the z_of_interest points
ftim=[] #it keeps track if it is the first time that the slice is plotted
vxmin_arr=[]
vymin_arr=[]
vzmin_arr=[]
enmin_arr=[]
rhomin_arr=[]
tempmin_arr=[]
mubmin_arr=[]
dens_hadmin_arr=[]
s_entr_densmin_arr=[]
rho_babmin_arr=[]
vxmax_arr=[]
vymax_arr=[]
vzmax_arr=[]
enmax_arr=[]
rhomax_arr=[]
tempmax_arr=[]
mubmax_arr=[]
dens_hadmax_arr=[]
s_entr_densmax_arr=[]
rho_babmax_arr=[]

nz=len(zz)
if(nz<2):
  zz_selected.append(0)
  ftim.append(True)
else:
  dz=zz[1]-zz[0]
  zmin=zz[0]-dz/2.
  zmax=zz[-1]+dz/2.
  for i in range(len(z_of_interest)):
    z_test=z_of_interest[i]
    if((z_test >= zmin) and (z_test<=zmax)):
      zz_selected.append(int(math.floor((z_test-zmin)/dz)))
      ftim.append(True)
 
nzsel=len(zz_selected)
if(nzsel==0):
  print("No selected points are in the available z-range")
  sys.exit(1)

if(use_aniso):
  print("Applying the anisotropic correction to the energy density distribution")
  mu=muA
  temp=tempA
  s_entr_dens=s_entr_densA
  for it in range(len(tt)):
      for iki in range(nzsel):
          ik=zz_selected[iki]
          for ii in range(len(xx)):
              for jj in range(len(yy)):
                  if(ppar[it,ii,jj,ik]==0):
                    en[it,ii,jj,ik]=0.
                  else:
                    aniso_ratio=ptra[it,ii,jj,ik]/ppar[it,ii,jj,ik]
                    x_aniso=aniso_ratio**(4./3.)
                    if(x_aniso<1):
                      r_aniso=1/2.*x_aniso**(-1./3.)*(1.+(x_aniso*np.arctanh(np.sqrt(1.-x_aniso)))/(np.sqrt(1.-x_aniso)))
                    elif(x_aniso>1):
                      r_aniso=1/2.*x_aniso**(-1./3.)*(1.+(x_aniso*np.arctan(np.sqrt(x_aniso-1.)))/(np.sqrt(x_aniso-1.)))
                    else: #x_aniso==1
                      r_aniso=1.
                    en[it,ii,jj,ik]=en[it,ii,jj,ik]/r_aniso
else:
  mu=muS
  temp=tempS
  s_entr_dens=s_entr_densS

fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 16
fig_size[1] = 8
plt.rcParams["figure.figsize"]=fig_size


for it in range(len(tt)):
  if(tt[it]<time_min):
    continue
  if(tt[it]>time_max):
    sys.exit(0)
  print("*****\n\nDoing timestep "+str(it)+", t="+'{:4.2f}'.format(tt[it]))
  for iki in range(nzsel):
    ik=zz_selected[iki]
    if(ftim[iki]==True):
      ftim[iki]==False
      testval=max(abs(np.amin(vx[:,:,:,ik])),abs(np.amax(vx[:,:,:,ik])))
      vxmin_arr.append(-testval)
      vxmax_arr.append(testval)
      testval=max(abs(np.amin(vy[:,:,:,ik])),abs(np.amax(vy[:,:,:,ik])))
      vymin_arr.append(-testval)
      vymax_arr.append(testval)
      testval=max(abs(np.amin(vz[:,:,:,ik])),abs(np.amax(vz[:,:,:,ik])))
      vzmin_arr.append(-testval)
      vzmax_arr.append(testval)
      enmin_arr.append(np.amin(en[:,:,:,ik]))
      rhomin_arr.append(np.amin(rho[:,:,:,ik]))
      mubmin_arr.append(np.amin(mu[:,:,:,ik]))
      tempmin_arr.append(np.amin(temp[:,:,:,ik]))
      enmax_arr.append(np.amax(en[:,:,:,ik]))
      rhomax_arr.append(np.amax(rho[:,:,:,ik]))
      mubmax_arr.append(np.amax(mu[:,:,:,ik]))
      tempmax_arr.append(np.amax(temp[:,:,:,ik]))
      dens_hadmax_arr.append(np.amax(dens_had[:,:,:,ik]))
      s_entr_densmax_arr.append(np.amax(s_entr_dens[:,:,:,ik]))
      rho_babmax_arr.append(np.amax(rho_bab[:,:,:,ik]))
      dens_hadmin_arr.append(np.amin(dens_had[:,:,:,ik]))
      s_entr_densmin_arr.append(np.amin(s_entr_dens[:,:,:,ik]))
      rho_babmin_arr.append(np.amin(rho_bab[:,:,:,ik]))

    if((np.amax(mu[it,:,:,ik])>0) or (np.amax(temp[it,:,:,ik])>0)):
      outdir=od+"/z_"+'{:+05.2f}'.format(zz[ik])
      if(not os.path.exists(outdir)):
        os.mkdir(outdir)
      print("Z slice: "+str(ik)+" ,z="+'{:4.2f}'.format(zz[ik]))

#     plots with local maxima and minima

      plt.suptitle("t="+'{:4.2f}'.format(tt[it])+" - z="+'{:4.2f}'.format(zz[ik]))

      maxvalue=np.amax(vx[it,:,:,ik])
      minvalue=np.amin(vx[it,:,:,ik])
      topval=max(abs(maxvalue),abs(minvalue))
      plt.subplot(241)
      plt.imshow(vx[it,:,:,ik].transpose(),extent=[xx[0], xx[-1], yy[0],yy[-1]], origin='lower',cmap='PiYG',vmin=-topval,vmax=topval)
      plt.title("V_x")
      plt.xlabel('x [fm]')
      plt.ylabel('y [fm]')
      ax=plt.gca()
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)
      plt.colorbar(label=r'$\mu$'+" [c units]",cax=cax,format="%6.3f")

      maxvalue=np.amax(vy[it,:,:,ik])
      minvalue=np.amin(vy[it,:,:,ik])
      topval=max(abs(maxvalue),abs(minvalue))
      plt.subplot(242)
      plt.imshow(vy[it,:,:,ik].transpose(),extent=[xx[0], xx[-1], yy[0],yy[-1]], origin='lower',cmap='PiYG',vmin=-topval,vmax=topval)
      plt.title("V_y")
      plt.xlabel('x [fm]')
      plt.ylabel('y [fm]')
      ax=plt.gca()
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)
      plt.colorbar(label=r'$\mu$'+" [c units]",cax=cax,format="%6.3f")

      maxvalue=np.amax(vz[it,:,:,ik])
      minvalue=np.amin(vz[it,:,:,ik])
      topval=max(abs(maxvalue),abs(minvalue))
      plt.subplot(243)
      plt.imshow(vz[it,:,:,ik].transpose(),extent=[xx[0], xx[-1], yy[0],yy[-1]], origin='lower',cmap='PiYG',vmin=-topval,vmax=topval)
      plt.title("V_z")
      plt.xlabel('x [fm]')
      plt.ylabel('y [fm]')
      ax=plt.gca()
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)
      plt.colorbar(label=r'$\mu$'+" [c units]",cax=cax,format="%6.3f")

      plt.subplot(245)
      plt.imshow(rho[it,:,:,ik].transpose()*1000,extent=[xx[0], xx[-1], yy[0],yy[-1]], origin='lower',cmap='gnuplot2')
      plt.title("Net baryon dens.")
      plt.xlabel('x [fm]')
      plt.ylabel('y [fm]')
      ax=plt.gca()
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)
      plt.colorbar(label=r'$\rho_B$'+" [1/fm^3]",cax=cax,format="%6.3f")

      plt.subplot(246)
      plt.imshow(en[it,:,:,ik].transpose()*1000,extent=[xx[0], xx[-1], yy[0],yy[-1]], origin='lower',cmap='gnuplot2')
      plt.title("Energy density")
      plt.xlabel('x [fm]')
      plt.ylabel('y [fm]')
      ax=plt.gca()
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)
      plt.colorbar(label=r'$\epsilon$'+" [MeV/fm^3]",cax=cax,format="%6.3f")

      plt.subplot(247)
      plt.imshow(mu[it,:,:,ik].transpose()*1000,extent=[xx[0], xx[-1], yy[0],yy[-1]], origin='lower',cmap='gnuplot2')
      plt.title("Baryon chem. pot.")
      plt.xlabel('x [fm]')
      plt.ylabel('y [fm]')
      ax=plt.gca()
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)
      plt.colorbar(label=r'$\mu$'+" [MeV]",cax=cax,format="%6.3f")

      plt.subplot(248)
      plt.imshow(temp[it,:,:,ik].transpose()*1000,extent=[xx[0], xx[-1], yy[0],yy[-1]], origin='lower',cmap='gnuplot2')
      plt.title("Temperature")
      plt.xlabel('x [fm]')
      plt.ylabel('y [fm]')
      ax=plt.gca()
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)
      plt.colorbar(label="T [MeV]",cax=cax,format="%6.3f")


      plt.tight_layout()

      plt.savefig(outdir+"/"+"t_"+'{:05.2f}'.format(tt[it])+".png",dpi=150,pad_inches=0.)
      plt.close('all')
    
#     plots with global maxima and minima, suitable to make animations

      plt.suptitle("t="+'{:4.2f}'.format(tt[it])+" - z="+'{:4.2f}'.format(zz[ik]))

      plt.subplot(241)
      plt.imshow(vx[it,:,:,ik].transpose(),extent=[xx[0], xx[-1], yy[0],yy[-1]], origin='lower',cmap='PiYG',vmin=vxmin_arr[iki],vmax=vxmax_arr[iki])
      plt.title("V_x")
      plt.xlabel('x [fm]')
      plt.ylabel('y [fm]')
      ax=plt.gca()
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)
      plt.colorbar(label=r'$\mu$'+" [c units]",cax=cax,format="%6.3f")

      plt.subplot(242)
      plt.imshow(vy[it,:,:,ik].transpose(),extent=[xx[0], xx[-1], yy[0],yy[-1]], origin='lower',cmap='PiYG',vmin=vymin_arr[iki],vmax=vymax_arr[iki])
      plt.title("V_y")
      plt.xlabel('x [fm]')
      plt.ylabel('y [fm]')
      ax=plt.gca()
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)
      plt.colorbar(label=r'$\mu$'+" [c units]",cax=cax,format="%6.3f")

      plt.subplot(243)
      plt.imshow(vz[it,:,:,ik].transpose(),extent=[xx[0], xx[-1], yy[0],yy[-1]], origin='lower',cmap='PiYG',vmin=vzmin_arr[iki],vmax=vzmax_arr[iki])
      plt.title("V_z")
      plt.xlabel('x [fm]')
      plt.ylabel('y [fm]')
      ax=plt.gca()
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)
      plt.colorbar(label=r'$\mu$'+" [c units]",cax=cax,format="%6.3f")

      plt.subplot(245)
      plt.imshow(np.log10(rho[it,:,:,ik].transpose()*1000+0.00001),extent=[xx[0], xx[-1], yy[0],yy[-1]], origin='lower',cmap='gnuplot2',vmin=-3,vmax=np.log10(rhomax_arr[iki]*1000))
      plt.title("Net baryon dens.")
      plt.xlabel('x [fm]')
      plt.ylabel('y [fm]')
      ax=plt.gca()
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)
      plt.colorbar(label=r'$\rho_B$'+" [1/fm^3]",cax=cax,format="%6.3f")

      plt.subplot(246)
      plt.imshow(np.log10(en[it,:,:,ik].transpose()*1000+0.00001),extent=[xx[0], xx[-1], yy[0],yy[-1]], origin='lower',cmap='gnuplot2',vmin=-2,vmax=np.log10(enmax_arr[iki]*1000))
      plt.title("Energy density")
      plt.xlabel('x [fm]')
      plt.ylabel('y [fm]')
      ax=plt.gca()
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)
      plt.colorbar(label=r'$\epsilon$'+" [MeV/fm^3]",cax=cax,format="%6.3f")

      plt.subplot(247)
      plt.imshow(np.log10(mu[it,:,:,ik].transpose()*1000+0.00001),extent=[xx[0], xx[-1], yy[0],yy[-1]], origin='lower',cmap='gnuplot2',vmin=0,vmax=3.)
      plt.title("Baryon chem. pot.")
      plt.xlabel('x [fm]')
      plt.ylabel('y [fm]')
      ax=plt.gca()
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)
      plt.colorbar(label="Log_{10}("+r'$\mu$'+" [MeV])",cax=cax,format="%6.3f")

      plt.subplot(248)
      plt.imshow(np.log10(temp[it,:,:,ik].transpose()*1000+0.00001),extent=[xx[0], xx[-1], yy[0],yy[-1]], origin='lower',cmap='gnuplot2',vmin=0,vmax=2.5)
      plt.title("Temperature")
      plt.xlabel('x [fm]')
      plt.ylabel('y [fm]')
      ax=plt.gca()
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)
      plt.colorbar(label="Log_{10}(T [MeV])",cax=cax,format="%6.3f")


      plt.tight_layout()

      plt.savefig(outdir+"/"+"t_"+'{:05.2f}'.format(tt[it])+"_global.png",dpi=150,pad_inches=0.)
      plt.close('all')

#     plots with local maxima and minima of e, n, e/n, T, s, s/T^3, nb, nbab

      plt.suptitle("t="+'{:4.2f}'.format(tt[it])+" - z="+'{:4.2f}'.format(zz[ik]))

      plt.subplot(241)
      plt.imshow(en[it,:,:,ik].transpose()*1000,extent=[xx[0], xx[-1], yy[0],yy[-1]], origin='lower',cmap='gnuplot2')
      plt.title("Energy density")
      plt.xlabel('x [fm]')
      plt.ylabel('y [fm]')
      ax=plt.gca()
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)
      plt.colorbar(label=r'$\epsilon$'+" [MeV/fm^3]",cax=cax,format="%6.3f")

      plt.subplot(242)
      plt.imshow(s_entr_dens[it,:,:,ik].transpose(),extent=[xx[0], xx[-1], yy[0],yy[-1]], origin='lower',cmap='gnuplot2')
      plt.title("Entropy density")
      plt.xlabel('x [fm]')
      plt.ylabel('y [fm]')
      ax=plt.gca()
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)
      plt.colorbar(label=r"s [1/fm^3]",cax=cax,format="%6.3f")

      plt.subplot(244)
      plt.imshow(temp[it,:,:,ik].transpose()*1000,extent=[xx[0], xx[-1], yy[0],yy[-1]], origin='lower',cmap='gnuplot2',vmin=0)
      plt.title("Temperature")
      plt.xlabel('x [fm]')
      plt.ylabel('y [fm]')
      ax=plt.gca()
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)
      plt.colorbar(label="T [MeV]",cax=cax,format="%6.3f")

      plt.subplot(243)
      plt.imshow(dens_had[it,:,:,ik].transpose(),extent=[xx[0], xx[-1], yy[0],yy[-1]], origin='lower',cmap='gnuplot2')
      plt.title("Particle density")
      plt.xlabel('x [fm]')
      plt.ylabel('y [fm]')
      ax=plt.gca()
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)
      plt.colorbar(label=r"n [1/fm^3]",cax=cax,format="%6.3f")

      plt.subplot(245)
      plt.imshow(rho[it,:,:,ik].transpose(),extent=[xx[0], xx[-1], yy[0],yy[-1]], origin='lower',cmap='gnuplot2')
      plt.title("Net baryon dens.")
      plt.xlabel('x [fm]')
      plt.ylabel('y [fm]')
      ax=plt.gca()
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)
      plt.colorbar(label=r'$\rho_B$'+" [1/fm^3]",cax=cax,format="%6.3f")

      plt.subplot(246)
      plt.imshow(rho_bab[it,:,:,ik].transpose(),extent=[xx[0], xx[-1], yy[0],yy[-1]], origin='lower',cmap='gnuplot2',vmin=0)
      plt.title("Baryon+antib. density")
      plt.xlabel('x [fm]')
      plt.ylabel('y [fm]')
      ax=plt.gca()
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)
      plt.colorbar(label=r'$\rho_{b+ab}$'+" [1/fm^3]",cax=cax,format="%6.3f")
   
      enratio=np.divide(en[it,:,:,ik],dens_had[it,:,:,ik],out=np.zeros_like(en[it,:,:,ik]), where=(dens_had[it,:,:,ik]!=0))

      plt.subplot(247)
      plt.imshow(enratio.transpose(),extent=[xx[0], xx[-1], yy[0],yy[-1]], origin='lower',cmap='gnuplot2',vmin=0,vmax=min(2,np.amax(enratio)))
      plt.title("Energy per particle")
      plt.xlabel('x [fm]')
      plt.ylabel('y [fm]')
      ax=plt.gca()
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)
      plt.colorbar(label="<E>/<N> [GeV]",cax=cax,format="%6.3f")

      st3ratio=np.divide(s_entr_dens[it,:,:,ik],(temp[it,:,:,ik]**3),out=np.zeros_like(s_entr_dens[it,:,:,ik]), where=(temp[it,:,:,ik]!=0))

      plt.subplot(248)
      plt.imshow(st3ratio.transpose()*hc3,extent=[xx[0], xx[-1], yy[0],yy[-1]], origin='lower',cmap='gnuplot2',vmax=min(25,np.amax(st3ratio*hc3)))
      plt.title("s/T^3")
      plt.xlabel('x [fm]')
      plt.ylabel('y [fm]')
      ax=plt.gca()
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)
      plt.colorbar(label="s/T^3",cax=cax,format="%6.3f")


      plt.tight_layout()

      plt.savefig(outdir+"/"+"medium_t_"+'{:05.2f}'.format(tt[it])+".png",dpi=150,pad_inches=0.)
      plt.close('all')

#     plots of s with and without the anisotropic correction

      plt.suptitle("t="+'{:4.2f}'.format(tt[it])+" - z="+'{:4.2f}'.format(zz[ik]))

      plt.subplot(121)
      plt.imshow(s_entr_densS[it,:,:,ik].transpose(),extent=[xx[0], xx[-1], yy[0],yy[-1]], origin='lower',cmap='gnuplot2')
      plt.title("Entropy density - no an. corr.")
      plt.xlabel('x [fm]')
      plt.ylabel('y [fm]')
      ax=plt.gca()
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)
      plt.colorbar(label=r"$s [fm^{-3}]$",cax=cax,format="%6.3f")

      plt.subplot(122)
      plt.imshow(s_entr_densA[it,:,:,ik].transpose(),extent=[xx[0], xx[-1], yy[0],yy[-1]], origin='lower',cmap='gnuplot2')
      plt.title("Entropy density - with an. corr.")
      plt.xlabel('x [fm]')
      plt.ylabel('y [fm]')
      ax=plt.gca()
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)
      plt.colorbar(label=r"$s [fm^{-3}]$",cax=cax,format="%6.3f")

      plt.tight_layout()

      plt.savefig(outdir+"/"+"entropy_density_t_"+'{:05.2f}'.format(tt[it])+".png",dpi=150,pad_inches=0.)
      plt.close('all')
