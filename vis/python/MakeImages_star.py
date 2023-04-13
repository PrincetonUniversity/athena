import matplotlib
matplotlib.use('Agg')


import sys
sys.settrace
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pylab import *
import struct
import array
import os
from scipy.interpolate import griddata


import h5py

from mpi4py import MPI

# setup latex fonts
rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True


matplotlib.rc('font', family='serif', serif='cm10')
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']


# Athena++ modules
import athena_read

grav0=16.4983
r0=12.9837;
time0=3.93232 # in hours

rhomin=1.e-3
rhomax=0.2

Ermin=1.e-3
Ermax=1.e2



def MakeRhoVSlice(data_cart, vx_cart, vy_cart, minval, maxval, vlim1, vlim2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,time,outputname,vel=0,figtype=0):
    if figtype == 0:
      plots, axes = plt.subplots(figsize=(12,6),dpi=300)
      plt.subplots_adjust(left=0.22,right=0.70,top=0.85,bottom=0.1)
    else:
      plots, axes = plt.subplots(figsize=(16,7.4),dpi=300)
      plt.subplots_adjust(left=0.15,right=0.8,top=0.85,bottom=0.12)

    if figtype == 0:
      plt.xlabel('$ r/r_g$', size = 30)
      plt.ylabel('$ z/r_g$', size = 30)
    else:
      plt.xlabel('$ x/r_0$', size = 30)
      plt.ylabel('$ y/r_0$', size = 30)      


    plt.title(time,size=25,y=1.12)

    
    if vel>0:
      speed=np.sqrt(vx_cart**2+vy_cart**2)
      if vlim2 < vlim1:
         vlim2=speed.max()
      speed=np.clip(speed,vlim1,vlim2)
      logspeed=np.log10(speed)
      vcolor=axes.imshow(speed,cmap='jet',norm=LogNorm(vmin=vlim1,vmax=vlim2),origin='lower', extent=[xmin,xmax,ymin,ymax])
      if figtype == 0:
        velaxes = plots.add_axes([0.22, 0.89, 0.5, 0.03])
      else:
        velaxes = plots.add_axes([0.24, 0.89, 0.5, 0.03])

      velcbar = plots.colorbar(vcolor, cax = velaxes, orientation = 'horizontal')
      velcbar.ax.tick_params(labelsize=20)
      if speed.max() > vlim1:
        velfield=axes.streamplot(xgrid,ygrid,vx_cart,vy_cart,cmap='jet',density=1,color=logspeed,arrowsize=3.0)
      else:
        velfield=axes.streamplot(xgrid,ygrid,vx_cart,vy_cart,cmap='jet',density=1,arrowsize=3.0)
        
      axes.set_ylim([ymin,ymax])

       
    im = axes.imshow(data_cart,cmap='RdGy_r', norm=LogNorm(vmin = minval, vmax=maxval), \
                     origin='lower', extent=[xmin,xmax,ymin,ymax])
    if figtype == 0:      
      cbaxes = plots.add_axes([0.72,0.15,0.03,0.6])
    else:
      cbaxes = plots.add_axes([0.82,0.15,0.03,0.6])

    cbar=plots.colorbar(im,cax=cbaxes)
    cbar.set_label(label1, size=30)
    cbar.ax.tick_params(labelsize=25)
    if figtype == 0:
      axes.set_xticks([0,10,20,30,40,50])
    axes.yaxis.set_tick_params(labelsize=20)
    axes.xaxis.set_tick_params(labelsize=20)
    
    axes.set_aspect('auto')
    plt.savefig(outputname)
    plt.close(plots)


def MakePolarPlot(data,xgrid,ygrid,minval,maxval,title,clabel,outputname,ctick,xtick=None,xticklabel=None,ytick=None,yticklabel=None,logscale=0,cmap='RdGy_r'):

    plots, axes = plt.subplots(figsize=(8,12),dpi=300)
    plt.subplots_adjust(left=0.15,right=0.8,top=0.9,bottom=0.1)

    plt.xlabel('$\\phi$', size = 30)
    plt.ylabel('$ r/r_{\\odot}$', size = 30)


    plt.title(title,size=25,y=1.02)
    if logscale ==1:
      im=axes.pcolormesh(xgrid,ygrid,data,norm=LogNorm(vmin=minval,vmax=maxval),cmap=cmap)
    else:
      norm=mpl.colors.Normalize(vmin=minval, vmax=maxval)
      im=axes.pcolormesh(xgrid,ygrid,data,norm=norm,cmap=cmap)   	

    cbaxes = plots.add_axes([0.81,0.2,0.03,0.6])

    cbar=plots.colorbar(im,cax=cbaxes,ticks=ctick)
    cbar.set_label(clabel, size=30)
    cbar.ax.tick_params(labelsize=25)
    if xtick is not None:
      axes.set_xticks(xtick)
      axes.set_xticklabels(xticklabel)
    if ytick is not None:
      axes.set_yticks(ytick)
      axes.set_yticklabels(yticklabel)    	
    axes.yaxis.set_tick_params(labelsize=20)
    axes.xaxis.set_tick_params(labelsize=20)
    
    axes.set_aspect('auto')
    plt.savefig(outputname)
    plt.close(plots)


# Main function
def main(n,input_base,rmax,quantities=None,leveln=None):

  # Go through list of files
    # Determine filenames
    input_filename = input_base+'{:05d}'.format(n)+'.athdf'
    print input_filename + '\n'
    # Read attributes and data
    with h5py.File(input_filename, 'r') as f:
      attributes = f.attrs.items()
      attrs = dict(attributes)
      level = f.attrs['MaxLevel']
      time = f.attrs['Time']
    subsample = False
    if leveln is not None:
      if level > leveln:
        subsample = True
      level = leveln
    data = athena_read.athdf(input_filename, quantities=quantities,
        level=level, subsample=subsample)

    # Determine new grid size
    nx1 = attrs['RootGridSize'][0] * 2**level
    nx2 = attrs['RootGridSize'][1] * 2**level
    nx3 = attrs['RootGridSize'][2] * 2**level


    rho=data['rho']
    vr=data['vel1']
    vphi=data['vel3']
    vtheta=data['vel2']
    Er=data['Er']

    #plot rho, v and Er, B
    # make plots in r-theta plane and mid-plane

    x1f=data['x1f']
    x2f=data['x2f']
    x3f=data['x3f']


    x1v=np.zeros(nx1)
    x2v=np.zeros(nx2)
    x3v=np.zeros(nx3)

    for i in range(nx1):
      x1v[i]=0.75*(x1f[i+1]**4.0 - x1f[i]**4.0)/(x1f[i+1]**3.0-x1f[i]**3.0);
  
    for j in range(nx2):
      x2v[j]=((np.sin(x2f[j+1]) - x2f[j+1] * np.cos(x2f[j+1])) \
            -(np.sin(x2f[j]) - x2f[j] * np.cos(x2f[j]))) \
              / (np.cos(x2f[j]) - np.cos(x2f[j+1]))

    for k in range(nx3):
      x3v[k] = 0.5 * (x3f[k+1]+x3f[k])
  

   
###########################################################################
    # now the mid-plane plots

    # find the position corresponding to phiplot position
    thetapos=np.abs(x2v-0.5*np.pi).argmin()
    rmaxpos=np.abs(x1v-rmax).argmin()
    

    rhoslice=rho[:,thetapos,:rmaxpos]
    Erslice=Er[:,thetapos,:rmaxpos]
    Erslice=Erslice**0.25
    x1f=x1f[:rmaxpos+1]

    logr=np.log10(x1f)
    xmesh,ymesh=np.meshgrid(x3f,logr)




    outputname='star.'+'{:05d}'.format(n)+'_rho.png'

    labelname='$\\rho/\\rho_0$'

    timelabel='${\\rm t}='+"%4.1f"%(time*time0)+'{\ {\\rm hours}}$'

    xticks=[0,4*np.pi/180,8*np.pi/180]
    xticklabels=['$0$','$4^{\circ}$','$8^{\circ}$']
    yticks=[np.log10(11),np.log10(11.5),np.log10(12),np.log10(12.5),np.log10(13),np.log10(13.5)]
    yticklabels=['$11$','$11.5$','$12$','$12.5$','$13$','$13.5$']
    minval=1.e-7
    maxval=1.e2
    cticks=[minval,100*minval,1.e4*minval,1.e6*minval,1.e8*minval]
    MakePolarPlot(np.transpose(rhoslice),xmesh,ymesh,minval,maxval,timelabel,labelname,outputname,cticks,xtick=xticks,xticklabel=xticklabels,ytick=yticks,yticklabel=yticklabels,logscale=1)
    outputname2='star.'+'{:05d}'.format(n)+'_Er.png'
    labelname='$T_r/T_0$'
    cticks=[0,1,2,3,4]
    MakePolarPlot(np.transpose(Erslice),xmesh,ymesh,0.1,4,timelabel,labelname,outputname2,cticks,xtick=xticks,xticklabel=xticklabels,ytick=yticks,yticklabel=yticklabels,logscale=0,cmap='inferno')


#    outputname='disk.'+'{:05d}'.format(n)+'_Er_z.png'

#    label1='$E_r/a_rT_0^4$'

#    label2='$B$'

#    vlim1=1.e-3
#    vlim2=10.0


#    MakeRhoVSlice(Er_cart, Bx_cart, By_cart, Ermin, Ermax, vlim1,vlim2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,timelabel,outputname,vel,figtype=1)




 # get cell centered coordinate values




# in order to calculate radiation viscosity
# we need radiation energy density, vphi, vr, rho, kappa

quantities=['rho','Er','vel1','vel2','vel3']

ni=4000
no=4000

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

for i in range(ni,no+1,nprocs*5):
  file_num=i+rank*5
#  print file_num, rank
  main(file_num,'./Data/star.out1.',14,quantities)



