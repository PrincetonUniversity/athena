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


# setup latex fonts
#rcParams['text.usetex']=True
#rcParams['text.latex.unicode']=True


#matplotlib.rc('font', family='serif', serif='cm10')
#matplotlib.rc('text', usetex=True)
#matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']


# Athena++ modules
import athena_read

crat=2546.78
#time0=2.0*np.pi/gm1**0.5
time0=1/crat
rhomin=1.e-5
rhomax=10.0

Ermin=1.e-5
Ermax=1.0



def MakeRhoVSlice(data_cart, vx_cart, vy_cart, minval, maxval, vlim1, vlim2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,time,outputname,vel=0,figtype=0):
    if figtype == 0:
      plots, axes = plt.subplots(figsize=(6,10),dpi=300)
      plt.subplots_adjust(left=0.22,right=0.70,top=0.85,bottom=0.1)
    else:
      plots, axes = plt.subplots(figsize=(10,8.4),dpi=300)
      plt.subplots_adjust(left=0.16,right=0.78,top=0.85,bottom=0.1)

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
      cbaxes = plots.add_axes([0.8,0.15,0.03,0.6])

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



# Main function
def main(n,input_base,quantities=None,leveln=None, phiplot=0,width=23.44,ywidth=120,vel=0,rphiwidth=141.605):

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
    Er=data['Er']
    vr=data['vel1']
    vphi=data['vel3']
    vtheta=data['vel2']

#    Br=data['Bcc1']
#    Btheta=data['Bcc2']
#    Bphi=data['Bcc3']

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
    radius=np.zeros((nx3,nx1))
    phiangle=np.zeros((nx3,nx1))

    for j in range(0,nx3):
      for i in range(0,nx1):
        radius[j,i]=x1v[i] # convert rs to rg
        phiangle[j,i]=x3v[j]
    

    rhoslice=rho[:,thetapos,:]


    r1D=radius.ravel()
    t1D=phiangle.ravel()

    Erslice=Er[:,thetapos,:].copy()
    
    vrslice=vr[:,thetapos,:].copy()
    vtslice=vphi[:,thetapos,:].copy()
    
#    Brslice=Br[:,thetapos,:].copy()
#    Btslice=Btheta[:,thetapos,:].copy()

    rho1D=rhoslice.ravel()
    Er1D=Erslice.ravel()
    vr1D=vrslice.ravel()
    vt1D=vtslice.ravel()
#    Br1D=Brslice.ravel()
#    Bt1D=Btslice.ravel()

    vx=vr1D*np.cos(t1D) - vt1D*np.sin(t1D)
    vy=vr1D*np.sin(t1D) + vt1D*np.cos(t1D)

#    Bx=Br1D*np.cos(t1D) - Bt1D*np.sin(t1D)
#    By=Br1D*np.sin(t1D) + Bt1D*np.cos(t1D)

    xcoord=r1D*np.cos(t1D)
    ycoord=r1D*np.sin(t1D)


#create grid for plotting
    nx=nx1*2
    ny=nx
    xmin=-rphiwidth
# xmax=np.max(rcoord)
    xmax=rphiwidth
    ymin=-rphiwidth
    ymax=rphiwidth
    xgrid=np.linspace(xmin, xmax, nx)
    ygrid=np.linspace(ymin, ymax, ny)

    xmesh,ymesh=np.meshgrid(xgrid,ygrid)

    rmesh=(xmesh**2.0+ymesh**2.0)**0.5

    rindix=rmesh < 5.0
    rindix2=rmesh > rphiwidth    

    rho_cart=griddata(np.c_[xcoord,ycoord],rho1D,(xmesh,ymesh),method='nearest')
    Er_cart=griddata(np.c_[xcoord,ycoord],Er1D,(xmesh,ymesh),method='nearest')

    rho_cart[rindix] = rhomin
    Er_cart[rindix] = Ermin

    rho_cart[rindix2] = rhomin
    Er_cart[rindix2] = Ermin


    vx_cart=griddata(np.c_[xcoord,ycoord],vx,(xmesh,ymesh),method='nearest')
    vy_cart=griddata(np.c_[xcoord,ycoord],vy,(xmesh,ymesh),method='nearest')

#    Bx_cart=griddata(np.c_[xcoord,ycoord],Bx,(xmesh,ymesh),method='nearest')
#    By_cart=griddata(np.c_[xcoord,ycoord],By,(xmesh,ymesh),method='nearest')


    vx_cart[rindix] = 0.0
    vy_cart[rindix] = 0.0

#    Bx_cart[rindix] = 0.0
#    By_cart[rindix] = 0.0

    vx_cart[rindix2] = 0.0
    vy_cart[rindix2] = 0.0

#    Bx_cart[rindix2] = 0.0
#    By_cart[rindix2] = 0.0



    outputname='disk.'+'{:05d}'.format(n)+'_rho_z.png'

    labelname='$\\rho/\\rho_0$'

    timelabel='${\\rm time}='+"%4.2f"%(time/time0)+'{\ t_0}$'

    label1='$\\rho/\\rho_0$'
    label2='$v/c$'

    vlim1=10.0
    vlim2=1.e4


    MakeRhoVSlice(rho_cart, vx_cart, vy_cart, rhomin, rhomax, vlim1,vlim2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,timelabel,outputname,vel,figtype=1)

    outputname='disk.'+'{:05d}'.format(n)+'_Er_z.png'

    label1='$E_r/a_rT_0^4$'

    label2='$v$'

    vlim1=1.e-3
    vlim2=10.0


    MakeRhoVSlice(Er_cart, vx_cart, vy_cart, Ermin, Ermax, vlim1,vlim2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,timelabel,outputname,vel,figtype=1)




 # get cell centered coordinate values




# in order to calculate radiation viscosity
# we need radiation energy density, vphi, vr, rho, kappa

quantities=['rho','Er','vel1','vel2','vel3']

ni=200
no=200


rank = 0
nprocs = 1

for i in range(ni,no+1,nprocs):
  file_num=i+rank
#  print file_num, rank
  main(file_num,'disk.out1.',quantities,vel=1)


