#!/usr/bin/env python
# coding: utf-8

# In[7]:


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
from astropy.convolution import convolve, Box1DKernel


import h5py
from mpi4py import MPI


plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})
# for Palatino and other serif fonts use:
#plt.rcParams.update({
#    "text.usetex": True,
#    "font.family": "serif",
#    "font.serif": ["Palatino"],
#})
from importlib import reload


# In[8]:


import athena_read
athena_read=reload(athena_read)


# In[9]:


def MakeRhoVSlice(data_cart, vx_cart, vy_cart, minval, maxval, vlim1, vlim2, width, xmin, xmax, ymin, ymax,
                  xgrid,ygrid,label1,label2,time,outputname,vel=0,figtype=0,colormap='RdGy_r'):
    if figtype == 0:
        plots, axes = plt.subplots(figsize=(11,10),dpi=300)
        plt.subplots_adjust(left=0.11,right=0.85,top=0.92,bottom=0.1)
    else:
        plots, axes = plt.subplots(figsize=(10,8.4),dpi=300)
        plt.subplots_adjust(left=0.16,right=0.78,top=0.85,bottom=0.1)

    if figtype == 0:
        plt.xlabel('$ x/r_0$', size = 30)
        plt.ylabel('$ z/r_0$', size = 30)
    else:
        plt.xlabel('$ x/r_0$', size = 30)
        plt.ylabel('$ y/r_0$', size = 30)      


    plt.title(time,size=25,y=1.02)

    
    if vel>0:
        speed=np.sqrt(vx_cart**2+vy_cart**2)
        if vlim2 < vlim1:
             vlim2=speed.max()
        speed=np.clip(speed,vlim1,vlim2)
        logspeed=np.log10(speed)
        vcolor=axes.imshow(speed,cmap='jet',norm=LogNorm(vmin=vlim1,vmax=vlim2),origin='lower', 
                         extent=[xmin,xmax,ymin,ymax])
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


    im = axes.imshow(data_cart,cmap=colormap, norm=LogNorm(vmin = minval, vmax=maxval),                      origin='lower', extent=[xmin,xmax,ymin,ymax])
    if figtype == 0:      
        cbaxes = plots.add_axes([0.86,0.18,0.03,0.6])
    else:
        cbaxes = plots.add_axes([0.8,0.15,0.03,0.6])

    cbar=plots.colorbar(im,cax=cbaxes)
    cbar.set_label(label1, size=30)
    cbar.ax.tick_params(labelsize=25)
#    if figtype == 0:
#      axes.set_xticks([0,10,20,30,40,50])
    axes.yaxis.set_tick_params(labelsize=20)
    axes.xaxis.set_tick_params(labelsize=20)
    
    axes.set_aspect('auto')
    plt.savefig(outputname)
    plt.close(plots)


# In[10]:


def main(n,input_base,quantities=None,leveln=None, vel=0,colormap='RdGy_r'):
    
    input_filename = input_base+'{:05d}'.format(n)+'.athdf'
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
    
    nx1 = attrs['RootGridSize'][0] * 2**level
    nx2 = attrs['RootGridSize'][1] * 2**level
    nx3 = attrs['RootGridSize'][2] * 2**level

    rho=data['rho']
    vr=data['vel1']
    vphi=data['vel3']
    vtheta=data['vel2']
    x1f=data['x1f']
    x2f=data['x2f']
    x3f=data['x3f']
 
    x1v=np.zeros(nx1)
    x2v=np.zeros(nx2)
    x3v=np.zeros(nx3)

    for i in range(nx1):
        x1v[i]=0.75*(x1f[i+1]**4.0 - x1f[i]**4.0)/(x1f[i+1]**3.0-x1f[i]**3.0)

    
    for j in range(nx2):
        x2v[j]=((np.sin(x2f[j+1]) - x2f[j+1] * np.cos(x2f[j+1]))             -(np.sin(x2f[j]) - x2f[j] * np.cos(x2f[j])))               / (np.cos(x2f[j]) - np.cos(x2f[j+1]))

    for k in range(nx3):
        x3v[k] = 0.5 * (x3f[k+1]+x3f[k])

    phipos1=np.abs(x3v-x3v[0]).argmin()
    phipos2=np.abs(x3v-x3v[0]-np.pi).argmin()
    phi_1=x3v[phipos1]
    phi_2=x3v[phipos2]
    radius=np.zeros((nx2,2*nx1))
    tangle=np.zeros((nx2,2*nx1))

    for j in range(0,nx2):
        for i in range(nx1,2*nx1):
            radius[j,i]=x1v[i-nx1]
            tangle[j,i]=x2v[j]
        for i in range(nx1-1,-1,-1):
            radius[j,i]=-x1v[nx1-1-i]
            tangle[j,i]=x2v[j]

    r1D=radius.ravel()
    t1D=tangle.ravel()
    xcoord=r1D*np.sin(t1D)
    ycoord=r1D*np.cos(t1D)
    rhoslice1=rho[phipos1,:,:]
    vrslice1=vr[phipos1,:,:]
    vtslice1=vtheta[phipos1,:,:]
    rhoslice2=rho[phipos2,:,:]
    vrslice2=vr[phipos2,:,:]
    vtslice2=vtheta[phipos2,:,:]
    rhoslice=np.zeros((nx2,2*nx1))
    vrslice=np.zeros((nx2,2*nx1))
    vtslice=np.zeros((nx2,2*nx1))
    for j in range(0,nx2):
        rhoslice[nx2-1-j,:nx1]=rhoslice2[j,::-1]
        rhoslice[j,nx1:]=rhoslice1[j,:]
        vrslice[nx2-1-j,:nx1]=vrslice2[j,::-1]
        vrslice[j,nx1:]=vrslice1[j,:]
        vtslice[nx2-1-j,:nx1]=vtslice2[j,::-1]
        vtslice[j,nx1:]=vtslice1[j,:]

        
    rho1D=rhoslice.ravel()
    vr1D=vrslice.ravel()
    vt1D=vtslice.ravel()

    vx=vr1D*np.cos(t1D) - vt1D*np.sin(t1D)
    vy=vr1D*np.sin(t1D) + vt1D*np.cos(t1D)

    nx=2*nx1
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
    rindix=rmesh < 1.0
    rindix2=rmesh > rphiwidth 
    
    rho_cart=griddata(np.c_[xcoord,ycoord],rho1D,(xmesh,ymesh),method='nearest')
    rho_cart[rindix]=rhomin
    rho_cart[rindix2]=rhomin


    
    vx_cart=griddata(np.c_[xcoord,ycoord],vx,(xmesh,ymesh),method='nearest')
    vy_cart=griddata(np.c_[xcoord,ycoord],vy,(xmesh,ymesh),method='nearest')
    vx_cart[rindix] = 0.0
    vy_cart[rindix] = 0.0
    vx_cart[rindix2] = 0.0
    vy_cart[rindix2] = 0.0
    
    
    outputname='disk.'+'{:05d}'.format(n)+'_rho.png'

    labelname='$\\rho/\\rho_0$'

    timelabel='$t='+"%4.2f"%(time/time0)+'{\ t_0}$'

    label1='$\\rho/\\rho_0$'
    label2='$v/c$'

    MakeRhoVSlice(rho_cart, vx_cart, vy_cart, rhomin, rhomax, vlim1,vlim2, width, 
              xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,timelabel,outputname,vel,figtype=0,colormap=colormap)
    


# In[ ]:


ni=430
no=430
colormap='inferno'
input_base='./disk.out1.'
quantities=['rho','vel1','vel2','vel3']
leveln=None
rphiwidth=23.44
width=23.44
rhomin=1.e-4
rhomax=1.
Ermin=1.
Ermax=1.e6
omega0=1.79925
time0=2*np.pi/omega0
vlim1=1.0
vlim2=1.e4
vel=0
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

for i in range(ni,no+1,nprocs):
    file_num=i+rank
    print(file_num, rank)
    main(file_num,input_base,quantities,vel=0,colormap=colormap)


# In[5]:




