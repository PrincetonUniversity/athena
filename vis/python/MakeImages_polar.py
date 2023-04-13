#!/usr/bin/env python
# coding: utf-8

# In[29]:


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
import glob
import h5py
from scipy.interpolate import griddata
import scipy.integrate as integrate


plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Helvetica",
    "mathtext.default": "regular",
    "text.latex.preamble": [r"""\usepackage{bm}"""],
})


# In[2]:

from mpi4py import MPI
import athena_read

def makeplot(data,rmax,vmin,vmax,clabel,timelabel,outputname,cmap='inferno'):
    plots, axes = plt.subplots(figsize=(4.5,3.5),dpi=300)
    im=axes.pcolormesh(x_grid,y_grid,data,norm=LogNorm(vmin=vmin, vmax=vmax),cmap=cmap)
    axes.set_aspect('equal')
    axes.set_xlabel(r'$\mathbf{r \cos\phi}$')
    axes.set_ylabel(r'$\mathbf{r \sin\phi}$')
    axes.set_ylim([-rmax, rmax])
    axes.set_xlim([-rmax,rmax])
    cbar=plots.colorbar(im, ax=axes)
    cbar.set_label(clabel,style='italic')
    plt.title(timelabel)
    plt.savefig(outputname)



# In[3]:


grav0=1.09465
r0=28.1427
time0=4.54868/24 # in days

rhomin=1.e-6
rhomax=10.0

Trmin=0.01
Trmax=10

rmax=800


# In[4]:


files=sorted(glob.glob('star*athdf'))
quantities=['rho','vel1','vel2','vel3','Er']

nfile=len(files)

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

for i in range(0,nfile,nprocs):
# In[5]:
#for input_filename in files:
    input_filename=files[i+rank]
    print(input_filename,rank)
    with h5py.File(input_filename, 'r') as f:
        attributes = f.attrs.items()
        attrs = dict(attributes)
        level = f.attrs['MaxLevel']
        datatime = f.attrs['Time']
        subsample = False

        data = athena_read.athdf(input_filename, quantities=quantities,
            level=level, subsample=subsample)


# In[6]:


    x1f=data['x1f']
    x2f=data['x2f']
    x3f=data['x3f']
    x1v=data['x1v']
    x2v=data['x2v']
    x3v=data['x3v']
    rho=data['rho']
    vel1=data['vel1']
    vel2=data['vel2']
    vel3=data['vel3']
    Er=data['Er']


# In[7]:


    nr=int(len(x1v))
    ntheta=int(len(x2v))
    nphi=int(len(x3v))
# extend radius, theta, phi to 3d grid for all cells
#cell size is d(r^3/3)(dcostheta)dphi
    thetapos=np.abs(x2v-0.5*np.pi).argmin()
    rmaxpos=np.abs(x1v-rmax).argmin()
    grid_phi,grid_r=np.meshgrid(x3v,x1v[:rmaxpos],indexing='ij')


# In[8]:


#get radial surface area
    rhoslice=rho[:,thetapos,:rmaxpos]
    Erslice=Er[:,thetapos,:rmaxpos]


# In[9]:


    grid_r, grid_phi = np.meshgrid(x1f[:rmaxpos+1], x3f) 
    x_grid = grid_r*np.cos(grid_phi)
    y_grid = grid_r*np.sin(grid_phi)


# In[21]:


    output_name1='YSGL_rho_'+input_filename[10:15]+'.png'
    timelabel='$t='+f"{datatime*time0:4.1f}"+r'{\ {\rm days}}$'
    output_name2='YSGL_Tr_'+input_filename[10:15]+'.png'


# In[34]:
    rholabel=r'$\mathbf{\rho/\rho_0}$'
    makeplot(rhoslice,rmax,rhomin,rhomax,rholabel,timelabel,output_name1,cmap='viridis')

    trlabel= r'$\mathbf{T_r/T_0}$'
    makeplot(Erslice**0.25,rmax,Trmin,Trmax,trlabel,timelabel,output_name2,cmap='inferno')




# In[ ]:




