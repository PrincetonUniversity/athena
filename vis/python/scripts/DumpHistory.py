#!/usr/bin/env python
# coding: utf-8

# In[1]:


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


# In[2]:


from ReduceSimpleData import *
from mpi4py import MPI

# In[3]:


files=sorted(glob.glob('disk.out1*athdf'))
num_file=len(files)


# In[4]:


nrloc=4
rloc=np.zeros(nrloc)
rloc[0]=120.0
rloc[1]=140.0
rloc[2]=160.0
rloc[3]=180.0
rindex=np.zeros(nrloc, dtype=np.uint32)


# In[7]:
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

count=0

for n in range(0,num_file,nprocs):
    fi=n+rank
    filename=files[fi]
    print(filename)
    data=ReduceData(filename)
    x1v=data['x1v']
    x2v=data['x2v']
    x2f=data['x2f']

    time=data['Time']
    rho=data['rho']  
    rhovr=data['rhovr']
    pgas=data['press']

    sigma_a=data['Sigma_a_0']
    sigma_p=data['Sigma_p_0']
    sigma_s=data['Sigma_s_0']
    
    Er=data['Er']
    Fr1=data['Fr1']
    Fr2=data['Fr2']
    Fr01=data['Fr01']
    Fr02=data['Fr02']

    Bcc1=data['Bcc1']
    Bcc2=data['Bcc2']
    Bcc3=data['Bcc3']
    PB1=data['PB1']
    PB2=data['PB2']
    PB3=data['PB3']
    PB=PB1+PB2+PB3
    lambda_a=data['lambda_a']
    dim=rho.shape
    ntheta=dim[0]
    nr=dim[1]
    radius=x1v
    
    if count == 0:
        ST_time=np.zeros(1)
        ST_sigma=np.zeros(nr)
        
        ST_rho_r1=np.zeros(ntheta)
        ST_rho_r2=np.zeros(ntheta)
        ST_rho_r3=np.zeros(ntheta)
        ST_rho_r4=np.zeros(ntheta)
        
        ST_Er_r1=np.zeros(ntheta)
        ST_Er_r2=np.zeros(ntheta)
        ST_Er_r3=np.zeros(ntheta)
        ST_Er_r4=np.zeros(ntheta)  
        
        ST_pg_r1=np.zeros(ntheta)
        ST_pg_r2=np.zeros(ntheta)
        ST_pg_r3=np.zeros(ntheta)
        ST_pg_r4=np.zeros(ntheta)  
        
        ST_kappa_r1=np.zeros(ntheta)
        ST_kappa_r2=np.zeros(ntheta)
        ST_kappa_r3=np.zeros(ntheta)
        ST_kappa_r4=np.zeros(ntheta)
        
        ST_B1_r1=np.zeros(ntheta)
        ST_B1_r2=np.zeros(ntheta)
        ST_B1_r3=np.zeros(ntheta)
        ST_B1_r4=np.zeros(ntheta)
        
        ST_B2_r1=np.zeros(ntheta)
        ST_B2_r2=np.zeros(ntheta)
        ST_B2_r3=np.zeros(ntheta)
        ST_B2_r4=np.zeros(ntheta)

        ST_B3_r1=np.zeros(ntheta)
        ST_B3_r2=np.zeros(ntheta)
        ST_B3_r3=np.zeros(ntheta)
        ST_B3_r4=np.zeros(ntheta)
        
        
        ST_PB_r1=np.zeros(ntheta)
        ST_PB_r2=np.zeros(ntheta)
        ST_PB_r3=np.zeros(ntheta)
        ST_PB_r4=np.zeros(ntheta)
        
        ST_lambda_r1=np.zeros(ntheta)
        ST_lambda_r2=np.zeros(ntheta)
        ST_lambda_r3=np.zeros(ntheta)
        ST_lambda_r4=np.zeros(ntheta)
        
        ST_rhovr=np.zeros((ntheta,nr))

 
        ST_Fr1_r1=np.zeros(ntheta)
        ST_Fr1_r2=np.zeros(ntheta)
        ST_Fr1_r3=np.zeros(ntheta)
        ST_Fr1_r4=np.zeros(ntheta)
        
        
        ST_Fr2_r1=np.zeros(ntheta)
        ST_Fr2_r2=np.zeros(ntheta)
        ST_Fr2_r3=np.zeros(ntheta)
        ST_Fr2_r4=np.zeros(ntheta)
        
        
        for i in range(nrloc):
            rindex[i]=np.abs(radius - rloc[i]).argmin()
            
    ST_time[0]=time
    for i in range(nr):
        ST_sigma[i] = np.sum(rho[:,i]*x1v[i]*np.abs(np.sin(x2f[1:])-np.sin(x2f[:-1])))
        
    ST_rho_r1[:]=rho[:,rindex[0]]
    ST_rho_r2[:]=rho[:,rindex[1]]
    ST_rho_r3[:]=rho[:,rindex[2]]
    ST_rho_r4[:]=rho[:,rindex[3]]

    ST_Er_r1[:]=Er[:,rindex[0]]
    ST_Er_r2[:]=Er[:,rindex[1]]
    ST_Er_r3[:]=Er[:,rindex[2]]
    ST_Er_r4[:]=Er[:,rindex[3]]
    
    kappa=(sigma_a+sigma_s)
    
    ST_kappa_r1[:]=kappa[:,rindex[0]]
    ST_kappa_r2[:]=kappa[:,rindex[1]]
    ST_kappa_r3[:]=kappa[:,rindex[2]]
    ST_kappa_r4[:]=kappa[:,rindex[3]]    

    ST_B1_r1[:]=Bcc1[:,rindex[0]]
    ST_B1_r2[:]=Bcc1[:,rindex[1]]
    ST_B1_r3[:]=Bcc1[:,rindex[2]]
    ST_B1_r4[:]=Bcc1[:,rindex[3]] 

    ST_B2_r1[:]=Bcc2[:,rindex[0]]
    ST_B2_r2[:]=Bcc2[:,rindex[1]]
    ST_B2_r3[:]=Bcc2[:,rindex[2]]
    ST_B2_r4[:]=Bcc2[:,rindex[3]] 
    
    ST_B3_r1[:]=Bcc3[:,rindex[0]]
    ST_B3_r2[:]=Bcc3[:,rindex[1]]
    ST_B3_r3[:]=Bcc3[:,rindex[2]]
    ST_B3_r4[:]=Bcc3[:,rindex[3]] 

    ST_pg_r1[:]=pgas[:,rindex[0]]
    ST_pg_r2[:]=pgas[:,rindex[1]]
    ST_pg_r3[:]=pgas[:,rindex[2]]
    ST_pg_r4[:]=pgas[:,rindex[3]] 

    ST_PB_r1[:]=PB[:,rindex[0]]
    ST_PB_r2[:]=PB[:,rindex[1]]
    ST_PB_r3[:]=PB[:,rindex[2]]
    ST_PB_r4[:]=PB[:,rindex[3]] 

    ST_rhovr[:,:]=rhovr[:,:]
 
    ST_lambda_r1[:]=lambda_a[:,rindex[0]]
    ST_lambda_r2[:]=lambda_a[:,rindex[1]]
    ST_lambda_r3[:]=lambda_a[:,rindex[2]]
    ST_lambda_r4[:]=lambda_a[:,rindex[3]] 
    
    ST_Fr1_r1[:]=Fr1[:,rindex[0]]
    ST_Fr1_r2[:]=Fr1[:,rindex[1]]
    ST_Fr1_r3[:]=Fr1[:,rindex[2]]
    ST_Fr1_r4[:]=Fr1[:,rindex[3]] 
    
    ST_Fr2_r1[:]=Fr2[:,rindex[0]]
    ST_Fr2_r2[:]=Fr2[:,rindex[1]]
    ST_Fr2_r3[:]=Fr2[:,rindex[2]]
    ST_Fr2_r4[:]=Fr2[:,rindex[3]] 
    
    count=count+1
        

    outfile='hist_'+files[fi][10:15]+'.npz'
# In[8]:

    np.savez(outfile,time=ST_time,radius=radius,theta=x2v,x1f=data['x1f'],x2f=x2f,location=rloc,surface_density=ST_sigma,
         rho1=ST_rho_r1,rho2=ST_rho_r2,rho3=ST_rho_r3,rho4=ST_rho_r4,
         Er1=ST_Er_r1,Er2=ST_Er_r2,Er3=ST_Er_r3,Er4=ST_Er_r4,
         kappa1=ST_kappa_r1,kappa2=ST_kappa_r2,kappa3=ST_kappa_r3,kappa4=ST_kappa_r4,
         B1_r1=ST_B1_r1,B1_r2=ST_B1_r2,B1_r3=ST_B1_r3,B1_r4=ST_B1_r4,
         B2_r1=ST_B2_r1,B2_r2=ST_B2_r2,B2_r3=ST_B2_r3,B2_r4=ST_B2_r4,
         B3_r1=ST_B3_r1,B3_r2=ST_B3_r2,B3_r3=ST_B3_r3,B3_r4=ST_B3_r4,
         PB_r1=ST_PB_r1,PB_r2=ST_PB_r2,PB_r3=ST_PB_r3,PB_r4=ST_PB_r4,
         pg_r1=ST_pg_r1,pg_r2=ST_pg_r2,pg_r3=ST_pg_r3,pg_r4=ST_pg_r4,
         rhovr=ST_rhovr)



# In[10]:





# In[ ]:




