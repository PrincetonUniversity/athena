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


from ReduceHisData import *


# In[3]:


crat=8053.39


# In[4]:


files=sorted(glob.glob('disk.out1*athdf'))
ni=0
no=200
files=files[ni:no]

num_file=len(files)


# In[5]:



#select file range to average


# In[6]:


count=0

for filename in files:
    print(filename)
    data0=ReduceData(filename,crat)
    if count == 0:
        data=data0
        data['Time']=np.array(data0['Time'])
    else:
        data['Time']=np.append(data['Time'],data0['Time'])
        for q in data['quantities']:
            data[q] = data[q] + data0[q]
    
    count=count+1
        


# In[33]:


outfile='ave_'+files[0][10:15]+'_'+files[num_file-1][10:15]+'.npz'


# In[34]:


np.savez(outfile,time=data['Time'],x1f=data['x1f'],x2f=data['x2f'],x1v=data['x1v'],x2v=data['x2v'],rho=data['rho'],
         press=data['press'],vel1=data['vel1'],vel2=data['vel2'],vel3=data['vel3'],Er=data['Er'],
         Fr1=data['Fr1'],Fr2=data['Fr2'],Fr3=data['Fr3'],Pr11=data['Pr11'],Pr22=data['Pr22'],Pr33=data['Pr33'],
         Pr12=data['Pr12'],Pr13=data['Pr13'],Pr23=data['Pr23'],Pr21=data['Pr21'],Pr31=data['Pr31'],
        Pr32=data['Pr32'],Er0=data['Er0'],Fr01=data['Fr01'],Fr02=data['Fr02'],Fr03=data['Fr03'],
         Sigma_s=data['Sigma_s_0'],sigma_a=data['Sigma_a_0'],Sigma_p=data['Sigma_p_0'],
         Bcc1=data['Bcc1'],Bcc2=data['Bcc2'],Bcc3=data['Bcc3'],Maxwell=data['Maxwell'],Reynolds=data['Reynolds'],
         BrBphi=data['BrBphi'],BthetaBphi=data['BthetaBphi'],
        rhovrvphi=data['rhovrvphi'],rhovthetavphi=data['rhovthetavphi'],vrEr=data['vrEr'],vthetaEr=data['vthetaEr'],
         rhoPB=data['rhoPB'],rhosq=data['rhosq'],PB1=data['PB1'],
        PB2=data['PB2'],PB3=data['PB3'],PB=data['PB'],PBsq=data['PBsq'],kappa_s=data['kappa_s'],kappa_a=data['kappa_a'],
         Radacc=data['Radacc'],RhoVr=data['RhoVr'],RhoVtheta=data['RhoVtheta'],
        RhoVphi=data['RhoVphi'],RhoVout=data['RhoVout'],RhoVin=data['RhoVin'],Ekin1=data['Ekin1'],Ekin2=data['Ekin2'],
         Ekin3=data['Ekin3'],tgas=data['tgas'],Fr01Sigma=data['Fr01Sigma'],
        Fr02Sigma=data['Fr02Sigma'],Fr01kappa=data['Fr01kappa'],Fr02kappa=data['Fr02kappa'],
         rhovrsq=data['rhovrsq'],vrsq=data['vrsq'],rhovphiout=data['rhovphiout'],
        rhovphiin=data['rhovphiin'],rhoinflow=data['rhoinflow'],rhooutflow=data['rhooutflow'],
         rhovrvphisq=data['rhovrvphisq'],Ersq=data['Ersq'],
         rhovrvphiEr=data['rhovrvphiEr'],rhovrvphiPB=data['rhovrvphiPB'],radvis=data['radvis'],
         velshear=data['velshear'],divvel=data['divvel'],
        raddivwork=data['raddivwork'],radshearwork=data['radshearwork'],curvtrho=data['curvtrho'],
         curvtrhosq=data['curvtrhosq'],curvgrads=data['curvgrads'], curvgradssq=data['curvgradssq'],
         Fr1new=data['Fr1new'],Fr2new=data['Fr2new'],Fr3new=data['Fr3new'],lambda_a=data['lambda_a'])


# In[ ]:




