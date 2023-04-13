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

from ReduceAthenaData3 import *

from mpi4py import MPI


# setup latex fonts
#rcParams['text.usetex']=True
#rcParams['text.latex.unicode']=True


#matplotlib.rc('font', family='serif', serif='cm10')
#matplotlib.rc('text', usetex=True)
#matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']


crat=8053.39
rmax=800

ni=479
no=759

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

for i in range(ni,no+1,nprocs):
  fi=i+rank
  print(fi, rank)
  filename='disk.out1.'+'{:05d}'.format(fi)+'.athdf'
  data=ReduceData(filename,rmax,crat)
  quantities = np.array([x.decode('ascii', 'replace') for x in data['VariableNames'][:]])
  coord_quantities = ('x1f', 'x2f', 'x1v', 'x2v', 'Time') 
  outputname='average_'+'{:05d}'.format(fi)+'.athdf'
  outf=h5py.File(outputname, 'w')

  for q in quantities:
    outf.create_dataset(q,data=data[q],dtype='>f8',shape=data[q].shape)

  for q in coord_quantities:
    outf.create_dataset(q,data=data[q],dtype='>f8',shape=data[q].shape)

  outf.close()

