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

from ReduceAthenaData import *

from mpi4py import MPI


# setup latex fonts
#rcParams['text.usetex']=True
#rcParams['text.latex.unicode']=True


#matplotlib.rc('font', family='serif', serif='cm10')
#matplotlib.rc('text', usetex=True)
#matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']


vol_func = lambda rm,rp,thetam,thetap,phim,phip: \
           (1.0/3.0)*(rp**3-rm**3) * abs(np.cos(thetam)-np.cos(thetap)) * (phip-phim)

ni=13884
no=13929

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()


for i in range(ni,no+1,nprocs):
  fi=i+rank
  print fi, rank
  filename='./disk.out1.'+'{:05d}'.format(fi)+'.athdf'
  data=ReduceData(filename)
  quantities=data.keys()
  outputname='average_'+'{:05d}'.format(fi)+'.athdf'
  outf=h5py.File(outputname, 'w')

  for q in quantities:
    outf.create_dataset(q,data=data[q],dtype='>f8',shape=data[q].shape)

  outf.close()

