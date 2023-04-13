#! /usr/bin/env python
import os

f = open('outlist', 'r')

for i in range(43):
  oriname='CVIsoB2_rho'+'{:04d}'.format(i)+'.png'
  outname=f.readline()
  command='mv '+oriname+' '+outname
  print command
  os.system(command)

f.close()




