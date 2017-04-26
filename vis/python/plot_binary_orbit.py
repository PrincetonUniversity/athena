import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import sys

#direname = '/tigress/jiming/athena/bin/binary/'
#targname = 'orbit' #m4-tapper'
#histname = 'e0.5_q0.1_i0.tab'
#histname = 'e0.3_q0.8_i0.2_2.tab'
#fname=direname+targname+'/'+histname

def plot_binary_orbit(fname):
  """ semi-major axis unity """
  torb = 2.*np.pi
  dtype1 = np.dtype([('time', 'd'), ('xs1', 'd'),('xs2','d'),('xs3','d'),\
                     ('xp1','d'),('xp2','d'),('xp3','d')])
  ahist = np.loadtxt(fname, dtype=dtype1, skiprows=0, usecols=(0,4,5,6,7,8,9))
  
  mpl.rcParams['legend.fontsize'] = 10
  fig = plt.figure()
  ax = fig.gca(projection='3d')
  # plot trajectories
  ax.plot(ahist['xs1'],ahist['xs2'],ahist['xs3'],'b-*',markevery=[-1],zdir='z', label='secondary')
  ax.plot(ahist['xp1'],ahist['xp2'],ahist['xp3'],'g-*',markevery=[-1],zdir='z', label='primary')
  ax.legend(fontsize=15)
  plt.xlabel('x')
  plt.ylabel('y')
  
  # # Create cubic bounding box to simulate equal aspect ratio
  xmax,ymax,zmax=max(max(ahist['xs1']),max(ahist['xp1'])),\
                 max(max(ahist['xs2']),max(ahist['xp2'])),\
                 max(max(ahist['xs3']),max(ahist['xp3']))
  xmin,ymin,zmin=min(min(ahist['xs1']),min(ahist['xp1'])),\
                 min(min(ahist['xs2']),min(ahist['xp2'])),\
                 min(min(ahist['xs3']),min(ahist['xp3']))
  max_range = np.array([xmax-xmin, ymax-ymin, zmax-zmin]).max()
  Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(xmax+xmin)
  Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(ymax+ymin)
  Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(zmax+zmin)
  # Comment or uncomment following both lines to test the fake bounding box:
  for xb, yb, zb in zip(Xb, Yb, Zb):
     ax.plot([xb], [yb], [zb], 'w')
  
  
  # mark the final location and center-of-mass
  ax.plot([0,0],[0,0],[0,0],'r+',markersize=10)
  ax.plot(1.2*np.linspace(xmin,xmax,100),np.zeros(100),np.zeros(100),'k--')
  ax.plot(np.zeros(100),1.2*np.linspace(ymin,ymax,100),np.zeros(100),'k--')
  ax.plot(np.zeros(100),np.zeros(100),6*np.linspace(zmin,zmax,100),'k--')
  plt.grid()
  plt.show()


if __name__=='__main__':
    if len( sys.argv ) < 2:
        print "Please specify orbit profile filename" 
        exit(  )
    fname = sys.argv[1]
    plot_binary_orbit(fname)

# end of the script
