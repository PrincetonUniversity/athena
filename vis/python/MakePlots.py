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
rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True


matplotlib.rc('font', family='serif', serif='cm10')
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']

vmax=100
sigma=1.e3

vol_func = lambda rm,rp,thetam,thetap: \
            1.0/3.0*(rp**3-rm**3) * abs(np.cos(thetam)-np.cos(thetap)) * 2.0*np.pi

rarea_func = lambda r,thetam,thetap: \
           r**2.0 * abs(np.cos(thetam)-np.cos(thetap)) * 2.0*np.pi

#the tabulated color
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]   

for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.) 

def MakeRhoVSlice(data_cart, vx_cart, vy_cart, minval, maxval, vlim1, vlim2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,outputname,vel=0,logscale=1,xticks=None,taux=None,tauy=None):
    plots, axes = plt.subplots(figsize=(5.5,10),dpi=300)
    plt.xlabel('$r\\sin\\theta/r_g$', size = 30)
    plt.ylabel('$r\\cos\\theta/r_g$', size = 30)
    plt.subplots_adjust(left=0.23,right=0.7,top=0.85,bottom=0.1)
#    plt.title(time,size=25,y=1.12)
    
    if vel>0:
      speed=np.sqrt(vx_cart**2+vy_cart**2)
      if vlim2 < vlim1:
         vlim2=speed.max()
      speed=np.clip(speed,vlim1,vlim2)
      logspeed=np.log10(speed)
      vcolor=axes.imshow(speed,cmap='jet',norm=LogNorm(vmin=vlim1,vmax=vlim2),origin='lower', extent=[xmin,xmax,ymin,ymax])
      velaxes = plots.add_axes([0.24, 0.9, 0.5, 0.03])
      velcbar = plots.colorbar(vcolor, cax = velaxes, orientation = 'horizontal')
      velcbar.set_label(label2, size=30, labelpad=-90)
      velcbar.ax.tick_params(labelsize=25)

      velfield=axes.streamplot(xgrid,ygrid,vx_cart,vy_cart,cmap='jet',arrowsize=3,density=0.6,linewidth=1.5,color=logspeed)
        
      axes.set_ylim([ymin,ymax])

    if taux is not None:
      axes.plot(taux,tauy,color='black',linestyle='dashed',linewidth=3.0)

    if logscale > 0:
      im = axes.imshow(data_cart,cmap='RdGy_r', norm=LogNorm(vmin = minval, vmax=maxval), \
                       origin='lower', extent=[xmin,xmax,ymin,ymax])
    else:
      im = axes.imshow(data_cart,cmap='RdGy_r', vmin=minval, vmax=maxval, \
                       origin='lower', extent=[xmin,xmax,ymin,ymax])
    
    
    
    cbaxes = plots.add_axes([0.72,0.15,0.03,0.6])
    cbar=plots.colorbar(im,cax=cbaxes)
    cbar.set_label(label1, size=30)
    cbar.ax.tick_params(labelsize=25)
    if xticks is not None:
      axes.set_xticks(xticks)
    axes.yaxis.set_tick_params(labelsize=20)
    axes.xaxis.set_tick_params(labelsize=20)
  

    axes.set_aspect('auto')

    plt.savefig(outputname)
    plt.close(plots)



def PlotHistory(data, col1, label1, ymin, ymax, ylabel, filename, col2=0, label2='', col3=0, label3='',xsticks=None):
    plots, axes = plt.subplots(figsize=(8,9),dpi=300)
    plt.xlabel('$ t/(r_g/c)$', size = 30)
    plt.ylabel(ylabel, size = 30)
    plt.subplots_adjust(left=0.2,right=0.88,top=0.9,bottom=0.1)
    plt.ylim([ymin,ymax])
#    plt.title(time,size=25,y=1.12)
    plt.plot(data[:,0],data[:,col1],color='black',label=label1,linewidth=2.0)
    if col2 >0:
      plt.plot(data[:,0],data[:,col2],color='r',label=label2,linewidth=2.0)
    if col3 > 0:
      plt.plot(data[:,0],data[:,col3],color='g',label=label3,linewidth=3.0)
    plt.legend(loc='upper left',frameon=False)
    axes.set_aspect('auto')
    axes.yaxis.set_tick_params(labelsize=20)
    if xsticks is not None:
      axes.set_xticks(xsticks)
    axes.xaxis.set_tick_params(labelsize=20)
    plt.savefig(filename)
    plt.close(plots)

def PlotProfile(datax, datay, xmin, xmax, ymin, ymax,  ylabel, label1, filename, xlabel='$r/r_g$', logscale=0, datay1_2=None, datay1_3=None, datax2=None, datay2=None, datay2_2=None, datay2_3=None, datax3=None, datay3=None, datay3_2=None, datay3_3=None, datax4=None, datay4=None, datax5=None, datay5=None, label2='', label3='', label4='', label5='',datay2line=None):
    plots, axes = plt.subplots(figsize=(9,11),dpi=300)
    plt.xlabel(xlabel, size = 30)
    plt.ylabel(ylabel, size = 30)
    plt.subplots_adjust(left=0.2,right=0.88,top=0.9,bottom=0.1)
    plt.ylim([ymin,ymax])
    plt.xlim([xmin,xmax])
    if logscale > 0:
      axes.set_yscale('log')
#    axes.set_xscale('log')
#    plt.title(time,size=25,y=1.12)
    if label1 != None:
      plt.plot(datax,datay,color='black',label=label1,linewidth=2.0)
    else:
      plt.plot(datax,datay,color='black',linewidth=2.0)      
    if datay1_2 is not None:
      plt.plot(datax,datay1_2,color=tableau20[19],linestyle='dashed',linewidth=2.0)
    if datay1_3 is not None:
      plt.plot(datax,datay1_3,color='black',linestyle='dashed',linewidth=4.0)
    if datay2 is not None:
      if datay2line is not None:
        plt.plot(datax2,datay2,color='red',label=label2,linewidth=2.0,linestyle=datay2line)
      else:
        plt.plot(datax2,datay2,color='red',label=label2,linewidth=2.0)
    if datay2_2 is not None:
      plt.plot(datax2,datay2_2,color=tableau20[19],linestyle='dashed',linewidth=2.0)
    if datay2_3 is not None:
      plt.plot(datax2,datay2_3,color='red',linestyle='dashed',linewidth=4.0)
    if datay3 is not None:
      plt.plot(datax3,datay3,color='green',label=label3,linewidth=2.0)
    if datay3_2 is not None:
      plt.plot(datax3,datay3_2,color=tableau20[19],linestyle='dashed',linewidth=2.0)
    if datay3_3 is not None:
      plt.plot(datax3,datay3_3,color='green',linestyle='dashed',linewidth=4.0)
    if datay4 is not None:
      plt.plot(datax3,datay3,color=tableau20[0],label=label4,linewidth=4.0)
    if datay5 is not None:
      plt.plot(datax3,datay3,color=tableau20[0],label=label5,linewidth=4.0)      
    if label1 is not None:
      plt.legend(loc='upper left',frameon=False)
    axes.set_aspect('auto')
    axes.yaxis.set_tick_params(labelsize=20)
    axes.xaxis.set_tick_params(labelsize=20)
    plt.savefig(filename)
    plt.close(plots)


def PlotTwoData(data, data2, xmin, xmax, ymin, ymax, labels, num, filename, xlabel='$r/r_g$'):
    plots, axes = plt.subplots(num,2,figsize=(12,10),dpi=300)
    plt.subplots_adjust(left=0.12,right=0.97,top=0.95,bottom=0.1,hspace=0.08,wspace=0.32)
    plt.xlabel(xlabel, size = 30)
 

    axes[0,0].plot(data[0],data[2],linewidth=2.0,color='black')
 #   axes[0,0].plot(data[0],data[1],linewidth=2.0,color='red',linestyle='dashed')
    axes[0,0].plot(data[0],data[6],linewidth=2.0,color='black',linestyle='dashed')
    axes[0,0].set_ylabel(labels[0], size=25.0)
    axes[0,0].set_ylim([ymin[0,0],ymax[0,0]])
    axes[0,0].yaxis.set_tick_params(labelsize=20)
#    axes[0,0].set_xlabel(xlabel, size=30)
    axes[0,0].set_xlim([xmin[0],xmax[0]])
    axes[0,0].set_xticklabels([])
    axes[0,0].annotate('$\\rm AGNM2$', xy=(0.08, 0.8), size=20,xycoords='axes fraction')

#    plt.legend(loc='upper left',frameon=False)
    
    if num > 1:
        axes[1,0].plot(data[0],data[3],linewidth=2.0,color='black')
 #       if logscale1 > 0:
 #         axes[1,0].set_yscale('log')
        axes[1,0].set_ylabel(labels[1], size=25)
        axes[1,0].set_ylim([ymin[1,0],ymax[1,0]])
        axes[1,0].yaxis.set_tick_params(labelsize=20)    
  #      axes[1,0].set_xlabel(xlabel, size=30)
        axes[1,0].set_xlim([xmin[0],xmax[0]]) 
  #      axes[1,0].xaxis.set_tick_params(labelsize=20) 
        axes[1,0].set_xticklabels([])

    if num > 2:
        axes[2,0].plot(data[0],data[4],linewidth=2.0,color='black')
#        if logscale1 > 0:
#          axes[2,0].set_yscale('log')
        axes[2,0].set_ylabel(labels[2], size=20)
        axes[2,0].set_ylim([ymin[2,0],ymax[2,0]])
        axes[2,0].yaxis.set_tick_params(labelsize=20) 
        axes[2,0].xaxis.set_tick_params(labelsize=20)  
 #       axes[2,0].set_xlabel(xlabel, size=30)
        axes[2,0].set_xlim([xmin[0],xmax[0]])
        axes[2,0].set_xticklabels([])

    if num > 3:
        axes[3,0].plot(data[0],data[5],linewidth=2.0,color='black')
#        if logscale1 > 0:
#          axes[2,0].set_yscale('log')
        axes[3,0].set_ylabel(labels[4], size=25)
        axes[3,0].set_ylim([ymin[3,0],ymax[3,0]])
        axes[3,0].yaxis.set_tick_params(labelsize=20)  
        axes[3,0].xaxis.set_tick_params(labelsize=20) 
        axes[3,0].set_xlabel(xlabel, size=30)
        axes[3,0].set_xlim([xmin[0],xmax[0]]) 

    
    axes[0,1].plot(data2[0],data2[2],linewidth=2.0,color='black')
    axes[0,1].plot(data2[0],data2[6],linewidth=2.0,color='black',linestyle='dashed')
    axes[0,1].set_ylabel(labels[0], size=25)
    axes[0,1].set_ylim([ymin[0,1],ymax[0,1]])
    axes[0,1].yaxis.set_tick_params(labelsize=20)
#    axes[0,1].set_xlabel(xlabel, size=30)
    axes[0,1].set_xlim([xmin[1],xmax[1]])
#    axes[0,1].xaxis.set_tick_params(labelsize=20) 
    axes[0,1].set_xticklabels([])
#    plt.legend(loc='upper left',frameon=False)
    axes[0,1].annotate('$\\rm AGNBM1$', xy=(0.08, 0.8), size=20,xycoords='axes fraction')
    
    if num > 1:
        axes[1,1].plot(data2[0],data2[3],linewidth=2.0,color='black')
 #       if logscale1 > 0:
 #         axes[1,0].set_yscale('log')
        axes[1,1].set_ylabel(labels[1], size=25)
        axes[1,1].set_ylim([ymin[1,1],ymax[1,1]])
        axes[1,1].yaxis.set_tick_params(labelsize=20)    
 #       axes[1,1].set_xlabel(xlabel, size=30)
        axes[1,1].set_xlim([xmin[1],xmax[1]]) 
#        axes[1,1].xaxis.set_tick_params(labelsize=20) 
        axes[1,1].set_xticklabels([])

    if num > 2:
        axes[2,1].plot(data2[0],data2[4],linewidth=2.0,color='black')
#        if logscale1 > 0:
#          axes[2,0].set_yscale('log')
        axes[2,1].set_ylabel(labels[3], size=25)
        axes[2,1].set_ylim([ymin[2,1],ymax[2,1]])
        axes[2,1].yaxis.set_tick_params(labelsize=20) 
        axes[2,1].xaxis.set_tick_params(labelsize=20)  
        axes[2,1].set_xlabel(xlabel, size=30)
 #       axes[2,1].set_xlabel(xlabel, size=30)
        axes[2,1].set_xlim([xmin[1],xmax[1]])
        axes[2,1].set_xticklabels([])

    if num > 3:
        axes[3,1].plot(data2[0],data2[5],linewidth=2.0,color='black')
#        if logscale1 > 0:
#          axes[2,0].set_yscale('log')
        axes[3,1].set_ylabel(labels[4], size=25)
        axes[3,1].set_ylim([ymin[3,1],ymax[3,1]])
        axes[3,1].yaxis.set_tick_params(labelsize=20)  
        axes[3,1].xaxis.set_tick_params(labelsize=20)  
        axes[3,1].set_xlabel(xlabel, size=30)
        axes[3,1].set_xlim([xmin[1],xmax[1]])

    plt.savefig(filename)
    plt.clf()

        
def PlotProfile3(datax,datay1,datay2,datay3,xmin,xmax,ymin,ymax,ylabel1,ylabel2,ylabel3,filename,logscale1=0,logscale2=0,logscale3=0,var4=None):
    plots, axes = plt.subplots(3,1,figsize=(9,11),dpi=300,sharex=True)
    plt.subplots_adjust(left=0.16,right=0.95,top=0.95,bottom=0.1,hspace=0.05)
    plt.xlabel('$ r/r_s $', size = 30)
    plt.xlim([xmin,xmax])
    axes[0].plot(datax,datay1,linewidth=2.0,color='black')
    if logscale1 > 0:
      axes[0].set_yscale('log')
    axes[0].set_ylabel(ylabel1, size=30)
    axes[0].set_ylim([ymin[0],ymax[0]])
    
    axes[0].yaxis.set_tick_params(labelsize=20)


    axes[1].plot(datax,datay2,linewidth=2.0,color='black')
    if logscale2 > 0:
      axes[1].set_yscale('log')
    axes[1].set_ylabel(ylabel2, size=30)

    axes[1].yaxis.set_tick_params(labelsize=20)

    axes[1].set_ylim([ymin[1],ymax[1]])

    axes[2].plot(datax,datay3,linewidth=2.0,color='black')
    if var4 is not None:
      axes[2].plot(datax,var4,linewidth=2.0,color='blue')
    if logscale3 > 0:
      axes[2].set_yscale('log')
    axes[2].set_ylabel(ylabel3, size=30)

    axes[2].set_ylim([ymin[2],ymax[2]])
    axes[2].yaxis.set_tick_params(labelsize=20)

    axes[2].xaxis.set_tick_params(labelsize=20)

    plt.savefig(filename)
    plt.clf()


def PlotProfile4(datax,datay1,datay1_2,datay2,datay2_2,datay3,datay3_2,datay4,datay4_2,xmin,xmax,ymin,ymax,ylabel1,ylabel2,ylabel3,ylabel4,filename,logscale1=0,logscale2=0,logscale3=0,logscale4=0):
    plots, axes = plt.subplots(4,1,figsize=(9,11),dpi=300,sharex=True)
    plt.subplots_adjust(left=0.16,right=0.95,top=0.95,bottom=0.1,hspace=0.05)
    plt.xlabel('$ r/r_s $', size = 30)
    plt.xlim([xmin,xmax])
    axes[0].set_xscale('log')
    axes[0].xaxis.set_tick_params(labelsize=20)
    axes[0].plot(datax,datay1,linewidth=2.0,color='black')
    axes[0].plot(datax,datay1_2,linewidth=2.0,color='blue')
    if logscale1 > 0:
      axes[0].set_yscale('log')
    axes[0].set_ylabel(ylabel1, size=30)
    axes[0].set_ylim([ymin[0],ymax[0]])
    
    axes[0].yaxis.set_tick_params(labelsize=20)


    axes[1].plot(datax,datay2,linewidth=2.0,color='black')
    axes[1].plot(datax,datay2_2,linewidth=2.0,color='blue')
    if logscale2 > 0:
      axes[1].set_yscale('log')
    axes[1].set_ylabel(ylabel2, size=30)

    axes[1].yaxis.set_tick_params(labelsize=20)

    axes[1].set_ylim([ymin[1],ymax[1]])

    axes[2].plot(datax,datay3,linewidth=2.0,color='black')
    axes[2].plot(datax,datay3_2,linewidth=2.0,color='blue')
    if logscale3 > 0:
      axes[2].set_yscale('log')
    axes[2].set_ylabel(ylabel3, size=30)

    axes[2].set_ylim([ymin[2],ymax[2]])
    axes[2].yaxis.set_tick_params(labelsize=20)

    axes[3].plot(datax,datay4,linewidth=2.0,color='black')
    axes[3].plot(datax,datay4_2,linewidth=2.0,color='blue')
    if logscale3 > 0:
      axes[3].set_yscale('log')
    axes[3].set_ylabel(ylabel4, size=30)

    axes[3].set_ylim([ymin[3],ymax[3]])
    axes[3].yaxis.set_tick_params(labelsize=20)


    plt.savefig(filename)
    plt.clf()


def Compare4Plot(datax1,datay1,datax2,datay2,datax3,datay3,datax4,datay4,xmin,xmax,ymin,ymax,ylabel1,ylabel2,ylabel3,ylabel4,filename,datay1_2=None,datay2_2=None,datay3_2=None,datay4_2=None,datay1_3=None,datay2_3=None,datay3_3=None,datay4_3=None,datay2label=None,datay2label2=None,datay2label3=None,logscale1=0,logscale2=0,logscale3=0,logscale4=0,xlabel='$ t/(10^4\\times  t_0) $',nameflag=1,namex=0.08,namey=0.2,sharex=False,title=None,leftedge=0.16):
    plots, axes = plt.subplots(4,1,figsize=(9,11),dpi=300,sharex=sharex)
    plt.subplots_adjust(left=leftedge,right=0.95,top=0.95,bottom=0.1,hspace=0.2)
    plt.xlabel(xlabel, size = 30)
    if title is not None:
      plt.title(title,size=25,y=4.65)

#    axes[0].set_xscale('log')
    axes[0].xaxis.set_tick_params(labelsize=20)
    axes[0].plot(datax1,datay1,linewidth=2.0,color='black')
    if datay1_2 is not None:
       axes[0].plot(datax1,datay1_2,linewidth=2.0,color='black',linestyle='dashed')
    if datay1_3 is not None:
       axes[0].plot(datax1,datay1_3,linewidth=2.0,color='black',linestyle='dotted')
    if logscale1 > 0:
      axes[0].set_yscale('log')
    axes[0].set_ylabel(ylabel1, size=30)
    axes[0].set_xlim([xmin[0],xmax[0]])
    axes[0].set_ylim([ymin[0],ymax[0]])
    axes[0].xaxis.set_tick_params(labelsize=20)
    axes[0].yaxis.set_tick_params(labelsize=20)
    if nameflag > 0:
      axes[0].annotate('$\\rm AGNM1$', xy=(namex, namey), size=20,xycoords='axes fraction')


    axes[1].plot(datax2,datay2,linewidth=2.0,color='black',label=datay2label)
    if datay2_2 is not None:
       axes[1].plot(datax2,datay2_2,linewidth=2.0,color='black',linestyle='dashed',label=datay2label2)
    if datay2_3 is not None:
       axes[1].plot(datax2,datay2_3,linewidth=2.0,color='black',linestyle='dotted',label=datay2label3)
    if datay2_2 is not None and datay2_3 is not None and datay1_3 is None:
       axes[1].legend(loc='center left', bbox_to_anchor=(0.005, 0.6),frameon=False)
    if logscale2 > 0:
      axes[1].set_yscale('log')
    axes[1].set_ylabel(ylabel2, size=30)
    axes[1].xaxis.set_tick_params(labelsize=20)
    axes[1].yaxis.set_tick_params(labelsize=20)
    axes[1].set_xlim([xmin[1],xmax[1]])
    axes[1].set_ylim([ymin[1],ymax[1]])
    if nameflag > 0:
      axes[1].annotate('$\\rm AGNM2$', xy=(namex, namey), size=20,xycoords='axes fraction')


    axes[2].plot(datax3,datay3,linewidth=2.0,color='black')
    if datay3_2 is not None:
      axes[2].plot(datax3,datay3_2,linewidth=2.0,color='black',linestyle='dashed')
    if datay3_3 is not None:
      axes[2].plot(datax3,datay3_3,linewidth=2.0,color='black',linestyle='dotted')
    if logscale3 > 0:
      axes[2].set_yscale('log')
    axes[2].set_ylabel(ylabel3, size=30)
    axes[2].set_xlim([xmin[2],xmax[2]])
    axes[2].set_ylim([ymin[2],ymax[2]])
    axes[2].xaxis.set_tick_params(labelsize=20)
    axes[2].yaxis.set_tick_params(labelsize=20)
    if nameflag > 0:
      axes[2].annotate('$\\rm AGNBM1$', xy=(namex, namey), size=20,xycoords='axes fraction')


    axes[3].plot(datax4,datay4,linewidth=2.0,color='black')
    if datay4_2 is not None:
       axes[3].plot(datax4,datay4_2,linewidth=2.0,color='black',linestyle='dashed')
    if datay4_3 is not None:
       axes[3].plot(datax4,datay4_3,linewidth=2.0,color='black',linestyle='dotted')
    if logscale4 > 0:
      axes[3].set_yscale('log')
    axes[3].set_ylabel(ylabel4, size=30)
    axes[3].set_xlim([xmin[3],xmax[3]])
    axes[3].set_ylim([ymin[3],ymax[3]])
    axes[3].xaxis.set_tick_params(labelsize=20)
    axes[3].yaxis.set_tick_params(labelsize=20)
    if nameflag > 0:
      axes[3].annotate('$\\rm AGNBM2$', xy=(namex, namey), size=20,xycoords='axes fraction')

    plt.savefig(filename)
    plt.clf()



    
def PlotSTdiag(data, num, time, height, min, max, varname, filename, scale=1, scale1=0, scale2=0, line=0, tautop=None, taubot=None,title=None, cticks1=None,cticks1_label=None,cticks2=None, cticks2_label=None, cticks3=None, cticks3_label=None):
    plots, axes = plt.subplots(num,1,figsize=(9,7),dpi=300,sharex=True) 
    plt.xlabel('$t/t_0$', size = 20)
    plt.subplots_adjust(left=0.1,right=0.85,top=0.95,bottom=0.1,hspace=0.1)
    if line > 0:
       axes[0].plot(time,tautop,color='r',linewidth=2.0)
       axes[0].plot(time,taubot,color='r',linewidth=2.0)
    if title is not None:
       axes[0].set_title(title)
    xmin=np.min(time)
    xmax=np.max(time)
    ymin=np.min(np.log10(height))
    ymax=np.max(np.log10(height))

#    The first panel
    if scale > 0:
       im = axes[0].imshow(data[0], cmap='RdGy_r', norm=LogNorm(vmin = min[0], vmax = max[0]), \
            origin='lower', extent=[xmin, xmax, ymin, ymax])
    else:
       im = axes[0].imshow(data[0], cmap='RdGy_r', vmin = min[0], vmax = max[0], \
            origin='lower', extent=[xmin, xmax, ymin, ymax])
            
    axes[0].set_aspect('auto')
    axes[0].set_ylabel('$\\log(r/r_{\\odot})$', size = 20)
    axes[0].yaxis.set_tick_params(labelsize=15)
            
    cbaxes = plots.add_axes([0.86, 0.69, 0.03, 0.25])

    if cticks1 is not None:
      cbar = plots.colorbar(im, cax = cbaxes, ticks=cticks1)
    else:
      cbar = plots.colorbar(im, cax = cbaxes)

    cbar.set_label(varname[0], size = 20)

    if cticks1_label is not None:
      cbar.ax.set_yticklabels(cticks1_label)

# The second panel
    if num > 1:
       if scale1 > 0:
          im = axes[1].imshow(data[1], cmap='RdGy_r', norm=LogNorm(vmin = min[1], vmax = max[1]), \
               origin='lower', extent=[xmin, xmax, ymin, ymax])
       else:
          im = axes[1].imshow(data[1], cmap='RdGy_r', vmin = min[1], vmax = max[1], \
               origin='lower', extent=[xmin, xmax, ymin, ymax])
       axes[1].set_aspect('auto')
       axes[1].set_ylabel('$\\log(r/r_{\\odot})$', size = 20)
       axes[1].yaxis.set_tick_params(labelsize=15)
       axes[1].xaxis.set_tick_params(labelsize=15)
            
       cbaxes = plots.add_axes([0.86, 0.4, 0.03, 0.25])

       if cticks2 is not None:
         cbar = plots.colorbar(im, cax = cbaxes, ticks=cticks2)
       else:
         cbar = plots.colorbar(im, cax = cbaxes)

       cbar.set_label(varname[1], size = 20)

       if cticks2_label is not None:
         cbar.ax.set_yticklabels(cticks2_label)

# The third panel
    if num > 2:
       if scale2 > 0:
          im = axes[2].imshow(data[2], cmap='RdGy_r', norm=LogNorm(vmin = min[2], vmax = max[2]), \
               origin='lower', extent=[xmin, xmax, ymin, ymax])
       else:
          im = axes[2].imshow(data[2], cmap='RdGy_r', vmin = min[2], vmax = max[2], \
               origin='lower', extent=[xmin, xmax, ymin, ymax])
       axes[2].set_aspect('auto')
       axes[2].set_ylabel('$\\log(r/r_{\\odot})$', size = 20)
       axes[2].yaxis.set_tick_params(labelsize=15)
       axes[2].xaxis.set_tick_params(labelsize=15)
            
       cbaxes = plots.add_axes([0.86, 0.11 , 0.03, 0.25])
       
       if cticks3 is not None:
         cbar = plots.colorbar(im, cax = cbaxes, ticks=cticks3)
       else:
         cbar = plots.colorbar(im, cax = cbaxes)

       cbar.set_label(varname[2], size = 20)

       if cticks3_label is not None:
         cbar.ax.set_yticklabels(cticks3_label)

    
   
    plt.savefig(filename)
    plt.clf()



def PlotAveStructure(filename,rhomin,rhomax,Ermin,Ermax,PBmin,PBmax,Bin,Bmax,Frmin,Frmax,outputname1,outputname2,outputname3=None,outputname4=None,vecname1='Fr1',vecname2='Fr2',vecname3='RhoVr',vecname4='RhoVtheta'):
    f=h5py.File(filename, 'r')



    rho=f['rho'].value
    x1v=f['x1v'].value
    x2v=f['x2v'].value
    x1f=f['x1f'].value
    x2f=f['x2f'].value
    vel1=f[vecname3].value
    vel2=f[vecname4].value
    Fr1=f[vecname1].value
    Fr2=f[vecname2].value
    Er=f['Er'].value
    B1=f['B1'].value
    B2=f['B2'].value
    PB=f['PB'].value
    sigma_s=f['Sigma_s'].value
    sigma_a=f['Sigma_a'].value

    sigma=sigma_s+sigma_a

    kappa=sigma/(rho*kappaes)

    if vecname3=='RhoVr':
      vel1=vel1/rho
    
    if vecname4=='RhoVtheta':
      vel2=vel2/rho

    #Get location of tau=1 from rotation axis

    #convert to 1D
    nr=len(x1v)
    ntheta=len(x2v)
    dx1=np.zeros(nr)
    cellvol=np.zeros((ntheta,nr))
    rho_1D=np.zeros(nr*ntheta)
    vel1_1D=np.zeros(nr*ntheta)
    vel2_1D=np.zeros(nr*ntheta)
    Fr1_1D=np.zeros(nr*ntheta)
    Fr2_1D=np.zeros(nr*ntheta)
    Er_1D=np.zeros(nr*ntheta)
    radius=np.zeros(nr*ntheta)
    theta=np.zeros(nr*ntheta)
    B1_1D=np.zeros(nr*ntheta)
    B2_1D=np.zeros(nr*ntheta)
    PB_1D=np.zeros(nr*ntheta)
    sigma_1D=np.zeros(nr*ntheta)
    kappa_1D=np.zeros(nr*ntheta)
    dx_1D=np.zeros(nr*ntheta)

    for i in range(nr):
       dx1[i]=x1f[i+1]-x1f[i]

    for j in range(ntheta):
      for i in range(nr):
        rho_1D[j*nr+i]=rho[j,i]
        vel1_1D[j*nr+i]=vel1[j,i]
        vel2_1D[j*nr+i]=vel2[j,i]
        Fr1_1D[j*nr+i]=Fr1[j,i]
        Fr2_1D[j*nr+i]=Fr2[j,i]
        Er_1D[j*nr+i]=Er[j,i]
        B1_1D[j*nr+i]=B1[j,i]
        B2_1D[j*nr+i]=B2[j,i]
        PB_1D[j*nr+i]=PB[j,i]
        sigma_1D[j*nr+i]=sigma[j,i]
        kappa_1D[j*nr+i]=kappa[j,i]
        cellvol[j,i]=vol_func(x1f[i],x1f[i+1],x2f[j],x2f[j+1])
        radius[j*nr+i]=2*x1v[i]
        theta[j*nr+i]=x2v[j]
        dx_1D[j*nr+i]=dx1[i]*np.sin(x2v[j])




    vx=vel1_1D*np.sin(theta)+vel2_1D*np.cos(theta)
    vy=vel1_1D*np.cos(theta)-vel2_1D*np.sin(theta)

    Bx=B1_1D*np.sin(theta)+B2_1D*np.cos(theta)
    By=B1_1D*np.cos(theta)-B2_1D*np.sin(theta)


    Frx=Fr1_1D*np.sin(theta)+Fr2_1D*np.cos(theta)
    Fry=Fr1_1D*np.cos(theta)-Fr2_1D*np.sin(theta)

    # convert to cartesian coordinate
    xcoord=radius*np.sin(theta)
    ycoord=radius*np.cos(theta)


    width=50
    height=100
    nx=256
    ny=int(2*height*nx/width)

    xmin=0
    rmin=4
    xmax=width
    ymin=-height
    ymax=height

    xgrid=np.linspace(xmin, xmax, nx)
    ygrid=np.linspace(ymin, ymax, ny)

    xmesh,ymesh=np.meshgrid(xgrid,ygrid)

    rho_cart=griddata(np.c_[xcoord,ycoord],rho_1D,(xmesh,ymesh),method='nearest')
    sigma_cart=griddata(np.c_[xcoord,ycoord],sigma_1D,(xmesh,ymesh),method='nearest')
    kappa_cart=griddata(np.c_[xcoord,ycoord],kappa_1D,(xmesh,ymesh),method='nearest')   
    dx_cart=griddata(np.c_[xcoord,ycoord],dx_1D,(xmesh,ymesh),method='nearest')

    vx_cart=griddata(np.c_[xcoord,ycoord],vx,(xmesh,ymesh),method='nearest')
    vy_cart=griddata(np.c_[xcoord,ycoord],vy,(xmesh,ymesh),method='nearest')

    vx_cart=vx_cart/crat
    vy_cart=vy_cart/crat

    PB_cart=griddata(np.c_[xcoord,ycoord],PB_1D,(xmesh,ymesh),method='nearest')

    Bx_cart=griddata(np.c_[xcoord,ycoord],Bx,(xmesh,ymesh),method='nearest')
    By_cart=griddata(np.c_[xcoord,ycoord],By,(xmesh,ymesh),method='nearest')

    Er_cart=griddata(np.c_[xcoord,ycoord],Er_1D,(xmesh,ymesh),method='nearest')

    Frx_cart=griddata(np.c_[xcoord,ycoord],Frx,(xmesh,ymesh),method='nearest')
    Fry_cart=griddata(np.c_[xcoord,ycoord],Fry,(xmesh,ymesh),method='nearest')

#calculate the optical depth from rotation axis
    tau=np.zeros((ny,nx))
    xtaupos=np.zeros(ny)
    for j in range(ny):
       xtaupos[j]=xmin
       tau[j,0] = 0.5*dx_cart[j,0]*sigma_cart[j,0]
       if tau[j,0] < 1:
          xtaupos[j] = xgrid[0]
       for i in range(1,nx):
           tau[j,i] = tau[j,i-1] + 0.5*dx_cart[j,i]*sigma_cart[j,i]
           if tau[j,i] < 1:
              xtaupos[j] = xgrid[i]


    rmesh=(xmesh**2.0+ymesh**2.0)**0.5

    rindix=rmesh < rmin

    rho_cart[rindix]=rhomin
    vx_cart[rindix]=0.0
    vy_cart[rindix]=0.0
    Er_cart[rindix]=Ermin
    Bx_cart[rindix]=0.0
    By_cart[rindix]=0.0
    Frx_cart[rindix]=0.0
    Fry_cart[rindix]=0.0
    kappa_cart[rindix]=0.0

    vlim1=1.e-3
    vlim2=1

    label1='$\\rho/\\rho_0$'
    label2='$v/c$'

 

    MakeRhoVSlice(rho_cart, vx_cart, vy_cart, rhomin, rhomax, vlim1,vlim2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,outputname1,1,taux=xtaupos,tauy=ygrid)



    label1='$E_r/a_rT_0^4$'
    label2='$B/B_0$'


    vlim1=Bmin
    vlim2=Bmax

    MakeRhoVSlice(Er_cart, Bx_cart, By_cart, Ermin,Ermax, vlim1,vlim2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,outputname2,1)

    if outputname3 is not None:

      label1='$E_r/a_rT_0^4$'
      label2='$F_r/ca_rT_0^4$'


      vlim1=Frmin
      vlim2=Frmax

      MakeRhoVSlice(Er_cart, Frx_cart, Fry_cart, Ermin,Ermax, vlim1,vlim2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,outputname3,1)
      
    if outputname4 is not None:

      label1='$\\kappa/\kappa_0$'
      label2='$v/c$'
      
      MakeRhoVSlice(kappa_cart, vx_cart, vy_cart, 0.2,3, 1.e-3, 1, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,outputname4,1,logscale=0)


  #  outputname='AGNB2_Br.png'

  #  label1='$P_B/P_0$'
  #  label2='$B/B_0$'


  #  vlim1=1.e-5
  #  vlim2=0.01

  #  MakeRhoVSlice(PB_cart, Bx_cart, By_cart, PBmin,PBmax, vlim1,vlim2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,outputname,1)

def ReadAvedata(filename,thetalim1=None,thetalim2=None):


    file=h5py.File(filename,'r')


    x1v=file['x1v'].value # in unit of r_s
    x1v=x1v*2
    x2v=file['x2v'].value
    x1f=file['x1f'].value  # in unit of r_s
    x1f=x1f*2
    x2f=file['x2f'].value

    # find the location of theta limit
    if thetalim1 is not None:
      if thetalim1 > 0.0:
        nt1=np.abs(x2v-thetalim1).argmin()
      else:
        nt1=np.abs(x2v+thetalim1).argmin()        

    if thetalim2 is not None:
      if thetalim2 > 0.0:
        nt2=np.abs(x2v-thetalim2).argmin()
      else:
        nt2=np.abs(x2v+thetalim2).argmin()

    gamma=5.0/3.0

    #
    nr=len(x1v)
    ntheta=len(x2v)

    cellvol=np.zeros((ntheta,nr))
    x1area=np.zeros((ntheta,nr))
    radius=np.zeros((ntheta,nr))
    theta=np.zeros((ntheta,nr))

    for j in range(0,ntheta):
      for i in range(0,nr):
        x1area[j,i]=rarea_func(x1f[i+1],x2f[j],x2f[j+1])
        cellvol[j,i]=vol_func(x1f[i],x1f[i+1],x2f[j],x2f[j+1])
        radius[j,i]=x1v[i]
        theta[j,i]=x2v[j]


    rhovrvphi=file['rhovrvphi'].value
    BrBphi=file['BrBphi'].value
    Pr13=prat*file['Pr13'].value
    meanrhovrvphi=file['meanrhovrvphi'].value
    rhovr=file['RhoVr'].value
    vphi_rho=file['meanvphi'].value
    drhovrvphi=rhovrvphi-rhovr*vphi_rho

    
    if thetalim1 is not None:
      rhovrvphi=rhovrvphi[nt1:nt2,:]
      BrBphi=BrBphi[nt1:nt2,:]
      Pr13=Pr13[nt1:nt2,:]
      meanrhovrvphi=meanrhovrvphi[nt1:nt2,:]
      rhovr=rhovr[nt1:nt2,:]
      vphi_rho=vphi_rho[nt1:nt2,:]
      drhovrvphi=drhovrvphi[nt1:nt2,:]
      cellvol=cellvol[nt1:nt2,:]
      theta=theta[nt1:nt2,:]

    hydroflux=np.average(rhovrvphi*np.sin(theta),axis=0,weights=cellvol)*(x1v**3.0)
    maxflux=np.average(BrBphi*np.sin(theta),axis=0,weights=cellvol)*(x1v**3.0)
    Prflux=np.average(Pr13*np.sin(theta),axis=0,weights=cellvol)*(x1v**3.0)
    Reynflux=np.average(drhovrvphi*np.sin(theta),axis=0,weights=cellvol)*(x1v**3.0)
    meanlflux=np.average(rhovr*vphi_rho*np.sin(theta),axis=0,weights=cellvol)*(x1v**3.0)
    sumflux=hydroflux+maxflux+Prflux

    rho=file['rho'].value
    Er=prat*file['Er'].value
    Pr=Er/3.0
    pgas=file['pgas'].value
    PB=file['PB'].value
    rhovin=file['RhoVin'].value
    rhovout=file['RhoVout'].value
    Sgas=np.log(pgas/rho**gamma)/(gamma-1.0)
    Br=file['B1'].value
    Btheta=file['B2'].value
    Bphi=file['B3'].value

    if thetalim1 is not None:
      rho=rho[nt1:nt2,:]
      Pr=Pr[nt1:nt2,:]
      pgas=pgas[nt1:nt2,:]
      PB=PB[nt1:nt2,:]
      rhovin=rhovin[nt1:nt2,:]
      rhovout=rhovout[nt1:nt2,:]
      Sgas=Sgas[nt1:nt2,:]

   # the area is in unit of r_g^2, need to convert to r_s^2 to match the unit of MEdd

    if thetalim1 is not None:
      massout=np.sum(rhovout*x1area[nt1:nt2,:],axis=0)*0.25/Medd
      massin=np.sum(rhovin*x1area[nt1:nt2,:],axis=0)*0.25/Medd
      totmass=np.sum(rhovr*x1area[nt1:nt2,:],axis=0)*0.25/Medd
    else:
      massout=np.sum(rhovout*x1area,axis=0)*0.25/Medd
      massin=np.sum(rhovin*x1area,axis=0)*0.25/Medd
      totmass=np.sum(rhovr*x1area,axis=0)*0.25/Medd     

    
    meanrho=np.average(rho,axis=0,weights=cellvol)
    meanPr=np.average(Pr,axis=0,weights=cellvol)
    meanpgas=np.average(pgas,axis=0,weights=cellvol)
    meanPB=np.average(PB,axis=0,weights=cellvol)
    meanrhovin=np.average(rhovin,axis=0,weights=cellvol)
    meanrhovout=np.average(rhovout,axis=0,weights=cellvol)
    meanrhov=np.average(rhovr,axis=0,weights=cellvol)
    meanReyn=np.average(drhovrvphi,axis=0,weights=cellvol)
    meanmax=np.average(-BrBphi,axis=0,weights=cellvol)
    meanPrphi=np.average(Pr13,axis=0,weights=cellvol)
    meanSgas=np.average(Sgas,axis=0,weights=cellvol*rho)
    meanBr=np.average(Br,axis=0,weights=cellvol)
    meanBtheta=np.average(Btheta,axis=0,weights=cellvol)
    meanBphi=np.average(Bphi,axis=0,weights=cellvol)

    returndata=[x1v,sumflux,hydroflux,Reynflux,maxflux,Prflux,meanlflux,meanrho,meanPr, \
                 meanpgas,meanPB,meanrhovin,meanrhovout,meanReyn,meanmax,meanPrphi,meanSgas,\
                 massin, massout,totmass,meanBr,meanBtheta,meanBphi]
    return returndata



def ReadVerticaldata(filename,rpos1=10,rpos2=20,rpos3=30,Frflag=0):


    file=h5py.File(filename,'r')


    x1v=file['x1v'].value # in unit of r_s
    x1v=x1v*2
    x2v=file['x2v'].value
    x1f=file['x1f'].value  # in unit of r_s
    x1f=x1f*2
    x2f=file['x2f'].value

    npros=np.zeros(3)
    radius=np.zeros(3)

    # find the location of radial position
    if rpos1 is not None:
        npros[0]=np.abs(x1v-rpos1).argmin()

    if rpos2 is not None:
        npros[1]=np.abs(x1v-rpos2).argmin()

    if rpos3 is not None:
        npros[2]=np.abs(x1v-rpos3).argmin()

    radius[0]=rpos1
    radius[1]=rpos2
    radius[2]=rpos3


    gamma=5.0/3.0

    #
    nr=len(x1v)
    ntheta=len(x2v)


    rho=file['rho'].value
    Er=prat*file['Er'].value
    Pr=Er/3.0
    pgas=file['pgas'].value
    PB=file['PB'].value
    PB2=file['PB2'].value
    BrBphi=file['BrBphi'].value
    Pr13=file['Pr13'].value
    rhovrvphi=file['rhovrvphi'].value
    meanrhovrvphi=file['meanrhovrvphi'].value
    rhovr=file['RhoVr'].value
    vphi_rho=file['meanvphi'].value
    Sgas=np.log(pgas/rho**gamma)/(gamma-1.0)
    drhovrvphi=rhovrvphi-rhovr*vphi_rho
    kappa=(file['Sigma_a'].value+file['Sigma_s'].value)/rho
#    MRIlambda=file['lambda_a'].value

    Fr01=file['Fr01'].value
    Fr02=file['Fr02'].value

    Fr01sigma=file['Fr01Sigma'].value
    Fr02sigma=file['Fr02Sigma'].value

    if Frflag > 0:
      Fr1=file['Fr1new'].value
      Fr2=file['Fr2new'].value
    else:
      Fr1=file['Fr1'].value
      Fr2=file['Fr2'].value


    rho_pos=np.zeros((3,ntheta))
    Pr_pos=np.zeros((3,ntheta))
    pgas_pos=np.zeros((3,ntheta))
    PB_pos=np.zeros((3,ntheta))

    BrBphi_pos=np.zeros((3,ntheta))
    Pr13_pos=np.zeros((3,ntheta))
    rhovrvphi_pos=np.zeros((3,ntheta))
    meanrhovrvphi_pos=np.zeros((3,ntheta))
    rhovr_pos=np.zeros((3,ntheta))
    vphi_rho_pos=np.zeros((3,ntheta))
    Sgas_pos=np.zeros((3,ntheta))
    drhovrvphi_pos=np.zeros((3,ntheta))
    lambdaA_pos=np.zeros((3,ntheta))
    ag_pos=np.zeros((3,ntheta))
    arad_pos=np.zeros((3,ntheta))
    arad2_pos=np.zeros((3,ntheta))
    MRIlambda_pos=np.zeros((3,ntheta))
    Frz0_pos=np.zeros((3,ntheta))
    Frz_pos=np.zeros((3,ntheta))



    for num in range(3):
      frz0=Fr01[:,int(npros[num])]*np.cos(x2v)-Fr02[:,int(npros[num])]*np.sin(x2v)
      Frz0_pos[num,:]=frz0
      Frz_pos[num,:]=Fr1[:,int(npros[num])]*np.cos(x2v)-Fr2[:,int(npros[num])]*np.sin(x2v)
      frz0sigma=Fr01sigma[:,int(npros[num])]*np.cos(x2v)-Fr02sigma[:,int(npros[num])]*np.sin(x2v)
      rho_pos[num,:]=rho[:,int(npros[num])]
      Pr_pos[num,:]=Pr[:,int(npros[num])]
      pgas_pos[num,:]=pgas[:,int(npros[num])]
      PB_pos[num,:]=PB[:,int(npros[num])]
      BrBphi_pos[num,:]=BrBphi[:,int(npros[num])]
      Pr13_pos[num,:]=Pr13[:,int(npros[num])]
      rhovrvphi_pos[num,:]=rhovrvphi[:,int(npros[num])]
      meanrhovrvphi_pos[num,:]=meanrhovrvphi[:,int(npros[num])]
      rhovr_pos[num,:]=rhovr[:,int(npros[num])]
      vphi_rho_pos[num,:]=vphi_rho[:,int(npros[num])]
      Sgas_pos[num,:]=Sgas[:,int(npros[num])]
      drhovrvphi_pos[num,:]=drhovrvphi[:,int(npros[num])]
      lambdaA_pos[num,:]=(2.0*PB2[:,int(npros[num])]/rho_pos[num,:])**0.5*2.0*np.pi/vphi_rho_pos[num,:]
      ag_pos[num,:]=(np.cos(x2v)*gm/(radius[num]*0.5-1.0)**2.0)*2.0/(crat**2.0)
      arad_pos[num,:]=kappa[:,int(npros[num])]*prat*frz0*2.0/(crat**2.0)
      arad2_pos[num,:]=prat*frz0sigma*2.0/(rho_pos[num,:]*crat**2.0)
   #   MRIlambda_pos[num,:]=MRIlambda[:,int(npros[num])]*2.0*np.pi/vphi_rho_pos[num,:]



    returndata=[x2v,rho_pos,Pr_pos,pgas_pos,PB_pos,BrBphi_pos,Pr13_pos,rhovrvphi_pos,meanrhovrvphi_pos, \
                rhovr_pos,vphi_rho_pos,Sgas_pos,drhovrvphi_pos,lambdaA_pos,ag_pos,arad_pos,arad2_pos,MRIlambda_pos,Frz_pos,Frz0_pos]

    return returndata



def RadialGradient(radius,data):
  n=radius.size
  ddata=np.zeros(n)

  
  ddata[0]=(data[1]-data[0])/(radius[1]-radius[0])
  ddata[n-1]=(data[n-1]-data[n-2])/(radius[n-1]-radius[n-2])
  for i in range(1,n-1):
    ddata[i]=(data[i+1]-data[i-1])/(radius[i+1]-radius[i-1])

  return ddata



###################################################################
# The constants

# The diffusion equation is 
#\partial T/\partial t= (gamma-1)kappa/\rho Grad^2 T
#initial profile T(x,0)=exp(-40x^2)
#the solution is T(x,t)=1/(1+160t D)exp(-40x^2/(1+160D))
#where D=(gamma-1)kappa/\rho


datafile0='bin/tc.out1.00000.athdf'
datafile1='bin/tc.out1.00030.athdf'
datafile2='bin/tc.out1.00060.athdf'

data0=h5py.File(datafile0,'r')
x=data0['x1v'][0]
rho0=data0['prim'][0,0,0,4,:]
pgas0=data0['prim'][1,0,0,4,:]
t0=data0.attrs['Time']
Tg0=pgas0/rho0


data1=h5py.File(datafile1,'r')
rho1=data1['prim'][0,0,0,4,:]
pgas1=data1['prim'][1,0,0,4,:]
t1=data1.attrs['Time']
Tg1=pgas1/rho1

data2=h5py.File(datafile2,'r')
rho2=data2['prim'][0,0,0,4,:]
pgas2=data2['prim'][1,0,0,4,:]
t2=data2.attrs['Time']
Tg2=pgas2/rho2


filename='TC_diffusion.pdf'
ylabel='$T$'
# PlotProfile(datax, datay, xmin, xmax, ymin, ymax,  ylabel, label1, filename, xlabel='$r/r_g$', logscale=0, datay1_2=None, datay1_3=None, datax2=None, datay2=None, datay2_2=None, datay2_3=None, datax3=None, datay3=None, datay3_2=None, datay3_3=None, datax4=None, datay4=None, datax5=None, datay5=None, label2='', label3='', label4='', label5='',datay2line=None):
# the analytical diffusion solution is
#u(x,t)=1/(4Dt nu^2+1)exp(-nu^2 x^2/(4Dt nu^2+1))

#the initial profile is 
#u(x,0)=exp(-40x^2)
#nu^2=40

gamma=5.0/3.0
kappa=0.01
rho=1


D=(gamma-1.0)*kappa/rho
Tsol0=exp(-40*x**2)
Tsol1=exp(-40*x**2/(4*D*t1*40.0+1))/(4*D*t1*40+1)**0.5
Tsol2=exp(-40*x**2/(4*D*t2*40.0+1))/(4*D*t2*40+1)**0.5

label1='$t='+"%4.2f"%(t0)+'$'
label2='$t='+"%4.2f"%(t1)+'$'
label3='$t='+"%4.2f"%(t2)+'$'

PlotProfile(x,Tg0, -1, 1, 0, 1.0,  ylabel, label1, filename, xlabel='$x$', logscale=0,datay1_2=Tsol0,datax2=x,datay2=Tg1,datay2_2=Tsol1,datax3=x,datay3=Tg2,datay3_2=Tsol2,label2=label2,label3=label3)


