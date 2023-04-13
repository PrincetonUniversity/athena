import matplotlib
matplotlib.use('Agg')

import numpy as np

import h5py



# setup latex fonts
#rcParams['text.usetex']=True
#rcParams['text.latex.unicode']=True


#matplotlib.rc('font', family='serif', serif='cm10')
#matplotlib.rc('text', usetex=True)
#matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']

# The full quantities' names are:
#[[0]'rho' [1]'press' [2]'vel1' [3]'vel2' [4]'vel3' [5]'Bcc1' [6]'Bcc2' [7]'Bcc3' [8]'Er' [9]'Fr1' [10]'Fr2'
# [11]'Fr3' [12]'Pr11' [13]'Pr22' [14]'Pr33' [15]'Pr12' [16]'Pr13' [17]'Pr23' [18]'Pr21' [19]'Pr31' [20]'Pr32' 
# [21]'Er0'
# [22]'Fr01' [23]'Fr02' [24]'Fr03' [25]'Sigma_s' [26]'Sigma_a' [27]'Sigma_p' [28]'ir_0' [29]'ir_1' [30]'ir_2'
# [31]'ir_3' [32]'ir_4' [33]'ir_5' [34]'ir_6' [35]'ir_7' [36]'ir_8' [37]'ir_9']



def ReduceData(filename,rmax,crat):

# Read attributes and data
  f=h5py.File(filename, 'r')
  level = f.attrs['MaxLevel']
  block_size = f.attrs['MeshBlockSize']
  root_grid_size = f.attrs['RootGridSize']
  levels = f['Levels'][:]
  logical_locations = f['LogicalLocations'][:]
  nx1 = root_grid_size[0] * 2**level
  nx2 = root_grid_size[1] * 2**level if root_grid_size[1] > 1 else 1
  nx3 = root_grid_size[2] * 2**level if root_grid_size[2] > 1 else 1

  if nx3 > 1:
    dim = 3
  elif nx2 > 1:
    dim = 2
  else:
    dim = 1

  data={}

# the coordinates
  for d in range(1,4):
    nx = (nx1,nx2,nx3)[d-1]
    xmin = f.attrs['RootGridX'+repr(d)][0]
    xmax = f.attrs['RootGridX'+repr(d)][1]
    xrat_root = f.attrs['RootGridX'+repr(d)][2]
    if (xrat_root == 1.0):
      data['x'+repr(d)+'f'] = np.linspace(xmin, xmax, nx+1)
    else:
      xrat = xrat_root ** (1.0 / 2**level)
      data['x'+repr(d)+'f'] = \
         xmin + (1.0-xrat**np.arange(nx+1)) / (1.0-xrat**nx) * (xmax-xmin)


  var_quantities = np.array([x.decode('ascii', 'replace') for x in f.attrs['VariableNames'][:]])
  coord_quantities = ('x1f', 'x2f', 'x3f', 'x1v', 'x2v', 'x3v')
  attr_quantities = [key for key in f.attrs]
  other_quantities = ('Levels',)    
  quantities = var_quantities
  quantities = [str(q) for q in quantities if q not in coord_quantities 
              and q not in attr_quantities and q not in other_quantities]
  for key in attr_quantities:
    data[str(key)] = f.attrs[key]
  num_blocks = f.attrs['NumMeshBlocks']
  dataset_names = np.array([x.decode('ascii', 'replace')
                                  for x in f.attrs['DatasetNames'][:]])
  dataset_sizes = f.attrs['NumVariables'][:]
  dataset_sizes_cumulative = np.cumsum(dataset_sizes)
  variable_names = np.array([x.decode('ascii', 'replace')
                                   for x in f.attrs['VariableNames'][:]])
  quantity_datasets = []
  quantity_indices = []
  spec_datasets = []
  spec_indices = []
  for q in quantities:
    var_num = np.where(variable_names == q)[0][0]
    dataset_num = np.where(dataset_sizes_cumulative > var_num)[0][0]
    if dataset_num == 0:
          dataset_index = var_num
    else:
          dataset_index = var_num - dataset_sizes_cumulative[dataset_num-1]
    quantity_datasets.append(dataset_names[dataset_num])
    quantity_indices.append(dataset_index)
    if q == 'rho' or q=='vel1' or q=='vel2' or q=='vel3' or q=='press' or q=='Bcc1' or \
                q=='Bcc2' or q=='Bcc3' or q=='Er' or q=='Sigma_s_0' or q=='Sigma_a_0' or q=='Fr01' or q=='Fr02' \
                or q == 'Fr03' or q == 'Er0':
      spec_datasets.append(dataset_names[dataset_num])
      spec_indices.append(dataset_index)



# get rho, v1, v2, v3, B1, B2, B3, kappa_s, kappa_a


# now add the derived quantities
  quantities.append('Maxwell')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('Reynolds')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('BrBphi')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('BthetaBphi')
  quantity_datasets.append('None')
  quantity_indices.append(0)
  
  quantities.append('rhovrvphi')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('rhovthetavphi')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('vrEr')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('vthetaEr')
  quantity_datasets.append('None')
  quantity_indices.append(0)
  
  quantities.append('rhoPB')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('rhosq')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('PB1')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('PB2')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('PB3')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('PB')
  quantity_datasets.append('None')
  quantity_indices.append(0)


  quantities.append('PBsq')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('kappa_s')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('kappa_a')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('Radacc')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('RhoVr')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('RhoVtheta')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('RhoVphi')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('RhoVout')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('RhoVin')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('Ekin1')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('Ekin2')
  quantity_datasets.append('None')
  quantity_indices.append(0)
  
  quantities.append('Ekin3')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('tgas')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('Fr01Sigma')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('Fr02Sigma')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('Fr01kappa')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('Fr02kappa')
  quantity_datasets.append('None')
  quantity_indices.append(0)
  
  quantities.append('rhovrsq')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('vrsq')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('rhovphiout')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('rhovphiin')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('rhoinflow')
  quantity_datasets.append('None')
  quantity_indices.append(0)
  
  quantities.append('rhooutflow')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('rhovrvphisq')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('Ersq')
  quantity_datasets.append('None')
  quantity_indices.append(0)
  
  quantities.append('rhovrvphiEr')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('rhovrvphiPB')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('radvis')
  quantity_datasets.append('None')
  quantity_indices.append(0)
  
  quantities.append('velshear')
  quantity_datasets.append('None')
  quantity_indices.append(0)
  
  quantities.append('divvel')
  quantity_datasets.append('None')
  quantity_indices.append(0)
  
  quantities.append('raddivwork')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('radshearwork')
  quantity_datasets.append('None')
  quantity_indices.append(0)


  quantities.append('curvtrho')
  quantity_datasets.append('None')
  quantity_indices.append(0)
    

  quantities.append('curvtrhosq')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('curvgrads')
  quantity_datasets.append('None')
  quantity_indices.append(0)


  quantities.append('curvgradssq')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('Fr1new')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('Fr2new')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('Fr3new')
  quantity_datasets.append('None')
  quantity_indices.append(0)


  quantities.append('lambda_a')
  quantity_datasets.append('None')
  quantity_indices.append(0)



  n_rmax=np.argmin(np.fabs(data['x1f']-rmax))
  for q in quantities:
    data[q] = np.zeros((nx2,n_rmax))
  vol_phi=np.zeros((nx2,n_rmax))


# Go through blocks in data file
  for block_num in range(num_blocks):
    # Extract location information
    block_level = levels[block_num]
    block_location = logical_locations[block_num,:]
   # Calculate scale (number of copies per dimension)
    s= 2**(level-block_level)
  
  # the block size
    block_phi_size=f['x3f'][block_num][block_size[2]]-f['x3f'][block_num][0]
    radius=f['x1f'][block_num]
    theta=f['x2f'][block_num]
    phi=f['x3f'][block_num]

    if np.max(radius) < rmax:
  
  # get cell center coordinates
      x1v=np.zeros(block_size[0])
      x2v=np.zeros(block_size[1])
      x3v=np.zeros(block_size[2])
    
      for ni in range(block_size[0]):
        x1v[ni]=0.75*(radius[ni+1]**4.0 - radius[ni]**4.0)/(radius[ni+1]**3.0-radius[ni]**3.0);
    
      for nj in range(block_size[1]):
        x2v[nj]=((np.sin(theta[nj+1]) - theta[nj+1] * np.cos(theta[nj+1])) \
              -(np.sin(theta[nj]) - theta[nj] * np.cos(theta[nj]))) \
                / (np.cos(theta[nj]) - np.cos(theta[nj+1]));

      for nk in range(block_size[2]):
        x3v[nk] = 0.5 * (phi[nk+1]+phi[nk])
    
      grid_phi,grid_theta,grid_r=np.meshgrid(x3v,x2v,x1v,indexing='ij')
    

  # Calculate fine-level begin indices
      il = block_location[0] * block_size[0] * s
      jl = block_location[1] * block_size[1] * s if dim >= 2 else 0

  # Calculate fine-level end indices
      iu = il + block_size[0] * s
      ju = jl + block_size[1] * s if dim >= 2 else 1

  # Calculate fine-level offsets
      io_vals = range(s)
      jo_vals = range(s) if dim >= 2 else (0,)

            
      rho_data=f[spec_datasets[0]][spec_indices[0],block_num,:]
      pgas_data=f[spec_datasets[1]][spec_indices[1],block_num,:]
      v1_data=f[spec_datasets[2]][spec_indices[2],block_num,:]
      vout_data=v1_data.clip(min=0.0)
      vin_data=v1_data.clip(max=0.0)
      v2_data=f[spec_datasets[3]][spec_indices[3],block_num,:]
      v3_data=f[spec_datasets[4]][spec_indices[4],block_num,:]
      B1_data=f[spec_datasets[5]][spec_indices[5],block_num,:]
      B2_data=f[spec_datasets[6]][spec_indices[6],block_num,:]
      B3_data=f[spec_datasets[7]][spec_indices[7],block_num,:]
      Er_data=f[spec_datasets[8]][spec_indices[8],block_num,:]
      Er0_data=f[spec_datasets[9]][spec_indices[9],block_num,:]
      Fr01_data=f[spec_datasets[10]][spec_indices[10],block_num,:]
      Fr02_data=f[spec_datasets[11]][spec_indices[11],block_num,:]
      Fr03_data=f[spec_datasets[12]][spec_indices[12],block_num,:]
      sigma_s_data=f[spec_datasets[13]][spec_indices[13],block_num,:]
      sigma_a_data=f[spec_datasets[14]][spec_indices[14],block_num,:]
      PB_data=0.5*(np.multiply(B1_data,B1_data)+np.multiply(B2_data,B2_data)+np.multiply(B3_data,B3_data))
      PB1_data=0.5*np.multiply(B1_data,B1_data)   
      PB2_data=0.5*np.multiply(B2_data,B2_data)
      PB3_data=0.5*np.multiply(B3_data,B3_data)
      rhovphi_data=np.multiply(rho_data,v3_data)
      rhovrvphi_data=np.multiply(v1_data,rhovphi_data)
      tr_data=(Er_data)**0.25
      entropy_data=tr_data**3.0/rho_data

      velmag=np.sqrt(v1_data**2+v2_data**2+v3_data**2)
      betafactor=velmag/crat
      lorentz=1.0/np.sqrt(1.0-betafactor**2)
      
  #calculate velocity shear and divergence
      divvel=np.zeros((block_size[2],block_size[1],block_size[0]))
      velshear=np.zeros((block_size[2],block_size[1],block_size[0]))
      
      curvtrho=np.zeros((block_size[2],block_size[1],block_size[0]))
      curvtrhosq=np.zeros((block_size[2],block_size[1],block_size[0]))
      curvgrads=np.zeros((block_size[2],block_size[1],block_size[0]))
      curvgradssq=np.zeros((block_size[2],block_size[1],block_size[0]))
      
      kappatot=np.divide(sigma_a_data+sigma_s_data,rho_data)
      
      
      for nk in range(0,block_size[2]):
        for nj in range(0,block_size[1]):
          for ni in range(0,block_size[0]):
             #get cell center position
             
             dvrdr=0
             dvtdt=0
             dvpdp=0
             
             dvpdt=0.0
             dvtdp=0.0
             dvrdp=0.0
             dvpdr=0.0
             dvtdr=0.0
             dvrdt=0.0
             
             dsdr=0.0
             dsdt=0.0
             dsdp=0.0
             
             if ni==0:
               dvrdr=(x1v[ni+1]**2*v1_data[nk,nj,ni+1]-x1v[ni]**2*v1_data[nk,nj,ni])/(x1v[ni+1]-x1v[ni])
               dvpdr=(x1v[ni+1]*v3_data[nk,nj,ni+1]-x1v[ni]*v3_data[nk,nj,ni])/(x1v[ni+1]-x1v[ni])
               dvtdr=(x1v[ni+1]*v2_data[nk,nj,ni+1]-x1v[ni]*v2_data[nk,nj,ni])/(x1v[ni+1]-x1v[ni])
               dsdr=(entropy_data[nk,nj,ni+1]-entropy_data[nk,nj,ni])/(x1v[ni+1]-x1v[ni])
             elif ni==block_size[0]-1:
               dvrdr=(x1v[ni]**2*v1_data[nk,nj,ni]-x1v[ni-1]**2*v1_data[nk,nj,ni-1])/(x1v[ni]-x1v[ni-1])
               dvpdr=(x1v[ni]*v3_data[nk,nj,ni]-x1v[ni]*v3_data[nk,nj,ni-1])/(x1v[ni]-x1v[ni-1])
               dvtdr=(x1v[ni]*v2_data[nk,nj,ni]-x1v[ni]*v2_data[nk,nj,ni-1])/(x1v[ni]-x1v[ni-1])
               dsdr=(entropy_data[nk,nj,ni]-entropy_data[nk,nj,ni-1])/(x1v[ni]-x1v[ni-1])
             else:
               dvrdr=(x1v[ni+1]**2*v1_data[nk,nj,ni+1]-x1v[ni-1]**2*v1_data[nk,nj,ni-1])/(x1v[ni+1]-x1v[ni-1])
               dvpdr=(x1v[ni+1]*v3_data[nk,nj,ni+1]-x1v[ni-1]*v3_data[nk,nj,ni-1])/(x1v[ni+1]-x1v[ni-1])
               dvtdr=(x1v[ni+1]*v2_data[nk,nj,ni+1]-x1v[ni-1]*v2_data[nk,nj,ni-1])/(x1v[ni+1]-x1v[ni-1])
               dsdr=(entropy_data[nk,nj,ni+1]-entropy_data[nk,nj,ni-1])/(x1v[ni+1]-x1v[ni-1])
             
             dvpdr=dvpdr/x1v[ni]
             dvtdr=dvtdr/x1v[ni]
             
             dvrdr=dvrdr/(x1v[ni]**2)

             if nj==0:
               dvtdt=(np.sin(x2v[nj+1])*v2_data[nk,nj+1,ni]-np.sin(x2v[nj])*v2_data[nk,nj,ni])/(x2v[nj+1]-x2v[nj])
               dvpdt=(np.sin(x2v[nj+1])*v3_data[nk,nj+1,ni]-np.sin(x2v[nj])*v3_data[nk,nj,ni])/(x2v[nj+1]-x2v[nj])
               dvrdt=(v1_data[nk,nj+1,ni]-v1_data[nk,nj,ni])/(x2v[nj+1]-x2v[nj])
               dsdt=(entropy_data[nk,nj+1,ni]-entropy_data[nk,nj,ni])/(x2v[nj+1]-x2v[nj])
             elif nj==block_size[1]-1:
               dvtdt=(np.sin(x2v[nj])*v2_data[nk,nj,ni]-np.sin(x2v[nj-1])*v2_data[nk,nj-1,ni])/(x2v[nj]-x2v[nj-1])
               dvpdt=(np.sin(x2v[nj])*v3_data[nk,nj,ni]-np.sin(x2v[nj-1])*v3_data[nk,nj-1,ni])/(x2v[nj]-x2v[nj-1])
               dvrdt=(v1_data[nk,nj,ni]-v1_data[nk,nj-1,ni])/(x2v[nj]-x2v[nj-1])
               dsdt=(entropy_data[nk,nj,ni]-entropy_data[nk,nj-1,ni])/(x2v[nj]-x2v[nj-1])
             else:
               dvtdt=(np.sin(x2v[nj+1])*v2_data[nk,nj+1,ni]-np.sin(x2v[nj-1])*v2_data[nk,nj-1,ni])/(x2v[nj+1]-x2v[nj-1])
               dvpdt=(np.sin(x2v[nj+1])*v3_data[nk,nj+1,ni]-np.sin(x2v[nj-1])*v3_data[nk,nj-1,ni])/(x2v[nj+1]-x2v[nj-1])
               dvrdt=(v1_data[nk,nj+1,ni]-v1_data[nk,nj-1,ni])/(x2v[nj+1]-x2v[nj-1])
               dsdt=(entropy_data[nk,nj+1,ni]-entropy_data[nk,nj-1,ni])/(x2v[nj+1]-x2v[nj-1])

             dvtdt=dvtdt/(x1v[ni]*np.sin(x2v[nj]))
             
             dvpdt=dvpdt/(x1v[ni]*np.sin(x2v[nj]))
             dvrdt=dvrdt/x1v[ni]
             dsdt=dsdt/x1v[ni]

             if nk==0:
               dvpdp=(v3_data[nk+1,nj,ni]-v3_data[nk,nj,ni])/(x3v[nk+1]-x3v[nk])
               dvrdp=(v1_data[nk+1,nj,ni]-v1_data[nk,nj,ni])/(x3v[nk+1]-x3v[nk])
               dvtdp=(v2_data[nk+1,nj,ni]-v2_data[nk,nj,ni])/(x3v[nk+1]-x3v[nk])
               dsdp=(entropy_data[nk+1,nj,ni]-entropy_data[nk,nj,ni])/(x3v[nk+1]-x3v[nk])
             elif nk==block_size[2]-1:
               dvpdp=(v3_data[nk,nj,ni]-v3_data[nk-1,nj,ni])/(x3v[nk]-x3v[nk-1])
               dvrdp=(v1_data[nk,nj,ni]-v1_data[nk-1,nj,ni])/(x3v[nk]-x3v[nk-1])
               dvtdp=(v2_data[nk,nj,ni]-v2_data[nk-1,nj,ni])/(x3v[nk]-x3v[nk-1])
               dsdp=(entropy_data[nk,nj,ni]-entropy_data[nk-1,nj,ni])/(x3v[nk]-x3v[nk-1])
             else:
               dvpdp=(v3_data[nk+1,nj,ni]-v3_data[nk-1,nj,ni])/(x3v[nk+1]-x3v[nk-1])
               dvrdp=(v1_data[nk+1,nj,ni]-v1_data[nk-1,nj,ni])/(x3v[nk+1]-x3v[nk-1])
               dvtdp=(v2_data[nk+1,nj,ni]-v2_data[nk-1,nj,ni])/(x3v[nk+1]-x3v[nk-1])
               dsdp=(entropy_data[nk+1,nj,ni]-entropy_data[nk-1,nj,ni])/(x3v[nk+1]-x3v[nk-1])
               
             dvpdp=dvpdp/(x1v[ni]*np.sin(x2v[nj]))
             dvrdp=dvrdp/(x1v[ni]*np.sin(x2v[nj]))
             dvtdp=dvtdp/(x1v[ni]*np.sin(x2v[nj]))
             dsdp=dsdp/(x1v[ni])*np.sin(x2v[nj])

             divvel[nk,nj,ni]=dvrdr+dvtdt+dvpdp
             
             curvtrho[nk,nj,ni]=(dvrdp-dvpdr)/rho_data[nk,nj,ni]
             curvtrhosq[nk,nj,ni]=curvtrho[nk,nj,ni]*curvtrho[nk,nj,ni]
             curvgrads[nk,nj,ni]=((dvpdt-dvtdp)*dsdr+(dvrdp-dvpdr)*dsdt+(dvtdr-dvrdt)*dsdp)/rho_data[nk,nj,ni]
             curvgradssq[nk,nj,ni]=curvgrads[nk,nj,ni]*curvgrads[nk,nj,ni]
             
             # the shear stress



             dvphidr=0

             if ni==0:
               dvphidr=(v3_data[nk,nj,ni+1]/x1v[ni+1]-v3_data[nk,nj,ni]/x1v[ni])/(x1v[ni+1]-x1v[ni])
             elif ni==block_size[0]-1:
               dvphidr=(v3_data[nk,nj,ni]/x1v[ni]-v3_data[nk,nj,ni-1]/x1v[ni-1])/(x1v[ni]-x1v[ni-1])
             else:
               dvphidr=(v3_data[nk,nj,ni+1]/x1v[ni+1]-v3_data[nk,nj,ni-1]/x1v[ni-1])/(x1v[ni+1]-x1v[ni-1])

             dvphidr=dvphidr*x1v[ni]

             dvrdphi=0

             if nk==0:
               dvrdphi=(v1_data[nk+1,nj,ni]-v1_data[nk,nj,ni])/(x3v[nk+1]-x3v[nk])
             elif nk==block_size[2]-1:
               dvrdphi=(v1_data[nk,nj,ni]-v1_data[nk-1,nj,ni])/(x3v[nk]-x3v[nk-1])
             else:
               dvrdphi=(v1_data[nk+1,nj,ni]-v1_data[nk-1,nj,ni])/(x3v[nk+1]-x3v[nk-1])

             dvrdphi=dvrdphi/(x1v[ni]*np.sin(x2v[nj]))
             
             velshear[nk,nj,ni]=dvphidr+dvrdphi


      radviscoef=(8.0/27.0)* Er_data
      radviscoef=np.multiply(radviscoef,velshear)
      radviscoef=np.divide(radviscoef,rho_data)
      radviscoef=np.divide(radviscoef,kappatot)
             
             
             
      
      inflowindex=v1_data < 0.0
      outflowindex = v1_data > 0.0
      rho_inflow=np.copy(rho_data)
      rho_inflow[outflowindex] = 1.e-14
      rho_outflow=np.copy(rho_data)
      rho_outflow[inflowindex] = 1.e-14
      
      vphi_inflow = np.copy(v3_data)
      vphi_inflow[outflowindex] = 0.0
      vphi_outflow = np.copy(v3_data)
      vphi_outflow[inflowindex] = 0.0
      
      
      

      for jo in jo_vals:
        for io in io_vals:
            vol_phi[jl+jo:ju+jo:s,il+io:iu+io:s] \
                = vol_phi[jl+jo:ju+jo:s,il+io:iu+io:s] + block_phi_size


  # Assign values
      for q,dataset,index in zip(quantities,quantity_datasets,quantity_indices):
        if q=='rhosq':
           oridata=np.multiply(rho_data,rho_data)
        elif q=='PB':
           oridata=PB_data
        elif q=='PB1':
           oridata=PB1_data
        elif q=='PB2':
           oridata=PB2_data
        elif q=='PB3':
           oridata=PB3_data
        elif q=='PBsq':
           oridata=np.multiply(PB_data,PB_data)
        elif q=='kappa_s':
           oridata=np.divide(sigma_s_data,rho_data)
        elif q=='kappa_a':
           oridata=np.divide(sigma_a_data,rho_data)
        elif q=='Maxwell':
           bxdata=np.multiply(B2_data,np.cos(grid_theta))+np.multiply(B1_data,np.sin(grid_theta))
           oridata=-np.multiply(bxdata,B3_data)
        elif q=='Reynolds':
           vxdata=np.multiply(v2_data,np.cos(grid_theta))+np.multiply(v1_data,np.sin(grid_theta))
           oridata=np.multiply(vxdata,np.multiply(v3_data,rho_data))
        elif q=='Radacc':
           oridata=np.divide(np.multiply(Fr01_data,(sigma_s_data+sigma_a_data)),rho_data)
        elif q=='RhoVr':
           oridata=np.multiply(v1_data,rho_data)
        elif q=='RhoVtheta':
           oridata=np.multiply(v2_data,rho_data)
        elif q=='RhoVphi':
           oridata=np.multiply(v3_data,rho_data)
        elif q=='RhoVout':
           oridata=np.multiply(vout_data,rho_data)
        elif q=='RhoVin':
           oridata=np.multiply(vin_data,rho_data)
        elif q=='Ekin1':
           oridata=0.5*np.multiply(v1_data, v1_data)
           oridata=np.multiply(oridata,rho_data)
        elif q=='Ekin2':
           oridata=0.5*np.multiply(v2_data, v2_data)
           oridata=np.multiply(oridata,rho_data)
        elif q=='Ekin3':
           oridata=0.5*np.multiply(v3_data, v3_data)
           oridata=np.multiply(oridata,rho_data)
        elif q=='BrBphi':
           oridata=np.multiply(B1_data,B3_data)
        elif q=='BthetaBphi':
           oridata=np.multiply(B2_data,B3_data)
        elif q=='rhovrvphi':
           oridata=rhovrvphi_data
        elif q=='rhovthetavphi':
           oridata=np.multiply(v2_data,rhovphi_data)
        elif q=='vrEr':
           oridata=np.multiply(v1_data,Er_data)
        elif q=='vthetaEr':
           oridata=np.multiply(v2_data,Er_data)
        elif q=='rhoPB':
           oridata=np.multiply(rho_data,PB_data)
        elif q=='tgas':
           oridata=np.divide(pgas_data,rho_data)
        elif q=='Fr01Sigma':
           oridata=np.multiply(Fr01_data,(sigma_s_data+sigma_a_data))
        elif q=='Fr02Sigma':
           oridata=np.multiply(Fr02_data,(sigma_s_data+sigma_a_data))
        elif q=='Fr01kappa':
           oridata=np.divide(np.multiply(Fr01_data,(sigma_s_data+sigma_a_data)),rho_data)
        elif q=='Fr02kappa':
           oridata=np.divide(np.multiply(Fr02_data,(sigma_s_data+sigma_a_data)),rho_data)
        elif q=='rhovrsq':
           rhovr=np.multiply(v1_data,rho_data)
           oridata=np.multiply(rhovr,rhovr)
        elif q=='vrsq':
           oridata=np.multiply(v1_data,v1_data)
        elif q=='rhovphiout':
           oridata=np.multiply(vphi_outflow,rho_outflow)
        elif q=='rhovphiin':
           oridata=np.multiply(vphi_inflow,rho_inflow)
        elif q=='rhoinflow':
           oridata=rho_inflow
        elif q=='rhooutflow':
           oridata=rho_outflow
        elif q=='rhovrvphisq':
           oridata=np.multiply(rhovrvphi_data,rhovrvphi_data)
        elif q=='Ersq':
           oridata=np.multiply(Er_data,Er_data)
        elif q=='rhovrvphiEr':
           oridata=np.multiply(rhovrvphi_data,Er_data)
        elif q=='rhovrvphiPB':
           oridata=np.multiply(rhovrvphi_data,PB_data)
        elif q=='radvis':
           oridata=radviscoef
        elif q=='velshear':
           oridata=velshear
        elif q=='divvel':
           oridata=divvel
        elif q=='raddivwork':
           oridata=np.multiply(Er_data,divvel)
        elif q=='radshearwork':
           oridata=np.multiply(radviscoef,0.5*velshear)
        elif q=='curvtrho':
           oridata=curvtrho
        elif q=='curvtrhosq':
           oridata=curvtrhosq
        elif q=='curvgrads':
           oridata=curvgrads
        elif q=='curvgradssq':
           oridata=curvgradssq
        elif q=='Fr1new':
           oridata=(lorentz**2)*(Fr01_data+v1_data*4.0*Er0_data/(crat*3.0))
        elif q=='Fr2new':
           oridata=(lorentz**2)*(Fr02_data+v2_data*4.0*Er0_data/(crat*3.0))
        elif q=='Fr3new':
           oridata=(lorentz**2)*(Fr03_data+v3_data*4.0*Er0_data/(crat*3.0))
        elif q=='lambda_a':
           oridata=(PB2_data*2.0/rho_data)**0.5
        else:
           oridata=f[dataset][index,block_num,:]

        
        for jo in jo_vals:
          for io in io_vals:
            data[q][jl+jo:ju+jo:s,il+io:iu+io:s] \
                     = data[q][jl+jo:ju+jo:s,il+io:iu+io:s] + \
                       block_phi_size * np.mean(oridata,axis=0)
        


# divide by the total phi volume

  for q in quantities:
      data[q][:,:] = data[q][:,:]/vol_phi[:,:]
      data[q][:,:] = np.nan_to_num(data[q][:,:])


# get the cell center coordinates
  # get cell center coordinates
  data['x1v']=np.zeros(nx1)
  data['x2v']=np.zeros(nx2)
  
  for ni in range(nx1):
    data['x1v'][ni]=0.75*(data['x1f'][ni+1]**4.0 - data['x1f'][ni]**4.0)/(data['x1f'][ni+1]**3.0-data['x1f'][ni]**3.0);
  
  for nj in range(nx2):
    data['x2v'][nj]=((np.sin(data['x2f'][nj+1]) - data['x2f'][nj+1] * np.cos(data['x2f'][nj+1])) \
           -(np.sin(data['x2f'][nj]) - data['x2f'][nj] * np.cos(data['x2f'][nj]))) \
              / (np.cos(data['x2f'][nj]) - np.cos(data['x2f'][nj+1]));

  #add time
  data['Time']=f.attrs['Time']

  f.close()

  return data


