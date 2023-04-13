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



def ReduceData(filename):

# Read attributes and data
  f=h5py.File(filename, 'r')
  level = f.attrs['MaxLevel']
  block_size = f.attrs['MeshBlockSize']
  root_grid_size = f.attrs['RootGridSize']
  nx1 = root_grid_size[0] * 2**level
  nx2 = root_grid_size[1] * 2**level if root_grid_size[1] > 1 else 1
  nx3 = root_grid_size[2] * 2**level if root_grid_size[2] > 1 else 1

  if nx3 > 1:
    dim = 3
  elif nx2 > 1:
    dim = 2
  else:
    dim = 1


  quantities = f.attrs['VariableNames'][:]
  quantities = [str(q) for q in quantities \
        if q != 'x1f' and q != 'x2f' and q != 'x3f']

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

# Get metadata describing file layout
  num_blocks = f.attrs['NumMeshBlocks']
  dataset_names = f.attrs['DatasetNames'][:]
  dataset_sizes = f.attrs['NumVariables'][:]
  dataset_sizes_cumulative = np.cumsum(dataset_sizes)
  variable_names = f.attrs['VariableNames'][:]
  levels = f['Levels'][:]
  logical_locations = f['LogicalLocations'][:]
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
                q=='Bcc2' or q=='Bcc3':
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


  quantities.append('rhovrvphiPB')
  quantity_datasets.append('None')
  quantity_indices.append(0)



  quantities.append('lambda_a')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('lz')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('lx')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('ly')
  quantity_datasets.append('None')
  quantity_indices.append(0)  

# the quantities
# get the azimuthal averaged data

  for q in quantities:
    data[q] = np.zeros((nx2,nx1))


# the total volume summed along phi direction
  vol_phi=np.zeros((nx2,nx1))


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
    PB_data=0.5*(np.multiply(B1_data,B1_data)+np.multiply(B2_data,B2_data)+np.multiply(B3_data,B3_data))
    PB1_data=0.5*np.multiply(B1_data,B1_data)   
    PB2_data=0.5*np.multiply(B2_data,B2_data)
    PB3_data=0.5*np.multiply(B3_data,B3_data)
    rhovphi_data=np.multiply(rho_data,v3_data)
    rhovrvphi_data=np.multiply(v1_data,rhovphi_data)


    velmag=np.sqrt(v1_data**2+v2_data**2+v3_data**2)

    
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

    r_3D=np.zeros((block_size[2],block_size[1],block_size[0]))
    t_3D=np.zeros((block_size[2],block_size[1],block_size[0]))    
    p_3D=np.zeros((block_size[2],block_size[1],block_size[0]))


    for k in range(block_size[2]):
      for j in range(block_size[1]):
        r_3D[k,j,:]=x1v

    for k in range(block_size[2]):
      for i in range(block_size[0]):
        t_3D[k,:,i]=x2v

    for j in range(block_size[1]):
      for i in range(block_size[0]):
        p_3D[:,j,i]=x3v


    lz = rho_data * r_3D * v3_data * np.sin(t_3D)
    lx = -rho_data * r_3D * (v3_data * np.cos(t_3D) * np.cos(p_3D) + v2_data * np.sin(p_3D))
    ly = -rho_data * r_3D * (v3_data * np.cos(t_3D) * np.sin(p_3D) + v2_data * np.cos(p_3D))
    
    
    

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
      elif q=='Maxwell':
         bxdata=np.multiply(B2_data,np.cos(grid_theta))+np.multiply(B1_data,np.sin(grid_theta))
         oridata=-np.multiply(bxdata,B3_data)
      elif q=='Reynolds':
         vxdata=np.multiply(v2_data,np.cos(grid_theta))+np.multiply(v1_data,np.sin(grid_theta))
         oridata=np.multiply(vxdata,np.multiply(v3_data,rho_data))
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
      elif q=='rhoPB':
         oridata=np.multiply(rho_data,PB_data)
      elif q=='tgas':
         oridata=np.divide(pgas_data,rho_data)
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
      elif q=='rhovrvphiPB':
         oridata=np.multiply(rhovrvphi_data,PB_data)
      elif q=='lambda_a':
         oridata=(PB2_data*2.0/rho_data)**0.5
      elif q=='lz':
         oridata=lz
      elif q=='lx':
         oridata=lx
      elif q=='ly':
         oridata=ly
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


