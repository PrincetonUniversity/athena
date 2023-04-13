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




  quantities.append('PB1')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('PB2')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('PB3')
  quantity_datasets.append('None')
  quantity_indices.append(0)


  quantities.append('rhovr')
  quantity_datasets.append('None')
  quantity_indices.append(0)

  quantities.append('lambda_a')
  quantity_datasets.append('None')
  quantity_indices.append(0)



  n_rmax=nx1
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
    vr_data=f[spec_datasets[2]][spec_indices[2],block_num,:]
    B1_data=f[spec_datasets[5]][spec_indices[5],block_num,:]
    B2_data=f[spec_datasets[6]][spec_indices[6],block_num,:]
    B3_data=f[spec_datasets[7]][spec_indices[7],block_num,:]

    PB_data=0.5*(np.multiply(B1_data,B1_data)+np.multiply(B2_data,B2_data)+np.multiply(B3_data,B3_data))
    PB1_data=0.5*np.multiply(B1_data,B1_data)   
    PB2_data=0.5*np.multiply(B2_data,B2_data)
    PB3_data=0.5*np.multiply(B3_data,B3_data)



    
    

  

    for jo in jo_vals:
      for io in io_vals:
          vol_phi[jl+jo:ju+jo:s,il+io:iu+io:s] \
              = vol_phi[jl+jo:ju+jo:s,il+io:iu+io:s] + block_phi_size


  # Assign values
    for q,dataset,index in zip(quantities,quantity_datasets,quantity_indices):
      if q == 'PB1':
         oridata = PB1_data
      elif q == 'PB2':
         oridata = PB2_data
      elif q == 'PB3':
         oridata = PB3_data
      elif q == 'rhovr':
         oridata = np.multiply(rho_data,vr_data)
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
  data['quantities']=quantities
  f.close()

  return data


