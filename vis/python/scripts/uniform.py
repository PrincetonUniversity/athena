"""
Read .athdf data files and write new ones as single block at constant refinement level.

Note: Requires h5py.

Note: Only works for 3D data.
"""

# Python modules
import h5py
import numpy as np

from mpi4py import MPI



# Athena++ modules
import athena_read

# Main function
def main(n,input_base,output_base,quantities=None,leveln=None,xdmf=1):

  # Go through list of files
    # Determine filenames
    input_filename = input_base+'{:05d}'.format(n)+'.athdf'
    output_filename = output_base+'{:05d}'.format(n)+'.athdf'

    # Read attributes and data
    with h5py.File(input_filename, 'r') as f:
      attributes = f.attrs.items()
      attrs = dict(attributes)
      level = f.attrs['MaxLevel']
      time = f.attrs['Time']
    subsample = False
    if leveln is not None:
      if level > leveln:
        subsample = True
      level = leveln
    data = athena_read.athdf(input_filename, quantities=quantities,
        level=level, subsample=subsample)

    # Determine new grid size
    nx1 = attrs['RootGridSize'][0] * 2**level
    nx2 = attrs['RootGridSize'][1] * 2**level
    nx3 = attrs['RootGridSize'][2] * 2**level

    # Create new HDF5 file
    with h5py.File(output_filename, 'w') as f:

      # Write attributes
      for key,val in attributes:
        if key == 'RootGridX1' or key == 'RootGridX2' or key == 'RootGridX3':
          value = [val[0], val[1], val[2]**(1.0/2.0**level)]
        elif key == 'RootGridSize':
          value = [nx1, nx2, nx3]
        elif key == 'NumMeshBlocks':
          value = 1
        elif key == 'MeshBlockSize':
          value = [nx1, nx2, nx3]
        elif key == 'MaxLevel':
          value = 0
        elif key == 'NumVariables' and quantities is not None:
          value = [len(quantities)]
        elif key == 'DatasetNames' and quantities is not None:
          value = ['quantities']
        elif key == 'VariableNames' and quantities is not None:
          value = quantities
        else:
          value = val
        f.attrs.create(key, value, dtype=val.dtype)

      # Write datasets
      f.create_dataset('Levels', data=[0], dtype='>i4')
      f.create_dataset('LogicalLocations', data=[0,0,0], dtype='>i8', shape=(1,3))
      f.create_dataset('x1f', data=data['x1f'], dtype='>f4', shape=(1,nx1+1))
      f.create_dataset('x2f', data=data['x2f'], dtype='>f4', shape=(1,nx2+1))
      f.create_dataset('x3f', data=data['x3f'], dtype='>f4', shape=(1,nx3+1))
      var_offset = 0
      for dataset_name,num_vars in zip(f.attrs['DatasetNames'],f.attrs['NumVariables']):
        f.create_dataset(dataset_name, dtype='>f4', shape=(num_vars,1,nx3,nx2,nx1))
        for var_num in range(num_vars):
          variable_name = f.attrs['VariableNames'][var_num+var_offset]
          f[dataset_name][var_num,0,:,:,:] = data[variable_name]
        var_offset += num_vars

    # Create new XDMF file
    if xdmf > 0:
      with open(output_filename+'.xdmf', 'w') as f:
        f.write('<?xml version="1.0" ?>\n')
        f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
        f.write('<Xdmf Version="2.0">\n')
        f.write('<Domain>\n')
        f.write('<Grid Name="Mesh" GridType="Collection">\n')
        f.write('  <Grid Name="MeshBlock0" GridType="Uniform">\n')
        f.write('     <Time Type="Single" Value="{0}"/>\n'.format(time))
        f.write(('    <Topology TopologyType="3DRectMesh"' \
            + ' NumberOfElements="{0} {1} {2}"/>\n').format(nx3+1,nx2+1,nx1+1))
        f.write('    <Geometry GeometryType="VXVYVZ">\n')
        for nx,xf_string in zip((nx1,nx2,nx3),('x1f','x2f','x3f')):
          f.write('      <DataItem ItemType="HyperSlab" Dimensions="{0}">\n'\
              .format(nx+1))
          f.write(('        <DataItem Dimensions="3 2" NumberType="Int">' \
              + ' 0 0 1 1 1 {0} </DataItem>\n').format(nx+1))
          f.write(('        <DataItem Dimensions="1 {0}" Format="HDF">' \
              + ' {1}:/{2} </DataItem>\n').format(nx+1,output_filename,xf_string))
          f.write('      </DataItem>\n')
        f.write('    </Geometry>\n')
        if quantities is None:
          num_variables = attrs['NumVariables']
          dataset_names = attrs['DatasetNames']
          variable_names = attrs['VariableNames']
        else:
          num_variables = [len(quantities)]
          dataset_names = ['quantities']
          variable_names = quantities
        var_offset = 0
        for dataset_name,num_vars in zip(dataset_names,num_variables):
          for var_num in range(num_vars):
            variable_name = variable_names[var_num+var_offset]
            f.write('    <Attribute Name="{0}" Center="Cell">\n'.format(variable_name))
            f.write('      <DataItem ItemType="HyperSlab" Dimensions="{0} {1} {2}">\n'\
                .format(nx3,nx2,nx1))
            f.write(('        <DataItem Dimensions="3 5" NumberType="Int">' \
                + ' {0} 0 0 0 0 1 1 1 1 1 1 1 {1} {2} {3} </DataItem>\n')\
                .format(var_num,nx3,nx2,nx1))
            f.write(('        <DataItem Dimensions="{0} 1 {1} {2} {3}" Format="HDF">' \
                + ' {4}:/{5} </DataItem>\n')\
                .format(num_vars,nx3,nx2,nx1,output_filename,dataset_name))
            f.write('      </DataItem>\n')
            f.write('    </Attribute>\n')
          var_offset += num_vars
        f.write('  </Grid>\n')
        f.write('</Grid>\n')
        f.write('</Domain>\n')
        f.write('</Xdmf>')

quantities=['rho','Er']

ni=7270
no=8745

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

for i in range(ni,no+1,nprocs*10):
  file_num=i+rank*10
  print file_num, rank
  main(file_num,'disk.out1.','vis.',quantities)


