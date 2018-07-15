# Test script for initializing problem with preexisting array

# Modules
import numpy as np
import sys
import h5py
import scripts.utils.athena as athena
import scripts.utils.comparison as comparison
sys.path.insert(0, '../../vis/python')
import athena_read # noqa

# Parameters
filename_input = 'initial_data.hdf5'
filename_output = 'from_array.cons.00000.athdf'
dataset_name = 'cons'
num_blocks = 2
nx1 = 8
nx2 = 6
nx3 = 4
gamma = 5.0/3.0

# Prepare Athena++
def prepare(**kwargs):

    # Configure and compile code
    athena.configure('hdf5',
                     prob='from_array',
                     **kwargs)
    athena.make()

    # Write file to be loaded
    cons_density = np.reshape(np.arange(384), (1, num_blocks, nx3, nx2, nx1))
    cons_momentum = np.zeros((3, num_blocks, nx3, nx2, nx1))
    cons_energy = np.ones((1, num_blocks, nx3, nx2, nx1)) / (gamma - 1.0)
    cons_input = np.vstack((cons_density, cons_momentum, cons_energy))
    with h5py.File('bin/{0}'.format(filename_input), 'w') as f:
        f.create_dataset(dataset_name, data=cons_input)


# Run Athena++
def run(**kwargs):
    arguments = ['time/tlim=0',
                 'time/ncycle_out=0',
                 'problem/input_filename={0}'.format(filename_input)]
    athena.run('mhd/athinput.from_array', arguments)


# Analyze outputs
def analyze():
    with h5py.File('bin/{0}'.format(filename_input), 'r') as f:
        cons_input = f[dataset_name][:]
    with h5py.File('bin/{0}'.format(filename_output), 'r') as f:
        cons_output = f[dataset_name][:]
    print(cons_input.shape)
    print(cons_input[0,0,1,1,1])
    print(cons_output.shape)
    print(cons_output[0,0,1,1,1])
    for n in range(5):
      for m in range(2):
        for k in range(4):
          for j in range(6):
            for i in range(8):
              if cons_input[n,m,k,j,i] != cons_output[n,m,k,j,i]:
                print('{0},{1},{2},{3},{4}: {5}, {6}'.format(n, m, k, j, i, cons_input[n,m,k,j,i], cons_output[n,m,k,j,i]))
    if np.all(cons_input == cons_output):
        return True
    return False
