# Test script for initializing problem with preexisting array

# Modules
import numpy as np
import sys
import h5py
import scripts.utils.athena as athena
sys.path.insert(0, '../../vis/python')
import athena_read  # noqa

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
    num_cells = num_blocks * nx1 * nx2 * nx3
    cons_density = np.reshape(np.arange(1, num_cells+1), (1, num_blocks, nx3, nx2, nx1))
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

    # Read input and output data
    with h5py.File('bin/{0}'.format(filename_input), 'r') as f:
        cons_input = f[dataset_name][:]
    with h5py.File('bin/{0}'.format(filename_output), 'r') as f:
        output_variables = f.attrs['VariableNames']
        cons_output = f[dataset_name][:]

    # Order output data to match inputs
    dens = cons_output[np.where(output_variables == 'dens')[0], ...]
    mom1 = cons_output[np.where(output_variables == 'mom1')[0], ...]
    mom2 = cons_output[np.where(output_variables == 'mom2')[0], ...]
    mom3 = cons_output[np.where(output_variables == 'mom3')[0], ...]
    etot = cons_output[np.where(output_variables == 'Etot')[0], ...]
    cons_output = np.vstack((dens, mom1, mom2, mom3, etot))

    # Check that outputs match inputs
    np.set_printoptions(precision=17, floatmode='unique')
    if np.allclose(cons_output, cons_input, rtol=1.0e-15, atol=1.0e-15):
        return True
    return False
