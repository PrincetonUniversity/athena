# Serial test script for initializing problem with preexisting array

# Standard modules
import sys

# Other modules
import logging
import numpy as np
import h5py

# Athena modules
import scripts.utils.athena as athena
sys.path.insert(0, '../../vis/python')
import athena_read  # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module

# Parameters
filename_input = 'initial_data.hdf5'
filename_output = 'from_array.cons.00000.athdf'
dataset_cons = 'cons'
dataset_b1 = 'b1'
dataset_b2 = 'b2'
dataset_b3 = 'b3'
nb1 = 4
nx1 = 4
nx2 = 6
nx3 = 4
gamma = 5.0/3.0


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)

    # Configure and compile code
    athena.configure('b',
                     'hdf5', 'h5double',
                     prob='from_array',
                     **kwargs)
    athena.make()

    # Calculate initial field values
    b1 = np.empty((nx3, nx2, nb1 * nx1 + 1))
    b1[...] = np.arange(nx2)[None, :, None] - np.arange(nx3)[:, None, None]
    b1_input = np.empty((nb1, nx3, nx2, nx1 + 1))
    b2_input = np.zeros((nb1, nx3, nx2 + 1, nx1))
    b3_input = np.zeros((nb1, nx3 + 1, nx2, nx1))
    for n in range(nb1):
        b1_input[n, ...] = b1[:, :, n*nx1:(n+1)*nx1+1]
    # (second-order accurate assumption)
    b1v = 0.5 * (b1_input[:, :, :, :-1] + b1_input[:, :, :, 1:])

    # Calculate initial conserved values
    num_cells = nb1 * nx1 * nx2 * nx3
    density = np.reshape(np.arange(1, num_cells+1), (1, nb1, nx3, nx2, nx1))
    momentum = np.zeros((3, nb1, nx3, nx2, nx1))
    energy = np.ones((1, nb1, nx3, nx2, nx1)) / (gamma - 1.0) + 0.5 * b1v[None, ...] ** 2
    cons_input = np.vstack((density, momentum, energy))

    # Write file to be loaded
    with h5py.File('bin/{0}'.format(filename_input), 'w') as f:
        f.create_dataset(dataset_cons, data=cons_input)
        f.create_dataset(dataset_b1, data=b1_input)
        f.create_dataset(dataset_b2, data=b2_input)
        f.create_dataset(dataset_b3, data=b3_input)


# Run Athena++
def run(**kwargs):
    arguments = ['time/tlim=0',
                 'time/ncycle_out=0',
                 'mesh/nx1={0}'.format(nb1 * nx1),
                 'mesh/nx2={0}'.format(nx2),
                 'mesh/nx3={0}'.format(nx3),
                 'meshblock/nx1={0}'.format(nx1),
                 'meshblock/nx2={0}'.format(nx2),
                 'meshblock/nx3={0}'.format(nx3),
                 'problem/input_filename={0}'.format(filename_input)]
    athena.run('mhd/athinput.from_array', arguments)


# Analyze outputs
def analyze():
    analyze_status = True
    # Read input data
    with h5py.File('bin/{0}'.format(filename_input), 'r') as f:
        cons_input = f[dataset_cons][:]
        b1_input = f[dataset_b1][:]
        b2_input = f[dataset_b2][:]
        b3_input = f[dataset_b3][:]

    # Calculate cell-centered field inputs from face-centered values
    # (second-order accurate assumption)
    b1v = 0.5 * (b1_input[:, :, :, :-1] + b1_input[:, :, :, 1:])
    b2v = 0.5 * (b2_input[:, :, :-1, :] + b2_input[:, :, 1:, :])
    b3v = 0.5 * (b3_input[:, :-1, :, :] + b3_input[:, 1:, :, :])

    # Read output data
    with h5py.File('bin/{0}'.format(filename_output), 'r') as f:
        num_vars = f.attrs['NumVariables']
        dataset_names = f.attrs['DatasetNames'].astype('U')
        output_vars = f.attrs['VariableNames'].astype('U')
        cons_output = f['cons'][:]
        field_output = f['B'][:]

    # Order conserved output data to match inputs
    index_cons = np.where(dataset_names == 'cons')[0][0]
    num_vars_cons = num_vars[index_cons]
    num_vars_pre_cons = np.sum(num_vars[:index_cons])
    output_vars_cons = output_vars[num_vars_pre_cons:num_vars_pre_cons+num_vars_cons]
    dens = cons_output[np.where(output_vars_cons == 'dens')[0], ...]
    mom1 = cons_output[np.where(output_vars_cons == 'mom1')[0], ...]
    mom2 = cons_output[np.where(output_vars_cons == 'mom2')[0], ...]
    mom3 = cons_output[np.where(output_vars_cons == 'mom3')[0], ...]
    etot = cons_output[np.where(output_vars_cons == 'Etot')[0], ...]
    cons_output = np.vstack((dens, mom1, mom2, mom3, etot))

    # Order field output data to match inputs
    index_field = np.where(dataset_names == 'B')[0][0]
    num_vars_field = num_vars[index_field]
    num_vars_pre_field = np.sum(num_vars[:index_field])
    output_vars_field = output_vars[num_vars_pre_field:num_vars_pre_field+num_vars_field]
    b1_output = field_output[np.where(output_vars_field == 'Bcc1')[0][0], ...]
    b2_output = field_output[np.where(output_vars_field == 'Bcc2')[0][0], ...]
    b3_output = field_output[np.where(output_vars_field == 'Bcc3')[0][0], ...]

    # Check that outputs match inputs
    if not np.all(cons_output == cons_input):
        analyze_status = False
    if not np.all(b1_output == b1v):
        analyze_status = False
    if not np.all(b2_output == b2v):
        analyze_status = False
    if not np.all(b3_output == b3v):
        analyze_status = False
    return analyze_status
