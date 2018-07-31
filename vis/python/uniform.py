"""
Read .athdf data files and write new ones as single block at constant refinement
level.

Note: Requires h5py.

Note: Only works for 3D data.
"""

# Python standard modules
import argparse
import os

# Other Python modules
import h5py

# Athena++ modules
import athena_read

# Main function


def main(**kwargs):

    # Determine which files to process given possible MPI information
    file_nums = range(kwargs['start'], kwargs['end']+1, kwargs['stride'])
    if kwargs['m']:
        size = int(os.environ['OMPI_COMM_WORLD_SIZE'])
        rank = int(os.environ['OMPI_COMM_WORLD_RANK'])
        num_files = len(file_nums)
        num_files_per_rank = num_files/size
        num_files_extra = num_files % size
        num_files_list = ([num_files_per_rank+1] * num_files_extra
                          + [num_files_per_rank] * (size-num_files_extra))
        num_files_previous = sum(num_files_list[:rank])
        num_files_current = num_files_list[rank]
        file_nums_local = file_nums[num_files_previous:num_files_previous
                                    + num_files_current]
    else:
        file_nums_local = file_nums

    # Go through list of files
    for n in file_nums_local:

        # Determine filenames
        input_filename = '{0}.out{1}.{2:05d}.athdf'\
            .format(kwargs['input_filename'], kwargs['output_num'], n)
        output_filename = '{0}.out{1}.{2:05d}.athdf'\
            .format(kwargs['output_filename'], kwargs['output_num'], n)
        output_dir, output_base = os.path.split(output_filename)

        # Read attributes and data
        with h5py.File(input_filename, 'r') as f:
            attributes = f.attrs.items()
            attrs = dict(attributes)
            level = f.attrs['MaxLevel']
        subsample = False
        if kwargs['level'] is not None:
            if level > kwargs['level']:
                subsample = True
            level = kwargs['level']
        data = athena_read.athdf(input_filename, quantities=kwargs['quantities'],
                                 level=level, subsample=subsample)

        # Determine new grid size
        nx1 = attrs['RootGridSize'][0] * 2**level
        nx2 = attrs['RootGridSize'][1] * 2**level
        nx3 = attrs['RootGridSize'][2] * 2**level

        # Create new HDF5 file
        with h5py.File(output_filename, 'w') as f:

            # Write attributes
            for key, val in attributes:
                if key == 'RootGridX1' or key == 'RootGridX2' or key == 'RootGridX3':
                    if val[2] > 0.0:
                        value = [val[0], val[1], val[2]**(1.0/2.0**level)]
                    else:
                        value = [val[0], val[1], val[2]]
                elif key == 'RootGridSize':
                    value = [nx1, nx2, nx3]
                elif key == 'NumMeshBlocks':
                    value = 1
                elif key == 'MeshBlockSize':
                    value = [nx1, nx2, nx3]
                elif key == 'MaxLevel':
                    value = 0
                elif key == 'NumVariables' and kwargs['quantities'] is not None:
                    value = [len(kwargs['quantities'])]
                elif key == 'DatasetNames' and kwargs['quantities'] is not None:
                    value = ['quantities']
                elif key == 'VariableNames' and kwargs['quantities'] is not None:
                    value = kwargs['quantities']
                else:
                    value = val
                f.attrs.create(key, value, dtype=val.dtype)

            # Write datasets
            f.create_dataset('Levels', data=[0], dtype='>i4')
            f.create_dataset(
                'LogicalLocations', data=[
                    0, 0, 0], dtype='>i8', shape=(
                    1, 3))
            f.create_dataset('x1f', data=data['x1f'], dtype='>f4', shape=(1, nx1 + 1))
            f.create_dataset('x2f', data=data['x2f'], dtype='>f4', shape=(1, nx2 + 1))
            f.create_dataset('x3f', data=data['x3f'], dtype='>f4', shape=(1, nx3 + 1))
            f.create_dataset('x1v', data=data['x1v'], dtype='>f4', shape=(1, nx1))
            f.create_dataset('x2v', data=data['x2v'], dtype='>f4', shape=(1, nx2))
            f.create_dataset('x3v', data=data['x3v'], dtype='>f4', shape=(1, nx3))
            var_offset = 0
            for dataset_name, num_vars in zip(
                    f.attrs['DatasetNames'], f.attrs['NumVariables']):
                f.create_dataset(
                    dataset_name, dtype='>f4', shape=(
                        num_vars, 1, nx3, nx2, nx1))
                for var_num in range(num_vars):
                    variable_name = f.attrs['VariableNames'][var_num + var_offset]
                    f[dataset_name][var_num, 0, :, :, :] = data[variable_name]
                var_offset += num_vars

        # Create new XDMF file
        if kwargs['x']:
            with open(output_filename + '.xdmf', 'w') as f:
                f.write('<?xml version="1.0" ?>\n')
                f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
                f.write('<Xdmf Version="2.0">\n')
                f.write('<Domain>\n')
                f.write('<Grid Name="Mesh" GridType="Collection">\n')
                f.write('  <Grid Name="MeshBlock0" GridType="Uniform">\n')
                f.write(('    <Topology TopologyType="3DRectMesh"' +
                         ' NumberOfElements="{0} {1} {2}"/>\n').format(nx3+1, nx2+1,
                                                                       nx1+1))
                f.write('    <Geometry GeometryType="VXVYVZ">\n')
                for nx, xf_string in zip((nx1, nx2, nx3), ('x1f', 'x2f', 'x3f')):
                    f.write(
                        '      <DataItem ItemType="HyperSlab" Dimensions="{0}">\n'.format(nx + 1)) # noqa
                    f.write(('        <DataItem Dimensions="3 2" NumberType="Int">'
                             + ' 0 0 1 1 1 {0} </DataItem>\n').format(nx + 1))
                    f.write(
                        ('        <DataItem Dimensions="1 {0}" Format="HDF">'
                         + ' {1}:/{2} </DataItem>\n').format(nx + 1, output_base, xf_string)) # noqa
                    f.write('      </DataItem>\n')
                f.write('    </Geometry>\n')
                if kwargs['quantities'] is None:
                    num_variables = attrs['NumVariables']
                    dataset_names = attrs['DatasetNames']
                    variable_names = attrs['VariableNames']
                else:
                    num_variables = [len(kwargs['quantities'])]
                    dataset_names = ['quantities']
                    variable_names = kwargs['quantities']
                var_offset = 0
                for dataset_name, num_vars in zip(dataset_names, num_variables):
                    for var_num in range(num_vars):
                        variable_name = variable_names[var_num + var_offset]
                        f.write(
                            '    <Attribute Name="{0}" Center="Cell">\n'.format(variable_name)) # noqa
                        f.write(
                            '      <DataItem ItemType="HyperSlab" Dimensions="{0} {1} {2}">\n' .format(nx3, nx2, nx1)) # noqa
                        f.write(('        <DataItem Dimensions="3 5" NumberType="Int">'
                                 + ' {0} 0 0 0 0 1 1 1 1 1 1 1 {1} {2} {3} </DataItem>\n')
                                .format(var_num, nx3, nx2, nx1))
                        f.write(
                            ('        <DataItem Dimensions="{0} 1 {1} {2} {3}" Format="HDF">' # noqa
                             + ' {4}:/{5} </DataItem>\n') .format(num_vars, nx3, nx2, nx1, output_base, dataset_name)) # noqa
                        f.write('      </DataItem>\n')
                        f.write('    </Attribute>\n')
                    var_offset += num_vars
                f.write('  </Grid>\n')
                f.write('</Grid>\n')
                f.write('</Domain>\n')
                f.write('</Xdmf>')


# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input_filename',
                        type=str,
                        help='base name of files to be converted, including directory')
    parser.add_argument('output_filename',
                        type=str,
                        help='base name of new files to be saved, including directory')
    parser.add_argument('output_num',
                        type=int,
                        help='number of output to convert')
    parser.add_argument('start',
                        type=int,
                        help='first file number to be converted')
    parser.add_argument('end',
                        type=int,
                        help='last file number to be converted')
    parser.add_argument('stride',
                        type=int,
                        default=0,
                        help='stride in file numbers to be converted')
    parser.add_argument('-m',
                        action='store_true',
                        help='flag indicating this is an MPI job')
    parser.add_argument('-x',
                        action='store_false',
                        help='flag indicating no XDMF file should be written')
    parser.add_argument('-l', '--level',
                        type=int,
                        help='refinement level to use, overriding max level present')
    parser.add_argument('-q', '--quantities',
                        type=str,
                        nargs='+',
                        help='names of quantities to extract')
    args = parser.parse_args()
    main(**vars(args))
