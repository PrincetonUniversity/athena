# Regression test for all output types
#
# Runs Orszag Tang vortex test, restarting the job twice and making history (hst),
# formatted table (.tab), VTK, and HDF5 (if available) outputs.  Then reads last
# version of each file to make sure output data is correct

# Modules
import logging
import numpy as np
import sys
import scripts.utils.athena as athena
sys.path.insert(0, '../../vis/python')
import athena_read  # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure('b', 'hdf5',
                     prob='orszag_tang',
                     flux='hlld', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    arguments = ['time/ncycle_out=0', 'time/nlim=80']
    athena.run('mhd/athinput.test_outputs', arguments)
    arguments = ['time/ncycle_out=0', 'time/nlim=330']
    athena.restart('TestOutputs.00001.rst', arguments)
    arguments = ['time/ncycle_out=0', 'time/nlim=-1']
    athena.restart('TestOutputs.00004.rst', arguments)


# Analyze outputs
def analyze():
    analyze_status = True
    # check density max and Vz and Bz components in tab slice
    slice_data = athena_read.tab(
        filename='bin/TestOutputs.block0.out2.00010.tab',
        raw=True,
        dimensions=1)
    if max(slice_data[1, :]) < 0.25:
        analyze_status = False
    if max(slice_data[5, :]) != 0.0:
        analyze_status = False
    if max(slice_data[8, :]) != 0.0:
        analyze_status = False

    # check density max and Vz and Bz components in tab sum
    sum_data = athena_read.tab(filename='bin/TestOutputs.block0.out3.00010.tab', raw=True,
                               dimensions=1)
    if max(sum_data[1, :]) < 15.0 and max(sum_data[:, 1]) > 20.0:
        analyze_status = False
    if max(sum_data[5, :]) != 0.0:
        analyze_status = False
    if max(sum_data[8, :]) != 0.0:
        analyze_status = False

    # assuming domain is 64x64 w/o ghost zones output, slice near central interface x2=0.0
    # check density max and Vz and Bz components in VTK dump
    xf, yf, _, vtk_data = athena_read.vtk(
        filename='bin/TestOutputs.block0.out4.00010.vtk')
    if max(vtk_data['rho'][0, 32, :]) < 0.25:
        analyze_status = False
    if max(vtk_data['vel'][0, 32, :, 2]) != 0.0:
        analyze_status = False
    if max(vtk_data['Bcc'][0, 32, :, 2]) != 0.0:
        analyze_status = False
    # if max(xf) != 0.5 and min(xf) != -0.5:
    #   analyze_status = False
    # if max(yf) != 0.5 and min(yf) != -0.5:
    #  analyze_status = False
    # logger.debug(str(vtk_data['rho'].shape))

    hdf5_data = athena_read.athdf('bin/TestOutputs.out5.00010.athdf', dtype=np.float32)
    if max(hdf5_data['rho'][0, 32, :]) < 0.25:
        analyze_status = False
    if max(hdf5_data['vel3'][0, 32, :]) != 0.0:
        analyze_status = False
    if max(hdf5_data['Bcc3'][0, 32, :]) != 0.0:
        analyze_status = False

    return analyze_status
