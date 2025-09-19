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

    # consistency check of all HDF5 outputs
    hdf5_data_5 = athena_read.athdf('bin/TestOutputs.out5.00010.athdf', dtype=np.float32)
    if max(hdf5_data_5['rho'][0, 32, :]) < 0.25:
        analyze_status = False
    if max(hdf5_data_5['vel3'][0, 32, :]) != 0.0:
        analyze_status = False
    if max(hdf5_data_5['Bcc3'][0, 32, :]) != 0.0:
        analyze_status = False

    out_types = {6: np.float32, 7: np.float64, 8: np.uint8, 9: np.uint16, 10: np.uint32,
                 11: np.uint64, 12: np.uint8}

    vmin, vmax = -2.0, 2.0  # set in athinput.test_outputs
    for idx, dtype in out_types.items():
        fn = f'bin/TestOutputs.out{idx}.00010.athdf'
        hdf5_data = athena_read.athdf(fn, dtype=np.float64)
        for key in ['rho', 'press', 'vel1', 'vel2', 'vel3', 'Bcc1', 'Bcc2', 'Bcc3']:
            if key not in hdf5_data:
                continue
            data = hdf5_data[key]
            if idx > 7:
                # de-normalize data for unsigned types
                bits = 2**(idx-5)
                if idx == 12:
                    bits = 8
                nmax = 2**bits - 1
                data = (data / nmax) * (vmax - vmin) + vmin
            atol = {np.uint8: 2e-2, np.uint16: 7e-5}.get(dtype, None)
            opt = {'atol': atol} if atol else {}
            if not np.allclose(data, hdf5_data_5[key], **opt):
                diff = data - hdf5_data_5[key]
                loc = np.unravel_index(np.argmax(np.abs(diff)), diff.shape)
                err = diff[loc]
                msg = f'Key {key} in out{idx} ({dtype}) does not match out5 data;'
                msg += f' max error = |{err:.3g}|.'
                logger.warning(msg)
                analyze_status = False

    return analyze_status
