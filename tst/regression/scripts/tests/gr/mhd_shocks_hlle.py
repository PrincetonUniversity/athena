# Test script for relativistic MHD shock tubes in GR with HLLE

# Modules
import logging
import numpy as np
import sys
import scripts.utils.athena as athena
import scripts.utils.comparison as comparison
sys.path.insert(0, '../../vis/python')
import athena_read  # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure('bgt',
                     prob='gr_shock_tube',
                     coord='minkowski',
                     flux='hlle', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    arguments = [
      'job/problem_id=',
      'output1/file_type=vtk',
      'output1/variable=cons',
      'output1/dt=',
      'time/tlim=',
      'mesh/nx1=',
      'time/ncycle_out=0']
    times = [0.4, 0.55, 0.5]
    zones = [400, 800, 800]
    for i, time, zone in zip([1, 2, 4], times, zones):
        arguments_copy = list(arguments)
        arguments_copy[0] += 'gr_mhd_shock' + repr(i)
        arguments_copy[3] += repr(time)
        arguments_copy[4] += repr(time)
        arguments_copy[5] += repr(zone)
        athena.run('mhd_sr/athinput.mub_'+repr(i), arguments_copy)


# Analyze outputs
def analyze():
    analyze_status = True
    headers_ref = [('dens',), ('Etot',), ('mom', 0), ('mom', 1), ('mom', 2), ('cc-B', 0),
                   ('cc-B', 1), ('cc-B', 2)]
    headers_new = [('dens',), ('Etot',), ('mom', 0), ('mom', 1), ('mom', 2), ('Bcc', 0),
                   ('Bcc', 1), ('Bcc', 2)]
    tol_sets = [[0.02,  0.02,  0.03,  0.08, 0.0,   0.0, 0.02,  0.0],
                [0.003, 0.002, 0.007, 0.01, 0.007, 0.0, 0.005, 0.007],
                [0.003, 0.002, 0.03,  0.04, 0.008, 0.0, 0.002, 0.005]]
    for i, tols in zip([1, 2, 4], tol_sets):
        x_ref, _, _, data_ref = athena_read.vtk('data/sr_mhd_shock{0}_hlld.vtk'.format(i))
        x_new, _, _, data_new = \
            athena_read.vtk('bin/gr_mhd_shock{0}.block0.out1.00001.vtk'.format(i))
        for header_ref, header_new, tol in zip(headers_ref, headers_new, tols):
            array_ref = data_ref[header_ref[0]]
            if len(header_ref) == 1:
                array_ref = array_ref[0, 0, :]
            else:
                array_ref = array_ref[0, 0, :, header_ref[1]]
            array_new = data_new[header_new[0]]
            if len(header_new) == 1:
                array_new = array_new[0, 0, :]
            else:
                array_new = array_new[0, 0, :, header_new[1]]
            if header_new[0] == 'Etot':
                array_new = -array_new   # sign difference between SR and GR
            eps = comparison.l1_diff(x_ref, array_ref, x_new, array_new)
            if tol == 0.0:
                if eps > 0.0:
                    analyze_status = False
            else:
                eps /= comparison.l1_norm(x_ref, array_ref)
                if eps > tol or np.isnan(eps):
                    analyze_status = False
    return analyze_status
