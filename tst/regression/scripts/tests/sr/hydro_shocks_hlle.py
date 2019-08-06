# Test script for relativistic hydro shock tubes with HLLE

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
    athena.configure('s',
                     prob='gr_shock_tube',
                     coord='cartesian',
                     flux='hlle', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    arguments = [
      '',
      'output1/file_type=vtk',
      'output1/variable=cons',
      'output1/dt=0.4',
      'time/cfl_number=0.4',
      'time/tlim=0.4',
      'mesh/nx1=400',
      'time/ncycle_out=100']
    for i in range(1, 5):
        arguments[0] = 'job/problem_id=sr_hydro_shock' + repr(i)
        athena.run('hydro_sr/athinput.mb_'+repr(i), arguments)


# Analyze outputs
def analyze():
    analyze_status = True
    headers = [('dens',), ('Etot',), ('mom', 0)]
    tols = [[0.02, 0.01, 0.01], [0.01, 0.01, 0.02], [0.01, 0.01, 0.02], [0.5, 0.01, 0.02]]
    for i in range(1, 5):
        x_ref, _, _, data_ref = athena_read.vtk(
            'data/sr_hydro_shock{0}_hllc.vtk'.format(i))
        x_new, _, _, data_new = athena_read.vtk(
            'bin/sr_hydro_shock{0}.block0.out1.00001.vtk'.format(i))
        tols_particular = tols[i-1]
        for header, tol in zip(headers, tols_particular):
            array_ref = data_ref[header[0]]
            array_ref = array_ref[0, 0, :] if len(
                header) == 1 else array_ref[0, 0, :, header[1]]
            array_new = data_new[header[0]]
            array_new = array_new[0, 0, :] if len(
                header) == 1 else array_new[0, 0, :, header[1]]
            eps = comparison.l1_diff(x_ref, array_ref, x_new, array_new)
            eps /= comparison.l1_norm(x_ref, array_ref)
            if eps > tol or np.isnan(eps):
                analyze_status = False
    return analyze_status
