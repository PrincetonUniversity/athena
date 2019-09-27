# Test script for relativistic hydro shock tubes with passive scalars

# Modules
import logging
import numpy as np
import sys
import scripts.utils.athena as athena
sys.path.insert(0, '../../vis/python')
import athena_read  # noqa
athena_read.check_nan_flag = True
# set logger name based on module
logger = logging.getLogger('athena' + __name__[7:])


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure('s', prob='gr_shock_tube', coord='cartesian', flux='hllc',
                     nscalars='1', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    arguments = ['output1/variable=prim', 'time/ncycle_out=100']
    for n in range(1, 5):
        athena.run('hydro_sr/athinput.mb_'+repr(n), arguments)


# Analyze outputs
def analyze():

    # Approximate range of contact discontinuities
    x_lims = ((0.05, 0.15), (-0.15, -0.05), (0.225, 0.325), (0.3, 0.4))

    # Check that fractional concentration is correct on either side of contact
    analyze_status = True
    for n, x_lim in zip(range(1, 5), x_lims):
        data = athena_read.tab('bin/hydro_shock_rel_{0}.block0.out1.00001.tab'.format(n))
        x = data['x1v']
        chi = data['r0']
        chi_left = chi[np.where(x <= x_lim[0])[0]]
        chi_right = chi[np.where(x >= x_lim[1])[0]]
        if not np.allclose(chi_left, 0.0):
            analyze_status = False
        if not np.allclose(chi_right, 1.0):
            analyze_status = False
    return analyze_status
