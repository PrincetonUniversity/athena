# Regression test based on 2D advection of a weak magnetic field loop, xorder=4

# Advects the field loop diagonally across a periodic 2D domain with the RK3 + PPM +
# Laplacian flux correction + UCT4 (Felker & Stone 2018, JCP paper sec. 5.3) solver
# configuration. Checks that: the out-of-plane field energy remains identically zero
# (a strict test of the corner EMF construction: v3 = B3 = 0 implies E1 = E2 = 0),
# the magnetic energy of the loop is well preserved (low numerical dissipation of the
# fourth-order scheme), and no NaN appear.

# Modules
import logging
import scripts.utils.athena as athena
import sys
import numpy as np
sys.path.insert(0, '../../vis/python')
import athena_read                             # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module

# minimum acceptable fraction of the initial magnetic energy retained at t=2
# (one full diagonal crossing); measured on 2026-07: 0.8983 (RK3+UCT4, 128x64, hlld)
_me_retention_tol = 0.85
# out-of-plane field must remain at the roundoff floor; measured max 3-ME ~ 9e-35
_me3_tol = 1e-30


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure('b',
                     nghost=6,  # required for fourth-order MHD configurations
                     prob='field_loop',
                     coord='cartesian',
                     flux='hlld', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    arguments = ['time/ncycle_out=0',
                 'time/xorder=4', 'time/integrator=rk3',
                 'time/tlim=2.0',
                 'output1/dt=0.05',   # hst
                 'output2/dt=-1']
    athena.run('mhd/athinput.field_loop', arguments)


# Analyze outputs
def analyze():
    analyze_status = True
    hst = athena_read.hst('bin/Loop.hst')

    # out-of-plane magnetic energy must remain at the roundoff floor
    me3_max = np.abs(hst['3-ME']).max()
    logger.info('max out-of-plane magnetic energy: %e', me3_max)
    if me3_max > _me3_tol:
        logger.warning('out-of-plane magnetic energy grew beyond roundoff: {}'.format(
            me3_max))
        analyze_status = False

    # in-plane magnetic energy retention of the advected loop
    me0 = hst['1-ME'][0] + hst['2-ME'][0]
    mef = hst['1-ME'][-1] + hst['2-ME'][-1]
    retention = mef/me0
    logger.info('magnetic energy retention ME(t=2)/ME(0) = %f', retention)
    if retention < _me_retention_tol or retention > 1.0 + 1e-10:
        logger.warning('magnetic energy retention {} outside ({}, 1.0]'.format(
            retention, _me_retention_tol))
        analyze_status = False

    return analyze_status
