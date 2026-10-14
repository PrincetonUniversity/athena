# Regression test based on the 2D Orszag-Tang vortex, xorder=4

# Runs the Orszag-Tang vortex with the RK3 + PPM + Laplacian flux correction + UCT4
# (Felker & Stone 2018, JCP paper sec. 5.6) solver configuration and checks robustness
# invariants: the run completes without NaN, pressure and density stay positive, and
# total mass and energy are conserved (fully periodic domain).

# Modules
import logging
import scripts.utils.athena as athena
import sys
import numpy as np
sys.path.insert(0, '../../vis/python')
import athena_read                             # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure('b',
                     nghost=6,  # required for fourth-order MHD configurations
                     prob='orszag_tang',
                     coord='cartesian',
                     flux='hlld', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    arguments = ['time/ncycle_out=0',
                 'time/xorder=4', 'time/integrator=rk3',
                 'time/tlim=1.0',
                 'mesh/nx1=128', 'mesh/nx2=128',
                 'meshblock/nx1=64', 'meshblock/nx2=64',
                 'output1/dt=0.01',   # hst
                 'output2/dt=1.0']    # vtk prim
    athena.run('mhd/athinput.orszag-tang', arguments)


# Analyze outputs
def analyze():
    analyze_status = True

    # conservation of total mass and energy (fully periodic domain); the .hst output
    # format limits the achievable precision of this comparison
    hst = athena_read.hst('bin/OrszagTang.hst')
    for q in ('mass', 'tot-E'):
        rel = np.abs(hst[q][-1]/hst[q][0] - 1.0)
        logger.info('%s relative change over run: %e', q, rel)
        if rel > 1e-4:
            logger.warning('total {} not conserved: relative change {}'.format(q, rel))
            analyze_status = False

    # positivity and NaN check of the final state (NaNs raise via check_nan_flag)
    for blk in (0, 1, 2, 3):
        xf, yf, zf, data = athena_read.vtk(
            'bin/OrszagTang.block{}.out2.00001.vtk'.format(blk))
        if data['rho'].min() <= 0.0 or data['press'].min() <= 0.0:
            logger.warning('non-positive rho/press in final state: min rho {}, '
                           'min press {}'.format(data['rho'].min(),
                                                 data['press'].min()))
            analyze_status = False

    return analyze_status
