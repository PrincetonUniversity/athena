# Regression test based on the MHD Kelvin-Helmholtz instability of Lecoanet et al.
# (2016) with a weak aligned field, xorder=4 (thesis sec. 5.3.3)

# Runs the smooth double-shear-layer KHI at reduced resolution with the RK3 + PPM +
# Laplacian flux correction + UCT4 solver configuration, including explicit viscosity,
# resistivity, thermal conduction, and passive scalar (dye) diffusion at Re = 10^5, and
# the fourth-order accurate initial condition (analytic face-averaged B1, correct_ic).
# Checks robustness invariants: run completes without NaN, conservation of total mass
# and energy (fully periodic domain), and boundedness of the dye concentration.

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
                     prob='kh',
                     nscalars='1',
                     coord='cartesian',
                     flux='hlld', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    arguments = ['time/ncycle_out=0',
                 'time/xorder=4', 'time/integrator=rk3',
                 'time/correct_ic=true',
                 'time/tlim=2.0',
                 'mesh/nx1=64', 'mesh/nx2=128',
                 'meshblock/nx1=32', 'meshblock/nx2=64',
                 'output1/dt=0.05',   # hst
                 'output2/dt=2.0', 'output2/file_type=vtk',  # vtk prim
                 'output3/dt=-1']     # disable rst
    athena.run('mhd/athinput.kh-shear-lecoanet', arguments)


# Analyze outputs
def analyze():
    analyze_status = True

    hst = athena_read.hst('bin/kh-shear-lecoanet.hst')
    for q in ('mass', 'tot-E'):
        rel = np.abs(hst[q][-1]/hst[q][0] - 1.0)
        logger.info('%s relative change over run: %e', q, rel)
        if rel > 1e-4:
            logger.warning('total {} not conserved: relative change {}'.format(q, rel))
            analyze_status = False

    # dye concentration boundedness in the final state (r0 in [0,1] up to small
    # over/undershoots from the unlimited 4th-order corrections)
    for blk in range(4):
        xf, yf, zf, data = athena_read.vtk(
            'bin/kh-shear-lecoanet.block{}.out2.00001.vtk'.format(blk))
        if data['rho'].min() <= 0.0 or data['press'].min() <= 0.0:
            logger.warning('non-positive rho/press in final state')
            analyze_status = False
        if 'r0' in data:
            cmin, cmax = data['r0'].min(), data['r0'].max()
            if cmin < -0.01 or cmax > 1.01:
                logger.warning('dye concentration out of bounds: [{}, {}]'.format(
                    cmin, cmax))
                analyze_status = False

    return analyze_status
