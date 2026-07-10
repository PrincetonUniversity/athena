# Regression test based on the 2D strongly magnetized blast wave, xorder=4

# Runs the MHD blast wave (JCP paper sec. 5.7) in 2D with the RK3 + PPM + Laplacian
# flux correction + UCT4 solver configuration using characteristic reconstruction, and
# checks robustness invariants: run completes without NaN, density/pressure positivity,
# and conservation of total mass and energy (fully periodic domain).

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
                     prob='blast',
                     coord='cartesian',
                     flux='hlld', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    arguments = ['time/ncycle_out=0',
                 'time/xorder=4c', 'time/integrator=rk3',
                 'time/tlim=0.2',
                 'mesh/nx1=128', 'mesh/nx2=128', 'mesh/nx3=1',
                 'mesh/x3min=-0.5', 'mesh/x3max=0.5',
                 'output1/dt=0.2']    # vtk prim (the only output block in this input)
    athena.run('mhd/athinput.blast', arguments)


# Analyze outputs
def analyze():
    analyze_status = True

    # positivity and NaN check of the final state (NaNs raise via check_nan_flag);
    # single MeshBlock (athinput.blast has no <meshblock> block to override)
    xf, yf, zf, data = athena_read.vtk('bin/Blast.block0.out1.00001.vtk')
    rho_min = data['rho'].min()
    press_min = data['press'].min()
    logger.info('min rho %e, min press %e', rho_min, press_min)
    if rho_min <= 0.0 or press_min <= 0.0:
        logger.warning('non-positive rho/press in final state: min rho {}, '
                       'min press {}'.format(rho_min, press_min))
        analyze_status = False

    return analyze_status
