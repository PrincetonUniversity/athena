# Regression test of user-defined boundary functions with all reconstruction orders
#
# A user-defined boundary function receives the primitive variables and fills the ghost
# cells. A user-defined zero-gradient condition must give exactly the same result as the
# built-in outflow condition. With fourth-order reconstruction (time/xorder=4), the
# boundary conditions are applied to the primitives computed in the current stage before
# the fourth-order conversion from the conserved variables.
#
# A 1D pressure pulse leaves the domain through both boundaries, and is run with the
# built-in outflow condition and with the user-defined one for the combinations of time
# integrator and time/xorder below. The final table outputs must be identical.

# Modules
import logging
import numpy as np
import scripts.utils.athena as athena
import sys
sys.path.insert(0, '../../vis/python')
import athena_read                             # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module

# time/integrator and time/xorder combinations to test
solvers = [('rk2', '2'), ('rk3', '3'), ('rk4', '4')]

# largest tolerated difference between the two boundary implementations: identical
# results print the same numbers
tolerance = 1.0e-12


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure(
        nghost=4,  # required for fourth-order configurations
        prob='user_outflow_test',
        coord='cartesian',
        flux='hllc', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    for (torder, xorder) in solvers:
        common = ['time/integrator=' + torder, 'time/xorder=' + xorder]
        athena.run('hydro/athinput.user_outflow_test',
                   common + ['job/problem_id=Builtin_' + torder])
        athena.run('hydro/athinput.user_outflow_test',
                   common + ['job/problem_id=User_' + torder,
                             'mesh/ix1_bc=user', 'mesh/ox1_bc=user', 'problem/user_bc=1'])


# Analyze outputs
def analyze():
    analyze_status = True
    for (torder, xorder) in solvers:
        builtin = athena_read.tab('bin/Builtin_{}.block0.out1.00001.tab'.format(torder),
                                  raw=True, dimensions=1)
        user = athena_read.tab('bin/User_{}.block0.out1.00001.tab'.format(torder),
                               raw=True, dimensions=1)
        diff = np.max(np.abs(builtin - user))
        logger.info('{} + time/xorder={}: max |user - built-in outflow| = {:.3e}'.format(
            torder.upper(), xorder, diff))
        if not diff <= tolerance:
            logger.warning('{} + time/xorder={}: the user-defined boundary function '
                           'differs from the built-in outflow condition'.format(
                               torder.upper(), xorder))
            analyze_status = False
    return analyze_status
