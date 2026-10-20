# Regression test of the time integration of source terms
#
# A uniform gas at rest has no flux divergence, so it is driven only by source terms whose
# time integral is known exactly (problem generator source_time_test).
#
# Constant acceleration g: the momentum density after time t is exactly rho*g*t for every
# time integrator, because the source does not depend on time or on the state. This
# checks that all registers of the integrators receive the source terms with the right
# total weight (for SSPRK(5,4), including the register of its extra F(u^(3)) term).

# Modules
import logging
import scripts.utils.athena as athena
import sys
sys.path.insert(0, '../../vis/python')
import athena_read                             # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module

integrators = ['vl2', 'rk2', 'rk3', 'rk4', 'ssprk5_4']

# constant acceleration, and the largest tolerated relative error of the momentum
grav_acc1 = 1.0
tolerance_momentum = 1.0e-12


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure(prob='source_time_test', coord='cartesian', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    # constant acceleration
    for integrator in integrators:
        arguments = ['job/problem_id=Grav_' + integrator,
                     'time/integrator=' + integrator,
                     'hydro/grav_acc1=' + repr(grav_acc1)]
        athena.run('hydro/athinput.source_time_test', arguments)


# Analyze outputs
def analyze():
    analyze_status = True
    # constant acceleration: 1-momentum = rho*g*t, with rho=1 and a domain of length 1
    for integrator in integrators:
        data = athena_read.hst('bin/Grav_{}.hst'.format(integrator))
        time, momentum = data['time'][-1], data['1-mom'][-1]
        exact = grav_acc1*time
        error = abs(momentum - exact)/exact
        logger.info('{} constant acceleration: relative momentum error {:.2e} '
                    '(tolerance {:.1e})'.format(integrator, error, tolerance_momentum))
        if not error <= tolerance_momentum:
            logger.warning('{}: momentum {} instead of {}'.format(
                integrator, momentum, exact))
            analyze_status = False
    return analyze_status
