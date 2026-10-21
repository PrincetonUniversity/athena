# Regression test of the time integration of source terms
#
# A uniform gas at rest has no flux divergence, so it is driven only by source terms whose
# time integral is known exactly (problem generator source_time_test).
#
# 1. Constant acceleration g: the momentum density after time t is exactly rho*g*t for
#    every time integrator, because the source does not depend on time or on the state.
#    This checks that all registers of the integrators receive the source terms with the
#    right total weight (for SSPRK(5,4), including the register of its extra F(u^(3))
#    term).
#
# 2. User-defined energy source cos(omega*t): the total energy after time t is
#    E(0) + sin(omega*t)/omega. The source terms are evaluated at the time of the start
#    of each stage, so the error converges at the order of the time integrator only if
#    the stage times are right. The time step is varied through the CFL number.

# Modules
import logging
from math import log, sin
import scripts.utils.athena as athena
import sys
sys.path.insert(0, '../../vis/python')
import athena_read                             # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module

integrators = ['vl2', 'rk2', 'rk3', 'rk4', 'ssprk5_4']

# 1. constant acceleration, and the largest tolerated relative error of the momentum
grav_acc1 = 1.0
tolerance_momentum = 1.0e-12

# 2. time-dependent energy source: integrators with fourth-order accuracy, the amplitude
# and angular frequency of the source, the end time, the CFL numbers used to vary the
# time step, and the smallest tolerated convergence rate
order4_integrators = ['rk4', 'ssprk5_4']
src_amp = 1.0
src_omega = 6.0
src_tlim = 1.0
cfl_numbers = [0.8, 0.4, 0.2]
min_rate = 3.8


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
    # time-dependent energy source
    for integrator in order4_integrators:
        for cfl in cfl_numbers:
            arguments = ['job/problem_id=Cos_{}_{}'.format(integrator, cfl),
                         'time/integrator=' + integrator, 'time/cfl_number=' + repr(cfl),
                         'time/tlim=' + repr(src_tlim),
                         'problem/src_amp=' + repr(src_amp),
                         'problem/src_omega=' + repr(src_omega)]
            athena.run('hydro/athinput.source_time_test', arguments)


# Analyze outputs
def analyze():
    analyze_status = True
    # 1. constant acceleration: 1-momentum = rho*g*t, with rho=1 and a domain of length 1
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

    # 2. time-dependent energy source: total energy = E(0) + src_amp*sin(omega*t)/omega
    for integrator in order4_integrators:
        errors = []
        for cfl in cfl_numbers:
            data = athena_read.hst('bin/Cos_{}_{}.hst'.format(integrator, cfl))
            time, energy = data['time'][-1], data['tot-E'][-1]
            exact = data['tot-E'][0] + src_amp*sin(src_omega*time)/src_omega
            errors.append(abs(energy - exact))
        logger.info('{} time-dependent source: energy errors {} for CFL numbers '
                    '{}'.format(integrator, ['{:.2e}'.format(e) for e in errors],
                                cfl_numbers))
        for i in range(1, len(cfl_numbers)):
            rate = log(errors[i-1]/errors[i])/log(cfl_numbers[i-1]/cfl_numbers[i])
            logger.info('{}: convergence rate {:.2f} (minimum {})'.format(
                integrator, rate, min_rate))
            if not rate >= min_rate:
                logger.warning('{}: the energy error converges at the rate {:.2f}, '
                               'slower than {}'.format(integrator, rate, min_rate))
                analyze_status = False
    return analyze_status
