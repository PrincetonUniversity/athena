# Regression test that a restart reproduces the time steps of an uninterrupted run
#
# The time step may grow by at most a factor of two from one cycle to the next
# (Mesh::NewTimeStep). A restart file stores the time step that was already chosen for
# the next cycle, and the time step limits are recomputed when the code restarts. The
# first step after a restart must be the same as in the uninterrupted run, also when
# the step is limited by its growth.
#
# The problem generator dt_growth_test limits the time step by a user-defined function
# to a small value first and then lets it grow by a factor of two per cycle up to the CFL
# limit. The problem is run for 14 cycles in one go ("Cont"), and in two pieces ("Rst"):
# the first stops after 8 cycles, in the middle of the growth, and the second restarts
# from its restart file. The time and time step of the two runs must agree.

# Modules
import logging
import scripts.utils.athena as athena
import sys
sys.path.insert(0, '../../vis/python')
import athena_read                             # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module

# largest tolerated relative difference of the time and the time step between the runs
tolerance = 1.0e-12


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure(prob='dt_growth_test', coord='cartesian', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    # uninterrupted run
    athena.run('hydro/athinput.dt_growth_test', ['job/problem_id=Cont', 'time/nlim=14'])
    # first piece: stop after 8 cycles; the restart file is written at the end of the run
    athena.run('hydro/athinput.dt_growth_test', ['job/problem_id=Rst', 'time/nlim=8'])
    # second piece: continue from the restart file written at the end of the first piece
    athena.restart('Rst.final.rst', ['time/nlim=14'])


# Analyze outputs
def analyze():
    analyze_status = True
    cont = athena_read.hst('bin/Cont.hst')
    rest = athena_read.hst('bin/Rst.hst')
    # Compare the last history row of each run: after 14 cycles the time is the sum of all
    # steps, and dt is the step the run would take next.
    for key in ('time', 'dt'):
        a, b = cont[key][-1], rest[key][-1]
        rel = abs(a - b)/abs(a)
        logger.info('{}: uninterrupted {:.14e}, restarted {:.14e}, relative difference '
                    '{:.2e}'.format(key, a, b, rel))
        if not rel <= tolerance:
            logger.warning('The restarted run differs from the uninterrupted run in '
                           '{} after 14 cycles'.format(key))
            analyze_status = False
    return analyze_status
