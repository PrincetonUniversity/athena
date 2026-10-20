# Regression test that a restart reproduces an uninterrupted run with time/correct_ic
#
# With time/correct_ic=true, Mesh::Initialize converts the midpoint (cell-centered point)
# values set by the problem generator into cell averages. A restart file already stores
# cell averages, so this conversion is done only when the problem generator has run.
# The test restarts a run and checks that the result is the same as without a restart.
#
# A large-amplitude 1D sound wave is run with RK4 + PPM (xorder=4) and correct_ic=true to
# t=0.2 in one go ("Cont"), and in two pieces ("Rst"): the first stops after 12 cycles,
# and the second restarts from its restart file. The first piece is limited by the
# number of cycles, not by time, so that none of its steps is shortened to reach an end
# time. The final table output of the two runs must agree.

# Modules
import glob
import logging
import numpy as np
import scripts.utils.athena as athena
import sys
sys.path.insert(0, '../../vis/python')
import athena_read                             # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module

# Arguments shared by both runs. A large-amplitude wave on a coarse grid makes the effect
# of the correction (relative size ~1e-3) much larger than the 5-digit precision of the
# table output.
_common = ['time/ncycle_out=0', 'time/xorder=4', 'time/integrator=rk4',
           'time/correct_ic=true', 'time/correct_err=false',
           'problem/wave_flag=0', 'problem/amp=0.5', 'problem/compute_error=false',
           'mesh/nx1=32', 'output2/file_type=tab', 'output2/dt=0.2']

# largest tolerated difference between the two runs: identical runs print identical
# numbers, while a change of the cell averages by the correction shows at the 1e-4 to
# 1e-3 level
tolerance = 1.0e-6


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure(
        nghost=4,  # required for fourth-order configurations
        prob='linear_wave',
        coord='cartesian',
        flux='hllc', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    # uninterrupted run
    arguments = _common + ['job/problem_id=Cont', 'time/tlim=0.2',
                           'output1/dt=-1']
    athena.run('hydro/athinput.linear_wave1d', arguments)
    # first piece: stop after 12 cycles; a restart file is written at the end of every
    # run (output1 is the history output by default, so use it for a restart output)
    arguments = _common + ['job/problem_id=Rst', 'time/tlim=0.2', 'time/nlim=12',
                           'output1/file_type=rst', 'output1/dt=1.0']
    athena.run('hydro/athinput.linear_wave1d', arguments)
    # second piece: continue from the restart file written at the end of the first piece
    athena.restart('Rst.final.rst', ['time/nlim=-1'])


def _last_table(problem_id):
    """Return the data of the last table output of the given run and its time."""
    filename = sorted(glob.glob('bin/{}.block0.out2.*.tab'.format(problem_id)))[-1]
    with open(filename) as table:
        time = float(table.readline().split('time=')[1].split()[0])
    return athena_read.tab(filename, raw=True, dimensions=1), time


# Analyze outputs
def analyze():
    analyze_status = True
    cont, t_cont = _last_table('Cont')
    rest, t_rest = _last_table('Rst')
    if not (abs(t_cont - 0.2) < 1e-12 and abs(t_rest - 0.2) < 1e-12):
        logger.warning('The last table outputs are at t={} and t={}, not 0.2'.format(
            t_cont, t_rest))
        return False
    diff = np.max(np.abs(cont - rest))
    logger.info('max |restarted - uninterrupted| = {:.3e} (tolerance {:.1e})'.format(
        diff, tolerance))
    if not diff <= tolerance:
        logger.warning('The restarted run differs from the uninterrupted run')
        analyze_status = False
    return analyze_status
