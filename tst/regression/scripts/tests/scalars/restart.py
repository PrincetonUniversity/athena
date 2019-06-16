# Regression test based on Newtonian hydro Sod shock tube problem
#
# Runs the Sod shock tube in x1, x2, and x3 directions successively, and checks errors
# against the analytic solution (which are computed by the executable automatically and
# stored in the temporary file shock_errors.dat)

# Modules
import logging
import scripts.utils.athena as athena
import sys
import os
from shutil import move
sys.path.insert(0, '../../vis/python')
import athena_read                             # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module
_exec = os.path.join('bin', 'athena')
_nscalars = [2]


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    global _nscalars
    for i in athena.global_config_args:
        tmp = i.split('=')
        if tmp[0] == '--nscalars' and len(tmp) == 2:
            _nscalars = [tmp[1]]
    for n in _nscalars:
        athena.configure(prob='shock_tube', coord='cartesian', nscalars=str(n), **kwargs)
        athena.make()
        move(_exec, _exec + '_' + str(n))
        os.system('cp -r obj obj_' + str(n))
    # Resuing obj/*.o may cause issues with Lcov, due to files without any coverage data.
    # Monitor this. E.g. obj_roe/ will contain hlle.o and hllc.o, but only roe.o is linked
    os.system('rm -rf obj')


# Run Athena++
def run(**kwargs):
    args = {'time/tlim': '0', 'output2/file_type': 'rst'}
    for n in _nscalars:
        move(_exec + '_' + str(n), _exec)
        os.system('mv obj_' + str(n) + ' obj')
        athena.run('hydro/athinput.sod', [i + '=' + args[i] for i in args],
                   lcov_test_suffix=str(n))
        athena.restart('Sod.00000.rst', [i + '=' + args[i] for i in args])
        os.system('rm -rf obj')
    return 'skip_lcov'


# Analyze outputs
def analyze():
    analyze_status = True
    return analyze_status
