# Regression test based on the Newtonian MHD linear wave convergence problem,
# 2D grid-aligned (horizontal, along x1, and vertical, along x2)

# Confirms fourth-order convergence rates of the RK3 + PPM + Laplacian flux correction
# + UCT4 (Felker & Stone 2018) solver configuration for grid-aligned MHD linear waves.
# The vertical direction exercises the x2-longitudinal code paths (UCT4x2, the x2 field
# conversions) that oblique and horizontal waves leave subdominant. Square cells
# (dx1f == dx2f) are REQUIRED by the 4th-order scheme and enforced by the code, so the
# horizontal runs use (nx1, nx2) = (N, N/2) and the vertical runs (2N, N) on the
# 1.0 x 0.5 domain.

# Modules
import logging
import scripts.utils.athena as athena
from math import log
import numpy as np
import sys
sys.path.insert(0, '../../vis/python')
import athena_read                             # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module

# Upper bounds on RMS-L1 errors for resolutions [32, 64, 128] (excluding lowest N=16)
# and lower bound on the convergence rate at N=64 (errors approach the ~1e-12 floor by
# N=128), per direction and wave mode.
# Tolerances are ~2x the errors measured on 2026-07 (macOS/AppleClang, -O3):
#   horizontal: fast 5.28e-10, 3.93e-11, 4.54e-12 / Alfven 2.80e-10, 1.80e-11, 1.16e-12
#               slow 3.05e-10, 1.93e-11, 3.91e-12 / entropy 2.86e-10, 1.84e-11, 1.17e-12
#   vertical:   fast 1.05e-09, 7.88e-11, 9.03e-12 (2 periods vs. 1 for horizontal);
#               Alfven/slow/entropy identical to horizontal to 4 significant digits
#               (exact x1/x2 symmetry of the scheme)
error_tols = {
    'horizontal': ((1.1e-9, 8.0e-11, 9.5e-12), (5.7e-10, 3.7e-11, 2.5e-12),
                   (6.2e-10, 4.0e-11, 8.0e-12), (5.8e-10, 3.8e-11, 2.5e-12)),
    'vertical':   ((2.2e-9, 1.6e-10, 1.9e-11), (5.7e-10, 3.7e-11, 2.5e-12),
                   (6.2e-10, 4.0e-11, 8.0e-12), (5.8e-10, 3.8e-11, 2.5e-12)),
}
# rates approach the error floor by N=128, so the rate check is at N=64
rate_tols = {
    'horizontal': (3.5, 3.5, 3.0, 3.5),
    'vertical':   (3.5, 3.5, 3.0, 3.5),
}

resolution_range = [16, 32, 64, 128]
num_nx1 = len(resolution_range)
wave_mode_names = ['fast', 'Alfven', 'slow', 'entropy']
# tlim = integer number of wave periods; lambda = 1.0 for horizontal (along x1) and
# 0.5 for vertical (along x2), so the periods lambda/c with c = 2, 1, 1/2, 1 for
# fast/Alfven/slow/entropy(vflow=1) differ per direction:
wave_mode_tlims = {'horizontal': [0.5, 1.0, 2.0, 1.0],
                   'vertical':   [0.5, 0.5, 1.0, 0.5]}
wave_mode_vflow = [0.0, 0.0, 0.0, 1.0]
directions = ['horizontal', 'vertical']


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure('b',
                     nghost=6,  # required for fourth-order MHD configurations
                     prob='linear_wave',
                     coord='cartesian',
                     flux='hlld', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    for direction in directions:
        for w in range(4):
            wf = w if w < 3 else 3
            for res in resolution_range:
                if direction == 'horizontal':
                    nx1, nx2 = res, res//2
                    extra = []
                else:
                    nx1, nx2 = 2*res, res
                    extra = ['problem/ang_3_vert=true']
                arguments = ['time/ncycle_out=0',
                             'time/xorder=4', 'time/integrator=rk3',
                             'problem/wave_flag={}'.format(wf),
                             'problem/vflow={}'.format(wave_mode_vflow[w]),
                             'mesh/nx1={}'.format(nx1), 'mesh/nx2={}'.format(nx2),
                             'time/tlim={}'.format(wave_mode_tlims[direction][w]),
                             'time/correct_err=true', 'time/correct_ic=true',
                             'output1/dt=-1', 'output2/dt=-1',
                             'problem/compute_error=true'] + extra
                athena.run('mhd/athinput.linear_wave2d_aligned', arguments)


# Analyze outputs
def analyze():
    analyze_status = True
    filename = 'bin/linearwave-errors.dat'
    data = np.array(athena_read.error_dat(filename))

    row = 0
    for direction in directions:
        err_tol = error_tols[direction]
        rate_tol = rate_tols[direction]
        logger.info('------------------------------')
        logger.info('RK3 + time/xorder=4, %s grid-aligned waves', direction)
        logger.info('------------------------------')
        for w, wave_mode in enumerate(wave_mode_names):
            logger.info("{} wave error convergence:".format(wave_mode.capitalize()))
            logger.info("N   |   rate   |   RMS-L1")
            rms_errs = data[row:row+num_nx1, 4]
            # the resolution along the wave direction:
            nres = (data[row:row+num_nx1, 0] if direction == 'horizontal'
                    else data[row:row+num_nx1, 1])
            row += num_nx1
            for i in range(1, num_nx1):
                rate = log(rms_errs[i-1]/rms_errs[i])/log(nres[i]/nres[i-1])
                logger.info("%d %g %g", int(nres[i]), rate, rms_errs[i])
                if (nres[i] == 64 and rate < rate_tol[w]):
                    logger.warning(
                        "{} {} wave converging at rate {} slower than {}".format(
                            direction, wave_mode, rate, rate_tol[w]))
                    analyze_status = False
                if (rms_errs[i] > err_tol[w][i-1]):
                    logger.warning(
                        "{} {} wave error {} is larger than tolerance {}".format(
                            direction, wave_mode, rms_errs[i], err_tol[w][i-1]))
                    analyze_status = False

    return analyze_status
