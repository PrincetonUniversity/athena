# Regression test based on the swing amplification test
# Convergence of L2 norm of the density error is tested

# Modules
import logging
import scripts.utils.athena as athena
import numpy as np
import sys
sys.path.insert(0, '../../vis/python')
import athena_read  # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module

resolution_range = [32, 64]

# Set dimensionless input parameters
Lx = 1.0        # box size
Omega0 = 1.0    # angular freqency
qshear = 1.0    # shear rate
Q = 2           # Toomre Q
nJ = 2.5        # Jeans number
nwx = -3        # number of waves in x
nwy = 1         # number of waves in y
amp = 1e-6      # perturbation amplitude

# derived parameters
kappa = np.sqrt(4-2*qshear)*Omega0
cs = np.sqrt(4.0-2.0*qshear)/np.pi/nJ/Q
cs2 = cs**2
kx = 2*np.pi*nwx/Lx
ky = 2*np.pi*nwy/Lx
gconst = nJ*cs2


def prepare(*args, **kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure('mpi', 'fft', prob='msa',
                     grav='blockfft', *args, **kwargs)
    athena.make()


def run(**kwargs):
    for n in resolution_range:
        arguments = ['job/problem_id=SwingAmplification_' + repr(n),
                     'mesh/nx1=' + repr(n),
                     'mesh/x1min={}'.format(-Lx/2.),
                     'mesh/x1max={}'.format(Lx/2.),
                     'mesh/nx2=' + repr(n),
                     'mesh/x2min={}'.format(-Lx/2.),
                     'mesh/x2max={}'.format(Lx/2.),
                     'meshblock/nx1=' + repr(n/2),
                     'meshblock/nx2=' + repr(n/2),
                     'hydro/iso_sound_speed={}'.format(cs),
                     'orbital_advection/qshear={}'.format(qshear),
                     'orbital_advection/Omega0={}'.format(Omega0),
                     'problem/Q={}'.format(Q),
                     'problem/nJ={}'.format(nJ),
                     'problem/amp={}'.format(amp),
                     'problem/nwx={}'.format(nwx),
                     'problem/nwy={}'.format(nwy),
                     'problem/compute_error=true']
        # athena.run('mhd/athinput.msa', arguments)
        athena.mpirun(kwargs['mpirun_cmd'], kwargs['mpirun_opts'],
                      4, 'mhd/athinput.msa', arguments)


def analyze():
    # read data from error file
    filename = 'bin/msa-errors.dat'
    data = athena_read.error_dat(filename)
    analyze_status = True

    # Check absolute error and convergence rate lower bounds of all waves
    # Asymptotic second-order convergence should have ratio <= 0.25 below
    if data[1][4] > 2.3e-7:
        logger.warning("error in the shearing wave too large %g",
                       data[1][4])
        analyze_status = False
    if data[1][4]/data[0][4] > 0.25:
        logger.warning("not converging for shearing wave %g %g",
                       data[0][4], data[1][4])
        analyze_status = False

    return analyze_status
