# Regression test based on Newtonian hydro linear wave convergence problem

# Run grid-aligned linear waves in (x1) 1D, (x2) 2D, (x3) 3D and confirm symmetry in
# errors. Uniform square grid, midpoint approximation used for initialization and in error
# calculations. Setting CFL=0.3 for all runs.

# Modules
import logging
import scripts.utils.athena as athena
import numpy as np
import sys
sys.path.insert(0, '../../vis/python')
import athena_read                             # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module


# List of time/integrator and time/xorder combinations to test:
solvers = [('vl2', '2'), ('rk3', '4c'), ('ssprk5_4', '4')]

# Grid cells in the wavevector direction
resolution_range = [32]
num_nx1 = len(resolution_range)
# Number of times Athena++ is run for each above configuration:
nrows_per_solver = 3*num_nx1


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
    for (torder, xorder) in solvers:
        # L-going sound wave
        for i in resolution_range:
            # 1D domain with Lx1=1.0, wavevector points along x1:
            arguments = ['time/ncycle_out=0', 'time/cfl_number=0.3',
                         'time/xorder=' + xorder, 'time/integrator=' + torder,
                         'problem/wave_flag=0', 'problem/vflow=0.0',
                         'mesh/nx1=' + repr(i),
                         'mesh/x1max=1.0',
                         'output2/dt=-1', 'time/tlim=1.0',
                         'problem/compute_error=true']
            athena.run('hydro/athinput.linear_wave1d', arguments)

            # 2D domain, 2.0 x 1.0, wavevector points along x2:
            arguments = ['time/ncycle_out=0', 'time/cfl_number=0.3',
                         'time/xorder=' + xorder, 'time/integrator=' + torder,
                         'problem/wave_flag=0', 'problem/vflow=0.0',
                         'mesh/nx1=' + repr(i*2), 'mesh/nx2=' + repr(i),
                         'mesh/x1max=2.0', 'mesh/x2max=1.0',
                         'problem/ang_3_vert=true',
                         'output2/dt=-1', 'time/tlim=1.0',
                         'problem/compute_error=true']
            athena.run('hydro/athinput.linear_wave2d', arguments)

            # 3D domain, 2.0 x 1.0 x 1.0, wave vector points along x3:
            arguments = ['time/ncycle_out=0', 'time/cfl_number=0.3',
                         'time/xorder=' + xorder, 'time/integrator=' + torder,
                         # 'time/correct_ic=true', 'time/correct_err=true',
                         'problem/wave_flag=0', 'problem/vflow=0.0',
                         'mesh/nx1=' + repr(i*2), 'mesh/nx2=' + repr(i),
                         'mesh/nx3=' + repr(i),
                         'meshblock/nx1=' + repr(i*2), 'meshblock/nx2=' + repr(i),
                         'meshblock/nx3=' + repr(i),
                         'mesh/x1max=2.0', 'mesh/x2max=1.0', 'mesh/x3max=1.0',
                         'problem/ang_2_vert=true',
                         'output2/dt=-1', 'time/tlim=1.0',
                         'problem/compute_error=true']
            athena.run('hydro/athinput.linear_wave3d', arguments)


# Analyze outputs
def analyze():
    analyze_status = True
    # read data from error file
    filename = 'bin/linearwave-errors.dat'
    data = athena_read.error_dat(filename)

    for (torder, xorder) in solvers:
        # effectively list.pop() range of rows for this solver configuration
        solver_results = np.array(data[0:nrows_per_solver])
        # del data[0:nrows_per_solver]

        data = np.delete(data, np.s_[0:nrows_per_solver], 0)

        logger.debug('{} + {}'.format(torder.upper(), xorder))
        # L-going sound wave: Ncycle, RMS-L1, d_L1, E_L1, d_max, E_max
        indices = [3, 4, 5, 9, 11, 15]
        results_1D = np.take(np.squeeze(solver_results[0, :]), indices)
        results_2D = np.take(np.squeeze(solver_results[1, :]), indices)
        results_3D = np.take(np.squeeze(solver_results[2, :]), indices)

        # Check that results are nearly identical for sound waves along each dir/ndims
        # (Differences in d_max can exceed 2e-15, but other differences are 2e-16 at most)
        atol = 5e-15
        rtol = 1e-8

        # Useful optional diagnostics for determining if differences are meaningful in FP:
        logger.debug(str(np.allclose(results_1D, results_2D, atol=atol, rtol=rtol)))
        logger.debug("numpy tolerance = {} %g".format(atol + 1e-10*abs(results_2D)))
        logger.debug(np.allclose(results_2D, results_3D, atol=atol, rtol=rtol))
        logger.debug("numpy tolerance = {} %g".format(atol + 1e-10*abs(results_3D)))
        logger.debug(" ".join(map(str, [results_1D, results_2D])))
        logger.debug(str(results_1D - results_2D))

        if not (np.allclose(results_1D, results_2D, atol=atol, rtol=rtol)
                and np.allclose(results_2D, results_3D, atol=atol, rtol=rtol)):
            logger.warning("1D/2D/3D grid-aligned sound wave results:")
            logger.warning("Ncycle  RMS-L1-Error  d_L1  E_L1  d_max  E_max")
            logger.warning(str(results_1D))
            logger.warning(str(results_2D))
            logger.warning(str(results_3D))
            logger.warning("Exhibit differences that are not close to round-off")
            analyze_status = False

    return analyze_status
