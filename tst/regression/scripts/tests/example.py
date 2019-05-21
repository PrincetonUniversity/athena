"""
Example test script.

This is a complete, working example that can be run as part of the test suite. It does a
simple test of a relativistic shock tube using the GR framework. There are many comments
in order to make this file self-explanatory, but the actual working code is only 28
lines long.

There are three functions defined here:
    prepare(**kwargs)
    run(**kwargs)
    analyze()
All three must be defined with the same names and no required inputs in order to make a
working script. They are called in sequence from the main test script run_tests.py.
Additional support functions can be defined here, to be called by the three primary fns.

Heavy use is made of support utilities defined in scripts/utils/athena.py. These are
general-purpose Python scripts that interact with Athena++. They should be used whenever
possible, since they work together to compile and run Athena++ and read the output data.
In particular, proper use of them will result in all files outside tst/regression/ being
in the same state after the test as they were before (including whatever configured
version of Athena++ existed in athena/bin/), as well as cleaning up any new files
produced by the test.
"""

# Modules
import numpy as np                             # standard Python module for numerics
import sys                                     # standard Python module to change path
import scripts.utils.athena as athena          # utilities for running Athena++
import scripts.utils.comparison as comparison  # more utilities explicitly for testing
sys.path.insert(0, '../../vis/python')         # insert path to Python read scripts
import athena_read                             # utilities for reading Athena++ data # noqa
athena_read.check_nan_flag = True              # raise exception when encountering NaNs


def prepare(**kwargs):
    """
    Configure and make the executable.

    This function is called first. It is responsible for calling the configure script and
    make to create an executable. It takes no inputs and produces no outputs.
    """

    # Configure as though we ran
    #     python configure.py -g -t --prob=shock_tube_gr --coord=minkowski
    # from the athena/ directory. Note that additional -<flag> command-line arguments can
    # be specified as additional '<flag>' arguments before the <key>='<value>' arguments
    # to athena.configure(). Any number of --<key>=<value> command-line arguments can also
    # be supplied. Note athena.configure() expects the values only to be quoted, e.g.
    # --<key>='<value>'.
    athena.configure('g', 't',
                     prob='gr_shock_tube',
                     coord='minkowski',
                     **kwargs)

    # Call make as though we ran
    #     make clean
    #     make
    # from the athena/ directory.
    athena.make()


def run(**kwargs):
    """
    Run the executable.

    This function is called second. It is responsible for calling the Athena++ binary in
    such a way as to produce testable output. It takes no inputs and produces no outputs.
    """

    # Create list of runtime arguments to override the athinput file. Each element in the
    # list is simply a string of the form '<block>/<field>=<value>', where the contents of
    # the string are exactly what one would type on the command line run running Athena++.
    arguments = ['time/ncycle_out=0',
                 'job/problem_id=gr_shock_tube',
                 'output1/file_type=vtk',
                 'output1/variable=cons',
                 'output1/dt=0.4',
                 'time/cfl_number=0.4',
                 'time/tlim=0.4',
                 'mesh/nx1=400']

    # Run Athena++ as though we called
    #     ./athena -i ../inputs/hydro_sr/athinput.mb_1 job/problem_id=gr_shock_tube <...>
    # from the bin/ directory. Note we omit the leading '../inputs/' below when specifying
    # the athinput file.)
    athena.run('hydro_sr/athinput.mb_1', arguments)
    # No return statement/value is ever required from run(), but returning anything other
    # than default None will cause run_tests.py to skip executing the optional Lcov cmd
    # immediately after this module.run() finishes, e.g. if Lcov was already invoked by:
    # athena.run('hydro_sr/athinput.mb_1', arguments, lcov_test_suffix='mb_1')
    return 'skip_lcov'


def analyze():
    """
    Analyze the output and determine if the test passes.

    This function is called third; nothing from this file is called after it. It is
    responsible for reading whatever data it needs and making a judgment about whether or
    not the test passes. It takes no inputs. Output should be True (test passes) or False
    (test fails).
    """

    # Read in reference data. The tst/regression/data/ directory has reference runs for
    # comparing future output of the code. We only need to specify file names starting
    # with "data/". Now athena_read.vtk() returns four objects: the x-interface locations,
    # the y-interface locations, the z-interface locations, and the values of the
    # variables themselves. In a 1D problem we ignore the second and third returned
    # values, assigning them to the _ variable as is typical Python style.
    x_ref, _, _, data_ref = athena_read.vtk('data/sr_hydro_shock1_hlle.vtk')

    # Read in the data produced during this test. This will usually be stored in the
    # tst/regression/bin/ directory, but again we omit the first part of the path. Note
    # the file name is what we expect based on the job/problem_id field supplied in run().
    x_new, _, _, data_new = athena_read.vtk('bin/gr_shock_tube.block0.out1.00001.vtk')

    # Extract the quantities of interest. Suppose we want to check that the total energy
    # and the x-momentum are the same as those given in the reference dataset. The fourth
    # object returned by athena_read.vtk() is a dictionary of 3D (scalars) or 4D (vectors)
    # NumPy arrays, whose keys ('Etot' and 'mom' in this case) are exactly the names of
    # the arrays as stored in the vtk file. Here we extract the reference values, where
    # the fourth index specifies which component of the vector quantity to extract. The
    # choice of slicing will give us 1D arrays without any singleton dimensions.
    e_ref = data_ref['Etot'][0, 0, :]
    mx_ref = data_ref['mom'][0, 0, :, 0]

    # Similarly, we extract the newly created values.
    e_new = -data_new['Etot'][0, 0, :]   # sign flip between SR and GR definitions
    mx_new = data_new['mom'][0, 0, :, 0]

    # Next we compute the differences between the reference arrays and the newly created
    # ones in the L^1 sense. That is, given functions f and g, we want
    #     \int |f(x)-g(x)| dx.
    # The utility script comparison.l1_diff() does this exactly, conveniently taking N+1
    # interface locations and N volume-averaged quantities. The two datasets can have
    # different values of N.
    error_abs_e = comparison.l1_diff(x_ref, e_ref, x_new, e_new)
    error_abs_mx = comparison.l1_diff(x_ref, mx_ref, x_new, mx_new)

    # The errors are more meaningful if we account for the length of the domain and the
    # typical magnitude of the function itself. Fortunately, comparison.l1_norm() computes
    #     \int |f(x)| dx.
    # (Note neither comparison.l1_diff() nor comparison.l1_norm() divides by the length of
    # the domain.)
    error_rel_e = error_abs_e / comparison.l1_norm(x_ref, e_ref)
    error_rel_mx = error_abs_mx / comparison.l1_norm(x_ref, mx_ref)

    # Finally, we test that the relative errors in the two quantities are no more than 1%.
    # If they are, we return False at the very end of the function and file; otherwise
    # we return True. NumPy provides a way of checking if the error is NaN, which also
    # indicates something went wrong. The same check can (and should) be enabled
    # automatically at the point of reading the input files via the athena_read.py
    # functions by setting "athena_read.check_nan_flag=True" (as done at the top of this
    # file). Regression test authors should keep in mind the caveats of floating-point
    # calculations and perform multiple checks for NaNs when necessary.

    # The main test script will record the result and delete both tst/regression/bin/ and
    # obj/ folders before proceeding on to the next test.
    analyze_status = True
    if error_rel_e > 0.01 or np.isnan(error_rel_e):
        analyze_status = False
    if error_rel_mx > 0.01 or np.isnan(error_rel_mx):
        analyze_status = False

    # Note, if the problem generator in question outputs a unique CSV file containing
    # quantitative error measurements (e.g. --prob=linear_wave outputs
    # linearwave-errors.dat when problem/compute_error=true at runtime), then these values
    # can also be input and used in this analyze() function. It is recommended to use:
    # athena_read.error_dat('bin/linearwave-errors.dat')
    # This wrapper function to np.loadtxt() can automatically check for the presence of
    # NaN values as in the other athena_read.py functions.

    return analyze_status
