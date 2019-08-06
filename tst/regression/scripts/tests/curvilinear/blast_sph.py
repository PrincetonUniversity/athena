# Regression test to check whether blast wave remains spherical in spherical_polar coords

# Modules
import logging
import scripts.utils.athena as athena
import sys
sys.path.insert(0, '../../vis/python')
import athena_read  # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure(
        prob='blast',
        coord='spherical_polar', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    arguments = ['time/ncycle_out=0', 'problem/compute_error=true']
    athena.run('hydro/athinput.blast_sph', arguments)


# Analyze output
def analyze():
    analyze_status = True
    # read data from error file
    filename = 'bin/blastwave-shape.dat'
    data = athena_read.error_dat(filename)

    # check blast is spherical
    if data[0][3] > 1.0:
        logger.warning("Distortion of blast wave in spherical coords too large %g",
                       data[0][3])
        analyze_status = False

    return analyze_status
