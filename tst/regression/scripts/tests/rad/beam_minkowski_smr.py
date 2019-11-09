# Test script for GR radiation - beam in Minkowski spacetime with SMR

# Modules
import logging
import numpy as np
import sys
import scripts.utils.athena as athena
import scripts.utils.comparison as comparison
sys.path.insert(0, '../../vis/python')
import athena_read  # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure('g', 'r', 'hdf5', prob='gr_rad_beam', coord='minkowski', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    input_file = 'rad_gr/athinput.beam_cart'
    job_str = 'job/problem_id={0}'
    ref_str = 'mesh/refinement={0}'
    athena.run(input_file, [job_str.format('uniform'), ref_str.format('none')])
    athena.run(input_file, [job_str.format('refined'), ref_str.format('static')])


# Analyze outputs
def analyze():

    # Read data
    data_uniform = athena_read.athdf('bin/uniform.rad.00010.athdf')
    data_refined = athena_read.athdf('bin/refined.rad.00010.athdf', level=0)
    energy_uniform = data_uniform['R00']
    energy_refined = data_refined['R00']

    # Check extrema
    if np.min(energy_uniform) > 0.0 or np.min(energy_refined) > 0.0:
        logger.warning('did not find cells with no radiation')
        return False
    if not np.isclose(np.max(energy_uniform), 0.82753056, rtol=0.01):
        logger.warning('uniform grid has wrong maximum energy')
        return False
    if not np.isclose(np.max(energy_refined), 0.5032955, rtol=0.01):
        logger.warning('refined grid has wrong maximum energy')
        return False
    if np.argmax(energy_uniform) not in (1860, 2116):
        logger.warning('uniform grid has maximum energy in wrong location')
        return False
    if np.argmax(energy_refined) not in (1860, 2116):
        logger.warning('refined grid has maximum energy in wrong location')
        return False

    # Check that refined results look like uniform results
    if not np.allclose(energy_uniform, energy_refined, rtol=0.4, atol=0.2):
        logger.warning('refined results differ too much from uniform results')
        return False
    return True
