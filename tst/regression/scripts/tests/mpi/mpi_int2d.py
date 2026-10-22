# Regression test for the int12, int13, and int23 outputs with MPI
#
# Runs the spherical-polar blast wave with AMR from outputs/int2d_curvilinear on four
# MPI ranks, so that the integrals are summed over MeshBlocks owned by different ranks,
# and checks them against integrals computed from full-precision formatted table dumps
# of the same data.

# Modules
import logging
import scripts.utils.athena as athena
import scripts.utils.int2d_check as int2d_check
import scripts.tests.outputs.int2d_curvilinear as int2d_curvilinear

logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module

_coord = 'spherical_polar'
_label = 'mpi_' + _coord + '_amr'
_mesh, _arguments = int2d_curvilinear._coords[_coord]


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure('mpi', prob='blast', coord=_coord, **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    athena.mpirun(kwargs['mpirun_cmd'], kwargs['mpirun_opts'], 4,
                  int2d_curvilinear._input,
                  ['job/problem_id=' + _label]
                  + int2d_curvilinear.mesh_arguments(_mesh) + _arguments
                  + int2d_curvilinear.amr_arguments)


# Analyze outputs
def analyze():
    analyze_status = True
    for ext in ('int12', 'int13', 'int23'):
        if not int2d_check.check(_label, ext, _mesh, _coord, _label):
            analyze_status = False
    return analyze_status
