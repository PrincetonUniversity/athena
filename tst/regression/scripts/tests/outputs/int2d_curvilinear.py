# Regression test for the int12, int13, and int23 outputs in curvilinear coordinates
#
# Runs a blast wave in cylindrical and spherical-polar coordinates on 3D meshes of
# several MeshBlocks, with and without static refinement, writing the 2D integrals of
# the primitive variables together with full-precision formatted table dumps of the
# same data. Checks that the integrals match those computed from the dumps with the
# exact areas of the coordinate surfaces.

# Modules
import logging
import numpy as np
import scripts.utils.athena as athena
import scripts.utils.int2d_check as int2d_check

logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module

_input = 'hydro/athinput.blast_int2d'

# coordinate system: (root mesh, command-line arguments for the mesh and the blast)
_coords = {
    'cylindrical': (
        [(0.5, 1.5, 16), (0.0, 0.5 * np.pi, 16), (-0.5, 0.5, 16)],
        ['problem/x1_0=1.0', 'problem/x2_0=0.8', 'problem/x3_0=0.05',
         'problem/radius=0.3',
         'refinement1/x1min=1.0', 'refinement1/x1max=1.2',
         'refinement1/x2min=0.0', 'refinement1/x2max=0.3',
         'refinement1/x3min=0.0', 'refinement1/x3max=0.2']),
    'spherical_polar': (
        [(1.0, 2.0, 16), (0.25 * np.pi, 0.75 * np.pi, 16),
         (-0.25 * np.pi, 0.25 * np.pi, 16)],
        ['problem/x1_0=1.5', 'problem/x2_0=1.5', 'problem/x3_0=0.1',
         'problem/radius=0.3',
         'refinement1/x1min=1.5', 'refinement1/x1max=1.7',
         'refinement1/x2min=1.6', 'refinement1/x2max=1.9',
         'refinement1/x3min=0.0', 'refinement1/x3max=0.3']),
}

# (label suffix, extra command-line arguments)
_refinements = [('uniform', []), ('smr', ['mesh/refinement=static'])]

# configure options passed to prepare(), reused when building each coordinate system
_config_kwargs = {}


def mesh_arguments(mesh):
    arguments = []
    for axis, (xmin, xmax, nx) in enumerate(mesh):
        arguments += ['mesh/x{}min={}'.format(axis + 1, repr(xmin)),
                      'mesh/x{}max={}'.format(axis + 1, repr(xmax)),
                      'mesh/nx{}={}'.format(axis + 1, nx)]
    return arguments


# Prepare Athena++ (each coordinate system is built in run())
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    _config_kwargs.update(kwargs)


# Run Athena++
def run(**kwargs):
    for coord, (mesh, arguments) in _coords.items():
        athena.configure(prob='blast', coord=coord, **_config_kwargs)
        athena.make()
        for suffix, refinement in _refinements:
            label = coord + '_' + suffix
            athena.run(_input, ['job/problem_id=' + label] + mesh_arguments(mesh)
                       + arguments + refinement)


# Analyze outputs
def analyze():
    analyze_status = True
    for coord, (mesh, _) in _coords.items():
        for suffix, _ in _refinements:
            label = coord + '_' + suffix
            for ext in ('int12', 'int13', 'int23'):
                if not int2d_check.check(label, ext, mesh, coord, label):
                    analyze_status = False
    return analyze_status
