# Regression test for the int12, int13, and int23 outputs in Cartesian coordinates
#
# Runs a blast wave on 2D and 3D meshes of several MeshBlocks with no refinement, with
# static refinement, and with adaptive refinement, writing the 2D integrals of the
# primitive variables together with full-precision formatted table dumps of the same
# data. Checks that the integrals match those computed from the dumps.

# Modules
import logging
import scripts.utils.athena as athena
import scripts.utils.int2d_check as int2d_check

logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module

_input = 'hydro/athinput.blast_int2d'
_mesh_3d = [(-0.5, 0.5, 16), (-0.5, 0.5, 16), (-0.5, 0.5, 16)]
_mesh_2d = [(-0.5, 0.5, 16), (-0.5, 0.5, 16), (-0.5, 0.5, 1)]
_args_2d = ['mesh/nx3=1', 'meshblock/nx3=1']

# (label, extra command-line arguments, root mesh)
_cases = [
    ('3d_uniform', [], _mesh_3d),
    ('3d_smr', ['mesh/refinement=static'], _mesh_3d),
    ('3d_amr', ['mesh/refinement=adaptive'], _mesh_3d),
    ('2d_uniform', _args_2d, _mesh_2d),
    ('2d_smr', _args_2d + ['mesh/refinement=static'], _mesh_2d),
    ('2d_amr', _args_2d + ['mesh/refinement=adaptive'], _mesh_2d),
]


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure(prob='blast', coord='cartesian', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    for label, arguments, _ in _cases:
        athena.run(_input, ['job/problem_id=' + label] + arguments)


# Analyze outputs
def analyze():
    analyze_status = True
    for label, _, mesh in _cases:
        for ext in ('int12', 'int13', 'int23'):
            if not int2d_check.check(label, ext, mesh, 'cartesian', label):
                analyze_status = False
    return analyze_status
