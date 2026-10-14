# Regression test to check that reflection symmetry survives high-order reconstruction

# A spherical blast wave on an even, origin-centered grid has mirror-symmetric initial
# conditions under x->-x, y->-y (and z->-z in 3D), with the symmetry planes on cell
# faces. Next to those faces, and in the exactly uniform ambient medium next to the
# blast edge, neighboring cell averages are bitwise equal, so the branch taken by the
# extremum-preserving PPM limiter (xorder=3,4) sits exactly on its tie. A 1-ulp
# asymmetry (e.g. from FMA contraction, or from a stencil that is not evaluated in a
# mirror-symmetric order) would then send mirror-image cells down different branches
# and grow to a truncation-error-sized asymmetry within a few steps. This test runs the
# blast in 2D and 3D with xorder=3 and 4 and requires the solution to stay symmetric to
# roundoff.

# Modules
import logging
import numpy as np
import scripts.utils.athena as athena
import sys
sys.path.insert(0, '../../vis/python')
import athena_read                             # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module

# (label, extra arguments). The blast has a smooth (log-linear) ramp between r=0.05 and
# r=0.3 and is exactly uniform inside and outside, so neighboring cell averages are
# bitwise equal next to the kinks and across the symmetry planes. Only a few steps are
# taken: a limiter branch flip shows up immediately (~1e-5), while roundoff seeds
# (e.g. FMA contraction) that get amplified at generic, non-degenerate ties need many
# more steps to grow.
d3 = ['mesh/nx1=32', 'mesh/nx2=32', 'mesh/nx3=32',
      'meshblock/nx1=32', 'meshblock/nx2=32', 'meshblock/nx3=32', 'time/cfl_number=0.15']
cases = [('2d_x4_rk4', ['time/xorder=4', 'time/integrator=rk4']),
         ('2d_x3_rk3', ['time/xorder=3', 'time/integrator=rk3']),
         ('3d_x4_rk4', ['time/xorder=4', 'time/integrator=rk4'] + d3),
         ('3d_x3_rk3', ['time/xorder=3', 'time/integrator=rk3'] + d3)]

# largest tolerated asymmetry relative to max|q|: a limiter branch flip gives >= 1e-5
# within a step or two, while roundoff seeds stay below ~1e-12 over these few steps
tolerance = 1.0e-10


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure(
        nghost=4,  # required for fourth-order configurations
        prob='blast',
        coord='cartesian',
        flux='hllc', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    for label, extra in cases:
        arguments = ['job/problem_id=' + label, 'time/ncycle_out=0'] + extra
        athena.run('hydro/athinput.blast_mirror', arguments)


def load(label):
    """Read the last frame into dict name -> array[z, y, x] (one MeshBlock)."""
    filename = 'bin/' + label + '.block0.out1.00001.tab'
    with open(filename) as f:
        f.readline()
        names = f.readline().lstrip('#').split()
    data = np.loadtxt(filename)
    col = {name: data[:, c] for c, name in enumerate(names)}
    axes = []
    for coord in ('x3v', 'x2v', 'x1v'):
        if coord in col:
            axes.append(np.searchsorted(np.unique(col[coord]), col[coord]))
        else:
            axes.append(np.zeros(len(data), dtype=int))
    shape = tuple(a.max() + 1 for a in axes)
    fields = {}
    for name in names:
        if name in ('i', 'j', 'k', 'x1v', 'x2v', 'x3v'):
            continue
        a = np.empty(shape)
        a[axes[0], axes[1], axes[2]] = col[name]
        fields[name] = a
    return fields


# Analyze outputs
def analyze():
    analyze_status = True
    # (axis to flip, velocity component that is odd under that reflection)
    reflections = [(2, 'vel1'), (1, 'vel2'), (0, 'vel3')]
    for label, extra in cases:
        fields = load(label)
        for axis, odd in reflections:
            if fields['rho'].shape[axis] == 1:
                continue
            worst = 0.0
            for name, q in fields.items():
                mirror = np.flip(q, axis=axis)
                asym = np.abs(q + mirror if name == odd else q - mirror).max()
                worst = max(worst, asym / np.abs(q).max())
            logger.info('{}: max relative asymmetry under reflection of axis {} = {:.3e}'
                        .format(label, axis, worst))
            if worst > tolerance:
                logger.warning('{}: reflection symmetry about axis {} is broken ({:.3e})'
                               .format(label, axis, worst))
                analyze_status = False
    return analyze_status
