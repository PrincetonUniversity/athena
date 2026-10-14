# Regression test to check that reflection symmetry survives the SR Riemann solvers

# The special-relativistic version of symmetry/mhd_blast_mirror: a blast wave at rest in
# a uniform field on an even, origin-centered grid is mirror symmetric under x->-x, y->-y
# (and z->-z in 3D), and with the field along a diagonal the 2D problem is symmetric
# under the 180-degree rotation (x,y)->(-x,-y). The SR HLLE/HLLC/LLF solvers and the
# relativistic wave-speed functions used to break the symmetry at the 1-ulp level (#165):
# the HLL state summed its flux difference in a left-to-right order, the magnetosonic
# speeds took lambda_plus and -lambda_minus from different root formulas when the normal
# velocity was exactly zero, and HLLC always took the L* flux on a stationary contact.
# The extremum-preserving PPM limiter amplifies such seeds to O(1e-3) within a few tens
# of steps. This test runs the blast with HLLE and LLF (MHD) and HLLC (hydro), with PLM
# and with extremum-preserving PPM, in 2D and 3D, and requires the solution to stay
# bit-exactly symmetric. FMA contraction is disabled for the build because it seeds the
# same kind of asymmetry on its own.

# Modules
import logging
import numpy as np
import os
import scripts.utils.athena as athena
import sys
from shutil import move
sys.path.insert(0, '../../vis/python')
import athena_read                             # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module

# (Riemann solver, whether to enable magnetic fields)
_configs = [('hlle', True), ('llf', True), ('hllc', False)]
_exec = os.path.join('bin', 'athena')

# (label, symmetry, extra arguments). PLM keeps a 1-ulp seed at the roundoff level for
# hundreds of steps, so the PLM cases mostly guard the bit-exactness itself; the PPM
# cases are the ones that grow a seed to O(1e-3) by these times.
d3 = ['mesh/nx1=32', 'mesh/nx2=32', 'mesh/nx3=32',
      'meshblock/nx1=32', 'meshblock/nx2=32', 'meshblock/nx3=32', 'time/cfl_number=0.15',
      'time/tlim=0.15', 'output1/dt=0.15']
cases = [('2d_x2_vl2', 'mirror', ['time/xorder=2', 'time/integrator=vl2']),
         ('2d_x3_rk3', 'mirror', ['time/xorder=3', 'time/integrator=rk3']),
         ('3d_x2_vl2', 'mirror', ['time/xorder=2', 'time/integrator=vl2'] + d3),
         ('3d_x3_rk3', 'mirror', ['time/xorder=3', 'time/integrator=rk3'] + d3),
         ('2d_x3_rk3_diag', 'point', ['time/xorder=3', 'time/integrator=rk3',
                                      'problem/angle=45.0'])]

# largest tolerated asymmetry relative to max|q|: the asymmetry is exactly zero without
# FMA contraction; a single reordered operation shows up as ~1e-14 within a few steps
# and grows from there with PPM
tolerance = 1.0e-14


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    for flux, mhd in _configs:
        args = ['s', 'b'] if mhd else ['s']
        athena.configure(*args,
                         nghost=3,
                         prob='blast',
                         coord='cartesian',
                         flux=flux,
                         cflag='-ffp-contract=off',
                         **kwargs)
        # the hydro and MHD builds cannot share object files, so build each from scratch
        athena.make()
        move(_exec, _exec + '_' + flux)


# Run Athena++
def run(**kwargs):
    for flux, mhd in _configs:
        move(_exec + '_' + flux, _exec)
        for label, symmetry, extra in cases:
            arguments = ['job/problem_id=' + flux + '_' + label,
                         'time/ncycle_out=0'] + extra
            athena.run('mhd_sr/athinput.blast_mirror', arguments)


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
    # Fields that are odd under each transformation, using the invariance of the MHD
    # equations under B->-B to keep the initial B0 x-hat (or diagonal) field fixed:
    # (axes to flip, odd fields). B is an axial vector, so under a mirror reflection the
    # component normal to the mirror is even and the in-plane components are odd; the
    # y and z mirrors then need B->-B to keep B0 x-hat, which flips the parity of Bcc1.
    transformations = {'mirror': [((2,), ['vel1', 'Bcc2', 'Bcc3']),
                                  ((1,), ['vel2', 'Bcc2']),
                                  ((0,), ['vel3', 'Bcc3'])],
                       'point': [((1, 2), ['vel1', 'vel2', 'Bcc3'])]}
    for flux, mhd in _configs:
        for case_label, symmetry, extra in cases:
            label = flux + '_' + case_label
            fields = load(label)
            for axes, odd in transformations[symmetry]:
                if any(fields['rho'].shape[axis] == 1 for axis in axes):
                    continue
                worst = 0.0
                for name, q in fields.items():
                    mirror = np.flip(q, axis=axes)
                    asym = np.abs(q + mirror if name in odd else q - mirror).max()
                    scale = np.abs(q).max()
                    worst = max(worst, asym / scale if scale > 0.0 else asym)
                logger.info('{}: max relative asymmetry under flipping axes {} = {:.3e}'
                            .format(label, axes, worst))
                if worst > tolerance:
                    logger.warning('{}: symmetry under flipping axes {} is broken '
                                   '({:.3e})'.format(label, axes, worst))
                    analyze_status = False
    return analyze_status
