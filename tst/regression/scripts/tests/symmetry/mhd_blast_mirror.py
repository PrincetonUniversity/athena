# Regression test to check that reflection symmetry survives constrained transport

# An MHD blast wave in a uniform field on an even, origin-centered grid has
# mirror-symmetric initial conditions under x->-x, y->-y (and z->-z in 3D), with the
# symmetry planes on cell faces; with the field along a diagonal the 2D problem is
# instead symmetric under the 180-degree rotation (x,y)->(-x,-y). Only operations that
# are evaluated in a reflection-symmetric order keep that symmetry to the bit: the
# upwind weighting of the corner EMFs, the order in which the eight contributions to
# each corner EMF are summed, and the HLLD solver's choice of flux on a stationary
# contact all used to break it at the 1-ulp level (#165), and the extremum-preserving
# PPM limiter then amplifies the seed to a truncation-error-sized asymmetry within a
# few tens of steps. This test runs the blast with HLLE, HLLD and Roe, with PLM and with
# extremum-preserving PPM, in 2D and 3D, and requires the solution to stay bit-exactly
# symmetric. FMA contraction is disabled for the build because it seeds the same kind
# of asymmetry on its own (a*b+c is contracted or not depending on operand order).

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

_fluxes = ['hlle', 'hlld', 'roe']
_exec = os.path.join('bin', 'athena')

# (label, symmetry, extra arguments). PLM keeps a 1-ulp seed at the roundoff level for
# hundreds of steps, so the PLM cases mostly guard the bit-exactness itself; the PPM
# cases are the ones that grow an asymmetric corner EMF to O(1e-3) by these times.
d3 = ['mesh/nx1=32', 'mesh/nx2=32', 'mesh/nx3=32',
      'meshblock/nx1=32', 'meshblock/nx2=32', 'meshblock/nx3=32', 'time/cfl_number=0.15']
cases = [('2d_x2_vl2', 'mirror', ['time/xorder=2', 'time/integrator=vl2',
                                  'time/tlim=0.5', 'output1/dt=0.5']),
         ('2d_x3_rk3', 'mirror', ['time/xorder=3', 'time/integrator=rk3',
                                  'time/tlim=0.5', 'output1/dt=0.5']),
         ('3d_x2_vl2', 'mirror', ['time/xorder=2', 'time/integrator=vl2',
                                  'time/tlim=0.3', 'output1/dt=0.3'] + d3),
         ('3d_x3_rk3', 'mirror', ['time/xorder=3', 'time/integrator=rk3',
                                  'time/tlim=0.3', 'output1/dt=0.3'] + d3),
         ('2d_x3_rk3_diag', 'point', ['time/xorder=3', 'time/integrator=rk3',
                                      'time/tlim=0.5', 'output1/dt=0.5',
                                      'problem/angle=45.0'])]

# largest tolerated asymmetry relative to max|q|: the corner EMF is now evaluated in a
# reflection-symmetric order, so the asymmetry is exactly zero without FMA contraction;
# a single reordered operation shows up as ~1e-14 within a few steps and grows from
# there with PPM
tolerance = 1.0e-14


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    global _fluxes
    for i in athena.global_config_args:
        tmp = i.split('=')
        if tmp[0] == '--flux' and len(tmp) == 2:
            _fluxes = [tmp[1]]
    for flux in _fluxes:
        athena.configure('b',
                         nghost=3,
                         prob='blast',
                         coord='cartesian',
                         flux=flux,
                         cflag='-ffp-contract=off',
                         **kwargs)
        # to save time, reuse compiled .o files for all executables created in this test:
        athena.make(clean_first=False)
        move(_exec, _exec + '_' + flux)
        os.system('cp -r obj obj_' + flux)
    os.system('rm -rf obj')


# Run Athena++
def run(**kwargs):
    for flux in _fluxes:
        move(_exec + '_' + flux, _exec)
        os.system('mv obj_' + flux + ' obj')
        for label, symmetry, extra in cases:
            arguments = ['job/problem_id=' + flux + '_' + label,
                         'time/ncycle_out=0'] + extra
            athena.run('mhd/athinput.blast_mirror', arguments)
        os.system('rm -rf obj')


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
    for flux in _fluxes:
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
