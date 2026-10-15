# Helpers for the int12/int13/int23 output regression tests
#
# Computes the expected 2D integrals independently from full-precision formatted table
# (.tab) dumps of the same data, one file per MeshBlock, and compares them with the
# values that athena_read.int2d() reads from the integral output files.
#
# For each cell, the integrand is multiplied by the area of the coordinate surface
# through the cell center (normal to the third dimension) restricted to the cell, and
# the product is added to every bin of the third dimension that the cell covers.
# Cell faces are reconstructed from the uniform mesh spacing at the level of the
# MeshBlock, so the meshes used with these helpers must be uniform in every direction.

# Modules
import glob
import logging
import re
import sys
import numpy as np
sys.path.insert(0, '../../vis/python')
import athena_read  # noqa

logger = logging.getLogger('athena.int2d_check')

# third dimension (0-based axis) integrated out by each file type
_third_axis = {'int12': 2, 'int13': 1, 'int23': 0}
_rtol = 1.0e-12


def read_tab(filename):
    """Return (time, dict of columns) from a .tab file."""
    with open(filename, 'r') as f:
        time = float(re.search(r'time=(\S+)', f.readline()).group(1))
        names = f.readline().split()[1:]
    data = np.atleast_2d(np.loadtxt(filename))
    return time, {name: data[:, n] for n, name in enumerate(names)}


def cell_geometry(columns, mesh):
    """Return per-cell centers and faces along each axis for one MeshBlock.

    mesh is a list of (xmin, xmax, nx) of the root grid for x1, x2, x3."""
    ncell = len(columns['rho'])
    centers, lower, upper = [], [], []
    for axis, (xmin, xmax, nx) in enumerate(mesh):
        name = 'x{}v'.format(axis + 1)
        if nx == 1:
            centers.append(np.full(ncell, 0.5 * (xmin + xmax)))
            lower.append(np.full(ncell, xmin))
            upper.append(np.full(ncell, xmax))
            continue
        x = columns[name]
        spacing = np.median(np.diff(np.unique(x)))
        dx_root = (xmax - xmin) / nx
        level = int(round(np.log2(dx_root / spacing)))
        dx = dx_root / 2**level
        index = np.floor((x - xmin) / dx)
        centers.append(x)
        lower.append(xmin + index * dx)
        upper.append(xmin + (index + 1) * dx)
    return centers, lower, upper


def cell_area(coord, third, centers, lower, upper):
    """Area of the coordinate surface through each cell center, normal to axis third."""
    width = [u - lo for lo, u in zip(lower, upper)]
    if coord == 'cartesian':
        others = [a for a in range(3) if a != third]
        return width[others[0]] * width[others[1]]
    half_r2 = 0.5 * (upper[0]**2 - lower[0]**2)
    if coord == 'cylindrical':
        if third == 2:
            return half_r2 * width[1]
        if third == 1:
            return width[0] * width[2]
        return centers[0] * width[1] * width[2]
    if coord == 'spherical_polar':
        if third == 2:
            return half_r2 * width[1]
        if third == 1:
            return np.sin(centers[1]) * half_r2 * width[2]
        return centers[0]**2 * (np.cos(lower[1]) - np.cos(upper[1])) * width[2]
    raise ValueError('unknown coordinate system ' + coord)


def expected_integrals(tab_files, mesh, coord, ext, xf, varnames):
    """Integrals over all MeshBlocks, binned on the faces xf of the third dimension."""
    third = _third_axis[ext]
    bin_centers = 0.5 * (xf[:-1] + xf[1:])
    result = {name: np.zeros(len(bin_centers)) for name in varnames}
    for filename in tab_files:
        _, columns = read_tab(filename)
        centers, lower, upper = cell_geometry(columns, mesh)
        area = cell_area(coord, third, centers, lower, upper)
        first = np.searchsorted(bin_centers, lower[third])
        last = np.searchsorted(bin_centers, upper[third])
        for name in varnames:
            weighted = columns[name] * area
            for n in range(len(weighted)):
                result[name][first[n]:last[n]] += weighted[n]
    return result


def check(basename, ext, mesh, coord, label):
    """Compare basename.ext with integrals computed from basename.block*.out1.*.tab."""
    status = True
    out = athena_read.int2d('bin/{}.{}'.format(basename, ext))
    varnames = list(out._fields[2:])
    time = np.atleast_1d(out.time)
    xf = np.asarray(out.xf)

    # bin faces must lie on the uniform mesh
    xmin, xmax, _ = mesh[_third_axis[ext]]
    nbins = len(xf) - 1
    faces = xmin + (xmax - xmin) * np.arange(nbins + 1) / nbins
    if not np.allclose(xf, faces, rtol=_rtol, atol=_rtol * abs(xmax - xmin)):
        logger.warning('%s %s: bin faces do not match the mesh', label, ext)
        status = False

    numbers = sorted({re.search(r'out1\.(\d+)\.tab', f).group(1) for f in
                      glob.glob('bin/{}.block*.out1.*.tab'.format(basename))})
    if len(numbers) != len(time):
        logger.warning('%s %s: %d integral records but %d tab dumps', label, ext,
                       len(time), len(numbers))
        return False

    for n, number in enumerate(numbers):
        pattern = 'bin/{}.block*.out1.{}.tab'.format(basename, number)
        tab_files = sorted(glob.glob(pattern))
        tab_time, _ = read_tab(tab_files[0])
        if abs(time[n] - tab_time) > 1.0e-5 * max(abs(tab_time), 1.0e-5):
            logger.warning('%s %s: record %d time %g, tab dump time %g', label, ext, n,
                           time[n], tab_time)
            status = False
        ref = expected_integrals(tab_files, mesh, coord, ext, xf, varnames)
        for name in varnames:
            got = np.atleast_2d(getattr(out, name))[n]
            scale = max(np.abs(ref[name]).max(), 1.0e-300)
            err = np.abs(got - ref[name]).max() / scale
            if err > _rtol:
                logger.warning('%s %s: %s at t=%g differs by %.3e (relative)',
                               label, ext, name, time[n], err)
                status = False
    return status
