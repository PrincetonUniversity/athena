#!/usr/bin/env python3
"""Analysis tools for the MHD Kelvin-Helmholtz problem of thesis section 5.3.3.

Completes the quantitative analysis that the thesis lists as future work for the
MHD extension of the Lecoanet et al. (2016) KHI comparison:
  1. dye entropy production S(t) = -Integral( c*ln(c) ) dV  (thesis eq. 5.22)
  2. L2 error of the dye concentration field (and |B|) relative to a reference
     solution at a higher resolution and/or from a higher-accuracy solver,
     following the thesis L2(D_ref, c) metric (thesis eq. 5.23)

Requires HDF5 (.athdf) outputs of the primitive variables and the passive scalar,
e.g. from inputs/mhd/athinput.kh-shear-lecoanet with -hdf5 -h5double.

Example usage for a resolution study (from the run directory):
    # dye entropy time series for one run:
    python khi_mhd_analysis.py entropy 'kh-shear-lecoanet.out2.*.athdf'

    # L2 error of coarse run vs reference run at matched output time:
    python khi_mhd_analysis.py l2 coarse/kh-shear-lecoanet.out2.00120.athdf \
                                  ref/kh-shear-lecoanet.out2.00120.athdf
"""

import argparse
import glob
import sys

import numpy as np

import athena_read


def _dye_concentration(data):
    """Return the dye concentration c from an athdf dict (specific scalar r0)."""
    if 'r0' in data:
        return data['r0']
    if 's0' in data:  # conserved dye density -> specific concentration
        return data['s0']/data['rho']
    raise RuntimeError('no passive scalar (r0/s0) found in output; '
                       'configure with --nscalars=1 and output prim variables')


def dye_entropy(filename):
    """Compute S = -Int c*ln(c) dV for a single .athdf output (thesis eq. 5.22)."""
    data = athena_read.athdf(filename)
    c = np.clip(_dye_concentration(data), 1e-300, None)
    # uniform Cartesian mesh assumed:
    dx = data['x1f'][1] - data['x1f'][0]
    dy = data['x2f'][1] - data['x2f'][0]
    dv = dx*dy
    if len(data['x3f']) > 2:
        dv *= data['x3f'][1] - data['x3f'][0]
    return data['Time'], -np.sum(c*np.log(c))*dv


def block_average(q, factor):
    """Volume average a 2D/3D array over factor^ndim blocks (restriction to the
    coarse mesh), for comparing solutions at different resolutions."""
    if factor == 1:
        return q
    if q.ndim == 3 and q.shape[0] == 1:  # 2D data stored as (1, ny, nx)
        q = q[0]
    if q.ndim == 2:
        ny, nx = q.shape
        return q.reshape(ny//factor, factor, nx//factor, factor).mean(axis=(1, 3))
    nz, ny, nx = q.shape
    return q.reshape(nz//factor, factor, ny//factor, factor,
                     nx//factor, factor).mean(axis=(1, 3, 5))


def l2_error(file_coarse, file_ref, quantity='dye'):
    """L2(ref, q) = sqrt( Int (q - q_ref)^2 dV ) with the reference solution
    restricted (volume-averaged) onto the coarse mesh (thesis eq. 5.23)."""
    dc = athena_read.athdf(file_coarse)
    dr = athena_read.athdf(file_ref)
    if abs(dc['Time'] - dr['Time']) > 1e-10:
        print('WARNING: comparing outputs at different times: '
              f"{dc['Time']} vs {dr['Time']}", file=sys.stderr)

    if quantity == 'dye':
        qc, qr = _dye_concentration(dc), _dye_concentration(dr)
    elif quantity == 'B':
        qc = np.sqrt(dc['Bcc1']**2 + dc['Bcc2']**2 + dc['Bcc3']**2)
        qr = np.sqrt(dr['Bcc1']**2 + dr['Bcc2']**2 + dr['Bcc3']**2)
    else:
        qc, qr = dc[quantity], dr[quantity]

    nxc = len(dc['x1v'])
    nxr = len(dr['x1v'])
    if nxr % nxc != 0:
        raise RuntimeError(f'reference nx1={nxr} is not a multiple of coarse '
                           f'nx1={nxc}')
    qr_restricted = block_average(np.squeeze(qr), nxr//nxc)
    qc = np.squeeze(qc)

    dx = dc['x1f'][1] - dc['x1f'][0]
    dy = dc['x2f'][1] - dc['x2f'][0]
    dv = dx*dy
    if len(dc['x3f']) > 2:
        dv *= dc['x3f'][1] - dc['x3f'][0]
    return np.sqrt(np.sum((qc - qr_restricted)**2)*dv)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='mode', required=True)

    p_ent = sub.add_parser('entropy', help='dye entropy time series')
    p_ent.add_argument('pattern', help='glob pattern of .athdf outputs')

    p_l2 = sub.add_parser('l2', help='L2 error vs reference solution')
    p_l2.add_argument('coarse', help='.athdf output of the run under test')
    p_l2.add_argument('reference', help='.athdf output of the reference run')
    p_l2.add_argument('--quantity', default='dye',
                      help="'dye' (default), 'B', or any athdf dataset name")

    args = parser.parse_args()
    if args.mode == 'entropy':
        files = sorted(glob.glob(args.pattern))
        if not files:
            sys.exit(f'no files match {args.pattern}')
        print('# time  dye_entropy')
        for f in files:
            t, s = dye_entropy(f)
            print(f'{t:.6e}  {s:.10e}')
    else:
        err = l2_error(args.coarse, args.reference, args.quantity)
        print(f'L2({args.quantity}) = {err:.10e}')


if __name__ == '__main__':
    main()
