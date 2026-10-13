# Regression test to check whether an axisymmetric disk stays axisymmetric to the bit
# when the mesh is split into MeshBlocks

# The Keplerian disk of pgen/disk.cpp is built from r alone, so every phi column follows
# identical arithmetic and the run is axisymmetric to the last bit. Runs it on one block,
# on 16 blocks along phi, and on 4 blocks with x2max written to 12 digits, and requires
# every run to be axisymmetric to the bit and the multi-block runs to reproduce the
# one-block run. Catches a uniform cell width formed per block from its rounded edges.

# Modules
import logging
import glob
import numpy as np
import scripts.utils.athena as athena
import sys
sys.path.insert(0, '../../vis/python')
import athena_read                             # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module

# (label, meshblock nx2, x2max) - nx2 = 256 so 256/nx2 blocks along phi
cases = [('one_block', 256, '6.2831853071795862'),
         ('16_phi_blocks', 16, '6.2831853071795862'),
         ('4_phi_blocks_12digit_2pi', 64, '6.28318530718')]


def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure(prob='disk', coord='cylindrical', **kwargs)
    athena.make()


def run(**kwargs):
    for label, nx2, x2max in cases:
        arguments = ['job/problem_id=' + label,
                     'time/ncycle_out=0', 'time/tlim=3.14',
                     'output1/dt=3.14', 'output2/dt=-1',
                     'meshblock/nx2=' + repr(nx2),
                     'mesh/x2max=' + x2max]
        athena.run('hydro/athinput.disk_cyl_meshblock', arguments)


def load(label):
    """Every MeshBlock's last frame pooled: rows of (x1v, rho, press, vel1, vel2)."""
    files = sorted(glob.glob('bin/' + label + '.block*.out1.00001.tab'))
    data = np.concatenate([np.loadtxt(f) for f in files])
    return files, data


def analyze():
    analyze_status = True
    reference = None
    for label, nx2, x2max in cases:
        files, data = load(label)
        nblocks = len(files)
        x1v = data[:, 1]
        # axisymmetry: at every radius, every variable takes one value over phi
        spread = 0.0
        for col in (4, 5, 6, 7):                # rho, press, vel1, vel2
            for r in np.unique(x1v):
                v = data[x1v == r, col]
                spread = max(spread, (v.max() - v.min()) / max(abs(v).max(), 1e-300))
        logger.info('{}: {} MeshBlock(s), largest relative spread over phi = {:.3e}'
                    .format(label, nblocks, spread))
        if spread != 0.0:
            logger.warning('{}: the run is not axisymmetric to the bit'.format(label))
            analyze_status = False
        # the multi-block runs must reproduce the one-block run
        prof = np.array([[r] + [data[x1v == r, c][0] for c in (4, 5, 6, 7)]
                         for r in np.unique(x1v)])
        if label == 'one_block':
            reference = prof
        elif x2max == cases[0][2] and reference is not None:
            if not np.array_equal(prof, reference):
                diff = np.abs(prof - reference).max()
                logger.warning('{}: differs from the one-block run by up to {:.3e}'
                               .format(label, diff))
                analyze_status = False
    return analyze_status
