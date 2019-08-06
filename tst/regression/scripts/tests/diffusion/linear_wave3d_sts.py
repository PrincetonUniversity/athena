# Regression test based on the decaying linear wave due to viscosity,
# Ohmic resistivity and thermal conduction. The decay rate is fit and
# then compared with analytic solution.  This test employs STS.

# Modules
# (needed for global variables modified in run_tests.py, even w/o athena.run(), etc.)
import scripts.utils.athena as athena  # noqa
import scripts.tests.diffusion.linear_wave3d as linear_wave3d
import logging

linear_wave3d.method = 'STS'
# Override analyze() paramaters from non-STS diffusion module
# linear_wave3d.error_rel_tols = [0.22, 0.05]

# lower bound on convergence rate at final (Nx1=64) asymptotic convergence regime
linear_wave3d.rate_tols = [1.0]
linear_wave3d.logger = logging.getLogger('athena' + __name__[7:])


def prepare(*args, **kwargs):
    linear_wave3d.prepare('sts', *args, **kwargs)


def run(**kwargs):
    return linear_wave3d.run(**kwargs)


def analyze():
    return linear_wave3d.analyze()
