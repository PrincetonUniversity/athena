# Regression test based on the diffusion of a Gaussian
# scalar distribution.  Convergence of L1 norm of the error
# in r0 is tested.

# Modules
# (needed for global variables modified in run_tests.py, even w/o athena.run(), etc.)
import scripts.utils.athena as athena  # noqa
import scripts.tests.diffusion.scalar_diffusion as scalar_diffusion
import logging

scalar_diffusion.sts_integrators = ['rkl1', 'rkl2']
scalar_diffusion.rate_tols = [-0.99, -1.99]
scalar_diffusion.logger = logging.getLogger('athena' + __name__[7:])


def prepare(*args, **kwargs):
    return scalar_diffusion.prepare('sts', *args, **kwargs)


def run(**kwargs):
    return scalar_diffusion.run(**kwargs)


def analyze():
    return scalar_diffusion.analyze()
