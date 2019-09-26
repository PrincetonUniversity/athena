# Regression test based on the diffusion of a Gaussian
# velocity field.  Convergence of L1 norm of the error
# in v is tested.

# Modules
# (needed for global variables modified in run_tests.py, even w/o athena.run(), etc.)
import scripts.utils.athena as athena  # noqa
import scripts.tests.diffusion.viscous_diffusion as viscous_diffusion
import logging

viscous_diffusion.sts_integrators = ['rkl1', 'rkl2']
viscous_diffusion.rate_tols = [-0.99, -1.99]
viscous_diffusion.logger = logging.getLogger('athena' + __name__[7:])


def prepare(*args, **kwargs):
    return viscous_diffusion.prepare('sts', *args, **kwargs)


def run(**kwargs):
    return viscous_diffusion.run(**kwargs)


def analyze():
    return viscous_diffusion.analyze()
