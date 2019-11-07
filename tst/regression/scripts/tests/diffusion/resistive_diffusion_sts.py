# Regression test based on the diffusion of a Gaussian
# magnetic field.  Convergence of L1 norm of the error
# in bcc is tested.

# Modules
# (needed for global variables modified in run_tests.py, even w/o athena.run(), etc.)
import scripts.utils.athena as athena  # noqa
import scripts.tests.diffusion.resistive_diffusion as resistive_diffusion
import logging

resistive_diffusion.sts_integrators = ['rkl1', 'rkl2']
resistive_diffusion.rate_tols = [-0.99, -1.99]
resistive_diffusion.logger = logging.getLogger('athena' + __name__[7:])


def prepare(*args, **kwargs):
    return resistive_diffusion.prepare('sts', *args, **kwargs)


def run(**kwargs):
    return resistive_diffusion.run(**kwargs)


def analyze():
    return resistive_diffusion.analyze()
