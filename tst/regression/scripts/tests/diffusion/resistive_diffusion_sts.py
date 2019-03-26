# Regression test based on the diffusion of a Gaussian
# magnetic field.  Convergence of L1 norm of the error
# in b is tested.  Expected 1st order conv. for STS.

# Modules
# (needed for global variables modified in run_tests.py, even w/o athena.run(), etc.)
import scripts.utils.athena as athena # noqa
import scripts.tests.diffusion.resistive_diffusion as resistive_diffusion

resistive_diffusion.method = 'STS'
resistive_diffusion.rate_tols = [-0.99]


def prepare(*args, **kwargs):
    return resistive_diffusion.prepare('sts', *args, **kwargs)


def run(**kwargs):
    return resistive_diffusion.run(**kwargs)


def analyze():
    return resistive_diffusion.analyze()
