# Regression test based on the diffusion of a Gaussian
# velocity field.  Convergence of L1 norm of the error
# in v is tested.  Expected 1st order conv. for STS.

# Modules
# (needed for global variables modified in run_tests.py, even w/o athena.run(), etc.)
import scripts.utils.athena as athena  # noqa
import scripts.tests.diffusion.viscous_diffusion as viscous_diffusion

viscous_diffusion.method = 'STS'
viscous_diffusion.rate_tols = [-0.99]


def prepare(*args, **kwargs):
    return viscous_diffusion.prepare('sts', *args, **kwargs)


def run(**kwargs):
    return viscous_diffusion.run(**kwargs)


def analyze():
    return viscous_diffusion.analyze()
