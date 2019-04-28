# Regression test based on the decaying linear wave due to thermal
# conduction. The decay rate is fit and then compared with analytic
# solution.  This test employs STS.

# Modules
# (needed for global variables modified in run_tests.py, even w/o athena.run(), etc.)
import scripts.utils.athena as athena  # noqa
import scripts.tests.diffusion.thermal_attenuation as thermal_attenuation
import logging

thermal_attenuation.method = 'STS'
thermal_attenuation.logger = logging.getLogger('athena' + __name__[7:])


def prepare(*args, **kwargs):
    thermal_attenuation.prepare('sts', *args, **kwargs)


def run(**kwargs):
    return thermal_attenuation.run(**kwargs)


def analyze():
    return thermal_attenuation.analyze()
