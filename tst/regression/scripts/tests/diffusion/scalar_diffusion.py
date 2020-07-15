# Regression test based on the diffusion of a Gaussian
# scalar distribution.  Convergence of L1 norm of the error
# in r0 is tested.

# Modules
import logging
import scripts.utils.athena as athena
import numpy as np
import sys
sys.path.insert(0, '../../vis/python')
import athena_read  # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module

_amp = 1.e-6
_nu = 0.25
_t0 = 0.5
_tf = 2.0
_Lx1 = 12.0

resolution_range = [256, 512]
sts_integrators = ['Explicit']  # no STS, explicit integration applied
rate_tols = [-1.99]


def prepare(*args, **kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure(prob='scalar_diff', *args,
                     eos='isothermal', flux='roe',
                     nscalars=1, **kwargs)
    athena.make()


def run(**kwargs):
    for integrator in sts_integrators:
        for n in resolution_range:
            arguments = ['job/problem_id=ScalarDiffusion_'
                         + repr(n) + '_' + integrator,
                         'output2/file_type=tab', 'output2/variable=r0',
                         'output2/data_format=%24.16e', 'output2/dt={}'.format(_tf),
                         'time/cfl_number=0.8',
                         'time/tlim={}'.format(_tf), 'time/nlim=10000',
                         'time/sts_integrator={}'.format(integrator),
                         'time/ncycle_out=0',
                         'mesh/nx1=' + repr(n),
                         'mesh/x1min={}'.format(-_Lx1/2.),
                         'mesh/x1max={}'.format(_Lx1/2.),
                         'mesh/ix1_bc=outflow', 'mesh/ox1_bc=outflow',
                         'mesh/nx2=1', 'mesh/x2min=-1.0', 'mesh/x2max=1.0',
                         'mesh/ix2_bc=periodic', 'mesh/ox2_bc=periodic',
                         'mesh/nx3=1', 'mesh/x3min=-1.0', 'mesh/x3max=1.0',
                         'mesh/ix3_bc=periodic', 'mesh/ox3_bc=periodic',
                         'hydro/iso_sound_speed=1.0',
                         'problem/amp={}'.format(_amp), 'problem/iprob=0',
                         'problem/t0={}'.format(_t0),
                         'problem/nu_scalar_iso={}'.format(_nu)]
            athena.run('hydro/athinput.scalar_diff', arguments)


def analyze():
    l1ERROR = [[] for err in range(0, len(sts_integrators))]
    conv = []

    for i in range(len(sts_integrators)):
        for n in resolution_range:
            x1v, r0 = athena_read.tab('bin/ScalarDiffusion_' + str(n) + '_'
                                      + sts_integrators[i] + '.block0.out2.00001.tab',
                                      raw=True, dimensions=1)
            dx1 = _Lx1/len(x1v)
            analytic = (_amp/np.sqrt(4.*np.pi*_nu*(_t0+_tf))
                        * np.exp(-(x1v**2.)/(4.*_nu*(_t0+_tf))))
            l1ERROR[i].append(sum(np.absolute(r0-analytic)*dx1))

    # estimate L1 convergence
    analyze_status = True
    for i in range(len(sts_integrators)):
        method = sts_integrators[i].upper()
        conv.append(np.diff(np.log(np.array(l1ERROR[i])))
                    / np.diff(np.log(np.array(resolution_range))))
        logger.info('[Scalar Diffusion {}]: Convergence order = {}'
                    .format(method, conv[i]))

        if conv[i] > rate_tols[i]:
            logger.warning('[Scalar Diffusion {}]: '
                           'Scheme NOT converging at expected order.'.format(method))
            analyze_status = False
        else:
            logger.info('[Scalar Diffusion {}]: '
                        'Scheme converging at expected order.'.format(method))

    return analyze_status
