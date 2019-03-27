# Regression test based on the diffusion of a Gaussian
# magnetic field.  Convergence of L1 norm of the error
# in b is tested.  Expected 2nd order conv. for explicit.

# Modules
import scripts.utils.athena as athena
import numpy as np
import sys
sys.path.insert(0, '../../vis/python')
import athena_read  # noqa
athena_read.check_nan_flag = True

_amp = 1.e-6
_eta = 0.4
_t0 = 0.5
_Lx1 = 8.0

resolution_range = [64, 128]
method = 'Explicit'
rate_tols = [-1.99]


def prepare(*args, **kwargs):
    athena.configure('b', *args,
                     prob='resist',
                     eos='isothermal', **kwargs)
    athena.make()


def run(**kwargs):
    for i in resolution_range:
        arguments0 = ['output2/file_type=tab', 'output2/variable=bcc2',
                      'output2/data_format=%24.16e', 'output2/dt={}'.format(_t0),
                      'time/cfl_number=0.8',
                      'time/tlim={}'.format(_t0), 'time/nlim=800',
                      'time/ncycle_out=0',
                      'mesh/nx1=' + repr(i),
                      'mesh/x1min={}'.format(-_Lx1/2.),
                      'mesh/x1max={}'.format(_Lx1/2.),
                      'mesh/ix1_bc=outflow', 'mesh/ox1_bc=outflow',
                      'mesh/nx2=1', 'mesh/x2min=-1.0', 'mesh/x2max=1.0',
                      'mesh/ix2_bc=periodic', 'mesh/ox2_bc=periodic',
                      'mesh/nx3=1', 'mesh/x3min=-1.0', 'mesh/x3max=1.0',
                      'mesh/ix3_bc=periodic', 'mesh/ox3_bc=periodic',
                      'hydro/iso_sound_speed=1.0',
                      'problem/amp={}'.format(_amp), 'problem/iprob=0',
                      'problem/eta_ohm={}'.format(_eta)]
        arguments = arguments0+['job/problem_id=res' + repr(i)]
        athena.run('mhd/athinput.resist', arguments)


def analyze():
    l1ERROR = []

    for n in resolution_range:
        x1v, bcc2 = athena_read.tab("bin/res"+str(n)+".block0.out2.00001.tab", raw=True,
                                    dimensions=1)
        sigma = np.sqrt(2.*_eta*_t0)
        dx1 = _Lx1/len(x1v)
        analytic = ((_amp/np.sqrt(2.*np.pi*sigma**2.))
                    * (1./np.sqrt(1.+(2.*_eta*_t0/sigma**2.)))
                    * np.exp(-(x1v**2.)
                    / (2.*sigma**2.*(1.+(2.*_eta*_t0/sigma**2.)))))
        l1ERROR.append(sum(np.absolute(bcc2-analytic)*dx1))

    # estimate L1 convergence
    conv = np.diff(np.log(np.array(l1ERROR)))/np.diff(np.log(np.array(resolution_range)))
    print('[Resistive Diffusion {}]: Convergence order = {}'.format(method, conv))

    analyze_status = True
    if conv > rate_tols[-1]:
        print('[Resistive Diffusion {}]: '
              'Scheme NOT converging at expected order.'.format(method))
        analyze_status = False
    else:
        print('[Resistive Diffusion {}]: '
              'Scheme converging at expected order.'.format(method))

    return analyze_status
