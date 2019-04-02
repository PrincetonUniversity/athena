# Regression test based on the decaying linear wave due to viscosity,
# Ohmic resistivity and thermal conduction. The decay rate is fit and
# then compared with analytic solution

# Modules
import logging
import numpy as np
import sys
import scripts.utils.athena as athena
sys.path.insert(0, '../../vis/python')
import athena_read                             # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module


def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure('b',
                     prob='linear_wave',
                     flux='hlld',
                     eos='adiabatic')
    athena.make()


def run(**kwargs):
    arguments0 = ['output2/dt=0.03',
                  'time/tlim=3.0',
                  'time/ncycle_out=0',
                  # L-going slow wave
                  'problem/wave_flag=2',
                  'problem/amp=1.0e-4', 'problem/vflow=0.0',
                  'problem/nu_iso=0.01', 'problem/eta_ohm=0.02', 'problem/kappa_iso=0.02']
    arguments = arguments0+['job/problem_id=DecayLinWave']
    athena.run('mhd/athinput.linear_wave3d', arguments)


def analyze():
    ksqr = (2.0*np.pi)**2
    # fast mode decay rate = (19\nu/4+3\eta+3(\gamma-1)^2*kappa/gamma/4)*(2/15)*k^2
    # slow mode decay rate = (4\nu/+3\eta/4+3(\gamma-1)^2*kappa/gamma/4)*(2/15)*k^2
    rate = 2.0*(4.0*0.01+3.0*0.02/4.0+0.02*4.0/5.0)/15.0*ksqr

    basename = 'bin/DecayLinWave.block0.out2.'
    nframe = 100
    dumprate = 0.03
    max_vy = np.zeros(nframe)
    tt = np.zeros(nframe)
    for i in range(nframe):
        x1f, x2f, x3f, data = athena_read.vtk(basename+str(i).zfill(5)+'.vtk')
        max_vy[i] = np.max(data['vel'][..., 1])
        tt[i] = i*dumprate

    # estimate the decay rate from simulation
    new_rate, coeff = np.polyfit(tt, np.log(np.abs(max_vy)), 1, w=np.sqrt(np.abs(max_vy)))
    new_rate = -new_rate
    print('[Decaying Linear Wave-3D]: Ref(decay_rate) = {}'.format(rate))
    print('[Decaying Linear Wave-3D]: New(decay_rate) = {}'.format(new_rate))

    flag = True
    error_rel = np.fabs(rate/new_rate-1.0)
    if error_rel > 0.1:
        print('[Decaying Linear Wave-3D]: decay rate is off by 10 percent')
        flag = False
    else:
        print('[Decaying Linear Wave-3D]: decay rate is within 10 percent precision')

    return flag
