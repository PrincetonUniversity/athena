# Regression test based on the diffusion of a Gaussian
# magnetic field.  Convergence of L1 norm of the error
# in b is tested.  Expected 2nd order conv. for explicit.

# Modules
import logging
import scripts.utils.athena as athena
import numpy as np
import sys
sys.path.insert(0, '../../vis/python')
import athena_read  # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module


def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure('b',
                     prob='resist',
                     eos='isothermal')
    athena.make()


def run(**kwargs):
    for i in (64, 128):
        arguments0 = ['output2/file_type=tab', 'output2/variable=bcc2',
                      'output2/data_format=%24.16e', 'output2/dt=0.5',
                      'time/cfl_number=0.8', 'time/tlim=0.5', 'time/nlim=800',
                      'time/ncycle_out=0',
                      'mesh/nx1=' + repr(i), 'mesh/x1min=-4.0', 'mesh/x1max=4.0',
                      'mesh/ix1_bc=outflow', 'mesh/ox1_bc=outflow',
                      'mesh/nx2=1', 'mesh/x2min=-1.0', 'mesh/x2max=1.0',
                      'mesh/ix2_bc=periodic', 'mesh/ox2_bc=periodic',
                      'mesh/nx3=1', 'mesh/x3min=-1.0', 'mesh/x3max=1.0',
                      'mesh/ix3_bc=periodic', 'mesh/ox3_bc=periodic',
                      'hydro/iso_sound_speed=1.0',
                      'problem/amp=1.0e-6', 'problem/iprob=0',
                      'problem/eta_ohm=0.4']
        arguments = arguments0+['job/problem_id=res' + repr(i)]
        athena.run('mhd/athinput.resist', arguments)


def analyze():
    res = [64, 128]
    l1ERROR = []

    for n in res:
        x1v, bcc2 = athena_read.tab("bin/res"+str(n)+".block0.out2.00001.tab", raw=True,
                                    dimensions=1)

        # initial conditions
        amp = 1.e-6
        eta = 0.4
        t0 = 0.5
        sigma = np.sqrt(2.*eta*t0)

        dx1 = 8./len(x1v)
        analytic = ((amp/np.sqrt(2.*np.pi*sigma**2.))
                    * (1./np.sqrt(1.+(2.*eta*t0/sigma**2.)))
                    * np.exp(-(x1v**2.)
                    / (2.*sigma**2.*(1.+(2.*eta*t0/sigma**2.)))))
        l1ERROR.append(sum(np.absolute(bcc2-analytic)*dx1))

    # estimate L1 convergence
    conv = np.diff(np.log(np.array(l1ERROR)))/np.diff(np.log(np.array(res)))
    logger.info('[Resistive Diffusion Explicit]: Convergence order = {}'.format(conv))

    flag = True
    if conv > -1.99:
        logger.warning(
            '[Resistive Diffusion Explicit]: Scheme NOT Converging at ~2nd order.')
        flag = False
    else:
        logger.info('[Resistive Diffusion Explicit]: Scheme Converging at ~2nd order.')

    return flag
