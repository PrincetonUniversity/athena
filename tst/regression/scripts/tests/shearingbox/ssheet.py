# Regression test of shearing box and orbital advection with 2d hydro shwave.

# Modules
import logging
import cmath
import math
import mpmath as mp
import scripts.utils.athena as athena
import sys
sys.path.insert(0, '../../vis/python')
import athena_read  # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure(prob='ssheet', flux='hlle',
                     eos='isothermal', **kwargs)
    athena.make()


# Run Athena++ w/wo Orbital Advection
def run(**kwargs):
    # HD shwave w/o Orbital Advection
    arguments = [
        'job/problem_id=SSHEET_SHWAVE',
        'output1/file_type=hst', 'output1/dt=10.0',
        'output2/file_type=vtk', 'output2/variable=prim',
        'output2/dt=-1', 'time/cfl_number=0.4',
        'time/tlim=8000.0', 'time/nlim=2000',
        'time/xorder=2', 'time/integrator=vl2',
        'mesh/nx1=64', 'mesh/x1min=-2.0', 'mesh/x1max=2.0',
        'mesh/ix1_bc=shear_periodic', 'mesh/ox1_bc=shear_periodic',
        'mesh/nx2=64', 'mesh/x2min=-2.0', 'mesh/x2max=2.0',
        'mesh/ix2_bc=periodic', 'mesh/ox2_bc=periodic',
        'mesh/nx3=1', 'mesh/x3min=-0.5', 'mesh/x3max=0.5',
        'hydro/iso_sound_speed=0.001', 'problem/ipert=3',
        'problem/amp=4.0e-4', 'problem/nwx=-4', 'problem/nwy=1',
        'orbital_advection/Omega0=1.0e-3', 'orbital_advection/qshear=1.5',
        'orbital_advection/shboxcoord=1', 'orbital_advection/OAorder=0',
        'problem/error_output=true', 'time/ncycle_out=0']
    athena.run('hydro/athinput.ssheet', arguments)

    # HD shwave w/  Orbital Advection
    arguments = [
        'job/problem_id=SSHEET_SHWAVE_ORB',
        'output1/file_type=hst', 'output1/dt=10.0',
        'output2/file_type=vtk', 'output2/variable=prim',
        'output2/dt=-1', 'time/cfl_number=0.4',
        'time/tlim=8000.0', 'time/nlim=2000',
        'time/xorder=2', 'time/integrator=vl2',
        'mesh/nx1=64', 'mesh/x1min=-2.0', 'mesh/x1max=2.0',
        'mesh/ix1_bc=shear_periodic', 'mesh/ox1_bc=shear_periodic',
        'mesh/nx2=64', 'mesh/x2min=-2.0', 'mesh/x2max=2.0',
        'mesh/ix2_bc=periodic', 'mesh/ox2_bc=periodic',
        'mesh/nx3=1', 'mesh/x3min=-0.5', 'mesh/x3max=0.5',
        'hydro/iso_sound_speed=0.001', 'problem/ipert=3',
        'problem/amp=4.0e-4', 'problem/nwx=-4', 'problem/nwy=1',
        'orbital_advection/Omega0=1.0e-3', 'orbital_advection/qshear=1.5',
        'orbital_advection/shboxcoord=1', 'orbital_advection/OAorder=2',
        'problem/error_output=true', 'time/ncycle_out=0']
    athena.run('hydro/athinput.ssheet', arguments)


# Analyze outputs
def analyze():
    # set parameters
    ky = 0.5*math.pi
    kx0 = -2.0*math.pi
    Omega0 = 0.001
    qshear = 1.5
    kappa2 = 2.0*(2.0-qshear)*Omega0**2.0
    cs = 0.001
    constC = 0.5*(cs**2.0*ky**2.0+kappa2)/(qshear*Omega0*cs*ky)
    c1 = -1.82085106245
    c2 = -0.8208057072217868

    # read results of SSEET_SHWAVE
    fname = 'bin/SSHEET_SHWAVE.hst'
    a = athena_read.hst(fname)
    time1 = a['time']
    dvyc1 = a['dvyc']
    dvys1 = a['dvys']
    nf1 = len(time1)
    norm1c = 0.0
    norm1s = 0.0
    for n in range(nf1):
        tau_ = qshear*Omega0*time1[n]+kx0/ky
        T_ = 1.0j * cmath.sqrt(2.0*cs*ky/(qshear*Omega0))*tau_
        exp_ = cmath.exp(-0.25j*T_*T_)
        fterm_ = exp_*mp.hyp1f1(0.25-0.5j*constC, 0.5, 0.5j*T_*T_)
        first_ = fterm_.real
        exp_ = cmath.exp(0.25j*(math.pi-T_*T_))
        sterm_ = exp_*T_*mp.hyp1f1(0.75-0.5j*constC, 1.5, 0.5j*T_*T_)
        second_ = sterm_.real
        advy = c1*first_+c2*second_
        norm1c += abs(dvyc1[n]-advy)
        norm1s += abs(dvys1[n])
    norm1c /= nf1
    norm1s /= nf1

    # read results of SSEET_SHWAVE_ORB
    fname = 'bin/SSHEET_SHWAVE_ORB.hst'
    b = athena_read.hst(fname)
    time2 = b['time']
    dvyc2 = b['dvyc']
    dvys2 = b['dvys']
    nf2 = len(time2)
    norm2c = 0.0
    norm2s = 0.0
    for n in range(nf2):
        tau_ = qshear*Omega0*time2[n]+kx0/ky
        T_ = 1.0j * cmath.sqrt(2.0*cs*ky/(qshear*Omega0))*tau_
        exp_ = cmath.exp(-0.25j*T_*T_)
        fterm_ = exp_*mp.hyp1f1(0.25-0.5j*constC, 0.5, 0.5j*T_*T_)
        first_ = fterm_.real
        exp_ = cmath.exp(0.25j*(math.pi-T_*T_))
        sterm_ = exp_*T_*mp.hyp1f1(0.75-0.5j*constC, 1.5, 0.5j*T_*T_)
        second_ = sterm_.real
        advy = c1*first_+c2*second_
        norm2c += abs(dvyc2[n]-advy)
        norm2s += abs(dvys2[n])
    norm2c /= nf2
    norm2s /= nf2

    logger.warning('[SSHEET_SHWAVE]: L1-Norm Errors')
    msg = '[SSHEET_SHWAVE]: epsilon_c = {}, epsilon_s = {}'
    flag = True
    logger.warning('[SSHEET_SHWAVE]: w/o Orbital Advection')
    logger.warning(msg.format(norm1c, norm1s))
    if norm1c > 0.2:
        logger.warning('[SSHEET_SHWAVE]: dvyc Error is more than 20%.')
        flag = False
    if norm1s > 0.01:
        logger.warning('[SSHEET_SHWAVE]: dvys Error is more than 1%.')
        flag = False

    logger.warning('[SSHEET_SHWAVE]: w/  Orbital Advection')
    logger.warning(msg.format(norm2c, norm2s))
    if norm1c > 0.2:
        logger.warning('[SSHEET_SHWAVE]: dvyc Error is more than 20%.')
        flag = False
    if norm1s > 0.01:
        logger.warning('[SSHEET_SHWAVE]: dvys Error is more than 1%.')
        flag = False

    return flag
