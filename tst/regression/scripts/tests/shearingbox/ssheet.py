# Regression test of shearing box and orbital advection with 2d hydro shwave.

# Modules
import logging
import cmath
import math
import mpmath as mp
import os
import scripts.utils.athena as athena
import sys
sys.path.insert(0, '../../vis/python')
import athena_read  # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure(prob='ssheet', nscalars=1,
                     eos='isothermal', **kwargs)
    athena.make()
    os.system('mv bin/athena bin/athena_scalar')
    os.system('mv obj obj_scalar')
    athena.configure(prob='ssheet',
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
        'problem/Omega0=1.0e-3', 'problem/qshear=1.5', 'problem/shboxcoord=1',
        'problem/orbital_advection=false',
        'time/ncycle_out=0']
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
        'problem/Omega0=1.0e-3', 'problem/qshear=1.5', 'problem/shboxcoord=1',
        'problem/orbital_advection=true', 'problem/orbital_splitting=2',
        'time/ncycle_out=0']
    athena.run('hydro/athinput.ssheet', arguments)

    os.system('rm -rf obj')
    os.system('mv obj_scalar obj')
    os.system('mv bin/athena_scalar bin/athena')
    # passive scalar in shearingbox w/o Orbital Advection
    arguments = [
        'job/problem_id=SSHEET_SCALAR',
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
        'hydro/iso_sound_speed=0.001', 'problem/ipert=1',
        'problem/amp=4.0e-4', 'problem/nwx=-4', 'problem/nwy=1',
        'problem/Omega0=1.0e-3', 'problem/qshear=1.5', 'problem/shboxcoord=1',
        'problem/orbital_advection=false',
        'time/ncycle_out=0']
    athena.run('hydro/athinput.ssheet', arguments)

    # passive scalar in shearingbox w/  Orbital Advection
    arguments = [
        'job/problem_id=SSHEET_SCALAR_ORB',
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
        'hydro/iso_sound_speed=0.001', 'problem/ipert=1',
        'problem/amp=4.0e-4', 'problem/nwx=-4', 'problem/nwy=1',
        'problem/Omega0=1.0e-3', 'problem/qshear=1.5', 'problem/shboxcoord=1',
        'problem/orbital_advection=true', 'problem/orbital_splitting=2',
        'time/ncycle_out=0']
    athena.run('hydro/athinput.ssheet', arguments)


# Analyze outputs
def analyze():
    # HD SHWAVE
    # set parameters
    ky = 0.5*math.pi
    kx0 = -2.0*math.pi
    Omega0 = 0.001
    qshear = 1.5
    kappa2 = 2.0*(2.0-qshear)*Omega0**2.0
    cs = 0.001
    constC = 0.5*(cs**2.0*ky**2.0+kappa2)/(qshear*Omega0*cs*ky)

    c1 = -1.82088e-07
    c2 = -8.20766e-08
    dvy0 = 1.0e-7

    # read results of SSEET_SHWAVE
    fname = 'bin/SSHEET_SHWAVE.hst'
    a = athena_read.hst(fname)
    time1 = a['time']
    dvyc1 = a['dvyc']
    nf1 = len(time1)
    norm1 = 0.0
    for n in xrange(nf1):
        tau_ = qshear*Omega0*time1[n]+kx0/ky
        T_ = 1.0j * cmath.sqrt(2.0*cs*ky/(qshear*Omega0))*tau_
        exp_ = cmath.exp(-0.25j*T_*T_)
        fterm_ = exp_*mp.hyp1f1(0.25-0.5j*constC, 0.5, 0.5j*T_*T_)
        first_ = fterm_.real
        exp_ = cmath.exp(0.25j*(math.pi-T_*T_))
        sterm_ = exp_*T_*mp.hyp1f1(0.75-0.5j*constC, 1.5, 0.5j*T_*T_)
        second_ = sterm_.real
        advy = c1*first_+c2*second_
        norm1 += abs(dvyc1[n]-advy)/dvy0
    norm1 /= nf1

    # read results of SSEET_SHWAVE_ORB
    fname = 'bin/SSHEET_SHWAVE_ORB.hst'
    b = athena_read.hst(fname)
    time2 = b['time']
    dvyc2 = b['dvyc']
    nf2 = len(time2)
    norm2 = 0.0
    for n in xrange(nf2):
        tau_ = qshear*Omega0*time2[n]+kx0/ky
        T_ = 1.0j * cmath.sqrt(2.0*cs*ky/(qshear*Omega0))*tau_
        exp_ = cmath.exp(-0.25j*T_*T_)
        fterm_ = exp_*mp.hyp1f1(0.25-0.5j*constC, 0.5, 0.5j*T_*T_)
        first_ = fterm_.real
        exp_ = cmath.exp(0.25j*(math.pi-T_*T_))
        sterm_ = exp_*T_*mp.hyp1f1(0.75-0.5j*constC, 1.5, 0.5j*T_*T_)
        second_ = sterm_.real
        advy = c1*first_+c2*second_
        norm2 += abs(dvyc2[n]-advy)/dvy0
    norm2 /= nf2

    logger.warning('[SSHEET_SHWAVE]: check L1 norm of the dvyc deviation')
    msg = '[SSHEET_SHWAVE]: L1 Norm {} = {}'
    logger.warning(msg.format('w/o Orbital Advection', norm1))
    logger.warning(msg.format('w/  Orbital Advection', norm2))
    flag = True
    if norm1 > 0.2:
        logger.warning('[SSHEET_SHWAVE]: the deviation is more than 20%')
        logger.warning('                 w/o Orbital Advection')
        flag = False
    if norm2 > 0.2:
        logger.warning('[SSHEET_SHWAVE]: the deviation is more than 20%')
        logger.warning('                 w/  Orbital Advection')
        flag = False

    # passive in shearingbox
    amp = 4.0e-4
    # read results of SSHEET_SCALAR
    fname = 'bin/SSHEET_SCALAR.hst'
    c = athena_read.hst(fname)
    time3 = c['time']
    scalar3 = c['ghost_scalar']
    nf3 = len(time3)
    norm3 = scalar3[nf3-1]/amp

    # read results of SSHEET_SCALAR_ORB
    fname = 'bin/SSHEET_SCALAR_ORB.hst'
    d = athena_read.hst(fname)
    time4 = d['time']
    scalar4 = d['ghost_scalar']
    nf4 = len(time4)
    norm4 = scalar4[nf4-1]/amp

    logger.warning('[SSHEET_SCALAR] check L1 norm of the scalar deviation')
    msg = '[SSHEET_SCALAR]: L1 Norm {} = {}'
    logger.warning(msg.format('w/o Orbital Advection', norm3))
    logger.warning(msg.format('w/  Orbital Advection', norm4))
    if norm3 > 0.1:
        logger.warning('[SSHEET_SCALAR]: the deviation is more than 10%')
        logger.warning('                 w/o Orbital Advection')
        flag = False
    if norm4 > 0.1:
        logger.warning('[SSHEET_SCALAR]: the deviation is more than 10%')
        logger.warning('                 w/  Orbital Advection')
        flag = False

    return flag
