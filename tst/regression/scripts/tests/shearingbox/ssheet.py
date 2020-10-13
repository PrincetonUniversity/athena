# Regression test of shearing box and orbital advection with 2d hydro shwave.

# Modules
import logging
import math, cmath
import numpy as np
import os
import scripts.utils.athena as athena
import sys
sys.path.insert(0, '../../vis/python')
import athena_read  # noqa
from mpmath import *
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
    # w/o Orbital Advection for ipert=1 w/o passive scalars
    arguments = [
        'job/problem_id=ssheet_ipert1_na',
        'output1/file_type=hst', 'output1/dt=10.0',
        'output2/file_type=vtk', 'output2/variable=prim',
        'output2/dt=-1', 'time/cfl_number=0.4',
        'time/tlim=8000.0', 'time/nlim=2000',
        'time/xorder=2', 'time/integrator=vl2', 'time/ncycle_out=10',
        'mesh/nx1=64', 'mesh/x1min=-2.0', 'mesh/x1max=2.0',
        'mesh/ix1_bc=shear_periodic', 'mesh/ox1_bc=shear_periodic',
        'mesh/nx2=64', 'mesh/x2min=-2.0', 'mesh/x2max=2.0',
        'mesh/ix2_bc=periodic', 'mesh/ox2_bc=periodic',
        'mesh/nx3=1', 'mesh/x3min=-0.5', 'mesh/x3max=0.5',
        'hydro/iso_sound_speed=0.001', 'problem/ipert=1',
        'problem/amp=4.0e-4', 'problem/nwx=-4', 'problem/nwy=1', 'problem/nwz=0',
        'problem/Omega0=1.0e-3', 'problem/qshear=1.5', 'problem/shboxcoord=1',
        'problem/orbital_advection=false', 'time/ncycle_out=0']
    athena.run('hydro/athinput.ssheet', arguments)

    # w/  Orbital Advection for ipert=1 w/o passive scalars
    arguments = [
        'job/problem_id=ssheet_ipert1_oa',
        'output1/file_type=hst', 'output1/dt=10.0',
        'output2/file_type=vtk', 'output2/variable=prim',
        'output2/dt=-1', 'time/cfl_number=0.4',
        'time/tlim=8000.0', 'time/nlim=500',
        'time/xorder=2', 'time/integrator=vl2', 'time/ncycle_out=10',
        'mesh/nx1=64', 'mesh/x1min=-2.0', 'mesh/x1max=2.0',
        'mesh/ix1_bc=shear_periodic', 'mesh/ox1_bc=shear_periodic',
        'mesh/nx2=64', 'mesh/x2min=-2.0', 'mesh/x2max=2.0',
        'mesh/ix2_bc=periodic', 'mesh/ox2_bc=periodic',
        'mesh/nx3=1', 'mesh/x3min=-0.5', 'mesh/x3max=0.5',
        'hydro/iso_sound_speed=0.001', 'problem/ipert=1',
        'problem/amp=4.0e-4', 'problem/nwx=-4', 'problem/nwy=1', 'problem/nwz=0',
        'problem/Omega0=1.0e-3', 'problem/qshear=1.5', 'problem/shboxcoord=1',
        'problem/orbital_advection=true', 'time/ncycle_out=0']
    athena.run('hydro/athinput.ssheet', arguments)

    os.system('rm -rf obj')
    os.system('mv obj_scalar obj')
    os.system('mv bin/athena_scalar bin/athena')
    # w/o Orbital Advection for ipert=0 w/  passive scalars
    arguments = [
        'job/problem_id=ssheet_ipert0_na',
        'output1/file_type=hst', 'output1/dt=10.0',
        'output2/file_type=vtk', 'output2/variable=prim',
        'output2/dt=-1', 'time/cfl_number=0.4',
        'time/tlim=8000.0', 'time/nlim=2000',
        'time/xorder=2', 'time/integrator=vl2', 'time/ncycle_out=10',
        'mesh/nx1=64', 'mesh/x1min=-2.0', 'mesh/x1max=2.0',
        'mesh/ix1_bc=shear_periodic', 'mesh/ox1_bc=shear_periodic',
        'mesh/nx2=64', 'mesh/x2min=-2.0', 'mesh/x2max=2.0',
        'mesh/ix2_bc=periodic', 'mesh/ox2_bc=periodic',
        'mesh/nx3=1', 'mesh/x3min=-0.5', 'mesh/x3max=0.5',
        'hydro/iso_sound_speed=0.001', 'problem/ipert=0',
        'problem/amp=4.0e-4', 'problem/nwx=-4', 'problem/nwy=1', 'problem/nwz=0',
        'problem/Omega0=1.0e-3', 'problem/qshear=1.5', 'problem/shboxcoord=1',
        'problem/orbital_advection=false', 'time/ncycle_out=0']
    athena.run('hydro/athinput.ssheet', arguments)

    # w/  Orbital Advection for ipert=0 w/  passive scalars
    arguments = [
        'job/problem_id=ssheet_ipert0_oa',
        'output1/file_type=hst', 'output1/dt=10.0',
        'output2/file_type=vtk', 'output2/variable=prim',
        'output2/dt=-1', 'time/cfl_number=0.4',
        'time/tlim=8000.0', 'time/nlim=500',
        'time/xorder=2', 'time/integrator=vl2', 'time/ncycle_out=10',
        'mesh/nx1=64', 'mesh/x1min=-2.0', 'mesh/x1max=2.0',
        'mesh/ix1_bc=shear_periodic', 'mesh/ox1_bc=shear_periodic',
        'mesh/nx2=64', 'mesh/x2min=-2.0', 'mesh/x2max=2.0',
        'mesh/ix2_bc=periodic', 'mesh/ox2_bc=periodic',
        'mesh/nx3=1', 'mesh/x3min=-0.5', 'mesh/x3max=0.5',
        'hydro/iso_sound_speed=0.001', 'problem/ipert=0',
        'problem/amp=4.0e-4', 'problem/nwx=-4', 'problem/nwy=1', 'problem/nwz=0',
        'problem/Omega0=1.0e-3', 'problem/qshear=1.5', 'problem/shboxcoord=1',
        'problem/orbital_advection=true', 'time/ncycle_out=0']
    athena.run('hydro/athinput.ssheet', arguments)


# Analyze outputs
def analyze():
    # ipert=1
    # set parameters 
    ky     = 0.5*pi
    kx0    = -2.0*pi
    Omega0 = 0.001
    qshear = 1.5
    kappa2 = 2.0*(2.0-qshear)*Omega0**2.0
    cs     = 0.001
    constC = 0.5*(cs**2.0*ky**2.0+kappa2)/(qshear*Omega0*cs*ky)

    c1 = -1.82088e-07
    c2 = -8.20766e-08
    dvy0 = 1.0e-7

    # read results w/o Orbital Advection
    fname = 'bin/ssheet_ipert1_na.hst'
    a    = athena_read.hst(fname)
    time1 = a['time']
    dvyc1 = a['dvyc']
    nf1   = len(time1)
    norm1 = 0.0
    for n in xrange(nf1):
        tau_    = qshear*Omega0*time1[n]+kx0/ky
        T_      = 1.0j * cmath.sqrt(2.0*cs*ky/(qshear*Omega0))*tau_
        exp_    = cmath.exp(-0.25j*T_*T_)
        fterm_  = exp_*hyp1f1(0.25-0.5j*constC, 0.5, 0.5j*T_*T_)
        first_  = fterm_.real
        exp_    = cmath.exp(0.25j*(pi-T_*T_))
        sterm_  = exp_*T_*hyp1f1(0.75-0.5j*constC, 1.5, 0.5j*T_*T_)
        second_ = sterm_.real
        advy = c1*first_+c2*second_
        norm1 += abs(dvyc1[n]-advy)/dvy0
    norm1 /= nf1

    # read results w/  Orbital Advection
    fname = 'bin/ssheet_ipert1_oa.hst'
    b    = athena_read.hst(fname)
    time2 = b['time']
    dvyc2 = b['dvyc']
    nf2   = len(time2)
    norm2 = 0.0
    for n in xrange(nf2):
        tau_    = qshear*Omega0*time2[n]+kx0/ky
        T_      = 1.0j * cmath.sqrt(2.0*cs*ky/(qshear*Omega0))*tau_
        exp_    = cmath.exp(-0.25j*T_*T_)
        fterm_  = exp_*hyp1f1(0.25-0.5j*constC, 0.5, 0.5j*T_*T_)
        first_  = fterm_.real
        exp_    = cmath.exp(0.25j*(pi-T_*T_))
        sterm_  = exp_*T_*hyp1f1(0.75-0.5j*constC, 1.5, 0.5j*T_*T_)
        second_ = sterm_.real
        advy = c1*first_+c2*second_
        norm2 += abs(dvyc2[n]-advy)/dvy0
    norm2 /= nf2

    logger.warning('[ssheet]: check L1 norm of the vyc deviation for ipert=1')
    msg = '[ssheet]: L1 Norm {} = {}'
    logger.warning(msg.format('w/o Orbital Advection', norm1))
    logger.warning(msg.format('w/  Orbital Advection', norm2))
    flag = True
    if norm1 > 0.2:
        logger.warning('[ssheet]: the deviation is more than 20% w/o Orbital Advection')
        flag = False
    if norm2 > 0.2:
        logger.warning('[ssheet]: the deviation is more than 20% w/  Orbital Advection')
        flag = False

    # ipert=0 w/o orbital advection
    amp   = 4.0e-4
    fname = 'bin/ssheet_ipert0_na.hst'
    c     = athena_read.hst(fname)
    time3   = c['time']
    scalar3 = c['ghost_scalar']
    nf3     = len(time3)
    norm3   = scalar3[nf3-1]/amp

    # ipert=0 w/  orbital advection
    fname = 'bin/ssheet_ipert0_oa.hst'
    d     = athena_read.hst(fname)
    time4   = d['time']
    scalar4 = d['ghost_scalar']
    nf4     = len(time4)
    norm4   = scalar4[nf4-1]/amp

    logger.warning('[ssheet] check L1 norm of the scalar deviation for ipert=0')
    logger.warning(msg.format('w/o Orbital Advection', norm3))
    logger.warning(msg.format('w/  Orbital Advection', norm4))
    if norm3 > 0.1:
        logger.warning('[ssheet]: the deviation is more than 10% w/o Orbital Advection')
        flag = False
    if norm4 > 0.1:
        logger.warning('[ssheet]: the deviation is more than 10% w/  Orbital Advection')
        flag = False

    return flag
