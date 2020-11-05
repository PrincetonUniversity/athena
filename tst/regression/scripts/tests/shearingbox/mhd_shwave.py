# Regression test of shearing box and orbital advection with 3D MHD shwave.

# Modules
import logging
import math
import scripts.utils.athena as athena
import sys
sys.path.insert(0, '../../vis/python')
import athena_read  # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure('b', prob='jgg', flux='hlld',
                     eos='isothermal', **kwargs)
    athena.make()


# Run Athena++ w/wo Orbital Advection
def run(**kwargs):
    # w/o Orbital Advection
    arguments = [
        'job/problem_id=JGG',
        'output1/file_type=hst', 'output1/dt=0.01',
        'output1/data_format=%1.16f',
        'output2/file_type=vtk', 'output2/variable=prim',
        'output2/dt=1.0', 'time/cfl_number=0.3',
        'time/tlim=3.0', 'time/nlim=1000',
        'time/xorder=2', 'time/integrator=vl2',
        'mesh/nx1=32', 'mesh/x1min=-0.25', 'mesh/x1max=0.25',
        'mesh/ix1_bc=shear_periodic', 'mesh/ox1_bc=shear_periodic',
        'mesh/nx2=16', 'mesh/x2min=-0.25', 'mesh/x2max=0.25',
        'mesh/ix2_bc=periodic', 'mesh/ox2_bc=periodic',
        'mesh/nx3=16', 'mesh/x3min=-0.25', 'mesh/x3max=0.25',
        'mesh/ix3_bc=periodic', 'mesh/ox3_bc=periodic',
        'hydro/iso_sound_speed=1.0', 'problem/d0=1.0',
        'problem/amp=1.0e-6', 'problem/beta=20.0', 'problem/ipert=2',
        'problem/nwx=-2', 'problem/nwy=1', 'problem/nwz=1',
        'problem/Omega0=1.0', 'problem/qshear=1.5',
        'problem/orbital_advection=true',
        'time/ncycle_out=0']
    athena.run('mhd/athinput.jgg', arguments)

    # w/  Orbital Advection
    arguments = [
        'job/problem_id=JGG_ORB',
        'output1/file_type=hst', 'output1/dt=0.01',
        'output1/data_format=%1.16f',
        'output2/file_type=vtk', 'output2/variable=prim',
        'output2/dt=1.0', 'time/cfl_number=0.3',
        'time/tlim=3.0', 'time/nlim=1000',
        'time/xorder=2', 'time/integrator=vl2',
        'mesh/nx1=32', 'mesh/x1min=-0.25', 'mesh/x1max=0.25',
        'mesh/ix1_bc=shear_periodic', 'mesh/ox1_bc=shear_periodic',
        'mesh/nx2=16', 'mesh/x2min=-0.25', 'mesh/x2max=0.25',
        'mesh/ix2_bc=periodic', 'mesh/ox2_bc=periodic',
        'mesh/nx3=16', 'mesh/x3min=-0.25', 'mesh/x3max=0.25',
        'mesh/ix3_bc=periodic', 'mesh/ox3_bc=periodic',
        'hydro/iso_sound_speed=1.0', 'problem/d0=1.0',
        'problem/amp=1.0e-6', 'problem/beta=20.0', 'problem/ipert=2',
        'problem/nwx=-2', 'problem/nwy=1', 'problem/nwz=1',
        'problem/Omega0=1.0', 'problem/qshear=1.5',
        'problem/orbital_advection=true', 'problem/orbital_splitting=2',
        'time/ncycle_out=0']
    athena.run('mhd/athinput.jgg', arguments)


# Analyze outputs
def analyze():
    # set parameters
    nx = -2
    ny = 1
    nz = 1
    Lx = 0.5
    Ly = 0.5
    Lz = 0.5
    Omega0 = 1.0
    qshear = 1.5
    iso_cs = 1.0
    rho0 = 1.0
    epsilon = 1.0e-6
    beta = 20.0
    kx0 = nx*2.0*math.pi/Lx
    ky = ny*2.0*math.pi/Ly
    kz = nz*2.0*math.pi/Lz
    vA2 = iso_cs*iso_cs/beta
    vAx = math.sqrt(vA2/(kx0*kx0+ky*ky))*ky
    sch = Omega0/iso_cs

    # read results w/o Orbital Advection
    fname = 'bin/JGG.hst'
    a = athena_read.hst(fname)
    time1 = a['time']
    dby1 = a['dBy']
    nf1 = len(time1)
    norm1 = 0.0
    CS = 0.0
    for i in xrange(nf1):
        if i != 0:
            dt = time1[i]-time1[i-1]
            for j in xrange(10):
                time_temp = time1[i-1]+0.1*dt*(j+1)
                kx = kx0+qshear*Omega0*time_temp*2.0*math.pi/Lx
                k = math.sqrt(kx*kx+ky*ky+kz*kz)
                vAy = -math.sqrt(vA2/(kx0*kx0+ky*ky))*kx
                By = math.sqrt(rho0)*vAy
                omega_over_k2 = (iso_cs*iso_cs+vAx*vAx+vAy*vAy)
                CS += math.sqrt(omega_over_k2)*k*0.1*dt
        kx = kx0+qshear*Omega0*time1[i]*2.0*math.pi/Lx
        k = math.sqrt(kx*kx+ky*ky+kz*kz)
        vAy = -math.sqrt(vA2/(kx0*kx0+ky*ky))*kx
        By = math.sqrt(rho0)*vAy
        omega_over_k2 = (iso_cs*iso_cs+vAx*vAx+vAy*vAy)
        dBa = By*epsilon*math.sqrt(sch*k*math.sqrt(beta/(1.0+beta)))*math.cos(CS)
        if i == 0:
            dBy0 = dBa
        norm1 += abs(dby1[i]-dBa)/dBy0
    norm1 /= nf1

    # read results w/  Orbital Advection
    fname = 'bin/JGG_ORB.hst'
    b = athena_read.hst(fname)
    time2 = b['time']
    dby2 = b['dBy']
    nf2 = len(time2)
    norm2 = 0.0
    CS = 0.0
    for i in xrange(nf2):
        if i != 0:
            dt = time2[i]-time2[i-1]
            for j in xrange(10):
                time_temp = time2[i-1]+0.1*dt*(j+1)
                kx = kx0+qshear*Omega0*time_temp*2.0*math.pi/Lx
                k = math.sqrt(kx*kx+ky*ky+kz*kz)
                vAy = -math.sqrt(vA2/(kx0*kx0+ky*ky))*kx
                By = math.sqrt(rho0)*vAy
                omega_over_k2 = (iso_cs*iso_cs+vAx*vAx+vAy*vAy)
                CS += math.sqrt(omega_over_k2)*k*0.1*dt
        kx = kx0+qshear*Omega0*time2[i]*2.0*math.pi/Lx
        k = math.sqrt(kx*kx+ky*ky+kz*kz)
        vAy = -math.sqrt(vA2/(kx0*kx0+ky*ky))*kx
        By = math.sqrt(rho0)*vAy
        omega_over_k2 = (iso_cs*iso_cs+vAx*vAx+vAy*vAy)
        dBa = By*epsilon*math.sqrt(sch*k*math.sqrt(beta/(1.0+beta)))*math.cos(CS)
        if i == 0:
            dBy0 = dBa
        norm2 += abs(dby2[i]-dBa)/dBy0
    norm2 /= nf2

    msg = '[MHD_SHWAVE]: L1 Norm of b2 deviation {} = {}'
    logger.warning(msg.format('w/o Orbital Advection', norm1))
    logger.warning(msg.format('w/  Orbital Advection', norm2))
    flag = True
    if norm1 > 0.2:
        logger.warning('[MHD_SHWAVE]: L1 Norm is off by 20% w/o Orbital Advection')
        flag = False
    if norm2 > 0.2:
        logger.warning('[MHD_SHWAVE]: L1 Norm is off by 20% w/  Orbital Advection')
        flag = False

    return flag
