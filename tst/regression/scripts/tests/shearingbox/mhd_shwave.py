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
        'orbital_advection/Omega0=1.0', 'orbital_advection/qshear=1.5',
        'orbital_advection/order=0', 'problem/error_output=true',
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
        'orbital_advection/Omega0=1.0', 'orbital_advection/qshear=1.5',
        'orbital_advection/order=2', 'problem/error_output=true',
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

    Bx = 0.1*iso_cs*math.sqrt(rho0)
    By0 = 0.2*iso_cs*math.sqrt(rho0)

    kx0 = float(nx)*2.0*math.pi/Lx
    ky = float(ny)*2.0*math.pi/Ly
    kz = float(nz)*2.0*math.pi/Lz

    # read results w/o Orbital Advection
    fname = 'bin/JGG.hst'
    a = athena_read.hst(fname)
    time1 = a['time']
    dByc1 = a['dByc']
    nf1 = len(time1)

    # initialize parameters
    n_ex = 8
    norm1 = 0.0
    for i in xrange(nf1):
        dst_time = time1[i]
        if (i == 0):
            dt = 0.0
            k = math.sqrt(kx0**2+ky**2+kz**2)
            dBya = epsilon*By0
            dBya *= math.sqrt(iso_cs*k/Omega0*math.sqrt(beta/(1.0+beta)))
            omega_tot = 0.0
            dBy0 = dBya
            dBya = 1.0
        else:
            # calculate omega_tot
            dt = dst_time/float(n_ex*i)
            omega = 0.0
            omega_tot = 0.0
            for j in xrange(n_ex*i+1):
                present_time = float(j)*dt
                kx = kx0+qshear*Omega0*present_time*ky
                k = math.sqrt(kx**2+ky**2+kz**2)
                By = By0-qshear*Omega0*present_time*Bx
                omega = math.sqrt(iso_cs**2+Bx**2+By**2)*k
                if (j == 0):
                    omega_tot = omega
                elif (j == n_ex*i):
                    omega_tot += omega
                elif ((j % 2) == 0):
                    omega_tot += 2.0*omega
                else:
                    omega_tot += 4.0*omega
            omega_tot *= dt/3.0
            # calculate dBy_analytic
            kx = kx0+qshear*Omega0*dst_time*ky
            k = math.sqrt(kx**2+ky**2+kz**2)
            By = By0-qshear*Omega0*dst_time*Bx
            dBya = epsilon*By*math.cos(omega_tot)/dBy0
            dBya *= math.sqrt(iso_cs*k/Omega0*math.sqrt(beta/(1.0+beta)))
        norm1 += abs(dByc1[i]-dBya)
    norm1 /= nf1

    # read results w/  Orbital Advection
    fname = 'bin/JGG_ORB.hst'
    b = athena_read.hst(fname)
    time2 = b['time']
    dByc2 = b['dByc']
    nf2 = len(time2)

    # initialize parameters
    n_ex = 8
    norm2 = 0.0
    for i in xrange(nf1):
        dst_time = time2[i]
        if (i == 0):
            dt = 0.0
            k = math.sqrt(kx0**2+ky**2+kz**2)
            dBya = epsilon*By0
            dBya *= math.sqrt(iso_cs*k/Omega0*math.sqrt(beta/(1.0+beta)))
            omega_tot = 0.0
            dBy0 = dBya
            dBya = 1.0
        else:
            # calculate omega_tot
            dt = dst_time/float(n_ex*i)
            omega = 0.0
            omega_tot = 0.0
            for j in xrange(n_ex*i+1):
                present_time = float(j)*dt
                kx = kx0+qshear*Omega0*present_time*ky
                k = math.sqrt(kx**2+ky**2+kz**2)
                By = By0-qshear*Omega0*present_time*Bx
                omega = math.sqrt(iso_cs**2+Bx**2+By**2)*k
                if (j == 0):
                    omega_tot = omega
                elif (j == n_ex*i):
                    omega_tot += omega
                elif ((j % 2) == 0):
                    omega_tot += 2.0*omega
                else:
                    omega_tot += 4.0*omega
            omega_tot *= dt/3.0
            # calculate dBy_analytic
            kx = kx0+qshear*Omega0*dst_time*ky
            k = math.sqrt(kx**2+ky**2+kz**2)
            By = By0-qshear*Omega0*dst_time*Bx
            dBya = epsilon*By*math.cos(omega_tot)/dBy0
            dBya *= math.sqrt(iso_cs*k/Omega0*math.sqrt(beta/(1.0+beta)))
        norm2 += abs(dByc2[i]-dBya)
    norm2 /= nf2

    logger.warning('[MHD_SHWAVE]: L1-Norm Errors')
    msg = '[MHD_SHWAVE]: epsilon_c = {}'
    flag = True
    logger.warning('[MHD_SHWAVE]: w/o Orbital Advection')
    logger.warning(msg.format(norm1))
    if norm1 > 0.20:
        logger.warning('[MHD_SHWAVE]: dByc Error is more than 20%.')
        flag = False
    logger.warning('[MHD_SHWAVE]: w/  Orbital Advection')
    logger.warning(msg.format(norm2))
    if norm2 > 0.20:
        logger.warning('[MHD_SHWAVE]: dByc Error is more than 20%.')
        flag = False

    return flag
