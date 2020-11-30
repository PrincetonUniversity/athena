# Regression test of shearing box with 3d MRI

# Modules
import logging
import numpy as np
import scripts.utils.athena as athena
import sys
sys.path.insert(0, '../../vis/python')
import athena_read  # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure('b', prob='hgb', flux='hlld',
                     eos='isothermal', **kwargs)
    athena.make()


# Run Athena++ w/wo Orbital Advection
def run(**kwargs):
    arguments = [
        'job/problem_id=HGB',
        'output1/file_type=hst', 'output1/dt=0.062831853',
        'output2/file_type=vtk', 'output2/variable=prim', 'output2/dt=31.4616',
        'time/cfl_number=0.3', 'time/tlim=62.83185', 'time/nlim=10000',
        'time/xorder=2', 'time/integrator=vl2', 'time/ncycle_out=10',
        'mesh/nx1=32', 'mesh/x1min=-0.5', 'mesh/x1max=0.5',
        'mesh/ix1_bc=shear_periodic', 'mesh/ox1_bc=shear_periodic',
        'mesh/nx2=24', 'mesh/x2min=-1.57079632679', 'mesh/x2max=1.57079632679',
        'mesh/ix2_bc=periodic', 'mesh/ox2_bc=periodic',
        'mesh/nx3=32', 'mesh/x3min=-0.5', 'mesh/x3max=0.5',
        'mesh/ix3_bc=periodic', 'mesh/ox3_bc=periodic',
        'meshblock/nx1=32', 'meshblock/nx2=24', 'meshblock/nx3=32',
        'hydro/iso_sound_speed=1.0', 'problem/beta=100',
        'problem/d0=1.0', 'problem/amp=0.025', 'problem/nwx=-2', 'problem/nwy=1',
        'problem/ipert=1', 'problem/ifield=1',
        'orbital_advection/Omega0=1.0', 'orbital_advection/qshear=1.5',
        'orbital_advection/OAorder=0', 'time/ncycle_out=0']
    athena.run('mhd/athinput.hgb', arguments)


# Analyze outputs
def analyze():
    # omg  = 1.0e-3 # unused
    rho0 = 1.0
    cs = 1.0
    pres = rho0 * cs**2
    vol = 1.0 * np.pi * 1.0
    index = -500

    # reference
    fname = 'data/mhd_mri_3d.hst'
    a = athena_read.hst(fname)
    me = (a['1-ME'] + a['2-ME'] + a['3-ME'])
    ref_stress = np.average(a['-BxBy'][index:] / vol / pres)
    ref_me = np.average(me[index:] / vol / pres)
    ref_ratio = ref_me / ref_stress

    # resutls
    fname = 'bin/HGB.hst'
    b = athena_read.hst(fname)
    me = (b['1-ME'] + b['2-ME'] + b['3-ME'])
    new_stress = np.average(b['-BxBy'][index:] / vol / pres)
    new_me = np.average(me[index:] / vol / pres)
    new_ratio = new_me / new_stress

    msg = '[MRI-3D]: {}(stress,ME,ratio) = {} {} {}'
    logger.warning(msg.format('Ref', ref_stress, ref_me, ref_ratio))
    logger.warning(msg.format('New', new_stress, new_me, new_ratio))
    flag = True
    error_rel = np.fabs((new_stress / ref_stress) - 1.0)
    if error_rel > 0.5:
        logger.warning('[MRI-3D]: averaged stress is off by a factor > 2')
        flag = False
    error_rel = np.fabs((new_me / ref_me) - 1.0)
    if error_rel > 0.5:
        logger.warning('[MRI-3D]: averaged magnetic energy is off by a factor > 2')
        flag = False
    error_rel = np.fabs(new_ratio - ref_ratio)
    if error_rel > 1.0:
        logger.warning('[MRI-3D]: energy-to-stress ratio is off by an amount > 1.0')
        flag = False

    return flag
