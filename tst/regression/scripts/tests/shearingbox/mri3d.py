# Regression test of shearing box with 3d MRI

# Modules
import numpy as np
import scripts.utils.athena as athena


def prepare(**kwargs):
    athena.configure('b', 'shear',
                     prob='hgb',
                     flux='hlld',
                     eos='isothermal', **kwargs)
    athena.make()


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
        'hydro/iso_sound_speed=1.0', 'problem/beta=100', 'problem/amp=0.025',
        'problem/ipert=1', 'problem/ifield=1',
        'problem/Omega0=1.0', 'problem/qshear=1.5', 'time/ncycle_out=0']

    athena.run('mhd/athinput.hgb', arguments)


def analyze():

    # omg  = 1.0e-3 # unused
    rho0 = 1.0
    cs = 1.0
    pres = rho0 * cs**2
    vol = 1.0 * np.pi * 1.0
    index = -500
    dtype_array = np.dtype([('me1', 'f8'), ('me2', 'f8'), ('me3', 'f8'),
                            ('stress', 'f8')])
    col_indices = (9, 10, 11, 12)

    fname = 'data/mhd_mri_3d.hst'
    a = np.loadtxt(fname, dtype=dtype_array, skiprows=2, usecols=col_indices)
    me = (a['me1'] + a['me2'] + a['me3'])
    ref_stress = np.average(a['stress'][index:] / vol / pres)
    ref_me = np.average(me[index:] / vol / pres)
    ref_ratio = ref_me / ref_stress

    fname = 'bin/HGB.hst'
    b = np.loadtxt(fname, dtype=dtype_array, skiprows=2, usecols=col_indices)
    me = (b['me1'] + b['me2'] + b['me3'])
    new_stress = np.average(b['stress'][index:] / vol / pres)
    new_me = np.average(me[index:] / vol / pres)
    new_ratio = new_me / new_stress

    print('[MRI-3D]: Ref(stress,ME,ratio) = {} {} {}'.format(ref_stress,
                                                             ref_me,
                                                             ref_ratio))
    print('[MRI-3D]: New(stress,ME,ratio) = {} {} {}'.format(new_stress,
                                                             new_me,
                                                             new_ratio))
    flag = True
    error_rel = np.fabs((new_stress / ref_stress) - 1.0)
    if error_rel > 0.5:
        print('[MRI-3D]: averaged stress is off by a factor > 2')
        flag = False
    error_rel = np.fabs((new_me / ref_me) - 1.0)
    if error_rel > 0.5:
        print('[MRI-3D]: averaged magnetic energy is off by a factor > 2')
        flag = False
    error_rel = np.fabs(new_ratio - ref_ratio)
    if error_rel > 1.0:
        print('[MRI-3D]: energy-to-stress ratio is off by an amount > 1.0')
        flag = False

    return flag
