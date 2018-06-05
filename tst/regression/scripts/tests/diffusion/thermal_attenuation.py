# Regression test based on the decaying linear wave due to thermal
# conduction. The decay rate is fit and then compared with analytic
# solution

# Modules
import numpy as np
import sys
import scripts.utils.athena as athena
sys.path.insert(0, '../../vis/python')
import athena_read                            # noqa


def prepare(**kwargs):
    athena.configure(
        prob='linear_wave',
        flux='hllc',
        eos='adiabatic')
    athena.make()


def run(**kwargs):
    arguments0 = ['output1/file_type=hst', 'output1/dt=0.01',
                  'output2/file_type=vtk', 'output2/variable=prim', 'output2/dt=0.03',
                  'time/cfl_number=0.3', 'time/tlim=3.0', 'time/nlim=-1',
                  'time/xorder=2', 'time/integrator=vl2', 'time/ncycle_out=0',
                  'mesh/nx1=64', 'mesh/x1min=0.0', 'mesh/x1max=3.0',
                  'mesh/ix1_bc=periodic', 'mesh/ox1_bc=periodic',
                  'mesh/nx2=32', 'mesh/x2min=0.0', 'mesh/x2max=1.5',
                  'mesh/ix2_bc=periodic', 'mesh/ox2_bc=periodic',
                  'mesh/nx3=32', 'mesh/x3min=0.0', 'mesh/x3max=1.5',
                  'mesh/ix3_bc=periodic', 'mesh/ox3_bc=periodic',
                  'mesh/num_threads=1', 'mesh/refinement=none',
                  'meshblock/nx1=64', 'meshblock/nx2=32', 'meshblock/nx3=32',
                  'hydro/gamma=1.666666666666667', 'hydro/iso_sound_speed=1.0',
                  # L-going fast wave
                  'problem/compute_error=false', 'problem/wave_flag=0',
                  'problem/amp=1.0e-4', 'problem/vflow=0.0',
                  'problem/nu_iso=0.00', 'problem/kappa_iso=0.04']
    arguments = arguments0+['job/problem_id=DecayLinWave']
    athena.run('hydro/athinput.linear_wave3d', arguments)


def analyze():
    ksqr = (2.0*np.pi)**2
    # decay rate = (19\nu/4+3\eta+3(\gamma-1)^2*kappa/gamma/4)*(2/15)*k^2
    # (4nu/3+(gamma-1)^2/gamma*kappa)*k^2/2 decay rate of sound wave w/ thermal conduction
    rate = 2.0*0.04/15.0*ksqr

    basename = 'bin/DecayLinWave.block0.out2.'
    nframe = 101
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
