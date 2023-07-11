# Regression test based on thermalize relaxation problem
# for full radiation hydro equations, using the
# implicit multi-group radiation module

# Modules
import numpy as np
import sys
import scripts.utils.athena as athena
sys.path.insert(0, '../../vis/python')

# Prepare Athena++


def prepare(**kwargs):
    athena.configure('implicit_radiation',
                     prob='thermal_multigroup',
                     coord='cartesian',
                     flux='hllc')
    athena.make()

# Run Athena++


def run(**kwargs):
    # case 1
    arguments = [
        'problem/er_1=10.0',
        'problem/er_2=20.0',
        'problem/er_3=30.0',
        'problem/tgas=1.0',
        'problem/sigma_1=100.0',
        'problem/sigma_2=200.0',
        'problem/sigma_1=300.0',
        'radiation/prat=1.0',
        'radiation/crat=10.0',
        'radiation/error_limit=1.e-12',
        'radiation/n_frequency=3',
        'radiation/frequency_min=-4',
        'radiation/frequency_max=-8',
        'time/ncycle_out=100']
    athena.run('radiation/athinput.thermal_multigroup', arguments)

# Analyze outputs


def analyze():
    filename = 'bin/Averaged_quantity.dat'
    data = []
    with open(filename, 'r') as f:
        raw_data = f.readlines()
        for line in raw_data:
            if line.split()[0][0] == '#':
                continue
            data.append([float(val) for val in line.split()])

    # check absolute error and convergence of all three waves
    ref_sol = np.array([2.752165, 5.044411, 1.632063e+01, 3.600671e+01])
    max_err = np.abs(data[0][0]-ref_sol[0])
    for i in range(1, 4):
        err = np.abs(data[0][i]-ref_sol[i])
        if err > max_err:
            max_err = err

    if max_err > 1.e-5:
        print("error in case 1: T, Er0, Er1, Er2:", data[0], data[1], data[2], data[3])
        return False

    return True
