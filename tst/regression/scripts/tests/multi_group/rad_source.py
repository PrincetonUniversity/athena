# Regression test based on Newtonian hydro linear wave convergence problem
#
# Runs a linear wave convergence test in 3D including SMR and checks L1 errors (which
# are computed by the executable automatically and stored in the temporary file
# linearwave_errors.dat)

# Modules
import numpy as np
import sys
import scripts.utils.athena as athena
import os
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
        'radiation/Prat=1.0',
        'radiation/Crat=10.0',
        'radiation/error_limit=1.e-12',
        'radiation/n_frequency=3',
        'radiation/frequency_min=-4',
        'radiation/frequency_max=-8',
        'radiation/Compton=0',
        'radiation/T_unit=1.e4',
        'radiation/unit=0',
        'radiation/density_unit=1.e-10',
        'radiation/length_unit=3.e14',
        'radiation/cfl_rad=1',
        'radiation/Split_compton=0',
        'time/ncycle_out=100']
    athena.run('radiation/athinput.thermal_multigroup', arguments)
    bashcommand = "mv bin/Averaged_quantity.dat bin/Averaged_quantity1.dat"
    os.system(bashcommand)
    # case 2, with Compton
    arguments = [
        'problem/er_1=1.0',
        'problem/er_2=1.0',
        'problem/er_3=1.0',
        'problem/tgas=100.0',
        'problem/sigma_1=0.0',
        'problem/sigma_2=0.0',
        'problem/sigma_3=0.0',
        'radiation/Prat=0.546',
        'radiation/Crat=2.5485e4',
        'radiation/error_limit=1.e-12',
        'radiation/n_frequency=150',
        'radiation/frequency_min=-0.01',
        'radiation/frequency_max=-500',
        'radiation/Compton=1',
        'radiation/T_unit=1.e4',
        'radiation/unit=1',
        'radiation/density_unit=1.e-10',
        'radiation/length_unit=3.e14',
        'radiation/cfl_rad=0.01',
        'radiation/Split_compton=1',
        'problem/black_body=1',
        'problem/kappa_es=100',
        'time/tlim=1.e-2',
        'radiation/nmu=1',
        'time/ncycle_out=100']
    athena.run('radiation/athinput.thermal_multigroup', arguments)
    bashcommand = "mv bin/Averaged_quantity.dat bin/Averaged_quantity2.dat"
    os.system(bashcommand)

# Analyze outputs


def analyze():
    filename = 'bin/Averaged_quantity1.dat'
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

    filename = 'bin/Averaged_quantity2.dat'
    data = []
    with open(filename, 'r') as f:
        raw_data = f.readlines()
        for line in raw_data:
            if line.split()[0][0] == '#':
                continue
            data.append([float(val) for val in line.split()])

    ref_sol = np.array([7.170660e+01,  1.241716e-11,  5.827892e-12,
                        7.341434e-12,  9.225259e-12,  1.156634e-11,
                        1.447147e-11,  1.807164e-11,  2.252744e-11,
                        2.803555e-11,  3.483668e-11,  4.322521e-11,
                        5.356089e-11,  6.628308e-11,  8.192794e-11,
                        1.011494e-10,  1.247444e-10,  1.536836e-10,
                        1.891487e-10,  2.325771e-10,  2.857165e-10,
                        3.506902e-10,  4.300771e-10,  5.270068e-10,
                        6.452757e-10,  7.894864e-10,  9.652159e-10,
                        1.179218e-09,  1.439669e-09,  1.756458e-09,
                        2.141547e-09,  2.609389e-09,  3.177449e-09,
                        3.866812e-09,  4.702927e-09,  5.716495e-09,
                        6.944537e-09,  8.431675e-09,  1.023167e-08,
                        1.240928e-08,  1.504244e-08,  1.822497e-08,
                        2.206972e-08,  2.671239e-08,  3.231613e-08,
                        3.907702e-08,  4.723065e-08,  5.705999e-08,
                        6.890486e-08,  8.317327e-08,  1.003550e-07,
                        1.210378e-07,  1.459274e-07,  1.758701e-07,
                        2.118820e-07,  2.551820e-07,  3.072332e-07,
                        3.697914e-07,  4.449643e-07,  5.352827e-07,
                        6.437861e-07,  7.741264e-07,  9.306935e-07,
                        1.118767e-06,  1.344700e-06,  1.616146e-06,
                        1.942325e-06,  2.334365e-06,  2.805704e-06,
                        3.372587e-06,  4.054686e-06,  4.875845e-06,
                        5.865022e-06,  7.057435e-06,  8.495998e-06,
                        1.023310e-05,  1.233284e-05,  1.487380e-05,
                        1.795255e-05,  2.168804e-05,  2.622716e-05,
                        3.175177e-05,  3.848761e-05,  4.671561e-05,
                        5.678632e-05,  6.913835e-05,  8.432186e-05,
                        1.030288e-04,  1.261314e-04,  1.547320e-04,
                        1.902267e-04,  2.343873e-04,  2.894671e-04,
                        3.583360e-04,  4.446556e-04,  5.531036e-04,
                        6.896617e-04,  8.619869e-04,  1.079886e-03,
                        1.355922e-03,  1.706194e-03,  2.151325e-03,
                        2.717725e-03,  3.439181e-03,  4.358880e-03,
                        5.531945e-03,  7.028607e-03,  8.938171e-03,
                        1.137391e-02,  1.447910e-02,  1.843436e-02,
                        2.346657e-02,  2.985949e-02,  3.796633e-02,
                        4.822444e-02,  6.117204e-02,  7.746696e-02,
                        9.790711e-02,  1.234518e-01,  1.552431e-01,
                        1.946250e-01,  2.431586e-01,  3.026295e-01,
                        3.750426e-01,  4.626004e-01,  5.676557e-01,
                        6.926345e-01,  8.399184e-01,  1.011682e+00,
                        1.209678e+00,  1.434972e+00,  1.687621e+00,
                        1.966327e+00,  2.268066e+00,  2.587746e+00,
                        2.917934e+00,  3.248718e+00,  3.567774e+00,
                        3.860696e+00,  4.111657e+00,  4.304421e+00,
                        4.423673e+00,  4.456591e+00,  4.394496e+00,
                        4.234355e+00,  3.979875e+00,  3.641916e+00,
                        3.238031e+00,  2.791056e+00,  2.326840e+00,
                        6.649287e+00])
    max_err = np.abs(data[0][0] - ref_sol[0])/ref_sol[0]
    for i in range(1, 150):
        err = np.abs(data[0][i]-ref_sol[i])/ref_sol[i]
        if err > max_err:
            max_err = err
    # check absolute error and convergence of all three waves
    if max_err > 1.e-5:
        print("relative error too large:", data[0][0], data[0][1], max_err)

    return True
