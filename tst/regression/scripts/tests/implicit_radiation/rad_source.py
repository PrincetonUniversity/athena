# Regression test based on thermalize relaxation problem
# for full radiation hydro equations, using the
# implicit radiation hydro module

# Modules
import numpy as np
import sys
import scripts.utils.athena as athena
import os
sys.path.insert(0, '../../vis/python')

# Prepare Athena++


def prepare(**kwargs):
    athena.configure('implicit_radiation',
                     prob='thermal_relaxation',
                     coord='cartesian',
                     flux='hllc')
    athena.make()

# Run Athena++


def run(**kwargs):
    # case 1
    arguments = ['problem/er=10.0', 'problem/tgas=1.0', 'problem/sigma=1.0',
                 'radiation/prat=0.01', 'radiation/crat=10.0',
                 'radiation/error_limit=1.e-12',
                 'time/ncycle_out=100']
    athena.run('radiation/athinput.thermal_relaxation', arguments)
    bashcommand = "mv bin/thermal.hst bin/thermal1.hst"
    os.system(bashcommand)
    # case 2
    arguments = ['problem/er=10.0', 'problem/tgas=1.0', 'problem/sigma=100.0',
                 'radiation/prat=100.0', 'radiation/crat=10.0',
                 'radiation/error_limit=1.e-12',
                 'time/ncycle_out=100']
    athena.run('radiation/athinput.thermal_relaxation', arguments)
    bashcommand = "mv bin/thermal.hst bin/thermal2.hst"
    os.system(bashcommand)
    # case 3
    arguments = ['problem/er=1.0', 'problem/tgas=10.0', 'problem/sigma=100.0',
                 'radiation/prat=1.0', 'radiation/crat=10.0',
                 'radiation/error_limit=1.e-12',
                 'time/ncycle_out=100']
    athena.run('radiation/athinput.thermal_relaxation', arguments)
    bashcommand = "mv bin/thermal.hst bin/thermal3.hst"
    os.system(bashcommand)
# Analyze outputs


def analyze():
    filename = 'bin/thermal1.hst'
    data = []
    with open(filename, 'r') as f:
        raw_data = f.readlines()
        for line in raw_data:
            if line.split()[0][0] == '#':
                continue
            data.append([float(val) for val in line.split()])
    dim = np.shape(data)

    # check absolute error and convergence of all three waves
    if np.abs(data[dim[0]-1][9]-1.58746) > 1.0e-5 or \
            + np.abs(data[dim[0]-1][10]-1.25442) > 1.0e-5:
        print("error in case 1: tgas or Er", data[dim[0]-1][9], data[dim[0]-1][10])
        return False

    filename = 'bin/thermal2.hst'
    data = []
    with open(filename, 'r') as f:
        raw_data = f.readlines()
        for line in raw_data:
            if line.split()[0][0] == '#':
                continue
            data.append([float(val) for val in line.split()])
    dim = np.shape(data)

    # check absolute error and convergence of all three waves
    if np.abs(data[dim[0]-1][9]-2.66664) > 1.0e-3 or \
            + np.abs(data[dim[0]-1][10]-9.98833) > 1.0e-3:
        print("error in case 2: tgas or Er", data[dim[0]-1][9], data[dim[0]-1][10])
        return False

    filename = 'bin/thermal3.hst'
    data = []
    with open(filename, 'r') as f:
        raw_data = f.readlines()
        for line in raw_data:
            if line.split()[0][0] == '#':
                continue
            data.append([float(val) for val in line.split()])
    dim = np.shape(data)

    # check absolute error and convergence of all three waves
    if np.abs(data[dim[0]-1][9]-2.85609) > 1.0e-5 or \
            + np.abs(data[dim[0]-1][10]-13.1439) > 1.0e-4:
        print("error in case 3: tgas or Er", data[dim[0]-1][9], data[dim[0]-1][10])
        return False

    return True
