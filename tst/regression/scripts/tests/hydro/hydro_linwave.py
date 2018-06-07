# Regression test based on Newtonian hydro linear wave convergence problem
#
# Runs a linear wave convergence test in 3D including SMR and checks L1 errors (which
# are computed by the executable automatically and stored in the temporary file
# linearwave_errors.dat)

# Modules
import scripts.utils.athena as athena


# Prepare Athena++
def prepare(**kwargs):
    athena.configure(
        prob='linear_wave',
        coord='cartesian',
        flux='hllc', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    # L-going sound wave
    for i in (32, 64):
        arguments = ['time/ncycle_out=0',
                     'problem/wave_flag=0',
                     'problem/vflow=0.0',
                     'mesh/refinement=static',
                     'mesh/nx1=' + repr(i),
                     'mesh/nx2=' + repr(i/2),
                     'mesh/nx3=' + repr(i/2),
                     'meshblock/nx1=' + repr(i/4),
                     'meshblock/nx2=' + repr(i/4),
                     'meshblock/nx3=' + repr(i/4),
                     'output2/dt=-1',
                     'time/tlim=1.0',
                     'problem/compute_error=true']
        athena.run('hydro/athinput.linear_wave3d', arguments)
    # entropy wave
    for i in (32, 64):
        arguments = ['time/ncycle_out=0',
                     'problem/wave_flag=3',
                     'problem/vflow=1.0',
                     'mesh/refinement=static',
                     'mesh/nx1=' + repr(i),
                     'mesh/nx2=' + repr(i/2),
                     'mesh/nx3=' + repr(i/2),
                     'meshblock/nx1=' + repr(i/4),
                     'meshblock/nx2=' + repr(i/4),
                     'meshblock/nx3=' + repr(i/4),
                     'output2/dt=-1',
                     'time/tlim=1.0',
                     'problem/compute_error=true']
        athena.run('hydro/athinput.linear_wave3d', arguments)
    # L/R-going sound wave, no SMR
    for w in (0, 4):
        arguments = ['time/ncycle_out=0',
                     'problem/wave_flag=' + repr(w),
                     'output2/dt=-1', 'time/tlim=1.0', 'problem/compute_error=true']
        athena.run('hydro/athinput.linear_wave3d', arguments)


# Analyze outputs
def analyze():
    # read data from error file
    filename = 'bin/linearwave-errors.dat'
    data = []
    with open(filename, 'r') as f:
        raw_data = f.readlines()
        for line in raw_data:
            if line.split()[0][0] == '#':
                continue
            data.append([float(val) for val in line.split()])

    # check absolute error and convergence of all three waves
    if data[1][4] > 3.7e-8:
        print("error in L-going sound wave too large", data[1][4])
        return False
    if data[1][4]/data[0][4] > 0.325:
        print("not converging for L-going sound wave", data[0][4], data[1][4])
        return False

    if data[3][4] > 2.7e-8:
        print("error in entropy wave too large", data[3][4])
        return False
    if data[3][4]/data[2][4] > 0.33:
        print("not converging for entropy wave", data[2][4], data[3][4])
        return False

    # check error identical for waves in each direction
    if data[4][4] != data[5][4]:
        print("error in L/R-going sound waves not equal", data[4][4], data[5][4])
        return False

    return True
