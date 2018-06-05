# Regression test based on MHD circularly polarized Alfven wave convergence problem
#
# Runs an isothermal cpaw convergence test in 2D including SMR and checks L1 errors (which
# are computed by the executable automatically and stored in the temporary file
# linearwave_errors.dat)

# Modules
import scripts.utils.athena as athena


# Prepare Athena++
def prepare(**kwargs):
    athena.configure('b', prob='cpaw', eos='isothermal', flux='hlld', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    # run R-going wave at two resolutions
    for i in (128, 256):
        arguments = ['time/ncycle_out=0',
                     'mesh/refinement=static',
                     'mesh/nx1=' + repr(i), 'mesh/nx2=' + repr(i/2),
                     'meshblock/nx1=' + repr(i/4), 'meshblock/nx2=' + repr(i/8),
                     'output2/dt=-1', 'time/tlim=1.0', 'problem/compute_error=true']
        athena.run('mhd/athinput.cpaw2d', arguments)
    # run L-going wave
    arguments = [
        'time/ncycle_out=0',
        'mesh/refinement=static',
        'mesh/nx1=256',
        'mesh/nx2=128',
        'meshblock/nx1=64',
        'meshblock/nx2=32',
        'output2/dt=-1',
        'time/tlim=1.0',
        'problem/compute_error=true',
        'problem/dir=2']
    athena.run('mhd/athinput.cpaw2d', arguments)


# Analyze outputs
def analyze():
    # read data from error file
    filename = 'bin/cpaw-errors.dat'
    data = []
    with open(filename, 'r') as f:
        raw_data = f.readlines()
        for line in raw_data:
            if line.split()[0][0] == '#':
                continue
            data.append([float(val) for val in line.split()])

    print(data[0][4])
    print(data[1][4])
    print(data[2][4])

    # check absolute error and convergence
    if data[1][4] > 2.0e-4:
        print("error in L-going fast wave too large", data[1][4])
        return False
    if data[1][4]/data[0][4] > 0.3:
        print("not converging for L-going fast wave", data[0][4], data[1][4])
        return False

    # check error identical for waves in each direction
    if abs(data[2][4] - data[1][4]) > 2.0e-6:
        print("error in L/R-going Alfven waves not equal", data[2][4], data[0][4])
        return False

    return True
