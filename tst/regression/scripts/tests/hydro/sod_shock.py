# Regression test based on Newtonian hydro Sod shock tube problem
#
# Runs the Sod shock tube in x1, x2, and x3 directions successively, and checks errors
# against the analytic solution (which are computed by the executable automatically and
# stored in the temporary file shock_errors.dat)

# Modules
import scripts.utils.athena as athena


# Prepare Athena++
def prepare(**kwargs):
    athena.configure(prob='shock_tube', coord='cartesian', flux='hllc', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    # run in X1 direction
    for i in (128, 256):
        arguments = [
            'time/ncycle_out=0',
            'mesh/nx1=' + repr(i),
            'mesh/nx2=1',
            'mesh/nx3=1',
            'mesh/ix1_bc=outflow',
            'mesh/ox1_bc=outflow',
            'mesh/ix2_bc=periodic',
            'mesh/ox2_bc=periodic',
            'mesh/ix3_bc=periodic',
            'mesh/ox3_bc=periodic',
            'time/cfl_number=0.3',
            'problem/shock_dir=1',
            'problem/compute_error=true']
        athena.run('hydro/athinput.sod', arguments)
    # run in X2 direction
    for i in (128, 256):
        arguments = [
            'time/ncycle_out=0',
            'mesh/nx1=4',
            'mesh/nx2=' + repr(i),
            'mesh/nx3=1',
            'mesh/ix1_bc=periodic',
            'mesh/ox1_bc=periodic',
            'mesh/ix2_bc=outflow',
            'mesh/ox2_bc=outflow',
            'mesh/ix3_bc=periodic',
            'mesh/ox3_bc=periodic',
            'time/cfl_number=0.3',
            'problem/shock_dir=2',
            'problem/compute_error=true']
        athena.run('hydro/athinput.sod', arguments)
    # run in X3 direction
    for i in (128, 256):
        arguments = [
            'time/ncycle_out=0',
            'mesh/nx1=4',
            'mesh/nx2=4',
            'mesh/nx3=' + repr(i),
            'mesh/ix1_bc=periodic',
            'mesh/ox1_bc=periodic',
            'mesh/ix2_bc=periodic',
            'mesh/ox2_bc=periodic',
            'mesh/ix3_bc=outflow',
            'mesh/ox3_bc=outflow',
            'time/cfl_number=0.3',
            'problem/shock_dir=3',
            'problem/compute_error=true']
        athena.run('hydro/athinput.sod', arguments)


# Analyze outputs
def analyze():
    # read data from error file
    filename = 'bin/shock-errors.dat'
    data = []
    with open(filename, 'r') as f:
        raw_data = f.readlines()
        for line in raw_data:
            if line.split()[0][0] == '#':
                continue
            data.append([float(val) for val in line.split()])

    # check Ncycles same for each direction
    if data[1][3] != data[3][3]:
        print("Ncycles in x1 not equal to Ncycles in x2", data[1][3], data[3][3])
        return False
    if data[1][3] != data[5][3]:
        print("Ncycles in x1 not equal to Ncycles in x3", data[1][3], data[5][3])
        return False

    # check absolute error and convergence in x1
    if data[0][4] > 0.011:
        print("error in x1 too large", data[0][4])
        return False
    if data[1][4] / data[0][4] > 0.6:
        print("not converging in x1", data[0][4], data[1][4])
        return False

    # check absolute error and convergence in x2
    if data[2][4] > 0.011:
        print("error in x2 too large", data[2][4])
        return False
    if data[3][4] / data[2][4] > 0.6:
        print("not converging in x2", data[2][4], data[3][4])
        return False

    # check absolute error and convergence in x3
    if data[4][4] > 0.011:
        print("error in x3 too large", data[4][4])
        return False
    if data[5][4] / data[4][4] > 0.6:
        print("not converging in x3", data[4][4], data[5][4])
        return False

    return True
