# Regression test based on linear wave convergence problem
# for full radiation hydro equations, using the
# implicit radiation hydro module

# Modules
import sys
import scripts.utils.athena as athena
sys.path.insert(0, '../../vis/python')

# Prepare Athena++


def prepare(**kwargs):
    athena.configure('implicit_radiation',
                     prob='rad_linearwave',
                     coord='cartesian',
                     flux='hllc')
    athena.make()

# Run Athena++


def run(**kwargs):
    # case 1
    for i in (32, 64, 128, 256):
        arguments = ['problem/regime=1', 'radiation/prat=0.01', 'radiation/crat=10.0',
                     'radiation/error_limit=1.e-15',
                     'time/tlim=0.7745966144169111', 'problem/compute_error=true',
                     'mesh/nx1=' + repr(i), 'mesh/nx2=8', 'mesh/nx3=1',
                     'meshblock/nx1=' + repr(i/2),
                     'meshblock/nx2=8',
                     'meshblock/nx3=1',
                     'time/ncycle_out=100']
        athena.run('radiation/athinput.rad_linearwave', arguments)
#
#  bashcommand="cp bin/linearwave-errors.dat ~/linearwave-errors_regime1.dat"
#  os.system(bashcommand)

    # case 2
    for i in (32, 64, 128, 256):
        arguments = ['problem/regime=2', 'radiation/prat=0.01', 'radiation/crat=10.0',
                     'radiation/error_limit=1.e-15',
                     'time/tlim=0.7745911961524788', 'problem/compute_error=true',
                     'mesh/nx1=' + repr(i), 'mesh/nx2=8', 'mesh/nx3=1',
                     'meshblock/nx1=' + repr(i/2),
                     'meshblock/nx2=8',
                     'meshblock/nx3=1',
                     'time/ncycle_out=100']
        athena.run('radiation/athinput.rad_linearwave', arguments)

    # case 3
    for i in (32, 64, 128, 256):
        arguments = ['problem/regime=3', 'radiation/prat=0.01', 'radiation/crat=10.0',
                     'radiation/error_limit=1.e-15',
                     'time/tlim=0.7741319089038714', 'problem/compute_error=true',
                     'mesh/nx1=' + repr(i), 'mesh/nx2=8', 'mesh/nx3=1',
                     'meshblock/nx1=' + repr(i/2),
                     'meshblock/nx2=8',
                     'meshblock/nx3=1',
                     'time/ncycle_out=100']
        athena.run('radiation/athinput.rad_linearwave', arguments)

    # case 4
    for i in (32, 64, 128, 256):
        arguments = ['problem/regime=4', 'radiation/prat=0.01', 'radiation/crat=10.0',
                     'radiation/error_limit=1.e-15',
                     'time/tlim=0.7739734114015662', 'problem/compute_error=true',
                     'mesh/nx1=' + repr(i), 'mesh/nx2=8', 'mesh/nx3=1',
                     'meshblock/nx1=' + repr(i/2),
                     'meshblock/nx2=8',
                     'meshblock/nx3=1',
                     'time/ncycle_out=100']
        athena.run('radiation/athinput.rad_linearwave', arguments)

    # case 5
    for i in (32, 64, 128, 256):
        arguments = ['problem/regime=5', 'radiation/prat=0.01', 'radiation/crat=10.0',
                     'radiation/error_limit=1.e-15', 'radiation/taucell=15',
                     'time/tlim=0.7746132569285668', 'problem/compute_error=true',
                     'mesh/nx1=' + repr(i), 'mesh/nx2=8', 'mesh/nx3=1',
                     'meshblock/nx1=' + repr(i/2),
                     'meshblock/nx2=8',
                     'meshblock/nx3=1']
        athena.run('radiation/athinput.rad_linearwave', arguments)

    # case 6
    for i in (32, 64, 128, 256):
        arguments = ['problem/regime=6', 'radiation/prat=1.0', 'radiation/crat=10.0',
                     'radiation/error_limit=1.e-15',
                     'time/tlim=0.7747068066908429', 'problem/compute_error=true',
                     'mesh/nx1=' + repr(i), 'mesh/nx2=8', 'mesh/nx3=1',
                     'meshblock/nx1=' + repr(i/2),
                     'meshblock/nx2=8',
                     'meshblock/nx3=1',
                     'time/ncycle_out=100']
        athena.run('radiation/athinput.rad_linearwave', arguments)

    # case 7
    for i in (32, 64, 128, 256):
        arguments = ['problem/regime=7', 'radiation/prat=1.0', 'radiation/crat=10.0',
                     'radiation/error_limit=1.e-15',
                     'time/tlim=0.7856204848539599', 'problem/compute_error=true',
                     'mesh/nx1=' + repr(i), 'mesh/nx2=8', 'mesh/nx3=1',
                     'meshblock/nx1=' + repr(i/2),
                     'meshblock/nx2=8',
                     'meshblock/nx3=1',
                     'time/ncycle_out=100']
        athena.run('radiation/athinput.rad_linearwave', arguments)

    # case 8
    for i in (32, 64, 128, 256):
        arguments = ['problem/regime=8', 'radiation/prat=1.0', 'radiation/crat=10.0',
                     'radiation/error_limit=1.e-15',
                     'time/tlim=0.9819810298708518', 'problem/compute_error=true',
                     'mesh/nx1=' + repr(i), 'mesh/nx2=8', 'mesh/nx3=1',
                     'meshblock/nx1=' + repr(i/2),
                     'meshblock/nx2=8',
                     'meshblock/nx3=1',
                     'time/ncycle_out=100']
        athena.run('radiation/athinput.rad_linearwave', arguments)

#  bashcommand="cp bin/linearwave-errors.dat ~/linearwave-errors.dat"
#  os.system(bashcommand)

# Analyze outputs


def analyze():
    filename = 'bin/linearwave-errors.dat'
    data = []
    with open(filename, 'r') as f:
        raw_data = f.readlines()
        for line in raw_data:
            if line.split()[0][0] == '#':
                continue
            data.append([float(val) for val in line.split()])

    # check absolute error and convergence of all three waves
    # regime1
    if data[3][4] > 8.0e-10:
        print("error in regime 1: ", data[3][4])
        return False
    if data[3][4]/data[2][4] > 0.55:
        print("not converging for regime 1: ", data[3][4], data[2][4])
        return False

        # regime2
    if data[7][4] > 1.2e-9:
        print("error in regime 2: ", data[7][4])
        return False
    if data[7][4]/data[6][4] > 0.55:
        print("not converging for regime 2: ", data[7][4], data[6][4])
        return False

    # regime3
    if data[11][4] > 6.6e-9:
        print("error in regime 3: ", data[11][4])
        return False
    if data[11][4]/data[10][4] > 0.55:
        print("not converging for regime 3: ", data[11][4], data[10][4])
        return False

    # regime4
    if data[15][4] > 8.1e-9:
        print("error in regime 4: ", data[15][4])
        return False
    if data[15][4]/data[14][4] > 0.55:
        print("not converging for regime 4: ", data[15][4], data[14][4])
        return False

    # regime5
    if data[19][4] > 5.0e-9:
        print("error in regime 5: ", data[19][4])
        return False
    if data[19][4]/data[18][4] > 0.55:
        print("not converging for regime 5: ", data[19][4], data[18][4])
        return False

    # regime6
    if data[23][4] > 8.4e-10:
        print("error in regime 6: ", data[23][4])
        return False
    if data[23][4]/data[22][4] > 0.55:
        print("not converging for regime 6: ", data[23][4], data[22][4])
        return False

    # regime7
    if data[27][4] > 2.11e-9:
        print("error in regime 7: ", data[27][4])
        return False
    if data[27][4]/data[26][4] > 0.55:
        print("not converging for regime 7: ", data[27][4], data[26][4])
        return False

    # regime8
    if data[31][4] > 2.1e-9:
        print("error in regime 8: ", data[31][4])
        return False
    if data[31][4]/data[30][4] > 0.55:
        print("not converging for regime 8: ", data[31][4], data[30][4])
        return False

    return True
