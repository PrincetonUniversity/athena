# Regression test based on cosmic ray diffusion test problem

# Modules
import sys
import scripts.utils.athena as athena
sys.path.insert(0, '../../vis/python')

# Prepare Athena++


def prepare(**kwargs):
    athena.configure('cr',
                     prob='cr_diffusion',
                     coord='cartesian',
                     flux='hllc')
    athena.make()

# Run Athena++


def run(**kwargs):
    # case 1: static diffusion along x direction
    arguments = ['mesh/nx1=256', 'mesh/ix1_bc=outflow', 'mesh/ox1_bc=outflow',
                 'mesh/nx2=4', 'mesh/ix2_bc=periodic', 'mesh/ox2_bc=periodic',
                 'meshblock/nx1=32', 'meshblock/nx2=4', 'problem/direction=0',
                 'problem/v0=0', 'time/ncycle_out=100']
    athena.run('cosmic_ray/athinput.cr_diffusion', arguments)

    # case 2: dynamic static diffusion along x direction
    arguments = ['mesh/nx1=256', 'mesh/ix1_bc=outflow', 'mesh/ox1_bc=outflow',
                 'mesh/nx2=4', 'mesh/ix2_bc=periodic', 'mesh/ox2_bc=periodic',
                 'meshblock/nx1=32', 'meshblock/nx2=4', 'problem/direction=0',
                 'problem/v0=1', 'time/ncycle_out=100']
    athena.run('cosmic_ray/athinput.cr_diffusion', arguments)

    # case 3: static diffusion along y direction
    arguments = ['mesh/nx1=4', 'mesh/ix1_bc=periodic', 'mesh/ox1_bc=periodic',
                 'mesh/nx2=256', 'mesh/ix2_bc=outflow', 'mesh/ox2_bc=outflow',
                 'meshblock/nx1=4', 'meshblock/nx2=32' 'problem/direction=0',
                 'problem/v0=0', 'time/ncycle_out=100']
    athena.run('cosmic_ray/athinput.cr_diffusion', arguments)

    # case 4: dynamic static diffusion along y direction
    arguments = ['mesh/nx1=4', 'mesh/ix1_bc=periodic', 'mesh/ox1_bc=periodic',
                 'mesh/nx2=256', 'mesh/ix2_bc=outflow', 'mesh/ox2_bc=outflow',
                 'meshblock/nx1=4', 'meshblock/nx2=32', 'problem/direction=1',
                 'problem/v0=1', 'time/ncycle_out=100']
    athena.run('cosmic_ray/athinput.cr_diffusion', arguments)

    # case 5: static diffusion along z direction
    arguments = ['mesh/nx1=4', 'mesh/ix1_bc=periodic', 'mesh/ox1_bc=periodic',
                 'mesh/nx2=4', 'mesh/ix2_bc=periodic', 'mesh/ox2_bc=periodic',
                 'mesh/nx3=256', 'mesh/ix3_bc=outflow', 'mesh/ox3_bc=outflow',
                 'meshblock/nx1=4', 'meshblock/nx2=4', 'meshblock/nx3=32',
                 'problem/direction=2', 'problem/v0=0', 'time/ncycle_out=100']
    athena.run('cosmic_ray/athinput.cr_diffusion', arguments)

    # case 6: dynamic static diffusion along z direction
    arguments = ['mesh/nx1=4', 'mesh/ix1_bc=periodic', 'mesh/ox1_bc=periodic',
                 'mesh/nx2=4', 'mesh/ix2_bc=periodic', 'mesh/ox2_bc=periodic',
                 'mesh/nx3=256', 'mesh/ix3_bc=outflow', 'mesh/ox3_bc=outflow',
                 'meshblock/nx1=4', 'meshblock/nx2=4', 'meshblock/nx3=32',
                 'problem/direction=2', 'problem/v0=1', 'time/ncycle_out=100']
    athena.run('cosmic_ray/athinput.cr_diffusion', arguments)


# Analyze outputs


def analyze():
    filename = 'bin/diffusion_error.dat'
    data = []
    with open(filename, 'r') as f:
        raw_data = f.readlines()
        for line in raw_data:
            if line.split()[0][0] == '#':
                continue
            data.append([float(val) for val in line.split()])

    # check absolute error and convergence of all three waves
    # regime1
    if data[0][8] > 8e-5:
        print("error in static diffusion along x: ", data[0][8])
        return False

        # regime2
    if data[1][8] > 9e-5:
        print("error in dynamic diffusion along x: ", data[1][8])
        return False

    if data[2][8] > 8e-5:
        print("error in static diffusion along y: ", data[2][8])
        return False

        # regime2
    if data[3][8] > 9e-5:
        print("error in dynamic diffusion along y: ", data[3][8])
        return False

    if data[4][8] > 8e-5:
        print("error in static diffusion along z: ", data[4][8])
        return False

        # regime2
    if data[5][8] > 9e-5:
        print("error in dynamic diffusion along z: ", data[5][8])
        return False

    return True
