# Regression test based on Newtonian hydro linear wave convergence problem
#
# Runs a linear wave convergence test in 3D including SMR and checks L1 errors (which
# are computed by the executable automatically and stored in the temporary file
# linearwave_errors.dat)

# Modules
import os
import scripts.utils.athena as athena


# Prepare Athena++
def prepare(**kwargs):
    athena.configure('mpi', 'fft',
                     prob='jeans',
                     grav='fft', **kwargs)
    athena.make()
    os.system('mv bin/athena bin/athena_mpi_fft')

    # athena.configure('mpi',
    #                  prob='jeans',
    #                  grav='mg', **kwargs)
    # athena.make()
    # os.system('mv bin/athena bin/athena_mpi_mg')

    athena.configure('fft',
                     prob='jeans',
                     grav='fft',
                     **kwargs)
    athena.make()
    os.system('mv bin/athena bin/athena_fft')

    # athena.configure(prob='jeans',
    #                  grav='mg',
    #                  **kwargs)
    # athena.make()


# Run Athena++
def run(**kwargs):
    arguments = ['time/ncycle_out=10',
                 'mesh/nx1=64', 'mesh/nx2=32', 'mesh/nx3=32',
                 'meshblock/nx1=16',
                 'meshblock/nx2=16',
                 'meshblock/nx3=16',
                 'output2/dt=-1', 'time/tlim=1.0', 'problem/compute_error=true']
    # athena.run('hydro/athinput.jeans_3d', arguments)

    os.system('mv bin/athena_fft bin/athena')
    athena.run('hydro/athinput.jeans_3d', arguments)

    # os.system('mv bin/athena_mpi_mg bin/athena')
    # athena.mpirun(kwargs['mpirun_cmd'], kwargs['mpirun_opts'],
    #               1, 'hydro/athinput.jeans_3d', arguments)
    # athena.mpirun(kwargs['mpirun_cmd'], kwargs['mpirun_opts'],
    #               2, 'hydro/athinput.jeans_3d', arguments)
    # athena.mpirun(kwargs['mpirun_cmd'], kwargs['mpirun_opts'],
    #               4, 'hydro/athinput.jeans_3d', arguments)

    os.system('mv bin/athena_mpi_fft bin/athena')
    athena.mpirun(kwargs['mpirun_cmd'], kwargs['mpirun_opts'],
                  1, 'hydro/athinput.jeans_3d', arguments)
    athena.mpirun(kwargs['mpirun_cmd'], kwargs['mpirun_opts'],
                  2, 'hydro/athinput.jeans_3d', arguments)
    athena.mpirun(kwargs['mpirun_cmd'], kwargs['mpirun_opts'],
                  4, 'hydro/athinput.jeans_3d', arguments)


# Analyze outputs
def analyze():
    # read data from error file
    filename = 'bin/jeans-errors.dat'
    data = []
    with open(filename, 'r') as f:
        raw_data = f.readlines()
        for line in raw_data:
            if line.split()[0][0] == '#':
                continue
            data.append([float(val) for val in line.split()])

    print(data[0][4], data[1][4])
    print(data[2][4], data[3][4])  # , data[4][4])
    # print(data[5][4], data[6][4], data[7][4])

    # check errors between runs w/wo MPI and different numbers of cores
    # if data[0][4] != data[2][4]:
    #     print(
    #         "Linear wave error with one core w/wo MPI not identical for MG gravity",
    #         data[0][4],
    #         data[2][4])
    #     return False
    # if data[0][4] > 1.e-7:
    #     print("Linear wave error is too large for MG gravity", data[0][4])
    #     return False
    if data[0][4] > 1.e-7:
        print("Linear wave error is too large for FFT gravity", data[0][4])
        return False
    if data[0][4] != data[1][4]:
        print(
            "Linear wave error with one core w/wo MPI not identical for FFT gravity",
            data[0][4],
            data[1][4])
        return False
    # if abs(data[3][4]-data[0][4]) > 5.0e-4:
    #     print(
    #         "Linear wave error between 2 and 1 cores too large for MG gravity",
    #         data[3][4],
    #         data[0][4])
    #     return False
    # if abs(data[4][4]-data[0][4]) > 5.0e-4:
    #     print(
    #         "Linear wave error between 4 and 1 cores too large for MG gravity",
    #         data[4][4],
    #         data[0][4])
    #     return False
    if abs(data[2][4]-data[0][4]) > 5.0e-4:
        print(
            "Linear wave error between 2 and 1 cores too large for FFT gravity",
            data[2][4],
            data[0][4])
        return False
    if abs(data[3][4]-data[0][4]) > 5.0e-4:
        print(
            "Linear wave error between 4 and 1 cores too large for FFT gravity",
            data[3][4],
            data[0][4])
        return False

    return True
