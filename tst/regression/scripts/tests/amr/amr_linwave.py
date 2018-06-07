# Regression test based on Newtonian 2D MHD linear wave test problem with AMR
#
# Runs a 2D linear wave test with AMR, using a refinement condition that tracks the
# velocity maxima.  Then checks L1 and L_infty (max) error.  This test is very sensitive
# to finding errors in AMR prolongation/restriction/boundaries

# Modules
import scripts.utils.athena as athena


# Prepare Athena++
def prepare(**kwargs):
    athena.configure('b',
                     prob='linear_wave',
                     coord='cartesian',
                     flux='hlld', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    # L-going fast wave (set by default in input)
    arguments = ['time/ncycle_out=10', 'output1/dt=-1']
    athena.run('mhd/athinput.linear_wave2d_amr', arguments)


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

    if data[0][4] > 2.0e-8:
        print("RMS error in L-going fast wave too large", data[0][4])
        return False
    if data[0][13] > 5.5:
        print("maximum relative error in L-going fast wave too large", data[0][13])
        return False

    return True
