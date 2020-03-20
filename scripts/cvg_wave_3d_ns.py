#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 ,-*
(_) Created on <Sa Okt 19 2019> @ 18:44:15

@authors: Boris Daszuta
@function: Convergence testing

"""
import utils.asset_tools as uat
import numpy as np
import matplotlib.pyplot as plt

###############################################################################
verbose = False
N_nx1 = 2 ** np.arange(5, 6)
cfl_number = 0.1
tlim = 10
integrator = 'rk4'
NGHOST = 2
num_threads = 4

# use threadnumber to control meshblock partition
meshblock_thread_partition = True

compile_flags='''--prob=cvg_wave_3d_trig -w --nghost={NGHOST}
--cxx=g++ -omp
-hdf5 -h5double --hdf5_path=$usr'''.format(NGHOST=NGHOST)

compile_flags='''--prob=cvg_wave_3d_trig -w --nghost={NGHOST}
--cxx=g++ -omp'''.format(NGHOST=NGHOST)

# nested structure within the output base directory
dir_outputs_nested_structure = ('wave_3d', )

###############################################################################

if meshblock_thread_partition:
    if np.remainder(N_nx1, num_threads).sum() != 0:
        raise ValueError("num_threads must divide for meshblock partition")

# for storing data
extracted_data = []

# parse compile flags
pcf = uat.problem_compile_flags(compile_flags=compile_flags)

for ix_nx1 in range(len(N_nx1)):

    #
    # Construct problem specification [this controls the .inp generator]
    #
    ps = uat.problem_specification()

    # <mesh>
    ps.set_value(section="mesh",
                 num_threads=num_threads,
                 refinement="static",
                 nx1=N_nx1[ix_nx1],
                 x1min=-1.0, x1max=1.0,
                 ix1_bc="periodic", ox1_bc="periodic",
                 nx2=N_nx1[ix_nx1],
                 x2min=-1.0, x2max=1.0,
                 ix2_bc="periodic", ox2_bc="periodic",
                 nx3=N_nx1[ix_nx1],
                 x3min=-1.0, x3max=1.0,
                 ix3_bc="periodic", ox3_bc="periodic")

    # <meshblock>
    mb_nx1 = N_nx1[ix_nx1]
    if meshblock_thread_partition:
        mb_nx1 = mb_nx1 // num_threads
    ps.set_value(section="meshblock",
                 nx1=mb_nx1, nx2=mb_nx1, nx3=mb_nx1)

    # <time>
    ps.set_value(section="time",
                 cfl_number=cfl_number,
                 tlim=tlim, integrator=integrator)

    # add <output1>
    # ps.add_section(
    #     section="output1",
    #     file_type="hdf5", variable="wave", dt=0.1,
    #     data_format="%.16e")

    # # add <output2> for restarts
    # ps.add_section(section="output2", file_type="rst", dt=0.1)

    # add <wave>
    ps.add_section(section="wave", c=1.0, use_Sommerfeld="false")

    # add <refinement1> section
    # ps.add_section(
    #     section="refinement1",
    #     x1min=0.0, x1max=1.0,
    #     x2min=-1.0, x2max=1.0,
    #     x3min=-1.0, x3max=1.0,
    #     level=5)

    # <job>
    ps.set_value(section="job",
                 problem_id="cvg_wave",
                 comment="dummy_comment",
                 compile_flags=pcf.compile_flags)


    #
    # we now init. the IO / compiler handler
    #
    pIOCh = uat.problem_IOC_handler(
        problem_specification=ps,
        problem_compile_flags=pcf,
        dir_outputs_nested_structure=dir_outputs_nested_structure,
        verbose=verbose,
        make_threads=num_threads)


    # monolithic directory structure / prob inp / compiler
    pIOCh.prepare_athena()

    # perform run if required
    pIOCh.execute_local()

#
# :D
#
