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
N_nx1 = 2 ** np.arange(4, 14)
cfl_number = 0.5
tlim = 1
integrator = 'rk4'
NGHOST = 4
num_threads = 4

# use threadnumber to control meshblock partition
meshblock_thread_partition = False

compile_flags='''--prob=cvg_wave_1d_trig -w --nghost={NGHOST}
--cxx=g++ -omp
-hdf5 -h5double --hdf5_path=$usr'''.format(NGHOST=NGHOST)

# nested structure within the output base directory
dir_outputs_nested_structure = ('wave_1d', )

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
                 nx2=1,
                 x2min=-1.0, x2max=1.0,
                 ix2_bc="periodic", ox2_bc="periodic",
                 nx3=1,
                 x3min=-1.0, x3max=1.0,
                 ix3_bc="periodic", ox3_bc="periodic")

    # <meshblock>
    mb_nx1 = N_nx1[ix_nx1]
    if meshblock_thread_partition:
        mb_nx1 = mb_nx1 // num_threads
    ps.set_value(section="meshblock",
                 nx1=mb_nx1, nx2=1, nx3=1)

    # <time>
    ps.set_value(section="time",
                 cfl_number=cfl_number,
                 tlim=tlim, integrator=integrator)

    # add <output1>
    ps.add_section(
        section="output1",
        file_type="hdf5", variable="wave", dt=0.1,
        data_format="%.16e")

    # add <output2> for restarts
    ps.add_section(section="output2", file_type="rst", dt=0.1)

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


###############################################################################
# parse all directories for known data
###############################################################################

# need to extract known times and file-names

# we only consider single-block, rk4
filter_problem_input = (
    lambda _: _['mesh']['nx1'] != _['meshblock']['nx1'],
#    lambda _: _['time']['integrator'] != 'rk4',
)
filter_problem_input = None

# throw away data at some times
filter_dataset = (
    #    lambda _: (_['Time'] < 0.5) or (_['Time'] > 0.8),
    lambda _: (_['Time'] != 0.5) and (_['Time'] != 1.0),
)

# format of returned info
key_format = (('Time',),
              ('mesh', 'nx1'),
              ('mesh', 'nx1'),
              ('meshblock', 'num_threads'),
              ('compile_flags', '--nghost'),
              ('time', 'integrator'),
)

data_file_info = pIOCh.dataset_parser(
    filter_problem_input=filter_problem_input,
    filter_dataset=filter_dataset,
    key_format=key_format)

data = pIOCh.data_load(data_file_info, quantities=['wError'])



# compute L^2, L_\infty norms as we have exact and numerical solution

def parse_errors():
    # compute max_abs_err
    # approx. discrete L^2 norm
    known_times = {}
    for ix_nx1, el in enumerate(extracted_data):
        for ix_t, t_el in enumerate(el):
            T = t_el['Time']
            ixs = (ix_nx1, ix_t)
            if T in known_times:
                known_times[T].append(ixs)
            else:
                known_times[T] = [ixs]

    err_info = {}
    for k, v in known_times.items():
        c_len = len(v)  # total number of index tuples with vals at this time
        err_arr = np.zeros((c_len, 2))
        for ix_v, ix_tup in enumerate(v):
            nx1 = N_nx1[ix_tup[0]]
            err_arr[ix_v, 0] = nx1
            err_arr[ix_v, 1] = np.abs(
                extracted_data[ix_tup[0]][ix_tup[1]]['wError']).max()

        err_info[k] = err_arr
    return err_info


# out = parse_errors()
# plt.loglog(out[0.5][:,0], out[0.5][:,1], '-xr')
# plt.loglog(out[1][:,0], out[1][:,1], '-ob')

def data_extract(data, Time=None, nghost=None, integrator=None,
                 meshblock_thread_partition=False):
    nx1_lst = []
    err_lst = []

    for k, v in data.items():
        _Time, _nx1, _nx1_mb, _num_threads, _nghost, _integrator = k

        flag_1 = ((str(Time) == _Time) and
                  (str(nghost) == _nghost) and
                  (str(integrator) == _integrator))
        flag_2 = (int(_num_threads) * int(_nx1_mb) == int(_nx1))

        if meshblock_thread_partition:
            flag = flag_1 and flag_2
        else:
            flag = flag_1

        if flag:
            nx1_lst.append(int(_nx1))
            err_lst.append(abs(v['wError']).max())

    arr_nx1 = np.array(nx1_lst)
    ix_sort = arr_nx1.argsort()
    arr_nx1 = arr_nx1[ix_sort]
    arr_err = np.array(err_lst)[ix_sort]

    return arr_nx1, arr_err


def plot_max_abs(data):
    plt.figure(1)

    nx1, err = data_extract(data, Time='0.5', nghost=2, integrator='rk4',
                            meshblock_thread_partition=False)
    plt.loglog(nx1, err, 'or')

    nx1, err = data_extract(data, Time='0.5', nghost=4, integrator='rk4',
                            meshblock_thread_partition=False)
    plt.loglog(nx1, err, 'xr')


    nx1, err = data_extract(data, Time='0.5', nghost=2, integrator='rk3',
                            meshblock_thread_partition=False)
    plt.loglog(nx1, err, 'og')

    nx1, err = data_extract(data, Time='0.5', nghost=4, integrator='rk3',
                            meshblock_thread_partition=False)
    plt.loglog(nx1, err, 'xg')



    nx1, err = data_extract(data, Time='0.5', nghost=2, integrator='rk4',
                            meshblock_thread_partition=True)
    plt.loglog(nx1, err, '-or')

    nx1, err = data_extract(data, Time='0.5', nghost=4, integrator='rk4',
                            meshblock_thread_partition=True)
    plt.loglog(nx1, err, '-xr')


    nx1, err = data_extract(data, Time='0.5', nghost=2, integrator='rk3',
                            meshblock_thread_partition=True)
    plt.loglog(nx1, err, '-og')

    nx1, err = data_extract(data, Time='0.5', nghost=4, integrator='rk3',
                            meshblock_thread_partition=True)
    plt.loglog(nx1, err, '-xg')


    # fit to internal points
    # nx1, err = data_extract(data, Time='0.5', nghost=2, integrator='rk4')

    # ln_nx1, ln_err = np.log(nx1)[2:-3], np.log(err)[2:-3]

    # pf = np.polyfit(ln_nx1, ln_err, 1)

    # superpose -4 fit


    plt.show()

#
# :D
#
