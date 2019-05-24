The executable binary file `bin/athena` is created after compilation. In order to start a simulation, you need to specify an input file (see also [[The Input File]]) corresponding to the problem generator specified in configuration.

    > athena -i athinput.example

By default, all the output files are created in the current working directory. The code will output simulation progress and some diagnosis messages to stdout.

For parallel simulations, see [[Using MPI and OpenMP]].

### Command-Line Options
A variety of command-line options have been implemented in Athena. The list of the options as well as the code configuration is given by the `-h` switch:

    > athena -h

* `-h` : show the help message
* `-i <file>` : specify an input file (default = athinput)
* `-r <file>` : restart from a restarting file
* `-d <directory>` : specify a run dir (default = current dir)
* `-n` : test an input file and quit
* `-c` : show configuration and quit
* `-m <nproc>` : output mesh structure and quit (see [[Static Mesh Refinement]])
* `-t hh:mm:ss` : set wall time limit for final output

The `-t` option is convenient for running simulations on supercomputing systems (see also [[Restart File|Outputs#restart-file]]).  Allow some time for the final output (a restart output must be specified in the input file/options) as the code makes the final output after the specified time limit is reached. Alternatively, you can send a signal (SIGTERM or SIGINT) to stop Athena++ gently. Some job schedulers such as Moab, SLURM and NQSII allow users to set the time interval when a SIGTERM signal is sent to a job before the job is actually killed. For details, consult the documentation of your system.

### Overriding Input Parameters
Some input parameters can be overridden by command-line parameters with the `<parameter_block>/<parameter>=<value>` format. For example, if you want to extend the simulation time limit when restarting,

    > athena -r example.out1.00010.rst time/tlim=100

This feature must be used with caution, because it can cause inconsistency. Also, you cannot change some parameters, such as the sizes of Mesh and MeshBlock. Changing the time or number limit of the simulation or the output parameters is generally safe.
