### Input File Format
Athena++ requires an input file containing runtime parameters. Usually this file is given the name `athinput.<problem-name>`, where `<problem-name>` is an identifier of the problem. Sample input files are provided in `inputs/`.

Within the input file, parameters are grouped into named blocks, with the name of each block appearing on a single line within angle brackets, for example
```
    <time>
    cfl_number = 0.4       # The Courant, Friedrichs, & Lewy (CFL) Number
    nlim       = -1        # cycle limit
    tlim       = 1.0       # time limit
```
Block names must always appear in angle brackets on a separate line (blank lines above and below the block names are not required, but can be used for clarity).

Below each block name is a list of parameters, with syntax
```
    parameter = value [# comments]
```
White space after the parameter name, after the `=`, and before the `#` character is ignored. Everything after (and including) the `#` character is also ignored. Only one parameter value can appear per line. Comment lines (i.e. lines beginning with `#`) are allowed for documentation purposes. Both block names and parameter names are case sensitive.

### List of Input Blocks and Parameters
* `<comment>` (optional) : for documentation purposes
  * `problem` (optional) : problem description
  * `reference` (optional) : journal reference (if any)
  * `configure` (optional) : suggested configuration parameters
* `<job>` (**mandatory**)
  * `problem_id` (**mandatory**) : used in the output file names
* `<output[n]>` (optional) : output information (`[n]` is an integer)
  * `file_type` (**mandatory**) : file type (`vtk`, `hdf5`, `rst`, etc.; see [[Outputs]])
  * `dt` (**mandatory**) : output interval in computing time
  * `variable` (depends) : variable(s) to be output (see [[Outputs]])
  * `data_format` (optional; used only for `.tab` and `.hst`) : format specifier string used for writing data (e.g. `%12.5e`)
  * `id` (optional) : output ID used in file names (default = `out[n]`)
  * `x?_slice` (optional) :  sliced output in orthogonal directions at the specified position (`? = 1, 2, or 3`)
  * `ghost_zones` (optional) : include boundary ghost cells in output (default = `false`)
  * `cartesian_vector` (optional) : add vector variables converted into Cartesian coordinates (default = `false`)
  * `next_time` (optional)  : override the time for the end first output interval `dt`
* `<time>` (**mandatory**) : the CFL number and limit of the simulation
  * `cfl_number` (**mandatory**) : the Courant, Friedrichs, & Lewy (CFL) Number
  * `nlim` (**mandatory**) : time step limit (-1 = infinity)
  * `tlim` (**mandatory**) : time limit in computing time
  * `start_time` (optional) : time at the beginning of new simulation (default = 0)
  * `integrator` (optional) : time-integration scheme. Choices:
    * `vl2` (*default*) : second-order accurate van Leer predictor-corrector scheme
    * `rk2` : second-order accurate Runge-Kutta/Heun's method
    * `rk3` : third-order accurate Strong Stability Preserving (SSP) variant
    * `rk4` : fourth-order accurate, four-stage, 2 storage register method `RK4()4[2S]` from Table 2 of Ketcheson (2010)
    * `ssprk5_4` : fourth-order accurate, five-stage, 3 storage register, Strong Stability Preserving SSPRK(5,4) method
  * `xorder` (optional) : method for spatial reconstruction. Choices:
    * `2` (*default*) : Piecewise Linear Method (PLM) applied to primitive variables
    * `2c` : PLM applied to characteristic variables
    * `3`: Piecewise Parabolic Method (PPM) applied to primitive variables
    * `3c`: PPM applied to characteristic variables
    * `4` or `4c`: Enable transverse terms. PPM applied to primitive or characteristic variables, respectively. See [[High-Order Methods]] for more details. 
  * `ncycle_out` (optional) : interval for writing summary info to stdout. Choices:
    * `1` (*default*) : write out every cycle
    * `0` : suppresses all except final cycle summary info
    * `n` = positive integer : write out once per `n` cycles
  * `dt_diagnostics` (optional) : output extra timestep diagnostic info to stdout (at every `ncycle_out` interval). Choices:
    * `-1` (*default*) : disabled
    * `0` : list `dt_hyperbolic`, `dt_parabolic`, and `dt_user` (if active) and their respective ratios to the overall main integrator `dt`. If the code was configured with [[Super-Time-Stepping]], the following metric is also output: `stage_ratio = (dt*nstages)/(dt_parabolic*(nstages + nstages_sts))`, where `nstages` corresponds to the main multistage explicit time-integrator. 
    * `n` = positive integer : same behavior as `0`, but output extra diagnostics at every STS stage for proof-of-progress if `-sts` was used during [[Configuring]].
  * `sts_max_dt_ratio` (optional) : prevent [[Super-Time-Stepping]] from operating if the ratio of overall `dt` to `dt_parabolic` exceeds this value (-1 = infinity = default)
* `<mesh>` (**mandatory**) : grid configuration
  * `nx1,nx2,nx3` (**mandatory**) : the number of cells in the x1, x2, x3 directions, respectively (`nx3=1` means 2D, `nx2=1` implies 1D)
  * `x1min,x1max`, etc. (**mandatory**) : positions of the minimum and maximum surfaces (i.e., the box size)
  * `x1rat`, etc. (optional) : size ratio between neighboring cells in the given direction (see [[Coordinate Systems and Meshes]]). Choices:
    * `1.0` (*default*) : uniform coordinate spacing
    * `-1.0` : user-defined mesh generation function
    * `r`= positive Real  : size ratio fixed to `r`. Solver will raise a warning if `0.9 < r < 1.1` due to a possible loss in accuracy.
  * `ix1_bc,ox1_bc`, etc. (**mandatory**) : boundary conditions (see [[Boundary Conditions]])
  * `num_threads` (optional) : maximum number of OpenMP threads (default = 1)
  * `refinement` (optional) : enabling adaptive mesh refinement (default = none, see [[Adaptive Mesh Refinement]])
  * `numlevel` (optional) : the number of AMR refinement levels (default = 1, see [[Adaptive Mesh Refinement]])
  * `derefine_count` (optional) : the number of timesteps required before derefinement (default = 1, see [[Adaptive Mesh Refinement]])
* `<meshblock>` (optional) : domain decomposition unit (see [[Using MPI and OpenMP]])
  * `nx1,nx2,nx3` (mandatory) : the number of the cells per MeshBlock (decomposition unit) in x1, x2, x3, respectively
* `<refinement[n]>` (optional) : (static) refinement regions (`[n]` is an integer, see [[Static Mesh Refinement]])
  * `x1min,x1max,` etc. (mandatory) : positions of the refined regions
  * `level` (mandatory) : refinement level (root level = 0)
* `<loadbalancing>` (optional) : parameters for load balancing (see [[Load Balancing]])
  * `balancer` (optional) : load balancing method (default/manual/automatic)
  * `interval` (optional) : interval between load balancing (default = 10)
  * `tolerance` (optional) : acceptable load imbalance (default = 0.5 = 50%)

* `<hydro>` (**mandatory**) : parameters of hydrodynamics
  * `iso_sound_speed` (**mandatory** for isothermal EOS) : sounds speed in the isothermal EOS
  * `gamma` (**mandatory** for adiabatic EOS) : adiabatic index
  * `active` (optional) : `true`= dynamically evolve the fluid, `background`= the fluid is defined but the hydrodynamic variables are constant in time, `false` (under construction) = no hydrodynamic variables are defined across the domain (default = `true`)
  * `dfloor, pfloor` (optional) : density and pressure floors
  * `gamma_max` (optional) : maximum Lorentz factor in SR and GR (default = 1000)
  * `rho_min, rho_pow, u_min, u_pow` (optional) : additional controls for density and internal energy floors used in GR
  * `rho_pmag_min, u_pmag_min` (optional) : limits on plasma sigma and beta used in MHD in SR and GR
* `<coord>` (depends) : parameters for GR coordinate system
  * `m` (depends) : mass of black hole
  * `a` (depends) : spin of black hole
* `<problem>` (optional) : problem-specific parameters
  * `GM` (optional) : enables point source gravity; gravitational constant x mass of the point source - this works only in 3D spherical polar coordinates or 2D cylindrical coordinates

#### Changes from old Athena
Although the structure of the input file is almost the same, the names of the blocks and parameters are different. The input files used with Athena must be rewritten accordingly.
