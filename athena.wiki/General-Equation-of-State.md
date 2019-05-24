### General Equation of State

In addition to isothermal and adiabatic equations of state (EOS), Athena++ has the capability of utilizing more general/realistic EOS (MHD is not currently supported). Currently, two general EOS options have been implemented `general/hydrogen` and `general/eos_table`. These options are specified at configuration, e.g.:

    > python configure.py --prob shock_tube --eos general/hydrogen

### General EOS framework

The file `src/eos/general/general_hydro.cpp` establishes the framework for all general EOS (i.e. all the `general\<...>` EOS options). Each new general EOS must define the following functions:
* `Real EquationOfState::PresFromRhoEg(Real rho, Real egas)` gas pressure as a function of mass density and internal energy density
* `Real EquationOfState::EgasFromRhoP(Real rho, Real pres)` internal energy density as a function of mass density and gas pressure
* `Real EquationOfState::AsqFromRhoP(Real rho, Real pres)` adiabatic sound speed squared as a function of mass density and gas pressure
* `Real EquationOfState::RiemannAsq(Real rho, Real hint)` adiabatic sound speed squared as a function of mass density and specific internal enthalpy ((egas + pres)/rho)

The function `PresFromRhoEg` is used to convert conservatives to primitives, and `EgasFromRhoP` is used to convert primitives to conservative. Problem generators intended to be used with a general EOS must utilized these EOS functions for the corresponding conversions.

### EOS Table

The `general/eos_table` option allows for the user to utilize a specified EOS table. This requires four data fields to implement the above EOS functions. These data fields are:
1. Log(p/e(e/rho,rho))
2. Log(e/p(p/rho,rho))
3. Log(asq*rho/p(p/rho,rho))
4. Log(asq*rho/h(h/rho,rho))

where "Log" is implicitly base 10 and the second independent variable (rho=density) is the fastest indexing variable. The variables p, e, asq are gas pressure, internal energy density and adiabatic sound speed squared respectively, and h=e+p. The tables need to be rectangular in log space. The file containing the EOS data is specified in the input file along with some other options:

    <hydro>
    eos_file_name   = gamma_is_1.100.tab # Specify EOS table filename
    eos_file_type   = ascii              # Specify EOS table file type [ascii,binary,hdf5]
    eos_read_ratios = true               # Table file specifies ratios between different
                                         # variables, e.g e/rho vs p/rho (default=true)
    eos_rho_unit    = 1.0                # Table unit/sim unit for mass density (default=1.0)
    eos_egas_unit   = 1.0                # Table unit/sim unit for internal energy density 
                                         # (default=1.0)

There are three file types that can be read in: ascii table, raw binary, and hdf5. The `eos_read_ratios` option allows for the different first independent variables (e/rho, p/rho, h/rho) to have different normalisations compared to each other. If the default/specified range for the first independent variable is `[min, max]` then the range for field `i` is `[ratio[i]*min, ratio[i]*max]`. The primary inclusion for this option is to allow the different fields to cover the same parameter space for a table that is asymptotically ideal. When `eos_read_ratios` is `false` then the ratios are set to 1, making the min/max values for `e/rho`, `p/rho`, `h/rho` the same.

Here is a simple ascii (`ascii`) table which implements an ideal/adiabatic EOS with gamma=1.1:

    # Entries must be space separated.
    #  n_var,  n_espec, n_rho
    # (fields) (rows)  (columns)
    4 2 3
    # Log espec lim (specific internal energy e/rho)
    -1.0000e+01 2.0000e+01
    # Log rho lim
    -2.4000e+01 4.0000e+00
    # Ratios = 1, eint/pres, eint/pres, eint/h
    # This line is required iff EOS_read_ratios
    1.0000e+00 1.0000e+01 1.0000e+01 9.0909e-01
    # Log p/e(e/rho,rho)
    -1.0000e+00 -1.0000e+00 -1.0000e+00
    -1.0000e+00 -1.0000e+00 -1.0000e+00
    # Log e/p(p/rho,rho)
    1.0000e+00 1.0000e+00 1.0000e+00
    1.0000e+00 1.0000e+00 1.0000e+00
    # Log asq*rho/p(p/rho,rho)
    4.1393e-02 4.1393e-02 4.1393e-02
    4.1393e-02 4.1393e-02 4.1393e-02
    # Log asq*rho/h(h/rho,rho)
    -1.0000e+00 -1.0000e+00 -1.0000e+00
    -1.0000e+00 -1.0000e+00 -1.0000e+00

The `binary` file type is a raw binary file where the size (n_var, n_espec, n_rho) is given in 32 bit integers, and all floats have the same precision as `Real`. Data is written in the order as the above ascii table example.

The `hdf5` option requires HDF5 to be enabled (`-hdf5`). By default, the code assumes that the following datasets are defined:
* `LogDensLim` = [Log min(rho), Log max(rho)]
* `LogEspecLim` = [Log min(e/rho), Log max(e/rho)]
* `ratios` shape = [4] (only read if `EOS_read_ratios`)
* `p/e(e/rho,rho)` shape = [n_espec, n_rho], Log p/e(e/rho,rho)
* `e/p(p/rho,rho)` shape = [n_espec, n_rho], Log e/p(p/rho,rho)
* `asq*rho/p(p/rho,rho)` shape = [n_espec, n_rho], Log asq*rho/p(p/rho,rho)
* `asq*rho/h(h/rho,rho)` shape = [n_espec, n_rho], Log asq*rho/h(h/rho,rho)

### Hydrogen EOS

This EOS assumes that there are three constituents in the plasma: neutral hydrogen, ionized hydrogen, and free electrons. This EOS also assumes that the relative density of these species are governed by local thermal equilibrium, i.e. the Saha equation describing hydrogen recombination. Note only the ground state of neutral hydrogen is taken into account for computing partition functions.

The assumed units for this EOS are as follows:
* density: 0.253384 g/cc
* pressure: 3.302272e12 erg/cc
* speed: 3.6100785e6 cm/s

### EOS Test Problem Generator

There is a `eos_test.cpp` problem generator with the sole purpose of testing and debugging the desired EOS. There is an associated default input file `inputs/hydro/athinput.eos_test`, the `<problem>` parameters of which are:

    <problem>
    print_table         = true  # Print whole EOS table
    exponentiate_table  = true  # Print 10^<table value>
    eos_loop            = true  # Run EOS user input test loop

To test an ASCII table (such as the above example) one should compile with the command

    > python configure.py --prob eos_test --eos general/eos_table

Here is the output of the EOS test with the the above example table saved as `test.tab`:

    > ./athena -i athinput.eos_test
    Equation of state (EOS) diagnostics:
    General EOS enabled.
    Using table 'test.tab'.
    Shape (nVar, nEgas, nRho): 4, 2, 3
    logEgasMin, logEgasMax: -10, 20
    logRhoMin, logRhoMax: -24, 4
    eUnit, rhoUnit, hUnit: 1, 1, 1
    EosRatios: 1, 10, 10, 0.90909, 
    
    var = 0
    0.1 0.1 0.1 
    0.1 0.1 0.1 
    
    var = 1
    10 10 10 
    10 10 10 
    
    var = 2
    1.1 1.1 1.1 
    1.1 1.1 1.1 
    
    var = 3
    0.1 0.1 0.1 
    0.1 0.1 0.1 
    
    
    Input fluid parameters and retrieve EOS parameters.
    Non-positive inputs will exit loop.
    Input density (mass/volume): 1.0e0
    Input internal energy (energy/volume): 1.0e0
    Density, internal energy: 1, 1
    P(d, e)    , h(d, e)    , ASq(d, P)  , PErr   , ASqErr
    0.1, 1.1, 0.11, 0, 7.2495e-07
    
    Input density (mass/volume): 0
    Input internal energy (energy/volume): 0
    
    
    Setup complete, entering main loop...
    
    cycle=0 time=0 dt=0.0953462
    
    Terminating on cycle limit
    time=0 cycle=0
    tlim=0 nlim=0
    
    zone-cycles = 0
    cpu time used  = 9e-05
    zone-cycles/cpu_second = 0

where we have inputted `1.0e0` via `stdin` for both density and internal energy for the first input loop. `PErr` is the relative error of `SimpleEgas(rho, SimplePres(rho, egas))` and `ASqErr` is the relative difference between `SimpleAsq` and `RiemannAsq`.