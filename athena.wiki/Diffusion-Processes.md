There are several general categories of diffusion processes that Athena++ can model:
* Hydrodynamic diffusion such as fluid viscosity and thermal conduction
* Field diffusion such as Ohmic dissipation and ambipolar diffusion
* [[Passive Scalars]] diffusion of the dye concentration 
 
Each specific diffusion process has been natively implemented in the code.

### Configuration
Users do not need to manually enable the diffusion solver during the configuration stage. These capabilities are automatically invoked when the runtime parameters specify nonzero coefficients. 

However, certain diffusion processes depend on specific physics and will be ignored if the code is configured without the requisite options. For example, in order to calculate Ohmic dissipation and ambipolar diffusion, the magnetic field has to be enabled during configuration with `-b`. Similarly, thermal conduction is calculated only when a non-barotropic equation of state is enabled when [[Configuring]] the code. And passive scalar diffusion will only be enabled if `--nscalars=1` or greater is used when configuring the solver. 

### Input File
The additional microphysics are enrolled automatically through the following optional input parameters: 
```
    <problem>
    nu_iso      = 0.01       # isotropic viscosity coefficient
    nu_aniso    = 0.0        # anisotropic viscosity coefficient
    kappa_iso   = 0.01       # isotropic thermal conduction coefficient
    kappa_aniso = 0.0        # anisotropic thermal conduction coefficient
    eta_ohm     = 0.01       # Ohmic resistivity coefficient
    eta_ad      = 0.0        # Ambipolar diffusion coefficient
    nu_scalar_iso  = 0.01    # isotropic passive scalar diffusion coefficient
```

Positive coefficients invoke the calculation of the corresponding physical processes. Note, the anisotropic viscosity and thermal conduction are **not** currently implemented in Athena++. However, users may define such processes in the [[Problem Generators]] and enroll functions via the `Mesh::EnrollViscosityCoefficient()`, `Mesh::EnrollConductionCoefficient()`, `Mesh::EnrollFieldDiffusivity()` methods. 

### Algorithms
Currently all diffusion processes are calculated from the primitive variables and update the conservative quantities; the resulting diffusion fluxes are directly added to the total hydro and/or field fluxes in an unsplit fashion (no temporal operator splitting is used). In order to calculate the viscosity stress tensor, we make use of the covariant derivative of vectors, with the help of geometric scale factors and their derivatives. Please see [[Covariant Expressions for Vector and Tensor Operators]] for a summary of the general expressions and a table of the scale factors used in the code.

### Units
Note, the `kappa_iso` and `kappa_aniso` coefficients correspond to diffusivities, not conductivities. Also note that the current implementation uses a dimensionless system of units in that the factor (m&#773;/<span>k<sub>B</sub></span>) is not included in calculating the temperature (instead, <var>T</var>=<var>P</var>/<var>d</var> is used). That is, the energy flux is set to - <var>kappa</var>*<var>d</var> d(<var>P</var>/<var>d</var>)/dx, and kappa must have dimensions L<sup>2</sup>/t.  For an energy flux of the form -<var>kappa</var> d<var>T</var>/dx (e.g. with density-independent Spitzer conductivity), the value assigned to <var> kappa</var> must include the factor (<span>k<sub>B</sub></span>/m&#773;) and must include a factor 1/<var> d</var>.
