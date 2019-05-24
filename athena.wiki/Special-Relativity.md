### Governing equations

Athena++ supports special-relativistic fluid dynamics. In this mode the primitive variables are

* comoving rest mass density <var>&rho;</var>,
* gas pressure <var>p</var><sub>gas</sub>, and
* 3-velocity **v** with components <var>v</var><sup>&thinsp;<var>i</var></sup>.

In MHD we also have a magnetic field **B** with components <var>B</var><sup>&thinsp;<var>i</var></sup>.

The conserved variables are

* coordinate-frame density <var>D</var>&#8287;=&#8287;<var>&gamma;</var><var>&rho;</var>,
* total energy density <var>E</var>&#8287;=&#8287;<var>&gamma;</var><sup>&thinsp;2</sup><var>w</var><sub>tot</sub>&#8287;&#8722;&#8287;<var>p</var><sub>tot</sub>&#8287;&#8722;&#8287;<var>b</var><sup>&thinsp;<var>t</var></sup><var>b</var><sup>&thinsp;<var>t</var></sup>, and
* momentum density **M** with components <var>M</var><sup>&thinsp;<var>i</var></sup>&#8287;=&#8287;<var>&gamma;</var><sup>&thinsp;2</sup><var>w</var><sub>tot</sub><var>v</var><sup>&thinsp;<var>i</var></sup>&#8287;&#8722;&#8287;<var>b</var><sup>&thinsp;<var>t</var></sup><var>b</var><sup>&thinsp;<var>i</var></sup>.

Here we define the following useful quantities:

* Lorentz factor <var>&gamma;</var>&#8287;=&#8287;(1&#8722;**v**<sup>2</sup>)<sup>-1/2</sup>,
* projected 4-field with components <var>b</var><sup>&thinsp;<var>t</var></sup>&#8287;=&#8287;<var>&gamma;&thinsp;</var>**v&#8901;B** and <var>b</var><sup>&thinsp;<var>i</var></sup>&#8287;=&#8287;<var>B</var><sup>&thinsp;<var>i</var></sup>/<var>&gamma;</var>&#8287;+&#8287;<var>b</var><sup>&thinsp;<var>t</var></sup><var>v</var><sup>&thinsp;<var>i</var></sup>,
* magnetic pressure <var>p</var><sub>mag</sub>&#8287;=&#8287;<var>&eta;</var><sub><var>&mu;</var><var>&nu;</var></sub><var>b</var><sup>&thinsp;<var>&mu;</var></sup><var>b</var><sup>&thinsp;<var>&nu;</var></sup>/2,
* total pressure <var>p</var><sub>tot</sub>&#8287;=&#8287;<var>p</var><sub>gas</sub>&#8287;+&#8287;<var>p</var><sub>mag</sub>,
* adiabatic index &Gamma;,
* comoving gas enthalpy <var>w</var><sub>gas</sub>&#8287;=&#8287;<var>&rho;</var>&#8287;+&#8287;&Gamma;/(&Gamma;&#8722;1)&#8287;&times;&#8287;<var>p</var><sub>gas</sub>,
* comoving magnetic enthalpy <var>w</var><sub>mag</sub>&#8287;=&#8287;2<var>p</var><sub>mag</sub>,
* comoving total enthalpy <var>w</var><sub>tot</sub>&#8287;=&#8287;<var>w</var><sub>gas</sub>&#8287;+&#8287;<var>w</var><sub>mag</sub>.

The evolution equations are then

* &#8706;<sub><var>t</var></sub><var>D</var>&#8287;+&#8287;&#8706;<sub><var>j</var></sub>&thinsp;(<var>&gamma;</var><var>&rho;</var><var>v</var><sup>&thinsp;<var>j</var></sup>)&#8287;=&#8287;0,
* &#8706;<sub><var>t</var></sub><var>E</var>&#8287;+&#8287;&#8706;<sub><var>j</var></sub><var>M</var><sup>&thinsp;<var>j</var></sup>&#8287;=&#8287;0,
* &#8706;<sub><var>t</var></sub><var>M</var><sup>&thinsp;<var>i</var></sup>&#8287;+&#8287;&#8706;<sub><var>j</var></sub>&thinsp;(<var>&gamma;</var><sup>&thinsp;2</sup><var>w</var><sub>tot</sub><var>v</var><sup>&thinsp;<var>i</var></sup><var>v</var><sup>&thinsp;<var>j</var></sup>&#8287;+&#8287;<var>p</var><sub>tot</sub><var>&delta;</var><sup>&thinsp;<var>i</var><var>j</var></sup>&#8287;&#8722;&#8287;<var>b</var><sup>&thinsp;<var>i</var></sup><var>b</var><sup>&thinsp;<var>j</var></sup>)&#8287;=&#8287;0, and
* &#8706;<sub><var>t</var></sub><var>B</var><sup>&thinsp;<var>i</var></sup>&#8287;+&#8287;&#8706;<sub><var>j</var></sub>&thinsp;(<var>&gamma;</var><var>v</var><sup>&thinsp;<var>j</var></sup><var>b</var><sup>&thinsp;<var>i</var></sup>&#8287;&#8722;&#8287;<var>&gamma;</var><var>v</var><sup>&thinsp;<var>i</var></sup><var>b</var><sup>&thinsp;<var>j</var></sup>)&#8287;=&#8287;0.

The momentum equation might have source terms on the right-hand side in non-Cartesian coordinates.

### General code considerations

To enable special relativity, add the `-s` flag when configuring. This option will automatically select an appropriate Riemann solver and equation of state. The Riemann solver can be local Lax-Friedrichs, HLLE, HLLC (hydrodynamics only), or HLLD (MHD only). Note the Roe solver is unavailable for relativity. The equation of state must be adiabatic (i.e. gamma-law).

SR is compatible with the same set of coordinates as Newtonian physics: Cartesian, cylindrical, and spherical-polar. Note the coordinates labeled "Minkowski" are for use with _general_ relativity only. For those with GR experience, it is worth noting that cylindrical and spherical-polar coordinates are implemented in a locally orthonormal way for SR, rather than in the natural coordinate way. For example, a purely radial velocity field that is everywhere half the speed of light will have constant components (0.5, 0, 0) in spherical coordinates, rather than spatially varying components and understood basis vectors that also vary in magnitude spatially. This is reflected in the output data.

When calculating timesteps in SR, Athena++ assumes a signal speed of the speed of light (i.e. unity). Therefore a fixed grid will have fixed timesteps given by the CFL number multiplied by the smallest cell width.

SR is fully compatible with static and adaptive mesh refinement.

### Floors and ceilings on variables

Caution is advised when simulating extreme conditions (very high Lorentz factors, high magnetizations, etc.). Due to the extra nonlinear coupling in SR compared to Newtonian physics, there are more ways in which a simulation can go wrong, often manifesting in recovering superluminal velocities in the conserved-to-primitive variable inversion.

In order to help robustness, various limits are available in the `<hydro>` section of the input file, described below:

<table>
  <tr> <th>Name</th> <th>Default Value</th> <th>Description</th> </tr>
  <tr> <td><code>dfloor</code></td> <td>&asymp;10<sup>-35</sup></td> <td>floor on <var>&rho;</var></td> </tr>
  <tr> <td><code>pfloor</code></td> <td>&asymp;10<sup>-35</sup></td> <td>floor on <var>p</var><sub>gas</sub></td> </tr>
  <tr> <td><code>gamma_max</code></td> <td>1000</td> <td>ceiling on <var>&gamma;</var></td> </tr>
  <tr> <td><code>sigma_max</code></td> <td>0</td> <td>ceiling on 2<var>p</var><sub>mag</sub>/<var>&rho;</var>; only used if positive</td> </tr>
  <tr> <td><code>beta_min</code></td> <td>0</td> <td>floor on <var>p</var><sub>gas</sub>/<var>p</var><sub>mag</sub>; only used if positive</td> </tr>
</table>

If Lorentz factors above <var>&gamma;</var><sub>max</sub> are found, the 3-velocity components are all rescaled by the same factor in order to make <var>&gamma;</var>&#8287;=&#8287;<var>&gamma;</var><sub>max</sub>. For the limits on plasma <var>&sigma;</var> and <var>&beta;</var>, <var>&rho;</var> and <var>p</var><sub>gas</sub> are adjusted and the magnetic field is left alone in order to preserve &nabla;&#8901;**B**&#8287;=&#8287;0.