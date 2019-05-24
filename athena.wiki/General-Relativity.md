### Governing equations

Let <var>g</var><sub><var>&mu;</var><var>&nu;</var></sub> be the components of a stationary metric with time coordinate <var>x</var><sup>&thinsp;0</sup>, so &#8706;<sub>0</sub><var>g</var><sub><var>&mu;</var><var>&nu;</var></sub>&#8287;=&#8287;0. We can define the lapse <var>&alpha;</var>&#8287;=&#8287;(-<var>g</var><sup>&thinsp;00</sup>)<sup>-1/2</sup>. Then the unit future-pointing timelike normal to surfaces of constant <var>x</var><sup>&thinsp;0</sup> has components <var>n</var><sub><var>&mu;</var></sub>&#8287;=&#8287;-<var>&alpha;</var><var>&delta;</var><sub>0<var>&mu;</var></sub>. Projection onto the spatial part of the tangent space can be accomplished with the components <var>j</var><sub><var>&mu;</var><var>&nu;</var></sub>&#8287;=&#8287;<var>g</var><sub><var>&mu;</var><var>&nu;</var></sub>&#8287;+&#8287;<var>n</var><sub><var>&mu;</var></sub><var>n</var><sub><var>&nu;</var></sub>.

When working with general relativity, the primitive variables in Athena++ are taken to be

* comoving rest mass density <var>&rho;</var>,
* gas pressure <var>p</var><sub>gas</sub>, and
* spatial components of the projected 4-velocity <var>u&#771;</var><sup>&thinsp;<var>i</var></sup>&#8287;=&#8287;<var>j</var><sup><var>i</var></sup><sub><var>&mu;</var></sub>u<sup>&thinsp;<var>&mu;</var></sup>.

In MHD there will also be the magnetic field components <var>B</var><sup>&thinsp;<var>i</var></sup>&#8287;=&#8287;(&#9733;<var>F</var>)<sup><var>i</var>0</sup>, where <var>F</var> is the electromagnetic field tensor. Note that these are not the projected variables &#8492;<sup><var>i</var></sup>&#8287;=&#8287;-<var>j</var><sup>&thinsp;<var>i</var></sup><sub><var>&mu;</var></sub><var>n</var><sub><var>&nu;</var></sub>(&#9733;<var>F</var>)<sup><var>&mu;</var><var>&nu;</var></sup>&#8287;=&#8287;<var>&alpha;</var><var>B</var><sup>&thinsp;<var>i</var></sup>.

The following quantities are useful to define:

* projected field components <var>b</var><sup>&thinsp;<var>&mu;</var></sup>&#8287;=&#8287;<var>u</var><sub>&nu;</sub>(&#9733;<var>F</var>)<sup><var>&nu;</var><var>&mu;</var></sup>, in particular
  * <var>b</var><sup>&thinsp;0</sup>&#8287;=&#8287;<var>g</var><sub><var>i</var><var>&mu;</var></sub><var>B</var><sup>&thinsp;<var>i</var></sup><var>u</var><sup>&thinsp;<var>&mu;</var></sup>,
  * <var>b</var><sup>&thinsp;<var>i</var></sup>&#8287;=&#8287;(<var>B</var><sup>&thinsp;<var>i</var></sup>&#8287;+&#8287;<var>b</var><sup>&thinsp;0</sup><var>u</var><sup>&thinsp;<var>i</var></sup>)/<var>u</var><sup>&thinsp;0</sup>;
* magnetic pressure <var>p</var><sub>mag</sub>&#8287;=&#8287;<var>b</var><sub><var>&mu;</var></sub><var>b</var><sup>&thinsp;<var>&mu;</var></sup>/2;
* total pressure <var>p</var><sub>tot</sub>&#8287;=&#8287;<var>p</var><sub>gas</sub>&#8287;+&#8287;<var>p</var><sub>mag</sub>;
* adiabatic index &Gamma;;
* gas enthalpy <var>w</var><sub>gas</sub>&#8287;=&#8287;<var>&rho;</var>&#8287;+&#8287;&Gamma;/(&Gamma;&#8722;1)&#8287;&times;&#8287;<var>p</var><sub>gas</sub>;
* magnetic enthalpy <var>w</var><sub>mag</sub>&#8287;=&#8287;2<var>p</var><sub>mag</sub>;
* total enthalpy <var>w</var><sub>tot</sub>&#8287;=&#8287;<var>w</var><sub>gas</sub>&#8287;+&#8287;<var>w</var><sub>mag</sub>;
* stress-energy components <var>T</var><sup>&thinsp;<var>&mu;</var><var>&nu;</var></sup>&#8287;=&#8287;<var>w</var><sub>tot</sub><var>u</var><sup>&thinsp;<var>&mu;</var></sup><var>u</var><sup>&thinsp;<var>&nu;</var></sup>&#8287;+&#8287;<var>p</var><sub>tot</sub><var>g</var><sup>&thinsp;<var>&mu;</var><var>&nu;</var></sup>&#8287;&#8722;&#8287;<var>b</var><sup>&thinsp;<var>&mu;</var></sup><var>b</var><sup>&thinsp;<var>&nu;</var></sup>;
* metric determinant <var>g</var>;
* connection coefficients &Gamma;<sup><var>&sigma;</var></sup><sub><var>&mu;</var><var>&nu;</var></sub>&#8287;=&#8287;1/2&#8287;&times;&#8287;<var>g</var><sup>&thinsp;<var>&sigma;</var><var>&lambda;</var></sup>(&#8706;<sub><var>&mu;</var></sub><var>g</var><sub><var>&nu;</var><var>&lambda;</var></sub>&#8287;+&#8287;&#8706;<sub><var>&nu;</var></sub><var>g</var><sub><var>&mu;</var><var>&lambda;</var></sub>&#8287;&#8722;&#8287;&#8706;<sub><var>&lambda;</var></sub><var>g</var><sub><var>&mu;</var><var>&nu;</var></sub>).

The evolution equations are

* &#8706;<sub>0</sub>&thinsp;((-<var>g</var>)<sup>1/2</sup><var>&rho;</var><var>u</var><sup>&thinsp;0</sup>)&#8287;+&#8287;&#8706;<sub><var>j</var></sub>&thinsp;((-<var>g</var>)<sup>1/2</sup><var>&rho;</var><var>u</var><sup>&thinsp;<var>j</var></sup>)&#8287;=&#8287;0,
* &#8706;<sub>0</sub>&thinsp;((-<var>g</var>)<sup>1/2</sup><var>T</var><sup>&thinsp;0</sup><sub><var>&mu;</var></sub>)&#8287;+&#8287;&#8706;<sub><var>j</var></sub>((-<var>g</var>)<sup>1/2</sup><var>T</var><sup>&thinsp;<var>j</var></sup><sub><var>&mu;</var></sub>)&#8287;=&#8287;(-<var>g</var>)<sup>1/2</sup><var>T</var><sup>&thinsp;<var>&nu;</var></sup><sub><var>&sigma;</var></sub>&Gamma;<sup><var>&sigma;</var></sup><sub><var>&mu;</var><var>&nu;</var></sub>,
* &#8706;<sub>0</sub>&thinsp;((-<var>g</var>)<sup>1/2</sup>(&#9733;<var>F</var>)<sup><var>i</var>0</sup>)&#8287;+&#8287;&#8706;<sub><var>j</var></sub>&thinsp;((-<var>g</var>)<sup>1/2</sup>(&#9733;<var>F</var>)<sup><var>i</var><var>j</var></sup>)&#8287;=&#8287;0.

Athena++ considers the conserved variables to be

* conserved density <var>D</var>&#8287;=&#8287;<var>&rho;</var><var>u</var><sup>&thinsp;0</sup>,
* energy density <var>E</var>&#8287;=&#8287;<var>T</var><sup>&thinsp;0</sup><sub>0</sub>, and
* momentum density <var>M</var><sub><var>i</var></sub>&#8287;=&#8287;<var>T</var><sup>&thinsp;0</sup><sub><var>i</var></sub>,

as well as the magnetic fields <var>B</var><sup>&thinsp;<var>i</var></sup>. In particular, the factors of (-<var>g</var>)<sup>1/2</sup> are not included in the output data, nor should they be included when setting conserved variables in problem generators. Also, these are coordinate-frame quantities, not quantities in the normal frame often seen in 3+1 decompositions. This is the difference between components <var>A</var><sup>&thinsp;<var>&mu;</var></sup> (coordinate frame) and -<var>n</var><sub><var>&mu;</var></sub><var>A</var><sup>&thinsp;<var>&mu;</var></sup> and <var>j</var><sup>&thinsp;<var>i</var></sup><sub><var>&mu;</var></sub><var>A</var><sup>&thinsp;<var>&mu;</var></sup> (normal frame).

### Considerations for GR in Athena++

To enable GR, add the `-g` flag when configuring. The `-s` flag should **not** be used, as this is reserved for SR without GR. There are several important choices one must make, detailed here.

##### Coordinates (`--coord` option)

Coordinate systems for GR are entirely distinct from those used in SR or Newtonian simulations. In particular, `cartesian`, `cylindrical`, and `spherical_polar` coordinates do not function in GR. Instead choices include:

* `minkowski` (<var>t</var>,<var>x</var>,<var>y</var>,<var>z</var>): <var>ds</var><sup>&thinsp;2</sup>&#8287;=&#8287;-d<var>t</var><sup>&thinsp;2</sup>&#8287;+&#8287;d<var>x</var><sup>&thinsp;2</sup>&#8287;+&#8287;d<var>y</var><sup>&thinsp;2</sup>&#8287;+&#8287;d<var>z</var><sup>&thinsp;2</sup>;
* `schwarzschild` (<var>t</var>,<var>r</var>,<var>&theta;</var>,<var>&#981;</var>): <var>ds</var><sup>&thinsp;2</sup>&#8287;=&#8287;-(1&#8722;2<var>M</var>/<var>r</var>)d<var>t</var><sup>&thinsp;2</sup>&#8287;+&#8287;(1&#8722;2<var>M</var>/<var>r</var>)<sup>-1</sup>d<var>r</var><sup>&thinsp;2</sup>&#8287;+&#8287;<var>r</var><sup>&thinsp;2</sup>d<var>&theta;</var><sup>&thinsp;2</sup>&#8287;+&#8287;<var>r</var><sup>&thinsp;2</sup>sin<sup>2</sup>(<var>&theta;</var>)d<var>&#981;</var><sup>&thinsp;2</sup>;
* `kerr-schild` (<var>t</var>,<var>r</var>,<var>&theta;</var>,<var>&#981;</var>): <var>ds</var><sup>&thinsp;2</sup>&#8287;=&#8287;<var>g</var><sub><var>&mu;</var><var>&nu;</var></sub>d<var>x</var><sup>&thinsp;<var>&mu;</var></sup>d<var>x</var><sup>&thinsp;<var>&nu;</var></sup>, with:
  * &Sigma;&#8287;=&#8287;<var>r</var><sup>&thinsp;2</sup>&#8287;+&#8287;<var>a</var><sup>&thinsp;2</sup>cos<sup>2</sup>(<var>&theta;</var>),
  * <var>g</var><sub>00</sub>&#8287;=&#8287;-(1&#8722;2<var>M</var><var>r</var>/&Sigma;),
  * <var>g</var><sub>01</sub>&#8287;=&#8287;2<var>M</var><var>r</var>/&Sigma;,
  * <var>g</var><sub>03</sub>&#8287;=&#8287;-(2<var>M</var><var>a</var><var>r</var>/&Sigma;)sin<sup>2</sup>(<var>&theta;</var>),
  * <var>g</var><sub>11</sub>&#8287;=&#8287;1&#8287;+&#8287;2<var>M</var><var>r</var>/&Sigma;,
  * <var>g</var><sub>13</sub>&#8287;=&#8287;-(1&#8287;+&#8287;2<var>M</var><var>r</var>/&Sigma;)<var>a</var>&thinsp;sin<sup>2</sup>(<var>&theta;</var>),
  * <var>g</var><sub>22</sub>&#8287;=&#8287;&Sigma;,
  * <var>g</var><sub>33</sub>&#8287;=&#8287;(<var>r</var><sup>&thinsp;2</sup>&#8287;+&#8287;<var>a</var><sup>&thinsp;2</sup>&#8287;+&#8287;(2<var>M</var><var>a</var><sup>&thinsp;2</sup><var>r</var>/&Sigma;)sin<sup>2</sup>(<var>&theta;</var>))sin<sup>2</sup>(<var>&theta;</var>),
  * <var>g</var><sub>02</sub>&#8287;=&#8287;<var>g</var><sub>12</sub>&#8287;=&#8287;<var>g</var><sub>23</sub>&#8287;=&#8287;0;
* `gr_user`: [[user-defined coordinates|Arbitrary Coordinates]].

Kerr-Schild coordinates are a standard choice for black holes, since they are horizon penetrating.

Note that all primitive and conserved variables are in the coordinate basis. The natural basis vectors are not normalized, even in cases where it would be easy to do so such as Schwarzschild.

##### Equation of state (`--eos` option)

GR in Athena++ only works with an adiabatic (gamma-law) equation of state.

##### Riemann solver (`--flux` option)

In order of decreasing numerical diffusion, the Riemann solvers compatible with GR are `llf` (local Lax-Friedrichs), `hlle`, `hllc`, and `hlld`. HLLC only works for hydrodynamics, and HLLD only works for MHD. In addition, HLLC and HLLD can only work using frame transformations at interfaces to use SR-compatible algorithms (see the [Athena++ GRMHD method paper][GRMHD]), and thus they require the `-t` option. Frame transformations can be either on or off when using LLF or HLLE.

At the polar axis coordinate singularity in Schwarzschild and Kerr-Schild coordinates, frame transformations become singular. The Riemann solver will use the appropriate non-transforming method (LLF in that case, HLLE in all other cases) at those interfaces automatically.

##### Mesh refinement

Athena++'s mesh refinement, both [[static|Static-Mesh-Refinement]] and [[adaptive|Adaptive-Mesh-Refinement]], is designed to be fully compatible with GR. Note however that GR simulations often involve spherical coordinate systems, and [[polar boundaries require careful consideration|Boundary-Conditions#polar-boundary-conditions]].

##### Timestep size

The size of a timestep is set by the shortest light crossing time in the grid. It does not change even if the sound (hydrodynamics) or fast magnetosonic (MHD) speeds are considerably slower than unity (i.e. the speed of light). The <var>x</var><sup>&thinsp;1</sup>-width of cells is defined to be &int;(<var>g</var><sub>11</sub>)<sup>1/2</sup>d<var>x</var><sup>&thinsp;1</sup> along the line of constant <var>x</var><sup>&thinsp;2</sup> and <var>x</var><sup>&thinsp;3</sup> from one face to the other passing through the cell center. The <var>x</var><sup>&thinsp;2</sup>- and <var>x</var><sup>&thinsp;3</sup>-widths are defined likewise. In coordinate systems where the integral has no simple antiderivative, lower bounds are used instead. Also, in complex coordinate systems the cell center may be defined as the arithmetic mean of the face coordinates rather than the halfway point in volume.

##### Floors and ceilings on variables

As with SR, the complexity of the equations of GR can prove troublesome for robustness. There are a number of ways in which a simulation can produce problematic values, often manifesting in recovering negative densities, negative pressures, or superluminal velocities in the conserved-to-primitive variable inversion.

In order to help robustness, various limits are available in the `<hydro>` section of the input file, described below:

<table>
  <tr> <th>Name</th> <th>Default Value</th> <th>Description</th> </tr>
  <tr> <td><code>dfloor</code></td> <td>&asymp;10<sup>-35</sup></td> <td rowspan="3"><var>&rho;</var><sub>floor</sub>&#8287;=&#8287;max(<code>dfloor</code>,&#8287;<code>rho_min</code>&#8287;&times;&#8287;(<var>x</var><sup>&thinsp;1</sup>)<sup><code>rho_pow</code></sup>)</td> </tr>
  <tr> <td><code>rho_min</code></td> <td><code>dfloor</code></td> </tr>
  <tr> <td><code>rho_pow</code></td> <td>0</td> </tr>
  <tr> <td><code>pfloor</code></td> <td>&asymp;10<sup>-35</sup></td> <td rowspan="3"><var>p</var><sub>gas,floor</sub>&#8287;=&#8287;max(<code>pfloor</code>,&#8287;<code>pgas_min</code>&#8287;&times;&#8287;(<var>x</var><sup>&thinsp;1</sup>)<sup><code>pgas_pow</code></sup>)</td> </tr>
  <tr> <td><code>pgas_min</code></td> <td><code>pfloor</code></td> </tr>
  <tr> <td><code>pgas_pow</code></td> <td>0</td> </tr>
  <tr> <td><code>gamma_max</code></td> <td>1000</td> <td>ceiling on <var>&gamma;</var></td> </tr>
  <tr> <td><code>sigma_max</code></td> <td>0</td> <td>ceiling on 2<var>p</var><sub>mag</sub>/<var>&rho;</var>; only used if positive</td> </tr>
  <tr> <td><code>beta_min</code></td> <td>0</td> <td>floor on <var>p</var><sub>gas</sub>/<var>p</var><sub>mag</sub>; only used if positive</td> </tr>
</table>

The power-law scalings of <var>&rho;</var><sub>floor</sub> and <var>p</var><sub>gas,floor</sub> are meant to be used in spherical-like coordinates, where <var>x</var><sup>&thinsp;1</sup> is the radial coordinate. These scaling components of the floors are ignored if the exponents are 0.

The Lorentz factor is in the normal frame: <var>&gamma;</var>&#8287;=&#8287;-<var>n</var><sub><var>&mu;</var></sub><var>u</var><sup>&thinsp;<var>&mu;</var></sup>&#8287;=&#8287;<var>&alpha;</var><var>u</var><sup>&thinsp;0</sup>&#8287;=&#8287;(1&#8287;+&#8287;<var>g</var><sub><var>i</var><var>j</var></sub><var>u&#771;</var><sup>&thinsp;<var>i</var></sup><var>u&#771;</var><sup>&thinsp;<var>j</var></sup>)<sup>1/2</sup>. If a value <var>&gamma;</var>&#8287;&gt;&#8287;<var>&gamma;</var><sub>max</sub> is found, the velocity components <var>u&#771;</var><sup>&thinsp;<var>i</var></sup> are each multiplied by ((<var>&gamma;</var><sub>max</sub><sup>2</sup>&#8722;1)/(<var>&gamma;</var><sup>&thinsp;2</sup>&#8722;1))<sup>1/2</sup> to make <var>&gamma;</var>&#8287;=&#8287;<var>&gamma;</var><sub>max</sub>.

For the limits on plasma <var>&sigma;</var> and <var>&beta;</var>, <var>&rho;</var> and <var>p</var><sub>gas</sub> are adjusted and the magnetic field is left alone in order to preserve &#8706;<sub><var>i</var></sub>&thinsp;((-<var>g</var>)<sup>1/2</sup><var>B</var><sup>&thinsp;<var>i</var></sup>)&#8287;=&#8287;0.

  [GRMHD]: http://adsabs.harvard.edu/abs/2016ApJS..225...22W