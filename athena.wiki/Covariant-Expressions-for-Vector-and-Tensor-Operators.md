This page summarizes the covariant expressions for vector and tensor operators used in Athena++. These operators are mostly used for computing viscosity and thermal conduction terms in [[Diffusion Processes]]. Please see Appendix A of [Stone & Norman (1992)](http://adsabs.harvard.edu/abs/1992ApJS...80..753S) for details. 

The metric tensor for orthogonal coordinate systems reads

![metric_tensor](images/g_matrix.png)

where <span>h<sub>1</sub></span>, <span>h<sub>2</sub></span>, and <span>h<sub>3</sub></span> are scale factors or Lamé coefficients. They take different forms depending on the specific coordinate system (see Table 1 below).

<div align="center">
<p> <b>Table 1</b>: Scale factors and their derivatives in various coordinate systems. </p>

| scale factors | Cartesian | Cylindrical | Spherical |
| :---:         |    :---:  |   :---:     |  :---:    |
| <span>h<sub>1</sub></span>   | 1     | 1   | 1 |
| <span>h<sub>2</sub></span>   | 1     | <span>r </span>   | <span>r </span> |
| <span>h<sub>3</sub> = h<sub>31</sub> &middot;  h<sub>32</sub></span>   | 1     | <span>1 </span>   | <span>r &middot; sin&theta; </span>  |
| <span>h<sub>31</sub></span>   | 1     | <span>1 </span>   | <span>r </span>  |
| <span>h<sub>32</sub></span>   | 1     | <span>1 </span>   | <span>sin&theta; </span>  |
| <span>&#8706;h<sub>2</sub> / &#8706;x<sub>1</sub></span> |  0 |  1  | 1  |
| <span>&#8706;h<sub>31</sub> / &#8706;x<sub>1</sub></span> |  0 |  0  | 1  |
| <span>&#8706;h<sub>32</sub> / &#8706;x<sub>2</sub></span> |  0 |  0  | <span>cos&theta;</span>  |

</div>

### Gradient of a scalar:
<span>&nabla; &Phi; = &#8706;&Phi; / (h<sub>i</sub> &#8706;x<sub>i</sub>) <b>&ecirc;<sub>i</sub></b> </span>,

with <span>i &isin; [x1,x2,x3]</span> and assuming Einstein summation (the repeated indices are implicitly summed over). 

### Divergence of a vector:
<span>&nabla; &middot; <b>F</b> = g<sup>-&frac12;</sup> &#8706;(g<sup>&frac12;</sup> F<sub>i</sub> / h<sub>i</sub>) / &#8706;x<sub>i</sub></span>,

where <span><b>F</b> = F<sub>i</sub> <b>&ecirc;<sub>i</sub></b> </span>, and 
<span>g<sup>&frac12;</sup> = h<sub>1</sub>h<sub>2</sub>h<sub>3</sub></span>.

### Curl of a vector:
<span>&nabla; &times; <b>F</b> = h<sub>k</sub> g<sup>-&frac12;</sup>; &varepsilon;<sub>kji</sub> &#8706;(h<sub>i</sub> F<sub>i</sub>) / &#8706;x<sub>j</sub> <b>&ecirc;<sub>i</sub></b></span>.

### Gradient of a vector (covariant derivative):
<span>&nabla; <b>F</b> = [&#8706;(h<sub>i</sub>F<sub>i</sub>) / &#8706;x<sub>j</sub> - &Gamma;<sup>k</sup><sub>i j</sub> h<sub>k</sub> F<sub>k</sub>] (h<sub>i</sub> h<sub>j</sub>)<sup>-1</sup>  <b>&ecirc;<sub>i</sub></b> <b>&ecirc;<sub>j</sub></b> </span>,

where <span>&Gamma;<sup>k</sup><sub>i j</sub></span> is the Christoffel symbol:

<span>&Gamma;<sup>i</sup><sub>i i</sub> = (2h<sub>i</sub><sup>2</sup>)<sup>-1</sup> &#8706; h<sub>i</sub><sup>2</sup> / &#8706; x<sub>i</sub></span>

<span>&Gamma;<sup>i</sup><sub>i j</sub> = &Gamma;<sup>i</sup><sub>j i</sub> = (2h<sub>i</sub><sup>2</sup>)<sup>-1</sup> &#8706; h<sub>i</sub><sup>2</sup> / &#8706; x<sub>j</sub></span>

<span>&Gamma;<sup>i</sup><sub>j j</sub> = -(2h<sub>i</sub><sup>2</sup>)<sup>-1</sup> &#8706; h<sub>j</sub><sup>2</sup> / &#8706; x<sub>i</sub></span>

<span>&Gamma;<sup>i</sup><sub>j k</sub> = 0</span>

### Example
For instance, the viscous stress tensor (neglecting the bulk viscosity) can be written as

<span>&sigma;</sup><sub>i j</sub> = -&mu;( v<sub>i; j</sub> + v<sub>j; i</sub> -&frac23; v<sup>k</sup><sub>;k</sub> g<sub>i j</sub> )<span>,

where <span>&mu; = &rho; &nu;</span> is the dynamic viscosity. The expression consists of two gradient terms and one term containing the divergence of the velocity. All of these terms can be calculated in any of the supported [[Coordinate Systems and Meshes]] using the above formulae.
