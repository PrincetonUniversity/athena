### Performance of base Newtonian solver configurations
Typical performance figures for Athena++ are presented in Tables 1 and 2 for serial single-core and MPI-parallelized multicore results on several target architectures. The default second-order temporally-accurate predictor-corrector scheme (`time/integrator=vl2`) was used in all cases, but the reconstruction method of the corrector step was set to either PLM (`time/xorder=2`) or PPM (`time/xorder=3`). Various Riemann solvers were tested.

These figures can be considered as the benchmark values for the basic hydrodynamics and magnetohydrodynamics solver capabilities of the code; optional functionality such as non-Cartesian and/or nonuniform [[Coordinate Systems and Meshes]], [[Special Relativity]], [[General Relativity]], [[Shearing Box]], [[Self Gravity with FFT]], etc. should all be measured relative to the below values. 

<table border="0" class="dataframe">
<caption><b>Table 1</b>: Single-core performance</caption>
  <thead>
    <tr>
      <th></th>
      <th></th>
      <th></th>
      <th colspan="3" halign="left">MZone-cycles/sec</th>
    </tr>
    <tr>
      <th></th>
      <th></th>
      <th></th>
      <th>Xeon Phi<br>KNL 7210</th>
      <th>Broadwell<br>E5-2680 v4</th>
      <th>Skylake-SP<br>Gold 6148</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="6" valign="top">Hydro Sod</th>
      <th rowspan="3" valign="top">PLM</th>
      <th>HLLC</th>
      <td>1.533</td>
      <td>2.730</td>
      <td>4.503</td>
    </tr>
    <tr>
      <th>HLLE</th>
      <td>1.618</td>
      <td>2.868</td>
      <td>4.880</td>
    </tr>
    <tr>
      <th>Roe</th>
      <td>1.555</td>
      <td>2.872</td>
      <td>4.654</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">PPM</th>
      <th>HLLC</th>
      <td>0.752</td>
      <td>1.336</td>
      <td>2.411</td>
    </tr>
    <tr>
      <th>HLLE</th>
      <td>0.762</td>
      <td>1.365</td>
      <td>2.528</td>
    </tr>
    <tr>
      <th>Roe</th>
      <td>0.762</td>
      <td>1.361</td>
      <td>2.424</td>
    </tr>
    <tr>
      <th rowspan="6" valign="top">MHD Brio-Wu</th>
      <th rowspan="3" valign="top">PLM</th>
      <th>HLLD</th>
      <td>0.705</td>
      <td>1.340</td>
      <td>2.403</td>
    </tr>
    <tr>
      <th>HLLE</th>
      <td>0.803</td>
      <td>1.406</td>
      <td>2.307</td>
    </tr>
    <tr>
      <th>Roe</th>
      <td>0.649</td>
      <td>1.143</td>
      <td>1.921</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">PPM</th>
      <th>HLLD</th>
      <td>0.392</td>
      <td>0.719</td>
      <td>1.291</td>
    </tr>
    <tr>
      <th>HLLE</th>
      <td>0.419</td>
      <td>0.749</td>
      <td>1.259</td>
    </tr>
    <tr>
      <th>Roe</th>
      <td>0.373</td>
      <td>0.666</td>
      <td>1.119</td>
    </tr>
  </tbody>
</table><table border="0" class="dataframe">
<caption><b>Table 2</b>: Full node performance</caption>
  <thead>
    <tr>
      <th></th>
      <th></th>
      <th></th>
      <th colspan="3" halign="left">MZone-cycles/sec</th>
    </tr>
    <tr>
      <th></th>
      <th></th>
      <th></th>
      <th>Xeon Phi<br>KNL 7210</th>
      <th>(2x) Broadwell<br>E5-2680 v4</th>
      <th>(2x) Skylake-SP<br>Gold 6148</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="6" valign="top">Hydro Sod</th>
      <th rowspan="3" valign="top">PLM</th>
      <th>HLLC</th>
      <td>66.908</td>
      <td>29.444</td>
      <td>40.750</td>
    </tr>
    <tr>
      <th>HLLE</th>
      <td>67.405</td>
      <td>29.440</td>
      <td>40.764</td>
    </tr>
    <tr>
      <th>Roe</th>
      <td>67.094</td>
      <td>29.418</td>
      <td>40.820</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">PPM</th>
      <th>HLLC</th>
      <td>40.279</td>
      <td>20.269</td>
      <td>32.248</td>
    </tr>
    <tr>
      <th>HLLE</th>
      <td>40.182</td>
      <td>20.328</td>
      <td>32.403</td>
    </tr>
    <tr>
      <th>Roe</th>
      <td>40.196</td>
      <td>20.336</td>
      <td>32.234</td>
    </tr>
    <tr>
      <th rowspan="6" valign="top">MHD Brio-Wu</th>
      <th rowspan="3" valign="top">PLM</th>
      <th>HLLD</th>
      <td>30.886</td>
      <td>16.244</td>
      <td>22.711</td>
    </tr>
    <tr>
      <th>HLLE</th>
      <td>32.526</td>
      <td>16.483</td>
      <td>22.757</td>
    </tr>
    <tr>
      <th>Roe</th>
      <td>29.145</td>
      <td>15.140</td>
      <td>22.673</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">PPM</th>
      <th>HLLD</th>
      <td>19.378</td>
      <td>11.123</td>
      <td>17.733</td>
    </tr>
    <tr>
      <th>HLLE</th>
      <td>20.430</td>
      <td>11.313</td>
      <td>17.684</td>
    </tr>
    <tr>
      <th>Roe</th>
      <td>18.972</td>
      <td>10.623</td>
      <td>17.495</td>
    </tr>
  </tbody>
</table>

<!-- Version 545be021343e1a0a663d6a6f2c76412ac78be7e4 w/ 64^3 MeshBlock. Table generated on 2018-07-26 18:07:29.858757 -->

**Notes on methodology:**
* Both benchmark problems are 3D shock tube problems using the adiabatic equation of state. 
* Each table entry represents the mean of 20 trials of independent, exclusive compute node [Slurm](https://www.schedmd.com/) allocations on clusters managed by [Princeton Research Computing](https://researchcomputing.princeton.edu/)
* The solver was configured with `--nghost=2` for all PLM tests and `--nghost=4` for all PPM tests. 
* The Intel C++ Compiler version 18.0.3 was used to generate all of these results. The only compiler flags used are those defined by the latest version's `--cxx=icc` [Configuring] option.
  * Similarly, the 2018 Revision 3 Intel MPI library was the only MPI library used for the multicore study
* Table 2 uses the same problem size *per-core* as the single-core tests in Table 1. Flat MPI is used to parallelize the problem with 1 rank assigned per physical core. 
  * However, the multicore tests on the KNL were the only set to use hybrid OpenMP+MPI parallelization. Assigning 4 OpenMP threads (each assigned a 64x32x32 MeshBlock) per MPI rank achieved high performance utilizing the 4-way hyperthreading of the 64x physical cores (256 logical cores) on these nodes. See the discussion in [[Using MPI and OpenMP]].

*KNL-specific details:*
- Flat memory mode. Cache memory mode was simulated by prepending the binary call with `numactl -p 1 ./athena ...` to prefer that Athena++ used the ~16 GB of MCDRAM. 
- Quadrant clustering mode.

### Performance cost of optional code features
Under construction.

<!-- TO DO:
- Add additional architecture info to a new table, e.g. # cores
- Attach athinput files
-->