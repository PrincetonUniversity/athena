### Visualization Software
Athena++ can output formatted table, VTK, and HDF5 data (see also [[Outputs]]).

The formatted table (`.tab`) format is intended for relatively small simulations, in particular 1D and 2D ones. This is just a simple text file which can be plotted by your favorite plotting software. For example, we recommend [[gnuplot|http://www.gnuplot.info/]] as it is simple and free. The meaning of each column is described in the file header.

VTK and HDF5 are standard formats commonly used in numerical simulations. Because VTK is non-parallel IO, it is intended for relatively small 2D and 3D simulations. You can merge the VTK files using `vis/vtk/join_vtk++.c`, but currently this code works only for simulations without mesh refinement. For massively parallel simulations and/or simulations with mesh refinement, HDF5 is strongly recommended.

To visualize VTK or HDF5 files, there are some options publicly available, such as [[VisIt|https://visit.llnl.gov/]] and [[ParaView|http://www.paraview.org/]]. [[yt|http://yt-project.org/]] also supports the Athena++ HDF5 output with some minor limitations. For usage of these programs, please consult their documentation.

To read HDF5 data with VisIt or ParaView, load the corresponding `.xdmf` (eXtensible Data Model and Format) files, which tell the software how to read the `.athdf` files. All the variables are stored as scalars. For convenience, we provide typical variable expression files for VisIt in the `vis/vtk` directory. You can load these files from "Controls → Expressions".

For very large simulations (with thousands of MeshBlocks), VisIt can be slow to read in data. It may be helpful to first [[resample|Resampling-HDF5-Outputs]] the HDF5 files in this case.

### Python Plotting and Analysis
Several simple plotting scripts are included. These are described in [[Plotting Scripts]]. For more detailed analysis, see [[Reading Data into Python]] for how to work with Athena++ outputs in a Python script.

### Reading .athdf Files
The format of `.athdf` files is HDF5, and can be read from your favorite programming language using the HDF5 library/package/etc. The data structure is pretty self-explanatory, and you can see it using [[HDFVIEW|https://www.hdfgroup.org/products/java/hdfview/]] provided by the HDF Group. See also [[HDF5-Format]].
