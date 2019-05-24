The Python script [`athena/vis/python/uniform.py`][1] can be used to rewrite `.athdf` HDF5 files with a single, uniform MeshBlock. This is useful for reading into [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/) (whose performance suffers when thousands of MeshBlocks are present), or for any other analysis tool that would be better served without dealing with the MeshBlock structure. 

*Note*: It is a known issue in VisIt that the number of `MeshBlocks` (stored as HDF5 data containers) in an open database is not updated when the timeslice is changed. This may result in missing/empty patches of the domain when AMR derefines a `MeshBlock`; reopening the database will fix the missing data.

The script is used as
```
python uniform.py <input_filename> <output_filename> <start> <end> <stride> [options]
```

The necessary arguments are:
- `input_filename`: base name of files to be rewritten, including the path but excluding the 5-digit number. Files are not modified unless `output_filename` is the same.
- `output_filename`: base name of output files to be written. This can be same as `input_filename`.
- `start`: integer of first file to rewrite, corresponding to 5-digit identifier (does not need to include leading zeros).
- `end`: integer of last file to rewrite, corresponding to 5-digit identifier (does not need to include leading zeros).
- `stride`: stride in file identifiers, useful for only rewriting some of them.

Optional arguments are:
- `-m`: run script in parallel in an MPI environment. Files are distributed equally among the processes.
- `-x`: do not write corresponding XDMF (`.athdf.xdmf`) files
- `-l <level>`, `--level <level>`: force refinement level of output to be this, relative to root grid in input. Must be nonnegative integer. If omitted, each file is upsampled to its own maximum refinement level.
- `-q <quantity_1> [<quantity_2> ...]`, `--quantities <quantity_1> [<quantity_2> ...]`: list of quantities to include in rewriting. These should be variable names, not dataset names. If omitted, all cell-centered quantities are included.

For example, to create 10 new files that include only density and the first component of magnetic field, one might use the script as
```
python uniform.py data/old.prim data/new.prim 0 9 1 --quantities rho Bcc1
```
This would use the data in `data/old.prim.0000{0..9}.athdf` to create `data/new.prim.0000{0..9}.athdf` and `data/new.prim.0000{0..9}.athdf.xdmf`.

If a refinement level is forced with `--level` to be less than the maximum refinement of a file, the `subsample=True` option will be used in the internal call to the [Python HDF5 data reader][2].

The new HDF5 file will internally consist of a single MeshBlock covering the entire domain. The format is the same as the [[HDF5 format produced by Athena++|HDF5-Format]], with the exception that datasets are not preserved if the `--quantities` option is used. In this case the only cell-centered dataset present will be called `"quantities"`, and all variables will be inside this dataset.

  [1]: https://github.com/PrincetonUniversity/athena-public-version/blob/master/vis/python/uniform.py
  [2]: Reading-data-into-Python#athdffilename-datanone-quantitiesnone-level0-subsamplefalse-fast_restrictfalse-vol_funcnone-vol_paramsnone