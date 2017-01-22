#! /usr/bin/env python

"""
Script for plotting vertical (r,theta) or midplane (r,phi) slices of data in
spherical coordinates.

Run "python azimuthal_average.py -h" to see description of inputs.

See documentation on athena_read.athdf() for important notes about reading files
with mesh refinement.

The -c <colormap> option is highly recommended to change the default. Consider
"RdBu_r" or "gist_heat" for example.

Users are encouraged to make their own versions of this script for improved
results by adjusting figure size, spacings, tick locations, axes labels, etc.
The script must also be modified to plot any functions of the quantities in the
file, including combinations of multiple quantities.
"""

# Python modules
import numpy as np
import argparse
import h5py

# Python plotting modules
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# Athena++ modules
import athena_read

# Main function
def main(**kwargs):

  # Determine refinement level to use
  if kwargs['level'] is not None:
    level = kwargs['level']
  else:
    with h5py.File(kwargs['data_file'], 'r') as f:
      level = f.attrs['MaxLevel']

  # Read data
  data = athena_read.athdf(kwargs['data_file'], quantities=(kwargs['quantity'],), \
      level=level)

  # Extract basic coordinate information
  r = data['x1v']
  theta = data['x2v']
  phi = data['x3v']
  nx2 = len(theta)
  nx3 = len(phi)

  # Set radial extent
  if kwargs['r_max'] is not None:
    r_max = kwargs['r_max']
  else:
    r_max = data['x1f'][-1]

  # Perform slicing/averaging
  if kwargs['midplane']:
    if nx2%2 == 0:
      vals = np.mean(data[kwargs['quantity']][:,nx2/2-1:nx2/2+1,:], axis=1)
    else:
      vals = data[kwargs['quantity']][:,nx2/2,:]
    if kwargs['average']:
      vals = np.repeat(np.mean(vals, axis=0, keepdims=True), nx3, axis=0)
  else:
    if kwargs['average']:
      vals_right = np.mean(data[kwargs['quantity']], axis=0)
      vals_left = vals_right
    else:
      vals_right = \
          0.5 * (data[kwargs['quantity']][-1,:,:] + data[kwargs['quantity']][0,:,:])
      vals_left = 0.5 \
          * (data[kwargs['quantity']][nx3/2-1,:,:] + data[kwargs['quantity']][nx3/2,:,:])

  # Join data through boundaries
  if kwargs['midplane']:
    vals_extended = np.vstack((vals[-1:,:], vals, vals[:1,:]))
  else:
    vals_extended = \
        np.vstack((vals_left[0:1,:], vals_right, vals_right[::-1,:], vals_left[0:1,:]))

  # Create grids
  if kwargs['midplane']:
    phi_extended = \
        np.concatenate((phi[-1:]-2.0*np.pi, phi, phi[:1]+2.0*np.pi))
    r_grid,phi_grid = np.meshgrid(r, phi_extended)
    x_grid = r_grid * np.cos(phi_grid)
    y_grid = r_grid * np.sin(phi_grid)
  else:
    theta_extended = np.concatenate((-theta[0:1], theta, 2.0*np.pi-theta[::-1], \
        2.0*np.pi+theta[0:1]))
    r_grid,theta_grid = np.meshgrid(r, theta_extended)
    x_grid = r_grid * np.sin(theta_grid)
    y_grid = r_grid * np.cos(theta_grid)

  # Determine colormapping properties
  cmap = plt.get_cmap(kwargs['colormap'])
  vmin = kwargs['vmin']
  vmax = kwargs['vmax']
  if kwargs['log']:
    norm = colors.LogNorm()
  else:
    norm = colors.Normalize()

  # Make plot
  plt.figure()
  im = plt.pcolormesh(x_grid, y_grid, vals_extended, cmap=cmap, vmin=vmin, vmax=vmax, \
      norm=norm)
  plt.gca().set_aspect('equal')
  plt.xlim((-r_max, r_max))
  plt.ylim((-r_max, r_max))
  plt.colorbar(im)
  plt.savefig(kwargs['output_file'])

# Execute main function
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('data_file',
      type=str,
      help='name of input file, possibly including path')
  parser.add_argument('quantity',
      type=str,
      help='name of quantity to be plotted')
  parser.add_argument('output_file',
      type=str,
      help='name of output to be (over)written, possibly including path')
  parser.add_argument('-m', '--midplane',
      action='store_true',
      help='flag indicating plot should be midplane (r,phi) rather than (r,theta)')
  parser.add_argument('-a', '--average',
      action='store_true',
      help='flag indicating phi-averaging should be done')
  parser.add_argument('-l', '--level',
      type=int,
      default=None,
      help='refinement level to be used in plotting (default: max level in file)')
  parser.add_argument('-r', '--r_max',
      type=float,
      default=None,
      help='maximum radial extent of plot')
  parser.add_argument('-c', '--colormap',
      type=str,
      default=None,
      help='name of Matplotlib colormap to use instead of default; highly recommended; \
          try "RdBu_r" or "gist_heat" if looking for suggestions')
  parser.add_argument('--vmin',
      type=float,
      default=None,
      help='data value to correspond to colormap minimum; use --vmin=<val> if <val> has \
          negative sign')
  parser.add_argument('--vmax',
      type=float,
      default=None,
      help='data value to correspond to colormap maximum; use --vmax=<val> if <val> has \
          negative sign')
  parser.add_argument('--log',
      action='store_true',
      help='flag indicating data should be colormapped logarithmically')
  args = parser.parse_args()
  main(**vars(args))
