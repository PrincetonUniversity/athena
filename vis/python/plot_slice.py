#! /usr/bin/env python

"""
Script for plotting 2D data or 2D slices of 3D data, intended primarily for
Cartesian grids.

Run "plot_slice.py -h" to see description of inputs.

See documentation on athena_read.athdf() for important notes about reading files
with mesh refinement.

Users are encouraged to make their own versions of this script for improved
results by adjusting figure size, spacings, tick locations, axes labels, etc.
The script must also be modified to plot any functions of the quantities in the
file, including combinations of multiple quantities.
"""

# Python standard modules
import argparse
import warnings

# Other Python modules
import numpy as np

# Athena++ modules
import athena_read


# Main function
def main(**kwargs):

    # Load Python plotting modules
    if kwargs['output_file'] != 'show':
        import matplotlib
        matplotlib.use('agg')
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors

    # Verify user inputs
    slice_erroneously_specified = False
    if kwargs['slice_location'] is not None:
        if kwargs['stream'] is not None:
            if (kwargs['average'] or kwargs['sum']) and kwargs['stream_average']:
                slice_erroneously_specified = True
        else:
            if kwargs['average'] or kwargs['sum']:
                slice_erroneously_specified = True
    if slice_erroneously_specified:
        raise RuntimeError('Slice location specified but all quantities are to be '
                           'averaged or summed')

    # Set default slice location (even if averaging or summing)
    if kwargs['slice_location'] is None:
        kwargs['slice_location'] = 0.0

    # Determine refinement level to use
    if kwargs['level'] is not None:
        level = kwargs['level']
    else:
        level = None

    # Determine if vector quantities should be read
    quantities = [kwargs['quantity']]
    if kwargs['stream'] is not None:
        if kwargs['direction'] == 1:
            quantities.append(kwargs['stream'] + '2')
            quantities.append(kwargs['stream'] + '3')
        elif kwargs['direction'] == 2:
            quantities.append(kwargs['stream'] + '1')
            quantities.append(kwargs['stream'] + '3')
        else:
            quantities.append(kwargs['stream'] + '1')
            quantities.append(kwargs['stream'] + '2')

    # Read data
    if quantities[0] == 'Levels':
        data = athena_read.athdf(kwargs['data_file'], quantities=quantities[1:],
                                 level=level, return_levels=True, num_ghost=kwargs['num_ghost'])
    else:
        data = athena_read.athdf(kwargs['data_file'], quantities=quantities, level=level, num_ghost=kwargs['num_ghost'])

    # Check that coordinates work with user choices
    coordinates = data['Coordinates'].decode('ascii', 'replace')
    ave_or_sum = kwargs['average'] or kwargs['sum'] or kwargs['stream_average']
    warn_projection = False
    warn_vector = False
    projection_type = None
    if coordinates in ('cartesian', 'minkowski'):
        pass
    elif coordinates == 'cylindrical':
        if ave_or_sum and kwargs['direction'] == 1:
            warn_projection = True
        if kwargs['stream'] and kwargs['direction'] in (1, 3):
            warn_vector = True
        if kwargs['direction'] == 3:
            # pass
            projection_type = "polar"
    elif coordinates in ('spherical_polar', 'schwarzschild', 'kerr-schild'):
        if ave_or_sum and kwargs['direction'] in (1, 2):
            warn_projection = True
        if kwargs['stream']:
            warn_vector = True
    else:
        warnings.warn('Coordinates not recognized; results may be misleading')
    if warn_projection:
        warnings.warn('Sums/slices are not computed with correct volumes')
    if warn_vector:
        warnings.warn('Vector plot may be misleading')

    # Name coordinates
    if coordinates in ('cartesian', 'minkowski'):
        coord_labels = (r'$x$', r'$y$', r'$z$')
    elif coordinates == 'cylindrical':
        coord_labels = (r'$R$', r'$\phi$', r'$z$')
    elif coordinates in ('spherical_polar', 'schwarzschild', 'kerr-schild'):
        coord_labels = (r'$r$', r'$\theta$', r'$\phi$')
    else:
        coord_labels = (r'$x^1$', r'$x^2$', r'$x^3$')

    # Extract basic coordinate information
    if kwargs['direction'] == 1:
        xf = data['x2f']
        xv = data['x2v']
        yf = data['x3f']
        yv = data['x3v']
        zf = data['x1f']
        x_label = coord_labels[1]
        y_label = coord_labels[2]
    elif kwargs['direction'] == 2:
        xf = data['x1f']
        xv = data['x1v']
        yf = data['x3f']
        yv = data['x3v']
        zf = data['x2f']
        x_label = coord_labels[0]
        y_label = coord_labels[2]
    if kwargs['direction'] == 3:
        xf = data['x1f']
        xv = data['x1v']
        yf = data['x2f']
        yv = data['x2v']
        zf = data['x3f']
        x_label = coord_labels[0]
        y_label = coord_labels[1]

    # Create grids
    x_grid, y_grid = np.meshgrid(xf, yf)
    x_stream, y_stream = np.meshgrid(xv, yv)

    # Extract scalar data
    vals = data[kwargs['quantity']]
    if kwargs['average'] or kwargs['sum']:
        vals = np.sum(vals, axis=3-kwargs['direction'])
        if kwargs['sum']:
            vals *= zf[-1] - zf[0]
        vals /= len(zf) - 1
    else:
        if kwargs['slice_location'] < zf[0]:
            index = 0
        elif kwargs['slice_location'] >= zf[-1]:
            index = -1
        else:
            index = np.where(zf <= kwargs['slice_location'])[0][-1]
        if kwargs['direction'] == 1:
            vals = vals[:, :, index]
        elif kwargs['direction'] == 2:
            vals = vals[:, index, :]
        else:
            vals = vals[index, :, :]

    # Extract vector data
    if kwargs['stream'] is not None:
        if kwargs['direction'] == 1:
            vals_x = data[kwargs['stream'] + '2']
            vals_y = data[kwargs['stream'] + '3']
        elif kwargs['direction'] == 2:
            vals_x = data[kwargs['stream'] + '1']
            vals_y = data[kwargs['stream'] + '3']
        else:
            vals_x = data[kwargs['stream'] + '1']
            vals_y = data[kwargs['stream'] + '2']
        if kwargs['stream_average']:
            vals_x = np.sum(vals_x, axis=3-kwargs['direction'])
            vals_y = np.sum(vals_y, axis=3-kwargs['direction'])
            vals_x /= len(zf) - 1
            vals_y /= len(zf) - 1
        else:
            if kwargs['slice_location'] < zf[0]:
                index = 0
            elif kwargs['slice_location'] >= zf[-1]:
                index = -1
            else:
                index = np.where(zf <= kwargs['slice_location'])[0][-1]
            if kwargs['direction'] == 1:
                vals_x = vals_x[:, :, index]
                vals_y = vals_y[:, :, index]
            elif kwargs['direction'] == 2:
                vals_x = vals_x[:, index, :]
                vals_y = vals_y[:, index, :]
            else:
                vals_x = vals_x[index, :, :]
                vals_y = vals_y[index, :, :]

    # Determine plot limits
    x_min = kwargs['x_min'] if kwargs['x_min'] is not None else xf[0]
    x_max = kwargs['x_max'] if kwargs['x_max'] is not None else xf[-1]
    y_min = kwargs['y_min'] if kwargs['y_min'] is not None else yf[0]
    y_max = kwargs['y_max'] if kwargs['y_max'] is not None else yf[-1]
    v_min = kwargs['vmin'] if kwargs['vmin'] is not None else vals.min()
    v_max = kwargs['vmax'] if kwargs['vmax'] is not None else vals.max()

    # Determine colormap norm
    if kwargs['logc']:
        norm = colors.LogNorm(v_min, v_max)
    else:
        norm = colors.Normalize(v_min, v_max)

    # Make plot
    # should make the size editable?
    fig = plt.figure(1, figsize=(12, 12))
    ax = fig.add_subplot(1,1,1,projection=projection_type)
    if projection_type == 'polar':
        # switch axis for radial and azimuthal
        im = ax.pcolormesh(y_grid, x_grid, vals, cmap=kwargs['colormap'], norm=norm)
    else:
        im = ax.pcolormesh(x_grid, y_grid, vals, cmap=kwargs['colormap'], norm=norm)

    if kwargs['stream'] is not None:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                'ignore',
                'invalid value encountered in greater_equal',
                RuntimeWarning,
                'numpy')
        if projection_type == 'polar':
            # switch axis for radial and azimuthal and Transpose
            ax.streamplot(y_stream.T, x_stream.T, vals_y.T, vals_x.T,
                        density=kwargs['stream_density'], color='k')
        else:
            ax.streamplot(x_stream, y_stream, vals_x, vals_y,
                        density=kwargs['stream_density'], color='k')

    if projection_type == 'polar':
        ax.set_rmin(x_min)
        ax.set_rmax(x_max)
    else:
        ax.set_xlim((x_min, x_max))
        ax.set_ylim((y_min, y_max))
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)

    fig.colorbar(im)

    if not kwargs['fill']:
        ax.set_aspect('equal')

    if kwargs['output_file'] == 'show':
        fig.show()
    else:
        fig.savefig(kwargs['output_file'], bbox_inches='tight')


# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('data_file',
                        help='name of input file, possibly including path')
    parser.add_argument('quantity',
                        help='name of quantity to be plotted')
    parser.add_argument('output_file',
                        help=('name of output to be (over)written, possibly including '
                              'path; use "show" to show interactive plot instead'))
    parser.add_argument('-d', '--direction',
                        type=int,
                        choices=(1, 2, 3),
                        default=3,
                        help=('direction orthogonal to slice for 3D data'))
    parser.add_argument('--slice_location',
                        type=float,
                        default=None,
                        help=('coordinate value along which slice is to be taken '
                              '(default: 0)'))
    parser.add_argument('-a', '--average',
                        action='store_true',
                        help=('flag indicating averaging should be done in orthogonal '
                              'direction for 3D data'))
    parser.add_argument('-s', '--sum',
                        action='store_true',
                        help=('flag indicating summation should be done in orthogonal '
                              'direction for 3D data'))
    parser.add_argument('-l',
                        '--level',
                        type=int,
                        default=None,
                        help=('refinement level to be used in plotting (default: max '
                              'level in file)'))
    parser.add_argument('--x_min',
                        type=float,
                        default=None,
                        help='minimum extent of plot in first plotted direction')
    parser.add_argument('--x_max',
                        type=float,
                        default=None,
                        help='maximum extent of plot in first plotted direction')
    parser.add_argument('--y_min',
                        type=float,
                        default=None,
                        help='minimum extent of plot in second plotted direction')
    parser.add_argument('--y_max',
                        type=float,
                        default=None,
                        help='maximum extent of plot in second plotted direction')
    parser.add_argument('-f', '--fill',
                        action='store_true',
                        help='flag indicating image should fill plot area, even if this '
                             'distorts the aspect ratio')
    parser.add_argument('-c',
                        '--colormap',
                        default=None,
                        help=('name of Matplotlib colormap to use instead of default'))
    parser.add_argument('--vmin',
                        type=float,
                        default=None,
                        help=('data value to correspond to colormap minimum; use '
                              '--vmin=<val> if <val> has negative sign'))
    parser.add_argument('--vmax',
                        type=float,
                        default=None,
                        help=('data value to correspond to colormap maximum; use '
                              '--vmax=<val> if <val> has negative sign'))
    parser.add_argument('--logc',
                        action='store_true',
                        help='flag indicating data should be colormapped logarithmically')
    parser.add_argument('--stream',
                        default=None,
                        help='name of vector quantity to use to make stream plot')
    parser.add_argument('--stream_average',
                        action='store_true',
                        help='flag indicating stream plot should be averaged in '
                             'orthogonal direction for 3D data')
    parser.add_argument('--stream_density',
                        type=float,
                        default=1.0,
                        help='density of stream lines')
    parser.add_argument('--num_ghost',
                        type=int,
                        default=0,
                        help=('Include number of ghost cells in each direction'))
    args = parser.parse_args()
    main(**vars(args))
