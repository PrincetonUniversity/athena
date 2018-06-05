#! /usr/bin/env python

"""
Script for finding optimal static mesh refinement grid in spherical coordinates.

Requires matplotlib if output other than "show" is specified. Requires scipy if
r_ratio is not specified.
"""

# Python standard modules
import argparse
import math

# Other Python modules
import numpy as np


# Main function
def main(**kwargs):

    # Extract inputs
    r_min = kwargs['r_min']
    r_max = kwargs['r_max']
    theta_min = kwargs['theta_min']
    theta_max = np.pi - theta_min
    num_r = kwargs['num_r']
    num_theta = kwargs['num_theta']
    num_phi = kwargs['num_phi']
    num_r_block = kwargs['num_r_block']
    num_theta_block = kwargs['num_theta_block']
    num_phi_block = kwargs['num_phi_block']
    max_levels = kwargs['max_levels']
    r_ratio = kwargs['r_ratio']
    metric = kwargs['metric']
    parameters = kwargs['parameters']
    theta_compress = kwargs['theta_compress']
    minimum_width = kwargs['minimum_width']
    output = kwargs['output']
    colormap = kwargs['colormap']
    grid_refined = kwargs['grid_refined']
    log = kwargs['log']

    # Verify inputs
    if not np.isfinite(r_min) or not np.isfinite(r_max) or r_min < 0.0 or r_max <= r_min:
        raise RuntimeError('r_min < r_max must be nonnegative numbers')
    if not np.isfinite(theta_min) or theta_min < 0.0 or theta_min >= np.pi/2.0:
        raise RuntimeError('must have 0 <= theta_min < pi/2')
    if (num_r % num_r_block != 0 or num_theta % num_theta_block != 0
            or num_phi % num_phi_block != 0):
        raise RuntimeError('blocks must evenly divide root grid in all dimensions')
    if num_r_block < 4 or num_theta_block < 4 or num_phi_block < 4:
        raise RuntimeError('blocks must have at least 4 cells in all dimensions')
    if max_levels < 0:
        raise RuntimeError('max_levels must be nonnegative')
    if (max_levels > 0 and
            (num_r_block % 2 != 0 or num_theta_block % 2 != 0 or num_phi_block % 2 != 0)):
        raise RuntimeError('blocks must have even number of cells in all dimensions to \
        support refinement')
    if r_ratio is not None and (not np.isfinite(r_ratio) or r_ratio <= 0.0):
        raise RuntimeError('r_ratio must be positive if specified')
    if r_ratio is not None and r_ratio < 1.0:
        raise RuntimeError('r_ratio < 1 not supported at this time')
    if metric is None and parameters is not None:
        raise RuntimeError('parameters cannot be specified without metric')
    if metric is not None and parameters is None:
        raise RuntimeError('must specify parameters for metric')
    if metric == 'schwarzschild':
        if len(parameters) != 1:
            raise RuntimeError('must specify 1 parameter (mass) for this metric')
        try:
            float(parameters[0])
        except ValueError:
            raise RuntimeError('invalid parameter')
    if metric == 'boyer-lindquist' or metric == 'kerr-schild':
        if len(parameters) != 2:
            raise RuntimeError(
                'must specify 2 parameters (mass and spin) for this metric')
        try:
            float(parameters[0])
            float(parameters[1])
        except ValueError:
            raise RuntimeError('invalid parameters')
    if theta_compress <= 0.0 or theta_compress > 1.0:
        raise RuntimeError('must have 0 < theta_compress <= 1')
    if minimum_width is not None and minimum_width <= 0.0:
        raise RuntimeError('must have minimum_width > 0')

    # Calculate minimum width
    if r_ratio is None:
        r_ratio = log_ratio(r_max/r_min, num_r)
    delta_phi = 2.0*np.pi / num_phi
    if minimum_width is None:
        r1 = r_min
        r2 = pos_face(r_min, r_max, r_ratio, num_r, 1)
        theta1 = theta_adjust(theta_min, theta_compress)
        theta2 = theta_adjust(
            pos_face(
                theta_min,
                theta_max,
                1.0,
                num_theta,
                1),
            theta_compress)
        w_r_min, w_theta_min, w_phi_min = widths(r1, r2, theta1, theta2,
                                                 delta_phi, metric, parameters)
        width_min = min(w_r_min, w_theta_min, w_phi_min)
    else:
        width_min = minimum_width

    # Determine refinement
    refinement = []
    r_bounds = []
    theta_bounds = []
    for L in range(max_levels+1):

        # Prepare description of all blocks at current level
        num_blocks_r = num_r/num_r_block * 2**L
        num_blocks_theta = num_theta/num_theta_block * 2**L
        r_bounds.append(np.empty(num_blocks_r+1))
        theta_bounds.append(np.empty(num_blocks_theta+1))
        for i in range(num_blocks_r+1):
            r_bounds[L][i] = pos_face(r_min, r_max, r_ratio**(1.0/2**L),
                                      num_r*2**L, i*num_r_block)
        for j in range(num_blocks_theta+1):
            theta_unadjusted = pos_face(theta_min, theta_max, 1.0,
                                        num_theta*2**L, j*num_theta_block)
            theta_bounds[L][j] = theta_adjust(theta_unadjusted, theta_compress)

        # Record which blocks might be refined (level 0) or might exist (otherwise)
        refinement.append(np.ones((num_blocks_r, num_blocks_theta), dtype=bool))
        if L == 0:
            continue
        for i in range(num_blocks_r/2):
            for j in range(num_blocks_theta/2):
                if not refinement[L-1][i, j]:
                    refinement[L][i*2:(i+1)*2, j*2:(j+1)*2] = False

        # Go through all blocks to see if they can exist based on cell widths
        for i in range(num_blocks_r):
            r1 = r_bounds[L][i]
            r2 = pos_face(r_bounds[L][i], r_bounds[L][i+1],
                          r_ratio**(1.0/2**L), num_r_block, 1)
            for j in range(num_blocks_theta):
                if not refinement[L][i, j]:
                    continue
                if j < (num_blocks_theta+1)/2:
                    theta1 = theta_adjust(theta_bounds[L][j], theta_compress)
                    theta2_unadjusted = pos_face(
                        theta_bounds[L][j], theta_bounds[L][j+1], 1.0, num_theta_block, 1)
                    theta2 = theta_adjust(theta2_unadjusted, theta_compress)
                else:
                    theta1_unadjusted = pos_face(theta_bounds[L][j], theta_bounds[L][j+1],
                                                 1.0, num_theta_block, num_theta_block-1)
                    theta1 = theta_adjust(theta1_unadjusted, theta_compress)
                    theta2 = theta_adjust(theta_bounds[L][j+1], theta_compress)
                w_r, w_theta, w_phi = widths(r1, r2, theta1, theta2, delta_phi/2**L,
                                             metric, parameters)
                if min(w_r, w_theta, w_phi) < width_min:
                    refinement[L][i, j] = False
                    refinement[L-1][i/2, j/2] = False

        # Make sure only entire blocks from previous level are refined
        for i in range(num_blocks_r/2):
            for j in range(num_blocks_theta/2):
                if not refinement[L-1][i, j]:
                    refinement[L][i*2:(i+1)*2, j*2:(j+1)*2] = False

        # Make sure blocks adjacent to unrefined ones cannot be refined
        refinement_copy = np.copy(refinement[L])
        for i in range(num_blocks_r):
            for j in range(num_blocks_theta):
                if not refinement_copy[i, j]:
                    continue
                ii_min = max(i-1, 0)
                ii_max = min(i+1, num_blocks_r-1)
                jj_min = max(j-1, 0)
                jj_max = min(j+1, num_blocks_theta-1)
                for ii in range(ii_min, ii_max+1):
                    for jj in range(jj_min, jj_max+1):
                        if not refinement_copy[ii, jj]:
                            refinement[L][i, j] = False

    # Determine refinement regions
    refinement_regions = []
    num_blocks = np.zeros(max_levels+1, dtype=int)
    num_blocks[0] = (
        (num_r/num_r_block) * (num_theta/num_theta_block) * (num_phi/num_phi_block))
    for L in range(max_levels):

        # Calculate number of blocks in this level
        num_blocks_r = num_r/num_r_block * 2**L
        num_blocks_theta = num_theta/num_theta_block * 2**L
        num_blocks_phi = num_phi/num_phi_block * 2**L

        # Find block limit of region to be refined
        j_lims = np.empty(num_blocks_r, dtype=int)
        for i in range(num_blocks_r):
            if not refinement[L][i, (num_blocks_theta-1)/2]:
                j_lims[i] = -1
                continue
            j_lim = 0
            for j in range((num_blocks_theta-1)/2, -1, -1):
                if not refinement[L][i, j]:
                    j_lim = j+1
                    break
            j_lims[i] = j_lim

        # Divide region into chunks
        previous_j_lim = -1
        i_start = None
        for i in range(num_blocks_r+1):
            if i == num_blocks_r:
                j_lim = -1
            else:
                j_lim = j_lims[i]
            if j_lim != previous_j_lim:
                if previous_j_lim >= 0:
                    refinement_regions.append((L, i_start, i-1, previous_j_lim))
                    num_refined_r = i - i_start
                    num_refined_theta = num_blocks_theta - 2 * previous_j_lim
                    num_refined_phi = num_blocks_phi
                    num_blocks_refined = num_refined_r*num_refined_theta*num_refined_phi
                    num_blocks[L] -= num_blocks_refined
                    num_blocks[L+1] += 8 * num_blocks_refined
                if j_lim >= 0:
                    i_start = i
            previous_j_lim = j_lim

    # Report refinement regions
    num_refinement = len(refinement_regions)
    if num_refinement == 0:
        print('\nNo refinement possible')
    elif num_refinement == 1:
        print('\nThe following input block can be used to specify maximal refinement:')
    else:
        print('\nThe following {} input blocks can be used to specify maximal refinement:'
              .format(num_refinement))
    for refinement_num, (L, i_start, i_end, j_lim) in enumerate(refinement_regions):
        num_blocks_r = num_r/num_r_block * 2**L
        num_blocks_theta = num_theta/num_theta_block * 2**L
        if i_start == 0:
            r1 = r_min
        else:
            r1 = pos_face(r_bounds[L][i_start], r_bounds[L][i_start+1],
                          r_ratio**(1.0/2**L), num_r_block, num_r_block/2)
        if i_end == num_blocks_r-1:
            r2 = r_max
        else:
            r2 = pos_face(r_bounds[L][i_end], r_bounds[L][i_end+1], r_ratio**(1.0/2**L),
                          num_r_block, num_r_block/2)
        if j_lim == 0:
            theta1 = theta_adjust(theta_min, theta_compress)
        else:
            theta1_unadjusted = pos_face(theta_bounds[L][j_lim], theta_bounds[L][j_lim+1],
                                         1.0, num_theta_block, num_theta_block/2)
            theta1 = theta_adjust(theta1_unadjusted, theta_compress)
        theta2 = np.pi - theta1
        print('\n<refinement{0}>'.format(refinement_num+1))
        print('level = ' + repr(L+1))
        print('x1min = ' + repr(r1))
        print('x1max = ' + repr(r2))
        print('x2min = ' + repr(theta1))
        print('x2max = ' + repr(theta2))
        print('x3min = ' + repr(0.0))
        print('x3max = ' + repr(2.0*np.pi))

    # Report number of blocks
    max_digits_level = len(repr(max_levels))
    max_digits_count = len(repr(sum(num_blocks)))
    print('\nNumber of blocks used:')
    for L in range(max_levels+1):
        digits_level = len(repr(L))
        # digits_count = len(repr(num_blocks[L]))  # unused
        string = ('    Level {0}:' + ' '*(max_digits_level-digits_level) + ' {1:'
                  + repr(max_digits_count) + 'd}')
        print(string.format(L, num_blocks[L]))
    print('    Total:  ' + ' '*max_digits_level + repr(sum(num_blocks)))

    # Report cell width
    print('\nLimiting width: {0:.3e}'.format(width_min))
    if minimum_width is None and width_min != w_phi_min:
        print('Note: phi-width not smallest width of this cell')

    # Report ratio information
    r_ratio_optimal = log_ratio(r_max/r_min, num_r)
    print('\nGeometric ratio used: x1rat = ' + repr(r_ratio))
    if r_ratio != r_ratio_optimal:
        print('Suggested ratio (so that r / Delta r is constant): '
              + repr(r_ratio_optimal))
    print('')

    # Create image of grid
    if output is not None:
        plot_grid(refinement, r_bounds, theta_bounds, output, colormap, grid_refined, log)


# Function for calculating geometric ratio closest to logarithmic spacing
def log_ratio(f, n):
    from scipy.optimize import brentq

    def res(ratio): return (
        1.0 + 1.0/(f-1.0) * (1.0 - f/ratio**(n-1)) * (ratio**n - 1.0) / (ratio - 1.0)
    )
    ratio_min = 1.0
    ratio_max = f**(1.0/(n-1))
    ratio1 = ratio_min + 0.5*(ratio_max-ratio_min)
    while res(ratio1) > 0.0:
        ratio1 -= 0.5*(ratio1-ratio_min)
    ratio2 = ratio_max
    ratio = brentq(res, ratio1, ratio2, xtol=1.0e-12)
    return ratio


# Function for calculating interface positions
def pos_face(x1, x2, ratio, n, n_face):
    ratio_powers = ratio ** np.arange(n)
    x = x1 + (x2-x1) * math.fsum(ratio_powers[:n_face]) / math.fsum(ratio_powers)
    return x


# Function for compressing theta-position
def theta_adjust(theta_unadjusted, theta_compress):
    theta_adjusted = (
        theta_unadjusted + (1.0-theta_compress)/2.0 * np.sin(2.0*theta_unadjusted))
    return theta_adjusted


# Function for calculating cell widths
def widths(r1, r2, theta1, theta2, delta_phi, metric, parameters):
    if metric == 'schwarzschild':
        m = float(parameters[0])
        r = (0.5 * (r1**3 + r2**3)) ** (1.0/3.0)
        theta = np.arccos(0.5 * (np.cos(theta1) + np.cos(theta2)))
        alpha1 = (1.0 - 2.0*m/r1) ** 0.5
        alpha2 = (1.0 - 2.0*m/r2) ** 0.5
        w_r = (r2 * alpha2 - r1 * alpha1
               + m * np.log((r2 * (1.0+alpha2) - m) / (r1 * (1.0+alpha1) - m)))
        w_theta = r * (theta2 - theta1)
        w_phi = r * np.sin(theta) * delta_phi
    if metric == 'boyer-lindquist':
        m = float(parameters[0])
        a = float(parameters[1])
        r = 0.5 * (r1 + r2)
        theta = 0.5 * (theta1 + theta2)
        delta1 = r1**2 - 2.0*m*r + a**2
        delta2 = r2**2 - 2.0*m*r + a**2
        sigma = r**2 + a**2 * np.cos(theta)**2
        w_r = (delta2**0.5 - delta1**0.5
               + m * np.log((r2 + delta2**0.5 - m) / (r1 + delta1**0.5 - m)))
        w_theta = r * (theta2 - theta1)
        w_phi = (np.sin(theta) * delta_phi
                 * (r**2 + a**2 + 2.0*m*a**2*r/sigma * np.sin(theta)**2) ** 0.5)
    elif metric == 'kerr-schild':
        m = float(parameters[0])
        a = float(parameters[1])
        r = 0.5 * (r1 + r2)
        theta = 0.5 * (theta1 + theta2)
        sigma = r**2 + a**2 * np.cos(theta)**2
        w_r = ((r2**2 + m**2) ** 0.5 - (r1**2 + m**2) ** 0.5
               + m * np.log(((r2**2 + m**2) ** 0.5 + r2) / ((r1**2 + m**2) ** 0.5 + r1)))
        w_theta = r * (theta2 - theta1)
        w_phi = (np.sin(theta) * delta_phi
                 * (r**2 + a**2 + 2.0*m*a**2*r/sigma * np.sin(theta)**2) ** 0.5)
    else:
        r = (0.5 * (r1**3 + r2**3)) ** (1.0/3.0)
        theta = np.arccos(0.5 * (np.cos(theta1) + np.cos(theta2)))
        w_r = r2 - r1
        w_theta = r * (theta2 - theta1)
        w_phi = r * np.sin(theta) * delta_phi
    return (w_r, w_theta, w_phi)


# Function for making image of grid
def plot_grid(refinement, r_bounds, theta_bounds, output, colormap, grid_refined, log):

    # Load Python plotting modules
    import matplotlib
    if output != 'show':
        matplotlib.use('agg')
    matplotlib.rc('font', family='serif', size=14)
    import matplotlib.cm as cm
    import matplotlib.colors as colors
    import matplotlib.patches as patches
    import matplotlib.pyplot as plt

    # Locate corners of cells
    eps = 1.0e-6
    r = r_bounds[-1]
    delta_r = np.ediff1d(r)
    r_in = r - eps * np.concatenate((delta_r[:1], delta_r))
    r_out = r + eps * np.concatenate((delta_r, delta_r[-1:]))
    r_vals = np.vstack((r_in, r_out)).T.flatten()
    if log:
        r_vals = np.log10(r_vals)
    theta = theta_bounds[-1]
    delta_theta = np.ediff1d(theta)
    if theta_bounds[0][0] == 0.0:
        theta_north = theta[:-1] + eps * delta_theta
        theta_south = theta[1:] - eps * delta_theta
        theta_vals = np.vstack((theta_north, theta_south)).T.flatten()
        theta_vals = np.concatenate(
            (-theta_vals[0:1], theta_vals,
             2.0*np.pi-theta_vals[::-1], 2.0*np.pi+theta_vals[0:1])
        )
    else:
        theta_north = theta[1:] - eps * delta_theta
        theta_south = theta[:-1] + eps * delta_theta
        theta_north_end = max(theta[0] - eps * delta_theta[0], 0.0)
        theta_south_end = min(theta[-1] + eps * delta_theta[-1], np.pi)
        theta_north = np.concatenate(([theta_north_end], theta_north))
        theta_south = np.concatenate((theta_south, [theta_south_end]))
        theta_vals = np.vstack((theta_north, theta_south)).T.flatten()
        theta_vals = np.concatenate((theta_vals, 2.0*np.pi-theta_vals[::-1]))
    r_grid, theta_grid = np.meshgrid(r_vals, theta_vals)
    x_grid = r_grid * np.sin(theta_grid)
    y_grid = r_grid * np.cos(theta_grid)

    # Determine levels as a function of position
    max_level = len(refinement) - 1
    levels = np.zeros((2*len(r), 2*len(delta_theta)), dtype=int)
    for l, refinement_current in enumerate(refinement[:-1]):
        block_size = 2 ** (max_level-l+1)
        for (i, j), refined in np.ndenumerate(refinement_current):
            if refined:
                levels[i*block_size+1:(i+1)*block_size+1,
                       j*block_size:(j+1)*block_size] = l + 1
    levels[0, :] = -1
    levels[-1, :] = -1
    levels = levels.T
    if theta_bounds[0][0] == 0.0:
        levels = np.vstack((levels[0:1, :], levels, levels[::-1, :], levels[0:1, :]))
    else:
        boundary_row = -1 * np.ones(2*len(r), dtype=int)
        levels = np.vstack((boundary_row, levels, boundary_row, boundary_row,
                            levels[::-1, :], boundary_row))

    # Set discrete colormap
    bounds = range(max_level+2)
    norm = colors.BoundaryNorm(bounds, 256)
    if max_level == 0:
        cmap = cm.get_cmap(colormap)
        color = cmap(0.0)
        colormap, norm = colors.from_levels_and_colors(range(3), (color, color))

    # Plot colors to indicate level
    plt.figure()
    ax = plt.gca()
    im = plt.pcolormesh(x_grid, y_grid, levels, cmap=colormap, norm=norm)
    im.cmap.set_under('w')
    ax.set_aspect('equal')
    rlim = r[-1]
    if log:
        rlim = np.log10(rlim)
    plt.xlim((-rlim, rlim))
    plt.ylim((-rlim, rlim))
    if log:
        plt.xlabel(r'$\log_{10}(r)\ \sin(\theta)$')
        plt.ylabel(r'$\log_{10}(r)\ \cos(\theta)$')
    else:
        plt.xlabel(r'$r\ \sin(\theta)$')
        plt.ylabel(r'$r\ \cos(\theta)$')

    # Make colorbar
    ticks = np.arange(max_level+1) + 0.5
    if max_level == 0:
        ticks = (1,)
    labels = [repr(n) for n in range(max_level+1)]
    cbar = plt.colorbar(im, ticks=ticks)
    cax = cbar.ax
    cax.set_yticklabels(labels)
    cax.tick_params(length=0)

    # Draw root grid block boundaries
    if grid_refined < 0:
        r_vals_grid = (r_bounds[0][0], r_bounds[0][-1])
        theta_vals_grid = (theta_bounds[0][0], theta_bounds[0][-1])
    else:
        r_vals_grid = r_bounds[0]
        theta_vals_grid = theta_bounds[0]
    for r_val in r_vals_grid:
        if log:
            r_val_plot = np.log10(r_val)
        else:
            r_val_plot = r_val
        theta1 = 90.0 - theta_bounds[0][0] * 180.0/np.pi
        theta2 = 90.0 - theta_bounds[0][-1] * 180.0/np.pi
        theta3 = 180.0 - theta1
        theta4 = 180.0 - theta2
        arc = patches.Arc((0, 0), r_val_plot*2, r_val_plot *
                          2, theta1=theta2, theta2=theta1)
        ax.add_artist(arc)
        arc = patches.Arc((0, 0), r_val_plot*2, r_val_plot *
                          2, theta1=theta3, theta2=theta4)
        ax.add_artist(arc)
    for theta_val in theta_vals_grid:
        r_inner = r_bounds[0][0]
        r_outer = r_bounds[0][-1]
        if log:
            r_inner = np.log10(r_inner)
            r_outer = np.log10(r_outer)
        x1 = r_inner * np.sin(theta_val)
        y1 = r_inner * np.cos(theta_val)
        x2 = r_outer * np.sin(theta_val)
        y2 = r_outer * np.cos(theta_val)
        plt.plot((x1, x2), (y1, y2), 'k')
        plt.plot((-x1, -x2), (y1, y2), 'k')

    # Draw refined boundaries
    for l in range(min(max_level, grid_refined)):
        for (i, j), refined in np.ndenumerate(refinement[l]):
            if refined:
                r_val = r_bounds[l+1][2*i+1]
                if log:
                    r_val = np.log10(r_val)
                theta1 = 90.0 - theta_bounds[l][j] * 180.0/np.pi
                theta2 = 90.0 - theta_bounds[l][j+1] * 180.0/np.pi
                theta3 = 180.0 - theta1
                theta4 = 180.0 - theta2
                arc = patches.Arc((0, 0), r_val*2, r_val*2, theta1=theta2, theta2=theta1)
                ax.add_artist(arc)
                arc = patches.Arc((0, 0), r_val*2, r_val*2, theta1=theta3, theta2=theta4)
                ax.add_artist(arc)
                theta_val = theta_bounds[l+1][2*j+1]
                r1 = r_bounds[l][i]
                r2 = r_bounds[l][i+1]
                if log:
                    r1 = np.log10(r1)
                    r2 = np.log10(r2)
                x1 = r1 * np.sin(theta_val)
                y1 = r1 * np.cos(theta_val)
                x2 = r2 * np.sin(theta_val)
                y2 = r2 * np.cos(theta_val)
                plt.plot((x1, x2), (y1, y2), 'k')
                plt.plot((-x1, -x2), (y1, y2), 'k')

    # Show or save figure
    if output == 'show':
        plt.show()
    else:
        plt.savefig(output, bbox_inches='tight')


# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('r_min',
                        type=float,
                        help='minimum radius')
    parser.add_argument('r_max',
                        type=float,
                        help='maximum radius')
    parser.add_argument(
        'theta_min',
        type=float,
        help=('minimum polar angle (maximum is assumed to be pi-complement); polar '
              'boundaries assumed if this is 0')
    )
    parser.add_argument('num_r',
                        type=int,
                        help='number of cells in radial direction')
    parser.add_argument('num_theta',
                        type=int,
                        help='number of cells in polar direction')
    parser.add_argument('num_phi',
                        type=int,
                        help='number of cells in azimuthal direction')
    parser.add_argument('num_r_block',
                        type=int,
                        help='number of cells in radial direction in one block')
    parser.add_argument('num_theta_block',
                        type=int,
                        help='number of cells in polar direction in one block')
    parser.add_argument('num_phi_block',
                        type=int,
                        help='number of cells in azimuthal direction in one block')
    parser.add_argument('max_levels',
                        type=int,
                        help='maximum number of mesh refinement levels to consider')
    parser.add_argument(
        '-r',
        '--r_ratio',
        type=float,
        help='ratio of adjacent separations in radius (optimal value chosen if omitted)')
    parser.add_argument('-m', '--metric',
                        choices=('schwarzschild', 'boyer-lindquist', 'kerr-schild',),
                        default=None,
                        help='metric to assume if in GR')
    parser.add_argument(
        '-p', '--parameters',
        nargs='+',
        help=('parameters (mass M, possibly spin a, 0 <= a < M) to be used if metric is '
              'specified')
    )
    parser.add_argument(
        '--theta_compress',
        type=float,
        default=1.0,
        help='parameter h governing midplane compression of theta-surfaces')
    parser.add_argument('--minimum_width',
                        type=float,
                        help='override for smallest allowed cell width')
    parser.add_argument(
        '-o', '--output',
        help=('name of image file to write showing grid; use "show" to show interactive '
              'plot instead')
    )
    parser.add_argument(
        '-c', '--colormap',
        default='cool',
        help='name of colormap')
    parser.add_argument(
        '-g',
        '--grid_refined',
        type=int,
        default=0,
        help='maximum refinement level at which grid of block boundaries should be drawn')
    parser.add_argument(
        '-l',
        '--log',
        action='store_true',
        help='flag indicating output image should show radius logarithmically')
    args = parser.parse_args()
    main(**vars(args))
