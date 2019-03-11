#! /usr/bin/env python

"""
Script for plotting 1D data from .athdf, .hst, or .tab files.

Run "plot_lines.py -h" to see full description of inputs.

Multiple lines can be plotted by having any of the first three arguments be
comma-separated lists. If one list runs out before another, its last entry will
be repeated as necessary.

Use "show" for the 4th argument to show interactive plot instead of saving to
file.
"""

# Python modules
import argparse

# Athena++ modules
import athena_read


# Main function
def main(**kwargs):

    # Extract inputs
    data_files = kwargs['data_files'].split(',')
    x_names = kwargs['x_names'].split(',')
    y_names = kwargs['y_names'].split(',')
    output_file = kwargs['output_file']
    styles = kwargs['styles'].split(',')
    colors = kwargs['colors']
    labels = kwargs['labels']
    x_log = kwargs['x_log']
    y_log = kwargs['y_log']
    x_min = kwargs['x_min']
    x_max = kwargs['x_max']
    y_min = kwargs['y_min']
    y_max = kwargs['y_max']
    x_label = kwargs['x_label']
    y_label = kwargs['y_label']

    # Verify inputs
    num_lines = max(len(data_files), len(x_names), len(y_names))
    if data_files[0] == '':
        raise RuntimeError('First entry in data_files must be nonempty')
    if x_names[0] == '':
        raise RuntimeError('First entry in x_names must be nonempty')
    if y_names[0] == '':
        raise RuntimeError('First entry in y_names must be nonempty')
    if len(data_files) < num_lines:
        data_files += data_files[-1:] * (num_lines - len(data_files))
    if len(x_names) < num_lines:
        x_names += x_names[-1:] * (num_lines - len(x_names))
    if len(y_names) < num_lines:
        y_names += y_names[-1:] * (num_lines - len(y_names))
    for n in range(num_lines):
        if data_files[n] == '':
            data_files[n] = data_files[n-1]
        if x_names[n] == '':
            x_names[n] = x_names[n-1]
        if y_names[n] == '':
            y_names[n] = y_names[n-1]
    for data_file in data_files:
        valid_file = (data_file[-6:] == '.athdf' or data_file[-4:] == '.hst'
                      or data_file[-4:] == '.tab')
        if not valid_file:
            raise RuntimeError('Files must have .athdf, .hst, or .tab extension')
    if len(styles) < num_lines:
        styles += styles[-1:] * (num_lines - len(styles))
    for n in range(num_lines):
        styles[n] = styles[n].lstrip()
        if styles[n] == '':
            styles[n] = '-'
    if colors is None:
        colors = [None] * num_lines
    else:
        colors = colors.split(',')
        if len(colors) < num_lines:
            colors += colors[-1:] * (num_lines - len(colors))
        for n in range(num_lines):
            if colors[n] == '':
                colors[n] = None
    if num_lines == 1 and colors[0] is None:
        colors[0] = 'k'
    if labels is None:
        labels = [None] * num_lines
    else:
        labels = labels.split(',')
        if len(labels) < num_lines:
            labels += [None] * (num_lines - len(labels))
        for n in range(num_lines):
            if labels[n] == '':
                labels[n] = None
    labels_used = False
    for n in range(num_lines):
        if labels[n] is not None:
            labels_used = True
            break

    # Load Python plotting modules
    if output_file != 'show':
        import matplotlib
        matplotlib.use('agg')
    import matplotlib.pyplot as plt

    # Read data
    x_vals = []
    y_vals = []
    for n in range(num_lines):
        if data_files[n][-6:] == '.athdf':
            data = athena_read.athdf(data_files[n])
        elif data_files[n][-4:] == '.hst':
            data = athena_read.hst(data_files[n])
        else:
            data = athena_read.tab(data_files[n])
        x_vals.append(data[x_names[n]].flatten())
        y_vals.append(data[y_names[n]].flatten())

    # Plot data
    plt.figure()
    for n in range(num_lines):
        plt.plot(x_vals[n], y_vals[n], styles[n], color=colors[n], label=labels[n])
    if x_log:
        plt.xscale('log')
    if y_log:
        plt.yscale('log')
    plt.xlim((x_min, x_max))
    plt.ylim((y_min, y_max))
    if x_label is not None:
        plt.xlabel(x_label)
    if y_label is not None:
        plt.ylabel(y_label)
    if labels_used:
        plt.legend(loc='best')
    if output_file == 'show':
        plt.show()
    else:
        plt.savefig(output_file, bbox_inches='tight')


# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
      'data_files',
      help=('comma-separated list of input files; empty strings repeat previous entries; '
            'list is extended if x_names or y_names is longer')
    )
    parser.add_argument(
        'x_names',
        help=('comma-separated list of abscissas; empty strings repeat previous entries; '
              'list is extended if data_files or y_names is longer')
    )
    parser.add_argument(
        'y_names',
        help=('comma-separated list of ordinates; empty strings repeat previous entries; '
              'list is extended if data_files or x_names is longer')
    )
    parser.add_argument(
        'output_file',
        help=('name of output to be (over)written; use "show" to show interactive plot '
              'instead')
    )
    parser.add_argument(
      '-s', '--styles',
      default='-',
      help=('comma-separated list of line or marker styles, such as "-" or "o"; use the '
            '" -s=..." form of the argument if the first entry begins with a dash; empty '
            'strings are interpreted as solid lines; last entry is repeated as necessary')
    )
    parser.add_argument(
      '-c', '--colors',
      help=('comma-separated list of color codes, such as "k", "blue", or "#123abc"; '
            'empty strings result in black (single line) or default color cycling '
            '(multiple lines); last entry is repeated as necessary')
    )
    parser.add_argument(
      '-l', '--labels',
      help=('comma-separated list of labels for legend; empty strings are not added to '
            'legend; strings can include mathematical notation inside $...$ (e.g. "-l '
            '\'$\\rho$\'")')
    )
    parser.add_argument(
        '--x_log',
        action='store_true',
        help='flag indicating x-axis should be log scaled')
    parser.add_argument(
        '--y_log',
        action='store_true',
        help='flag indicating y-axis should be log scaled')
    parser.add_argument('--x_min',
                        type=float,
                        help='minimum for x-axis')
    parser.add_argument('--x_max',
                        type=float,
                        help='maximum for x-axis')
    parser.add_argument('--y_min',
                        type=float,
                        help='minimum for y-axis')
    parser.add_argument('--y_max',
                        type=float,
                        help='maximum for y-axis')
    parser.add_argument('--x_label',
                        help='label to use for x-axis')
    parser.add_argument('--y_label',
                        help='label to use for y-axis')
    args = parser.parse_args()
    main(**vars(args))
