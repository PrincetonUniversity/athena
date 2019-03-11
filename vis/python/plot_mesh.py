#! /usr/bin/env python

"""
Script for plotting mesh structure in mesh_structure.dat (default) file produced
by running Athena++ with "-m <np>" argument.

Can optionally specify "-i <input_file>" and/or "-o <output_file>". Output
defaults to using "show()" command rather than saving to file.
"""

# Python modules
import argparse


# Main function
def main(**kwargs):

    # Extract inputs
    input_file = kwargs['input']
    output_file = kwargs['output']

    # Load Python plotting modules
    if output_file != 'show':
        import matplotlib
        matplotlib.use('agg')
    import matplotlib.pyplot as plt
    # not used explicitly, but required for 3D projections
    from mpl_toolkits.mplot3d import Axes3D  # noqa

    # Read and plot block edges
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    x = []
    y = []
    z = []
    with open(input_file) as f:
        for line in f:
            if line[0] != '\n' and line[0] != '#':
                numbers_str = line.split()
                x.append(float(numbers_str[0]))
                y.append(float(numbers_str[1]))
                # append zero if 2D
                if(len(numbers_str) > 2):
                    z.append(float(numbers_str[2]))
                else:
                    z.append(0.0)
            if line[0] == '\n' and len(x) != 0:
                ax.plot(x, y, z, 'k-')
                x = []
                y = []
                z = []
    if output_file == 'show':
        plt.show()
    else:
        plt.savefig(output_file, bbox_inches='tight')


# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input',
                        default='mesh_structure.dat',
                        help='name of mesh structure file')
    parser.add_argument('-o',
                        '--output',
                        default='show',
                        help=('name of output image file to create; omit to '
                              'display rather than save image'))
    args = parser.parse_args()
    main(**vars(args))
