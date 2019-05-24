import numpy as np
import os
import sys
import logging
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)
import matplotlib               # noqa
# matplotlib.use('agg')
matplotlib.use('pdf')
import matplotlib.pyplot as plt  # noqa
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter  # noqa
from matplotlib.lines import Line2D  # noqa
sys.path.insert(0, '../../vis/python')
import athena_read                             # noqa
athena_read.check_nan_flag = True

# TODO(felker): suppress "logging" module's DEBUG statemets from matplotlib
# TODO(felker): improve coupling with calling script, ../mignone_radial_1d.py
# Plotting wrapper fns should take parameters to avoid having to match below
# file-level scoped variables with the ones in../mignone_radial_1d.py

# Problem parameters: (Mignone formulation and Athena++ options)
alpha = 1
iprob = 1
integrator = 'rk3'
coords = ['cylindrical', 'spherical_polar']
cases = ['A', 'B']
coord_m = [1, 2]
a_params = [10.0, 16.0]
b_params = [0.0, 0.5]
case_parameters = ['{a=10, b=0}', '{a=16, b=0.5}']
xorders = [2, 3]
xorder_strs = ['PLM', r'$\mathrm{PPM}_{4}$']

nx1_profile = 64   # resolution of Athena++ output to compare w/ analytic solution
nsamples = 1000    # x1 samples for analytic solution


def InitialGaussianProfile(x1, a, b):
    return np.exp(-a**2*(x1-b)**2)


def EvolvedGaussianProfile(x1, a, b, m, t):
    x_initial = x1*np.exp(-alpha*t)
    q_initial = InitialGaussianProfile(x_initial, a, b)
    amp = np.exp(-(m + 1)*alpha*t)
    return amp*q_initial


# Plot appearance options:
figsize = (12.8, 9.6)
dpi_global = 300

case_xlims = [
    [0.0, 1.0],
    [0.82, 1.9],
]

coord_ylims = [
    [0.0, 0.14],
    [0.0, 0.055],
]

ylims = [
    [4e-12, 2e-3],
    [1e-7, 1e-1],
    [1e-14, 5e-4],
    [3e-8, 4e-2],
]

major_yticks = [
    [1e-10, 1e-8, 1e-6, 1e-4],
    [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1],
    [1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4],
    [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2],
]

minor_yticks = [
    [1e-11, 1e-9, 1e-7, 1e-5, 1e-3],
    None,
    [1e-13, 1e-11, 1e-9, 1e-7, 1e-5],
    None,
]

minor_yticks_auto = [
    False,
    True,
    False,
    False,
]

xorder_symbols = {
    2: "+",   # P = filled variant
    # weno3 : 'o',
    # ppm0 : '' asterisk?
    # ppm3 : '^',
    3: "s",
    # ppm5 : 'D',
}

xorder_colors = {
    2: "teal",  # 317465  # rgb(49, 116, 101)
    # weno3 : 'g',  # #3e8a28, rgb(62, 138, 40)
    # ppm0 : 'k' asterisk?
    # ppm3 : "# e5b03c",  # rgb(229, 176, 60),
    3: "#cb5b25",  # rgb(203, 91, 37)
    # ppm5 : 'r',    # #a62b17, rgb(166, 43, 23)
}

profile_legend_handles = [
    Line2D([0], [0], marker=xorder_symbols[xorder_], color='w',
           label=xorder_str_, markeredgecolor=xorder_colors[xorder_],
           markerfacecolor='w', fillstyle='none', markersize=8)
    for xorder_, xorder_str_ in zip(xorders, xorder_strs)]

profile_legend_handles += [Line2D([0], [0], color='k', label='Exact', linewidth=1.0)]


def figure2_profiles():
    fig = plt.figure(figsize=(1.0*figsize[0], 1.0*figsize[1]), dpi=dpi_global)
    axes = fig.subplots(2, 2, gridspec_kw={'hspace': 0.4})

    for coord_, ylims_, m_, axes_row_ in zip(coords, coord_ylims, coord_m, axes):
        for case_, xlims_, param_str_, a_, b_, ax in zip(cases, case_xlims,
                                                         case_parameters, a_params,
                                                         b_params, axes_row_):
            for xorder_, xorder_str_ in zip(xorders, xorder_strs):
                filename = os.path.join(
                    'bin', '{}_case_{}_{}_xorder_{}_nx1_{}.tab'.format(
                        coord_, case_, integrator, xorder_, nx1_profile))
                data = athena_read.tab(filename)
                x = data['x1v']
                y = data['r0']
                ax.plot(x, y, '{}'.format(xorder_symbols[xorder_]),
                        fillstyle='none', color=xorder_colors[xorder_], label=xorder_str_,
                        markersize=8)
                x_samples = np.linspace(0, 2, nsamples)
                y_samples = EvolvedGaussianProfile(x_samples, a_, b_, m_,  1.0)
                # EvolvedGaussianProfile(x_samples, a_, b_, m_,  0.0)  # Initial condition
                ax.plot(x_samples, y_samples, '-k', linewidth=1.0)

                # Drop "_polar" from "spherical_polar" Athena++ --coord choice for title
                coord_str = coord_.split('_')[0]
                ax.set_title('Radial advection ({})\n{}'.format(coord_str, param_str_))
                ax.set_xlabel(r'$\xi$')
                ax.set_ylabel(r'$Q$')
                # KGF: comment-out next 2x lines to autoscale axes limits
                ax.set_xlim(xlims_)
                ax.set_ylim(ylims_)
                ax.yaxis.set_minor_locator(AutoMinorLocator(4))
                ax.xaxis.set_minor_locator(AutoMinorLocator(4))
                ax.tick_params(direction='in', which='both', axis='both')
                # Hide the right and top spines / plot borders
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)
                # disable rounded edges:
                leg = ax.legend(handles=profile_legend_handles, fancybox=False,)
                leg.get_frame().set_edgecolor('k')

    output_name = 'athena_mignone_fig2'
    pdf_name = "{}.pdf".format(output_name)
    fig.savefig(pdf_name, bbox_inches='tight', dpi=dpi_global)


legend_handles = [
    Line2D([0], [0], marker=xorder_symbols[xorder_], color='w', label=xorder_str_,
           markeredgecolor=xorder_colors[xorder_],
           markerfacecolor='w', fillstyle='none', markersize=8)
    for xorder_, xorder_str_ in zip(xorders, xorder_strs)]

# y-position of first points on each example convergence line:
convergence_y0 = [
    3e-9,
    6e-5,
    5e-12,
    2e-5,
]

n3_xrange = np.array([32., 110.])
n2_xrange = np.array([70., 300.])


def figure3_convergence():
    fig = plt.figure(figsize=(0.75*figsize[0], 0.75*figsize[1]), dpi=dpi_global)
    axes = fig.subplots(2, 2, gridspec_kw={'hspace': 0.4})

    plot_id = 0
    for coord_, axes_row_ in zip(coords, axes):
        for case_, param_str_, ax in zip(cases, case_parameters, axes_row_):
            for xorder_, xorder_str_ in zip(xorders, xorder_strs):
                error_file = os.path.join(
                    'bin', 'errors_{}_case_{}_{}_xorder_{}.dat'.format(coord_, case_,
                                                                       integrator,
                                                                       xorder_))
                # read Athena++ data from error file
                data = athena_read.error_dat(error_file)
                ax.loglog(data[:, 0], data[:, 4], ':{}'.format(xorder_symbols[xorder_]),
                          fillstyle='none', color=xorder_colors[xorder_],
                          label=xorder_str_, markersize=8)

            # Add short example convergence lines in bottom-left corner of plot:
            # N^-3 from nx1=30 to 128:
            y_example = convergence_y0[plot_id]*(n3_xrange/n3_xrange[0])**(-3)
            [ln] = ax.loglog(n3_xrange, y_example, '--k', linewidth=0.75, dashes=(6.5, 5))
            x_text = 0.75*n3_xrange[-1]
            y_text = convergence_y0[plot_id]*(x_text/n3_xrange[0])**(-3)
            ax.annotate(r'$\propto \mathrm{N}^{-3}$', (x_text, y_text))
            #  print(ln.get_linewidth())
            # N^-2 from nx1=64 (or 70?) to around 300?
            y_example = convergence_y0[plot_id]*(n2_xrange/n2_xrange[0])**(-2)
            ax.loglog(n2_xrange, y_example, '--k', linewidth=0.75, dashes=(6.5, 5))
            x_text = 0.4*n2_xrange[-1]
            y_text = convergence_y0[plot_id]*(x_text/n2_xrange[0])**(-2)
            ax.annotate(r'$\propto \mathrm{N}^{-2}$', (x_text, y_text))

            # Drop "_polar" from "spherical_polar" Athena++ --coord choice for title
            coord_str = coord_.split('_')[0]
            ax.set_title('Radial advection ({})\n{}'.format(coord_str, param_str_))
            #  ax.set_xlabel(r'$N$')
            ax.set_xlabel(r'N')
            ax.set_ylabel(r'$\epsilon_1$', usetex=True)
            ax.tick_params(direction='in', which='both', axis='both')
            # Hide the right and top spines / plot borders
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.set_xlim([20, 4000])
            ax.set_ylim(ylims[plot_id])
            ax.set_yticks(major_yticks[plot_id], minor=False)
            if not minor_yticks_auto[plot_id]:
                if minor_yticks[plot_id]:
                    ax.set_yticks(minor_yticks[plot_id], minor=True)
                    ax.set_yticklabels([], minor=True)
                else:
                    ax.tick_params(axis='y', which='minor', left=False)

            ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
            plot_id += 1

            # Remove lines from legend entries, leaving only marker symbols:
            # disable rounded edges
            leg = ax.legend(handles=legend_handles, fancybox=False)
            leg.get_frame().set_edgecolor('k')

    output_name = 'athena_mignone_fig3'
    pdf_name = "{}.pdf".format(output_name)
    fig.savefig(pdf_name, bbox_inches='tight', dpi=dpi_global)
