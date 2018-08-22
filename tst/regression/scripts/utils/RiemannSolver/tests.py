#! /usr/bin/env python

import numpy as np
from .riemann import RiemannSol, StateVector
from .eos import Ideal, SimpleHydrogen


def sod_test(eos=None, gamma=None, plot=True):
    if eos is None:
        if gamma is None:
            gamma = 5. / 3.
        from .eos import Ideal
        eos = Ideal(gamma)

    l = StateVector(p=1, rho=1, u=0, eos=eos)
    l.complete()
    r = StateVector(p=.1, rho=.125, u=0, eos=eos)
    r.complete()
    rs = RiemannSol(l, r, eos)
    rs.gen_sol()

    if plot:
        rs.plot_sol()
        import matplotlib.pyplot as plt
        plt.show()

        rs.fan_plot()
        plt.show()

    return rs


def eos_test(plot=True):
    r0 = 1e-7
    p0 = 3e-8
    l = dict(rho=r0, p=p0, u=0.)
    r = dict(rho=.125 * r0, p=1e-9, u=0.)
    eoss = [SimpleHydrogen()]
    sols = [RiemannSol(StateVector(eos=eoss[-1], **l), StateVector(eos=eoss[-1], **r), eoss[-1])]
    sols[0].gen_sol()

    l['T'] = None
    l['p'] = sols[0].left.p
    r['T'] = None
    r['p'] = sols[0].right.p

    gs = np.array([eoss[0]._gamma1(w.rho, w.T) for w in sols[0].states])
    eoss.append(Ideal(gs.max()))
    sols.append(RiemannSol(StateVector(eos=eoss[-1], **l), StateVector(eos=eoss[-1], **r), eoss[-1]))
    eoss.append(Ideal(gs.min()))
    sols.append(RiemannSol(StateVector(eos=eoss[-1], **l), StateVector(eos=eoss[-1], **r), eoss[-1]))
    sols[1].gen_sol()
    sols[2].gen_sol()

    if plot:
        x0 = np.min([i.waves[0]['speed'][0] for i in sols])
        x1 = np.max([i.waves[2]['speed'][1] for i in sols])
        dx = x1 - x0
        x0 -= .1 * dx
        x1 += .1 * dx
        ax = None
        for sol in sols:
            ax = sol.plot_sol(ax=ax, speeds=False, xi_min=x0, xi_max=x1)

        lbls = ['H'] + [r'$\gamma={0:.3f}$'.format(i._g) for i in eoss[1:]]
        import matplotlib.pyplot as plt
        plt.legend(lbls)
        plt.show()

        sols[0].fan_plot()
        plt.show()

    return sols