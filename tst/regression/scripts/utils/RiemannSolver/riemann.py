import numpy as np
from scipy.optimize import brentq
from scipy.integrate import odeint
from scipy.interpolate import interp1d
from scripts.utils.EquationOfState.eos import parse_eos
from . import brent_opt, ode_opt


vector_state_names = ['rho', 'p', 'u']


def sanitize_lbl(label):
    """sanitize_lbl(label): Makes labels latex friendly."""
    import re
    out = label.split('$')
    for i in range(0, len(out), 2):
        out[i] = re.sub('_', r'\_', out[i])
    return '$'.join(out)


class StateVector(object):
    """Fluid state vector."""

    def __init__(self, p=None, rho=None, d=None, T=None, ei=None, u=None, eos=None):
        if d is not None:
            if rho is not None:
                if rho != d:
                    raise ValueError('Conflicting arguments "rho" and "d" specified.')
            rho = d
        self.eos = eos  # equation of state
        self.p = p      # pressure
        self.rho = rho  # density
        self.T = T      # temperature
        self.ei = ei    # internal energy density
        self.u = u      # fluid speed
        self._alt_names = {'press': 'p', 'dens': 'rho', 'vel1': 'u', 'vel': 'u'}

    @property
    def d(self):
        return self.rho

    def ram(self):
        """Computes ram pressure."""
        return self.p + self.rho * self.u ** 2

    def p_min(self):
        """Computes minimum pressure for Riemann solver."""
        out = self.p - self.rho * self.u ** 2
        if out <= 0:
            return .01 * self.p
        return out

    def complete(self, ignore_u=False):
        """Fills in empty parameters if possible."""
        if self.rho is None:
            raise ValueError('Density must be specified.')
        if (not ignore_u) and self.u is None:
            raise ValueError('Speed must be specified.')
        var = ['p', 'ei', 'T']  # need at least one of these to continue
        tmp = {i: getattr(self, i) for i in var if getattr(self, i, None)}
        if not tmp:
            raise ValueError('Insufficient information to complete state.')
        need = [i for i in var if i not in tmp]  # parameters we need to compute
        have = [i for i in var if i in tmp]      # parameters we have
        # is the eos ideal
        try:
            ideal = self.eos.ideal
        except AttributeError:
            ideal = False
        # determine independent variable of the EOS
        if ideal:
            indep = sorted(tmp.keys())[0]
        else:
            indep = self.eos.indep
        # get the independent variable of the EOS if we don't have it
        if indep in need:
            sol = None
            # see if the EOS has a function to get the independent variable
            for i in have:
                try:
                    sol = getattr(self.eos, indep + '_of_rho_' + i)(self.rho, self[i])
                    break
                except AttributeError:
                    pass
            # The EOS doesn't have this function so we need to invert the EOS
            if sol is None:
                eos_f = getattr(self.eos, have[0] + '_of_rho_' + indep)
                val = getattr(self, have[0])

                def f(y):
                    eos_f(self.rho, np.exp(y)) / val - 1.

                # do the root finding
                sol, r = brentq(f, -100, 100, **brent_opt)
                if not r.converged:
                    raise RuntimeError('Unable to converge on ' + indep + '.')
                sol = np.exp(sol)
            setattr(self, indep, sol)
            need.remove(indep)
        # We know the independent variable so use it to get the others
        y = getattr(self, indep)
        for i in need:
            setattr(self, i, getattr(self.eos, i + '_of_rho_' + indep)(self.rho, y))

    def flux(self):
        """Compute fluxes associated with state (mass, momentum, energy)."""
        ram = self.ram()
        return self.rho * self.u, ram + self.p, (.5 * ram + self.p + self.ei) * self.u

    def a(self):
        """Returns adiabatic sound speed."""
        return np.sqrt(self.eos.asq_of_rho_p(self.rho, self.p))

    def es(self):
        """Returns specific internal energy."""
        return self.ei / self.rho

    def __getitem__(self, item):
        """Support retrieving attributes as dictionary items."""
        try:
            out = getattr(self, item)
        except AttributeError:
            try:
                item = self._alt_names[item]
            except KeyError as err:
                try:
                    return self.eos_eval(item)
                except AttributeError:
                    raise err
            return self[item]
        try:
            return out()
        except TypeError:
            return out

    def show(self):
        """Print state variables."""
        msg = 'rho: {0:.3g}, pres: {1:,.3g}, vel: {2:.3g}'.format(
            float(self.rho), self.p, self.u)
        if self.T is not None:
            msg += ', temp: {0:.3g}'.format(self.T)
        print(msg)

    def eos_eval(self, func):
        indep = self[self.eos.indep]
        return getattr(self.eos, func)(self.rho, indep)

    def __repr__(self):
        data = np.array([self.rho, self.p, self.u, self.T], dtype=np.float64)
        data = ', '.join(['%.2g' % i for i in data])
        return '< SVec: ' + data + ' >'


class RiemannSol(object):
    """Specify a Riemann problem and solve it."""

    def __init__(self, left, right, eos):
        left.complete()
        right.complete()
        self.left = left    # left state vector
        self.right = right  # right state vector
        self.eos = eos      # equation of state
        # set some normalizations used in root finding
        self._norm = {'p': min(left.p, right.p)}  # pressure normalisation for root find
        a = max(abs(self.left.u), abs(self.right.u))
        b = .01 * np.sqrt(max(self.left.p / self.left.rho, self.right.p / self.right.rho))
        self._norm['speed'] = 1. / max(a, b)
        # set lower bound for middle pressure
        if left.u <= 0 and right.u >= 0:
            self._pl = min(left.p_min(), right.p_min())
        else:
            self._pl = .5 * min(left.p, right.p)
        self._pu = max(self.left.ram(), self.right.ram()) * \
            2  # upper bound for middle pressure
        self._rare_int_left = RareInt(self, left)    # left rarefaction ODE integrator
        self._rare_int_right = RareInt(self, right)  # right rarefaction ODE integrator
        self.lmid = None
        self.rmid = None
        self.waves = [{}, {}, {}]  # used to keep info on the waves in the problem

    def _shock_jump(self, edge, mid):
        """Find the middle state on the opposite side of a shock from a specified
        edge state and mid.p.

        Inputted mid state must have pressure specified

        Returns mid_state, flux across wave/shock"""
        dp = mid.p - edge.p
        psum = mid.p + edge.p
        if self.eos.ideal:
            # The solution for an ideal EOS is relatively simple
            psg = psum * self.eos.gamma()
            mid.u = np.sqrt(2. / ((dp + psg) * edge.rho)) * dp
            mid.rho = edge.rho * (dp + psg) / (psg - dp)
            mid.ei = self.eos.ei_of_rho_p(mid.rho, mid.p)
            f = np.sqrt(.5 * (dp + psg) * edge.rho)
            return mid, f
        else:
            # for a general EOS we need to root find for mass flux

            def f(fm):
                """Energy function (of mass flux) to find root of."""
                fsqr = fm**2
                rho = fsqr * edge.rho / (fsqr - dp * edge.rho)
                ei = (edge.ei * fsqr + .5 * dp * psum * edge.rho) / (fsqr - dp * edge.rho)
                out = self.eos.ei_of_rho_p(rho, mid.p) / ei - 1.
                return out

            # determine bounding fluxes (fa, fb) for root find
            eu = 2. * self._pu
            if dp > 0:
                fa = (1. + 1e-4) * (dp * edge.rho)**.5
            else:
                fa = 1e-99
            fb = 10. * np.sqrt(.5 * dp * (2. * eu + psum) * edge.rho / edge.ei)
            # find root
            fout, r = brentq(f, fa, fb, **brent_opt)
            if not r.converged:
                raise RuntimeError('Internal energy density not converged.')

            # use mass flux to compute mid state
            fsqr = fout**2
            mid.rho = fsqr * edge.rho / (fsqr - dp * edge.rho)
            mid.ei = self.eos.ei_of_rho_p(mid.rho, mid.p)
            mid.u = dp / fout
            return mid, fout  # return middle state, flux

    def _get_mids(self, p):
        """Determine the two middle states for a specified middle pressure (p*)."""
        lmid = StateVector(p=p, eos=self.eos)  # left middle state
        rmid = StateVector(p=p, eos=self.eos)  # right middle state
        # determine if the left/right waves are shock or rarefaction/simple waves
        if p > self.left.p:  # left wave is shock
            lmid, lf = self._shock_jump(self.left, lmid)
            lmid.u = self.left.u - lmid.u
        else:  # left wave is rarefaction/simple wave
            rho, u = self._rare_int_left(p)
            lmid.rho = rho
            lmid.u = self.left.u - u
            lf = np.nan
        if p > self.right.p:  # right wave is shock
            rmid, rf = self._shock_jump(self.right, rmid)
            rmid.u = self.right.u + rmid.u
        else:  # right wave is rarefaction/simple wave
            rho, u = self._rare_int_right(p)
            rmid.rho = rho
            rmid.u = self.right.u + u
            rf = np.nan
        return lmid, rmid, lf, rf

    def _du(self, pmult):
        """Returns difference in middle state speeds as a function of normalized
        middle pressure."""
        lmid, rmid, _, _ = self._get_mids(pmult * self._norm['p'])
        return (lmid.u - rmid.u) * self._norm['speed']

    def gen_sol(self):
        """Generate the solution to the specified Riemann problem."""
        pnorm = self._norm['p']
        # initialize rarefaction integrators
        self._rare_int_left.init_data(self._pl, self._pu)
        self._rare_int_right.init_data(self._pl, self._pu)
        # find the pressure that gives du=0
        p0, r = brentq(self._du, self._pl / pnorm, self._pu / pnorm, **brent_opt)
        if not r.converged:
            raise RuntimeError('Middle pressure not converged.')

        # save middle states, and keep left/right fluxes
        self.lmid, self.rmid, lf, rf = self._get_mids(p0 * pnorm)
        u = .5 * (self.lmid.u + self.rmid.u)
        self.lmid.u = u
        self.rmid.u = u
        self.lmid.complete()
        self.rmid.complete()
        p = self.lmid.p
        # save wave and mass flux info
        if p > self.left.p:
            tmp = {}
            tmp['kind'] = 'shock'
            tmp['flux'] = lf
            s = self.left.u - tmp['flux'] / self.left.rho
            tmp['speed'] = s, s
            tmp['left'] = self.left
            tmp['right'] = self.lmid
            self.waves[0].update(tmp)
        else:
            tmp = {}
            tmp['kind'] = 'simple'
            tmp['flux'] = None
            tmp['speed'] = self.left.u - self.left.a(), self.lmid.u - self.lmid.a()
            tmp['left'] = self.left
            tmp['right'] = self.lmid
            self.waves[0].update(tmp)
        if p > self.right.p:
            tmp = {}
            tmp['kind'] = 'shock'
            tmp['flux'] = rf
            s = self.right.u + tmp['flux'] / self.right.rho
            tmp['speed'] = s, s
            tmp['left'] = self.rmid
            tmp['right'] = self.right
            self.waves[2].update(tmp)
        else:
            tmp = {}
            tmp['kind'] = 'simple'
            tmp['flux'] = None
            tmp['speed'] = self.rmid.u + self.rmid.a(), self.right.u + self.right.a()
            tmp['left'] = self.rmid
            tmp['right'] = self.right
            self.waves[2].update(tmp)
        tmp = {}
        tmp['kind'] = 'contact'
        tmp['flux'] = 0.
        tmp['speed'] = self.lmid.u, self.lmid.u
        tmp['left'] = self.lmid
        tmp['right'] = self.rmid
        self.waves[1].update(tmp)

    def speeds(self):
        return sum([list(w['speed']) for w in self.waves], [])

    def vector_get_state(self, xi, add_var=None, inc_xi=False):
        """Return the state for a specified sorted vector of characteristic (xi=x/t)."""
        xi = np.atleast_1d(xi)
        names = vector_state_names[:]
        if add_var:
            names += add_var
            names = list(set(names))
        offset = 0
        if inc_xi:
            offset = 1
        data = np.ones((len(names) + offset,) + xi.shape) * np.nan

        indices = np.searchsorted(xi, self.speeds())
        for i, name in enumerate(names):
            data[i, :indices[0]] = self.left[name]
            data[i, indices[1]:indices[2]] = self.lmid[name]
            data[i, indices[3]:indices[4]] = self.rmid[name]
            data[i, indices[5]:] = self.right[name]
        for j in range(indices[0], indices[1]):
            state = self._rare_int_left.characteristic(xi[j])
            for i, name in enumerate(names):
                data[i, j] = state[name]
        for j in range(indices[4], indices[5]):
            state = self._rare_int_right.characteristic(xi[j])
            for i, name in enumerate(names):
                data[i, j] = state[name]
        if inc_xi:
            names.append('xi')
            data[-1] = xi
        return np.core.records.fromarrays(data, names=','.join(names))

    def get_state(self, xi):
        """Return the state for a specified characteristic (xi=x/t)."""
        try:
            return [self.get_state(i) for i in xi]
        except TypeError:
            # use the wave info to determine the part of the Riemann fan xi is in
            if xi < self.waves[0]['speed'][0]:
                return self.left
            if self.waves[0]['speed'][0] < xi < self.waves[0]['speed'][1]:
                return self._rare_int_left.characteristic(xi)
            if self.waves[0]['speed'][1] < xi < self.waves[1]['speed'][0]:
                return self.lmid
            if self.waves[1]['speed'][0] < xi < self.waves[1]['speed'][1]:
                raise ValueError('Should never get here.')
            if self.waves[1]['speed'][1] < xi < self.waves[2]['speed'][0]:
                return self.rmid
            if self.waves[2]['speed'][0] < xi < self.waves[2]['speed'][1]:
                return self._rare_int_right.characteristic(xi)
            if self.waves[2]['speed'][1] < xi:
                return self.right
            # If here we are at a discontinuity or wave edge
            args = ['p', 'u', 'rho']
            if self.waves[0]['speed'][0] == xi:
                if self.waves[0]['kind'] == 'simple':
                    return self.left
                tmp = {i: (self.left[i] if self.left[i]
                           == self.lmid[i] else np.nan) for i in args}
                return StateVector(eos=self.eos, **tmp)
            if self.waves[0]['speed'][1] == xi:
                return self.lmid
            if self.waves[1]['speed'][0] == xi:
                tmp = {i: (self.lmid[i] if self.lmid[i]
                           == self.rmid[i] else np.nan) for i in args}
                return StateVector(eos=self.eos, **tmp)
            if self.waves[2]['speed'][0] == xi:
                if self.waves[2]['kind'] == 'simple':
                    return self.rmid
                tmp = {i: (self.right[i] if self.rmid[i]
                           == self.right[i] else np.nan) for i in args}
                return StateVector(eos=self.eos, **tmp)
            if self.waves[2]['speed'][1] == xi:
                return self.right
            # Should not get here
            return None

    def data_array(self, xi, add_var=None):
        """Returns an array of state data corresponding to the specified array of
        characteristics."""
        xi = np.atleast_1d(xi)
        states = self.get_state(xi)
        names = ['xi', 'dens', 'press', 'vel1']
        if add_var is not None:
            names += add_var
        out = {n: np.array([w[n] for w in states], np.float64) for n in names[1:]}
        out['xi'] = xi
        for i in names[1:]:
            out[i] = np.array([w[i] for w in states])
        out['rho'] = out['dens']
        out['vel'] = out['vel1']
        return out

    def plot_sol(self, var=None, xi_min=None, xi_max=None, nsimp=100, speeds=True,
                 fig=True, ax=None, popt=None, discont=True, lbls=True, t=1, norm=None):
        """Plot solution to Riemann problem.

        Keyword arguments:
        var     -- variable(s) to be plotted (default ['rho', 'p', 'u'])
        xi_min  -- min xi to be plotted (default automatic)
        xi_max  -- max xi to be plotted (default automatic)
        nsimp   -- number of points to use in rarefaction/simple waves (default 100)
        speeds  -- choose whether to plot wave speeds (default True)
        fig     -- if True and ax is None create new figure (default True)
        ax      -- specify ax(es) to use for plotting (default None)
        popt    -- option dictionary to pass to plt.plot (default None)
        discont -- show discontinuities by not connecting lines (default True)
        lbls    -- show axis labels (default True)
        """
        import matplotlib.pyplot as plt

        ws = self.waves

        if popt is None:
            popt = {}
        if var is None:
            var = ['rho', 'p', 'u']
        var = np.atleast_1d(var)
        if norm is None:
            norm = [1 for i in var]
        if fig is True and ax is None:
            fig = plt.figure()

        if var.size > 1:
            n = var.size
            if ax is not None:
                axs = ax
            else:
                axs = [plt.subplot(n, 1, i + 1) for i in range(n)]
            for i, v in enumerate(var):
                ax = axs[i]
                self.plot_sol(v, xi_min=xi_min, xi_max=xi_max, nsimp=nsimp, speeds=speeds,
                              ax=ax, popt=popt, norm=norm[i], lbls=lbls)
            return axs
        else:
            var = var[0]
            try:
                norm = norm[0]
            except TypeError:
                pass
            if ax is None:
                ax = plt.gca()
            plt.sca(ax)

            if xi_min is None:
                xi_min = min(-.1, 1.1 * ws[0]['speed'][0])
            if xi_max is None:
                xi_max = max(.1, 1.1 * ws[2]['speed'][1])
            x = [xi_min, ws[0]['speed'][0]]
            y = [self.left[var]] * 2
            if ws[0]['kind'] == 'simple':
                xs = np.linspace(ws[0]['speed'][0], ws[0]['speed'][1], nsimp + 2)[1:-1]
                ys = [self._rare_int_left.characteristic(i)[var] for i in xs]
                x.extend(xs)
                y.extend(ys)
            elif discont:
                x.append(np.nan)
                y.append(np.nan)
            x.extend([ws[0]['speed'][1], ws[1]['speed'][0]])
            y.extend([self.lmid[var]] * 2)
            if discont:
                x.append(np.nan)
                y.append(np.nan)
            x.extend([ws[1]['speed'][1], ws[2]['speed'][0]])
            y.extend([self.rmid[var]] * 2)
            if ws[2]['kind'] == 'simple':
                xs = np.linspace(ws[2]['speed'][0], ws[2]['speed'][1], nsimp + 2)[1:-1]
                ys = [self._rare_int_right.characteristic(i)[var] for i in xs]
                x.extend(xs)
                y.extend(ys)
            elif discont:
                x.append(np.nan)
                y.append(np.nan)
            x.extend([ws[2]['speed'][1], xi_max])
            y.extend([self.right[var]] * 2)
            x = np.array(x)
            y = np.array(y)
            plt.plot(x * t, y * norm, **popt)
            if lbls:
                plt.ylabel(sanitize_lbl(var))
            plt.xlim(xi_min * t, xi_max * t)
            if speeds:
                waves = [i['speed'][0] for i in ws]
                waves += [i['speed'][1] for i in ws if i['speed'][1] not in waves]
                for i in waves:
                    plt.axvline(i, c='k', ls=':', lw=1)
            return ax

    def fan_plot(self, xlim=None, ylim=None, fig=True):
        """Plot the Riemann fan for this problem."""
        import matplotlib.pyplot as plt
        if xlim is None:
            xlim = [-1, 1]
        if ylim is None:
            ylim = [0, 1]
        xlim = np.atleast_1d(xlim)
        ylim = np.atleast_1d(ylim)
        dx = np.diff(xlim)[0]
        dy = np.diff(ylim)[0]
        y = np.array([0, ylim.max() + dy])
        ws = self.waves

        lws = {'shock': 3, 'contact': 1}
        if fig is True:
            fig = plt.figure()
        for w in ws:
            if w['kind'] == 'simple':
                plt.fill_betweenx(y, y * w['speed'][0], y * w['speed'][1], edgecolor='k',
                                  facecolor='.7', linestyle='-.')
            else:
                plt.plot(y * w['speed'][0], y, 'k-', lw=lws[w['kind']])
        plt.xlim(*xlim)
        plt.ylim(*ylim)
        plt.axis('off')
        plt.axes().set_aspect('equal')
        opt = dict(fontsize=12, ha='center', va='center')
        r = ylim.min() + .8 * dy
        th = .5 * np.pi + .5 * np.arctan2(1, ws[0]['speed'][0])
        plt.text(r * np.cos(th), r * np.sin(th), r'$\mathbf{U}_{L}$', **opt)
        th = .5 * np.arctan2(1, ws[0]['speed'][1]) + .5 * \
            np.arctan(1. / ws[1]['speed'][0])
        plt.text(r * np.cos(th), r * np.sin(th), r'$\mathbf{U}_{L*}$', **opt)
        th = .5 * np.arctan2(1, ws[1]['speed'][1]) + .5 * \
            np.arctan(1. / ws[2]['speed'][0])
        plt.text(r * np.cos(th), r * np.sin(th), r'$\mathbf{U}_{R*}$', **opt)
        th = .5 * np.arctan2(1, ws[2]['speed'][1])
        plt.text(r * np.cos(th), r * np.sin(th), r'$\mathbf{U}_{R}$', **opt)
        opt.pop('ha')
        opt.pop('va')
        aopt = dict(xycoords='data', textcoords='data',
                    arrowprops=dict(arrowstyle="<-", lw=2, color='.3'))
        plt.annotate("", xy=(xlim.min(), 0), xytext=(xlim.max(), 0), **aopt)
        plt.annotate("", xy=(0, ylim.min()), xytext=(0, ylim.max()), **aopt)
        plt.text(0, ylim.max() + .05 * dy, '$t$', ha='center', va='bottom', **opt)
        plt.text(xlim.max() + .025 * dx, 0, '$x$', ha='left', va='center', **opt)

        return fig

    def solve_plot(self, ax=None, popt=None):
        """Solve and plot Riemann problem (see self.gen_sol(), self.plot_sol())."""
        self.gen_sol()
        return self.plot_sol(ax=ax, popt=popt)

    def print_waves(self):
        """Print the wave information."""
        for w in self.waves:
            print(w['kind'], w['flux'], w['speed'])

    @property
    def states(self):
        return [self.left, self.lmid, self.rmid, self.right]

    def speed_row(self, sep=None, fmt='.7e'):
        """Format wave-speeds for printing in a table"""
        fmt = '{0:' + fmt + '}'
        speeds = [fmt.format(i) for i in self.speeds()]
        if speeds[0] == speeds[1]:
            speeds[1] = '-'
        if speeds[4] == speeds[5]:
            speeds[5] = '-'
        speeds.pop(3)
        if sep is None:
            return speeds
        return sep.join(speeds)

    def state_tbl(self, row_sep=None, col_sep=None, eol='', fmt='.7e', speeds=False):
        """Format all for state-vectors for printing in a table"""
        states = self.states
        var_names = ['rho', 'p', 'u']
        if self.eos.indep not in var_names:
            var_names.append(self.eos.indep)
        fmt = '{0:' + fmt + '}'
        out = [[fmt.format(s[v]) for v in var_names] for s in states]
        if speeds:
            speeds = self.speeds()
            out[0].extend([r'$-\infty$', fmt.format(speeds[0])])
            out[1].extend([fmt.format(speeds[1]), fmt.format(speeds[2])])
            out[2].extend([fmt.format(speeds[3]), fmt.format(speeds[4])])
            out[3].extend([fmt.format(speeds[5]), r'$\infty$'])
        if col_sep is not None:
            out = [col_sep.join(col) + eol for col in out]
        if row_sep is not None:
            out = row_sep.join(out)
        return out

    def rare_sol(self):
        out = dict()
        waves = self.waves
        add_var = [self.eos.indep]
        if waves[0]['kind'] == 'simple':
            xi = np.linspace(*waves[0]['speed'], num=101)
            tmp = self.vector_get_state(xi, add_var=add_var, inc_xi=True)
            for i in tmp.dtype.names:
                if i != 'xi':
                    tmp[i][0] = self.left[i]
                    tmp[i][-1] = self.lmid[i]
            out['left'] = tmp
        if waves[2]['kind'] == 'simple':
            xi = np.linspace(*waves[2]['speed'], num=101)
            tmp = self.vector_get_state(xi, add_var=add_var, inc_xi=True)
            for i in tmp.dtype.names:
                if i != 'xi':
                    tmp[i][-1] = self.right[i]
                    tmp[i][-0] = self.rmid[i]
            out['right'] = tmp
        return out

    @property
    def ic(self):
        var_names = ['d', 'u', self.eos.indep]
        out = dict()
        for side in ['left', 'right']:
            state = getattr(self, side)
            out.update({v + side[0]: state[v] for v in var_names})
        return out


class RareInt(object):
    """Rarefaction ODE integrator."""

    def __init__(self, owner, edge, sign=None, slow=False):
        self.owner = owner  # RiemannSol
        self.edge = edge
        self.eos = owner.eos
        try:
            self.ideal = owner.eos.ideal
        except AttributeError:
            self.ideal = False
        self._data = None
        self._rho = None
        self._u = None
        self._call = lambda x: None
        self._pmin = None
        self._pmax = None
        self._pa = None
        self.slow = slow
        if sign is None:
            if self.edge == self.owner.left:
                sign = -1
            if self.edge == self.owner.right:
                sign = 1
        self.sign = sign
        # if the eos is Ideal, there is an analytic solution to the ODE
        if self.ideal:
            def tmp(p):
                g = self.eos._g
                x = .5 * self.eos._gm1 / g
                c = np.sqrt(self.edge.p**(1. / g) / (g * self.edge.rho))
                return edge.rho * (p / edge.p)**(1/g), c * (p**x - self.edge.p**x) / x

            def rho(p):
                g = self.eos._g
                return edge.rho * (p / edge.p)**(1/g)

            def u(p):
                g = self.eos._g
                x = .5 * self.eos._gm1 / g
                c = np.sqrt(self.edge.p**(1. / g) / (g * self.edge.rho))
                return c * (p**x - self.edge.p**x) / x

            self._call = tmp
            self._rho = rho
            self._u = u

    def _rare_ode(self, y, p):
        """Rarefaction wave ODE. y=(density, velocity), p=pressure."""
        am2 = 1. / self.eos.asq_of_rho_p(y[0], p)
        return am2, np.sqrt(am2) / y[0]

    def characteristic(self, xi):
        """Return the state for a given characteristic."""

        # find the pressure corresponding to xi
        def f(p):
            rho = self._rho(p)
            u = self.sign * self._u(p)
            top = self.edge.u + u + self.sign * np.sqrt(self.eos.asq_of_rho_p(rho, p))
            return top / xi - 1.

        p, r = brentq(f, self._pmin, self._pmax, **brent_opt)  # root find
        if not r.converged:
            raise RuntimeError('Pressure within rarefaction wave not converged.')

        # use this pressure to determine the state
        u = self.edge.u + self.sign * self._u(p)
        out = StateVector(rho=self._rho(p), p=p, u=u, eos=self.eos)
        out.complete()
        return out

    def init_data(self, pmin, pmax, bp=12):
        """Initialize ODE integration."""
        if self._data is not None:
            raise RuntimeError('Data already initialized.')
        self._pmin = min(pmin, pmax)
        self._pmax = max(pmin, pmax)
        # if EOS is ideal we initialized in constructor
        if self.ideal:
            self._data = True
            return
        p0 = self.edge.p
        y0 = np.array([self.edge.rho, 0])
        pa = np.linspace(pmin, pmax, 2**bp + 1)
        self._pa = pa
        data = np.zeros((2**bp + 1, 2))
        loc = np.where(pa > p0)
        if loc[0].size > 0:
            tmp = np.concatenate(([p0], pa[loc]))
            tmp = odeint(self._rare_ode, y0, tmp, **ode_opt)
            data[loc] = tmp[1:]
        loc = np.where(pa < p0)
        if loc[0].size > 0:
            tmp = np.concatenate(([p0], pa[loc][::-1]))
            tmp = odeint(self._rare_ode, y0, tmp, **ode_opt)
            data[loc] = tmp[1:][::-1]
        self._data = data
        if self.slow:
            def f(p):
                idx = (np.abs(self._pa - p)).argmin()
                out = odeint(self._rare_ode, self._data[idx], np.array(
                    [self._pa[idx], p]), **ode_opt)[-1]
                return out
            self._call = f
            self._rho = lambda p: self._call(p)[0]
            self._u = lambda p: self._call(p)[1]
        else:
            self._rho = interp1d(pa, data[:, 0], kind='cubic')
            self._u = interp1d(pa, data[:, 1], kind='cubic')
            self._call = lambda p: (self._rho(p), self._u(p))

    def __call__(self, p):
        return self._call(p)


def riemann_problem(states, eos):
    """A RiemannSol wrapper to specify left/right state variables via a dict.
    Each key in state has form of [dpTu][lr].
    Example: state = dict(dl=1, pl=1, ul=0, dr=0.125, pr=1, ur=0)
    """
    eos = parse_eos(eos)
    left = {'eos': eos}
    right = {'eos': eos}
    for i in states:
        var = i[:-1]
        side = i[-1]
        if side == 'l':
            left[var] = float(states[i])
        elif side == 'r':
            right[var] = float(states[i])
        else:
            raise ValueError('Parameter "{0:}" not understood.'.format(i))
    left = StateVector(**left)
    right = StateVector(**right)
    out = RiemannSol(left, right, eos)
    out.gen_sol()
    return out
