import numpy as np
from scipy.optimize import brentq
from scipy.interpolate import RectBivariateSpline as RBS
import sys
from . import brent_opt


class EOS(object):
    """Parent class to implement equation of state functions"""
    def __init__(self):
        """Initialize EOS class"""
        self.ideal = False  # Whether this EOS is ideal
        self.indep = None   # independent variables not including density

    def valid(self):
        """Determine if this EOS is valid"""
        try:
            self.asq_of_rho_p(1, 1)
            self.ei_of_rho_p(1, 1)
        except NotImplementedError:
            return False
        return True

    def asq_of_rho_p(self, rho, p):
        """Adiabatic sound speed^2 as a function of density (rho) and pressure (p)"""
        raise NotImplementedError

    def ei_of_rho_p(self, rho, p):
        """Internal energy density as a function of density (rho) and pressure (p)"""
        raise NotImplementedError

    def es_of_rho_p(self, rho, p):
        """Specific internal energy as a function of density (rho) and pressure (p)"""
        return self.ei_of_rho_p(rho, p) / rho

    def p_of_rho_es(self, rho, es):
        """Pressure as a function of density (rho) and specific internal energy (es)"""
        raise NotImplementedError


class SimpleHydrogen(EOS):
    """Simple hydrogen equation of state"""
    def __init__(self):
        super(SimpleHydrogen, self).__init__()
        self.indep = 'T'  # Temperature is the independent variable other than density
        self.T_of_rho_ei = np.vectorize(self._T_of_rho_ei)
        self.T_of_rho_p = np.vectorize(self._T_of_rho_p)
        self.T_of_rho_h = np.vectorize(self._T_of_rho_h)

    def _phi(self, T):
        return np.exp(1. / T - 1.5 * np.log(T))

    def _x(self, rho, T):
        """Ionization fraction"""
        with np.errstate(over='ignore'):
            return 2. / (1 + np.sqrt(1 + 4.
                                     * np.exp(1. / T - 1.5 * np.log(T) + np.log(rho))))

    def _x_T(self, rho, T):
        """Temperature derivative of ionization fraction"""
        x = self._x(rho, T)
        return x**3 / (2. + x) * np.exp(1. / T - 3.5 * np.log(T)) * (1. + 1.5 * T) * rho

    def p_of_rho_T(self, rho, T):
        """Pressure as a function of density (rho) and temperature (T)"""
        return rho * T * (1. + self._x(rho, T))

    def ei_of_rho_T(self, rho, T):
        """Internal energy density as a function of density (rho) and temperature (T)"""
        return self._x(rho, T) * rho + 1.5 * self.p_of_rho_T(rho, T)

    def _b(self, rho, T):
        lt = np.log(T)
        c1 = np.exp(-1.25 * lt - .5 / T)
        c2 = np.exp(1.5 * lt - 1. / T)
        return 8. * rho * c1 / (np.sqrt(c2) + np.sqrt(c2 + 4. * rho))**3

    def gamma1(self, rho, T):
        """Gamma_1 as a function of density (rho) and temperature (T)"""
        x = self._x(rho, T)
        b = self._b(rho, T)
        return (b * (4. + 20. * T + 15. * T**2) + 10. * (2. + x - x**2)) /\
               (b * (2. + 3. * T)**2 + 6.*(2. + x - x**2))

    def asq_of_rho_T(self, rho, T):
        """Adiabatic sound speed^2 as a function of density (rho) and temperature (T)"""
        return T * (1. + self._x(rho, T)) * self.gamma1(rho, T)

    def asq_of_rho_p(self, rho, p):
        """Adiabatic sound speed^2 as a function of density (rho) and pressure (p)"""
        return p * self.gamma1(rho, self.T_of_rho_p(rho, p)) / rho

    def asq_of_rho_h(self, rho, h):
        """Adiabatic sound speed^2 function of density (rho) and specific enthalpy (h)"""
        return self.asq_of_rho_T(rho, self.T_of_rho_h(rho, h))

    def _T_of_rho_h(self, rho, h):
        """Temperature as a function of density (rho) and specific enthalpy (h)"""
        t1 = .4 * h * (1. + sys.float_info.epsilon)

        def f(y):
            return (self.p_of_rho_T(rho, y) + self.ei_of_rho_T(rho, y)) / (h * rho) - 1.

        T, r = brentq(f, .1 * t1, t1, **brent_opt)
        if not r.converged:
            raise RuntimeError('Unable to converge on temperature.')
        return T

    def _T_of_rho_p(self, rho, p):
        """Temperature as a function of density (rho) and pressure (p)"""
        t1 = p / rho * (1. + sys.float_info.epsilon)  # initial guess

        def f(y):  # function to find root of
            return self.p_of_rho_T(rho, y) / p - 1.

        try:
            T, r = brentq(f, .1 * t1, t1, **brent_opt)  # find root
        except ValueError:
            T, r = brentq(f, .05 * t1, 2 * t1, **brent_opt)
        if not r.converged:
            raise RuntimeError('Unable to converge on temperature.')
        return T

    def _T_of_rho_ei(self, rho, ei):
        """Temperature as a function of density (rho) and internal energy density (e)"""
        t1 = ei / rho * (1. + sys.float_info.epsilon)  # initial guess

        def f(y):   # function to find root of
            return self.ei_of_rho_T(rho, y) / ei - 1.

        T, r = brentq(f, .05 * t1, 2 * t1, **brent_opt)
        if not r.converged:
            raise RuntimeError('Unable to converge on temperature.')
        return T

    def ei_of_rho_p(self, rho, p):
        """Internal energy density as a function of density (rho) and pressure (p)"""
        return self.ei_of_rho_T(rho, self.T_of_rho_p(rho, p))

    def p_of_rho_es(self, rho, es):
        """Pressure as a function of density (rho) and specific internal energy (es)"""
        return self.p_of_rho_T(rho, self.T_of_rho_ei(rho, rho * es))


class Ideal(EOS):
    """Ideal equation of state class"""
    def __init__(self, gamma, R=1):
        """Adiabatic index "gamma" (>1) and ideal gas constant "R"."""
        if gamma <= 1:
            raise ValueError('The value for gamma must be larger than 1.')
        super(Ideal, self).__init__()
        self.ideal = True
        self._g = gamma
        self._gm1 = gamma - 1.
        self.R = R

    def gamma(self):
        """Returns gamma"""
        return self._g

    def asq_of_rho_p(self, rho, p):
        """Adiabatic sound speed^2 as a function of density (rho) and pressure (p)"""
        return self._g * p / rho

    def ei_of_rho_p(self, rho, p):
        """Internal energy density as a function of density (rho) and pressure (p)"""
        return p / self._gm1

    def T_of_rho_p(self, rho, p):
        """Temperature as a function of density (rho) and pressure (p)"""
        return p / (rho * self.R)

    def T_of_rho_ei(self, rho, ei):
        """Temperature as a function of density (rho) and internal energy density (ei)"""
        return ei * self._gm1 / (rho * self.R)

    def p_of_rho_ei(self, rho, ei):
        """Pressure as a function of density (rho) and internal energy density (ei)"""
        return ei * self._gm1

    def p_of_rho_es(self, rho, es):
        """Pressure as a function of density (rho) and specific internal energy (es)"""
        return rho * es * self._gm1

    def p_of_rho_T(self, rho, T):
        """Pressure as a function of density (rho) and temperature (T)"""
        return rho * T * self.R


class TestIdeal(Ideal):
    """Class to test if Riemann solver gives same answer as Ideal."""

    def __init__(self, gamma, R=1):
        super(TestIdeal, self).__init__(gamma, R=R)
        self.ideal = False
        self.indep = 'p'


class AthenaTable(EOS):
    def __init__(self, data, lrho, le, ratios=None, indep=None, dens_pow=-1, fn=None,
                 add_var=None):
        super(EOS, self).__init__()
        self.fn = fn
        if ratios is None:
            ratios = np.ones(data.shape[0])
        lr = np.log(ratios)
        self._lr = lr
        if indep is None:
            indep = 'ei'
        self.indep = indep
        self.data = data
        self.lrho = lrho
        self.le = le
        self.dens_pow = dens_pow
        var = ['p', 'e', 'asq_p']
        if add_var is not None:
            var.extend(add_var)
        d = {var[i]: RBS(lrho, le + lr[i], np.log10(data[i].T), kx=1, ky=1).ev
             for i in range(len(var))}
        self._interp_dict = d

    def _interp(self, rho, e, var):
        ld = np.log10(rho)
        return 10**self._interp_dict[var](ld, np.log10(e) + self.dens_pow * ld)

    def asq_of_rho_p(self, rho, p):
        """Adiabatic sound speed^2 as a function of density (rho) and pressure (p)"""
        return self._interp(rho, p, 'asq_p') * p / rho

    def ei_of_rho_p(self, rho, p):
        """Internal energy density as a function of density (rho) and pressure (p)"""
        return self._interp(rho, p, 'e') * p

    def es_of_rho_p(self, rho, p):
        """Specific internal energy as a function of density (rho) and pressure (p)"""
        return self._interp(rho, p, 'e') * p / rho

    def p_of_rho_ei(self, rho, ei):
        """Pressure as a function of density (rho) and internal energy density (ei)"""
        return self._interp(rho, ei, 'p') * ei

    def p_of_rho_es(self, rho, es):
        """Pressure as a function of density (rho) and specific internal energy (es)"""
        return self.p_of_rho_ei(rho, es / rho)


def parse_eos(eos):
    """Function to interpret input as an EOS"""
    if hasattr(eos, 'asq_of_rho_p'):
        return eos  # already is EOS class
    if eos == 'H' or eos == 'h':
        return SimpleHydrogen()
    try:
        return Ideal(float(eos))  # try parsing as a gamma value
    except ValueError:
        raise ValueError('Cannot parse EOS "{0:}".'.format(eos))
