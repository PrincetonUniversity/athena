import numpy as np
from scipy.optimize import brentq
import sys
from . import brent_opt


class EOS(object):
    def __init__(self):
        self.ideal = False
        self.indep = None
        pass

    def valid(self):
        try:
            self.asq_of_rho_p(1, 1)
            self.ei_of_rho_p(1, 1)
        except NotImplementedError:
            return False
        return True

    def asq_of_rho_p(self, rho, p):
        raise NotImplementedError

    def ei_of_rho_p(self, rho, p):
        raise NotImplementedError

    def es_of_rho_p(self, rho, p):
        return self.ei_of_rho_p(rho, p) / rho

    def p_of_rho_es(self, rho, p):
        raise NotImplementedError


class SimpleHydrogen(EOS):
    def __init__(self):
        super(SimpleHydrogen, self).__init__()
        self.indep = 'T'

    def _phi(self, T):
        return np.exp(1. / T - 1.5 * np.log(T))

    def _x(self, rho, T):
        return 2. / (1 + np.sqrt(1 + 4. * np.exp(1. / T - 1.5 * np.log(T) + np.log(rho))))

    def _x_T(self, rho, T):
        x = self._x(rho, T)
        return x**3 / (2. + x) * np.exp(1. / T - 3.5 * np.log(T)) * (1. + 1.5 * T) * rho

    def p_of_rho_T(self, rho, T):
        return rho * T * (1. + self._x(rho, T))

    def ei_of_rho_T(self, rho, T):
        return self._x(rho, T) * rho + 1.5 * self.p_of_rho_T(rho, T)

    def _b(self, rho, T):
        lt = np.log(T)
        c1 = np.exp(-1.25 * lt - .5 / T)
        c2 = np.exp(1.5 * lt - 1. / T)
        return 8. * rho * c1 / (np.sqrt(c2) + np.sqrt(c2 + 4. * rho))**3

    def gamma1(self, rho, T):
        x = self._x(rho, T)
        b = self._b(rho, T)
        return (b * (4. + 20. * T + 15. * T**2) + 10. * (2. + x - x**2)) /\
               (b * (2. + 3. * T)**2 + 6.*(2. + x - x**2))

    def asq_of_rho_T(self, rho, T):
        return T * (1. + self._x(rho, T)) * self.gamma1(rho, T)

    def asq_of_rho_p(self, rho, p):
        return p * self.gamma1(rho, self.T_of_rho_p(rho, p)) / rho

    def T_of_rho_p(self, rho, p):
        t1 = p / rho * (1. + sys.float_info.epsilon)

        def f(y):
            return self.p_of_rho_T(rho, y) / p - 1.

        T, r = brentq(f, .1 * t1, t1, **brent_opt)
        if not r.converged:
            raise RuntimeError('Unable to converge on temperature.')
        return T

    def T_of_rho_ei(self, rho, ei):
        t1 = ei / rho * (1. + sys.float_info.epsilon)

        def f(y):
            return self.ei_of_rho_T(rho, y) / ei - 1.

        T, r = brentq(f, .1 * t1, 10 * t1, **brent_opt)
        if not r.converged:
            raise RuntimeError('Unable to converge on temperature.')
        return T

    def ei_of_rho_p(self, rho, p):
        return self.ei_of_rho_T(rho, self.T_of_rho_p(rho, p))

    def p_of_rho_es(self, rho, es):
        return self.p_of_rho_T(rho, self.T_of_rho_ei(rho, rho * es))


class Ideal(EOS):
    def __init__(self, gamma, R=1):
        if gamma <= 1:
            raise ValueError('The value for gamma must be larger than 1.')
        super(Ideal, self).__init__()
        self.ideal = True
        self._g = gamma
        self._gm1 = gamma - 1.
        self.R = R

    def gamma(self):
        return self._g

    def asq_of_rho_p(self, rho, p):
        return self._g * p / rho

    def ei_of_rho_p(self, rho, p):
        return p / self._gm1

    def T_of_rho_p(self, rho, p):
        return p / (rho * self.R)

    def T_of_rho_ei(self, rho, ei):
        return ei * self._gm1 / (rho * self.R)

    def p_of_rho_ei(self, rho, ei):
        return ei * self._gm1

    def p_of_rho_es(self, rho, es):
        return rho * es * self._gm1


class TestIdeal(Ideal):
    """Class to test if Riemann solver gives same answer as Ideal."""

    def __init__(self, gamma, R=1):
        super(TestIdeal, self).__init__(gamma, R=R)
        self.ideal = False
        self.indep = 'p'


def parse_eos(eos):
    if hasattr(eos, 'asq_of_rho_p'):
        return eos
    if eos == 'H' or eos == 'h':
        return SimpleHydrogen()
    try:
        return Ideal(float(eos))
    except ValueError:
        raise ValueError('Cannot parse EOS "{0:}".'.format(eos))
