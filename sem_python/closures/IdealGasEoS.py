from numpy import sqrt, log, exp, vectorize

from .EoS import EoS, EoSParameters
from ..utilities.error_utilities import error


def assertNonNegativeSoundSpeedArgSingle(arg, p, v):
    if arg < 0:
        error(
            "Sound speed: negative x in sqrt(x), where x = gamma * p * v:\n" + "p = " + str(p) +
            "\nv = " + str(v))


assertNonNegativeSoundSpeedArg = vectorize(assertNonNegativeSoundSpeedArgSingle)


class IdealGasEoSParameters(EoSParameters):

    def __init__(self):
        EoSParameters.__init__(self)
        self.registerFloatParameter("gamma", "Ratio of specific heats")
        self.registerFloatParameter("R", "Specific gas constant")


class IdealGasEoS(EoS):

    def __init__(self, params):
        EoS.__init__(self)
        self.gamma = params.get("gamma")
        self.R = params.get("R")
        self.cp = self.gamma * self.R / (self.gamma - 1)
        self.cv = self.cp / self.gamma

    def rho(self, p, T):
        rho = p / ((self.gamma - 1) * self.cv * T)
        drho_dp = 1.0 / ((self.gamma - 1) * self.cv * T)
        drho_dT = -p * ((self.gamma - 1) * self.cv * T)**-2 * (self.gamma - 1) * self.cv

        return (rho, drho_dp, drho_dT)

    def rho_from_p_s(self, p, s):
        aux = (s + self.cv * log(p**(self.gamma - 1.0))) / self.cv
        daux_ds = 1.0 / self.cv
        daux_dp = 1.0 / p**(self.gamma - 1.0) * (self.gamma - 1.0) * p**(self.gamma - 2.0)

        T = exp(aux)**(1.0 / self.gamma)
        dT_daux = 1.0 / self.gamma * exp(aux)**(1.0 / self.gamma)
        dT_ds = dT_daux * daux_ds
        dT_dp = dT_daux * daux_dp

        rho, drho_dp_partial, drho_dT = self.rho(p, T)
        drho_dp = drho_dp_partial + drho_dT * dT_dp
        drho_ds = drho_dT * dT_ds

        return (rho, drho_dp, drho_ds)

    def e(self, v, p):
        e_value = 1.0 / (self.gamma - 1) * v * p
        de_dv = 1.0 / (self.gamma - 1) * p
        de_dp = 1.0 / (self.gamma - 1) * v
        return (e_value, de_dv, de_dp)

    def p(self, v, e):
        p_value = (self.gamma - 1) * e / v
        dp_dv = -(self.gamma - 1) * e / v / v
        dp_de = (self.gamma - 1) / v
        return (p_value, dp_dv, dp_de)

    def T(self, v, e):
        T_value = e / self.cv
        dT_dv = 0
        dT_de = 1.0 / self.cv
        return (T_value, dT_dv, dT_de)

    def c(self, v, p):
        # check for sqrt() of negative number
        arg = self.gamma * p * v
        assertNonNegativeSoundSpeedArg(arg, p, v)

        c_value = sqrt(self.gamma * p * v)
        dc_dv = 0.5 / sqrt(self.gamma * p * v) * self.gamma * p
        dc_dp = 0.5 / sqrt(self.gamma * p * v) * self.gamma * v
        return (c_value, dc_dv, dc_dp)

    def s(self, v, e):
        p, dp_dv, dp_de = self.p(v, e)
        T, dT_dv, dT_de = self.T(v, e)

        n = T**self.gamma * p**(1 - self.gamma)
        dn_dT = self.gamma * T**(self.gamma - 1) * p**(1 - self.gamma)
        dn_dp = T**self.gamma * (1 - self.gamma) * p**(-self.gamma)
        dn_dv = dn_dT * dT_dv + dn_dp * dp_dv
        dn_de = dn_dT * dT_de + dn_dp * dp_de

        s = self.cv * log(n)
        ds_dv = self.cv / n * dn_dv
        ds_de = self.cv / n * dn_de

        return (s, ds_dv, ds_de)

    def s_from_h_p(self, h, p):
        aux = p * (h / (self.gamma * self.cv))**(-self.gamma / (self.gamma - 1))
        daux_dh = p * (h / (self.gamma * self.cv))**(-self.gamma/(self.gamma - 1) - 1) \
          * (-self.gamma/(self.gamma - 1)) / (self.gamma * self.cv)
        daux_dp = (h / (self.gamma * self.cv))**(-self.gamma / (self.gamma - 1))

        s = -(self.gamma - 1) * self.cv * log(aux)
        ds_dh = -(self.gamma - 1) * self.cv / aux * daux_dh
        ds_dp = -(self.gamma - 1) * self.cv / aux * daux_dp

        return (s, ds_dh, ds_dp)

    def p_from_h_s(self, h, s):
        p = (h / (self.gamma * self.cv))**(self.gamma / (self.gamma - 1)) \
          * exp(-s / ((self.gamma - 1) * self.cv))
        dp_dh = (self.gamma / (self.gamma - 1)) * (h / (self.gamma * self.cv))**(self.gamma / (self.gamma - 1) - 1) \
          / (self.gamma * self.cv) * exp(-s / ((self.gamma - 1) * self.cv))
        dp_ds = (h / (self.gamma * self.cv))**(self.gamma / (self.gamma - 1)) \
          * exp(-s / ((self.gamma - 1) * self.cv)) / -((self.gamma - 1) * self.cv)

        return (p, dp_dh, dp_ds)

    def h(self, p, T):
        h = self.gamma * self.cv * T
        dh_dp = 0
        dh_dT = self.gamma * self.cv
        return (h, dh_dp, dh_dT)
