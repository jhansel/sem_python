from numpy import sqrt, vectorize, log, exp

from EoS import EoS, EoSParameters
from error_utilities import error

def assertNonNegativeSoundSpeedArgSingle(arg, p, v):
  if arg < 0:
    error("Sound speed: negative x in sqrt(x), where x = gamma * (p + p_inf) * v:\n" +
      "p = " + str(p) + "\nv = " + str(v))

assertNonNegativeSoundSpeedArg = vectorize(assertNonNegativeSoundSpeedArgSingle)

class StiffenedGasEoSParameters(EoSParameters):
  def __init__(self):
    EoSParameters.__init__(self)
    self.registerFloatParameter("gamma", "Ratio of specific heats")
    self.registerFloatParameter("cv", "Specific heat at constant volume")
    self.registerFloatParameter("q", "'q' parameter for stiffened gas")
    self.registerFloatParameter("p_inf", "'p_inf' parameter for stiffened gas")
    self.registerFloatParameter("q_prime", "'q_prime' parameter for stiffened gas")

class StiffenedGasEoS(EoS):
  def __init__(self, params):
    EoS.__init__(self)
    self.gamma = params.get("gamma")
    self.cv = params.get("cv")
    self.q = params.get("q")
    self.p_inf = params.get("p_inf")
    self.q_prime = params.get("q_prime")

  def rho(self, p, T):
    rho = (p + self.p_inf) / ((self.gamma - 1) * self.cv * T)
    drho_dp = 1.0 / ((self.gamma - 1) * self.cv * T)
    drho_dT = -(p + self.p_inf) * ((self.gamma - 1) * self.cv * T)**-2 * (self.gamma - 1) * self.cv

    return (rho, drho_dp, drho_dT)

  def rho_from_p_s(self, p, s):
    aux = (s - self.q_prime + self.cv * log((p + self.p_inf)**(self.gamma - 1.0))) / self.cv
    daux_ds = 1.0 / self.cv
    daux_dp = 1.0 / (p + self.p_inf)**(self.gamma - 1.0) * (self.gamma - 1.0) * (p + self.p_inf)**(self.gamma - 2.0)

    T = exp(aux)**(1.0 / self.gamma)
    dT_daux = 1.0 / self.gamma * exp(aux)**(1.0 / self.gamma)
    dT_ds = dT_daux * daux_ds
    dT_dp = dT_daux * daux_dp

    rho, drho_dp_partial, drho_dT = self.rho(p, T)
    drho_dp = drho_dp_partial + drho_dT * dT_dp
    drho_ds = drho_dT * dT_ds

    return (rho, drho_dp, drho_ds)

  def e(self, v, p):
    e_value = (p + self.gamma * self.p_inf) / (self.gamma - 1) * v + self.q
    de_dv = (p + self.gamma * self.p_inf) / (self.gamma - 1)
    de_dp = v / (self.gamma - 1)
    return (e_value, de_dv, de_dp)

  def p(self, v, e):
    p_value = (self.gamma - 1) * (e - self.q) / v - self.gamma * self.p_inf
    dp_dv = - (self.gamma - 1) * (e - self.q) / v / v
    dp_de = (self.gamma - 1) / v
    return (p_value, dp_dv, dp_de)

  def T(self, v, e):
    T_value = (e - self.q - v * self.p_inf) / self.cv
    dT_dv = - self.p_inf / self.cv
    dT_de = 1.0 / self.cv
    return (T_value, dT_dv, dT_de)

  def c(self, v, p):
    # check for sqrt() of negative number
    arg = self.gamma * (p + self.p_inf) * v
    assertNonNegativeSoundSpeedArg(arg, p, v)

    c_value = sqrt(arg)
    dc_dv = 0.5 / sqrt(self.gamma * (p + self.p_inf) * v) * self.gamma * (p + self.p_inf)
    dc_dp = 0.5 / sqrt(self.gamma * (p + self.p_inf) * v) * self.gamma * v
    return (c_value, dc_dv, dc_dp)

  def s(self, v, e):
    p, dp_dv, dp_de = self.p(v, e)
    T, dT_dv, dT_de = self.T(v, e)

    n = T**self.gamma * (p + self.p_inf)**(1 - self.gamma)
    dn_dT = self.gamma * T**(self.gamma - 1) * (p + self.p_inf)**(1 - self.gamma)
    dn_dp = T**self.gamma * (1 - self.gamma) * (p + self.p_inf)**(-self.gamma)
    dn_dv = dn_dT * dT_dv + dn_dp * dp_dv
    dn_de = dn_dT * dT_de + dn_dp * dp_de

    s = self.cv * log(n) + self.q_prime
    ds_dv = self.cv / n * dn_dv
    ds_de = self.cv / n * dn_de

    return (s, ds_dv, ds_de)

  def p_from_h_s(self, h, s):
    p = ((h - self.q) / (self.gamma * self.cv))**(self.gamma / (self.gamma - 1)) \
      * exp((self.q_prime - s) / ((self.gamma - 1) * self.cv)) - self.p_inf
    dp_dh = (self.gamma / (self.gamma - 1)) * ((h - self.q) / (self.gamma * self.cv))**(self.gamma / (self.gamma - 1) - 1) \
      / (self.gamma * self.cv) * exp((self.q_prime - s) / ((self.gamma - 1) * self.cv))
    dp_ds = ((h - self.q) / (self.gamma * self.cv))**(self.gamma / (self.gamma - 1)) \
      * exp((self.q_prime - s) / ((self.gamma - 1) * self.cv)) / -((self.gamma - 1) * self.cv)

    return (p, dp_dh, dp_ds)
