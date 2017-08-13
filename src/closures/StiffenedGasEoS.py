from numpy import sqrt, vectorize

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
    self.registerFloatParameter("q", "")
    self.registerFloatParameter("p_inf", "")

class StiffenedGasEoS(EoS):
  def __init__(self, params):
    EoS.__init__(self)
    self.gamma = params.get("gamma")
    self.cv = params.get("cv")
    self.q = params.get("q")
    self.p_inf = params.get("p_inf")

  def rho(self, p, T):
    return (p + self.p_inf) / ((self.gamma - 1) * self.cv * T)

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
