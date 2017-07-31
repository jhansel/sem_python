from numpy import sqrt

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/closures")
from EoS import EoS, EoSParameters

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
    return p / (self.gamma - 1) / self.cv / T

  def e(self, v, p):
    e_value = 1.0 / (self.gamma - 1) * v * p
    de_dv = 1.0 / (self.gamma - 1) * p
    de_dp = 1.0 / (self.gamma - 1) * v
    return (e_value, de_dv, de_dp)

  def p(self, v, e):
    p_value = (self.gamma - 1) * e / v
    dp_dv = - (self.gamma - 1) * e / v / v
    dp_de = (self.gamma - 1) / v
    return (p_value, dp_dv, dp_de)

  def T(self, v, e):
    T_value = e / self.cv
    dT_dv = 0
    dT_de = 1.0 / self.cv
    return (T_value, dT_dv, dT_de)

  def c(self, v, p):
    c_value = sqrt(self.gamma * p * v)
    dc_dv = 0.5 / sqrt(self.gamma * p * v) * self.gamma * p
    dc_dp = 0.5 / sqrt(self.gamma * p * v) * self.gamma * v
    return (c_value, dc_dv, dc_dp)
