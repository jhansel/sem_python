import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters

class BerryInterfacialAreaDensityParameters(AuxQuantity2PhaseParameters):
  def __init__(self):
    AuxQuantity2PhaseParameters.__init__(self)
    self.registerFloatParameter("a_int_min", "Minimum interfacial area density", 0)
    self.registerFloatParameter("a_int_max", "Maximum interfacial area density")

class BerryInterfacialAreaDensity(AuxQuantity2Phase):
  def __init__(self, params):
    AuxQuantity2Phase.__init__(self, params)
    self.a_int_min = params.get("a_int_min")
    self.a_int_max = params.get("a_int_max")
    self.coef_A = 4 * (self.a_int_min - self.a_int_max)
    self.coef_B = - self.coef_A
    self.coef_C = self.a_int_min

  def compute(self, data, der):
    vf1 = data["vf1"]
    data["a_int"] = self.coef_A * vf1**2 + self.coef_B * vf1 + self.coef_C

    da_int_dvf1 = 2 * self.coef_A * vf1 + self.coef_B

    der["a_int"] = {"vf1": da_int_dvf1}
