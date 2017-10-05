from numpy import vectorize

from .AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters
from ..utilities.error_utilities import error

def assertNonNegativePressureSingle(p, v, e):
  if p < 0:
    error("Negative pressure from p(v,e):\n" +
      "  p = " + str(p) + "\n" +
      "  v = " + str(v) + "\n" +
      "  e = " + str(e))

assertNonNegativePressure = vectorize(assertNonNegativePressureSingle)

class PressureParameters(AuxQuantity1PhaseParameters):
  def __init__(self):
    AuxQuantity1PhaseParameters.__init__(self)
    self.registerFunctionParameter("p_function", "Pressure function from EOS")

class Pressure(AuxQuantity1Phase):
  def __init__(self, params):
    AuxQuantity1Phase.__init__(self, params)
    self.name = self.p
    self.p_function = params.get("p_function")

  def compute(self, data, der):
    p, dp_dv, dp_de = self.p_function(data[self.v], data[self.e])
    data[self.name] = p

    assertNonNegativePressure(p, data[self.v], data[self.e])

    dp_daA1 = dp_dv * der[self.v]["aA1"]
    dp_darhoA = dp_dv * der[self.v][self.arhoA] + dp_de * der[self.e][self.arhoA]
    dp_darhouA = dp_de * der[self.e][self.arhouA]
    dp_darhoEA = dp_de * der[self.e][self.arhoEA]
    der[self.name]["aA1"] = dp_daA1
    der[self.name][self.arhoA] = dp_darhoA
    der[self.name][self.arhouA] = dp_darhouA
    der[self.name][self.arhoEA] = dp_darhoEA
