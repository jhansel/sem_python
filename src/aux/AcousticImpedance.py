from AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

class AcousticImpedanceParameters(AuxQuantity1PhaseParameters):
  def __init__(self):
    AuxQuantity1PhaseParameters.__init__(self)

class AcousticImpedance(AuxQuantity1Phase):
  def __init__(self, params):
    AuxQuantity1Phase.__init__(self, params)
    self.name = self.z

  def compute(self, data, der):
    data[self.name] = data[self.rho] * data[self.c]

    drhoc_dvf1 = der[self.rho]["vf1"] * data[self.c] + data[self.rho] * der[self.c]["vf1"]
    drhoc_darhoA = der[self.rho][self.arhoA] * data[self.c] + data[self.rho] * der[self.c][self.arhoA]
    drhoc_darhouA = data[self.rho] * der[self.c][self.arhouA]
    drhoc_darhoEA = data[self.rho] * der[self.c][self.arhoEA]

    der[self.name]["vf1"] = drhoc_dvf1
    der[self.name][self.arhoA] = drhoc_darhoA
    der[self.name][self.arhouA] = drhoc_darhouA
    der[self.name][self.arhoEA] = drhoc_darhoEA
