from AuxQuantity import AuxQuantity, AuxQuantityParameters

class AuxQuantity1PhaseParameters(AuxQuantityParameters):
  def __init__(self):
    AuxQuantityParameters.__init__(self)
    self.registerIntParameter("phase", "Phase index (0 or 1)")

class AuxQuantity1Phase(AuxQuantity):
  def __init__(self, params):
    AuxQuantity.__init__(self, params)
    self.phase_int = params.get("phase")
    self.phase = str(self.phase_int + 1)
    if self.phase_int == 0:
      self.sign = 1.0
      self.dgrad_aA_dgrad_aA1 = 1.0
    else:
      self.sign = -1.0
      self.dgrad_aA_dgrad_aA1 = -1.0

    self.aA = "aA" + self.phase
    self.arhoA = "arhoA" + self.phase
    self.arhouA = "arhouA" + self.phase
    self.arhoEA = "arhoEA" + self.phase
    self.vf = "vf" + self.phase
    self.rho = "rho" + self.phase
    self.rhoe = "rhoe" + self.phase
    self.u = "u" + self.phase
    self.E = "E" + self.phase
    self.e = "e" + self.phase
    self.v = "v" + self.phase
    self.p = "p" + self.phase
    self.T = "T" + self.phase
    self.c = "c" + self.phase
    self.z = "z" + self.phase

    self.grad_aA = "grad_" + self.aA
    self.grad_arhoA = "grad_" + self.arhoA
    self.grad_arhouA = "grad_" + self.arhouA
    self.grad_arhoEA = "grad_" + self.arhoEA
    self.grad_vf = "grad_" + self.vf
    self.grad_rho = "grad_" + self.rho
    self.grad_u = "grad_" + self.u
    self.grad_rhoe = "grad_" + self.rhoe

    self.visccoef_arhoA = "visccoef_" + self.arhoA
    self.visccoef_arhouA = "visccoef_" + self.arhouA
    self.visccoef_arhoEA = "visccoef_" + self.arhoEA

    self.viscflux_aA = "viscflux_" + self.aA
    self.viscflux_arhoA = "viscflux_" + self.arhoA
    self.viscflux_arhouA = "viscflux_" + self.arhouA
    self.viscflux_arhoEA = "viscflux_" + self.arhoEA
