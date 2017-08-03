import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters

class AmbrosoBetaParameters(AuxQuantity2PhaseParameters):
  def __init__(self):
    AuxQuantity2PhaseParameters.__init__(self)
    self.registerFloatParameter("chi", "Weighting coefficient")

class AmbrosoBeta(AuxQuantity2Phase):
  def __init__(self, params):
    AuxQuantity2Phase.__init__(self, params)
    self.chi = params.get("chi")

  def compute(self, data, der):
    denominator = self.chi * data["arho1"] + (1 - self.chi) * data["arho2"]
    beta = self.chi * data["arho1"] / denominator
    dbeta_darho1 = self.chi / denominator - self.chi * data["arho1"] / denominator / denominator * self.chi
    dbeta_darho2 = - self.chi * data["arho1"] / denominator / denominator * (1 - self.chi)

    data["beta"] = beta
    der["beta"] = {"arho1": dbeta_darho1, "arho2": dbeta_darho2}
