import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from AuxQuantity2Phase import AuxQuantity2Phase, AuxQuantity2PhaseParameters

class BetaParameters(AuxQuantity2PhaseParameters):
  def __init__(self):
    AuxQuantity2PhaseParameters.__init__(self)
    self.registerParameter("beta_function", "Function for computing beta")

class Beta(AuxQuantity2Phase):
  def __init__(self, params):
    AuxQuantity2Phase.__init__(self, params)
    self.beta_function = params.get("beta_function")

  def compute(self, data, der):
    data["beta"], dbeta_darho1, dbeta_darho2 = self.beta_function(data["arho1"], data["arho2"])
    der["beta"] = {"arho1": dbeta_darho1, "arho2": dbeta_darho2}
