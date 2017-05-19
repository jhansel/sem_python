import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

class VolumeFractionPhase1Parameters(AuxQuantity1PhaseParameters):
  def __init__(self):
    AuxQuantity1PhaseParameters.__init__(self)

class VolumeFractionPhase1(AuxQuantity1Phase):
  def __init__(self, params):
    AuxQuantity1Phase.__init__(self, params)

  def compute(self, data, der):
    data["vf"] = data["vf1"]

    der["vf"] = {"vf1" : 0 * data["vf1"] + 1.0}
