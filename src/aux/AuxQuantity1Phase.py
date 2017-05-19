import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from AuxQuantity import AuxQuantity, AuxQuantityParameters

class AuxQuantity1PhaseParameters(AuxQuantityParameters):
  def __init__(self):
    AuxQuantityParameters.__init__(self)
    self.registerFunctionParameter("p_function", "Pressure function from EOS")

class AuxQuantity1Phase(AuxQuantity):
  def __init__(self, params):
    AuxQuantity.__init__(self)
    self.p_function = params.get("p_function")
