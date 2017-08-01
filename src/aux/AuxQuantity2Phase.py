import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from AuxQuantity import AuxQuantity, AuxQuantityParameters

class AuxQuantity2PhaseParameters(AuxQuantityParameters):
  def __init__(self):
    AuxQuantityParameters.__init__(self)

class AuxQuantity2Phase(AuxQuantity):
  def __init__(self, params):
    AuxQuantity.__init__(self, params)
