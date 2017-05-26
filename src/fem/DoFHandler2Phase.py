import numpy as np

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/base")
from enums import ModelType, VariableName

sys.path.append(base_dir + "src/fem")
from DoFHandler import DoFHandler

class DoFHandler2Phase(DoFHandler):
  def __init__(self, mesh):
    DoFHandler.__init__(self, mesh)
    self.model_type = ModelType.TwoPhase
    self.setup()

  def getVolumeFraction(self, U, k):
    return U[self.i(k, self.vf1_index[0])]
