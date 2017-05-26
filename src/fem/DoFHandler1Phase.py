import numpy as np

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/base")
from enums import ModelType

sys.path.append(base_dir + "src/fem")
from DoFHandler import DoFHandler

class DoFHandler1Phase(DoFHandler):
  def __init__(self, mesh):
    DoFHandler.__init__(self, mesh)
    self.model_type = ModelType.OnePhase
    self.setup()

  def getVolumeFraction(self, U, k):
    return 1
