import numpy as np

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/base")
from enums import ModelType

sys.path.append(base_dir + "src/fem")
from DoFHandler import DoFHandler, DoFHandlerParameters

class DoFHandler2PhaseNonInteractingParameters(DoFHandlerParameters):
  def __init__(self):
    DoFHandlerParameters.__init__(self)
    self.registerFunctionParameter("initial_vf1", "Initial phase-1 volume fraction function")

class DoFHandler2PhaseNonInteracting(DoFHandler):
  def __init__(self, params):
    DoFHandler.__init__(self, params)
    initial_vf1 = params.get("initial_vf1")
    self.model_type = ModelType.TwoPhaseNonInteracting
    self.setup()

    # create array for volume fraction
    self.vf1 = np.zeros(self.n_node)
    for k in xrange(self.n_node):
      self.vf1[k] = initial_vf1(self.mesh.x[k])

  def getVolumeFraction(self, U, k):
    return self.vf1[k]
