import numpy as np

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/base")
from enums import ModelType

sys.path.append(base_dir + "src/fem")
from DoFHandler import DoFHandler

class DoFHandler2PhaseNonInteracting(DoFHandler):
  def __init__(self, mesh, ics):
    DoFHandler.__init__(self, mesh)
    self.model_type = ModelType.TwoPhaseNonInteracting
    self.setup()

    # create array for volume fraction
    self.vf1 = np.zeros(self.n_node)
    for k in xrange(self.n_node):
      self.vf1[k] = ics.vf1(mesh.x[k])

  def getVolumeFraction(self, U, k):
    return self.vf1[k]
