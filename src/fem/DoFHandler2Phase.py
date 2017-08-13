import numpy as np

from enums import ModelType, VariableName
from DoFHandler import DoFHandler, DoFHandlerParameters

class DoFHandler2PhaseParameters(DoFHandlerParameters):
  def __init__(self):
    DoFHandlerParameters.__init__(self)

class DoFHandler2Phase(DoFHandler):
  def __init__(self, params):
    DoFHandler.__init__(self, params)
    self.model_type = ModelType.TwoPhase
    self.setup()

  def getVolumeFraction(self, U, k):
    return U[self.i(k, self.vf1_index[0])]
