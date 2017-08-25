from Stabilization import Stabilization, StabilizationParameters

class NoStabilizationParameters(StabilizationParameters):
  def __init__(self):
    StabilizationParameters.__init__(self)

## No stabilization
class NoStabilization(Stabilization):
  def __init__(self, params):
    Stabilization.__init__(self, params)

  def needSolutionGradients(self):
    return False

  def createAuxQuantities(self):
    return []

  def createIndependentPhaseKernels(self, phase):
    return []

  def createPhaseInteractionKernels(self):
    return []
