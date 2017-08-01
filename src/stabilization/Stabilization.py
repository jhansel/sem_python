from abc import ABCMeta, abstractmethod

import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/input")
from Parameters import Parameters

class StabilizationParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)
    self.registerParameter("factory", "Factory")
    self.registerParameter("dof_handler", "Degree of freedom handler")

## Abstract base class for stabilization
class Stabilization(object):
  __metaclass__ = ABCMeta

  def __init__(self, params):
    self.factory = params.get("factory")
    self.dof_handler = params.get("dof_handler")

  @abstractmethod
  def needSolutionGradients(self):
    pass

  @abstractmethod
  def createIndependentPhaseAuxQuantities(self, phase):
    pass

  @abstractmethod
  def createPhaseInteractionAuxQuantities(self):
    pass

  @abstractmethod
  def createIndependentPhaseKernels(self, phase):
    pass

  @abstractmethod
  def createPhaseInteractionKernels(self):
    pass
