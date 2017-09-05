from enum import Enum

class ModelType(Enum):
  OnePhase = 1
  TwoPhaseNonInteracting = 2
  TwoPhase = 3

class VariableName(Enum):
  ARhoA = 1
  ARhoUA = 2
  ARhoEA = 3
  AA1 = 4
