from enum import Enum

class ModelType(Enum):
  OnePhase = 1
  TwoPhaseNonInteracting = 2
  TwoPhase = 3

class PhaseType(Enum):
  First = 1
  Second = 2

class VariableName(Enum):
  ARho = 1
  ARhoU = 2
  ARhoE = 3
  VF1 = 4
