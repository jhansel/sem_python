import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/base")
from enums import PhaseType

sys.path.append(base_dir + "src/input")
from Parameters import Parameters

class InitialConditions2PhaseParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)
    self.registerFunctionParameter("vf1", "Volume fraction of phase 1")
    self.registerFunctionParameter("p1", "Pressure of phase 1")
    self.registerFunctionParameter("T1", "Temperature of phase 1")
    self.registerFunctionParameter("u1", "Velocity of phase 1")
    self.registerFunctionParameter("p2", "Pressure of phase 2")
    self.registerFunctionParameter("T2", "Temperature of phase 2")
    self.registerFunctionParameter("u2", "Velocity of phase 2")

class InitialConditions2Phase(object):
  def __init__(self, params):
    self.vf0 = params.get("vf1")
    self.p0 = {PhaseType.First : params.get("p1"), PhaseType.Second : params.get("p2")}
    self.T0 = {PhaseType.First : params.get("T1"), PhaseType.Second : params.get("T2")}
    self.u0 = {PhaseType.First : params.get("u1"), PhaseType.Second : params.get("u2")}
