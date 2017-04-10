import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/base")
from enums import PhaseType

sys.path.append(base_dir + "src/input")
from Parameters import Parameters

class InitialConditions1PhaseParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)
    self.registerFunctionParameter("p", "Pressure")
    self.registerFunctionParameter("T", "Temperature")
    self.registerFunctionParameter("u", "Velocity")

class InitialConditions1Phase(object):
  def __init__(self, params):
    self.p0 = {PhaseType.First : params.get("p")}
    self.T0 = {PhaseType.First : params.get("T")}
    self.u0 = {PhaseType.First : params.get("u")}
