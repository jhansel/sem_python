import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/input")
from Parameters import Parameters

sys.path.append(base_dir + "src/utilities")
from error_utilities import error

class InitialConditions1PhaseParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)
    self.registerFunctionParameter("rho", "Density")
    self.registerFunctionParameter("p", "Pressure")
    self.registerFunctionParameter("T", "Temperature")
    self.registerFunctionParameter("u", "Velocity")

class InitialConditions1Phase(object):
  def __init__(self, params):
    self.p0 = [params.get("p")]
    self.u0 = [params.get("u")]

    # one may supply either rho or T, but not both
    if params.has("rho") and params.has("T"):
      error("ICs cannot specify both T and rho.")
    elif params.has("rho"):
      self.rho0 = [params.get("rho")]
      self.specified_rho = True
    elif params.has("T"):
      self.T0 = [params.get("T")]
      self.specified_rho = False
    else:
      error("Either 'rho' or 'T' must be specified for IC.")
