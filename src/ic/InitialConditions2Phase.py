import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/input")
from Parameters import Parameters

sys.path.append(base_dir + "src/utilities")
from error_utilities import error

class InitialConditions2PhaseParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)
    self.registerFunctionParameter("vf1", "Volume fraction of phase 1")
    self.registerFunctionParameter("rho1", "Density of phase 1")
    self.registerFunctionParameter("p1", "Pressure of phase 1")
    self.registerFunctionParameter("T1", "Temperature of phase 1")
    self.registerFunctionParameter("u1", "Velocity of phase 1")
    self.registerFunctionParameter("rho2", "Density of phase 2")
    self.registerFunctionParameter("p2", "Pressure of phase 2")
    self.registerFunctionParameter("T2", "Temperature of phase 2")
    self.registerFunctionParameter("u2", "Velocity of phase 2")

class InitialConditions2Phase(object):
  def __init__(self, params):
    self.vf0 = params.get("vf1")
    self.p0 = [params.get("p1"), params.get("p2")]
    self.u0 = [params.get("u1"), params.get("u2")]

    # one may supply either rho or T, but not both
    if params.has("rho1") and params.has("rho2"):
      has_rho = True
    else:
      has_rho = False
    if params.has("T1") and params.has("T2"):
      has_T = True
    else:
      has_T = False
    if (params.has("T1") and params.has("rho2")) or (params.has("rho1") and params.has("T2")):
      error("ICs cannot supply a mix of T and rho between the phases.")
    elif has_rho and has_T:
      error("ICs cannot specify both T and rho.")
    elif has_rho:
      self.rho0 = [params.get("rho1"), params.get("rho2")]
      self.specified_rho = True
    elif has_T:
      self.T0 = [params.get("T1"), params.get("T2")]
      self.specified_rho = False
    else:
      error("Either 'rho' or 'T' for each phase must be specified for IC.")
