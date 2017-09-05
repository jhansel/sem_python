from Parameters import Parameters
from error_utilities import error

class InitialConditions1PhaseParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)
    def one(x):
      return 1
    self.registerParsedFunctionParameter("A", "Cross-sectional area of flow channel", one)
    self.registerParsedFunctionParameter("rho", "Density")
    self.registerParsedFunctionParameter("p", "Pressure")
    self.registerParsedFunctionParameter("T", "Temperature")
    self.registerParsedFunctionParameter("u", "Velocity")

class InitialConditions1Phase(object):
  def __init__(self, params):
    self.A = params.get("A")
    self.p = [params.get("p")]
    self.u = [params.get("u")]

    # one may supply either rho or T, but not both
    if params.has("rho") and params.has("T"):
      error("ICs cannot specify both T and rho.")
    elif params.has("rho"):
      self.rho = [params.get("rho")]
      self.specified_rho = True
    elif params.has("T"):
      self.T = [params.get("T")]
      self.specified_rho = False
    else:
      error("Either 'rho' or 'T' must be specified for IC.")
