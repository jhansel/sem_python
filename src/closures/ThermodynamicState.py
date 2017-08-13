from Parameters import Parameters
from error_utilities import error

class ThermodynamicStateParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)
    self.registerFloatParameter("rho", "Density value")
    self.registerFloatParameter("p", "Pressure value")
    self.registerFloatParameter("T", "Temperature value")

## Defines a thermodynamic state
class ThermodynamicState(object):
  ## Constructor; determines provided values and computes rest
  def __init__(self, params):
    self.provided = dict()
    self.values = dict()
    self.getPropertyIfAvailable("rho", params)
    self.getPropertyIfAvailable("T", params)
    self.getPropertyIfAvailable("p", params)

    # count the number of provided properties
    n_provided_properties = 0
    for prop in self.provided:
      if self.provided[prop]:
        n_provided_properties += 1
    if n_provided_properties != 2:
      error("Exactly 2 thermodynamic properties must be provided.")

    # flag that all thermodynamic properties have been computed
    self.computed_all_properties = False

  ## Computes remaining thermodynamic properties
  # @param[in] eos  Equation of state object
  def computeRemainingProperties(self, eos):
    if self.provided["p"] and self.provided["T"]:
      self.values["rho"] = eos.rho(self.values["p"], self.values["T"])[0]
      self.values["v"] = 1.0 / self.values["rho"]
      self.values["e"] = eos.e(self.values["v"], self.values["p"])[0]
    elif self.provided["p"] and self.provided["rho"]:
      self.values["v"] = 1.0 / self.values["rho"]
      self.values["e"] = eos.e(self.values["v"], self.values["p"])[0]
      self.values["T"] = eos.T(self.values["v"], self.values["e"])[0]
    else:
      error("The provided combination of thermodynamic properties has not been implemented.")

    self.computed_all_properties = True

  ## Determines if a thermodynamic property was provided and gets it if it were
  # @param[in] prop  thermodyamic property to check
  # @param[in] params  dictionary of parameters defining thermodynamic state
  def getPropertyIfAvailable(self, prop, params):
    if params.has(prop):
      self.provided[prop] = True
      self.values[prop] = params.get(prop)
    else:
      self.provided[prop] = False

  ## Creates the thermodynamic state printout
  def __str__(self):
    # make sure that all properties have been computed
    if not self.computed_all_properties:
      error("The thermodynamic properties still need to be computed.")

    # create and return string
    printout = "\nFluid properties (SI Units):\n"
    printout += "  p   = " + str(self.values["p"]) + "\n"
    printout += "  T   = " + str(self.values["T"]) + "\n"
    printout += "  rho = " + str(self.values["rho"]) + "\n"
    printout += "  v   = " + str(self.values["v"]) + "\n"
    printout += "  e   = " + str(self.values["e"]) + "\n"
    return printout
