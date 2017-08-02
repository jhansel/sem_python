import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/base")
from enums import VariableName

sys.path.append(base_dir + "src/stabilization")
from Stabilization import Stabilization, StabilizationParameters

class LaxFriedrichsStabilizationParameters(StabilizationParameters):
  def __init__(self):
    StabilizationParameters.__init__(self)
    self.registerFloatParameter("mult", "Multiplier for the viscous coefficient", 1.0)
    self.registerBoolParameter("use_simple_dissipation", "Use simple form of dissipation", False)

## Lax-Friedrichs stabilization
class LaxFriedrichsStabilization(Stabilization):
  def __init__(self, params):
    Stabilization.__init__(self, params)
    self.mult = params.get("mult")
    self.use_simple_dissipation = params.get("use_simple_dissipation")

  def needSolutionGradients(self):
    return True

  def createIndependentPhaseAuxQuantities(self, phase):
    aux_list = list()
    var_list = ["vf", "arho", "arhou", "arhoE"]

    # add the viscous coefficients
    for var in var_list:
      params = {"phase": phase, "var": var, "mult": self.mult}
      aux_list.append(self.factory.createObject("LaxFriedrichsCoefficient", params))

    if not self.use_simple_dissipation:
      # internal energy density
      params = {"phase": phase}
      aux_list.append(self.factory.createObject("InternalEnergyDensity", params))

      # gradients
      aux_gradient_names = ["vf", "rho", "u", "rhoe"]
      for aux_gradient_name in aux_gradient_names:
        params = {"aux": aux_gradient_name + str(phase + 1)}
        aux_list.append(self.factory.createObject("AuxGradient", params))

      # add the viscous fluxes
      flux_classes = ["EntropyMinimumVolumeFractionFlux", "EntropyMinimumMassFlux",
        "EntropyMinimumMomentumFlux", "EntropyMinimumEnergyFlux"]
      params = {"phase": phase}
      for flux_class in flux_classes:
        aux_list.append(self.factory.createObject(flux_class, params))

    return aux_list

  def createPhaseInteractionAuxQuantities(self):
    return []

  def createIndependentPhaseKernels(self, phase):
    kernels = list()
    args = tuple([self.dof_handler])
    var_enums = [VariableName.ARho, VariableName.ARhoU, VariableName.ARhoE]
    for var_enum in var_enums:
      if self.use_simple_dissipation:
        params = {"phase": phase, "var": var_enum}
        kernels.append(self.factory.createObject("DissipationVariableGradient", params, args))
      else:
        var_name = self.dof_handler.variableEnumToName(var_enum, phase)
        flux_name = "viscflux_" + var_name
        params = {"phase": phase, "var": var_enum, "flux_name": flux_name}
        kernels.append(self.factory.createObject("DissipationAuxFlux", params, args))
    return kernels

  def createPhaseInteractionKernels(self):
    kernels = list()
    args = tuple([self.dof_handler])
    if self.use_simple_dissipation:
      params = {"phase": 0, "var": VariableName.VF1}
      kernels.append(self.factory.createObject("DissipationVariableGradient", params, args))
    else:
      params = {"phase": 0, "var": VariableName.VF1, "flux_name": "viscflux_vf1"}
      kernels.append(self.factory.createObject("DissipationAuxFlux", params, args))
    return kernels
