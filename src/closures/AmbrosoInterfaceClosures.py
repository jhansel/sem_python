import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/closures")
from InterfaceClosures import InterfaceClosures, InterfaceClosuresParameters

class AmbrosoInterfaceClosuresParameters(InterfaceClosuresParameters):
  def __init__(self):
    InterfaceClosuresParameters.__init__(self)
    self.registerFloatParameter("chi", "Weight fraction for phase 1", 0.5)
    self.registerFloatParameter("pressure_relaxation_time", "Relaxation time for pressures")

class AmbrosoInterfaceClosures(InterfaceClosures):
  def __init__(self, params):
    InterfaceClosures.__init__(self, params)
    self.pressure_relaxation_time = params.get("pressure_relaxation_time")
    self.chi = params.get("chi") # should be in (0,1)

  def createPhaseInteractionAuxQuantities(self):
    interaction_aux_names = ["AmbrosoBeta", "AmbrosoMu", "AmbrosoTheta", "AmbrosoInterfaceVelocity", "AmbrosoInterfacePressure"]
    interaction_aux = list()
    for aux_name in interaction_aux_names:
      params = dict()
      if aux_name == "AmbrosoBeta":
        params["chi"] = self.chi
      elif aux_name == "AmbrosoTheta":
        params["pressure_relaxation_time"] = self.pressure_relaxation_time
      interaction_aux.append(self.factory.createObject(aux_name, params))

    params = {"original_aux": "pI", "copy_aux": "pI_bar"}
    interaction_aux.append(self.factory.createObject("IdenticalAux", params))

    return interaction_aux
