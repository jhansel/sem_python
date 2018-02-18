from .enums import ModelType
from ..input.Parameters import Parameters


class ModelParameters(Parameters):

    def __init__(self, factory):
        Parameters.__init__(self, factory)
        self.registerStringSelectionParameter(
            "model", ["1phase", "2phase_noninteracting", "2phase"], "Identifier for the model")
        self.registerBoolParameter(
            "use_volume_fraction_pressure_relaxation",
            "Use the pressure relaxation term in the volume fraction equation", True)
        self.registerBoolParameter(
            "use_energy_pressure_relaxation",
            "Use the pressure relaxation term in the energy equation", True)


class Model(object):

    def __init__(self, params):
        model_string_to_type = {
            "1phase": ModelType.OnePhase,
            "2phase_noninteracting": ModelType.TwoPhaseNonInteracting,
            "2phase": ModelType.TwoPhase
        }
        self.model_type = model_string_to_type[params.get("model")]
        self.use_volume_fraction_pressure_relaxation = params.get(
            "use_volume_fraction_pressure_relaxation")
        self.use_energy_pressure_relaxation = params.get("use_energy_pressure_relaxation")

        if self.model_type == ModelType.OnePhase:
            self.phases = [0]
        else:
            self.phases = [0, 1]
        self.n_phases = len(self.phases)

        if self.model_type == ModelType.TwoPhase:
            self.phase_interaction = True
        else:
            self.phase_interaction = False
