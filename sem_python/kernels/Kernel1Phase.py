from ..base.enums import ModelType
from .Kernel import Kernel, KernelParameters


class Kernel1PhaseParameters(KernelParameters):

    def __init__(self):
        KernelParameters.__init__(self)
        self.registerIntParameter("phase", "Phase index (0 or 1)")
        self.registerParameter("var_enum", "Variable enumeration")


class Kernel1Phase(Kernel):

    def __init__(self, params):
        self.phase = params.get("phase")
        self.var_enum = params.get("var_enum")
        dof_handler = params.get("dof_handler")
        params.set("var_index", dof_handler.variable_index[self.var_enum][self.phase])
        Kernel.__init__(self, params)

        # create list of relevant variable indices
        self.arhoA_index = self.dof_handler.arhoA_index[self.phase]
        self.arhouA_index = self.dof_handler.arhouA_index[self.phase]
        self.arhoEA_index = self.dof_handler.arhoEA_index[self.phase]
        if self.dof_handler.model_type == ModelType.TwoPhase:
            self.aA1_index = self.dof_handler.aA1_index[0]
            self.var_indices = [
                self.aA1_index, self.arhoA_index, self.arhouA_index, self.arhoEA_index
            ]
        else:
            self.aA1_index = float("NaN")
            self.var_indices = [self.arhoA_index, self.arhouA_index, self.arhoEA_index]

        # create variable names
        phase_str = str(self.phase + 1)
        self.vf = "vf" + phase_str
        self.arhoA = "arhoA" + phase_str
        self.arhouA = "arhouA" + phase_str
        self.arhoEA = "arhoEA" + phase_str
        self.rho = "rho" + phase_str
        self.u = "u" + phase_str
        self.E = "E" + phase_str
        self.v = "v" + phase_str
        self.e = "e" + phase_str
        self.p = "p" + phase_str
        self.T = "T" + phase_str

        self.grad_arhoA = "grad_" + self.arhoA
        self.grad_arhouA = "grad_" + self.arhouA
        self.grad_arhoEA = "grad_" + self.arhoEA
