import numpy as np

from ..base.enums import ModelType, VariableName
from .FEMDoFHandler import FEMDoFHandler, FEMDoFHandlerParameters


class FEMDoFHandler2PhaseParameters(FEMDoFHandlerParameters):

    def __init__(self, factory):
        FEMDoFHandlerParameters.__init__(self, factory)


class FEMDoFHandler2Phase(FEMDoFHandler):

    def __init__(self, params):
        FEMDoFHandler.__init__(self, params)
        self.model_type = ModelType.TwoPhase
        self.n_phases = 2
        self.n_vf_equations = 1
        self.setup()

    def aA1(self, U, k):
        return U[self.i(k, self.aA1_index[0])]
