import numpy as np

from ..base.enums import ModelType
from .FEMDoFHandler import FEMDoFHandler, FEMDoFHandlerParameters


class FEMDoFHandler1PhaseParameters(FEMDoFHandlerParameters):

    def __init__(self, factory):
        FEMDoFHandlerParameters.__init__(self, factory)


class FEMDoFHandler1Phase(FEMDoFHandler):

    def __init__(self, params):
        FEMDoFHandler.__init__(self, params)
        self.model_type = ModelType.OnePhase
        self.n_phases = 1
        self.n_vf_equations = 0
        self.setup()

    def aA1(self, U, k):
        return self.A[k]
