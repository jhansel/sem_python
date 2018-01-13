import numpy as np

from ..base.enums import ModelType, VariableName
from .DoFHandler import DoFHandler, DoFHandlerParameters


class DoFHandler2PhaseParameters(DoFHandlerParameters):

    def __init__(self):
        DoFHandlerParameters.__init__(self)


class DoFHandler2Phase(DoFHandler):

    def __init__(self, params):
        DoFHandler.__init__(self, params)
        self.model_type = ModelType.TwoPhase
        self.n_phases = 2
        self.n_vf_equations = 1
        self.setup()

    def aA1(self, U, k):
        return U[self.i(k, self.aA1_index[0])]
