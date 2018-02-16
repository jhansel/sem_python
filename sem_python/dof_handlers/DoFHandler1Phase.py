import numpy as np

from ..base.enums import ModelType
from .DoFHandler import DoFHandler, DoFHandlerParameters


class DoFHandler1PhaseParameters(DoFHandlerParameters):

    def __init__(self, factory):
        DoFHandlerParameters.__init__(self, factory)


class DoFHandler1Phase(DoFHandler):

    def __init__(self, params):
        DoFHandler.__init__(self, params)
        self.model_type = ModelType.OnePhase
        self.n_phases = 1
        self.n_vf_equations = 0
        self.setup()

    def aA1(self, U, k):
        return self.A[k]
