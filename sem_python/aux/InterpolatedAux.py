import numpy as np

from .AuxQuantity import AuxQuantity, AuxQuantityParameters
from ..input.Parameters import Parameters


class InterpolatedAuxParameters(AuxQuantityParameters):

    def __init__(self, factory):
        AuxQuantityParameters.__init__(self, factory)
        self.registerParameter("variable", "Name of the variable to interpolate")
        self.registerParameter(
            "dependencies", "List of names of solution variables that this aux may depend upon")
        self.registerIntParameter("n_test_functions", "Number of test functions that have support on a cell")


class InterpolatedAux(AuxQuantity):

    def __init__(self, params):
        AuxQuantity.__init__(self, params)
        self.variable = params.get("variable")
        self.name = self.variable
        self.dependencies = params.get("dependencies")
        self.n_test_functions = params.get("n_test_functions")

    def compute(self, nodal_data, nodal_der, elem_data, elem_der):
        elem_data[self.variable] = np.zeros(self.n_test_functions)
        for k in range(self.n_test_functions):
            elem_data[self.variable] += nodal_data[self.variable][k] * elem_data["phi"][k, :]

        for var in self.dependencies:
            elem_der[self.variable][var] = np.zeros(self.n_test_functions)
            for k in range(self.n_test_functions):
                elem_der[self.variable][var] += nodal_der[self.variable][var][k] * elem_data["phi"][
                    k, :]
