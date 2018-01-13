import numpy as np

from ..input.Parameters import Parameters


class InterpolatedAuxParameters(Parameters):

    def __init__(self):
        Parameters.__init__(self)
        self.registerParameter("variable", "Name of the variable to interpolate")
        self.registerParameter(
            "dependencies", "List of names of solution variables that this aux may depend upon")
        self.registerIntParameter("n_dof_per_cell_per_var", "Number of DoFs per cell per variable")


class InterpolatedAux(object):

    def __init__(self, params):
        self.variable = params.get("variable")
        self.dependencies = params.get("dependencies")
        self.n_dof_per_cell_per_var = params.get("n_dof_per_cell_per_var")

    def compute(self, nodal_data, nodal_der, elem_data, elem_der):
        elem_data[self.variable] = 0 * elem_data["phi"][0, :]
        for k in range(self.n_dof_per_cell_per_var):
            elem_data[self.variable] += nodal_data[self.variable][k] * elem_data["phi"][k, :]

        for var in self.dependencies:
            elem_der[self.variable][var] = 0 * elem_data["phi"][0, :]
            for k in range(self.n_dof_per_cell_per_var):
                elem_der[self.variable][var] += nodal_der[self.variable][var][k] * elem_data["phi"][
                    k, :]
