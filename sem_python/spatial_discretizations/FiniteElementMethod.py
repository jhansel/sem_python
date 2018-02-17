from ..base.enums import ModelType
from .SpatialDiscretization import SpatialDiscretization, SpatialDiscretizationParameters


class FiniteElementMethodParameters(SpatialDiscretizationParameters):

    def __init__(self, factory):
        SpatialDiscretizationParameters.__init__(self, factory)
        self.registerNamedSubblock("Stabilization")
        self.registerBoolParameter("lump_mass_matrix", "Lump the mass matrix?", False)


## Finite element method spatial discretizations
class FiniteElementMethod(SpatialDiscretization):
    def __init__(self, params):
        SpatialDiscretization.__init__(self, params)

        if params.has("Stabilization"):
            self.has_stabilization = True
            self.stabilization_parameters = params.get("Stabilization")
        else:
            self.has_stabilization = False

        self.lump_mass_matrix = params.get("lump_mass_matrix")

    def createDoFHandler(self):
        dof_handler = self.factory.createObject("FEMDoFHandler")
        self.factory.storeObject(dof_handler, "dof_handler")

    def createAssemblyObjects(self):
        # quadrature
        quadrature = self.factory.createObject("Quadrature", {"n_q_points": 2})
        self.factory.storeObject(quadrature, "quadrature")

        # create stabilization, if any
        if self.has_stabilization:
            stabilization = self.factory.createObjectOfType(self.stabilization_parameters)
        else:
            stabilization = self.factory.createObject("NoStabilization")
        self.factory.storeObject(stabilization, "stabilization")

        # create assembly
        assembly_params = {"lump_mass_matrix": self.lump_mass_matrix}
        assembly = self.factory.createObject("FEMAssembly", assembly_params)
        self.factory.storeObject(assembly, "assembly")
