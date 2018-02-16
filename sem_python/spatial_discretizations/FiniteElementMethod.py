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

        # create stabilization, if any
        stabilization = None
        if params.has("Stabilization"):
            stabilization = self.factory.createObjectOfType(params.get("Stabilization"))
        else:
            stabilization = self.factory.createObject("NoStabilization")
        self.factory.storeObject(stabilization, "stabilization")

        # create assembly
        assembly_params = {"lump_mass_matrix": params.get("lump_mass_matrix")}
        assembly = self.factory.createObject("FEMAssembly", assembly_params)
        self.factory.storeObject(assembly, "assembly")
