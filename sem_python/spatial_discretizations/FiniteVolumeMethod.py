from .SpatialDiscretization import SpatialDiscretization, SpatialDiscretizationParameters


class FiniteVolumeMethodParameters(SpatialDiscretizationParameters):

    def __init__(self, factory):
        SpatialDiscretizationParameters.__init__(self, factory)


## Finite volume method spatial discretization
class FiniteVolumeMethod(SpatialDiscretization):
    def __init__(self, params):
        SpatialDiscretization.__init__(self, params)

        # create assembly
        assembly_params = {}
        assembly = self.factory.createObject("FVMAssembly", assembly_params)
        self.factory.storeObject(assembly, "assembly")
