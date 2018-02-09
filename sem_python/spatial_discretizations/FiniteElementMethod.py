from .SpatialDiscretization import SpatialDiscretization, SpatialDiscretizationParameters


class FiniteElementMethodParameters(SpatialDiscretizationParameters):

    def __init__(self, factory):
        SpatialDiscretizationParameters.__init__(self, factory)
        self.registerNamedSubblock("Stabilization")


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
