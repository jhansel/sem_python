from ..input.Parameters import Parameters


class SpatialDiscretizationParameters(Parameters):

    def __init__(self, factory):
        Parameters.__init__(self, factory)
        self.registerParameter("factory", "Factory")


## Base class for spatial discretizations
class SpatialDiscretization(object):
    def __init__(self, params):
        self.factory = params.get("factory")
