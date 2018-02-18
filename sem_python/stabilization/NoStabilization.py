from .Stabilization import Stabilization, StabilizationParameters


class NoStabilizationParameters(StabilizationParameters):

    def __init__(self, factory):
        StabilizationParameters.__init__(self, factory)


## No stabilization
class NoStabilization(Stabilization):

    def __init__(self, params):
        Stabilization.__init__(self, params)

    def needSolutionGradients(self):
        return False

    def createAuxQuantities(self):
        return []

    def createKernels(self):
        return []
