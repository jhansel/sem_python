from .AuxQuantity import AuxQuantity, AuxQuantityParameters


class AuxQuantity2PhaseParameters(AuxQuantityParameters):

    def __init__(self, factory):
        AuxQuantityParameters.__init__(self, factory)


class AuxQuantity2Phase(AuxQuantity):

    def __init__(self, params):
        AuxQuantity.__init__(self, params)
