from .AuxQuantity import AuxQuantity, AuxQuantityParameters


class IdenticalAuxParameters(AuxQuantityParameters):

    def __init__(self):
        AuxQuantityParameters.__init__(self)
        self.registerStringParameter("original_aux", "Name of original aux quantity")
        self.registerStringParameter("copy_aux", "Name of copy aux quantity")


class IdenticalAux(AuxQuantity):

    def __init__(self, params):
        AuxQuantity.__init__(self, params)
        self.original_aux = params.get("original_aux")
        self.name = params.get("copy_aux")

    def compute(self, data, der):
        data[self.name] = data[self.original_aux]
        der[self.name] = der[self.original_aux]
