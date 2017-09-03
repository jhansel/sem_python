from AuxQuantity import AuxQuantity, AuxQuantityParameters

class ConstantAuxParameters(AuxQuantityParameters):
  def __init__(self):
    AuxQuantityParameters.__init__(self)
    self.registerStringParameter("name", "Name of aux quantity")
    self.registerFloatParameter("value", "Constant value of aux quantity")

class ConstantAux(AuxQuantity):
  def __init__(self, params):
    AuxQuantity.__init__(self, params)
    self.name = params.get("name")
    self.value = params.get("value")

  def compute(self, data, der):
    data[self.name] = self.value
