from AuxQuantity import AuxQuantity, AuxQuantityParameters

from error_utilities import error

class TestAuxParameters(AuxQuantityParameters):
  def __init__(self):
    AuxQuantityParameters.__init__(self)
    self.registerStringParameter("var", "Name of this variable")
    self.registerListParameter("other_vars", "Variable names this quantity depends upon")
    self.registerListParameter("coefs", "Proportionality coefficients for each quantity")
    self.registerFloatParameter("b", "Constant shift", 1.0)

class TestAux(AuxQuantity):
  def __init__(self, params):
    AuxQuantity.__init__(self, params)
    self.name = params.get("var")
    self.other_vars = params.get("other_vars")
    self.coefs = params.get("coefs")
    self.b = params.get("b")
    if len(self.other_vars) != len(self.coefs):
      error("'other_vars' and 'coefs' list parameters must have same length.")
    self.n = len(self.other_vars)

  def compute(self, data, der):
    data[self.name] = self.b
    der[self.name] = dict()
    for i in xrange(self.n):
      data[self.name] += self.coefs[i] * data[self.other_vars[i]]
      der[self.name][self.other_vars[i]] = self.coefs[i]
