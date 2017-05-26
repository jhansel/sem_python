import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/aux")
from AuxQuantity import AuxQuantity, AuxQuantityParameters

sys.path.append(base_dir + "src/utilities")
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
    AuxQuantity.__init__(self)
    self.var = params.get("var")
    self.other_vars = params.get("other_vars")
    self.coefs = params.get("coefs")
    self.b = params.get("b")
    if len(self.other_vars) != len(self.coefs):
      error("'other_vars' and 'coefs' list parameters must have same length.")
    self.n = len(self.other_vars)

  def compute(self, data, der):
    data[self.var] = self.b
    der[self.var] = dict()
    for i in xrange(self.n):
      data[self.var] += self.coefs[i] * data[self.other_vars[i]]
      der[self.var][self.other_vars[i]] = self.coefs[i]
