import unittest

from AmbrosoPressureRelaxationCoef import AmbrosoPressureRelaxationCoef, AmbrosoPressureRelaxationCoefParameters
from TestAux import TestAux, TestAuxParameters
from AuxDerivativesTester import AuxDerivativesTester

# pressure relaxation aux
params = AmbrosoPressureRelaxationCoefParameters()
params.set("pressure_relaxation_time", 2.0)
test_aux = AmbrosoPressureRelaxationCoef(params)

# phase-1 volume fraction aux
params = TestAuxParameters()
params.set("var", "vf1")
params.set("other_vars", ["aA1"])
params.set("coefs", [1.3])
vf1_aux = TestAux(params)

# phase-1 pressure aux
params = TestAuxParameters()
params.set("var", "p1")
params.set("other_vars", ["aA1", "arhoA1", "arhouA1", "arhoEA1"])
params.set("coefs", [1.2, 2.2, 3.2, 4.2])
p1_aux = TestAux(params)

# phase-2 pressure aux
params = TestAuxParameters()
params.set("var", "p2")
params.set("other_vars", ["aA1", "arhoA2", "arhouA2", "arhoEA2"])
params.set("coefs", [1.5, 2.5, 3.5, 4.5])
p2_aux = TestAux(params)

other_aux = [vf1_aux, p1_aux, p2_aux]
root_vars = ["aA1", "arhoA1", "arhouA1", "arhoEA1", "arhoA2", "arhouA2", "arhoEA2"]

class AmbrosoPressureRelaxationCoefDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 5e-6)

if __name__ == "__main__":
  derivatives_tester = AuxDerivativesTester(True)
  _ = derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars)
