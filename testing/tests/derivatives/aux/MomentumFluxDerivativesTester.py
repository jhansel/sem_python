import unittest

from sem_python.aux.MomentumFlux import MomentumFlux, MomentumFluxParameters
from sem_python.aux.TestAux import TestAux, TestAuxParameters
from sem_python.aux.VolumeFractionPhase1 import VolumeFractionPhase1, VolumeFractionPhase1Parameters
from ....src.testers.AuxDerivativesTester import AuxDerivativesTester

# test aux
params = MomentumFluxParameters()
params.set("phase", 0)
test_aux = MomentumFlux(params)

# volume fraction aux
params = VolumeFractionPhase1Parameters()
params.set("phase", 0)
vf_aux = VolumeFractionPhase1(params)

# velocity aux
params = TestAuxParameters()
params.set("var", "u1")
params.set("other_vars", ["arhoA1", "arhouA1"])
params.set("coefs", [1.3, 2.2])
u_aux = TestAux(params)

# pressure aux
params = TestAuxParameters()
params.set("var", "p1")
params.set("other_vars", ["aA1", "arhoA1", "arhouA1", "arhoEA1"])
params.set("coefs", [1.4, 2.3, 2.1, 1.2])
p_aux = TestAux(params)

other_aux = [vf_aux, u_aux, p_aux]
root_vars = ["aA1", "arhoA1", "arhouA1", "arhoEA1"]

class MomentumFluxDerivativesTester(unittest.TestCase):
  def setUp(self):
    self.derivatives_tester = AuxDerivativesTester()

  def test(self):
    rel_diffs = self.derivatives_tester.checkDerivatives(test_aux, other_aux, root_vars, constant_data={"A": 0.3})
    for key in rel_diffs:
      self.assertLessEqual(rel_diffs[key], 1e-6)
