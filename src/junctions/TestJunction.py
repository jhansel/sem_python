from Junction import Junction, JunctionParameters

class TestJunctionParameters(JunctionParameters):
  def __init__(self):
    JunctionParameters.__init__(self)
    self.registerIntParameter("n_constraints", "Number of constraints to add for this test junction")

class TestJunction(Junction):
  def __init__(self, params):
    Junction.__init__(self, params)
    self.n_constraints += params.get("n_constraints")
