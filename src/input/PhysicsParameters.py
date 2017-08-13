from Parameters import Parameters

class PhysicsParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)
    self.registerFloatParameter("gravity", "X-component of gravitational acceleration")
