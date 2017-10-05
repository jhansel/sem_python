from ..input.Parameters import Parameters

class PhysicsParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)
    self.registerFloatListParameter("gravity", "3-D gravitational acceleration vector")
