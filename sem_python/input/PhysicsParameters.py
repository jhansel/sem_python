from ..input.Parameters import Parameters


class PhysicsParameters(Parameters):

    def __init__(self, factory):
        Parameters.__init__(self, factory)
        self.registerFloatListParameter("gravity", "3-D gravitational acceleration vector")
