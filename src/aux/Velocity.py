from AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

class VelocityParameters(AuxQuantity1PhaseParameters):
  def __init__(self):
    AuxQuantity1PhaseParameters.__init__(self)

class Velocity(AuxQuantity1Phase):
  def __init__(self, params):
    AuxQuantity1Phase.__init__(self, params)

  def compute(self, data, der):
    arho = data[self.arho]
    arhou = data[self.arhou]
    data[self.u] = arhou / arho

    du_darho = - arhou / arho / arho
    du_darhou = 1.0 / arho
    der[self.u] = {self.arho : du_darho, self.arhou : du_darhou}
