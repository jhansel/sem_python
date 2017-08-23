from AuxQuantity1Phase import AuxQuantity1Phase, AuxQuantity1PhaseParameters

class VelocityParameters(AuxQuantity1PhaseParameters):
  def __init__(self):
    AuxQuantity1PhaseParameters.__init__(self)

class Velocity(AuxQuantity1Phase):
  def __init__(self, params):
    AuxQuantity1Phase.__init__(self, params)
    self.name = self.u

  def compute(self, data, der):
    arho = data[self.arho]
    arhou = data[self.arhou]
    data[self.name] = arhou / arho

    du_darho = - arhou / arho / arho
    du_darhou = 1.0 / arho
    der[self.name] = {self.arho : du_darho, self.arhou : du_darhou}
