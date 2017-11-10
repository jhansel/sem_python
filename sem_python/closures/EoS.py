from abc import ABCMeta, abstractmethod

from ..input.Parameters import Parameters

class EoSParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)

## Abstract base class for equations of state
class EoS(object, metaclass=ABCMeta):
  @abstractmethod
  def rho(self, p, T):
    pass

  @abstractmethod
  def rho_from_p_s(self, p, s):
    pass

  @abstractmethod
  def e(self, v, p):
    pass

  @abstractmethod
  def p(self, v, e):
    pass

  @abstractmethod
  def T(self, v, e):
    pass

  @abstractmethod
  def c(self, v, e):
    pass

  @abstractmethod
  def s(self, v, e):
    pass

  @abstractmethod
  def s_from_h_p(self, h, p):
    pass

  @abstractmethod
  def p_from_h_s(self, h, s):
    pass

  @abstractmethod
  def h(self, p, T):
    pass
