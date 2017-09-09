from abc import ABCMeta, abstractmethod
from numpy.linalg import norm

from Parameters import Parameters

class MeshParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)
    self.registerStringParameter("name", "Name to assign to mesh", "Mesh")
    self.registerFloatListParameter("orientation", "3-D orientation vector of mesh", [1,0,0])

class Mesh(object):
  __metaclass__ = ABCMeta
  def __init__(self, params):
    self.name = params.get("name")

    # get orientation vector and normalize it
    orientation = params.get("orientation")
    if len(orientation) != 3:
      error("Mesh orientation vector must have 3 elements")
    magnitude = norm(orientation)
    self.orientation = [x / magnitude for x in orientation]

  def getMinimumCellWidth(self):
    pass
