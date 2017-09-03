from abc import ABCMeta, abstractmethod

from Parameters import Parameters

class MeshParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)
    self.registerStringParameter("name", "Name to assign to mesh", "Mesh")

class Mesh(object):
  __metaclass__ = ABCMeta
  def __init__(self, params):
    self.name = params.get("name")

  def getMinimumCellWidth(self):
    pass
