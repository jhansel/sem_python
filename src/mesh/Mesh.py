from abc import ABCMeta, abstractmethod

from Parameters import Parameters

class MeshParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)

class Mesh(object):
  __metaclass__ = ABCMeta
  def __init__(self, params):
    pass

  def getMinimumCellWidth(self):
    pass
