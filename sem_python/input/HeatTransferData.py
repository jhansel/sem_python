from .Parameters import Parameters

class HeatTransferDataParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)
    self.registerFloatParameter("htc_wall", "Wall heat transfer coefficient")
    self.registerFloatParameter("T_wall", "Wall temperature")
    self.registerFloatParameter("P_heat", "Heated perimeter")

class HeatTransferData(object):
  def __init__(self, params):
    self.htc_wall = params.get("htc_wall")
    self.T_wall = params.get("T_wall")
    self.P_heat = params.get("P_heat")
