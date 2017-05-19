import sys
sys.path.append("../../../src/base")
sys.path.append("../../../src/bc")
sys.path.append("../../../src/closures")
sys.path.append("../../../src/fem")

import unittest

from DoFHandler import DoFHandler
from enums import ModelType, PhaseType
from Factory import Factory

class FactoryTester(unittest.TestCase):
  def setUp(self):
    self.factory = Factory()

  ## Tests the createParametersObject() function without a type parameter
  def testCreateParametersObjectWithNoTypeParam(self):
    params = {"gamma" : "1.4", "R" : "200"}
    self.factory.createParametersObject("IdealGasEoS", params)

  ## Tests the createParametersObject() function with a type parameter
  def testCreateParametersObjectWithTypeParam(self):
    params = {"type" : "IdealGasEoS", "gamma" : "1.4", "R" : "200"}
    object_class = params["type"]
    self.factory.createParametersObject(object_class, params)

  ## Tests the createObject() function
  def testCreateObject(self):
    params = {"type" : "IdealGasEoS", "gamma" : "1.4", "R" : "200"}
    object_class = params["type"]
    self.factory.createObject(object_class, params)

  ## Tests the createObject() function when arguments in addition to input parameters are required
  def testCreateObjectWithExtraArgs(self):
    mesh_params = {"n_cell" : "5"}
    mesh = self.factory.createObject("UniformMesh", mesh_params)
    dof_handler = DoFHandler(mesh, ModelType.OnePhase, None)
    eos_params = {"gamma" : "1.4", "R" : "200"}
    eos = self.factory.createParametersObject("IdealGasEoS", eos_params)
    eos_map = {PhaseType.First : eos}
    bc_params = {"boundary" : "left", "p" : "1e5"}
    bc_params["phase"] = PhaseType.First
    args = (dof_handler, eos_map)
    self.factory.createObject("OutletBC", bc_params, args)
