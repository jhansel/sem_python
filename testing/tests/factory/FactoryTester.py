import unittest

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

  ## Tests the createObjectFromParametersObject() function
  def testCreateObjectFromParametersObject(self):
    object_class = "IdealGasEoS"
    parameters_object = self.factory.createParametersObject(object_class)
    parameters_object.set("gamma", 1.4)
    parameters_object.set("R", 200)
    self.factory.createObjectFromParametersObject(object_class, parameters_object)
