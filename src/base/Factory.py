import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

# bc
sys.path.append(base_dir + "src/bc")
from FreeBC import FreeBC, FreeBCParameters
from DirichletVolumeFractionBC import DirichletVolumeFractionBC, DirichletVolumeFractionBCParameters
from InletRhoUBC import InletRhoUBC, InletRhoUBCParameters
from OutletBC import OutletBC, OutletBCParameters
from SolidWallBC import SolidWallBC, SolidWallBCParameters

# closures
sys.path.append(base_dir + "src/closures")
from IdealGasEoS import IdealGasEoS, IdealGasEoSParameters
from StiffenedGasEoS import StiffenedGasEoS, StiffenedGasEoSParameters
from InterfaceClosures import InterfaceClosures, InterfaceClosuresParameters

# executioners
sys.path.append(base_dir + "src/executioners")
from ImplicitEulerExecutioner import ImplicitEulerExecutioner, ImplicitEulerExecutionerParameters
from SteadyStateExecutioner import SteadyStateExecutioner, SteadyStateExecutionerParameters

# ic
sys.path.append(base_dir + "src/ic")
from InitialConditions1Phase import InitialConditions1Phase, InitialConditions1PhaseParameters
from InitialConditions2Phase import InitialConditions2Phase, InitialConditions2PhaseParameters

# input
sys.path.append(base_dir + "src/input")
from ModelParameters import ModelParameters
from PhysicsParameters import PhysicsParameters

# mesh
sys.path.append(base_dir + "src/mesh")
from UniformMesh import UniformMesh, UniformMeshParameters

# output
sys.path.append(base_dir + "src/output")
from Postprocessor import Postprocessor, PostprocessorParameters

# solvers
sys.path.append(base_dir + "src/solvers")
from NonlinearSolver import NonlinearSolver, NonlinearSolverParameters

# utilities
sys.path.append(base_dir + "src/utilities")
from error_utilities import error

## Class for creating objects
class Factory(object):
  ## Creates a parameters object
  # @param object_class  class of object for which to create parameters object
  # @param params  dictionary of parameter names to their values as strings
  def createParametersObject(self, object_class, params):
    # parameters classes are always named as the object class plus "Parameters"
    parameters_class = object_class + "Parameters"

    # parameters classes should have no arguments to their constructors
    if parameters_class in globals():
      constructor = globals()[parameters_class]
    else:
      error("'" + parameters_class + "' is not a valid object type.")
    parameters_object = constructor()

    # set each of the parameters
    for param in params:
      if param != "type":
        parameters_object.set(param, params[param])

    return parameters_object

  ## Creates an object
  # @param object_class  class of object to create
  # @param params  dictionary of parameter names to their values as strings
  def createObject(self, object_class, params, args=None):
    # create the object's parameters object first
    parameters_object = self.createParametersObject(object_class, params)

    # create the object
    if object_class in globals():
      constructor = globals()[object_class]
    else:
      error("'" + object_class + "' is not a valid object type.")
    if args:
      the_object = constructor(parameters_object, *args)
    else:
      the_object = constructor(parameters_object)

    return the_object
