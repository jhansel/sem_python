# assembly
from ..assembly.FEMAssembly import FEMAssembly, FEMAssemblyParameters

# aux
from ..aux.AcousticImpedance import AcousticImpedance, AcousticImpedanceParameters
from ..aux.AuxGradient import AuxGradient, AuxGradientParameters
from ..aux.BerryInterfacePressure import BerryInterfacePressure, BerryInterfacePressureParameters
from ..aux.BerryInterfacePressureBar import BerryInterfacePressureBar, BerryInterfacePressureBarParameters
from ..aux.BerryInterfaceVelocity import BerryInterfaceVelocity, BerryInterfaceVelocityParameters
from ..aux.BerryInterfaceVelocityBar import BerryInterfaceVelocityBar, BerryInterfaceVelocityBarParameters
from ..aux.BerryInterfacialAreaDensity import BerryInterfacialAreaDensity, BerryInterfacialAreaDensityParameters
from ..aux.BerryPressureRelaxationCoef import BerryPressureRelaxationCoef, BerryPressureRelaxationCoefParameters
from ..aux.BerryVelocityRelaxationCoef import BerryVelocityRelaxationCoef, BerryVelocityRelaxationCoefParameters
from ..aux.ConstantAux import ConstantAux, ConstantAuxParameters
from ..aux.Density import Density, DensityParameters
from ..aux.EnergyFlux import EnergyFlux, EnergyFluxParameters
from ..aux.LaxFriedrichsCoefficient import LaxFriedrichsCoefficient, LaxFriedrichsCoefficientParameters
from ..aux.LaxFriedrichsCoefficientVolumeFraction import LaxFriedrichsCoefficientVolumeFraction, LaxFriedrichsCoefficientVolumeFractionParameters
from ..aux.EntropyMinimumMassFlux import EntropyMinimumMassFlux, EntropyMinimumMassFluxParameters
from ..aux.EntropyMinimumMomentumFlux import EntropyMinimumMomentumFlux, EntropyMinimumMomentumFluxParameters
from ..aux.EntropyMinimumEnergyFlux import EntropyMinimumEnergyFlux, EntropyMinimumEnergyFluxParameters
from ..aux.EntropyMinimumVolumeFractionFlux import EntropyMinimumVolumeFractionFlux, EntropyMinimumVolumeFractionFluxParameters
from ..aux.IdenticalAux import IdenticalAux, IdenticalAuxParameters
from ..aux.InternalEnergyDensity import InternalEnergyDensity, InternalEnergyDensityParameters
from ..aux.InterpolatedAux import InterpolatedAux, InterpolatedAuxParameters
from ..aux.MassFlux import MassFlux, MassFluxParameters
from ..aux.MomentumFlux import MomentumFlux, MomentumFluxParameters
from ..aux.Pressure import Pressure, PressureParameters
from ..aux.Temperature import Temperature, TemperatureParameters
from ..aux.SoundSpeed import SoundSpeed, SoundSpeedParameters
from ..aux.SpecificInternalEnergy import SpecificInternalEnergy, SpecificInternalEnergyParameters
from ..aux.SpecificTotalEnergy import SpecificTotalEnergy, SpecificTotalEnergyParameters
from ..aux.SpecificVolume import SpecificVolume, SpecificVolumeParameters
from ..aux.TestAux import TestAux, TestAuxParameters
from ..aux.Velocity import Velocity, VelocityParameters
from ..aux.VolumeFraction1Phase import VolumeFraction1Phase, VolumeFraction1PhaseParameters
from ..aux.VolumeFractionPhase1 import VolumeFractionPhase1, VolumeFractionPhase1Parameters
from ..aux.VolumeFractionPhase2 import VolumeFractionPhase2, VolumeFractionPhase2Parameters
from ..aux.AmbrosoBeta import AmbrosoBeta, AmbrosoBetaParameters
from ..aux.AmbrosoVelocityRelaxationCoef import AmbrosoVelocityRelaxationCoef, AmbrosoVelocityRelaxationCoefParameters
from ..aux.AmbrosoPressureRelaxationCoef import AmbrosoPressureRelaxationCoef, AmbrosoPressureRelaxationCoefParameters
from ..aux.AmbrosoInterfaceVelocity import AmbrosoInterfaceVelocity, AmbrosoInterfaceVelocityParameters
from ..aux.AmbrosoInterfacePressure import AmbrosoInterfacePressure, AmbrosoInterfacePressureParameters

# base
from .Model import Model, ModelParameters

# bc
from ..bc.FreeBC import FreeBC, FreeBCParameters
from ..bc.DirichletVolumeFractionBC import DirichletVolumeFractionBC, DirichletVolumeFractionBCParameters
from ..bc.InletP0T0BC import InletP0T0BC, InletP0T0BCParameters
from ..bc.InletRhoUBC import InletRhoUBC, InletRhoUBCParameters
from ..bc.OutletBC import OutletBC, OutletBCParameters
from ..bc.SolidWallBC import SolidWallBC, SolidWallBCParameters

# closures
from ..closures.IdealGasEoS import IdealGasEoS, IdealGasEoSParameters
from ..closures.StiffenedGasEoS import StiffenedGasEoS, StiffenedGasEoSParameters
from ..closures.TestEoS import TestEoS, TestEoSParameters
from ..closures.AmbrosoInterfaceClosures import AmbrosoInterfaceClosures, AmbrosoInterfaceClosuresParameters
from ..closures.BerryInterfaceClosures import BerryInterfaceClosures, BerryInterfaceClosuresParameters
from ..closures.ThermodynamicState import ThermodynamicState, ThermodynamicStateParameters

# executioners
from ..executioners.SteadyStateExecutioner import SteadyStateExecutioner, SteadyStateExecutionerParameters
from ..executioners.TransientExecutioner import TransientExecutioner, TransientExecutionerParameters

# fem
from ..fem.DoFHandler1Phase import DoFHandler1Phase, DoFHandler1PhaseParameters
from ..fem.DoFHandler2PhaseNonInteracting import DoFHandler2PhaseNonInteracting, DoFHandler2PhaseNonInteractingParameters
from ..fem.DoFHandler2Phase import DoFHandler2Phase, DoFHandler2PhaseParameters
from ..fem.FEValues import FEValues, FEValuesParameters
from ..fem.Quadrature import Quadrature, QuadratureParameters

# ic
from ..ic.InitialConditions1Phase import InitialConditions1Phase, InitialConditions1PhaseParameters
from ..ic.InitialConditions2Phase import InitialConditions2Phase, InitialConditions2PhaseParameters

# input
from ..input.HeatTransferData import HeatTransferData, HeatTransferDataParameters
from ..input.PhysicsParameters import PhysicsParameters

# junctions
from ..junctions.CloneJunction import CloneJunction, CloneJunctionParameters
from ..junctions.NewCompressibleJunction import NewCompressibleJunction, NewCompressibleJunctionParameters
from ..junctions.NewerCompressibleJunction import NewerCompressibleJunction, NewerCompressibleJunctionParameters
from ..junctions.NewestCompressibleJunction import NewestCompressibleJunction, NewestCompressibleJunctionParameters
from ..junctions.EqualFluxJunction import EqualFluxJunction, EqualFluxJunctionParameters
from ..junctions.EqualFluxLM1PhaseJunction import EqualFluxLM1PhaseJunction, EqualFluxLM1PhaseJunctionParameters
from ..junctions.EqualSolutionLM1PhaseJunction import EqualSolutionLM1PhaseJunction, EqualSolutionLM1PhaseJunctionParameters
from ..junctions.Junction1Phase import Junction1Phase, Junction1PhaseParameters
from ..junctions.TestJunction import TestJunction, TestJunctionParameters

# kernels
from ..kernels.InterpolatedAdvection import InterpolatedAdvection, InterpolatedAdvectionParameters
from ..kernels.DissipationAuxFlux import DissipationAuxFlux, DissipationAuxFluxParameters
from ..kernels.DissipationVariableGradient import DissipationVariableGradient, DissipationVariableGradientParameters
from ..kernels.VolumeFractionAdvection import VolumeFractionAdvection, VolumeFractionAdvectionParameters
from ..kernels.VolumeFractionPressureRelaxation import VolumeFractionPressureRelaxation, VolumeFractionPressureRelaxationParameters
from ..kernels.MassAdvection import MassAdvection, MassAdvectionParameters
from ..kernels.MomentumAdvection import MomentumAdvection, MomentumAdvectionParameters
from ..kernels.MomentumAreaGradient import MomentumAreaGradient, MomentumAreaGradientParameters
from ..kernels.MomentumGravity import MomentumGravity, MomentumGravityParameters
from ..kernels.MomentumVolumeFractionGradient import MomentumVolumeFractionGradient, MomentumVolumeFractionGradientParameters
from ..kernels.EnergyAdvection import EnergyAdvection, EnergyAdvectionParameters
from ..kernels.EnergyGravity import EnergyGravity, EnergyGravityParameters
from ..kernels.EnergyHeatTransfer import EnergyHeatTransfer, EnergyHeatTransferParameters
from ..kernels.EnergyPressureRelaxation import EnergyPressureRelaxation, EnergyPressureRelaxationParameters
from ..kernels.EnergyVolumeFractionGradient import EnergyVolumeFractionGradient, EnergyVolumeFractionGradientParameters

# mesh
from ..mesh.UniformMesh import UniformMesh, UniformMeshParameters

# output
from ..output.CSVOutput import CSVOutput, CSVOutputParameters
from ..output.PlotOutput import PlotOutput, PlotOutputParameters
from ..output.ScreenOutput import ScreenOutput, ScreenOutputParameters

# solvers
from ..solvers.NonlinearSolver import NonlinearSolver, NonlinearSolverParameters

# spatial discretizations
from ..spatial_discretizations.FiniteElementMethod import FiniteElementMethod, FiniteElementMethodParameters

# stabilization
from ..stabilization.NoStabilization import NoStabilization, NoStabilizationParameters
from ..stabilization.LaxFriedrichsStabilization import LaxFriedrichsStabilization, LaxFriedrichsStabilizationParameters

# time integrators
from ..time_integrators.ExplicitEuler import ExplicitEuler, ExplicitEulerParameters
from ..time_integrators.ImplicitEuler import ImplicitEuler, ImplicitEulerParameters

# time step sizers
from ..time_step_sizers.ConstantTimeStepSizer import ConstantTimeStepSizer, ConstantTimeStepSizerParameters
from ..time_step_sizers.CFLTimeStepSizer import CFLTimeStepSizer, CFLTimeStepSizerParameters

# utilities
from ..utilities.error_utilities import error


## Class for creating objects
class Factory(object):
    def __init__(self):
        self.stored_objects = {"factory": self}

    ##
    # Checks whether a parameter exists in factory warehouse
    #
    # @param[in] param   parameter to check
    # @returns true if parameter exists in factory warehouse
    #
    def hasParameter(self, param):
        return param in self.stored_objects

    ##
    # Gets a parameter from factory warehouse
    #
    # @param[in] param   parameter to get
    # @returns parameter from factory warehouse
    #
    def getParameter(self, param):
        return self.stored_objects[param]

    ## Creates a parameters object
    # @param object_class  class of object for which to create parameters object
    # @param params  dictionary of parameter names to their values as strings
    def createParametersObject(self, object_class, params=None):
        # parameters classes are always named as the object class plus "Parameters"
        parameters_class = object_class + "Parameters"

        # parameters classes should have no arguments to their constructors
        if parameters_class in globals():
            constructor = globals()[parameters_class]
        else:
            error("'" + parameters_class + "' is not a valid object type.")
        parameters_object = constructor(self)

        # set each of the parameters
        if params:
            for param in params:
                if param != "type":
                    parameters_object.set(param, params[param])

        return parameters_object

    ## Creates an object from its parameters object instead of a parameters dictionary
    # @param object_class  class of object to create
    # @param parameters_object  parameters object
    def createObjectFromParametersObject(self, object_class, parameters_object):
        # create the object
        if object_class in globals():
            constructor = globals()[object_class]
        else:
            error("'" + object_class + "' is not a valid object type.")

        the_object = constructor(parameters_object)

        return the_object

    ##
    # Creates an object, with class given as a parameter
    #
    # @param[in] object_class   class of object to create
    # @param[in] params   dictionary of parameter names to their values as strings
    # @returns created object
    #
    def createObject(self, object_class, params=None):
        # create the object's parameters object first
        parameters_object = self.createParametersObject(object_class, params)

        return self.createObjectFromParametersObject(object_class, parameters_object)

    ##
    # Creates an object, with class given by key "type"
    #
    # @param[in] params   dictionary of parameter names to their values as strings
    # @returns created object
    #
    def createObjectOfType(self, params):
        object_class = params["type"]
        return self.createObject(object_class, params)

    ## Creates and stores an object
    #
    # @param[in] the_object   the object to store
    # @param[in] object_name   string by which to identify the object in the warehouse
    #
    def storeObject(self, the_object, object_name):
        if object_name in self.stored_objects:
            error("There is already an object '" + object_name + "' stored in the factory")
        else:
            self.stored_objects[object_name] = the_object
