import numpy as np

from ..base.enums import ModelType, VariableName
from ..closures.thermodynamic_functions import computeVolumeFraction, \
  computeVelocity, computeDensity, computeSpecificVolume, \
  computeSpecificTotalEnergy, computeSpecificInternalEnergy
from ..input.Parameters import Parameters

class ExecutionerParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)
    self.registerParameter("model", "Model")
    self.registerParameter("ics", "Initial conditions")
    self.registerParameter("bcs", "Boundary conditions")
    self.registerParameter("junctions", "List of junctions")
    self.registerParameter("eos_list", "List of equations of state")
    self.registerParameter("interface_closures", "Interface closures")
    self.registerFloatParameter("gravity", "Acceleration due to gravity")
    self.registerParameter("dof_handler", "Degree of freedom handler")
    self.registerParameter("meshes", "List of meshes")
    self.registerParameter("nonlinear_solver_params", "Nonlinear solver parameters")
    self.registerParameter("stabilization", "Stabilization")
    self.registerParameter("factory", "Factory")
    self.registerBoolParameter("split_source", "Use source-splitting?", False)

class Executioner(object):
  def __init__(self, params):
    self.model = params.get("model")
    ics = params.get("ics")
    self.model_type = self.model.model_type
    self.bcs = params.get("bcs")
    self.junctions = params.get("junctions")
    self.eos_list = params.get("eos_list")
    interface_closures = params.get("interface_closures")
    self.gravity = params.get("gravity")
    self.dof_handler = params.get("dof_handler")
    self.meshes = params.get("meshes")
    self.nonlinear_solver_params = params.get("nonlinear_solver_params")
    self.factory = params.get("factory")
    stabilization = params.get("stabilization")
    self.split_source = params.get("split_source")
    self.need_solution_gradients = stabilization.needSolutionGradients()

    # quadrature
    quadrature_params = {}
    self.quadrature = self.factory.createObject("Quadrature", quadrature_params)

    # FE values
    fe_values_params = {"quadrature": self.quadrature, "dof_handler": self.dof_handler, "meshes": self.meshes}
    self.fe_values = self.factory.createObject("FEValues", fe_values_params)

    # initialize the solution
    self.U = np.zeros(self.dof_handler.n_dof)
    self.initializePhaseSolution(ics, 0)
    if self.model_type != ModelType.OnePhase:
      self.initializePhaseSolution(ics, 1)
    if self.model_type == ModelType.TwoPhase:
      self.initializeVolumeFractionSolution(ics)

    # set local solution update function
    if self.model_type == ModelType.OnePhase:
      self.computeLocalCellSolution = self.computeLocalCellSolutionOnePhase
      self.computeLocalNodeSolution = self.computeLocalNodeSolutionOnePhase
    else:
      self.computeLocalCellSolution = self.computeLocalCellSolutionTwoPhase
      self.computeLocalNodeSolution = self.computeLocalNodeSolutionTwoPhase

    # create aux quantities
    self.aux_list = self.createIndependentPhaseAuxQuantities(0)
    if self.model_type != ModelType.OnePhase:
      self.aux_list += self.createIndependentPhaseAuxQuantities(1)
    if self.model_type == ModelType.TwoPhase:
      self.aux_list += interface_closures.createAuxQuantities()
    self.aux_list += stabilization.createAuxQuantities()

    # get list of aux quantities
    self.aux_names = [aux.name for aux in self.aux_list]

    # create list of source kernels
    self.source_kernels = self.createIndependentPhaseSourceKernels(0)
    if self.model_type != ModelType.OnePhase:
      self.source_kernels += self.createIndependentPhaseSourceKernels(1)
    if self.model_type == ModelType.TwoPhase:
      self.source_kernels += self.createPhaseInteractionSourceKernels()

    # create advection kernels
    self.fem_kernels = self.createIndependentPhaseAdvectionKernels(0) + stabilization.createIndependentPhaseKernels(0)
    if self.model_type != ModelType.OnePhase:
      self.fem_kernels += self.createIndependentPhaseAdvectionKernels(1) + stabilization.createIndependentPhaseKernels(1)
    if self.model_type == ModelType.TwoPhase:
      self.fem_kernels += self.createPhaseInteractionAdvectionKernels() + stabilization.createPhaseInteractionKernels()

    # add source kernels to kernel lists if not using source-splitting
    if not self.split_source:
      self.fem_kernels += self.source_kernels

  def initializePhaseSolution(self, ics, phase):
    # get appropriate volume fraction function
    if self.model_type == ModelType.OnePhase:
      def initial_vf(x):
        return 1
    else:
      initial_vf1 = ics.vf1
      if phase == 0:
        initial_vf = initial_vf1
      else:
        def initial_vf(x):
          return 1 - initial_vf1(x)

    # get relevant IC functions
    initial_p = ics.p[phase]
    initial_u = ics.u[phase]
    if ics.specified_rho:
      initial_rho = ics.rho[phase]
    else:
      initial_T = ics.T[phase]

    # compute IC
    eos_phase = self.eos_list[phase]
    arho_index = self.dof_handler.variable_index[VariableName.ARho][phase]
    arhou_index = self.dof_handler.variable_index[VariableName.ARhoU][phase]
    arhoE_index = self.dof_handler.variable_index[VariableName.ARhoE][phase]
    for k in xrange(self.dof_handler.n_node):
      vf = initial_vf(self.dof_handler.x[k])
      p = initial_p(self.dof_handler.x[k])
      u = initial_u(self.dof_handler.x[k])
      if ics.specified_rho:
        rho = initial_rho(self.dof_handler.x[k])
      else:
        T = initial_T(self.dof_handler.x[k])
        rho = eos_phase.rho(p, T)
      e = eos_phase.e(1.0 / rho, p)[0]
      E = e + 0.5 * u * u
      self.U[self.dof_handler.i(k, arho_index)] = vf * rho
      self.U[self.dof_handler.i(k, arhou_index)] = vf * rho * u
      self.U[self.dof_handler.i(k, arhoE_index)] = vf * rho * E

  def initializeVolumeFractionSolution(self, ics):
    vf1_index = self.dof_handler.variable_index[VariableName.VF1][0]
    for k in xrange(self.dof_handler.n_node):
      self.U[self.dof_handler.i(k, vf1_index)] = ics.vf1(self.dof_handler.x[k])

  def createIndependentPhaseAuxQuantities(self, phase):
    # create list of aux quantities to create
    aux_names_phase = list()
    if phase == 0:
      if self.model_type == ModelType.OnePhase:
        aux_names_phase.append("VolumeFraction1Phase")
      else:
        aux_names_phase.append("VolumeFractionPhase1")
    else:
      aux_names_phase.append("VolumeFractionPhase2")
    aux_names_phase += ["Velocity", "SpecificTotalEnergy", "Density", \
      "SpecificVolume", "SpecificInternalEnergy", "Pressure", "Temperature", "SoundSpeed"]

    # create the aux quantities for this phase
    aux_list = list()
    for aux_name in aux_names_phase:
      params = {"phase": phase}
      if aux_name == "Pressure":
        params["p_function"] = self.eos_list[phase].p
      elif aux_name == "Temperature":
        params["T_function"] = self.eos_list[phase].T
      elif aux_name == "SoundSpeed":
        params["c_function"] = self.eos_list[phase].c
      aux_list.append(self.factory.createObject(aux_name, params))

    return aux_list

  def createIndependentPhaseAdvectionKernels(self, phase):
    kernels = list()
    params = {"phase": phase, "dof_handler": self.dof_handler}
    kernel_name_list = ["MassAdvection", "MomentumAdvection", "EnergyAdvection"]
    for kernel_name in kernel_name_list:
      kernels.append(self.factory.createObject(kernel_name, params))
    return kernels

  def createIndependentPhaseSourceKernels(self, phase):
    params = {"phase": phase, "dof_handler": self.dof_handler, "is_nodal": self.split_source}
    kernel_names = ["MomentumGravity", "EnergyGravity"]
    kernels = [self.factory.createObject(kernel_name, params) for kernel_name in kernel_names]
    return kernels

  def createPhaseInteractionAdvectionKernels(self):
    params1 = {"phase": 0, "dof_handler": self.dof_handler}
    params2 = {"phase": 1, "dof_handler": self.dof_handler}

    kernels = list()
    kernels.append(self.factory.createObject("VolumeFractionAdvection", params1))
    kernels.append(self.factory.createObject("MomentumVolumeFractionGradient", params1))
    kernels.append(self.factory.createObject("MomentumVolumeFractionGradient", params2))
    kernels.append(self.factory.createObject("EnergyVolumeFractionGradient", params1))
    kernels.append(self.factory.createObject("EnergyVolumeFractionGradient", params2))

    return kernels

  def createPhaseInteractionSourceKernels(self):
    params1 = {"phase": 0, "dof_handler": self.dof_handler, "is_nodal": self.split_source}
    params2 = {"phase": 1, "dof_handler": self.dof_handler, "is_nodal": self.split_source}

    kernels = list()
    kernels.append(self.factory.createObject("VolumeFractionPressureRelaxation", params1))
    kernels.append(self.factory.createObject("EnergyPressureRelaxation", params1))
    kernels.append(self.factory.createObject("EnergyPressureRelaxation", params2))

    return kernels

  # computes the steady-state residual and Jacobian without applying strong BC
  def assembleSteadyStateSystemWithoutStrongBC(self, U):
    r = np.zeros(self.dof_handler.n_dof)
    J = np.zeros(shape=(self.dof_handler.n_dof, self.dof_handler.n_dof))

    # volumetric terms
    self.addSteadyStateSystem(U, r, J)

    # BCs
    for bc in self.bcs:
      bc.applyWeakBC(U, r, J)

    # junctions
    for junction in self.junctions:
      junction.applyWeaklyToNonlinearSystem(U, self.U_old, r, J)

    return (r, J)

  # computes the full steady-state residual and Jacobian (strong BC applied)
  def assembleSteadyStateSystem(self, U):
    r, J = self.assembleSteadyStateSystemWithoutStrongBC(U)
    self.applyStrongConstraintsToNonlinearSystem(U, r, J)
    return (r, J)

  ## Applies strong constraints to a nonlinear system solved with Newton's method
  # @param[in] U  implicit solution vector
  # @param[in] r  nonlinear system residual vector
  # @param[in] J  nonlinear system Jacobian matrix
  def applyStrongConstraintsToNonlinearSystem(self, U, r, J):
    # BCs
    for bc in self.bcs:
      bc.applyStrongBCNonlinearSystem(U, r, J)

    # junctions
    for junction in self.junctions:
      junction.applyStronglyToNonlinearSystem(U, self.U_old, r, J)

  ## Applies strong constraints to a linear system matrix.
  #
  # This is separated from the corresponding RHS vector modification function
  # because the matrix needs to be modified only once; the RHS vector might
  # depend on time or the solution vector.
  #
  # @param[in] A      linear system matrix
  def applyStrongConstraintsToLinearSystemMatrix(self, A):
    # BCs
    for bc in self.bcs:
      bc.applyStrongBCLinearSystemMatrix(A)

    # junctions
    for junction in self.junctions:
      junction.applyStronglyToLinearSystemMatrix(A)

  ## Applies strong constraints to a linear system RHS vector.
  # @param[in] U_old  old solution, needed if Dirichlet values are solution-dependent
  # @param[in] b      linear system RHS vector
  def applyStrongConstraintsToLinearSystemRHSVector(self, U_old, b):
    # BCs
    for bc in self.bcs:
      bc.applyStrongBCLinearSystemRHSVector(U_old, b)

    # junctions
    for junction in self.junctions:
      junction.applyStronglyToLinearSystemRHSVector(U_old, b)

  ## Computes the steady-state residual and Jacobian
  def addSteadyStateSystem(self, U, r, J):
    data = dict()
    der = self.dof_handler.initializeDerivativeData(self.aux_names)

    data["phi"] = self.fe_values.get_phi()
    data["g"] = self.gravity

    for elem in xrange(self.dof_handler.n_cell):
      r_cell = np.zeros(self.dof_handler.n_dof_per_cell)
      J_cell = np.zeros(shape=(self.dof_handler.n_dof_per_cell, self.dof_handler.n_dof_per_cell))

      data["grad_phi"] = self.fe_values.get_grad_phi(elem)
      data["JxW"] = self.fe_values.get_JxW(elem)
      data["dx"] = self.dof_handler.h[elem]

      # compute solution
      self.computeLocalCellSolution(U, elem, data)

      # compute auxiliary quantities
      for aux in self.aux_list:
        aux.compute(data, der)

      # compute the local residual and Jacobian
      for kernel in self.fem_kernels:
        kernel.apply(data, der, r_cell, J_cell)

      # aggregate cell residual and matrix into global residual and matrix
      self.dof_handler.aggregateLocalCellVector(r, r_cell, elem)
      self.dof_handler.aggregateLocalCellMatrix(J, J_cell, elem)

  ## Computes the local cell solution and gradients for 1-phase flow
  def computeLocalCellSolutionOnePhase(self, U, elem, data):
    data["vf1"] = self.fe_values.computeLocalVolumeFractionSolution(U, elem)
    data["arho1"] = self.fe_values.computeLocalSolution(U, VariableName.ARho, 0, elem)
    data["arhou1"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoU, 0, elem)
    data["arhoE1"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoE, 0, elem)
    if self.need_solution_gradients:
      data["grad_vf1"] = self.fe_values.computeLocalVolumeFractionSolutionGradient(U, elem)
      data["grad_arho1"] = self.fe_values.computeLocalSolutionGradient(U, VariableName.ARho, 0, elem)
      data["grad_arhou1"] = self.fe_values.computeLocalSolutionGradient(U, VariableName.ARhoU, 0, elem)
      data["grad_arhoE1"] = self.fe_values.computeLocalSolutionGradient(U, VariableName.ARhoE, 0, elem)

  ## Computes the local node solution for 1-phase flow
  def computeLocalNodeSolutionOnePhase(self, U, k, data):
    arho1_index = self.dof_handler.arho_index[0]
    arhou1_index = self.dof_handler.arhou_index[0]
    arhoE1_index = self.dof_handler.arhoE_index[0]
    data["vf1"] = self.dof_handler.getVolumeFraction(U, k)
    data["arho1"] = U[self.dof_handler.i(k, arho1_index)]
    data["arhou1"] = U[self.dof_handler.i(k, arhou1_index)]
    data["arhoE1"] = U[self.dof_handler.i(k, arhoE1_index)]

  ## Computes the local cell solution and gradients for 2-phase flow
  def computeLocalCellSolutionTwoPhase(self, U, elem, data):
    data["vf1"] = self.fe_values.computeLocalVolumeFractionSolution(U, elem)
    data["arho1"] = self.fe_values.computeLocalSolution(U, VariableName.ARho, 0, elem)
    data["arhou1"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoU, 0, elem)
    data["arhoE1"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoE, 0, elem)
    data["arho2"] = self.fe_values.computeLocalSolution(U, VariableName.ARho, 1, elem)
    data["arhou2"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoU, 1, elem)
    data["arhoE2"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoE, 1, elem)
    data["grad_vf1"] = self.fe_values.computeLocalVolumeFractionSolutionGradient(U, elem)
    if self.need_solution_gradients:
      data["grad_arho1"] = self.fe_values.computeLocalSolutionGradient(U, VariableName.ARho, 0, elem)
      data["grad_arhou1"] = self.fe_values.computeLocalSolutionGradient(U, VariableName.ARhoU, 0, elem)
      data["grad_arhoE1"] = self.fe_values.computeLocalSolutionGradient(U, VariableName.ARhoE, 0, elem)
      data["grad_arho2"] = self.fe_values.computeLocalSolutionGradient(U, VariableName.ARho, 1, elem)
      data["grad_arhou2"] = self.fe_values.computeLocalSolutionGradient(U, VariableName.ARhoU, 1, elem)
      data["grad_arhoE2"] = self.fe_values.computeLocalSolutionGradient(U, VariableName.ARhoE, 1, elem)

  ## Computes the local node solution for 2-phase flow
  def computeLocalNodeSolutionTwoPhase(self, U, k, data):
    arho1_index = self.dof_handler.arho_index[0]
    arhou1_index = self.dof_handler.arhou_index[0]
    arhoE1_index = self.dof_handler.arhoE_index[0]
    arho2_index = self.dof_handler.arho_index[1]
    arhou2_index = self.dof_handler.arhou_index[1]
    arhoE2_index = self.dof_handler.arhoE_index[1]
    data["vf1"] = self.dof_handler.getVolumeFraction(U, k)
    data["arho1"] = U[self.dof_handler.i(k, arho1_index)]
    data["arhou1"] = U[self.dof_handler.i(k, arhou1_index)]
    data["arhoE1"] = U[self.dof_handler.i(k, arhoE1_index)]
    data["arho2"] = U[self.dof_handler.i(k, arho2_index)]
    data["arhou2"] = U[self.dof_handler.i(k, arhou2_index)]
    data["arhoE2"] = U[self.dof_handler.i(k, arhoE2_index)]

  def solve(self):
    self.nonlinear_solver.solve(self.U)
