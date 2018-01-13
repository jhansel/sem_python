import numpy as np

from ..base.enums import ModelType, VariableName
from ..closures.thermodynamic_functions import computeVelocity, \
  computeSpecificVolume, computeSpecificTotalEnergy, \
  computeSpecificInternalEnergy
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
    self.registerFloatListParameter("gravity", "3-D gravitational acceleration vector")
    self.registerParameter("ht_data", "List of HeatTransferData objects")
    self.registerParameter("dof_handler", "Degree of freedom handler")
    self.registerParameter("meshes", "List of meshes")
    self.registerParameter("nonlinear_solver_params", "Nonlinear solver parameters")
    self.registerParameter("stabilization", "Stabilization")
    self.registerParameter("factory", "Factory")
    self.registerBoolParameter("split_source", "Use source-splitting?", False)
    self.registerBoolParameter("group_fem", "Use group FEM?", False)
    self.registerBoolParameter("verbose", "Print execution information?", True)

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
    self.ht_data = params.get("ht_data")
    self.dof_handler = params.get("dof_handler")
    self.meshes = params.get("meshes")
    self.nonlinear_solver_params = params.get("nonlinear_solver_params")
    self.factory = params.get("factory")
    stabilization = params.get("stabilization")
    self.split_source = params.get("split_source")
    self.group_fem = params.get("group_fem")
    self.verbose = params.get("verbose")
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

    # initialize the junction constraint variables
    for junction in self.junctions:
      junction.initializeConstraintVariables(self.U)

    # set local solution update function
    if self.model_type == ModelType.OnePhase:
      self.computeLocalCellSolution = self.computeLocalCellSolutionOnePhase
      self.computeLocalNodeSolution = self.computeLocalNodeSolutionOnePhase
      self.getLocalNodalSolution = self.getLocalNodalSolutionOnePhase
    else:
      self.computeLocalCellSolution = self.computeLocalCellSolutionTwoPhase
      self.computeLocalNodeSolution = self.computeLocalNodeSolutionTwoPhase
      self.getLocalNodalSolution = self.getLocalNodalSolutionTwoPhase

    # create aux quantities
    self.aux_list = self.createIndependentPhaseAuxQuantities(0)
    if self.model_type != ModelType.OnePhase:
      self.aux_list += self.createIndependentPhaseAuxQuantities(1)
    if self.model_type == ModelType.TwoPhase:
      self.aux_list += self.createPhaseInteractionAuxQuantities() + interface_closures.createAuxQuantities()
    self.aux_list += stabilization.createAuxQuantities()

    # get list of aux quantities
    self.aux_names = [aux.name for aux in self.aux_list]

    # create nodal aux quantities if needed by group FEM
    if self.group_fem:
      self.nodal_aux_list = self.createNodalAuxQuantities(0)
      if self.model_type != ModelType.OnePhase:
        self.nodal_aux_list = self.createNodalAuxQuantities(1)

      self.group_fem_interp_aux_list = self.createInterpolatedAux(0)
      if self.model_type != ModelType.OnePhase:
        self.group_fem_interp_aux_list = self.createInterpolatedAux(1)
    else:
      self.nodal_aux_list = list()

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
    eos_phase = self.eos_list[phase]
    arhoA_index = self.dof_handler.variable_index[VariableName.ARhoA][phase]
    arhouA_index = self.dof_handler.variable_index[VariableName.ARhoUA][phase]
    arhoEA_index = self.dof_handler.variable_index[VariableName.ARhoEA][phase]

    for ic in ics:
      # get corresponding mesh
      mesh_name = ic.mesh_name
      i_mesh = self.dof_handler.mesh_name_to_mesh_index[mesh_name]
      mesh = self.meshes[i_mesh]

      # get appropriate volume fraction function
      if self.model_type == ModelType.OnePhase:
        def initial_vf(x):
          return 1
      else:
        initial_vf1 = ic.vf1
        if phase == 0:
          initial_vf = initial_vf1
        else:
          def initial_vf(x):
            return 1 - initial_vf1(x)

      # get relevant IC functions
      A_function = ic.A
      initial_p = ic.p[phase]
      initial_u = ic.u[phase]
      if ic.specified_rho:
        initial_rho = ic.rho[phase]
      else:
        initial_T = ic.T[phase]

      # compute IC
      for k_mesh in range(mesh.n_node):
        k = self.dof_handler.k_from_k_mesh(k_mesh, i_mesh)

        x = self.dof_handler.x[k]

        A = A_function(x)
        vf = initial_vf(x)
        p = initial_p(x)
        u = initial_u(x)
        if ic.specified_rho:
          rho = initial_rho(x)
        else:
          T = initial_T(x)
          rho, _, _ = eos_phase.rho(p, T)
        e = eos_phase.e(1.0 / rho, p)[0]
        E = e + 0.5 * u * u
        self.U[self.dof_handler.i(k, arhoA_index)] = vf * rho * A
        self.U[self.dof_handler.i(k, arhouA_index)] = vf * rho * u * A
        self.U[self.dof_handler.i(k, arhoEA_index)] = vf * rho * E * A

  def initializeVolumeFractionSolution(self, ics):
    aA1_index = self.dof_handler.variable_index[VariableName.AA1][0]

    for ic in ics:
      # get corresponding mesh
      mesh_name = ic.mesh_name
      i_mesh = self.dof_handler.mesh_name_to_mesh_index[mesh_name]
      mesh = self.meshes[i_mesh]

      for k_mesh in range(mesh.n_node):
        k = self.dof_handler.k_from_k_mesh(k_mesh, i_mesh)
        x = self.dof_handler.x[k]
        self.U[self.dof_handler.i(k, aA1_index)] = ic.vf1(x) * ic.A(x)

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

  def createPhaseInteractionAuxQuantities(self):
    aux_list = list()
    params = {"aux": "vf1", "variable_names": ["aA1", "A"]}
    aux_list.append(self.factory.createObject("AuxGradient", params))

    return aux_list

  def createNodalAuxQuantities(self, phase):
    # create list of aux quantities to create
    aux_names = list()
    if phase == 0:
      if self.model_type == ModelType.OnePhase:
        aux_names.append("VolumeFraction1Phase")
      else:
        aux_names.append("VolumeFractionPhase1")
    else:
      aux_names.append("VolumeFractionPhase2")
    aux_names += ["Velocity", "Density", "SpecificVolume", "SpecificTotalEnergy", "SpecificInternalEnergy", "Pressure"]

    # create the aux quantities for this phase
    aux_list = list()
    for aux_name in aux_names_phase:
      params = {"phase": phase}
      if aux_name == "Pressure":
        params["p_function"] = self.eos_list[phase].p
      aux_list.append(self.factory.createObject(aux_name, params))

    return aux_list

  def createInterpolatedAux(self, phase):
    aux_list = list()

    phase_str = str(phase + 1)
    arhoA = "arhoA" + phase_str
    arhouA = "arhouA" + phase_str
    arhoEA = "arhoEA" + phase_str

    # mass
    var = "inviscflux_arhoA" + phase_str
    params = {"phase": phase, "variable": var, "dependencies": [arhouA]}
    aux_list.append(self.factory.createObject(var, params))

    # momentum
    var = "inviscflux_arhouA" + phase_str
    params = {"phase": phase, "variable": var, "dependencies": ["aA1", arhoA, arhouA, arhoEA]}
    aux_list.append(self.factory.createObject(var, params))

    # energy
    var = "inviscflux_arhoEA" + phase_str
    params = {"phase": phase, "variable": var, "dependencies": ["aA1", arhoA, arhouA, arhoEA]}
    aux_list.append(self.factory.createObject(var, params))

    return aux_list

  def createIndependentPhaseAdvectionKernels(self, phase):
    kernels = list()
    params = {"phase": phase, "dof_handler": self.dof_handler}
    kernel_name_list = ["MassAdvection", "MomentumAdvection", "MomentumAreaGradient", "EnergyAdvection"]
    for kernel_name in kernel_name_list:
      kernels.append(self.factory.createObject(kernel_name, params))
    return kernels

  def createIndependentPhaseSourceKernels(self, phase):
    params = {"phase": phase, "dof_handler": self.dof_handler, "is_nodal": self.split_source}
    kernel_names = ["MomentumGravity", "EnergyGravity", "EnergyHeatTransfer"]
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
  def assembleSteadyStateSystemWithoutStrongConstraints(self, U):
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
    nodal_data = dict()
    nodal_der = self.dof_handler.initializeDerivativeData(self.aux_names)

    data["phi"] = self.fe_values.get_phi()
    for elem in range(self.dof_handler.n_cell):
      r_cell = np.zeros(self.dof_handler.n_dof_per_cell)
      J_cell = np.zeros(shape=(self.dof_handler.n_dof_per_cell, self.dof_handler.n_dof_per_cell))

      i_mesh = self.dof_handler.elem_to_mesh_index[elem]

      data["grad_phi"] = self.fe_values.get_grad_phi(elem)
      data["JxW"] = self.fe_values.get_JxW(elem)
      data["dx"] = self.dof_handler.h[elem]
      data["g"] = np.dot(self.meshes[i_mesh].orientation, self.gravity)
      data["T_wall"] = self.ht_data[i_mesh].T_wall
      data["htc_wall"] = self.ht_data[i_mesh].htc_wall
      data["P_heat"] = self.ht_data[i_mesh].P_heat

      # group FEM
      if self.group_fem:
        # compute nodal solution on element
        self.getLocalNodalSolution(U, elem, nodal_data)

        # compute nodal aux on element
        # NOTE: It would be approximately half as expensive (but more memory intensive)
        # to compute the nodal auxes as a global vector.
        for aux in self.nodal_aux_list:
          aux.compute(nodal_data, nodal_der)

        # compute flux auxes for group fem
        for aux in self.group_fem_interp_aux_list:
          aux.compute(nodal_data, nodal_der, data, der)

      # compute elemental solution
      self.computeLocalCellSolution(U, elem, data)

      # compute elemental auxiliary quantities
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
    data["A"] = self.fe_values.computeLocalArea(elem)
    data["aA1"] = self.fe_values.computeLocalVolumeFractionSolution(U, elem)
    data["arhoA1"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoA, 0, elem)
    data["arhouA1"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoUA, 0, elem)
    data["arhoEA1"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoEA, 0, elem)
    data["grad_A"] = self.fe_values.computeLocalAreaGradient(elem)
    if self.need_solution_gradients:
      data["grad_aA1"] = self.fe_values.computeLocalVolumeFractionSolutionGradient(U, elem)
      data["grad_arhoA1"] = self.fe_values.computeLocalSolutionGradient(U, VariableName.ARhoA, 0, elem)
      data["grad_arhouA1"] = self.fe_values.computeLocalSolutionGradient(U, VariableName.ARhoUA, 0, elem)
      data["grad_arhoEA1"] = self.fe_values.computeLocalSolutionGradient(U, VariableName.ARhoEA, 0, elem)

  ## Computes the local node solution for 1-phase flow
  def computeLocalNodeSolutionOnePhase(self, U, k, data):
    arhoA1_index = self.dof_handler.arhoA_index[0]
    arhouA1_index = self.dof_handler.arhouA_index[0]
    arhoEA1_index = self.dof_handler.arhoEA_index[0]
    data["A"] = self.dof_handler.A[k]
    data["aA1"] = self.dof_handler.aA1(U, k)
    data["arhoA1"] = U[self.dof_handler.i(k, arhoA1_index)]
    data["arhouA1"] = U[self.dof_handler.i(k, arhouA1_index)]
    data["arhoEA1"] = U[self.dof_handler.i(k, arhoEA1_index)]

  ## Computes the local cell solution and gradients for 2-phase flow
  def computeLocalCellSolutionTwoPhase(self, U, elem, data):
    data["A"] = self.fe_values.computeLocalArea(elem)
    data["aA1"] = self.fe_values.computeLocalVolumeFractionSolution(U, elem)
    data["arhoA1"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoA, 0, elem)
    data["arhouA1"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoUA, 0, elem)
    data["arhoEA1"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoEA, 0, elem)
    data["arhoA2"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoA, 1, elem)
    data["arhouA2"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoUA, 1, elem)
    data["arhoEA2"] = self.fe_values.computeLocalSolution(U, VariableName.ARhoEA, 1, elem)
    data["grad_A"] = self.fe_values.computeLocalAreaGradient(elem)
    data["grad_aA1"] = self.fe_values.computeLocalVolumeFractionSolutionGradient(U, elem)
    if self.need_solution_gradients:
      data["grad_arhoA1"] = self.fe_values.computeLocalSolutionGradient(U, VariableName.ARhoA, 0, elem)
      data["grad_arhouA1"] = self.fe_values.computeLocalSolutionGradient(U, VariableName.ARhoUA, 0, elem)
      data["grad_arhoEA1"] = self.fe_values.computeLocalSolutionGradient(U, VariableName.ARhoEA, 0, elem)
      data["grad_arhoA2"] = self.fe_values.computeLocalSolutionGradient(U, VariableName.ARhoA, 1, elem)
      data["grad_arhouA2"] = self.fe_values.computeLocalSolutionGradient(U, VariableName.ARhoUA, 1, elem)
      data["grad_arhoEA2"] = self.fe_values.computeLocalSolutionGradient(U, VariableName.ARhoEA, 1, elem)

  ## Computes the local node solution for 2-phase flow
  def computeLocalNodeSolutionTwoPhase(self, U, k, data):
    arhoA1_index = self.dof_handler.arhoA_index[0]
    arhouA1_index = self.dof_handler.arhouA_index[0]
    arhoEA1_index = self.dof_handler.arhoEA_index[0]
    arhoA2_index = self.dof_handler.arhoA_index[1]
    arhouA2_index = self.dof_handler.arhouA_index[1]
    arhoEA2_index = self.dof_handler.arhoEA_index[1]
    data["A"] = self.dof_handler.A[k]
    data["aA1"] = self.dof_handler.aA1(U, k)
    data["arhoA1"] = U[self.dof_handler.i(k, arhoA1_index)]
    data["arhouA1"] = U[self.dof_handler.i(k, arhouA1_index)]
    data["arhoEA1"] = U[self.dof_handler.i(k, arhoEA1_index)]
    data["arhoA2"] = U[self.dof_handler.i(k, arhoA2_index)]
    data["arhouA2"] = U[self.dof_handler.i(k, arhouA2_index)]
    data["arhoEA2"] = U[self.dof_handler.i(k, arhoEA2_index)]

  def getLocalNodalSolutionOnePhase(self, U, elem, data):
    data["A"] = self.fe_values.getLocalNodalArea(elem)
    data["arhoA1"] = self.fe_values.getLocalNodalSolution(U, VariableName.ARhoA, 0, elem)
    data["arhouA1"] = self.fe_values.getLocalNodalSolution(U, VariableName.ARhoUA, 0, elem)
    data["arhoEA1"] = self.fe_values.getLocalNodalSolution(U, VariableName.ARhoEA, 0, elem)

  def getLocalNodalSolutionTwoPhase(self, U, elem, data):
    data["A"] = self.fe_values.getLocalNodalArea(elem)
    data["aA1"] = self.fe_values.getLocalNodalVolumeFractionSolution(U, elem)
    data["arhoA1"] = self.fe_values.getLocalNodalSolution(U, VariableName.ARhoA, 0, elem)
    data["arhouA1"] = self.fe_values.getLocalNodalSolution(U, VariableName.ARhoUA, 0, elem)
    data["arhoEA1"] = self.fe_values.getLocalNodalSolution(U, VariableName.ARhoEA, 0, elem)
    data["arhoA2"] = self.fe_values.getLocalNodalSolution(U, VariableName.ARhoA, 1, elem)
    data["arhouA2"] = self.fe_values.getLocalNodalSolution(U, VariableName.ARhoUA, 1, elem)
    data["arhoEA2"] = self.fe_values.getLocalNodalSolution(U, VariableName.ARhoEA, 1, elem)

  def solve(self):
    self.nonlinear_solver.solve(self.U)
