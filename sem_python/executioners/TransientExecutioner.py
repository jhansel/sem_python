from copy import deepcopy
import numpy as np
from termcolor import colored

from ..base.enums import ModelType, VariableName
from .Executioner import Executioner, ExecutionerParameters
from ..utilities.error_utilities import error

class TransientExecutionerParameters(ExecutionerParameters):
  def __init__(self):
    ExecutionerParameters.__init__(self)
    self.registerFloatParameter("dt", "Nominal time step size")
    self.registerFloatParameter("cfl", "CFL number to compute time step size")
    self.registerFloatParameter("end_time", "End time")
    self.registerBoolParameter("lump_mass_matrix", "Lump the mass matrix?", False)
    self.registerFloatParameter("ss_tol", "Tolerance for steady-state check")

class TransientExecutioner(Executioner):
  def __init__(self, params):
    Executioner.__init__(self, params)

    # determine how time step size is determined
    if params.has("dt") and params.has("cfl"):
      error("The parameters 'dt' and 'cfl' cannot both be provided.")
    elif params.has("dt"):
      self.dt_nominal = params.get("dt")
      self.use_cfl_dt = False
    elif params.has("cfl"):
      self.cfl = params.get("cfl")
      self.use_cfl_dt = True
      self.cfl_dt_aux_list = self.createCFLAuxQuantities()
      self.cfl_aux_names = [aux.name for aux in self.cfl_dt_aux_list]
    else:
      error("Either parameter 'dt' or 'cfl' must be provided.")

    self.end_time = params.get("end_time")
    self.lump_mass_matrix = params.get("lump_mass_matrix")

    if params.has("ss_tol"):
      self.check_ss = True
      self.ss_tol = params.get("ss_tol")
    else:
      self.check_ss = False

    # tolerance to prevent small final time steps due to floating point precision error
    self.end_tolerance = 1e-12

    self.U_old = deepcopy(self.U)

    self.M = self.computeMassMatrix()

  def computeMassMatrix(self):
    M = np.zeros(shape=(self.dof_handler.n_dof, self.dof_handler.n_dof))

    self.addMassMatrixPhase(M, 0)
    if (self.model_type != ModelType.OnePhase):
      self.addMassMatrixPhase(M, 1)
    if (self.model_type == ModelType.TwoPhase):
      self.addMassMatrixVolumeFraction(M)

    return M

  def addMassMatrixPhase(self, M, phase):
    phi = self.fe_values.get_phi()

    arhoA_index = self.dof_handler.variable_index[VariableName.ARhoA][phase]
    arhouA_index = self.dof_handler.variable_index[VariableName.ARhoUA][phase]
    arhoEA_index = self.dof_handler.variable_index[VariableName.ARhoEA][phase]

    for e in xrange(self.dof_handler.n_cell):
      M_cell = np.zeros(shape=(self.dof_handler.n_dof_per_cell, self.dof_handler.n_dof_per_cell))

      JxW = self.fe_values.get_JxW(e)
      for q in xrange(self.quadrature.n_q):
        for k_local in xrange(self.dof_handler.n_dof_per_cell_per_var):
          i_arhoA = self.dof_handler.i(k_local, arhoA_index)
          i_arhouA = self.dof_handler.i(k_local, arhouA_index)
          i_arhoEA = self.dof_handler.i(k_local, arhoEA_index)
          for l_local in xrange(self.dof_handler.n_dof_per_cell_per_var):
            if self.lump_mass_matrix:
              M_cell[i_arhoA,i_arhoA] += phi[k_local,q] * phi[l_local,q] * JxW[q]
              M_cell[i_arhouA,i_arhouA] += phi[k_local,q] * phi[l_local,q] * JxW[q]
              M_cell[i_arhoEA,i_arhoEA] += phi[k_local,q] * phi[l_local,q] * JxW[q]
            else:
              j_arhoA = self.dof_handler.i(l_local, arhoA_index)
              j_arhouA = self.dof_handler.i(l_local, arhouA_index)
              j_arhoEA = self.dof_handler.i(l_local, arhoEA_index)

              M_cell[i_arhoA,j_arhoA] += phi[k_local,q] * phi[l_local,q] * JxW[q]
              M_cell[i_arhouA,j_arhouA] += phi[k_local,q] * phi[l_local,q] * JxW[q]
              M_cell[i_arhoEA,j_arhoEA] += phi[k_local,q] * phi[l_local,q] * JxW[q]

      # aggregate cell matrix into global matrix
      self.dof_handler.aggregateLocalCellMatrix(M, M_cell, e)

  def addMassMatrixVolumeFraction(self, M):
    phi = self.fe_values.get_phi()

    aA1_index = self.dof_handler.variable_index[VariableName.AA1][0]

    for e in xrange(self.dof_handler.n_cell):
      M_cell = np.zeros(shape=(self.dof_handler.n_dof_per_cell, self.dof_handler.n_dof_per_cell))

      JxW = self.fe_values.get_JxW(e)
      for q in xrange(self.quadrature.n_q):
        for k_local in xrange(self.dof_handler.n_dof_per_cell_per_var):
          i_aA1 = self.dof_handler.i(k_local, aA1_index)
          for l_local in xrange(self.dof_handler.n_dof_per_cell_per_var):
            if self.lump_mass_matrix:
              M_cell[i_aA1,i_aA1] += phi[k_local,q] * phi[l_local,q] * JxW[q]
            else:
              j_vf1 = self.dof_handler.i(l_local, aA1_index)
              M_cell[i_aA1,j_vf1] += phi[k_local,q] * phi[l_local,q] * JxW[q]

      # aggregate cell matrix into global matrix
      self.dof_handler.aggregateLocalCellMatrix(M, M_cell, e)

  def assembleTransientSystem(self, U):
    M_dU = np.matmul(self.M, U - self.U_old)
    return (M_dU, self.M)

  ## Integrates the source terms when source-splitting is used
  def takeSourceStep(self, U, dt):
    data = dict()
    der = self.dof_handler.initializeDerivativeData(self.aux_names)

    data["phi"] = np.zeros(shape=(self.dof_handler.n_dof_per_cell_per_var, 1))
    data["grad_phi"] = np.zeros(shape=(self.dof_handler.n_dof_per_cell_per_var, 1))
    data["phi"].fill(1.0)
    data["grad_phi"].fill(float("NaN"))
    data["JxW"] = 1.0

    n_var = self.dof_handler.n_var

    # loop over nodes
    for k in xrange(self.dof_handler.n_node):
      r_node = np.zeros(n_var)
      J_node = np.zeros(shape=(n_var, n_var))

      i_mesh = self.dof_handler.node_to_mesh_index[k]

      data["g"] = np.dot(self.meshes[i_mesh].orientation, self.gravity)
      data["T_wall"] = self.ht_data[i_mesh].T_wall
      data["htc_wall"] = self.ht_data[i_mesh].htc_wall
      data["P_heat"] = self.ht_data[i_mesh].P_heat

      # compute solution
      self.computeLocalNodeSolution(U, k, data)

      # compute auxiliary quantities
      for aux in self.aux_list:
        aux.compute(data, der)

      # compute the local residual and Jacobian
      for kernel in self.source_kernels:
        kernel.apply(data, der, r_node, J_node)

      # add integrated source (recall kernels assume LHS)
      self.dof_handler.aggregateLocalNodeVector(U, -dt * r_node, k)

  def createCFLAuxQuantities(self):
    aux_list = list()

    if self.model_type == ModelType.OnePhase:
      vf_classes = ["VolumeFraction1Phase"]
      self.n_phases = 1
    else:
      vf_classes = ["VolumeFractionPhase1", "VolumeFractionPhase2"]
      self.n_phases = 2

    for phase in xrange(self.n_phases):
      names = [vf_classes[phase]] + ["Density", "SpecificVolume", "Velocity", "SpecificTotalEnergy",
        "SpecificInternalEnergy", "Pressure", "SoundSpeed"]
      for name in names:
        if name == "Pressure":
          params = {"phase": phase, "p_function": self.eos_list[phase].p}
        elif name == "SoundSpeed":
          params = {"phase": phase, "c_function": self.eos_list[phase].c}
        else:
          params = {"phase": phase}
        aux_list.append(self.factory.createObject(name, params))

    return aux_list

  def computeTimeStepSizeFromCFL(self):
    # determine minimum cell width
    dx_min = self.meshes[0].getMinimumCellWidth()
    for mesh in self.meshes:
      dx_min = min(dx_min, mesh.getMinimumCellWidth())

    # determine maximum wave speed
    data = dict()
    der = self.dof_handler.initializeDerivativeData(self.cfl_aux_names)

    for phase in xrange(self.n_phases):
      phase_str = str(phase + 1)
      vf = "vf" + phase_str
      arhoA = "arhoA" + phase_str
      arhouA = "arhouA" + phase_str
      arhoEA = "arhoEA" + phase_str
      data[vf], data[arhoA], data[arhouA], data[arhoEA] = self.dof_handler.getPhaseSolution(self.U_old, phase)

    for aux in self.cfl_dt_aux_list:
      aux.compute(data, der)

    wave_speed_max = 0
    for phase in xrange(self.n_phases):
      phase_str = str(phase + 1)
      u = "u" + phase_str
      c = "c" + phase_str
      wave_speed = "wave_speed" + phase_str

      data[wave_speed] = abs(data[u]) + data[c]
      wave_speed_max = max(wave_speed_max, max(data[wave_speed]))

    # compute CFL time step size
    dt_CFL = dx_min / wave_speed_max

    return self.cfl * dt_CFL

  def run(self):
    transient_incomplete = True
    t = 0.0
    time_step = 1
    while (transient_incomplete):
      # compute time step size
      if self.use_cfl_dt:
        dt = self.computeTimeStepSizeFromCFL()
      else:
        dt = self.dt_nominal

      # shorten time step if at end of transient
      if (t + dt + self.end_tolerance >= self.end_time):
        transient_incomplete = False
        self.dt = self.end_time - t
      else:
        self.dt = dt

      # update time
      t += self.dt

      if self.verbose:
        print "\nTime step %i: t = %g, dt = %g" % (time_step, t, self.dt)

      # solve the time step
      self.solve()

      # perform source term integration
      if self.split_source:
        self.takeSourceStep(self.U, self.dt)

      # check for steady-state
      if self.check_ss:
        dU_dt = (self.U - self.U_old) / self.dt
        dU_dt_norm = np.linalg.norm(dU_dt, 2)
        U_old_norm = np.linalg.norm(self.U_old, 2)
        U_change_norm = dU_dt_norm / U_old_norm
        if self.verbose:
          print "Relative solution change: %e" % (U_change_norm)
        if U_change_norm < self.ss_tol:
          if self.verbose:
            print colored("\nConverged to steady-state!\n", "green")
          return self.U

      # store the old solution
      self.U_old = deepcopy(self.U)

      # save old solution and increment time step index
      self.U_old = deepcopy(self.U)
      time_step += 1

    if self.verbose:
      print ""

    return self.U

  def solve(self):
    self.nonlinear_solver.solve(self.U)
