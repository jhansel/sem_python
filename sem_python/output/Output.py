from collections import OrderedDict

from enums import ModelType
from thermodynamic_functions import computeDensity, computeVelocity, computeSpecificVolume, \
  computeSpecificTotalEnergy, computeSpecificInternalEnergy
from Parameters import Parameters
from Plotter import Plotter
from file_utilities import writeCSVFile

class OutputParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)
    self.registerBoolParameter("print_solution", "Option to print solution to console", False)
    self.registerBoolParameter("save_solution", "Option to save solution to a file", False)
    self.registerStringParameter("solution_file", "Name of solution output file", "solution.csv")
    self.registerIntParameter("output_precision", "Precision used in solution output file", 5)
    self.registerBoolParameter("plot_solution", "Option to plot solution", False)
    self.registerStringParameter("plot_file", "Name of plot file", "solution.pdf")
    self.registerBoolParameter("plot_area", "Option to plot area instead of temperature", False)

    self.registerParameter("model", "Model")
    self.registerParameter("eos_list", "List of equations of state")
    self.registerParameter("dof_handler", "Degree of freedom handler")
    self.registerParameter("meshes", "List of meshes")

class Output(object):
  def __init__(self, params):
    self.print_solution = params.get("print_solution")
    self.save_solution = params.get("save_solution")
    self.solution_file = params.get("solution_file")
    self.output_precision = params.get("output_precision")
    self.plot_solution = params.get("plot_solution")
    self.plot_file = params.get("plot_file")
    self.plot_area = params.get("plot_area")

    model = params.get("model")
    self.model_type = model.model_type
    self.eos_list = params.get("eos_list")
    self.dof_handler = params.get("dof_handler")
    self.meshes = params.get("meshes")

  def run(self, U):
    vf1, arhoA1, arhouA1, arhoEA1 = self.dof_handler.getPhaseSolution(U, 0)
    if (self.model_type != ModelType.OnePhase):
      vf2, arhoA2, arhouA2, arhoEA2 = self.dof_handler.getPhaseSolution(U, 1)

    # compute aux quantities
    n = self.dof_handler.n_node
    rho1 = computeDensity(vf1, arhoA1, self.dof_handler.A)[0]
    u1 = computeVelocity(arhoA1, arhouA1)[0]
    v1 = computeSpecificVolume(rho1)[0]
    E1 = computeSpecificTotalEnergy(arhoA1, arhoEA1)[0]
    e1 = computeSpecificInternalEnergy(u1, E1)[0]
    eos1 = self.eos_list[0]
    p1 = eos1.p(v1, e1)[0]
    T1 = eos1.T(v1, e1)[0]
    if (self.model_type != ModelType.OnePhase):
      rho2 = computeDensity(vf2, arhoA2, self.dof_handler.A)[0]
      u2 = computeVelocity(arhoA2, arhouA2)[0]
      v2 = computeSpecificVolume(rho2)[0]
      E2 = computeSpecificTotalEnergy(arhoA2, arhoEA2)[0]
      e2 = computeSpecificInternalEnergy(u2, E2)[0]
      eos2 = self.eos_list[1]
      p2 = eos2.p(v2, e2)[0]
      T2 = eos2.T(v2, e2)[0]

    # print solution
    if (self.print_solution):
      print
      if (self.model_type == ModelType.OnePhase):
        print ("%15s%15s%15s") % ("rho", "u", "p")
        for k in xrange(n):
          print ("%15g%15g%15g") % (rho1[k], u1[k], p1[k])
      elif (self.model_type == ModelType.TwoPhaseNonInteracting):
        print ("%15s%15s%15s%15s%15s%15s") % ("rho1", "u1", "p1", "rho2", "u2", "p2")
        for k in xrange(n):
          print ("%15g%15g%15g%15g%15g%15g") % (rho1[k], u1[k], p1[k], rho2[k], u2[k], p2[k])
      elif (self.model_type == ModelType.TwoPhase):
        print ("%15s%15s%15s%15s%15s%15s%15s") % ("vf1", "rho1", "u1", "p1", "rho2", "u2", "p2")
        for k in xrange(n):
          print ("%15g%15g%15g%15g%15g%15g%15g") % (vf1[k], rho1[k], u1[k], p1[k], rho2[k], u2[k], p2[k])

    # save solution to file
    if self.save_solution:
      # create the data dictionary
      data = OrderedDict()
      data["x"] = self.dof_handler.x
      if self.model_type == ModelType.OnePhase:
        data["rho"] = rho1
        data["u"] = u1
        data["p"] = p1
      else:
        if self.model_type == ModelType.TwoPhase:
          data["vf1"] = vf1
        data["rho1"] = rho1
        data["u1"] = u1
        data["p1"] = p1
        data["rho2"] = rho2
        data["u2"] = u2
        data["p2"] = p2

      # write to file
      writeCSVFile(data, self.solution_file, self.output_precision)

    # plot solution
    if (self.plot_solution):
      x_label = "Position, $x$"
      x = self.dof_handler.x
      x_min = min(x)
      x_max = max(x)

      def addSets(plotter, x, y, y_label, scale=1):
        for i_mesh in xrange(self.dof_handler.n_meshes):
          modified_y_label = y_label + ", Mesh " + str(i_mesh+1)
          plotter.addSet(x[i_mesh], y[i_mesh], modified_y_label, color=i_mesh+1, marker=i_mesh+1, scale=scale)

      if (self.model_type == ModelType.OnePhase):
        x_by_mesh = self.dof_handler.separateNodalQuantityByMesh(x)
        A_by_mesh = self.dof_handler.separateNodalQuantityByMesh(self.dof_handler.A)
        rho1_by_mesh = self.dof_handler.separateNodalQuantityByMesh(rho1)
        u1_by_mesh = self.dof_handler.separateNodalQuantityByMesh(u1)
        p1_by_mesh = self.dof_handler.separateNodalQuantityByMesh(p1)
        T1_by_mesh = self.dof_handler.separateNodalQuantityByMesh(T1)

        plotter = Plotter("", "Density, $\\rho$", (2,2))
        plotter.setXRange(x_min, x_max)
        addSets(plotter, x_by_mesh, rho1_by_mesh, "$\\rho$")
        plotter.fixNearConstantPlot()
        plotter.nextSubplot("", "Velocity, $u$")
        addSets(plotter, x_by_mesh, u1_by_mesh, "$u$")
        plotter.fixNearConstantPlot()
        plotter.nextSubplot(x_label, "Pressure, $p$ [kPa]")
        addSets(plotter, x_by_mesh, p1_by_mesh, "$p$", 1e-3)
        plotter.fixNearConstantPlot()
        if self.plot_area:
          plotter.nextSubplot(x_label, "Area, $A$ [m$^2$]")
          addSets(plotter, x_by_mesh, A_by_mesh, "$A$")
        else:
          plotter.nextSubplot(x_label, "Temperature, $T$ [K]")
          addSets(plotter, x_by_mesh, T1_by_mesh, "$T$")
        plotter.fixNearConstantPlot()
      elif (self.model_type == ModelType.TwoPhaseNonInteracting):
        plotter = Plotter(x_label, "Density, $\\rho$", (2,2))
        plotter.setXRange(x_min, x_max)
        plotter.addSet(x, rho1, "$\\rho_1$")
        plotter.fixNearConstantPlot()
        plotter.nextSubplot(x_label, "Density, $\\rho$")
        plotter.addSet(x, rho2, "$\\rho_2$")
        plotter.fixNearConstantPlot()
        plotter.nextSubplot(x_label, "Velocity, $u$")
        plotter.addSet(x, u1, "$u_1$")
        plotter.addSet(x, u2, "$u_2$")
        plotter.fixNearConstantPlot()
        plotter.nextSubplot(x_label, "Pressure, $p$ [kPa]")
        plotter.addSet(x, p1, "$p_1$", scale=1e-3)
        plotter.addSet(x, p2, "$p_2$", scale=1e-3)
        plotter.fixNearConstantPlot()
      elif (self.model_type == ModelType.TwoPhase):
        plotter = Plotter(x_label, "Volume Fraction, $\\alpha$", (2,2))
        plotter.setXRange(x_min, x_max)
        plotter.addSet(x, vf1, "$\\alpha_1$")
        plotter.fixNearConstantPlot()
        plotter.nextSubplot(x_label, "Density, $\\rho$")
        plotter.addSet(x, rho1, "$\\rho_1$")
        plotter.addSet(x, rho2, "$\\rho_2$")
        plotter.fixNearConstantPlot()
        plotter.nextSubplot(x_label, "Velocity, $u$")
        plotter.addSet(x, u1, "$u_1$")
        plotter.addSet(x, u2, "$u_2$")
        plotter.fixNearConstantPlot()
        plotter.nextSubplot(x_label, "Pressure, $p$ [kPa]")
        plotter.addSet(x, p1, "$p_1$", scale=1e-3)
        plotter.addSet(x, p2, "$p_2$", scale=1e-3)
        plotter.fixNearConstantPlot()

      # save plot
      plotter.save(self.plot_file)
