import os
import sys
base_dir = os.environ["SEM_PYTHON_DIR"]

sys.path.append(base_dir + "src/base")
from enums import ModelType, PhaseType

sys.path.append(base_dir + "src/closures")
from thermodynamic_functions import computeDensity, computeVelocity, computeSpecificVolume, \
  computeSpecificTotalEnergy, computeSpecificInternalEnergy

sys.path.append(base_dir + "src/input")
from Parameters import Parameters

sys.path.append(base_dir + "src/output")
from Plotter import Plotter

sys.path.append(base_dir + "src/utilities")
from file_utilities import writeCSVFile

class PostprocessorParameters(Parameters):
  def __init__(self):
    Parameters.__init__(self)
    self.registerBoolParameter("print_solution", "Option to print solution to console", False)
    self.registerBoolParameter("save_solution", "Option to save solution to a file", False)
    self.registerStringParameter("solution_file", "Name of solution output file", "solution.csv")
    self.registerIntParameter("output_precision", "Precision used in solution output file", 10)
    self.registerBoolParameter("plot_solution", "Option to plot solution", False)
    self.registerStringParameter("plot_file", "Name of plot file", "solution.pdf")

class Postprocessor(object):
  def __init__(self, params, model_type, eos_map, dof_handler, mesh):
    self.print_solution = params.get("print_solution")
    self.save_solution = params.get("save_solution")
    self.solution_file = params.get("solution_file")
    self.output_precision = params.get("output_precision")
    self.plot_solution = params.get("plot_solution")
    self.plot_file = params.get("plot_file")
    self.model_type = model_type
    self.eos_map = eos_map
    self.dof_handler = dof_handler
    self.mesh = mesh

  def run(self, U):
    vf1, arho1, arhou1, arhoE1 = self.dof_handler.getPhaseSolution(U, PhaseType.First)
    if (self.model_type != ModelType.OnePhase):
      vf2, arho2, arhou2, arhoE2 = self.dof_handler.getPhaseSolution(U, PhaseType.Second)

    # compute aux quantities
    n = self.dof_handler.n_node
    rho1 = computeDensity(vf1, arho1)[0]
    u1 = computeVelocity(arho1, arhou1)[0]
    v1 = computeSpecificVolume(rho1)[0]
    E1 = computeSpecificTotalEnergy(arho1, arhoE1)[0]
    e1 = computeSpecificInternalEnergy(u1, E1)[0]
    eos1 = self.eos_map[PhaseType.First]
    p1 = eos1.p(v1, e1)[0]
    if (self.model_type != ModelType.OnePhase):
      rho2 = [computeDensity(vf2[k], arho2[k])[0] for k in xrange(n)]
      u2 = [computeVelocity(arho2[k], arhou2[k])[0] for k in xrange(n)]
      v2 = computeSpecificVolume(rho2)[0]
      E2 = computeSpecificTotalEnergy(arho2, arhoE2)[0]
      e2 = computeSpecificInternalEnergy(u2, E2)[0]
      eos2 = self.eos_map[PhaseType.Second]
      p2 = eos2.p(v2, e2)[0]

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
      data = dict()
      data["x"] = self.mesh.x
      if self.model_type == ModelType.OnePhase:
        data["rho"] = rho1
        data["u"] = u1
        data["p"] = p1
      else:
        data["rho1"] = rho1
        data["u1"] = u1
        data["p1"] = p1
        data["rho2"] = rho2
        data["u2"] = u2
        data["p2"] = p2
        if self.model_type == ModelType.TwoPhase:
          data["vf1"] = vf1

      # write to file
      writeCSVFile(data, self.solution_file, self.output_precision)

    # plot solution
    if (self.plot_solution):
      x_label = "Position, $x$"
      if (self.model_type == ModelType.OnePhase):
        plotter = Plotter("", "Density, $\\rho$", 3)
        plotter.setXRange(self.mesh.x_min, self.mesh.x_max)
        plotter.addSet(self.mesh.x, rho1, "$\\rho$")
        plotter.fixNearConstantPlot()
        plotter.nextSubplot("", "Velocity, $u$")
        plotter.addSet(self.mesh.x, u1, "$u$")
        plotter.fixNearConstantPlot()
        plotter.nextSubplot(x_label, "Pressure, $p$ [kPa]")
        plotter.addSet(self.mesh.x, p1, "$p$", scale=1e-3)
        plotter.fixNearConstantPlot()
      elif (self.model_type == ModelType.TwoPhaseNonInteracting):
        plotter = Plotter(x_label, "Density, $\\rho$", (2,2))
        plotter.setXRange(self.mesh.x_min, self.mesh.x_max)
        plotter.addSet(self.mesh.x, rho1, "$\\rho_1$")
        plotter.fixNearConstantPlot()
        plotter.nextSubplot(x_label, "Density, $\\rho$")
        plotter.addSet(self.mesh.x, rho2, "$\\rho_2$")
        plotter.fixNearConstantPlot()
        plotter.nextSubplot(x_label, "Velocity, $u$")
        plotter.addSet(self.mesh.x, u1, "$u_1$")
        plotter.addSet(self.mesh.x, u2, "$u_2$")
        plotter.fixNearConstantPlot()
        plotter.nextSubplot(x_label, "Pressure, $p$ [kPa]")
        plotter.addSet(self.mesh.x, p1, "$p_1$", scale=1e-3)
        plotter.addSet(self.mesh.x, p2, "$p_2$", scale=1e-3)
        plotter.fixNearConstantPlot()
      elif (self.model_type == ModelType.TwoPhase):
        plotter = Plotter(x_label, "Volume Fraction, $\\alpha$", (2,2))
        plotter.setXRange(self.mesh.x_min, self.mesh.x_max)
        plotter.addSet(self.mesh.x, vf1, "$\\alpha_1$")
        plotter.fixNearConstantPlot()
        plotter.nextSubplot(x_label, "Density, $\\rho$")
        plotter.addSet(self.mesh.x, rho1, "$\\rho_1$")
        plotter.addSet(self.mesh.x, rho2, "$\\rho_2$")
        plotter.fixNearConstantPlot()
        plotter.nextSubplot(x_label, "Velocity, $u$")
        plotter.addSet(self.mesh.x, u1, "$u_1$")
        plotter.addSet(self.mesh.x, u2, "$u_2$")
        plotter.fixNearConstantPlot()
        plotter.nextSubplot(x_label, "Pressure, $p$ [kPa]")
        plotter.addSet(self.mesh.x, p1, "$p_1$", scale=1e-3)
        plotter.addSet(self.mesh.x, p2, "$p_2$", scale=1e-3)
        plotter.fixNearConstantPlot()

      # save plot
      plotter.save(self.plot_file)
