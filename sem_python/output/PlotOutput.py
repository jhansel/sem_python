from .Output import Output, OutputParameters
from .Plotter import Plotter
from ..utilities.error_utilities import error


class PlotOutputParameters(OutputParameters):

    def __init__(self, factory):
        OutputParameters.__init__(self, factory)
        self.registerStringParameter("file_name", "Name of solution output file", "solution.pdf")
        self.registerStringListParameter("plot_sets", "List of data to plot")
        self.registerBoolParameter(
            "plot_by_mesh", "Flag to plot each mesh separately on a plot", False)
        self.registerStringParameter("legend_location", "Location of legend", "upper right")
        self.registerFloatParameter("size_x", "Default x-length of figure in inches", 8)
        self.registerFloatParameter("size_y", "Default y-length of figure in inches", 6)
        self.registerFloatListParameter("y_bounds", "Bounds for y-axis, after scaling")


class PlotOutput(Output):

    def __init__(self, params):
        Output.__init__(self, params)
        self.file_name = params.get("file_name")
        plot_sets_raw = params.get("plot_sets")
        self.plot_by_mesh = params.get("plot_by_mesh")
        self.legend_location = params.get("legend_location")
        self.size_x = params.get("size_x")
        self.size_y = params.get("size_y")
        self.supplied_y_bounds = params.has("y_bounds")
        if self.supplied_y_bounds:
            self.y_bounds = params.get("y_bounds")

        self.n_subplots = len(plot_sets_raw)
        if self.n_subplots < 1:
            error("There must be at least one set specified by the 'plot_sets' parameter")
        elif self.n_subplots > 4:
            error("The maximum number of sub-plots is 4")

        if self.supplied_y_bounds and self.n_subplots > 1:
            error("Y-bounds can only be specified if there is only 1 sub-plot")

        # parse the plot set list in case there are multi-sets for a subplot
        # multi-sets must be of the format "[a,b,c]" (no spaces)
        self.plot_sets = list()
        for plot_set in plot_sets_raw:
            self.plot_sets.append(plot_set.rstrip("]").lstrip("[").split(","))

        # determine the sub-plot layout
        if self.n_subplots == 1:
            self.subplot_layout = (1, 1)
            self.print_x_label = [True]
        elif self.n_subplots == 2:
            self.subplot_layout = (1, 2)
            self.print_x_label = [False, True]
        elif self.n_subplots == 3:
            self.subplot_layout = (1, 3)
            self.print_x_label = [False, False, True]
        else:  # self.n_subplots == 4
            self.subplot_layout = (2, 2)
            self.print_x_label = [False, False, True, True]

        self.base_name_sets = [{"p", "p0"}, {"T", "T0"}]

        self.base_name_to_y_label_base_name = {"p": "p", "p0": "p", "T": "T", "T0": "T"}

        self.base_name_to_y_label = {
            "arhoA": "Mass equation variable",
            "arhouA": "Momentum equation variable",
            "arhoEA": "Energy equation variable",
            "rhoA": "Mass equation variable",
            "rhouA": "Momentum equation variable",
            "rhoEA": "Energy equation variable",
            "A": "Area",
            "vf": "Volume fraction",
            "p": "Pressure",
            "T": "Temperature",
            "u": "Velocity",
            "rho": "Density",
            "v": "Specific volume",
            "e": "Specific internal energy",
            "E": "Specific total energy",
            "h": "Specific enthalpy",
            "H": "Specific total enthalpy",
            "s": "Specific entropy",
            "p0": "Stagnation pressure",
            "T0": "Stagnation temperature"
        }

        self.base_name_to_scaling_factor = {"p": 1e-3, "p0": 1e-3}

        self.base_name_to_set_label = {
            "rhoA": "\\rho A",
            "rhouA": "\\rho u A",
            "rhoEA": "\\rho E A",
            "arhoA": "\\alpha\\rho A",
            "arhouA": "\\alpha\\rho u A",
            "arhoEA": "\\alpha\\rho E A",
            "A": "A",
            "vf": "\\alpha",
            "p": "p",
            "T": "T",
            "u": "u",
            "rho": "\\rho",
            "v": "v",
            "e": "e",
            "E": "E",
            "h": "h",
            "H": "H",
            "s": "s",
            "p0": "p_0",
            "T0": "T_0"
        }

        self.base_name_to_units = {
            "arhoA": "kg/m",
            "arhouA": "kg/s",
            "arhoEA": "J/m",
            "rhoA": "kg/m",
            "rhouA": "kg/s",
            "rhoEA": "J/m",
            "A": "m$^2$",
            "vf": "-",
            "p": "kPa",
            "T": "K",
            "u": "m/s",
            "rho": "kg/m$^3$",
            "v": "m$^3$/kg",
            "e": "J/kg",
            "E": "J/kg",
            "h": "J/kg",
            "H": "J/kg",
            "s": "J/kg",
            "p0": "kPa",
            "T0": "K"
        }

        self.base_names_with_products = ["rhoA", "rhouA", "rhoEA", "arhoA", "arhouA", "arhoEA"]

    def run(self, data):
        x_label = "Position, $x$"
        x = self.dof_handler.x
        x_min = min(x)
        x_max = max(x)
        if self.plot_by_mesh:
            x = self.dof_handler.separateNodalQuantityByMesh(x)

        # loop over sub-plots
        for i_subplot, plot_set in enumerate(self.plot_sets):
            y_label, set_labels, scaling_factor = self.determineLabelsAndScaling(plot_set)

            # only print lower-most x-label in each sub-plot column
            if self.print_x_label[i_subplot]:
                x_label_this_subplot = x_label
            else:
                x_label_this_subplot = ""

            # advance sub-plot
            if i_subplot == 0:
                plotter = Plotter(x_label_this_subplot, y_label, self.subplot_layout, \
                  default_size_x=self.size_x, default_size_y=self.size_y)
            else:
                plotter.nextSubplot(x_label_this_subplot, y_label)

            plotter.setXRange(x_min, x_max)

            if self.supplied_y_bounds:
                plotter.setYRange(self.y_bounds[0], self.y_bounds[1])

            # loop over quantities in sub-plot
            for i_name, name in enumerate(plot_set):
                if self.plot_by_mesh:
                    # plot each mesh separately
                    data_by_mesh = self.dof_handler.separateNodalQuantityByMesh(data[name])
                    for i_mesh in range(self.dof_handler.n_meshes):
                        modified_set_label = set_labels[i_name] + ", Mesh " + str(i_mesh + 1)
                        plotter.addSet(
                            x[i_mesh],
                            data_by_mesh[i_mesh],
                            modified_set_label,
                            color=i_mesh + 1,
                            marker=i_mesh + 1,
                            linetype=i_name,
                            scale=scaling_factor)
                else:
                    # plot all meshes together
                    plotter.addSet(
                        x,
                        data[name],
                        set_labels[i_name],
                        color=i_name + 1,
                        marker=i_name + 1,
                        linetype=0,
                        scale=scaling_factor)

                # don't plot "noise" in a near-constant solution
                plotter.fixNearConstantPlot()

                plotter.setLegendLocation(self.legend_location)

        # save plot
        plotter.save(self.file_name)

    def determineLabelsAndScaling(self, plot_sets):
        base_names = [plot_set.rstrip("1").rstrip("2") for plot_set in plot_sets]
        base_name_set = set(base_names)
        if len(base_name_set) > 1:
            if base_name_set in self.base_name_sets:
                y_label_base_name = self.base_name_to_y_label_base_name[base_names[0]]
            else:
                error("The set " + base_name_set + " is not a valid set for a sub-plot")
        else:
            y_label_base_name = base_names[0]

        # y-label
        y_label_set_label = self.base_name_to_set_label[y_label_base_name]
        units = self.base_name_to_units[y_label_base_name]
        y_label = self.base_name_to_y_label[y_label_base_name] + ", $" + y_label_set_label + "$ [" + units + "]"

        # set labels
        set_labels = list()
        for i, plot_set in enumerate(plot_sets):
            base_name = base_names[i]
            tentative_set_label = self.base_name_to_set_label[base_name]
            # add phase subscripts, if any
            if plot_set.endswith("1") or plot_set.endswith("2"):
                phase_number = plot_set[-1]
                if base_name in self.base_names_with_products:
                    set_label = "$(" + tentative_set_label + ")_" + phase_number + "$"
                else:
                    # prevent double subscript for stagnation quantities
                    if base_name.endswith("0"):
                        set_label = "$" + tentative_set_label.rstrip(
                            "0") + "{0," + phase_number + "}$"
                    else:
                        set_label = "$" + tentative_set_label + "_" + phase_number + "$"
            else:
                set_label = "$" + tentative_set_label + "$"
            set_labels.append(set_label)

        # scaling factor
        if y_label_base_name in self.base_name_to_scaling_factor:
            scaling_factor = self.base_name_to_scaling_factor[y_label_base_name]
        else:
            scaling_factor = 1

        return (y_label, set_labels, scaling_factor)
