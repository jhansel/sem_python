from sem_python.utilities.file_utilities import readCSVData
from sem_python.output.Plotter import Plotter

no_junction = readCSVData("no_junction.csv")
junction_1 = readCSVData("solution_1.csv")
junction_2 = readCSVData("solution_2.csv")

x_no_junction = no_junction["x"]
x_1 = junction_1["x"]
x_2 = junction_2["x"]


def makePlot(y_label, y_name, leg_loc):
    plotter = Plotter("$x$", y_label)
    plotter.addSet(x_no_junction, no_junction[y_name], "No Junction", color=0, linetype=0, marker=0)
    plotter.addSet(x_1, junction_1[y_name], "With Junction, Mesh 1", color=1, linetype=0, marker=1)
    plotter.addSet(x_2, junction_2[y_name], "With Junction, Mesh 2", color=2, linetype=0, marker=2)
    plotter.fixNearConstantPlot()
    plotter.setLegendLocation(leg_loc)
    plotter.save("equal_area_" + y_name + ".pdf")


y_labels = [
    "Density, $\\rho$ [kg/m$^3$]", "Velocity, $u$ [m/s]", "Pressure, $p$ [kPa]", "Area, $A$ [m$^2$]"
]
y_names = ["rho", "u", "p", "A"]
legend_locations = ["lower right", "upper left", "upper right", "upper right"]
for y_label, y_name, leg_loc in zip(y_labels, y_names, legend_locations):
    makePlot(y_label, y_name, leg_loc)
