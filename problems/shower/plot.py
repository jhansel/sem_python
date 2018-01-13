from sem_python.utilities.file_utilities import readCSVData
from sem_python.output.Plotter import Plotter

n_pipes = 3
problem_name = "shower"
x_max = 2
y_labels = ["Pressure, $p$ [kPa]", "Temperature, $T$ [K]"]
y_names = ["p", "T"]
legend_locations = ["upper right", "upper right"]

data = list()
for i_pipe in range(n_pipes):
    data.append(readCSVData(problem_name + "_" + str(i_pipe + 1) + ".csv"))


def makePlot(y_label, y_name, leg_loc):
    plotter = Plotter("$x$", y_label)
    for i_pipe in range(n_pipes):
        plotter.addSet(
            data[i_pipe]["x"],
            data[i_pipe][y_name],
            "Mesh " + str(i_pipe + 1),
            color=i_pipe + 1,
            linetype=0,
            marker=i_pipe + 1)
    plotter.setXRange(0, x_max)
    plotter.fixNearConstantPlot()
    plotter.setLegendLocation(leg_loc)
    plotter.save(problem_name + "_" + y_name + ".pdf")


for y_label, y_name, leg_loc in zip(y_labels, y_names, legend_locations):
    makePlot(y_label, y_name, leg_loc)
