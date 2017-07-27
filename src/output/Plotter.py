import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

# colors:
# 0: black
# 1: red
# 2: orange
# 3: green
# 4: blue
# 5: purple
# 6: hot pink
colors = ['k','indianred','orange','lightgreen','cornflowerblue','slateblue','orchid']
linetypes = ['-','--',':']

## Class for creating plots
class Plotter(object):
    def __init__(self, x_label, y_label, n_subplots=1, logscale_y=False, size_x=12, size_y=9):
        # set figure size
        plt.figure(figsize=(size_x, size_y))

        # use latex
        plt.rc('text', usetex=True)
        plt.rc('font', family='sans-serif')

        # get number of subplots in each dimension
        if isinstance(n_subplots, tuple):
          self.n_subplots_x, self.n_subplots_y = n_subplots
        else:
          self.n_subplots_y = n_subplots
          self.n_subplots_x = 1
        self.n_subplots = self.n_subplots_x * self.n_subplots_y
        self.current_subplot_index = 1

        # get axes of first subplot
        self.ax = plt.subplot(self.n_subplots_y, self.n_subplots_x, 1)

        # turn off offset
        self.ax.get_yaxis().get_major_formatter().set_useOffset(False)

        # log scale for first axis
        self.logscale_y = logscale_y

        # set the axis labels
        plt.xlabel(x_label)
        plt.ylabel(y_label)

        # flag for having set custom x range
        self.set_custom_x_range = False

        # legend
        self.put_legend_outside = False
        self.frame_legend = False
        self.legend_entries = list()

    def adjustLeftMargin(self, margin):
        plt.subplots_adjust(left=margin)

    def useLegendFrame(self):
        self.frame_legend = True

    def setYFormat(self, y_format):
        self.ax.yaxis.set_major_formatter(FormatStrFormatter(y_format))

    def setXRange(self,xmin,xmax):
        self.set_custom_x_range = True
        self.xmin = xmin
        self.xmax = xmax
        self.ax.set_xlim([self.xmin,self.xmax])

    def setYRange(self, ymin, ymax):
        self.ax.set_ylim([ymin, ymax])

    def fixNearConstantPlot(self, dy_min=1.0):
      ymin, ymax = self.ax.get_ylim()
      rel_diff = (ymax - ymin) / max(1e-15, abs(ymax))
      if (rel_diff < 1e-15):
        yavg = 0.5 * (ymin + ymax)
        ymin_new = yavg - 0.5 * dy_min
        ymax_new = yavg + 0.5 * dy_min
        self.setYRange(ymin_new, ymax_new)

    def putLegendOutside(self):
      self.put_legend_outside = True

      # shrink axis by 20% to allow room for legend outside of figure
      box = self.ax.get_position()
      self.ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    def nextSubplot(self, x_label, y_label, logscale_y=False):
        # update logscale
        self.logscale_y = logscale_y

        # create legend for previous subplot and then clear legend entries list
        if (self.put_legend_outside):
          self.ax.legend(self.legend_entries, loc='center left', frameon=self.frame_legend,
              bbox_to_anchor=(1,0.5), prop={'size':12})
        else:
          self.ax.legend(self.legend_entries, frameon=self.frame_legend, prop={'size':12})
        self.legend_entries = list()

        # create new subplot
        self.current_subplot_index += 1
        self.ax = plt.subplot(self.n_subplots_y, self.n_subplots_x, self.current_subplot_index)

        # turn off offset
        self.ax.get_yaxis().get_major_formatter().set_useOffset(False)

        if self.put_legend_outside:
          box = self.ax.get_position()
          self.ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        # set the x and y labels
        plt.xlabel(x_label)
        plt.ylabel(y_label)

        # set the x range if a custom x range was provided
        if (self.set_custom_x_range):
            self.ax.set_xlim([self.xmin,self.xmax])

    def addGrid(self):
        plt.grid(True)

    def addSet(self,x,y,set_name,color=-1,linetype=0,scale=1):
        # scale y values
        y = [yi * scale for yi in y]

        if (color >= 0):
            if (self.logscale_y):
              plt.semilogy(x,y,linetypes[linetype],color=colors[color])
            else:
              plt.plot(x,y,linetypes[linetype],color=colors[color])
        else:
            plt.plot(x,y,linetypes[linetype])
        self.legend_entries.append(set_name)

    def save(self, outputfile):
        # create legend for final subplot
        if (self.put_legend_outside):
          self.ax.legend(self.legend_entries, loc='center left', frameon=self.frame_legend,
              bbox_to_anchor=(1,0.5), prop={'size':12})
        else:
          self.ax.legend(self.legend_entries, frameon=self.frame_legend, prop={'size':12})

        # save the figure
        plt.savefig(outputfile, dpi=300)

if (__name__ == '__main__'):
    plotter = Plotter('x', 'y', 2)
    x = [1, 2, 3]
    y1 = [1, 2, 3]
    y2 = [2, 4, 6]
    z1 = [2, 3, 5]
    z2 = [0, 1, 4]
    plotter.addSet(x,y1,'y1',color=1)
    plotter.addSet(x,y2,'y2',color=2)
    plotter.nextSubplot('z')
    plotter.addSet(x,z1,'z1',color=3)
    plotter.addSet(x,z2,'z2')
    plotter.save('Plotter_test.pdf')
