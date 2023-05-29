import matplotlib.pyplot as plt
import mplhep as hep
from matplotlib.ticker import (FixedLocator, FixedFormatter)
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import numpy as np


def adjust_x_lim(axis, bins):
    if bins[-1] >= 1e3:
        axis.set_xscale("log")
        if bins[0] == 0:
            axis.set_xlim(0.2, bins[-1])
        else:
            axis.set_xlim(bins[0], bins[-1])
    else:
        axis.set_xlim(bins[0], bins[-1])


class Plotter:
    def __init__(self, rows, cols, **kwargs):
        plt.ioff()
        hep.style.use("CMS")
        plt.rcParams['axes.linewidth'] = 1
        plt.rcParams['hatch.linewidth'] = 0.1
        self.rows = rows
        self.cols = cols
        self.fig, self.axs = plt.subplots(self.rows, self.cols, **kwargs)
        plt.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.9, hspace=0.05)

    def get_axis(self, row, col):
        if self.rows == 1 and self.cols == 1:
            return self.axs
        elif self.rows == 1 or self.cols == 1:
            if self.rows == 1:
                index = col
            else:
                index = row
            return self.axs[index]
        else:
            return self.axs[row][col]

    def draw_hist(self, hist, row=0, col=0, **kwargs):
        # axis int or tuple
        axis = self.get_axis(row, col)
        artist = hep.histplot((hist.values, hist.bins), ax=axis, **kwargs)
        adjust_x_lim(axis, hist.bins)
        return artist

    def draw_stack(self, stack, row=0, col=0, use_mplhep=True, **kwargs):
        axis = self.get_axis(row, col)
        if use_mplhep:
            hep.histplot(stack.values_list, ax=axis, stack=True, bins=stack.bins, **kwargs)
        else:
            if "label" in kwargs:
                label = kwargs["label"]
            else:
                label = ["" for _ in range(len(stack.value_list))]
            bottom = 0
            width = stack.bins[1:] - stack.bins[0:-1]
            for index, values in enumerate(stack.values_list):
                axis.hist(stack.bins[:-1], stack.bins, weights=values/width, histtype='bar',
                          bottom=bottom, label=label[index])
                bottom += values/width
        adjust_x_lim(axis, stack.bins)

    def draw_error_boxes(self, hist, row=0, col=0, set_y_ones=False, **kwargs):
        axis = self.get_axis(row, col)
        center = hist.bins[:-1]+np.diff(hist.bins)/2.
        # down, up error for x axis (bin width)
        xerrs = np.expand_dims(np.diff(hist.bins)/2., axis=0)
        xerrs = np.append(xerrs, xerrs, axis=0)
        # down, up error for y axis
        yerrs = np.expand_dims(hist.errors, axis=0)
        yerrs = np.append(yerrs, yerrs, axis=0)

        if set_y_ones:
            data = np.ones(len(hist.values))
        else:
            data = hist.values
        # Loop over data points; create box from errors at each point
        errorboxes = [Rectangle((x - xe[0], y - ye[0]), xe.sum(), ye.sum(),
                                fill=False, edgecolor=None, facecolor=None)
                      for x, y, xe, ye in zip(center, data, xerrs.T, yerrs.T)]
        # Create patch collection with specified colour/alpha
        pc = PatchCollection(errorboxes, match_original=True, hatch="xx",
                             linewidth=0, zorder=100, **kwargs)
        axis.add_collection(pc)
        # Plot errorbars
        artists = axis.errorbar(center, data, fmt='none')
        # artist[0][0][0].get_data()
        # artist[0][0][2][0].get_segments()  # bin content
        # artist[0][0][2][1].get_segments()  # error boundary
        # Add collection to axes
        axis.add_collection(pc)
        return artists

    def set_labels(self, bins, labels, row=0, col=0):
        axis = self.get_axis(row, col)
        axis.minorticks_off()
        major_ticks = FixedLocator(bins[1:])
        axis.xaxis.set_major_locator(major_ticks)
        axis.xaxis.set_major_formatter(FixedFormatter(labels))

    def draw_isr_data_frame(self, data_frame, row=0, col=0, **kwargs):
        axis = self.get_axis(row, col)
        x_column_name = "mass"
        y_column_name = "pt"
        error_name = "stat_error"

        axis.errorbar(data_frame[x_column_name], data_frame[y_column_name],
                      xerr=data_frame[x_column_name+"_"+error_name],
                      yerr=data_frame[y_column_name+"_"+error_name],
                      **kwargs)

    def clear_all_axis(self):
        plt.rcdefaults()
        hep.style.use("CMS")
        plt.rcParams['axes.linewidth'] = 1
        plt.rcParams['hatch.linewidth'] = 0.1

        if self.rows == 1 and self.cols == 1:
            self.axs.cla()
        elif self.rows == 1 or self.cols == 1:
            if self.rows == 1:
                length = self.cols
            else:
                length = self.rows
            for index in range(length):
                self.axs[index].cla()
        else:
            for row in range(self.rows):
                for col in range(self.cols):
                    self.axs[row][col].cla()

    def save_plot(self, plot_name):
        self.fig.savefig(plot_name+".pdf")
