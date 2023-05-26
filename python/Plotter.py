import matplotlib.pyplot as plt
import mplhep as hep
from matplotlib.ticker import (FixedLocator, FixedFormatter)


def adjust_x_lim(axis, bins):
    if bins[-1] >= 1e3:
        axis.set_xscale("log")
        if bins[0] == 0:
            axis.set_xlim(0.2, bins[-1])
        else:
            axis.set_xlim(bins[0], bins[-1])
    else:
        axis.set_xlim(bins[0], bins[-1])


def set_labels(axis, bins, labels):
    axis.minorticks_off()
    major_ticks = FixedLocator(bins[1:])
    axis.xaxis.set_major_locator(major_ticks)
    axis.xaxis.set_major_formatter(FixedFormatter(labels))


class Plotter:
    def __init__(self, rows, cols, **kwargs):
        plt.ioff()
        hep.style.use("CMS")
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
        hep.histplot((hist.values, hist.bins), ax=axis, **kwargs)
        adjust_x_lim(axis, hist.bins)

    def draw_stack(self, stack, row=0, col=0, use_mplhep=True, **kwargs):
        # stack: list of named tuple [(),()]
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

    '''
    def draw_isr_means(self, axis, file_key, pt_hist_name, mass_hist_name,
                       pt_axis_steering="dipt[];dimass[UOC]", mass_axis_steering="dimass[UO];dipt[OC0]", **kwargs):

        pt_hist_path = self.channel + self.period + "/" + pt_hist_name
        mass_hist_path = self.channel + self.period + "/" + mass_hist_name
        isr_data = self.root_file_handle.get_isr_mean_dataframe(file_key, pt_hist_path, mass_hist_path,
                                                                pt_axis_steering, mass_axis_steering)

        axis.errorbar(
            isr_data["mean_mass"],
            isr_data["mean_pt"],
            xerr=isr_data["mean_mass_stat_error"],
            yerr=isr_data["mean_pt_stat_error"], **kwargs
        )
    '''

    def clear_all_axis(self):
        plt.rcdefaults()
        hep.style.use("CMS")

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
