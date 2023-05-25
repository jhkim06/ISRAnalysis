import numpy as np

import ROOTFiles as rf
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
    fig = None
    axs = None
    drawn_file_list = []
    denominator_index = 0

    def __init__(self, channel, period):
        self.channel = channel
        self.period = period
        self.path = "/Users/jhkim/cms_snu/data/Ultralegacy/"+channel+"/"+period
        self.root_file_handle = rf.ROOTFiles(self.channel, self.period)

    def create_subplots(self, rows=1, cols=1, **kwargs):
        self.clear_all_axis()
        plt.ioff()
        hep.style.use("CMS")
        self.fig, self.axs = plt.subplots(rows, cols, **kwargs)
        plt.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.9, hspace=0.05)
        return self.axs

    def draw_hist(self, axis, file_key, hist_name, axis_steering="", **kwargs):
        hist_path = self.channel+self.period+"/"+hist_name
        (values, bins), errors = self.root_file_handle.get_hist(file_key, hist_path,
                                                                axis_steering=axis_steering)
        if "yerr" in kwargs.keys():
            if not kwargs["yerr"]:
                errors = 0
                del kwargs['yerr']
        hep.histplot((values, bins), ax=axis, yerr=errors, xerr=True, **kwargs)
        adjust_x_lim(axis, bins)

        labels = self.root_file_handle.get_raw_labels(file_key, hist_path)
        if labels is not None:
            set_labels(axis, bins, labels)

        self.drawn_file_list.append(file_key)

    def draw_stack(self, axis, file_key_list, hist_name, axis_steering="", use_mplhep=True, **kwargs):
        hist_path = self.channel + self.period + "/" + hist_name
        (value_list, bins), error_list = self.root_file_handle.make_stack_list(file_key_list, hist_path, axis_steering)

        if use_mplhep:
            hep.histplot(value_list, ax=axis, stack=True, bins=bins, **kwargs)
        else:
            if "label" in kwargs:
                label = kwargs["label"]
            else:
                label = ["" for i in range(len(value_list))]
            width = bins[1:]-bins[0:-1]
            bottom = 0
            for index, value in enumerate(value_list):

                axis.hist(bins[:-1], bins, weights=value/width, histtype='bar', bottom=bottom, label=label[index])
                bottom += value/width

        adjust_x_lim(axis, bins)
        self.drawn_file_list.append(file_key_list)

    def set_denominator(self, drawn_file_index):
        self.denominator_index = drawn_file_index

    def draw_ratio(self, axis, nominator_file_list, nominator_hist_name, denominator_hist_name, axis_steering="",
                   **kwargs):
        nominator_hist_path = self.channel + self.period + "/" + nominator_hist_name
        denominator_hist_path = self.channel + self .period + "/" + denominator_hist_name

        if type(self.drawn_file_list[self.denominator_index]) == list:
            denominator_file_list = self.drawn_file_list[self.denominator_index]
        else:
            denominator_file_list = [self.drawn_file_list[self.denominator_index]]

        (values, bins), errors = self.root_file_handle.get_ratio_hist(nominator_file_list, denominator_file_list,
                                                                      nominator_hist_path, denominator_hist_path,
                                                                      axis_steering=axis_steering)
        # Draw
        hep.histplot((values, bins), ax=axis, yerr=errors, xerr=True, **kwargs)
        adjust_x_lim(axis, bins)

        labels = self.root_file_handle.get_raw_labels(nominator_file_list[0], nominator_hist_path)
        if labels is not None:
            set_labels(axis, bins, labels)

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

    def draw_isr_quantiles(self):
        pass

    def clear_all_axis(self):
        plt.rcdefaults()
        hep.style.use("CMS")
        self.drawn_file_list.clear()
        if self.axs is not None:
            if type(self.axs) == np.ndarray:
                for axis in self.axs:
                    axis.cla()
                    axis.tick_params(axis='both', reset=True)
            else:
                self.axs.cla()
                self.axs.tick_params(axis='both', reset=True)

    def save_plot(self, plot_name):
        self.fig.savefig(plot_name+".pdf")
