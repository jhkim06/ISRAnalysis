import ROOTFiles as rf
import matplotlib.pyplot as plt
import mplhep as hep


def adjust_x_lim(axis, bins):
    if bins[-1] - bins[0] >= 1e3:
        axis.set_xscale("log")
        axis.set_xlim(0.2, bins[-1])
    else:
        axis.set_xlim(bins[0], bins[-1])


class Plotter:
    fig = None
    axs = None

    def __init__(self, channel, period):
        self.channel = channel
        self.period = period
        self.path = "/Users/jhkim/cms_snu/data/Ultralegacy/"+channel+"/"+period
        self.root_file_handle = rf.ROOTFiles(self.path)

    def create_subplots(self, rows=1, cols=1, **kwargs):
        plt.ioff()
        hep.style.use("CMS")
        self.fig, self.axs = plt.subplots(rows, cols, **kwargs)
        plt.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.9, hspace=0.05)
        return self.axs

    def draw_hist(self, axis, file_key, hist_name, axis_steering="", **kwargs):
        hist_path = self.channel+self.period+"/"+hist_name
        (values, bins), errors = self.root_file_handle.get_hist(file_key, hist_path,
                                                                axis_steering=axis_steering)

        hep.histplot((values, bins), ax=axis, yerr=errors, **kwargs)
        adjust_x_lim(axis, bins)

    def draw_stack(self, axis, file_key_list, hist_name, axis_steering="", use_mplhep=True, **kwargs):
        hist_path = self.channel + self.period + "/" + hist_name
        (value_list, bins), error_list = self.root_file_handle.make_stack_list(file_key_list, hist_path, axis_steering)

        if use_mplhep:
            hep.histplot(value_list, ax=axis, stack=True, bins=bins, **kwargs)
        else:
            width = bins[1:]-bins[0:-1]
            bottom = 0
            for index, value in enumerate(value_list):
                axis.hist(bins[:-1], bins, weights=value/width, histtype='bar', bottom=bottom)
                bottom += value/width

    def clear_all_axis(self):
        if len(self.axs) > 1:
            for axis in self.axs:
                axis.cla()
        else:
            self.axs.cla()

    def save_plot(self, plot_name):
        self.fig.savefig(plot_name+".pdf")
