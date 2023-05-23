import ROOTFiles as rf
import matplotlib.pyplot as plt
import mplhep as hep


class Plotter:
    fig = None
    axs = None

    def __init__(self, channel, period):
        self.channel = channel
        self.period = period
        self.path = "/Users/jhkim/cms_snu/data/Ultralegacy/"+channel+"/"+period
        self.root_file_handle = rf.ROOTFiles(self.path)

    def create_subplots(self, rows=1, cols=1, **kwargs):
        hep.style.use("CMS")
        self.fig, self.axs = plt.subplots(rows, cols, **kwargs)
        plt.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.9, hspace=0.05)
        return self.axs

    def draw_hist(self, axis, file_key, hist_name, axis_steering="", **kwargs):
        hist_path = self.channel+self.period+"/"+hist_name
        (values, bins), errors = self.root_file_handle.get_hist(file_key, hist_path,
                                                                axis_steering=axis_steering)

        hep.histplot((values, bins), ax=axis, yerr=errors, **kwargs)

    def save_plot(self, plot_name):
        self.fig.savefig(plot_name+".pdf")
