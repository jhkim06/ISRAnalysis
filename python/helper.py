import ROOT as rt
import uproot as uprt
import numpy as np
import re
from collections import namedtuple

Hist = namedtuple('Hist', ['values', 'bins', 'errors'])


def root_to_numpy(hist):
    values = []
    bins = []
    error = []

    nbinsx = hist.GetNbinsX()

    for ibin in range(nbinsx):
        values.append(hist.GetBinContent(ibin + 1))
        error.append(hist.GetBinError(ibin + 1))
        bins.append(hist.GetXaxis().GetBinLowEdge(ibin + 1))

    bins.append(bins[-1] + hist.GetBinWidth(nbinsx))
    values = np.array(values)
    bins = np.array(bins)
    errors = np.array(error)

    h = Hist(values, bins, errors)
    return h


def return_hist(raw_hist, root_hist=False):
    if root_hist:
        return raw_hist
    else:
        return root_to_numpy(raw_hist)


def get_mass_window_edges(file_path, hist_path):
    file = rt.TFile.Open(file_path)
    # make function to get TUnfoldBinning
    directory = "/".join(hist_path.split("/")[:-1])
    hist_name = hist_path.split("/")[-1]   # assuming hists in the last subdirectory
    # get bin definition from file using hist name
    matches = re.findall(r'(\[[^]]*])', hist_name)
    bin_name = "[tunfold-bin]_" + re.sub(r'(\([^]]*\))', '', matches[1]) + "_" + "_".join(matches[2:])
    bin_path = directory + "/" + bin_name
    unfold_bin = file.Get(bin_path)

    edges = unfold_bin.GetDistributionBinning(1)
    edge_list = []
    for edge in edges:
        edge_list.append(edge)

    return edge_list


def get_raw_hist(file_path, hist_path, norm=False, stat_variation=0):

    rt.TH1.AddDirectory(False)
    file = rt.TFile.Open(file_path)
    hist = file.Get(hist_path)

    if stat_variation != 0:  # nominal +/- stat
        nbins = hist.GetNbinsX()
        for ibin in range(nbins+1):
            hist.SetBinContent(ibin, hist.GetBinContent(ibin) + stat_variation * hist.GetBinError(ibin))
    if norm:
        integral = hist.Integral()
        hist.Scale(1./integral)

    file.Close()
    return hist


def get_summary_statistics(hist, stat_type="mean", prob=0.6):
    Stat = namedtuple('stat', ['value', 'error'])

    if stat_type == "mean":
        value = hist.GetMean()
        error = hist.GetMeanError()
    elif stat_type == "quantile":
        hist_temp = hist.Clone("quantile")
        hist_temp.Scale(1, "width")

        length = 1
        q = np.zeros(length)
        prob_array = np.array([prob])
        hist_temp.GetQuantiles(length, q, prob_array)
        value = q[0]
        error = 0  # FIXME
    else:
        return None
    return Stat(value, error)


def get_raw_labels(file_path, hist_name):
    with uprt.open(file_path) as file:
        return file[hist_name].axis().edges(), file[hist_name].axis().labels()


def calculate_squared_root_sum(raw_df, reg_expression, new_col_name="total error"):

    selected_cols = raw_df.filter(regex=reg_expression)
    squared_sum = (selected_cols ** 2).sum(axis=1)
    root_of_squared_sum = np.sqrt(squared_sum)

    raw_df[new_col_name] = root_of_squared_sum
