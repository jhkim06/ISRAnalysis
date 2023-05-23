import ROOT as rt
import uproot as uprt
import numpy as np
import re


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
    error = np.array(error)

    return (values, bins), error


def get_mass_window_edges(file_path, hist_path):
    file = rt.TFile.Open(file_path)

    # make function to get TUnfoldBinning
    directory = hist_path.split("/")[0]
    hist_name = hist_path.split("/")[1]  # assuming no sub-directory
    # get bin definition from file using hist name
    matches = re.findall(r'(\[[^\]]*\])', hist_name)
    bin_name = "[tunfold-bin]_" + "_".join(matches[1:])
    bin_path = directory + "/" + bin_name
    unfold_bin = file.Get(bin_path)

    edges = unfold_bin.GetDistributionBinning(1)
    edge_list = []
    for edge in edges:
        edge_list.append(edge)

    return edge_list


def get_raw_hist(file_path, hist_path, axis_steering=""):
    rt.TH1.AddDirectory(False)
    file = rt.TFile.Open(file_path)
    hist = file.Get(hist_path)

    # option for TUnfoldBinning
    if axis_steering != "":
        # make function to get TUnfoldBinning
        directory = hist_path.split("/")[0]
        hist_name = hist_path.split("/")[1]  # assuming no sub-directory
        # get bin definition from file using hist name
        matches = re.findall(r'(\[[^\]]*\])', hist_name)
        bin_name = "[tunfold-bin]_" + "_".join(matches[1:])
        bin_path = directory + "/" + bin_name
        unfold_bin = file.Get(bin_path)

        original_axis_binning = True
        hist = unfold_bin.ExtractHistogram(hist_name + "_extracted_" + axis_steering,
                                           hist, 0, original_axis_binning, axis_steering)

    file.Close()
    return hist


def attach_hprefix_to_hist_name(file_name, hist_path):
    file_name_parsing = file_name.split(":")
    if len(file_name_parsing) == 2:  # hprefix set in config file
        directory = hist_path.split("/")[0]
        hist_name = hist_path.split("/")[1]
        hist_path = directory + "/" + file_name_parsing[1] + hist_name
        return hist_path
    else:
        return hist_path

def get_raw_labels(file_path, hist_name):
    with uprt.open(file_path) as file:
        return file[hist_name].axis().labels()
