import ROOT as rt
import uproot as uprt
import numpy as np


def root_to_numpy(hist):
    values = []
    bins = []
    error = []

    nbinsx = hist.GetNbinsX()

    for ibin in range(nbinsx):
        values.append(hist.GetBinContent(ibin+1))
        error.append(hist.GetBinError(ibin+1))
        bins.append(hist.GetXaxis().GetBinLowEdge(ibin+1))

    bins.append(bins[-1] + hist.GetBinWidth(nbinsx))
    values = np.array(values)
    bins = np.array(bins)
    error = np.array(error)

    return (values, bins), error


def get_raw_hist(file_path, hist_name):
    rt.TH1.AddDirectory(False)
    file = rt.TFile.Open(file_path)
    hist = file.Get(hist_name)
    file.Close()

    return hist


def get_raw_labels(file_path, hist_name):
    with uprt.open(file_path) as file:
        return file[hist_name].axis().labels()
