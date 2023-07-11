import json

import numpy as np
from numpy import ndarray
from helper import *
import pandas as pd
from collections import namedtuple


# ROOT histogram to numpy, pandas etc
class ROOTFiles:

    def __init__(self, path, channel, period, *args, sample_config="sample_config.json"):
        self.file_dir = {"data": path + "data/" + channel + "/" + period,
                         "mc": path + "mc/" + period}
        self.hist_dir = channel + period + "/"
        self.input_files = dict()
        self.set_input_files(path, channel+"_"+sample_config, *args)
        self.data_period = period

    def make_stack_list(self, file_key_list, hist_name, axis_steering="", file_sel=""):
        values_list: list[ndarray] = []
        bins = None
        errors_list = []
        for file_key in file_key_list:
            hist = self.get_hist(file_key, hist_name, axis_steering, file_sel=file_sel)
            values_list.append(hist.values)
            bins = hist.bins
            errors_list.append(hist.errors)

        Stack = namedtuple('Stack', ['values_list', 'bins', 'errors_list'])
        stack = Stack(values_list, bins, errors_list)
        return stack

    def get_file_directory(self, sample_type, file_sel):
        file_postfix = "/"
        if file_sel != "":
            file_postfix = "/" + file_sel.split(":")[0] + "/"
        return self.file_dir[sample_type] + file_postfix

    def get_hist_path(self, file_name, hist_name, file_sel=""):
        file_name_parsing = file_name.split(":")
        full_hist_dir = self.hist_dir
        if file_sel != "":
            full_hist_dir = self.hist_dir + file_sel.split(":")[1] + "/"
        if len(file_name_parsing) == 2:  # hprefix set in config file
            hist_path = full_hist_dir + file_name_parsing[1] + hist_name  # ex) tau_
            return hist_path
        else:
            hist_path = full_hist_dir + hist_name
            return hist_path

    # get values, bin, error as numpy array
    def get_hist(self, file_key, hist_name, axis_steering="", root_hist=False,
                 norm=False, binwnorm=False, file_sel="", stat_variation=0):
        sample_type, hist_label, file_name = file_key.split("/")
        raw_hist = None

        if sample_type == "data:bkg_subtracted":
            raw_hist = self.get_bkg_subtracted_data_hist(hist_name, axis_steering, root_hist=True,
                                                         file_sel=file_sel, bkg_stat=stat_variation)
        else:
            if file_name == "all":  # sum all histograms
                for index, file in enumerate(self.input_files[file_sel][sample_type][hist_label]):
                    file_path = self.get_file_directory(sample_type, file_sel) + \
                                self.input_files[file_sel][sample_type][hist_label][file]
                    hist_path = self.get_hist_path(file, hist_name, file_sel)
                    if index == 0:
                        raw_hist = get_raw_hist(file_path, hist_path, axis_steering, stat_variation=stat_variation)
                    else:
                        raw_hist.Add(get_raw_hist(file_path, hist_path, axis_steering, stat_variation=stat_variation))
            else:
                hist_path = self.get_hist_path(file_name, hist_name, file_sel)
                file_path = self.get_file_directory(sample_type, file_sel) + \
                            self.input_files[file_sel][sample_type][hist_label][file_name]
                raw_hist = get_raw_hist(file_path, hist_path, axis_steering, stat_variation=stat_variation)

        if norm:
            integral = raw_hist.Integral()
            raw_hist.Scale(1. / integral)
        if binwnorm:
            raw_hist.Scale(1., "width")

        if root_hist:
            return raw_hist
        else:
            return root_to_numpy(raw_hist)

    def get_combined_hist(self, file_key_list, hist_name, axis_steering="", root_hist=False,
                          norm=False, binwnorm=False, file_sel="", stat_variation=0):
        raw_hist = None
        for index, file_key in enumerate(file_key_list):
            hist = self.get_hist(file_key, hist_name, axis_steering, root_hist=True, file_sel=file_sel,
                                 stat_variation=stat_variation)
            if index == 0:
                raw_hist = hist
            else:
                raw_hist.Add(hist)
        if norm:
            integral = raw_hist.Integral()
            raw_hist.Scale(1. / integral)
        if binwnorm:
            raw_hist.Scale(1., "width")

        if root_hist:
            return raw_hist
        else:
            return root_to_numpy(raw_hist)

    def get_bkg_subtracted_data_hist(self, hist_name, axis_steering="", root_hist=False,
                                     norm=False, binwnorm=False, file_sel="", bkg_stat=0):
        data_hist = self.get_hist("data/" + self.data_period + "/all", hist_name, axis_steering, root_hist=True,
                                  file_sel=file_sel)
        # make file_key_list to background mc
        bkg_key_list = ["mc/" + bkg_key + "/all" for bkg_key in self.input_files[""]["mc"].keys()
                        if "background" in bkg_key]
        bkg_hist = self.get_combined_hist(bkg_key_list, hist_name, axis_steering, root_hist=True,
                                          file_sel=file_sel, stat_variation=bkg_stat)
        bkg_subtracted_data_hist = data_hist.Clone("data_bkg_subtracted")
        bkg_subtracted_data_hist.Add(bkg_hist, -1)
        if norm:
            integral = bkg_subtracted_data_hist.Integral()
            bkg_subtracted_data_hist.Scale(1. / integral)
        if binwnorm:
            bkg_subtracted_data_hist.Scale(1., "width")

        if root_hist:
            return bkg_subtracted_data_hist
        else:
            return root_to_numpy(bkg_subtracted_data_hist)

    def get_hist_systematic_bands(self, systematic_dict, file_key, hist_name, axis_steering="",
                                  norm=False, binwnorm=False, file_sel=""):
        nominal_hist = self.get_hist(file_key, hist_name, axis_steering, root_hist=True,
                                     norm=norm, binwnorm=binwnorm)
        # symmetric systematic
        systematic_list = []
        for systematic in systematic_dict.keys():
            delta_list = []
            for variation in systematic_dict[systematic]:
                variation_postfix = ""
                if variation != "":
                    variation_postfix += "_" + variation
                full_hist_name = hist_name + "_" + systematic + variation_postfix
                systematic_hist = self.get_hist(file_key, full_hist_name, axis_steering, root_hist=True, norm=norm,
                                                binwnorm=binwnorm, file_sel=file_sel)
                systematic_hist.Add(nominal_hist, -1)
                h = root_to_numpy(systematic_hist)
                delta_list.append(np.absolute(h.values))
            systematic_list.append(np.square(np.mean(delta_list, axis=0)))
        sqrt_sum = np.sqrt(np.sum(np.array(systematic_list), axis=0))  # square root sum
        return sqrt_sum

    def get_combined_hist_systematic_bands(self, systematic_dict, file_key_list, hist_name, axis_steering="",
                                           norm=False, binwnorm=False, file_sel=""):
        nominal_hist = self.get_combined_hist(file_key_list, hist_name, axis_steering=axis_steering, root_hist=True,
                                              norm=norm, binwnorm=binwnorm)
        # symmetric systematic
        systematic_list = []
        for systematic in systematic_dict.keys():
            delta_list = []
            for variation in systematic_dict[systematic]:
                variation_postfix = ""
                if variation != "":
                    variation_postfix += "_" + variation
                full_hist_name = hist_name + "_" + systematic + variation_postfix
                systematic_hist = self.get_combined_hist(file_key_list, full_hist_name, axis_steering=axis_steering,
                                                         root_hist=True,
                                                         norm=norm, binwnorm=binwnorm, file_sel=file_sel)
                systematic_hist.Add(nominal_hist, -1)
                h = root_to_numpy(systematic_hist)
                delta_list.append(np.absolute(h.values))
            systematic_list.append(np.square(np.mean(delta_list, axis=0)))
        sqrt_sum = np.sqrt(np.sum(np.array(systematic_list), axis=0))  # square root sum
        return sqrt_sum

    def get_scaled_hist(self, file_key_list, hist_name, axis_steering="", root_hist=False,
                        norm=False, binwnorm=False, file_sel="", systematic_dict=None):
        hist = self.get_combined_hist(file_key_list, hist_name, axis_steering, root_hist,
                                      norm, binwnorm, file_sel=file_sel)
        sys_error = None
        if systematic_dict is not None:
            sys_error = self.get_combined_hist_systematic_bands(systematic_dict, file_key_list, hist_name, axis_steering,
                                                                norm, binwnorm, file_sel)
            sys_error = np.nan_to_num(sys_error / hist.values)

        ones = np.nan_to_num(hist.values / hist.values)
        stat_error = np.nan_to_num(hist.errors / hist.values)
        return Hist(ones, hist.bins, stat_error), sys_error

    # TODO how to handle systematic error for division
    # currently considered statistical error assuming uncorrelated histograms
    def get_ratio_hist(self, file1_key_list, file2_key_list,
                       hist1_name, hist2_name, axis_steering="", root_hist=False,
                       other_instance_for_file2=None, norm=False, binwnorm=False,
                       file_sel1="", file_sel2=""):
        combined_nominator = self.get_combined_hist(file1_key_list, hist1_name,
                                                    axis_steering, root_hist=True, norm=norm,
                                                    binwnorm=binwnorm, file_sel=file_sel1)
        if other_instance_for_file2 is not None:
            combined_denominator = other_instance_for_file2.get_combined_hist(file2_key_list, hist2_name,
                                                                              axis_steering, root_hist=True,
                                                                              norm=norm, binwnorm=binwnorm,
                                                                              file_sel=file_sel2)
        else:
            combined_denominator = self.get_combined_hist(file2_key_list, hist2_name,
                                                          axis_steering, root_hist=True,
                                                          norm=norm, binwnorm=binwnorm, file_sel=file_sel2)
        ratio = combined_nominator.Clone("ratio")
        ratio.Divide(combined_denominator)

        if root_hist:
            return ratio
        else:
            return root_to_numpy(ratio)

    def get_ratio_hist2(self, file_key_list, hist_name, axis_steering="", nominator=None,
                        denominator=None, root_hist=False, norm=False, binwnorm=False, file_sel=""):
        if nominator is None and denominator is None:
            return
        if nominator is not None and denominator is not None:
            return

        if nominator is None:
            nominator = self.get_combined_hist(file_key_list, hist_name,
                                               axis_steering, root_hist=True, norm=norm,
                                               binwnorm=binwnorm, file_sel=file_sel)
            ratio = nominator.Clone("ratio")

        else:
            ratio = nominator.Clone("ratio")
            denominator = self.get_combined_hist(file_key_list, hist_name,
                                                 axis_steering, root_hist=True, norm=norm,
                                                 binwnorm=binwnorm, file_sel=file_sel)
        ratio.Divide(denominator)

        if root_hist:
            return ratio
        else:
            return root_to_numpy(ratio)

    # get ISR mean DataFrame from 2D pt-mass histogram
    def get_isr_dataframe(self, file_key, pt_hist_name, mass_hist_name, file_sel="", stat_type="mean", prob=0.5,
                          pt_steering="dipt[];dimass[UOC]", mass_steering="dimass[UO];dipt[OC0]", sys=False):

        file_path = self.get_file_path(file_key, file_sel=file_sel)
        # pt_hist_path = self.hist_dir + pt_hist_name
        sample_type, hist_label, file_name = file_key.split("/")
        # to get the edges of mass window
        if "data:bkg_subtracted" == sample_type:
            sample_type = "data"
        if file_name == "all":
            file_name = [*self.input_files[file_sel][sample_type][hist_label].keys()][0]  # use one of samples
        pt_hist_path = self.get_hist_path(file_name, pt_hist_name, file_sel)
        edge_list = get_mass_window_edges(file_path, pt_hist_path)

        mass_hist = self.get_hist(file_key, mass_hist_name, axis_steering=mass_steering, root_hist=True,
                                  file_sel=file_sel)
        dict_list = []
        matches = re.findall(r'\[([^\]]*)\]', pt_steering)
        for index, edge in enumerate(edge_list):
            if index < len(edge_list) - 1:
                pt_steering = "dipt[" + matches[0] + "];dimass[" + matches[1] + f"{index}]"
                # nominal
                pt_hist = self.get_hist(file_key, pt_hist_name, axis_steering=pt_steering, root_hist=True,
                                        file_sel=file_sel)
                mass_hist.GetXaxis().SetRangeUser(edge, edge_list[index + 1])

                stat = get_summary_statistics(pt_hist, stat_type, prob)
                temp_dict = {"mass_window": str(edge) + ":" + str(edge_list[index + 1]),
                             "pt": stat.value, "pt_stat error": stat.error,
                             "mass": mass_hist.GetMean(), "mass_stat error": mass_hist.GetMeanError()}

                # systematics
                # function to get systematic (symmetric)
                # pt_bkg stat_error
                # pt_eff sf_error
                if sys:
                    sys_name_dict = {"pt_bkg stat": [""]}
                    sys_dict = self.get_isr_sys(stat.value, sys_name_dict, file_key, pt_hist_name, pt_steering, file_sel)
                    temp_dict.update(sys_dict)

                dict_list.append(temp_dict)

        return pd.DataFrame(dict_list)

    def get_isr_sys(self, nominal, sys_dict, file_key, hist_name, axis_steering, file_sel,
                    stat_type="mean", prob=0.5):
        stat_up_hist = self.get_hist(file_key, hist_name, axis_steering=axis_steering, root_hist=True,
                                     file_sel=file_sel, stat_variation=1)
        stat_down_hist = self.get_hist(file_key, hist_name, axis_steering=axis_steering, root_hist=True,
                                       file_sel=file_sel, stat_variation=-1)

        stat_up = get_summary_statistics(stat_up_hist, stat_type, prob).value
        stat_down = get_summary_statistics(stat_down_hist, stat_type, prob).value
        sys = abs(nominal - stat_up) if abs(nominal - stat_up) > abs(nominal - stat_down) else abs(nominal - stat_down)
        return {"pt_bkg stat": sys}

    def get_raw_labels(self, file_key, hist_name):
        hist_path = self.hist_dir + hist_name
        file_path = self.get_file_path(file_key)

        return get_raw_labels(file_path, hist_path)

    def get_file_path(self, file_key, file_sel=""):
        sample_type, hist_label, file_name = file_key.split("/")
        if "data:bkg_subtracted" == sample_type:
            sample_type = "data"
        if file_name == "all":
            file_name = [*self.input_files[file_sel][sample_type][hist_label].keys()][0]  # use one of samples
        file_path = self.get_file_directory(sample_type, file_sel) + self.input_files[file_sel][sample_type][hist_label][file_name]
        return file_path

    def set_input_files(self, path, sample_config, *args):
        self.input_files[""] = dict()

        if sample_config:
            with open(path + sample_config, 'r') as config_file:
                config_json = json.load(config_file)
                self.input_files[""]["data"] = config_json["data"]
                self.input_files[""]["mc"] = config_json["mc"]
            for arg in args:
                with open(path + sample_config, 'r') as config_file:
                    config_json = json.load(config_file)
                    self.input_files[arg] = dict()
                    self.input_files[arg]["data"] = config_json["data"]
                    self.input_files[arg]["mc"] = config_json["mc"]

    def print_files(self):
        print(self.input_files)


if __name__ == "__main__":
    print("test ROOFiles")

    data_handle = ROOTFiles("/Users/jhkim/cms_snu/data/Ultralegacy/", "ee", "2016a")
    data_handle.print_files()
    print(data_handle.get_hist("data/2016prePFV/all", "ee2016a/cutflow"))
    print(data_handle.get_raw_labels("data/2016prePFV/all", "ee2016a/cutflow"))
