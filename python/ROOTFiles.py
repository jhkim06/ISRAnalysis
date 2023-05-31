import json

import numpy as np
from numpy import ndarray
from helper import *
import pandas as pd
from collections import namedtuple


# ROOT histogram to numpy, pandas etc
class ROOTFiles:

    def __init__(self, path, channel, period, *args, sample_config="sample_config.json"):
        self.file_dir = path+channel+"/"+period
        self.hist_dir = channel+period+"/"  # TODO case with additional sub directory?
        self.input_files = dict()
        self.set_input_files(sample_config, *args)
        self.data_period = period
        if period == "2016a":
            self.data_period = "2016prePFV"
        if period == "2016b":
            self.data_period = "2016postPFV"

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

    def get_file_directory(self, file_sel):
        file_postfix = "/"
        if file_sel != "":
            file_postfix = "/"+file_sel.split(":")[0] + "/"
        return self.file_dir+file_postfix

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
                 norm=False, binwnorm=False, file_sel=""):
        sample_type, hist_label, file_name = file_key.split("/")
        raw_hist = None

        if sample_type == "data:bkg_subtracted":
            raw_hist = self.get_bkg_subtracted_data_hist(hist_name, axis_steering, root_hist=True,
                                                         file_sel=file_sel)
        else:
            if file_name == "all":  # sum all histograms
                for index, file in enumerate(self.input_files[file_sel][sample_type][hist_label]):
                    file_path = self.get_file_directory(file_sel)+self.input_files[file_sel][sample_type][hist_label][file]
                    hist_path = self.get_hist_path(file, hist_name, file_sel)
                    if index == 0:
                        raw_hist = get_raw_hist(file_path, hist_path, axis_steering)
                    else:
                        raw_hist.Add(get_raw_hist(file_path, hist_path, axis_steering))
            else:
                hist_path = self.get_hist_path(file_name, hist_name, file_sel)
                file_path = self.get_file_directory(file_sel)+self.input_files[file_sel][sample_type][hist_label][file_name]
                raw_hist = get_raw_hist(file_path, hist_path, axis_steering)

        if norm:
            integral = raw_hist.Integral()
            raw_hist.Scale(1./integral)
        if binwnorm:
            raw_hist.Scale(1., "width")

        if root_hist:
            return raw_hist
        else:
            return root_to_numpy(raw_hist)

    def get_combined_hist(self, file_key_list, hist_name, axis_steering="", root_hist=False,
                          norm=False, binwnorm=False, file_sel=""):
        raw_hist = None
        for index, file_key in enumerate(file_key_list):
            hist = self.get_hist(file_key, hist_name, axis_steering, root_hist=True, file_sel=file_sel)
            if index == 0:
                raw_hist = hist
            else:
                raw_hist.Add(hist)
        if norm:
            integral = raw_hist.Integral()
            raw_hist.Scale(1./integral)
        if binwnorm:
            raw_hist.Scale(1., "width")

        if root_hist:
            return raw_hist
        else:
            return root_to_numpy(raw_hist)

    def get_bkg_subtracted_data_hist(self, hist_name, axis_steering="", root_hist=False,
                                     norm=False, binwnorm=False, file_sel=""):
        data_hist = self.get_hist("data/"+self.data_period+"/all", hist_name, axis_steering, root_hist=True,
                                  file_sel=file_sel)
        # make file_key_list to background mc
        bkg_key_list = ["mc/"+bkg_key+"/all" for bkg_key in self.input_files[""]["mc"].keys()
                        if "background" in bkg_key]
        bkg_hist = self.get_combined_hist(bkg_key_list, hist_name, axis_steering, root_hist=True,
                                          file_sel=file_sel)
        bkg_subtracted_data_hist = data_hist.Clone("data_bkg_subtracted")
        bkg_subtracted_data_hist.Add(bkg_hist, -1)
        if norm:
            integral = bkg_subtracted_data_hist.Integral()
            bkg_subtracted_data_hist.Scale(1./integral)
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
                    variation_postfix += "_"+variation
                full_hist_name = hist_name+"_"+systematic+variation_postfix
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
                    variation_postfix += "_"+variation
                full_hist_name = hist_name+"_"+systematic+variation_postfix
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
                                      norm, binwnorm)
        sys_error = self.get_combined_hist_systematic_bands(systematic_dict, file_key_list, hist_name, axis_steering,
                                           norm, binwnorm, file_sel)

        ones = np.nan_to_num(hist.values/hist.values)
        stat_error = np.nan_to_num(hist.errors/hist.values)
        return Hist(ones, hist.bins, stat_error), np.nan_to_num(sys_error/hist.values)

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

    # get ISR mean DataFrame from 2D pt-mass histogram
    def get_isr_dataframe(self, file_key, pt_hist_name, mass_hist_name, file_sel="", stat_type="mean", prob=0.5,
                          pt_steering="dipt[];dimass[UOC]", mass_steering="dimass[UO];dipt[OC0]"):
        file_path = self.get_file_path(file_key, file_sel=file_sel)
        pt_hist_path = self.hist_dir + pt_hist_name
        edge_list = get_mass_window_edges(file_path, pt_hist_path)

        mass_hist = self.get_hist(file_key, mass_hist_name, axis_steering=mass_steering, root_hist=True,
                                  file_sel=file_sel)

        dict_list = []
        matches = re.findall(r'\[([^\]]*)\]', pt_steering)
        for index, edge in enumerate(edge_list):
            if index < len(edge_list)-1:
                pt_steering = "dipt["+matches[0]+"];dimass["+matches[1]+f"{index}]"

                pt_hist = self.get_hist(file_key, pt_hist_name, axis_steering=pt_steering, root_hist=True,
                                        file_sel=file_sel)
                mass_hist.GetXaxis().SetRangeUser(edge, edge_list[index + 1])

                stat = get_summary_statistics(pt_hist, stat_type, prob)
                temp_dict = {"mass_window": str(edge)+":"+str(edge_list[index + 1]),
                             "pt": stat.value, "pt_stat_error": stat.error,
                             "mass": mass_hist.GetMean(), "mass_stat_error": mass_hist.GetMeanError()}
                dict_list.append(temp_dict)

        return pd.DataFrame(dict_list, columns=['mass_window',
                                                'pt', 'pt_stat_error',
                                                'mass', 'mass_stat_error'])

    def get_raw_labels(self, file_key, hist_name):
        hist_path = self.hist_dir+hist_name
        file_path = self.get_file_path(file_key)

        return get_raw_labels(file_path, hist_path)

    def get_file_path(self, file_key, file_sel=""):
        sample_type, hist_label, file_name = file_key.split("/")
        if "data:bkg_subtracted" == sample_type:
            sample_type = "data"
        if file_name == "all":
            file_name = [*self.input_files[file_sel][sample_type][hist_label].keys()][0]  # use one of samples
        file_path = self.get_file_directory(file_sel)+self.input_files[file_sel][sample_type][hist_label][file_name]
        return file_path

    def set_input_files(self, sample_config, *args):
        self.input_files[""] = dict()

        if sample_config:
            with open(self.file_dir+"/"+sample_config, 'r') as config_file:
                config_json = json.load(config_file)
                self.input_files[""]["data"] = config_json["data"]
                self.input_files[""]["mc"] = config_json["mc"]
            for arg in args:
                with open(self.file_dir + "/" + sample_config, 'r') as config_file:
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
