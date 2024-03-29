import json

from numpy import ndarray
from helper import *
from collections import namedtuple


# ROOT histogram to numpy, pandas etc
class ROOTFiles:  # FIXME maybe better to name this class as Experiments?

    def __init__(self, base_file_path, channel, period, *args, sample_config="sample_config.json"):

        self.file_path = {"data": base_file_path + "data/" + channel + "/" + period,
                          "mc": base_file_path + "mc/" + period}
        self.hist_path = channel + period + "/"
        self.input_files = dict()
        self.set_input_files(base_file_path, channel + "_" + sample_config, *args)
        self.data_period = period

        self.bin_width_norm = False
        self.additional_sel = ""  # nominal, systematic TODO additional_sel
        self.hist_name = ""

        self.mc_key_list = None  # TODO set default using sample_config
        self.sys_dict = {
            "bkg stat error": ("stat variation", [-1, 1], self.additional_sel),
            # "ID SF": {"hist name postfix": ("_PUweight_up", "_PUweight_down", "SYS__:")}
        }

    def set_bin_width_norm(self, bin_width_norm):
        self.bin_width_norm = bin_width_norm

    def set_additional_sel(self, additional_sel):
        self.additional_sel = additional_sel

    def set_hist_name(self, hist_name):
        self.hist_name = hist_name

    def set_mc_key_list(self, mc_key_list):
        self.mc_key_list = mc_key_list

    def set_input_files(self, path, sample_config, *args):
        with open(path + sample_config, 'r') as config_file:
            config_json = json.load(config_file)
            self.input_files["data"] = config_json["data"]
            self.input_files["mc"] = config_json["mc"]

    def get_file_directory(self, sample_type):
        file_postfix = "/"
        if self.additional_sel != "":
            file_postfix = "/" + self.additional_sel.split(":")[0] + "/"
        return self.file_path[sample_type] + file_postfix

    def get_file_path(self, hist_key):
        sample_type, hist_label, file_name = hist_key.split("/")
        if "data:bkg_subtracted" == sample_type:
            sample_type = "data"
        if file_name == "all":
            file_name = [*self.input_files[sample_type][hist_label].keys()][0]  # use one of samples
        file_path = self.get_file_directory(sample_type) + self.input_files[sample_type][hist_label][file_name]
        return file_path

    def get_obj_path(self, file_name, obj_name, obj_name_postfix=""):
        file_name_parsing = file_name.split(":")
        full_hist_dir = self.hist_path
        if self.additional_sel != "":
            full_hist_dir = self.hist_path + self.additional_sel.split(":")[1] + "/"
        if len(file_name_parsing) == 2:  # hprefix set in config file
            hist_path = full_hist_dir + file_name_parsing[1] + obj_name + obj_name_postfix  # ex) tau_
            return hist_path
        else:
            hist_path = full_hist_dir + obj_name + obj_name_postfix
            return hist_path

    def get_hist_path(self, file_name, hist_name_postfix=""):
        return self.get_obj_path(file_name, self.hist_name, hist_name_postfix)

    # get values, bin, error as numpy array
    def __get_hist(self, hist_key, norm=False, stat_variation=0, hist_name_postfix=""):
        sample_type, hist_label, file_name = hist_key.split("/")
        raw_hist = None

        if file_name == "all":  # sum all histograms in a file group
            for index, file in enumerate(self.input_files[sample_type][hist_label]):
                file_path = self.get_file_directory(sample_type) + \
                            self.input_files[sample_type][hist_label][file]
                hist_path = self.get_hist_path(file, hist_name_postfix)
                if index == 0:
                    raw_hist = get_raw_hist(file_path, hist_path, stat_variation=stat_variation)
                else:
                    raw_hist.Add(get_raw_hist(file_path, hist_path, stat_variation=stat_variation))
        else:
            file_path = self.get_file_directory(sample_type) + \
                        self.input_files[sample_type][hist_label][file_name]
            hist_path = self.get_hist_path(file_name, hist_name_postfix)
            raw_hist = get_raw_hist(file_path, hist_path, stat_variation=stat_variation)

        if norm:
            integral = raw_hist.Integral()
            raw_hist.Scale(1. / integral)

        if self.bin_width_norm:
            raw_hist.Scale(1., "width")

        return raw_hist

    def get_combined_hist(self, hist_key_list, root_hist=False, norm=False, stat_variation=0, hist_name_postfix=""):
        raw_hist = None
        if len(hist_key_list) == 1:
            sample_type, hist_label, file_name = hist_key_list[0].split("/")
            if "data:bkg_subtracted" == sample_type:
                raw_hist = self.__get_bkg_subtracted_data_hist(bkg_stat=stat_variation)
            else:
                raw_hist = self.__get_hist(hist_key_list[0], stat_variation=stat_variation)
        else:
            for index, hist_key in enumerate(hist_key_list):
                hist = self.__get_hist(hist_key, stat_variation=stat_variation,
                                       hist_name_postfix=hist_name_postfix)
                if index == 0:
                    raw_hist = hist
                else:
                    raw_hist.Add(hist)

        if norm:
            integral = raw_hist.Integral()
            raw_hist.Scale(1. / integral)

        return return_hist(raw_hist, root_hist)

    def get_hist(self, hist_name, hist_key_list, root_hist=False, norm=False, stat_variation=0, hist_name_postfix=""):
        self.set_hist_name(hist_name)
        return self.get_combined_hist(hist_key_list, root_hist, norm, stat_variation, hist_name_postfix)

    def make_stack_list(self, hist_key_list):
        values_list: list[ndarray] = []
        bins = None
        errors_list = []
        for hist_key in hist_key_list:
            hist = self.__get_hist(hist_key)
            hist = root_to_numpy(hist)
            values_list.append(hist.values)
            bins = hist.bins
            errors_list.append(hist.errors)

        Stack = namedtuple('Stack', ['values_list', 'bins', 'errors_list'])
        stack = Stack(values_list, bins, errors_list)
        return stack

    def get_expectation_hist(self, root_hist=False, norm=False, mc_stat=0, include_signal_mc=True):

        mc_key_list = [mc_key for mc_key in self.mc_key_list
                       if ("background" in mc_key) or ("signal" in mc_key and include_signal_mc)]
        mc_hist = self.get_combined_hist(mc_key_list, root_hist=True, stat_variation=mc_stat)

        if norm:
            integral = mc_hist.Integral()
            mc_hist.Scale(1. / integral)

        return return_hist(mc_hist, root_hist)

    # TODO hist_name_postfix
    def __get_bkg_subtracted_data_hist(self, norm=False, bkg_stat=0):
        data_hist = self.__get_hist("data/" + self.data_period + "/all")
        # make hist_key_list to background mc
        bkg_key_list = [bkg_key for bkg_key in self.mc_key_list if "background" in bkg_key]
        bkg_hist = self.get_combined_hist(bkg_key_list, root_hist=True, stat_variation=bkg_stat)
        bkg_subtracted_data_hist = data_hist.Clone("data_bkg_subtracted")
        bkg_subtracted_data_hist.Add(bkg_hist, -1)

        if norm:
            integral = bkg_subtracted_data_hist.Integral()
            bkg_subtracted_data_hist.Scale(1. / integral)

        return return_hist(bkg_subtracted_data_hist, root_hist=True)

    def get_scaled_hist(self, hist_key_list, root_hist=False, norm=False, systematic_dict=None):
        hist = self.get_combined_hist(hist_key_list, root_hist, norm)
        sys_error = None
        if systematic_dict is not None:
            sys_error = self.get_combined_hist_systematic_bands(systematic_dict, hist_key_list, norm)
            sys_error = np.nan_to_num(sys_error / hist.values)

        ones = np.nan_to_num(hist.values / hist.values)
        stat_error = np.nan_to_num(hist.errors / hist.values)
        return Hist(ones, hist.bins, stat_error), sys_error

    # TODO how to handle systematic error for division
    # currently considered statistical error assuming uncorrelated histograms
    def get_ratio_hist(self, file1_key_list, file2_key_list, root_hist=False, norm=False):
        combined_nominator = self.get_combined_hist(file1_key_list, root_hist=True, norm=norm)
        combined_denominator = self.get_combined_hist(file2_key_list, root_hist=True, norm=norm)
        ratio = combined_nominator.Clone("ratio")
        ratio.Divide(combined_denominator)

        return return_hist(ratio, root_hist)

    # either nominator or denominator is external TH histogram
    def get_ratio_hist2(self, hist_key_list, nominator=None,
                        denominator=None, root_hist=False, norm=False):
        # TODO error handling?
        if nominator is None and denominator is None:
            return
        if nominator is not None and denominator is not None:
            return

        if nominator is None:
            nominator = self.get_combined_hist(hist_key_list, root_hist=True, norm=norm)
            ratio = nominator.Clone("ratio")
        else:
            ratio = nominator.Clone("ratio")
            denominator = self.get_combined_hist(hist_key_list, root_hist=True, norm=norm)
        ratio.Divide(denominator)

        return return_hist(ratio, root_hist)

    def get_raw_labels(self, hist_key):
        hist_path = self.hist_path + self.hist_name
        file_path = self.get_file_path(hist_key)

        return get_raw_labels(file_path, hist_path)

    def get_hist_systematic_bands(self, systematic_dict, hist_key, norm=False):
        nominal_hist = self.__get_hist(hist_key, norm=norm)
        # symmetric systematic
        systematic_list = []
        for systematic in systematic_dict.keys():
            delta_list = []
            for variation in systematic_dict[systematic]:
                variation_postfix = ""
                if variation != "":
                    variation_postfix += "_" + variation
                hist_name_postfix = "_" + systematic + variation_postfix
                systematic_hist = self.__get_hist(hist_key, norm=norm,
                                                  hist_name_postfix=hist_name_postfix)
                systematic_hist.Add(nominal_hist, -1)
                h = root_to_numpy(systematic_hist)
                delta_list.append(np.absolute(h.values))
            systematic_list.append(np.square(np.mean(delta_list, axis=0)))
        sqrt_sum = np.sqrt(np.sum(np.array(systematic_list), axis=0))  # square root sum
        return sqrt_sum

    def get_combined_hist_systematic_bands(self, systematic_dict, hist_key_list, norm=False):
        nominal_hist = self.get_combined_hist(hist_key_list, root_hist=True, norm=norm)
        # symmetric systematic
        systematic_list = []
        for systematic in systematic_dict.keys():
            delta_list = []
            for variation in systematic_dict[systematic]:
                variation_postfix = ""
                if variation != "":
                    variation_postfix += "_" + variation
                hist_name_postfix = "_" + systematic + variation_postfix
                systematic_hist = self.get_combined_hist(hist_key_list, root_hist=True,
                                                         norm=norm, hist_name_postfix=hist_name_postfix)
                systematic_hist.Add(nominal_hist, -1)
                h = root_to_numpy(systematic_hist)
                delta_list.append(np.absolute(h.values))
            systematic_list.append(np.square(np.mean(delta_list, axis=0)))
        sqrt_sum = np.sqrt(np.sum(np.array(systematic_list), axis=0))  # square root sum
        return sqrt_sum


if __name__ == "__main__":
    print("test ROOFiles")

    data_handle = ROOTFiles("/Users/jhkim/cms_snu/data/Ultralegacy/", "ee", "2016a")
    data_handle.print_files()
