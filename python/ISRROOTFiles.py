from ROOTFiles import ROOTFiles

from helper import *
import pandas as pd
# from numpy import ndarray
# from collections import namedtuple


class ISRROOTFiles(ROOTFiles):
    def __init__(self, base_file_path, channel, period, *args, sample_config="sample_config.json"):
        super().__init__(base_file_path, channel, period, *args, sample_config=sample_config)

        self.file_key_list = None

        self.dipt_hist_names = None
        self.dimass_hist_names = ""

        self.dipt_bin_name = ""
        self.dimass_bin_name = ""

        self.dipt_axis_steering = ""
        self.dimass_axis_steering = ""

        self.dimass_window = []
        self.dipt_cut = None

        self.isr_stat_summary_type = "mean"
        self.isr_quantile_prob = 0.5

    def set_file_key_list(self, file_key_list):
        self.file_key_list = file_key_list

    def set_dipt_hist(self, hist_names, axis_steering=""):
        self.dipt_hist_names = hist_names
        self.dipt_axis_steering = axis_steering

    def set_dimass_hist(self, hist_name, axis_steering=""):
        self.dimass_hist_names = hist_name
        self.dimass_axis_steering = axis_steering

    def reset(self):
        self.set_isr_mass_window()

    def get_tunfold_bin(self, bin_name, file_key=""):
        if file_key == "":
            for mc_key in self.input_files["mc"]:
                if "signal" in mc_key:
                    file_key = "mc/" + mc_key + "/" + [*self.input_files["mc"][mc_key].keys()][0]
                    break

        sample_type, hist_label, file_name = file_key.split("/")
        file_path = super().get_file_path(file_key)
        bin_path = super().get_obj_path(file_name, bin_name)

        file = rt.TFile.Open(file_path)
        tunfold_bin = file.Get(bin_path).Clone("tunfold_bin")
        file.Close()

        return tunfold_bin

    def set_isr_mass_window(self):
        # todo error handle
        matches = re.findall(r'\[([^]]*)]', self.dipt_hist_names[0])
        bin_name = "[tunfold-bin]_" + "[" + matches[1] + "]_" + "[" + matches[2] + "]"

        unfold_bin = self.get_tunfold_bin(bin_name)
        edges = unfold_bin.GetDistributionBinning(1)

        for edge in edges:
            self.dimass_window.append(edge)

    def get_isr_hists(self, dipt_steering, dimass_steering, mass_index,
                      stat_variation=0, hist_name_postfix=""):

        # stat_variation=0, hist_name_postfix=""
        super().set_hist_name(self.dipt_hist_names[0])  # hist_postfix for systematic?
        super().set_axis_steering(dipt_steering)
        dipt_hist = super().get_combined_hist(self.file_key_list, root_hist=True, stat_variation=stat_variation,
                                              hist_name_postfix=hist_name_postfix)

        super().set_hist_name(self.dimass_hist_names)
        super().set_axis_steering(dimass_steering)
        dimass_hist = super().get_combined_hist(self.file_key_list, root_hist=True, stat_variation=stat_variation,
                                                hist_name_postfix=hist_name_postfix)
        dimass_hist.GetXaxis().SetRangeUser(self.dimass_window[mass_index], self.dimass_window[mass_index + 1])

        return dipt_hist, dimass_hist

    def get_isr_sys_dict(self, dipt_nominal, dimass_nominal, dipt_steering, dimass_steering, mass_index):

        dipt_sys_dict = {}
        dimass_sys_dict = {}
        for sys_name in self.sys_dict:
            dipt_delta_sum = 0
            dimass_delta_sum = 0
            for var in self.sys_dict[sys_name][1]:
                stat_variation = 0
                hist_name_postfix = ""

                if self.sys_dict[sys_name][0] == "stat variation":
                    stat_variation = var
                if self.sys_dict[sys_name][0] == "hist name postfix":
                    hist_name_postfix = var

                dimass_var, dipt_var = self.get_isr_summary(dipt_steering, dimass_steering, mass_index,
                                                            stat_variation, hist_name_postfix)
                dipt_delta_sum += abs(dipt_var.value-dipt_nominal)
                dimass_delta_sum += abs(dimass_var.value-dimass_nominal)
            dipt_sys_dict[sys_name] = dipt_delta_sum / len(self.sys_dict[sys_name][1])
            dimass_sys_dict[sys_name] = dimass_delta_sum / len(self.sys_dict[sys_name][1])

        return dimass_sys_dict, dipt_sys_dict

    def get_isr_summary(self, dipt_steering, dimass_steering, mass_index, stat_variation=0, hist_name_postfix=""):

        dipt_hist, dimass_hist = self.get_isr_hists(dipt_steering, dimass_steering, mass_index,
                                                    stat_variation, hist_name_postfix)

        dipt_stat = get_summary_statistics(dipt_hist, self.isr_stat_summary_type, self.isr_quantile_prob)
        dimass_stat = get_summary_statistics(dimass_hist, self.isr_stat_summary_type, self.isr_quantile_prob)

        return dimass_stat, dipt_stat

    # get ISR mean DataFrame from 2D pt-mass histogram
    def get_isr_dataframe(self, sys=False):

        dipt_matches = re.findall(r'\[([^]]*)]', self.dipt_axis_steering)

        dipt_dict_list = []
        dimass_dict_list = []
        for index, edge in enumerate(self.dimass_window):
            if index < len(self.dimass_window) - 1:
                dipt_steering = "dipt[" + dipt_matches[0] + "];dimass[" + dipt_matches[1] + f"{index}]"
                dimass_stat, dipt_stat = self.get_isr_summary(dipt_steering, self.dimass_axis_steering, index)

                dipt_dict = {"mass_window": str(edge) + ":" + str(self.dimass_window[index + 1]),
                             self.isr_stat_summary_type: dipt_stat.value, "stat error": dipt_stat.error}
                dimass_dict = {"mass_window": str(edge) + ":" + str(self.dimass_window[index + 1]),
                               self.isr_stat_summary_type: dimass_stat.value, "stat error": dimass_stat.error}

                if sys:
                    dimass_sys_dict, dipt_sys_dict = self.get_isr_sys_dict(dipt_stat.value, dimass_stat.value,
                                                                           dipt_steering, self.dimass_axis_steering,
                                                                           index)
                    dipt_dict.update(dipt_sys_dict)
                    dimass_dict.update(dimass_sys_dict)

                # systematics!
                dipt_dict_list.append(dipt_dict)
                dimass_dict_list.append(dimass_dict)

        dimass_df = pd.DataFrame(dimass_dict_list)
        dipt_df = pd.DataFrame(dipt_dict_list)
        reg_expression = r".*error$"
        calculate_squared_root_sum(dimass_df, reg_expression)
        calculate_squared_root_sum(dipt_df, reg_expression)
        return dimass_df, dipt_df
