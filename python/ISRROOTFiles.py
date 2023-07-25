from ROOTFiles import ROOTFiles
from ISRBinInfo import ISRBinInfo
from helper import *
import pandas as pd


class ISRROOTFiles(ROOTFiles):
    def __init__(self, base_file_path, channel, period, *args, sample_config="sample_config.json"):
        super().__init__(base_file_path, channel, period, *args, sample_config=sample_config)

        self.hist_key_list = None

        self.dipt_hist_name = ""
        self.dimass_hist_name = ""

        self.dipt_var_name = ""
        self.dimass_var_name = ""
        self.bin_info = None
        self.bins = None
        self.mass_windows = None

        self.isr_stat_summary_type = "mean"
        self.isr_quantile_prob = 0.5

    def set_hist_key_list(self, hist_key_list):
        self.hist_key_list = hist_key_list

    def set_bins(self):
        self.bins = {self.dipt_var_name: {"folded": None, "unfolded": None},
                     self.dimass_var_name: {"folded": None, "unfolded": None}}
        self.set_bin("dipt-dimass", "folded")
        self.set_bin("dipt-dimass", "unfolded")
        self.set_bin("dimass-dipt", "folded")
        self.set_bin("dimass-dipt", "unfolded")

    def set_bin(self, var_name, which_level):

        bin_name = self.bin_info.get_bin_name(var_name, which_level)
        unfold_bin = self.get_tunfold_bin(bin_name)
        self.bins[var_name][which_level] = unfold_bin

    def set_mass_window(self):
        edges = self.bins["dipt-dimass"]["unfolded"].GetDistributionBinning(1)
        self.mass_windows = [edge for edge in edges]

    def init_isr_hist_names(self, bin_info,
                            dipt_name="dipt-dimass", dimass_name="dimass-dipt"):

        self.dipt_var_name = dipt_name
        self.dimass_var_name = dimass_name

        self.bin_info = ISRBinInfo(bin_info)
        self.set_bins()
        self.set_mass_window()

        self.dipt_hist_name = self.bin_info.get_hist_name("hist", self.dipt_var_name, "folded")
        self.dimass_hist_name = self.bin_info.get_hist_name("hist", self.dimass_var_name, "folded")

    def get_tunfold_bin(self, bin_name, hist_key=""):
        if hist_key == "":
            for mc_key in self.input_files["mc"]:
                if "signal" in mc_key:
                    hist_key = "mc/" + mc_key + "/" + [*self.input_files["mc"][mc_key].keys()][0]
                    break

        sample_type, hist_label, file_name = hist_key.split("/")
        file_path = super().get_file_path(hist_key)
        bin_path = super().get_obj_path(file_name, bin_name)

        file = rt.TFile.Open(file_path)
        tunfold_bin = file.Get(bin_path).Clone("tunfold_bin")
        file.Close()

        return tunfold_bin

    def get_isr_hists(self, dipt_steering, dimass_steering, mass_index,
                      stat_variation=0, hist_name_postfix=""):

        # stat_variation=0, hist_name_postfix=""
        # get_hist(hist_name, hist_key_list, )
        dipt_hist = super().get_hist(self.dipt_hist_name, self.hist_key_list, root_hist=True,
                                     stat_variation=stat_variation,
                                     hist_name_postfix=hist_name_postfix)
        dipt_hist = self.bins[self.dipt_var_name]["folded"].ExtractHistogram("extracted_" + dipt_steering,
                                                                             dipt_hist, 0, True, dipt_steering)

        dimass_hist = super().get_hist(self.dimass_hist_name, self.hist_key_list, root_hist=True,
                                       stat_variation=stat_variation,
                                       hist_name_postfix=hist_name_postfix)

        dimass_hist = self.bins[self.dimass_var_name]["folded"].ExtractHistogram("extracted_" + dimass_steering,
                                                                                 dimass_hist, 0, True, dimass_steering)
        dimass_hist.GetXaxis().SetRangeUser(self.mass_windows[mass_index], self.mass_windows[mass_index + 1])

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
    def get_isr_dataframe(self, dimass_axis_steering="dimass[OU];dipt[OC0]",
                          dipt_axis_steering="dipt[O];dimass[UOC]",
                          sys=False):

        matches = re.findall(r'([a-z]+)\[([^]]*)]', dipt_axis_steering)

        dipt_dict_list = []
        dimass_dict_list = []
        for index, edge in enumerate(self.mass_windows):
            if index < len(self.mass_windows) - 1:

                dipt_steering = matches[0][0] + "[" + matches[0][1] + "];" + \
                                       matches[1][0] + "[" + matches[1][1] + f"{index}]"

                dimass_stat, dipt_stat = self.get_isr_summary(dipt_steering, dimass_axis_steering, index)

                dipt_dict = {"mass_window": str(edge) + ":" + str(self.mass_windows[index + 1]),
                             self.isr_stat_summary_type: dipt_stat.value, "stat error": dipt_stat.error}
                dimass_dict = {"mass_window": str(edge) + ":" + str(self.mass_windows[index + 1]),
                               self.isr_stat_summary_type: dimass_stat.value, "stat error": dimass_stat.error}

                if sys:
                    dimass_sys_dict, dipt_sys_dict = self.get_isr_sys_dict(dipt_stat.value, dimass_stat.value,
                                                                           dipt_steering, dimass_axis_steering,
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
