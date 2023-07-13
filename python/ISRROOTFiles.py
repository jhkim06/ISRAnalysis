from ROOTFiles import ROOTFiles

from numpy import ndarray
from helper import *
import pandas as pd
from collections import namedtuple


class ISRROOTFiles(ROOTFiles):
    def __init__(self, base_file_path, channel, period, *args, sample_config="sample_config.json"):
        super().__init__(base_file_path, channel, period, *args, sample_config=sample_config)

        self.dipt_hist_names = None
        self.dimass_hist_names = ""

        self.dipt_bin_name = ""
        self.dimass_bin_name = ""

        self.dipt_axis_steering = ""
        self.dimass_axis_steering = ""

        self.dimass_window = []
        self.dipt_cut = None

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
        print(matches)
        bin_name = "[tunfold-bin]_" + "[" + matches[1] + "]_" + "[" + matches[2] + "]"

        unfold_bin = self.get_tunfold_bin(bin_name)
        edges = unfold_bin.GetDistributionBinning(1)

        for edge in edges:
            self.dimass_window.append(edge)

    # get ISR mean DataFrame from 2D pt-mass histogram
    def get_isr_dataframe(self, file_key_list, stat_type="mean", prob=0.5):

        dipt_matches = re.findall(r'\[([^]]*)]', self.dipt_axis_steering)

        dipt_dict_list = []
        dimass_dict_list = []
        for index, edge in enumerate(self.dimass_window):
            if index < len(self.dimass_window) - 1:
                # for dipt update axis steering for each mass window

                dipt_steering = "dipt[" + dipt_matches[0] + "];dimass[" + dipt_matches[1] + f"{index}]"
                super().set_hist_name(self.dipt_hist_names[0])
                super().set_axis_steering(dipt_steering)
                dipt_hist = super().get_combined_hist(file_key_list, root_hist=True)

                super().set_hist_name(self.dimass_hist_names)
                super().set_axis_steering(self.dimass_axis_steering)
                dimass_hist = super().get_combined_hist(file_key_list, root_hist=True)
                dimass_hist.GetXaxis().SetRangeUser(edge, self.dimass_window[index + 1])

                dipt_stat = get_summary_statistics(dipt_hist, stat_type, prob)
                dimass_stat = get_summary_statistics(dimass_hist, stat_type, prob)
                dipt_dict = {"mass_window": str(edge) + ":" + str(self.dimass_window[index + 1]),
                             stat_type: dipt_stat.value, "stat error": dipt_stat.error}
                dimass_dict = {"mass_window": str(edge) + ":" + str(self.dimass_window[index + 1]),
                             stat_type: dimass_stat.value, "stat error": dimass_stat.error}

                dipt_dict_list.append(dipt_dict)
                dimass_dict_list.append(dimass_dict)

        return pd.DataFrame(dipt_dict_list), pd.DataFrame(dimass_dict_list)
