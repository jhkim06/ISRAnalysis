import ROOTFiles
import ROOT as rt
import pandas as pd
import ctypes
from helper import get_summary_statistics, root_to_numpy
import re


class ISRUnFolder:
    def __init__(self, data_handle, bin_names, pt_hist_name="dipt-dimass", mass_hist_name="dimass-dipt"):

        # assuming 2D unfolding as default
        self.data_handle = data_handle
        self.pt_hist_name = pt_hist_name
        self.mass_hist_name = mass_hist_name
        self.unfolders = {pt_hist_name: None, mass_hist_name: None}
        self.mass_windows = None
        self.bin_names = bin_names
        self.bins = {pt_hist_name: {"folded": None, "unfolded": None},
                     mass_hist_name: {"folded": None, "unfolded": None}}
        # self.axis1_label = self.bins.unfolded.GetDistributionAxisLabel(0)
        # self.axis2_label = self.bins.unfolded.GetDistributionAxisLabel(1)

    def create_unfolder(self, var_name, reg_mode, constraint_area, density_mode,
                        file_key="mc/signal:Drell-Yan/all", file_sel=""):

        response_matrix_name = self.make_response_matrix_name(var_name)
        self.response_matrix = self.data_handle.get_hist(file_key, response_matrix_name, root_hist=True, file_sel=file_sel)
        unfolded_bin = self.get_bin(var_name, "unfolded", file_key, file_sel=file_sel)
        folded_bin = self.get_bin(var_name, "folded", file_key, file_sel=file_sel)
        self.bins[var_name]["folded"] = folded_bin
        self.bins[var_name]["unfolded"] = unfolded_bin
        self.unfolders[var_name] = rt.TUnfoldDensity(self.response_matrix, rt.TUnfold.kHistMapOutputHoriz,
                                                     reg_mode,
                                                     constraint_area,
                                                     density_mode,
                                                     unfolded_bin, folded_bin)

        # set mass windows
        if "dipt-dimass" in var_name:
            edges = unfolded_bin.GetDistributionBinning(1)
            self.mass_windows = [edge for edge in edges]

    def set_unfold_input(self, var_name, file_key, subtract_signal_fake=False,
                         signal_file_key="mc/signal:Drell-Yan/all",  file_sel=""):
        axis1_var = var_name.split("(")[0].split("-")[0]
        axis2_var = var_name.split("(")[0].split("-")[1]

        hist_name = "[tunfold-hist]_" + \
                    "[" + var_name.split("(")[0]+"]_" + \
                    "[folded_" + self.bin_names[axis1_var]["folded"] + \
                    "-folded_" + self.bin_names[axis2_var]["folded"] + "]"
        input_hist = self.data_handle.get_hist(file_key, hist_name, root_hist=True, file_sel=file_sel)
        self.unfolders[var_name].SetInput(input_hist)

        if subtract_signal_fake:
            hist_name = "[tunfold-fake_hist]_" + \
                        "["+var_name+"]_[folded_" + self.bin_names[axis1_var]["folded"] + \
                        "-folded_" + self.bin_names[axis2_var]["folded"]+"]"
            signal_fake_hist = self.data_handle.get_hist(signal_file_key, hist_name, root_hist=True, file_sel=file_sel)
            self.unfolders[var_name].SubtractBackground(signal_fake_hist, "Drell-Yan fake")

        mc_hist = self.response_matrix.ProjectionY("projected_reco", 1, -1, "e")
        print(f"data/mc: {input_hist.Integral()/mc_hist.Integral()}")

    def set_unfold_input_external(self, ext_hist, var_name, subtract_signal_fake=False,
                                  signal_file_key="mc/signal:Drell-Yan/all",  file_sel=""):
        axis1_var = var_name.split("(")[0].split("-")[0]
        axis2_var = var_name.split("(")[0].split("-")[1]

        self.unfolders[var_name].SetInput(ext_hist)
        if subtract_signal_fake:
            hist_name = "[tunfold-fake_hist]_" + \
                        "["+var_name+"]_[folded_" + self.bin_names[axis1_var]["folded"] + \
                        "-folded_" + self.bin_names[axis2_var]["folded"]+"]"
            signal_fake_hist = self.data_handle.get_hist(signal_file_key, hist_name, root_hist=True, file_sel=file_sel)
            self.unfolders[var_name].SubtractBackground(signal_fake_hist, "Drell-Yan fake")

    def make_response_matrix_name(self, var_name):
        axis1_var = var_name.split("(")[0].split("-")[0]
        axis2_var = var_name.split("(")[0].split("-")[1]
        axis1_folded_bin_name = "folded_" + self.bin_names[axis1_var]["folded"]
        axis2_folded_bin_name = "folded_" + self.bin_names[axis2_var]["folded"]
        axis1_unfolded_bin_name = "unfolded_" + self.bin_names[axis1_var]["unfolded"]
        axis2_unfolded_bin_name = "unfolded_" + self.bin_names[axis2_var]["unfolded"]

        response_matrix_name = "[tunfold-matrix]_[" + var_name + "]" + "_" + \
                               "[" + axis1_folded_bin_name + "-" + axis2_folded_bin_name + "]" + "_" + \
                               "[" + axis1_unfolded_bin_name + "-" + axis2_unfolded_bin_name + "]"
        return response_matrix_name

    def get_bin(self, var_name, which_level, file_key="mc/signal:Drell-Yan/all", file_sel=""):
        axis1_var = var_name.split("(")[0].split("-")[0]
        axis2_var = var_name.split("(")[0].split("-")[1]

        sample_type, hist_label, file_name = file_key.split("/")
        if file_name == "all":
            file_name = [*self.data_handle.input_files[""][sample_type][hist_label].keys()][0]
            file_key = "/".join([sample_type, hist_label, file_name])
        bin_path = "[tunfold-bin]_[" + var_name.split("(")[0] + "]" + "_" + \
                   "[" + which_level + "_" + self.bin_names[axis1_var][which_level] + \
                   "-" + which_level + "_" + self.bin_names[axis2_var][which_level] + "]"
        unfold_bin = self.data_handle.get_hist(file_key, bin_path, root_hist=True, file_sel=file_sel)
        return unfold_bin

    def do_acceptance_correction(self, var_name, unfolded_hist, axis_steering,
                                 file_key="mc/signal:Drell-Yan/all", file_sel=""):
        axis1_var = var_name.split("(")[0].split("-")[0]
        axis2_var = var_name.split("(")[0].split("-")[1]

        hist_name = "[tunfold-fullphase_hist]_["+var_name+"]_" + \
                    "[unfolded_"+self.bin_names[axis1_var]["unfolded"] + \
                    "-unfolded_"+self.bin_names[axis2_var]["unfolded"]+"]"
        dy_full_phase_hist = self.data_handle.get_hist(file_key, hist_name, axis_steering,
                                                       root_hist=True, file_sel=file_sel)
        hist_name = "[tunfold-hist]_["+var_name+"]_" + \
                    "[unfolded_"+self.bin_names[axis1_var]["unfolded"] + \
                    "-unfolded_"+self.bin_names[axis2_var]["unfolded"]+"]"
        dy_unfolded_hist = self.data_handle.get_hist(file_key, hist_name, axis_steering,
                                                     root_hist=True, file_sel=file_sel)

        acceptance_th1 = dy_full_phase_hist.Clone("acceptance")
        acceptance_th1.Divide(dy_unfolded_hist)

        unfolded_hist.Multiply(acceptance_th1)
        return unfolded_hist

    def get_output_hist(self, var_name, axis_steering, use_axis_bin=True, do_accept=False, root_hist=False,
                        norm=False, binwnorm=False, file_key="mc/signal:Drell-Yan/all", file_sel=""):
        unfolded_hist = self.unfolders[var_name].GetOutput("unfolded_" + var_name,
                                                           ctypes.c_char_p(0),
                                                           ctypes.c_char_p(0),
                                                           axis_steering, use_axis_bin)
        if do_accept:
            unfolded_hist = self.do_acceptance_correction(var_name, unfolded_hist, axis_steering,
                                                          file_key=file_key, file_sel=file_sel)

        if norm:
            integral = unfolded_hist.Integral()
            unfolded_hist.Scale(1. / integral)
        if binwnorm:
            unfolded_hist.Scale(1., "width")
            
        if root_hist:
            return unfolded_hist
        else:
            return root_to_numpy(unfolded_hist)

    def get_isr_dataframe(self, mass_axis_steering="dimass[OU];dipt[UOC0]",
                          pt_axis_steering="dipt[O];dimass[UOC]", stat_type="mean", prob=0.5,
                          do_accept=False, file_sel="", pt_cut=None, file_key="mc/signal:Drell-Yan/all"):

        mass_unfolded_hist = self.unfolders[self.mass_hist_name].GetOutput("mass_unfolded",
                                                                           ctypes.c_char_p(0),
                                                                           ctypes.c_char_p(0),
                                                                           mass_axis_steering, True)
        if do_accept:
            mass_unfolded_hist = self.do_acceptance_correction(self.mass_hist_name, mass_unfolded_hist,
                                                               mass_axis_steering, file_sel=file_sel,
                                                               file_key=file_key)
        dict_list = []
        matches = re.findall(r'([a-z]+)\[([^\]]*)\]', pt_axis_steering)
        for index, edge in enumerate(self.mass_windows):
            if index < len(self.mass_windows) - 1:
                pt_steering_mass_bin = matches[0][0] + "[" + matches[0][1] + "];" + \
                                       matches[1][0] + "[" + matches[1][1] + f"{index}]"
                pt_unfolded_hist = self.unfolders[self.pt_hist_name].GetOutput("pt_unfolded",
                                                                               ctypes.c_char_p(0),
                                                                               ctypes.c_char_p(0),
                                                                               pt_steering_mass_bin, True)
                if do_accept:
                    if len(self.mass_windows) == 7 and index > 3:
                        pt_steering_mass_bin = matches[0][0] + "[" + matches[0][1] + "];" + \
                                               matches[1][0] + "[" + matches[1][1] + "45]"

                    pt_unfolded_hist = self.do_acceptance_correction(self.pt_hist_name, pt_unfolded_hist,
                                                                     pt_steering_mass_bin, file_sel=file_sel,
                                                                     file_key=file_key)
                if pt_cut is not None:
                    pt_unfolded_hist.GetXaxis().SetRangeUser(0, pt_cut)

                stat = get_summary_statistics(pt_unfolded_hist, stat_type, prob)

                mass_unfolded_hist.GetXaxis().SetRangeUser(edge, self.mass_windows[index + 1])
                temp_dict = {"mass_window": str(edge) + ":" + str(self.mass_windows[index + 1]),
                             "pt": stat.value, "pt_stat_error": stat.error,
                             "mass": mass_unfolded_hist.GetMean(),
                             "mass_stat_error": mass_unfolded_hist.GetMeanError()}
                dict_list.append(temp_dict)

        return pd.DataFrame(dict_list, columns=['mass_window',
                                                'pt', 'pt_stat_error',
                                                'mass', 'mass_stat_error'])