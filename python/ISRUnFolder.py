import pandas as pd
import ctypes
from helper import *


class ISRUnFolder:
    def __init__(self, data_handle):

        # assuming 2D unfolding as default
        self.data_handle = data_handle
        self.signal_hist_key = ""
        self.input_hist_key = ""

        self.dipt_var_name = data_handle.dipt_var_name
        self.dimass_var_name = data_handle.dimass_var_name
        self.bins = data_handle.bins
        self.mass_windows = data_handle.mass_windows
        self.bin_info = data_handle.bin_info

        self.sys_turned_on = False

        self.unfolders = {self.dipt_var_name: None, self.dimass_var_name: None}
        self.sys_unfolders = {self.dipt_var_name: dict(), self.dimass_var_name: dict()}

        # self.axis1_label = self.bins.unfolded.GetDistributionAxisLabel(0)
        # self.axis2_label = self.bins.unfolded.GetDistributionAxisLabel(1)
        # unfold configuration as dictionary
        # unfold_config[pt_hist_name] = {}
        # unfold_config[mass_hist_name] = {}
        # systematic unfolding self.sys_unfolders
        # type 1: input_hist variation
        # type 2: matrix variation
        # type 3: both input and matrix variation
        # type 4: unfolding parameter variation

    def set_input_hist_key(self, hist_key):
        self.input_hist_key = hist_key

    def set_signal_hist_key(self, hist_key):
        self.signal_hist_key = hist_key

    def create_unfolder(self, var_name, reg_mode, constraint_area, density_mode, sys=False):

        response_matrix_name = self.bin_info.get_response_matrix_name(var_name)
        response_matrix = self.data_handle.get_hist(response_matrix_name, [self.signal_hist_key], root_hist=True)

        self.unfolders[var_name] = rt.TUnfoldDensity(response_matrix,
                                                     rt.TUnfold.kHistMapOutputHoriz,
                                                     reg_mode,
                                                     constraint_area,
                                                     density_mode,
                                                     self.bins[var_name]["unfolded"],
                                                     self.bins[var_name]["folded"])
        if sys:
            self.sys_turned_on = True
            self._create_sys_unfolders(var_name, reg_mode, constraint_area, density_mode)  # TODO define properties
            # for unfolding conditions?

    def _create_sys_unfolders(self, var_name, reg_mode, constraint_area, density_mode):

        for sys_name in self.data_handle.sys_dict:
            self.sys_unfolders[var_name][sys_name] = dict()
            for var in self.data_handle.sys_dict[sys_name][1]:

                response_matrix_name = self.bin_info.get_response_matrix_name(var_name)
                if self.data_handle.sys_dict[sys_name][0] == "hist name postfix":
                    postfix = var
                    response_matrix_name += postfix

                response_matrix = self.data_handle.get_hist(response_matrix_name, [self.signal_hist_key],
                                                            root_hist=True)

                self.sys_unfolders[var_name][sys_name][var] = rt.TUnfoldDensity(response_matrix,
                                                                                rt.TUnfold.kHistMapOutputHoriz,
                                                                                reg_mode,
                                                                                constraint_area,
                                                                                density_mode,
                                                                                self.bins[var_name]["unfolded"],
                                                                                self.bins[var_name]["folded"])

    def set_unfold_input(self, var_name, subtract_signal_fake=False):

        hist_name = self.bin_info.get_hist_name("hist", var_name, "folded")
        input_hist = self.data_handle.get_hist(hist_name, [self.input_hist_key], root_hist=True)
        self.unfolders[var_name].SetInput(input_hist)

        if subtract_signal_fake:
            hist_name = self.bin_info.get_hist_name("fake_hist", var_name, "folded")
            signal_fake_hist = self.data_handle.get_hist(hist_name, [self.signal_hist_key], root_hist=True)
            self.unfolders[var_name].SubtractBackground(signal_fake_hist, "Drell-Yan fake")

        if self.sys_turned_on:
            for sys_name in self.data_handle.sys_dict:
                for var in self.data_handle.sys_dict[sys_name][1]:
                    stat_variation = 0
                    hist_name_postfix = ""

                    if self.data_handle.sys_dict[sys_name][0] == "stat variation":
                        stat_variation = var
                    if self.data_handle.sys_dict[sys_name][0] == "hist name postfix":
                        hist_name_postfix = var

                    hist_name = self.bin_info.get_hist_name("hist", var_name, "folded")
                    input_hist = self.data_handle.get_hist(hist_name, [self.input_hist_key], root_hist=True,
                                                           stat_variation=stat_variation,
                                                           hist_name_postfix=hist_name_postfix)

                    self.sys_unfolders[var_name][sys_name][var].SetInput(input_hist)

                    if subtract_signal_fake:
                        self.sys_unfolders[var_name][sys_name][var].SubtractBackground(signal_fake_hist,
                                                                                       "Drell-Yan fake")

    def set_unfold_input_external(self, ext_hist, var_name, subtract_signal_fake=False):

        self.unfolders[var_name].SetInput(ext_hist)
        if subtract_signal_fake:
            hist_name = self.bin_info.get_hist_name("fake_hist", var_name, "folded")
            signal_fake_hist = self.data_handle.get_hist(hist_name, [self.signal_hist_key], root_hist=True)
            self.unfolders[var_name].SubtractBackground(signal_fake_hist, "Drell-Yan fake")

    def unfold(self, var_name):
        self.unfolders[var_name].DoUnfold(0)

        if self.sys_turned_on:
            for sys_name in self.data_handle.sys_dict:
                for var in self.data_handle.sys_dict[sys_name][1]:
                    self.sys_unfolders[var_name][sys_name][var].DoUnfold(0)

    def do_acceptance_correction(self, var_name, unfolded_hist, axis_steering):

        hist_name = self.bin_info.get_hist_name("fullphase_hist", var_name, "unfolded")
        dy_full_phase_hist = self.data_handle.get_hist(hist_name, [self.signal_hist_key], root_hist=True)
        dy_full_phase_hist = self.bins[var_name]["unfolded"].ExtractHistogram("extracted_" + axis_steering,
                                                                              dy_full_phase_hist, 0, True,
                                                                              axis_steering)

        hist_name = self.bin_info.get_hist_name("hist", var_name, "unfolded")
        dy_unfolded_hist = self.data_handle.get_hist(hist_name, [self.signal_hist_key], root_hist=True)
        dy_unfolded_hist = self.bins[var_name]["unfolded"].ExtractHistogram("extracted_" + axis_steering,
                                                                            dy_unfolded_hist, 0, True,
                                                                            axis_steering)
        acceptance_th1 = dy_full_phase_hist.Clone("acceptance")
        acceptance_th1.Divide(dy_unfolded_hist)
        unfolded_hist.Multiply(acceptance_th1)

        return unfolded_hist

    def get_output_hist(self, var_name, axis_steering, use_axis_bin=True, do_accept=False, root_hist=False,
                        norm=False, binwnorm=False):

        unfolded_hist = self.unfolders[var_name].GetOutput("unfolded_" + var_name,
                                                           ctypes.c_char_p(0),
                                                           ctypes.c_char_p(0),
                                                           axis_steering, use_axis_bin)
        if do_accept:
            unfolded_hist = self.do_acceptance_correction(var_name, unfolded_hist, axis_steering)

        if norm:
            integral = unfolded_hist.Integral()
            unfolded_hist.Scale(1. / integral)
        if binwnorm:
            unfolded_hist.Scale(1., "width")
            
        if root_hist:
            return unfolded_hist
        else:
            return root_to_numpy(unfolded_hist)

    def get_isr_hists(self, dimass_steering, dipt_steering, mass_index, do_accept=False):

        mass_unfolded_hist = self.unfolders[self.dimass_var_name].GetOutput("mass_unfolded",
                                                                            ctypes.c_char_p(0),
                                                                            ctypes.c_char_p(0),
                                                                            dimass_steering, True)
        mass_unfolded_hist.GetXaxis().SetRangeUser(self.mass_windows[mass_index], self.mass_windows[mass_index + 1])

        pt_unfolded_hist = self.unfolders[self.dipt_var_name].GetOutput("pt_unfolded",
                                                                        ctypes.c_char_p(0),
                                                                        ctypes.c_char_p(0),
                                                                        dipt_steering, True)

        if do_accept:
            mass_unfolded_hist = self.do_acceptance_correction(self.dimass_var_name, mass_unfolded_hist,
                                                               dimass_steering)
            pt_unfolded_hist = self.do_acceptance_correction(self.dipt_var_name, pt_unfolded_hist,
                                                             dipt_steering)

        return mass_unfolded_hist, pt_unfolded_hist

    def get_isr_sys_hists(self, dimass_steering, dipt_steering, mass_index, sys_name, sys_variation,
                          do_accept=False):

        mass_unfolded_hist = self.sys_unfolders[self.dimass_var_name][sys_name][sys_variation].GetOutput(
            "mass_unfolded",
            ctypes.c_char_p(0),
            ctypes.c_char_p(0),
            dimass_steering, True
        )
        mass_unfolded_hist.GetXaxis().SetRangeUser(self.mass_windows[mass_index], self.mass_windows[mass_index + 1])

        pt_unfolded_hist = self.sys_unfolders[self.dipt_var_name][sys_name][sys_variation].GetOutput(
            "pt_unfolded",
            ctypes.c_char_p(0),
            ctypes.c_char_p(0),
            dipt_steering, True
        )
        if do_accept:
            mass_unfolded_hist = self.do_acceptance_correction(self.dimass_var_name, mass_unfolded_hist,
                                                               dimass_steering)
            pt_unfolded_hist = self.do_acceptance_correction(self.dipt_var_name, pt_unfolded_hist,
                                                             dipt_steering)

        return mass_unfolded_hist, pt_unfolded_hist

    def get_isr_sys_dict(self, dipt_nominal, dimass_nominal, dimass_steering, dipt_steering, mass_index,
                         do_accept=False):

        dipt_sys_dict = {}
        dimass_sys_dict = {}

        for sys_name in self.data_handle.sys_dict:
            dipt_delta_sum = 0
            dimass_delta_sum = 0
            for var in self.data_handle.sys_dict[sys_name][1]:
                dimass_var, dipt_var = self.get_isr_summary(dimass_steering, dipt_steering, mass_index,
                                                            sys_name, var, do_accept=do_accept)

                dipt_delta_sum += abs(dipt_var.value-dipt_nominal)
                dimass_delta_sum += abs(dimass_var.value-dimass_nominal)

            dipt_sys_dict[sys_name] = dipt_delta_sum / len(self.data_handle.sys_dict[sys_name][1])
            dimass_sys_dict[sys_name] = dimass_delta_sum / len(self.data_handle.sys_dict[sys_name][1])

        return dimass_sys_dict, dipt_sys_dict

    def get_isr_summary(self, dimass_steering, dipt_steering, mass_index, sys_name, sys_variation,
                        do_accept=False):

        dimass_hist, dipt_hist = self.get_isr_sys_hists(dimass_steering, dipt_steering, mass_index,
                                                        sys_name, sys_variation, do_accept=do_accept)

        dipt_stat = get_summary_statistics(dipt_hist, "mean", 0.5)
        dimass_stat = get_summary_statistics(dimass_hist, "mean", 0.5)

        return dimass_stat, dipt_stat

    def get_isr_dataframe(self, mass_axis_steering="dimass[OU];dipt[UOC0]",
                          pt_axis_steering="dipt[O];dimass[UOC]",
                          stat_type="mean", prob=0.5,
                          do_accept=False, pt_cut=None):

        dipt_dict_list = []
        dimass_dict_list = []

        matches = re.findall(r'([a-z]+)\[([^]]*)]', pt_axis_steering)
        for index, edge in enumerate(self.mass_windows):
            if index < len(self.mass_windows) - 1:

                pt_steering_mass_bin = matches[0][0] + "[" + matches[0][1] + "];" + \
                                       matches[1][0] + "[" + matches[1][1] + f"{index}]"

                mass_unfolded_hist, pt_unfolded_hist = self.get_isr_hists(mass_axis_steering, pt_steering_mass_bin,
                                                                          index, do_accept=do_accept)
                if pt_cut is not None:
                    pt_unfolded_hist.GetXaxis().SetRangeUser(0, pt_cut)

                stat = get_summary_statistics(pt_unfolded_hist, stat_type, prob)
                dimass_stat = get_summary_statistics(mass_unfolded_hist, stat_type, prob)

                dipt_dict = {"mass_window": str(edge) + ":" + str(self.mass_windows[index + 1]),
                             stat_type: stat.value, "stat error": stat.error}
                dimass_dict = {"mass_window": str(edge) + ":" + str(self.mass_windows[index + 1]),
                               stat_type: dimass_stat.value, "stat error": dimass_stat.error}

                if self.sys_turned_on:
                    dimass_sys_dict, dipt_sys_dict = self.get_isr_sys_dict(stat.value, dimass_stat.value,
                                                                           mass_axis_steering, pt_steering_mass_bin,
                                                                           index, do_accept=do_accept)
                    dipt_dict.update(dipt_sys_dict)
                    dimass_dict.update(dimass_sys_dict)

                dipt_dict_list.append(dipt_dict)
                dimass_dict_list.append(dimass_dict)

        dipt_df = pd.DataFrame(dipt_dict_list)
        dimass_df = pd.DataFrame(dimass_dict_list)
        reg_expression = r".*error$"
        calculate_squared_root_sum(dimass_df, reg_expression)
        calculate_squared_root_sum(dipt_df, reg_expression)

        return dimass_df, dipt_df
