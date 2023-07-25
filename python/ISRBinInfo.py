class ISRBinInfo:
    def __init__(self, bin_info):

        # di-lepton transverse momentum
        self.bin_info = dict()
        self._set_bin_info("dipt", bin_info)
        self._set_bin_info("dimass", bin_info)

    def _set_bin_info(self, var_name, bin_info):

        self.bin_info[var_name] = dict()
        self.bin_info[var_name]["bin"] = bin_info[var_name][0]
        self.bin_info[var_name]["folded"] = bin_info[var_name][1][0]
        self.bin_info[var_name]["unfolded"] = bin_info[var_name][1][1]
        self.bin_info[var_name]["as_second_bin"] = bin_info[var_name][1][2]
        self.bin_info[var_name]["under_overflow_info"] = bin_info[var_name][1][3]

    def get_bin_name(self, var_name, level):

        axis1_var = var_name.split("-")[0]
        axis2_var = var_name.split("-")[1]

        bin_name = "[tunfold-bin]_[" + var_name + "]" + "_" + \
                   "[" + self.bin_info[axis1_var][level] + "_" + self.bin_info[axis1_var]["bin"] + \
                   "-" + self.bin_info[axis2_var][level] + "_" + self.bin_info[axis2_var]["as_second_bin"] + "]"

        return bin_name

    def get_hist_name(self, hist_type, var_name, level):

        axis1_var = var_name.split("-")[0]
        axis2_var = var_name.split("-")[1]

        return "[tunfold-" + hist_type + "]_" + "[" + var_name + "]_" +\
               "[" + self.bin_info[axis1_var][level] + "_" + self.bin_info[axis1_var]["bin"] +\
               "-" + self.bin_info[axis2_var][level] + "_" + self.bin_info[axis2_var]["as_second_bin"] + "]"

    def get_response_matrix_name(self, var_name):

        axis1_var = var_name.split("-")[0]
        axis2_var = var_name.split("-")[1]

        axis1_folded_bin_name = self.bin_info[axis1_var]["folded"] + "_" + self.bin_info[axis1_var]["bin"]
        axis2_folded_bin_name = self.bin_info[axis2_var]["folded"] + "_" + self.bin_info[axis2_var]["as_second_bin"]
        axis1_unfolded_bin_name = self.bin_info[axis1_var]["unfolded"] + "_" + self.bin_info[axis1_var]["bin"]
        axis2_unfolded_bin_name = self.bin_info[axis2_var]["unfolded"] + "_" + self.bin_info[axis2_var]["as_second_bin"]

        response_matrix_name = "[tunfold-matrix]_[" + var_name + "]" + "_" + \
                               "[" + axis1_folded_bin_name + "-" + axis2_folded_bin_name + "]" + "_" + \
                               "[" + axis1_unfolded_bin_name + "-" + axis2_unfolded_bin_name + "]"

        return response_matrix_name
