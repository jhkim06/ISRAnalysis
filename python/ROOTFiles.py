import ROOT as rt
import numpy as np
import json
from helper import root_to_numpy, get_raw_hist, get_raw_labels


# ROOT histogram to numpy, pandas etc
class ROOTFiles:

    def __init__(self, file_path, sample_config=None, data=None, mc=None):
        self.path = file_path
        self.input_files = dict()
        self.set_input_files(sample_config, data, mc)

    def make_stack_list(self, stack_list):
        pass

    # get values, bin, error as numpy array
    def get_hist(self, file_key, hist_name):
        sample_type, hist_label, file_name = file_key.split("/")
        raw_hist = None
        if file_name == "all":  # sum of all histograms
            for index, file in enumerate(self.input_files[sample_type][hist_label]):
                file_path = self.path+"/"+self.input_files[sample_type][hist_label][file]
                if index == 0:
                    raw_hist = get_raw_hist(file_path, hist_name)
                else:
                    raw_hist.Add(get_raw_hist(file_path, hist_name))
        else:
            raw_hist = get_raw_hist(self.path+"/"+self.input_files[sample_type][hist_label][file_name], hist_name)

        return root_to_numpy(raw_hist)

    def get_raw_labels(self, file_key, hist_name):
        sample_type, hist_label, file_name = file_key.split("/")
        if file_name == "all":
            file_name = [*self.input_files[sample_type][hist_label].keys()][0]
        file_path = self.path+"/"+self.input_files[sample_type][hist_label][file_name]

        return get_raw_labels(file_path, hist_name)

    def set_input_files(self, sample_config, data, mc):
        self.input_files["data"] = dict()
        self.input_files["mc"] = dict()

        if sample_config:
            with open(self.path+"/"+sample_config, 'r') as config_file:
                config_json = json.load(config_file)
                self.input_files["data"] = config_json["data"]
                self.input_files["mc"] = config_json["mc"]

    def print_files(self):
        print(self.input_files)


if __name__ == "__main__":
    print("ROOFiles")
    channel = "ee"
    period = "2016a"
    path = "../../data/Ultralegacy/"+channel+"/"+period

    data_handle = ROOTFiles(path, "sample_config.json")
    data_handle.print_files()
    print(data_handle.get_hist("data/2016prePFV/all", "ee2016a/cutflow"))
    print(data_handle.get_raw_labels("data/2016prePFV/all", "ee2016a/cutflow"))
