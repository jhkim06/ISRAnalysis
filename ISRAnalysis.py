import os
import sys
from array import array
import ROOT as rt
import math

import pyScripts.unfoldUtil as unfoldutil
import pandas as pd

class HistType :

    unf_input_hist = 0
    unf_bkg_hist = 1
    unf_res_matrix = 2

class ISRAnalysis:
    
    def __init__(self, 
                 unfold_name_ = "Detector", 
                 year_ = "2016", 
                 channel_= "ee", 
                 channel_postfix_ = "", 
                 matrix_filekey_ = "DY",
                 bin_def_ = "", 
                 var_ = "Pt", 
                 regMode_ = 0, 
                 bias_ = 0., 
                 density_ = 0):
        
        self.unfold_name = unfold_name_
        self.channel = channel_
        self.year = year_
        self.matrix_filekey = matrix_filekey_ 
        self.mode = regMode_ 
        self.bin_def = bin_def_
        self.dir_path = channel_ + year_  

        folded_bin_name = self.bin_def[0] 
        unfolded_bin_name = self.bin_def[1]
        
        self.var = var_
        self.n_mass_bins = None
        
        self.out_dir_path     = "output/UltraLegacy/" +self.year+"/"+self.channel+"/"
        self.in_hist_path_txt = "inFiles/UltraLegacy/" +self.year+"/"+self.channel+"/fhist.txt"

        if channel_postfix_ != "" :
            self.out_dir_path     = "output/UltraLegacy/"+self.year+"/"+self.channel+"_"+channel_postfix_+"/"
            self.in_hist_path_txt = "inFiles/UltraLegacy/"+self.year+"/"+self.channel+"_"+channel_postfix_+"/fhist.txt"
    
        # make output directory
        if not os.path.exists(self.out_dir_path):
            os.makedirs(self.out_dir_path)
    
        # read text file including root file paths for unfolding
        self.file_paths = open(self.in_hist_path_txt, 'r')
        self.input_root_file = {}
        self.input_root_file["nominal"] = {}
        self.input_root_file["SYS"] = {}
        
        for path in self.file_paths:
            path = path.lstrip(' ').rstrip(' ').rstrip('\n')
            self.input_root_file[path.split()[0]][path.split()[1]] = path.split()[2]
        
        # Unfolding configuration
        self.bias = bias_
        self.density = density_
        
        # Create ISRUnfold object
        print("Create ISRUnfold objects...")

        # Create ISRUnfold object
        self.unfold  = rt.ISRUnfold(self.unfold_name, 
                                    self.channel, 
                                    self.year, 
                                    self.out_dir_path,   
                                    False, 
                                    self.var, 
                                    folded_bin_name,
                                    unfolded_bin_name,
                                    self.mode, 
                                    self.density)

        self.unfold.setBias(self.bias)
        
        # Set response matrix
        self.unfold.setNominalRM(self.input_root_file["nominal"][self.matrix_filekey], self.dir_path)

    def set_input_hist(self, use_MC_input = False, 
                    sys_type = "Type_0", 
                    sys_name = "", 
                    input_bin_def = None, 
                    useAccept=True) :
        
        hist_filekey_temp = "Data"

        nominal_or_SYS = "nominal"
        if sys_type == "Type_1" :
            nominal_or_SYS = "SYS"

        hist_postfix = self.make_sys_hist_postfix(HistType.unf_input_hist, sys_type, sys_name)

        self.unfold.setUnfInput(self.input_root_file[nominal_or_SYS][hist_filekey_temp], 
                                self.dir_path, 
                                sys_type, sys_name,
                                hist_postfix)

    def subFake(self, dir_name = "Detector_DY_Fake", sys_type = "Type_0", sys_name = "", isFSR = False):
            
        fakeList = ["DYJets", self.dy10to50HistName]
        hist_filekey_temp = "matrix" # FIXME save DY fake in the histogram file

        if isFSR :
            hist_filekey_temp = "fsr_matrix"

        histPostfix_temp = self.make_sys_hist_postfix(HistType.unf_bkg_hist, sys_type, sys_name)  
        #print(histPostfix_temp)

        for fake_name in fakeList:

            self.unfold.subBkgs(self.input_root_file[hist_filekey_temp], dir_name, fake_name, sys_type, sys_name, histPostfix_temp)
        
    def set_bkgs(self, sys_type = "Type_0", sys_name = ""):
   
        if self.channel == "ee" : 
            bkg_list = ["DY_tau", "ttbar", "top", "antitop", "ww", "wz", "zz"] 
        else :
            bkg_list = ["DY_tau", "ttbar", "top", "antitop", "ww", "wz", "zz"] 

        hist_postfix = self.make_sys_hist_postfix(HistType.unf_bkg_hist, sys_type, sys_name) 
        
        nominal_or_SYS = "nominal"
        if sys_type == "Type_1" or sys_type == "Type_2" or sys_type == "Type_3" :
            nominal_or_SYS = "SYS"

        for bkg_name in bkg_list:
            hist_filekey_temp = "hist"
            if "SingleTop" in bkg_name and "AlphaS" in sys_name : # FIXME
                hist_postfix = self.make_sys_hist_postfix(HistType.unf_bkg_hist, "Type_3", sys_name) 

            self.unfold.subBkgs(self.input_root_file[nominal_or_SYS][bkg_name], self.dir_path, 
                                bkg_name, 
                                sys_type, sys_name, hist_postfix)
            
    def set_systematics(self, sys_type, sys_name, isFSR = False):

        self.unfold.set_systematics(sys_name)

        matrix_filekey = self.matrix_filekey
        hist_postfix = self.make_sys_hist_postfix(HistType.unf_res_matrix, sys_type, sys_name)

        if sys_name == "UnfoldingModel" :
            matrix_filekey= "matrix_mg"
        if "fsrPHOTOS" in sys_name :
            matrix_filekey = "fsr_matrix_pythia"
        if "fsrPYTHIA" in sys_name :
            matrix_filekey = "fsr_matrix_photos"

        self.unfold.setSystematicRM(self.input_root_file["SYS"][matrix_filekey], 
                                    self.dir_path, 
                                    sys_name, 
                                    hist_postfix)

    def make_sys_hist_postfix(self, hist_type, sys_type, sys_name) :

        if hist_type == HistType.unf_input_hist : # Input histogram

            if sys_type == "Type_1":
                sys_hist_postfix = "_" + sys_name
            else :
                sys_hist_postfix = ""

            return sys_hist_postfix 

        elif hist_type == HistType.unf_bkg_hist : # Background histogram

            if sys_type == "Type_1" or sys_type == "Type_2":
                sys_hist_postfix = "_" + sys_name
            else :
                sys_hist_postfix = ""
            return sys_hist_postfix

        elif hist_type == HistType.unf_res_matrix : # Response matrix
            
            if sys_type == "Type_1" or sys_type == "Type_2" or sys_type == "Type_4" :
                sys_hist_postfix = "_" + sys_name
                if sys_name == "UnfoldingModel" or "fsr" in sys_name:
                    sys_hist_postfix = ""
            else :
                sys_hist_postfix = "" 

            return sys_hist_postfix

        else :
            print("Which histogram do you mean?")

    # Do unfold! 
    def do_unfold(self, partialReg=True):

        self.unfold.doISRUnfold(partialReg)

    def doAcceptance(self, isFSR = False) :

        if isFSR :
            self.unfold.doAcceptCorr(self.input_root_file["nominal"]['hist_accept_fullPhase'], self.dir_path)
         
        else : 
            self.unfold.doAcceptCorr(self.input_root_file["nominal"]['DY'], self.dir_path)

    def combineOutFiles(self, prefix = "") :

        pt_output_path = self.out_dir_path + self.unfold_name + "_" + self.channel + "_" + self.year + "_Pt.root" 
        mass_output_path = self.out_dir_path + self.unfold_name + "_" + self.channel + "_" + self.year + "_Mass.root" 
        target_file_path = self.out_dir_path + self.unfold_name + "_" + self.channel + "_" + self.year + ".root"   
        if prefix != "" :
            target_file_path = self.out_dir_path + self.unfold_name + "_" + self.channel + "_" + self.year + "_" + prefix + ".root"
        os.system('hadd -f ' + target_file_path + " " + pt_output_path + " " + mass_output_path)

    def closeOutFiles(self) :

        self.unfold.closeOutFile()
