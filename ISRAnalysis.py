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
                 regMode_ = 0, 
                 matrix_filekey_ = "DY",
                 bin_def_ = "", 
                 channel_postfix_ = "", 
                 var_ = "Pt", 
                 bias_ = 0., 
                 density_ = 0):
        
        # Initialize some variables 
        self.unfold_name = unfold_name_
        self.year = year_
        self.channel = channel_
        self.mode = regMode_ 
        self.matrix_filekey = matrix_filekey_ 
        self.bin_def = bin_def_
        self.dir_path = channel_ + year_  
        
        self.var = var_
        self.n_mass_bins = None
        
        self.out_dir_path    = "output/" +self.year+"/"+self.channel+"/"
        self.in_hist_path_txt = "inFiles/" +self.year+"/"+self.channel+"/fhist.txt"

        if channel_postfix_ != "" :
            self.out_dir_path    = "output/"+self.year+"/"+self.channel+"_"+channel_postfix_+"/"
            self.in_hist_path_txt = "inFiles/"+self.year+"/"+self.channel+"_"+channel_postfix_+"/fhist.txt"
    
        # Make output directory
        if not os.path.exists(self.out_dir_path):
            os.makedirs(self.out_dir_path)
    
        # Read text file including root file paths for unfolding
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
                                    int(self.year), 
                                    self.mode, 
                                    False, 
                                    self.var, 
                                    self.out_dir_path,   
                                    self.density)

        self.unfold.setBias(self.bias)
        
        # Set response matrix
        self.unfold.setNominalRM(self.input_root_file["nominal"][self.matrix_filekey], self.dir_path, self.bin_def)

    def setPartialReg2D(self, regMode = 0, startMass = 320., startPt = 0., endMass = 320., endPt = 99.) :
        pass
       
    def setInputHist(self, useMCInput = False, 
                    unfoldObj = None, 
                    sys_type = "Type_0", 
                    sys_name = "", 
                    isFSR = False, 
                    useMadgraph = False, 
                    inputBinDef = None, 
                    inputHistName="", 
                    useAccept=True) :
        
        hist_filekey_temp = "Data"
        hist_postfix = self.makeSysHistPostfix(HistType.unf_input_hist, sys_type, sys_name)

        if inputBinDef is not None :
            bin_def = inputBinDef
        else :
            bin_def = self.bin_def

        # TODO add selection key for nomianl or SYS
        self.unfold.setUnfInput(self.input_root_file["nominal"][hist_filekey_temp], 
                                self.dir_path, 
                                bin_def, sys_type, sys_name, hist_postfix)

    def setFromPreviousUnfold(self, unfoldObj) :

        self.unfold.setFromPrevUnfResult(unfoldObj.getISRUnfold(self.var), True)
            
    def subFake(self, dirName = "Detector_DY_Fake", sys_type = "Type_0", sys_name = "", isFSR = False, inputBinDef = None):
            
        fakeList = ["DYJets", self.dy10to50HistName]
        hist_filekey_temp = "matrix" # FIXME save DY fake in the histogram file

        if isFSR :
            hist_filekey_temp = "fsr_matrix"

        if inputBinDef is not None :
            bin_def = inputBinDef
           
        else :
            bin_def = self.bin_def
         

        histPostfix_temp = self.makeSysHistPostfix(HistType.unf_bkg_hist, sys_type, sys_name)  
        #print(histPostfix_temp)

        for fake_name in fakeList:

            self.unfold.subBkgs(self.input_root_file[hist_filekey_temp], dirName, bin_def, fake_name, sys_type, sys_name, histPostfix_temp)
        
    def setUnfoldBkgs(self, dirName = "ee2016", sys_type = "Type_0", sys_name = ""):
   
        if self.channel == "ee" : 
            bkgList = ["DY_tau", "ttbar", "top", "antitop", "ww", "wz", "zz"] 
        else :
            bkgList = ["DYJets10to50ToTauTau", "DYJetsToTauTau", "WW_pythia", "WZ_pythia", "ZZ_pythia", "TTLL_powheg", "SingleTop_tW_antitop_NoFullyHad", "SingleTop_tW_top_NoFullyHad"] 

        histPostfix_temp = self.makeSysHistPostfix(HistType.unf_bkg_hist, sys_type, sys_name) 
        #print(histPostfix_temp)
        
        for bkg_name in bkgList:
            hist_filekey_temp = "hist"
            if "SingleTop" in bkg_name and "AlphaS" in sys_name : # FIXME
                histPostfix_temp = self.makeSysHistPostfix(HistType.unf_bkg_hist, "Type_3", sys_name) 

            self.unfold.subBkgs(self.input_root_file["nominal"][bkg_name], dirName, self.bin_def, bkg_name, sys_type, sys_name, histPostfix_temp)
            
    def setSystematics(self, sys_type, sys_name, isFSR = False):

        #print("setSystematics")
        self.unfold.setSystematics(sys_name)

        matrix_filekey_temp = self.matrix_filekey
        matrix_dir_path_temp = self.matrix_dirPath
       
        hist_postfix_temp = self.makeSysHistPostfix(HistType.unf_res_matrix, sys_type, sys_name)

        if sys_name == "UnfoldingModel" :
            matrix_filekey_temp = "matrix_mg"
        if "fsrPHOTOS" in sys_name :
            matrix_filekey_temp = "fsr_matrix_pythia"
        if "fsrPYTHIA" in sys_name :
            matrix_filekey_temp = "fsr_matrix_photos"

        # make a function to make a full histogram name with systematic
        #print("setSystematics", sys_name, hist_postfix_temp)
        self.unfold.setSystematicRM(self.input_root_file["SYS"][matrix_filekey_temp], matrix_dir_path_temp, self.bin_def, sys_name, hist_postfix_temp)

    def makeSysHistPostfix(self, histType, sys_type, sys_name) :

        if histType == HistType.unf_input_hist : # Input histogram

            if sys_type == "Type_1":
                sysHistPostfix = "_" + sys_name
            else :
                sysHistPostfix = ""

            return sysHistPostfix 

        elif histType == HistType.unf_bkg_hist : # Background histogram

            if sys_type == "Type_1" or sys_type == "Type_2":
                sysHistPostfix = "_" + sys_name
            else :
                sysHistPostfix = ""
            return sysHistPostfix

        elif histType == HistType.unf_res_matrix : # Response matrix
            
            if sys_type == "Type_1" or sys_type == "Type_2" or sys_type == "Type_4" :
                sysHistPostfix = "_" + sys_name
                if sys_name == "UnfoldingModel" or "fsr" in sys_name:
                    sysHistPostfix = ""
            else :
                sysHistPostfix = "" 

            return sysHistPostfix

        else :
            print("Which histogram do you mean?")

    def checkIterEMUnfold(self):

        self.unfold.checkIterEMUnfold()

    # Do unfold! 
    def doUnfold(self, partialReg=True):

        self.unfold.doISRUnfold(partialReg)


    def doStatUnfold(self):

        self.unfold.doStatUnfold()
    

    def getISRUnfold(self):
        return self.unfold
    
    def doAcceptance(self, isFSR = False, outName = "", useMassBinned = False) :

        if isFSR :
        
            if useMassBinned :  
                self.unfold.doAcceptCorr(self.input_root_file["nominal"]['hist_accept_fullPhase'], 
                                         self.bin_def, self.input_root_file["nominal"]['hist_updated_accept'])
                
            else :
                self.unfold.doAcceptCorr(self.input_root_file["nominal"]['hist_accept_fullPhase'], self.bin_def)
                
         
        else : 
            if useMassBinned :
                self.unfold.doAcceptCorr(self.input_root_file["nominal"]['hist_accept_drp1'], 
                                         self.bin_def, self.input_root_file["nominal"]['hist_updated_accept'])
                
            else :
                self.unfold.doAcceptCorr(self.input_root_file["nominal"]['hist_accept_drp1'], self.bin_def)
      

    def drawCorrelation(self, var = "Mass", steering = None, useAxis = True, outName = ""):

        self.unfold.drawCorrelation(var, steering, useAxis, outName)

    def combineOutFiles(self, prefix = "") :

        pt_output_path = self.out_dir_path + self.unfold_name + "_" + self.channel + "_" + self.year + "_Pt.root" 
        mass_output_path = self.out_dir_path + self.unfold_name + "_" + self.channel + "_" + self.year + "_Mass.root" 
        target_file_path = self.out_dir_path + self.unfold_name + "_" + self.channel + "_" + self.year + ".root"   
        if prefix != "" :
            target_file_path = self.out_dir_path + self.unfold_name + "_" + self.channel + "_" + self.year + "_" + prefix + ".root"
        os.system('hadd -f ' + target_file_path + " " + pt_output_path + " " + mass_output_path)

    def closeOutFiles(self) :

        self.unfold.closeOutFile()
