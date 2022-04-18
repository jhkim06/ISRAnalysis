import os
import sys
from array import array
import ROOT as rt
import math

import pyScripts.unfoldUtil as unfoldutil
import pandas as pd


class ISRAnalysis:
    
    def __init__(self, unfold_name_ = "Detector", year_ = "2016", channel_= "electron", regMode_ = 0, doInputStat_ = False, doRMStat_ = False, ignoreBinZero_ = False, matrix_filekey_ = "matrix",
                 matrix_dirPath_ = "Detector_Dressed_DRp1_Fiducial",  binDef_ = ("",""), channel_postfix_ = "", doModelUnc_ = False):
        
        # Initialize some variables 
        self.unfold_name = unfold_name_
        self.matrix_filekey = matrix_filekey_
        self.matrix_dirPath  = matrix_dirPath_
        self.binDef = binDef_
        
        self.channel = channel_
        self.year = year_
        self.nMassBins = None
       
        # Set data and Drell-Yan input histogram names
        dataHistPrefix = "Double"
        if self.channel == "electron":
            dataHistPrefix = dataHistPrefix + "EG"
            if self.year == "2018":
                dataHistPrefix = "EGamma"
        if self.channel == "muon":
            dataHistPrefix = dataHistPrefix + "Muon"
        self.data_hist_name = "histo_"+dataHistPrefix
            
        self.dy10to50HistName = "DYJets10to50"
        MG_postfix = "_MG"
        if self.year != "2016":
            self.dy10to50HistName += MG_postfix
        
        self.outDirPath = "output/"+self.year+"/"+self.channel+"/"
        self.inHistPathTxt = "inFiles/"+self.year+"/"+self.channel+"/fhist.txt"

        if channel_postfix_ != "" :
            self.outDirPath = "output/"+self.year+"/"+self.channel+"_"+channel_postfix_+"/"
            self.inHistPathTxt = "inFiles/"+self.year+"/"+self.channel+"_"+channel_postfix_+"/fhist.txt"
    
        # Make output directory
        if not os.path.exists(self.outDirPath):
            os.makedirs(self.outDirPath)
    
        # Read text file including root file paths for unfolding
        self.filePaths = open(self.inHistPathTxt, 'r')
        self.inHistDic = {}
        
        for path in self.filePaths:
            modifiedPath = path.lstrip(' ').rstrip(' ').rstrip('\n')
            self.inHistDic[modifiedPath.split()[1]] = modifiedPath.split()[2]
        
        # Unfolding configuration
        self.bias = 1.
        self.mode = regMode_ 
        
        # Create ISRUnfold object
        # unfold_name : prefix for output plots
        # Make two ISRUnfold object for mass and pt
        print("Create ISRUnfold objects...")

        if doModelUnc_ :
            model_sys_file_path = self.inHistDic["matrix_reweightSF"]
            if "FSR" in self.unfold_name :
                model_sys_file_path = self.inHistDic["fsr_matrix_reweightSF"]
        else :
            model_sys_file_path = ""
            
        # Create ISRUnfold object
        self.unfold_pt   = rt.ISRUnfold(self.unfold_name, self.channel, int(self.year), int(self.mode), doInputStat_, doRMStat_, ignoreBinZero_, False, "Pt", self.outDirPath,   model_sys_file_path, doModelUnc_)
        self.unfold_mass = rt.ISRUnfold(self.unfold_name, self.channel, int(self.year), int(self.mode), doInputStat_, doRMStat_, ignoreBinZero_, False, "Mass", self.outDirPath, model_sys_file_path, doModelUnc_)

        self.unfold_pt.setBias(self.bias)
        self.unfold_mass.setBias(self.bias)
        
        # Set response matrix
        #print(self.binDef)
        self.unfold_pt.setNominalRM(self.inHistDic[self.matrix_filekey], self.matrix_dirPath, self.binDef[0])
        self.unfold_mass.setNominalRM(self.inHistDic[self.matrix_filekey], self.matrix_dirPath, self.binDef[1])
       
    def checkMatrixCond(self):
        self.unfold_mass.checkMatrixCond()
        self.unfold_pt.checkMatrixCond()

    def setInputHist(self, useMCInput = False, unfoldObj = None, dirName = "Detector", sys_type = "Type_0", sys_name = "", isFSR = False, useMadgraph = False, inputBinDef = None, inputHistName="", useAccept=True):
        
        input_hist_name = self.data_hist_name
        hist_filekey_temp = "hist"
        hist_postfix = self.makeSysHistPostfix(0, sys_type, sys_name)

        if inputBinDef is not None :
            pt_mass_bin_def = inputBinDef[0]
            mass_bin_def = inputBinDef[1]
        else :
            pt_mass_bin_def = self.binDef[0]
            mass_bin_def = self.binDef[1]

        # For closure test, use MC as unfolding input
        if useMCInput == True:

            if isFSR == False:
                if self.channel == "electron":
                    input_hist_name = "histo_DYJetsToEE"
                else :
                    input_hist_name = "histo_DYJetsToMuMu"

                if useMadgraph == True:
                    hist_filekey_temp = "hist_DYMG"

                if inputHistName != "" :
                    input_hist_name = inputHistName

            else :
                input_hist_name = "histo_DYJets"
                hist_filekey_temp = "hist_accept_drp1"

                if inputHistName != "" :
                    input_hist_name = inputHistName

            self.unfold_pt.setUnfInput(self.inHistDic[hist_filekey_temp], dirName, pt_mass_bin_def, input_hist_name, sys_type, sys_name, hist_postfix, isFSR)
            self.unfold_mass.setUnfInput(self.inHistDic[hist_filekey_temp], dirName, mass_bin_def, input_hist_name, sys_type, sys_name, hist_postfix, isFSR)

        else :

            if unfoldObj is None:

                self.unfold_pt.setUnfInput(self.inHistDic[hist_filekey_temp], dirName, pt_mass_bin_def, input_hist_name, sys_type, sys_name, hist_postfix)
                self.unfold_mass.setUnfInput(self.inHistDic[hist_filekey_temp], dirName, mass_bin_def, input_hist_name, sys_type, sys_name, hist_postfix)

            else:

                # Set unfold input from previous unfold 
                self.unfold_pt.setUnfInput(unfoldObj.getISRUnfold("Pt"),   sys_type, sys_name, useAccept) 
                self.unfold_mass.setUnfInput(unfoldObj.getISRUnfold("Mass"),sys_type, sys_name,useAccept)

    def setFromPreviousUnfold(self, unfoldObj) :

        self.unfold_pt.setFromPrevUnfResult(unfoldObj.getISRUnfold("Pt"), True)
        self.unfold_mass.setFromPrevUnfResult(unfoldObj.getISRUnfold("Mass"), True)
            
    def subFake(self, dirName = "Detector_DY_Fake", sys_type = "Type_0", sys_name = "", isFSR = False, inputBinDef = None):
            
        fakeList = ["DYJets", self.dy10to50HistName]
        hist_filekey_temp = "matrix" # FIXME save DY fake in the histogram file

        if isFSR :
            hist_filekey_temp = "fsr_matrix"

        if inputBinDef is not None :
            pt_mass_bin_def = inputBinDef[0]
            mass_bin_def = inputBinDef[1]
        else :
            pt_mass_bin_def = self.binDef[0]
            mass_bin_def = self.binDef[1]

        histPostfix_temp = self.makeSysHistPostfix(1, sys_type, sys_name)  
        #print(histPostfix_temp)

        for fake_name in fakeList:

            self.unfold_pt.subBkgs(self.inHistDic[hist_filekey_temp], dirName, pt_mass_bin_def, fake_name, sys_type, sys_name, histPostfix_temp)
            self.unfold_mass.subBkgs(self.inHistDic[hist_filekey_temp], dirName, mass_bin_def, fake_name, sys_type, sys_name, histPostfix_temp)
        
    def setUnfoldBkgs(self, dirName = "Detector", sys_type = "Type_0", sys_name = ""):
   
        if self.channel == "electron" : 
            bkgList = ["DYJets10to50ToTauTau", "DYJetsToTauTau", "WW_pythia", "WZ_pythia", "ZZ_pythia", "TTLL_powheg", "SingleTop_tW_antitop_NoFullyHad", "SingleTop_tW_top_NoFullyHad"] 
        else :
            bkgList = ["DYJets10to50ToTauTau", "DYJetsToTauTau", "WW_pythia", "WZ_pythia", "ZZ_pythia", "TTLL_powheg", "SingleTop_tW_antitop_NoFullyHad", "SingleTop_tW_top_NoFullyHad"] 

        histPostfix_temp = self.makeSysHistPostfix(1, sys_type, sys_name) 
        #print(histPostfix_temp)
        
        for bkg_name in bkgList:
            hist_filekey_temp = "hist"
            if "SingleTop" in bkg_name and "AlphaS" in sys_name : # FIXME
                histPostfix_temp = self.makeSysHistPostfix(1, "Type_3", sys_name) 

            self.unfold_pt.subBkgs(self.inHistDic[hist_filekey_temp], dirName, self.binDef[0], bkg_name, sys_type, sys_name, histPostfix_temp)
            self.unfold_mass.subBkgs(self.inHistDic[hist_filekey_temp], dirName, self.binDef[1], bkg_name, sys_type, sys_name, histPostfix_temp)
            
    def setSystematics(self, sys_type, sys_name, isFSR = False):

        #print("setSystematics")
        self.unfold_pt.setSystematics(sys_name)
        self.unfold_mass.setSystematics(sys_name)

        matrix_filekey_temp = self.matrix_filekey
        matrix_dir_path_temp = self.matrix_dirPath
       
        hist_postfix_temp = self.makeSysHistPostfix(2, sys_type, sys_name)

        if sys_name == "UnfoldingModel" :
            matrix_filekey_temp = "matrix_mg"
        if "fsrPHOTOS" in sys_name :
            matrix_filekey_temp = "fsr_matrix_pythia"
        if "fsrPYTHIA" in sys_name :
            matrix_filekey_temp = "fsr_matrix_photos"

        # make a function to make a full histogram name with systematic
        #print("setSystematics", sys_name, hist_postfix_temp)
        self.unfold_pt.setSystematicRM(self.inHistDic[matrix_filekey_temp], matrix_dir_path_temp, self.binDef[0], sys_name, hist_postfix_temp)
        self.unfold_mass.setSystematicRM(self.inHistDic[matrix_filekey_temp], matrix_dir_path_temp, self.binDef[1], sys_name, hist_postfix_temp)

    def makeSysHistPostfix(self, histType, sys_type, sys_name) :

        if histType == 0 : # Input histogram

            if sys_type == "Type_1":
                sysHistPostfix = "_" + sys_name
            else :
                sysHistPostfix = ""

            return sysHistPostfix 

        elif histType == 1 : # Background histogram

            if sys_type == "Type_1" or sys_type == "Type_2":
                sysHistPostfix = "_" + sys_name
            else :
                sysHistPostfix = ""
            return sysHistPostfix

        elif histType == 2 : # Response matrix
            
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

        self.unfold_pt.doISRUnfold(partialReg)
        self.unfold_mass.doISRUnfold(partialReg)

    def doStatUnfold(self):

        self.unfold_pt.doStatUnfold()
        self.unfold_mass.doStatUnfold()

    def getISRUnfold(self, var_ = "Mass"):
        
        if var_ == "Mass" :
            return self.unfold_mass
        else :
            return self.unfold_pt
    
    def doAcceptance(self, isFSR = False, outName = "", useMassBinned = False) :

        if isFSR :
        
            if useMassBinned :  
                self.unfold_pt.doAcceptCorr(self.inHistDic['hist_accept_fullPhase'], self.binDef[0], self.inHistDic['hist_updated_accept'])
                pass
            else :
                self.unfold_pt.doAcceptCorr(self.inHistDic['hist_accept_fullPhase'], self.binDef[0])
                pass
            self.unfold_mass.doAcceptCorr(self.inHistDic['hist_accept_fullPhase'], self.binDef[1])
        else : 
            if useMassBinned :
                self.unfold_pt.doAcceptCorr(self.inHistDic['hist_accept_drp1'], self.binDef[0], self.inHistDic['hist_updated_accept'])
                pass
            else :
                self.unfold_pt.doAcceptCorr(self.inHistDic['hist_accept_drp1'], self.binDef[0])
                pass
                
            self.unfold_mass.doAcceptCorr(self.inHistDic['hist_accept_drp1'], self.binDef[1])

    def drawCorrelation(self, var = "Mass", steering = None, useAxis = True, outName = ""):

        self.unfold.drawCorrelation(var, steering, useAxis, outName)

    def combineOutFiles(self, prefix = "") :

        pt_output_path = self.outDirPath + self.unfold_name + "_" + self.channel + "_" + self.year + "_Pt.root" 
        mass_output_path = self.outDirPath + self.unfold_name + "_" + self.channel + "_" + self.year + "_Mass.root" 
        target_file_path = self.outDirPath + self.unfold_name + "_" + self.channel + "_" + self.year + ".root"   
        if prefix != "" :
            target_file_path = self.outDirPath + self.unfold_name + "_" + self.channel + "_" + self.year + "_" + prefix + ".root"
        os.system('hadd -f ' + target_file_path + " " + pt_output_path + " " + mass_output_path)

    def closeOutFiles(self) :

        self.unfold_pt.closeOutFile()
        self.unfold_mass.closeOutFile()
