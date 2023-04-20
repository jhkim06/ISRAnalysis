import os
import sys
from array import array
import ROOT as rt
import math

import pyScripts.unfoldUtil as unfoldutil
import pandas as pd

class HistType :

    unfInputHist = 0
    unfBkgHist = 1
    unfResMatrix = 2

class ISRAnalysis:
    
    def __init__(self, unfold_name_ = "Detector", year_ = "2016", channel_= "ee", regMode_ = 0, doInputStat_ = False, doRMStat_ = False, ignoreBinZero_ = False, matrix_filekey_ = "DY",
                 binDef_ = "", channel_postfix_ = "", doModelUnc_ = False, var_ = "Pt", bias_ = 0., density_ = 0):
        
        # Initialize some variables 
        self.unfold_name = unfold_name_
        self.matrix_filekey = matrix_filekey_
        self.binDef = binDef_
        self.hist_path = channel_ + year_  # TODO use postfix
        
        self.channel = channel_
        self.year = year_
        self.var = var_
        self.nMassBins = None
        
        self.outDirPath    = "output/" +self.year+"/"+self.channel+"/"
        self.inHistPathTxt = "inFiles/" +self.year+"/"+self.channel+"/fhist.txt"

        if channel_postfix_ != "" :
            self.outDirPath    = "output/"+self.year+"/"+self.channel+"_"+channel_postfix_+"/"
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
        self.bias = bias_
        self.mode = regMode_ 
        self.density = density_
        
        # Create ISRUnfold object
        # unfold_name : prefix for output plots
        print("Create ISRUnfold objects...")

        if doModelUnc_ :
            model_sys_file_path = self.inHistDic["matrix_reweightSF"]
            if "FSR" in self.unfold_name :
                model_sys_file_path = self.inHistDic["fsr_matrix_reweightSF"]
        else :
            model_sys_file_path = ""
            
        # Create ISRUnfold object
        self.unfold  = rt.ISRUnfold(self.unfold_name, self.channel, int(self.year), self.mode, doInputStat_, doRMStat_, ignoreBinZero_, False, self.var, self.outDirPath,   model_sys_file_path, doModelUnc_, self.density)
        self.unfold.setBias(self.bias)
        
        # Set response matrix
        self.unfold.setNominalRM(self.inHistDic[self.matrix_filekey], self.hist_path, self.binDef)

    def setPartialReg2D(self, regMode = 0, startMass = 320., startPt = 0., endMass = 320., endPt = 99.) :
        pass
       
    def checkMatrixCond(self):
        self.unfold.checkMatrixCond()

    def setInputHistUnfSys(self) :

        self.unfold.setUnfInputUnfSys()

    def setInputHist(self, useMCInput = False, unfoldObj = None, sys_type = "Type_0", sys_name = "", isFSR = False, useMadgraph = False, inputBinDef = None, inputHistName="", useAccept=True):
        
        hist_filekey_temp = "Data"
        hist_postfix = self.makeSysHistPostfix(HistType.unfInputHist, sys_type, sys_name)

        if inputBinDef is not None :
            bin_def = inputBinDef
    
        else :
            bin_def = self.binDef

        if unfoldObj is None:
            self.unfold.setUnfInput(self.inHistDic[hist_filekey_temp], self.hist_path, bin_def, sys_type, sys_name, hist_postfix)

        else:
            # Set unfold input from previous unfold 
            self.unfold.setUnfInput(unfoldObj.getISRUnfold(self.var), sys_type, sys_name, useAccept)

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
            bin_def = self.binDef
         

        histPostfix_temp = self.makeSysHistPostfix(HistType.unfBkgHist, sys_type, sys_name)  
        #print(histPostfix_temp)

        for fake_name in fakeList:

            self.unfold.subBkgs(self.inHistDic[hist_filekey_temp], dirName, bin_def, fake_name, sys_type, sys_name, histPostfix_temp)
       
        
    def setUnfoldBkgs(self, dirName = "ee2016", sys_type = "Type_0", sys_name = ""):
   
        if self.channel == "ee" : 
            bkgList = ["DY_tau", "ttbar", "ww", "wz", "zz"] 
        else :
            bkgList = ["DYJets10to50ToTauTau", "DYJetsToTauTau", "WW_pythia", "WZ_pythia", "ZZ_pythia", "TTLL_powheg", "SingleTop_tW_antitop_NoFullyHad", "SingleTop_tW_top_NoFullyHad"] 

        histPostfix_temp = self.makeSysHistPostfix(HistType.unfBkgHist, sys_type, sys_name) 
        #print(histPostfix_temp)
        
        for bkg_name in bkgList:
            hist_filekey_temp = "hist"
            if "SingleTop" in bkg_name and "AlphaS" in sys_name : # FIXME
                histPostfix_temp = self.makeSysHistPostfix(HistType.unfBkgHist, "Type_3", sys_name) 

            self.unfold.subBkgs(self.inHistDic[bkg_name], dirName, self.binDef, bkg_name, sys_type, sys_name, histPostfix_temp)
            
    def setSystematics(self, sys_type, sys_name, isFSR = False):

        #print("setSystematics")
        self.unfold.setSystematics(sys_name)


        matrix_filekey_temp = self.matrix_filekey
        matrix_dir_path_temp = self.matrix_dirPath
       
        hist_postfix_temp = self.makeSysHistPostfix(HistType.unfResMatrix, sys_type, sys_name)

        if sys_name == "UnfoldingModel" :
            matrix_filekey_temp = "matrix_mg"
        if "fsrPHOTOS" in sys_name :
            matrix_filekey_temp = "fsr_matrix_pythia"
        if "fsrPYTHIA" in sys_name :
            matrix_filekey_temp = "fsr_matrix_photos"

        # make a function to make a full histogram name with systematic
        #print("setSystematics", sys_name, hist_postfix_temp)
        self.unfold.setSystematicRM(self.inHistDic[matrix_filekey_temp], matrix_dir_path_temp, self.binDef, sys_name, hist_postfix_temp)

    def makeSysHistPostfix(self, histType, sys_type, sys_name) :

        if histType == HistType.unfInputHist : # Input histogram

            if sys_type == "Type_1":
                sysHistPostfix = "_" + sys_name
            else :
                sysHistPostfix = ""

            return sysHistPostfix 

        elif histType == HistType.unfBkgHist : # Background histogram

            if sys_type == "Type_1" or sys_type == "Type_2":
                sysHistPostfix = "_" + sys_name
            else :
                sysHistPostfix = ""
            return sysHistPostfix

        elif histType == HistType.unfResMatrix : # Response matrix
            
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
                self.unfold.doAcceptCorr(self.inHistDic['hist_accept_fullPhase'], self.binDef, self.inHistDic['hist_updated_accept'])
                
            else :
                self.unfold.doAcceptCorr(self.inHistDic['hist_accept_fullPhase'], self.binDef)
                
         
        else : 
            if useMassBinned :
                self.unfold.doAcceptCorr(self.inHistDic['hist_accept_drp1'], self.binDef, self.inHistDic['hist_updated_accept'])
                
            else :
                self.unfold.doAcceptCorr(self.inHistDic['hist_accept_drp1'], self.binDef)
      

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

        self.unfold.closeOutFile()
