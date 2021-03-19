import os
import sys
from array import array
import ROOT as rt
import math

import pyScripts.unfoldUtil as unfoldutil
import pandas as pd

class ISRAnalysis:
    
    def __init__(self, unfold_name_ = "Detector", year_ = "2016", channel_= "electron", regMode_ = 0, doInputStat_ = False, doRMStat_ = False, matrix_filekey_ = "matrix",
                 matrix_dirPath_ = "Detector_Dressed_DRp1_Fiducial",  binDef_ = ("",""), channel_postfix_ = ""):
        
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
        self.bias = 1.0
        self.mode = regMode_ # 레귤라이제이션 모드
        
        # Create ISRUnfold object
        # unfold_name : prefix for output plots
        # Make two ISRUnfold object for mass and pt
        print("Creat ISRUnfold objects...")
        self.unfold_pt   = rt.ISRUnfold(self.unfold_name, self.channel, int(self.year), int(self.mode), doInputStat_, doRMStat_, False, "Pt", self.outDirPath)
        self.unfold_mass = rt.ISRUnfold(self.unfold_name, self.channel, int(self.year), int(self.mode), doInputStat_, doRMStat_, False, "Mass", self.outDirPath)

        self.unfold_pt.setBias(self.bias)
        self.unfold_mass.setBias(self.bias)
        
        # Set response matrix
        print(self.binDef)
        self.unfold_pt.setNominalRM(self.inHistDic[self.matrix_filekey], self.matrix_dirPath, self.binDef[0])
        self.unfold_mass.setNominalRM(self.inHistDic[self.matrix_filekey], self.matrix_dirPath, self.binDef[1])
       
    def checkMatrixCond(self):
        self.unfold_mass.checkMatrixCond()
        self.unfold_pt.checkMatrixCond()

    def setInputHist(self, useMCInput = False, unfoldObj = None, dirName = "Detector", sys_type = "Type_0", sys_name = "", isFSR = False, useMadgraph = False, inputBinDef = None):
        
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

            else :
                input_hist_name = "histo_DYJets"
                hist_filekey_temp = "hist_accept_drp1"

            self.unfold_pt.setUnfInput(self.inHistDic[hist_filekey_temp], dirName, pt_mass_bin_def, input_hist_name, sys_type, sys_name, hist_postfix, isFSR)
            self.unfold_mass.setUnfInput(self.inHistDic[hist_filekey_temp], dirName, mass_bin_def, input_hist_name, sys_type, sys_name, hist_postfix, isFSR)

        else :

            if unfoldObj is None:

                self.unfold_pt.setUnfInput(self.inHistDic[hist_filekey_temp], dirName, pt_mass_bin_def, input_hist_name, sys_type, sys_name, hist_postfix)
                self.unfold_mass.setUnfInput(self.inHistDic[hist_filekey_temp], dirName, mass_bin_def, input_hist_name, sys_type, sys_name, hist_postfix)

            else:

                # Set unfold input from previous unfold 
                self.unfold_pt.setUnfInput(unfoldObj.getISRUnfold("Pt"),   sys_type, sys_name, True) 
                self.unfold_mass.setUnfInput(unfoldObj.getISRUnfold("Mass"),sys_type, sys_name, True)

    def setFromPreviousUnfold(self, unfoldObj) :

        self.unfold_pt.setFromPrevUnfResult(unfoldObj.getISRUnfold("Pt"), True)
        self.unfold_mass.setFromPrevUnfResult(unfoldObj.getISRUnfold("Mass"), True)
            
    def subFake(self, dirName = "Detector_DY_Fake", sys_type = "Type_0", sys_name = "", isFSR = False, inputBinDef = None):
            
        fakeList = ["DYJets", self.dy10to50HistName]
        hist_filekey_temp = "matrix" # FIXME save DY fake in the histogram file

        if isFSR :
            hist_filekey_temp = "hist_accept_drp1"

        if inputBinDef is not None :
            pt_mass_bin_def = inputBinDef[0]
            mass_bin_def = inputBinDef[1]
        else :
            pt_mass_bin_def = self.binDef[0]
            mass_bin_def = self.binDef[1]

        histPostfix_temp = self.makeSysHistPostfix(1, sys_type, sys_name)  
        print(histPostfix_temp)

        for fake_name in fakeList:

            self.unfold_pt.subBkgs(self.inHistDic[hist_filekey_temp], dirName, pt_mass_bin_def, fake_name, sys_type, sys_name, histPostfix_temp)
            self.unfold_mass.subBkgs(self.inHistDic[hist_filekey_temp], dirName, mass_bin_def, fake_name, sys_type, sys_name, histPostfix_temp)
        
    def setUnfoldBkgs(self, dirName = "Detector", sys_type = "Type_0", sys_name = ""):
   
        bkgList = ["DYJets10to50ToTauTau", "DYJetsToTauTau", "WW_pythia", "WZ_pythia", "ZZ_pythia", "TTLL_powheg", "GamGamToLL"] 

        histPostfix_temp = self.makeSysHistPostfix(1, sys_type, sys_name) 
        print(histPostfix_temp)
        
        for bkg_name in bkgList:
            hist_filekey_temp = "hist"

            self.unfold_pt.subBkgs(self.inHistDic[hist_filekey_temp], dirName, self.binDef[0], bkg_name, sys_type, sys_name, histPostfix_temp)
            self.unfold_mass.subBkgs(self.inHistDic[hist_filekey_temp], dirName, self.binDef[1], bkg_name, sys_type, sys_name, histPostfix_temp)
            
    def setSystematics(self, sys_type, sys_name, isFSR = False):

        print("setSystematics")
        self.unfold_pt.setSystematics(sys_name)
        self.unfold_mass.setSystematics(sys_name)

        matrix_filekey_temp = self.matrix_filekey
        matrix_dir_path_temp = self.matrix_dirPath
       
        hist_postfix_temp = self.makeSysHistPostfix(2, sys_type, sys_name)

        if sys_name == "UnfoldingModel" :
            matrix_filekey_temp = "matrix_mg"

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
                if sys_name == "UnfoldingModel" :
                    sysHistPostfix = ""
            else :
                sysHistPostfix = "" 

            return sysHistPostfix

        else :
            print("Which histogram do you mean?")

    def checkIterEMUnfold(self):

        self.unfold.checkIterEMUnfold()

    # Do unfold! 
    def doUnfold(self):

        self.unfold_pt.doISRUnfold()
        self.unfold_mass.doISRUnfold()

    def doStatUnfold(self):

        self.unfold_pt.doStatUnfold()
        self.unfold_mass.doStatUnfold()

    def getISRUnfold(self, var_ = "Mass"):
        
        if var_ == "Mass" :
            return self.unfold_mass
        else :
            return self.unfold_pt
    
    def doAcceptance(self, isFSR = False, outName = "") :

        if isFSR :
            self.unfold_pt.doAcceptCorr(self.inHistDic['hist_accept_fullPhase'], self.binDef[0],  True)
            self.unfold_mass.doAcceptCorr(self.inHistDic['hist_accept_fullPhase'], self.binDef[1],  True)
        else : 
            self.unfold_pt.doAcceptCorr(self.inHistDic['hist_accept_drp1'], self.binDef[0])
            self.unfold_mass.doAcceptCorr(self.inHistDic['hist_accept_drp1'], self.binDef[1])

    def drawCorrelation(self, var = "Mass", steering = None, useAxis = True, outName = ""):

        self.unfold.drawCorrelation(var, steering, useAxis, outName)

    def combineOutFiles(self) :

        pt_output_path = self.outDirPath + self.unfold_name + "_" + self.channel + "_" + self.year + "_Pt.root" 
        mass_output_path = self.outDirPath + self.unfold_name + "_" + self.channel + "_" + self.year + "_Mass.root" 
        target_file_path = self.outDirPath + self.unfold_name + "_" + self.channel + "_" + self.year + ".root"   
        os.system('hadd -f ' + target_file_path + " " + pt_output_path + " " + mass_output_path)

    def closeOutFiles(self) :

        self.unfold_pt.closeOutFile()
        self.unfold_mass.closeOutFile()

    # Get histograms
    def getCDFPtVsMassTGraph(self, grTitle = ""):

        meanMass, meanPt = array('d'), array('d')
        meanMassStatErr, meanPtStatErr = array('d'), array('d')
        meanMassSysErr, meanPtSysErr = array('d'), array('d')

        # Muon
        cdf_muon_mass      = [47.72, 70.66, 90.99, 115.29, 243.33]
        cdf_muon_mass_stat = [0.05, 0.04, 0.01, 0.18, 1.63]
        cdf_muon_mass_sys  = [0.04, 0.07, 0.08, 0.14, 0.40]
        cdf_muon_mass_tot  = []

        cdf_electron_mass      = [47.83, 70.76, 90.98, 115.11, 245.46]
        cdf_electron_mass_stat = [0.05, 0.04, 0.01, 0.13, 1.29]
        cdf_electron_mass_sys  = [0.07, 0.04, 0.07, 0.14, 0.21]
        cdf_electron_mass_tot  = []
        
        # Electron
        cdf_muon_pt      = [9.12, 10.81, 11.84, 13.17, 16.18]
        cdf_muon_pt_stat = [0.09, 0.08, 0.03, 0.12, 0.61]
        cdf_muon_pt_sys  = [0.12, 0.14, 0.03, 0.12, 0.45]
        cdf_muon_pt_tot  = []

        cdf_electron_pt      = [9.10, 10.84, 11.79, 12.93, 16.41]
        cdf_electron_pt_stat = [0.13, 0.08, 0.02, 0.09, 0.56]
        cdf_electron_pt_sys  = [0.18, 0.10, 0.01, 0.09, 0.35]
        cdf_electron_pt_tot  = []

        for i in range(5):
            cdf_muon_mass_tot.append(math.sqrt(math.pow(cdf_muon_mass_sys[i], 2) + math.pow(cdf_muon_mass_stat[i], 2)))
            cdf_electron_mass_tot.append(math.sqrt(math.pow(cdf_electron_mass_sys[i], 2) + math.pow(cdf_electron_mass_stat[i], 2)))

            cdf_muon_pt_tot.append(math.sqrt(math.pow(cdf_muon_pt_stat[i], 2) + math.pow(cdf_muon_pt_sys[i], 2)))
            cdf_electron_pt_tot.append(math.sqrt(math.pow(cdf_electron_pt_stat[i], 2) + math.pow(cdf_electron_pt_sys[i], 2)))

        for ibin in range(5):
            if self.channel == "electron":
                    meanMass.append(cdf_electron_mass[ibin])
                    meanPt.append(cdf_electron_pt[ibin])
                    meanMassStatErr.append(cdf_electron_mass_stat[ibin])
                    meanPtStatErr.append(cdf_electron_pt_stat[ibin])

                    meanMassSysErr.append(cdf_electron_mass_tot[ibin])
                    meanPtSysErr.append(cdf_electron_pt_tot[ibin])
                                 
            if self.channel == "muon":
                    meanMass.append(cdf_muon_mass[ibin])
                    meanPt.append(cdf_muon_pt[ibin])
                    meanMassStatErr.append(cdf_muon_mass_stat[ibin])
                    meanPtStatErr.append(cdf_muon_pt_stat[ibin])

                    meanMassSysErr.append(cdf_muon_mass_tot[ibin])
                    meanPtSysErr.append(cdf_muon_pt_tot[ibin])

        gr = rt.TGraphErrors(self.nMassBins, meanMass, meanPt, meanMassSysErr, meanPtSysErr)
    
        gr.SetName(grTitle)
        return gr
