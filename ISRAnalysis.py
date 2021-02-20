import os
import sys
from array import array
import ROOT as rt
import math

import pyScripts.unfoldUtil as unfoldutil
import pandas as pd

class ISRAnalysis:
    
    def __init__(self, unfold_name_ = "Detector", year_ = "2016", channel_= "electron", regMode = 0, sys_ = False, matrix_filekey_ = "matrix",
                 matrix_dirPath_ = "Detector_Dressed_DRp1_Fiducial", matrix_histName_ = "Detector_Dressed_DRp1", binDef_ = "", channel_postfix_ = ""):
        
        # Initialize some variables 
        self.unfold_name = unfold_name_
        self.matrix_filekey = matrix_filekey_
        self.matrix_dirPath  = matrix_dirPath_
        self.matrix_histName = matrix_histName_
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
        self.dataHistName = "histo_"+dataHistPrefix
            
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
        
        # 언폴딩 컨피규레이션
        self.bias = 1.0
        self.mode = regMode # 레귤라이제이션 모드
        
        # Create ISRUnfold object
        # unfold_name : prefix for output plots
        # Make two ISRUnfold object for mass and pt
        self.unfold_pt   = rt.ISRUnfold(self.unfold_name, self.channel, int(self.year), int(self.mode), sys_, False, "Pt", self.outDirPath)
        self.unfold_mass = rt.ISRUnfold(self.unfold_name, self.channel, int(self.year), int(self.mode), sys_, False, "Mass", self.outDirPath)

        self.unfold_pt.setBias(self.bias)
        self.unfold_mass.setBias(self.bias)
        
        # Set response matrix
        self.unfold_pt.setNominalRM(self.inHistDic[self.matrix_filekey], self.matrix_dirPath, self.matrix_histName, self.binDef)
        self.unfold_mass.setNominalRM(self.inHistDic[self.matrix_filekey], self.matrix_dirPath, self.matrix_histName, self.binDef)
       
    def checkMatrixCond(self):
        self.unfold_mass.checkMatrixCond()
        self.unfold_pt.checkMatrixCond()

    def setInputHist(self, useMCInput = False, unfoldObj = None, dirName = "Detector", isSys = False, sysName = "nominal", sysPostfix = "", isFSR = False, useMadgraph = False):
        
        inputHistName = self.dataHistName
        hist_filekey_temp = "hist"

        # For closure test, use MC as unfolding input
        if useMCInput == True:

            if isFSR == False:
                if self.channel == "electron":
                    inputHistName = "histo_DYJetsToEE"
                else :
                    inputHistName = "histo_DYJetsToMuMu"

                if useMadgraph == True:
                    hist_filekey_temp = "hist_DYMG"
 
                self.unfold_pt.setUnfInput(self.binDef, self.inHistDic[hist_filekey_temp], dirName, inputHistName, isSys, sysName, sysPostfix)
                self.unfold_mass.setUnfInput(self.binDef, self.inHistDic[hist_filekey_temp], dirName, inputHistName, isSys, sysName, sysPostfix)
            else :
                inputHistName = "histo_DYJets"

                if useMadgraph == False:
                    hist_filekey_temp = "hist_accept_drp1"

                else:
                    hist_filekey_temp = "hist_accept_drp1_DYMG"

                self.unfold_pt.setUnfInput("Gen_FineCoarse", self.inHistDic[hist_filekey_temp], dirName, inputHistName, isSys, sysName, sysPostfix, isFSR)
                self.unfold_mass.setUnfInput("Gen_FineCoarse", self.inHistDic[hist_filekey_temp], dirName, inputHistName, isSys, sysName, sysPostfix, isFSR)

        else :
            if unfoldObj is None:

                self.unfold_pt.setUnfInput(self.binDef, self.inHistDic[hist_filekey_temp], dirName, inputHistName, isSys, sysName, sysPostfix)
                self.unfold_mass.setUnfInput(self.binDef, self.inHistDic[hist_filekey_temp], dirName, inputHistName, isSys, sysName, sysPostfix)
            else:
                # Set unfold input from previous unfold 
                self.unfold_pt.setUnfInput(unfoldObj.getISRUnfold("Pt"), isSys, sysName, sysPostfix, True)
                self.unfold_mass.setUnfInput(unfoldObj.getISRUnfold("Mass"), isSys, sysName, sysPostfix, True)

    def setFromPreviousUnfold(self, unfoldObj) :

        self.unfold_pt.setFromPrevUnfResult(unfoldObj.getISRUnfold("Pt"), True)
        self.unfold_mass.setFromPrevUnfResult(unfoldObj.getISRUnfold("Mass"), True)
            
    def subFake(self, isSys = False, dirName = "Detector_DY_Fake", sysName = "nominal", sysPostfix = "", histPostfix = "", isFSR = False):
            
        fakeList = {"DYJets": "DY", self.dy10to50HistName:"DY"}
       
        hist_filekey_temp = "matrix"
        if isFSR :
            hist_filekey_temp = "fsr_matrix"

        for fake in fakeList.items():

            if "ZptCorr" in sysName : 
                if sysPostfix != "Nominal" : 
                    hist_filekey_temp = 'matrix_zptcorr' 

            elif "FSR" in sysName : 
                histPostfix = ""

            elif "Background" in sysName : 
                histPostfix = ""

            self.unfold_pt.subBkgs(self.inHistDic[hist_filekey_temp], fake, isSys, self.binDef, dirName, sysName, sysPostfix, histPostfix)
            self.unfold_mass.subBkgs(self.inHistDic[hist_filekey_temp], fake, isSys, self.binDef, dirName, sysName, sysPostfix, histPostfix)
        
    def setUnfoldBkgs(self, isSys = False , dirName = "Detector", sysName = "nominal", sysPostfix = ""):
   
        bkgList = {}
        # 2016 데이터만 single top 샘플을 갖고 있다 
        if self.year == "2016" or self.year == "2017" or self.year == "2018":
            #bkgList = { "QCD": "Fake", "WJet": "Fake",
            bkgList = {
                        "WJets_MG": "WJets",
                        "DYJets10to50ToTauTau":"DY#rightarrow#tau#tau", "DYJetsToTauTau":"DY#rightarrow#tau#tau", 
                        "WW_pythia": "VV", "WZ_pythia": "VV", "ZZ_pythia": "VV", 
                        "TTLL_powheg": "t#bar{t}",
                      } 
        else :
            bkgList = {"WJets_MG": "WJets", 
                       "WW_pythia": "EWK", "WZ_pythia": "EWK", "ZZ_pythia": "EWK", 
                       "DYJets10to50ToTauTau":"EWK", "DYJetsToTauTau":"EWK", 
                       "TTLL_powheg": "Top"}
        
        for bkg in bkgList.items():
            hist_filekey_temp = "hist"
            histPostfix_temp = sysPostfix

            if "Unfolding" == sysName :
                histPostfix_temp = ""

            elif "Background" in sysName :
                histPostfix_temp = "" 
    
            _, histType = bkg
            if histType == "Fake" :
                histPostfix_temp = ""
                dirName = "Detector"

            self.unfold_pt.subBkgs(self.inHistDic[hist_filekey_temp], bkg, isSys, self.binDef, dirName, sysName, sysPostfix, histPostfix_temp)
            self.unfold_mass.subBkgs(self.inHistDic[hist_filekey_temp], bkg, isSys, self.binDef, dirName, sysName, sysPostfix, histPostfix_temp)
            
    def setSystematics(self, sysName, sysHistName, isFSR = False):

        self.unfold_pt.setSystematics(sysName, sysHistName)
        self.unfold_mass.setSystematics(sysName, sysHistName)

        matrix_filekey_temp = self.matrix_filekey
        dirPath_temp = self.matrix_dirPath
        histPostfix_temp = sysHistName

        if "Unfolding" in sysName :
            histPostfix_temp = ""

        elif "Unfold" in sysName :
            histPostfix_temp = ""
            if sysHistName != "Nominal" :
                matrix_filekey_temp = "fsr_matrix_DYMG"

        elif "ZptCorr" in sysName :
            histPostfix_temp = ""
            if sysHistName != "Nominal" :
                matrix_filekey_temp = "matrix_zptcorr"

        elif "FSR" in sysName :
            histPostfix_temp = ""

            if "PHOTOS" in sysHistName :
                matrix_filekey_temp = "fsr_matrix_powheg_photos"
            else :
                matrix_filekey_temp = "fsr_matrix_powheg_pythia"

        elif "Background" in sysName :
            histPostfix_temp = ""

        self.unfold_pt.setSystematicRM(self.inHistDic[matrix_filekey_temp], dirPath_temp, self.matrix_histName, sysName, sysHistName, histPostfix_temp, self.binDef)
        self.unfold_mass.setSystematicRM(self.inHistDic[matrix_filekey_temp], dirPath_temp, self.matrix_histName, sysName, sysHistName, histPostfix_temp, self.binDef)

    def drawResponseM(self, var = "Mass", sysName = "", sysPostfix = "", isDetector = True):

        self.unfold.drawResponseM(var, sysName, sysPostfix, isDetector)

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
            self.unfold_pt.doAcceptCorr(self.inHistDic['hist_accept_fullPhase'], self.binDef, outName, True)
            self.unfold_mass.doAcceptCorr(self.inHistDic['hist_accept_fullPhase'], self.binDef, outName, True)
        else : 
            self.unfold_pt.doAcceptCorr(self.inHistDic['hist_accept_drp1'], self.binDef, outName)
            self.unfold_mass.doAcceptCorr(self.inHistDic['hist_accept_drp1'], self.binDef, outName)

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
