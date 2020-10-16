import os
import sys
from array import array
import ROOT as rt
import math

import pyScripts.unfoldUtil as unfoldutil
import pandas as pd

class ISRAnalysis:
    
    def __init__(self, year_ = "2016", channel_= "electron", regMode = 0, sys_ = False, matrix_filekey_ = "matrix",
                 matrix_dirPath_ = "Detector_Dressed_DRp1_Fiducial", matrix_histName_ = "Detector_Dressed_DRp1", binDef_ = ""):
        
        # 디텍터 언폴딩 디렉토리 경로 & 매트릭스 이름 
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
        
        # 아웃풋 디렉토리
        self.outDirPath = "output/"+self.year+"/new_"+self.channel+"/"
        # 인풋 히스토그램 텍스트파일
        self.inHistPathTxt = "inFiles/"+self.year+"/"+self.channel+"/fhist.txt"
    
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
        self.unfold = rt.ISRUnfold(self.channel, int(self.year), int(self.mode), sys_)
        self.unfold.setOutputBaseDir(self.outDirPath)
        self.unfold.setBias(self.bias)
        
        # Set response matrix
        self.unfold.setNominalRM("Pt", self.inHistDic[self.matrix_filekey], self.matrix_dirPath, self.matrix_histName, self.binDef)
        self.unfold.setNominalRM("Mass", self.inHistDic[self.matrix_filekey], self.matrix_dirPath, self.matrix_histName, self.binDef)
       
    def checkMatrixCond(self, var = "Mass"):
        return self.unfold.checkMatrixCond(var)

    def setInputHist(self, useMCInput = False, useUnfoldOut = False, unfoldObj = None, dirName = "Detector", isSys = False, sysName = "nominal", sysPostfix = "", isFSR = False, useMadgraph = False):
        
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
 
                self.unfold.setUnfInput("Pt",   self.binDef, self.inHistDic[hist_filekey_temp], dirName, inputHistName, isSys, sysName, sysPostfix)
                self.unfold.setUnfInput("Mass", self.binDef, self.inHistDic[hist_filekey_temp], dirName, inputHistName, isSys, sysName, sysPostfix)
            else :
                inputHistName = "histo_DYJets"

                if useMadgraph == False:
                    hist_filekey_temp = "hist_accept_drp1"

                else:
                    hist_filekey_temp = "hist_accept_drp1_DYMG"

                self.unfold.setUnfInput("Pt",   "Gen_FineCoarse", self.inHistDic[hist_filekey_temp], dirName, inputHistName, isSys, sysName, sysPostfix, isFSR)
                self.unfold.setUnfInput("Mass", "Gen_FineCoarse", self.inHistDic[hist_filekey_temp], dirName, inputHistName, isSys, sysName, sysPostfix, isFSR)

        else :
            if useUnfoldOut == False:

                if "LepScale" in sysName :
                    hist_filekey_temp = "hist_LepScale"
                    dirName = dirName + "_" + sysPostfix  # FIXME put systematic postfix at the END, ex) Detector_Dressed_DRp1_Fiducial_LepScaleUp 

                elif "LepRes" in sysName :
                    hist_filekey_temp = "hist_LepRes"
                    dirName = dirName + "_" + sysPostfix  # FIXME put systematic postfix at the END, ex) Detector_Dressed_DRp1_Fiducial_LepScaleUp 
            
                self.unfold.setUnfInput("Pt",   self.binDef, self.inHistDic[hist_filekey_temp], dirName, inputHistName, isSys, sysName, sysPostfix)
                self.unfold.setUnfInput("Mass", self.binDef, self.inHistDic[hist_filekey_temp], dirName, inputHistName, isSys, sysName, sysPostfix)
            else:
                # Set unfold input from previous unfold 
                self.unfold.setUnfInput(unfoldObj, "Pt",   isSys, sysName, sysPostfix, True)
                self.unfold.setUnfInput(unfoldObj, "Mass", isSys, sysName, sysPostfix, True)

    def setFromPreviousUnfold(self, unfoldObj) :
        self.unfold.setFromPrevUnfResult(unfoldObj, True)
            
    def subFake(self, isSys = False, dirName = "Detector_DY_Fake", sysName = "nominal", sysPostfix = "", histPostfix = "", isFSR = False):
            
        fakeList = {"DYJets": "DY", self.dy10to50HistName:"DY"}
       
        hist_filekey_temp = "matrix"
        if isFSR :
            hist_filekey_temp = "fsr_matrix"

        for fake in fakeList.items():

            if "LepScale" in sysName and not isFSR:
                hist_filekey_temp ="matrix_LepScale"
                dirName = "Detector_DY_Fake" + "_" + sysPostfix

            elif "LepRes" in sysName and not isFSR:
                hist_filekey_temp ="matrix_LepRes"
                dirName = "Detector_DY_Fake" + "_" + sysPostfix

            elif "ZptCorr" in sysName : 
                if sysPostfix != "Nominal" : 
                    hist_filekey_temp = 'matrix_zptcorr' 

            elif "FSR" in sysName : 
                histPostfix = ""

            elif "PDF" in sysName :
                if isFSR :
                    hist_filekey_temp = "fsr_matrix_pdf"  
                else :
                    hist_filekey_temp = "matrix_pdf"

            elif "Background" in sysName : 
                histPostfix = ""

            self.unfold.subBkgs(self.inHistDic[hist_filekey_temp], fake, isSys, self.binDef, dirName, sysName, sysPostfix, histPostfix)
        
    def setUnfoldBkgs(self, isSys = False , dirName = "Detector", sysName = "nominal", sysPostfix = ""):
   
        bkgList = {}
        # 2016 데이터만 single top 샘플을 갖고 있다 
        if self.year == "2016" or self.year == "2017" or self.year == "2018":
            bkgList = { "QCD": "Fake", "WJet": "Fake",
            #bkgList = {
            #            "WJets_MG": "WJets",
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

            if "LepScale" in sysName :
                hist_filekey_temp = "hist_LepScale"
                dirName = "Detector" + "_" + sysPostfix
                histPostfix_temp = ""

            elif "LepRes" in sysName :
                hist_filekey_temp = "hist_LepRes"
                dirName = "Detector" + "_" + sysPostfix
                histPostfix_temp = ""

            elif "iterEM" == sysName :
                histPostfix_temp = ""

            elif "PDF" in sysName :
                hist_filekey_temp ="hist_pdf"

            elif "Background" in sysName :
                histPostfix_temp = "" 
    
            _, histType = bkg
            if histType == "Fake" :
                histPostfix_temp = ""
                dirName = "Detector"

            self.unfold.subBkgs(self.inHistDic[hist_filekey_temp], bkg, isSys, self.binDef, dirName, sysName, sysPostfix, histPostfix_temp)
            
    def setSystematics(self, sysName, sysHistName, isFSR = False):
        self.unfold.setSystematics(sysName, sysHistName)

        matrix_filekey_temp = self.matrix_filekey
        dirPath_temp = self.matrix_dirPath
        histPostfix_temp = sysHistName

        if "LepScale" in sysName :
            matrix_filekey_temp = "matrix_LepScale"
            dirPath_temp = "Detector_Dressed_DRp1_Fiducial" + "_" + sysHistName  # FIXME put systematic postfix at the END, ex) Detector_Dressed_DRp1_Fiducial_LepScaleUp 
            histPostfix_temp = ""

        elif "LepRes" in sysName :
            matrix_filekey_temp = "matrix_LepRes"
            dirPath_temp = "Detector_Dressed_DRp1_Fiducial" + "_" + sysHistName  # FIXME put systematic postfix at the END, ex) Detector_Dressed_DRp1_Fiducial_LepScaleUp 
            histPostfix_temp = ""

        elif "iterEM" in sysName :
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

            if sysHistName == "PHOTOS" :
                matrix_filekey_temp = "fsr_matrix_powheg_photos"
            else :
                matrix_filekey_temp = "fsr_matrix_powheg_pythia"

        elif "PDF" in sysName:
            if isFSR :
                matrix_filekey_temp = "fsr_matrix_pdf"
            else :
                matrix_filekey_temp = "matrix_pdf"

        elif "Background" in sysName :
            histPostfix_temp = ""

        self.unfold.setSystematicRM("Pt",   self.inHistDic[matrix_filekey_temp], dirPath_temp, self.matrix_histName, sysName, sysHistName, histPostfix_temp, self.binDef)
        self.unfold.setSystematicRM("Mass", self.inHistDic[matrix_filekey_temp], dirPath_temp, self.matrix_histName, sysName, sysHistName, histPostfix_temp, self.binDef)

    def getSystematics(self):
        self.unfold.printSystematics()

    def printMeanValues(self):
        self.unfold.printMeanValues(True)

    def printMeanValues_Accept(self):
        self.unfold.printMeanValues_Accept(True)

    def drawResponseM(self, var = "Mass", sysName = "", sysPostfix = "", isDetector = True):
        self.unfold.drawResponseM(var, sysName, sysPostfix, isDetector)

    def checkIterEMUnfold(self):
        self.unfold.checkIterEMUnfold()

    def drawDetPlot(self, var = "Mass", dirName = "Detector", steering = None, useAxis = True, sysName = "", outName = "", massBin = 0, binWidth = False, isBkgSubData = False):
        if "LepScale" in sysName :
            self.unfold.drawFoldedHists(var, self.inHistDic['hist'], dirName, steering, useAxis, sysName, outName, massBin, binWidth, self.inHistDic['hist_LepScale'], isBkgSubData)
        elif "LepRes" in sysName :
            self.unfold.drawFoldedHists(var, self.inHistDic['hist'], dirName, steering, useAxis, sysName, outName, massBin, binWidth, self.inHistDic['hist_LepRes'], isBkgSubData)
        elif "Scale" in sysName or "AlphaS" in sysName:
            self.unfold.drawFoldedHists(var, self.inHistDic['hist'], dirName, steering, useAxis, sysName, outName, massBin, binWidth, self.inHistDic['hist'], isBkgSubData)
        elif "PDF" in sysName :
            self.unfold.drawFoldedHists(var, self.inHistDic['hist'], dirName, steering, useAxis, sysName, outName, massBin, binWidth, self.inHistDic['hist_pdf'], isBkgSubData)
        else : 
            return self.unfold.drawFoldedHists(var, self.inHistDic['hist'], dirName, steering, useAxis, sysName, outName, massBin, binWidth, "", isBkgSubData)

    def drawPtDistributions(self) :
            self.unfold.drawPtDistributions(self.inHistDic['hist'])

    def drawPtBkgRatio(self) :
            self.unfold.drawPtBkgRatio(self.inHistDic['hist'])

    def drawUnfPlot(self, var = "Mass", steering = None, useAxis = True, sysName = "", outName = "", massBin = 0, binWidth = False, isType3Closure = False):
        self.unfold.drawUnfoldedHists(var, steering, useAxis, sysName, outName, massBin, binWidth, isType3Closure)

    def drawUnfVarPlot(self, var = "Mass", steering = None, useAxis = True, sysName = "", outName = "", massBin = 0, binWidth = False):
        self.unfold.drawUnfoldedVarHists(var, steering, useAxis, sysName, outName, massBin, binWidth)

    def drawAcceptVarPlot(self, var = "Mass", steering = None, useAxis = True, sysName = "", outName = "", massBin = 0, isAccept = False):
        self.unfold.drawAcceptVarHists(var, steering, useAxis, sysName, outName, massBin, isAccept)

    def drawSystematics(self, var = "Pt", isAccept = False, isHistStyle = False) :
        if isAccept:
            self.unfold.drawSystematics_Acceptance(var, isHistStyle)
        else :
            self.unfold.drawSystematics(var, isHistStyle) 
    # Do unfold! 
    def doUnfold(self, isSys = False):
        self.unfold.doISRUnfold(isSys)

    def doStatUnfold(self):
        print ("doStatUnfold")
        self.unfold.doStatUnfold()

    def setMeanValues(self, setDetector = False):
        # Set mean mass and pt values
        if setDetector :
            self.nMassBins = self.unfold.setMeanMass(self.inHistDic['hist'])
            self.unfold.setMeanPt(self.inHistDic['hist'])
        else :
            self.nMassBins = self.unfold.setMeanMass()
            self.unfold.setMeanPt()

    def setTheoryMeanValues(self, binDef = "_FineCoarse"):
        self.unfold.setTheoryMeanValues(self.inHistDic['hist_accept_fullPhase'], binDef)

    def setAcceptMeanValues(self):
        self.unfold.setMeanPt_Accept()
        self.unfold.setMeanMass_Accept()

    def setSysMeanValues(self):
        self.unfold.setSysMeanMass()
        self.unfold.setSysMeanPt()

    def setAcceptSysMeanValues(self):
        self.unfold.setSysMeanMass_Accept()
        self.unfold.setSysMeanPt_Accept()

    def setStatError(self):
        self.unfold.setStatError()

    def setSysError(self):
        self.unfold.setSysError()

    def setAcceptSysError(self):
        self.unfold.setSysError_Accept()

    def setTotSysError(self):
        self.unfold.setTotSysError()

    def setAcceptTotSysError(self):
        self.unfold.setTotSysError_Accept()
    
    def getISRUnfold(self):
        
        return self.unfold
    
    def drawStatVar(self, isPt = True):

        for ibin in range(self.nMassBins): 
            self.unfold.drawStatVariation(isPt, ibin)

    def drawPDFVar(self, isPt = True):

        for ibin in range(self.nMassBins): 
            self.unfold.drawPDFVariation(isPt, ibin)

    def drawSysVar(self, sysName, var = "Pt", isAccept = False):

        for ibin in range(self.nMassBins): 
            if isAccept == False:
                self.unfold.drawSysVariation(sysName, var, ibin)
            else :
                self.unfold.drawSysVariation_Accept(sysName, var, ibin)

    def doAcceptance(self, doSys = False, isFSR = False, outName = "") :
        if isFSR :
            self.unfold.doAcceptCorr(self.inHistDic['hist_accept_fullPhase'], "_FineCoarse", doSys, outName, True)
        else : 
            self.unfold.doAcceptCorr(self.inHistDic['hist_accept_drp1'], "_FineCoarse", doSys, outName)

    def drawAcceptPlot(self, var = "Mass", steering = None, useAxis = True, sysName = "", outName = "", massBin = 0, binWidth = False, isFSR = False): 
        if isFSR :
            self.unfold.drawAcceptCorrHists(var, self.inHistDic['hist_accept_fullPhase'], "_FineCoarse", steering, useAxis, sysName, outName, massBin, binWidth)  
        else :
            self.unfold.drawAcceptCorrHists(var, self.inHistDic['hist_accept_drp1'], "_FineCoarse", steering, useAxis, sysName, outName, massBin, binWidth)  

    def drawCorrelation(self, var = "Mass", steering = None, useAxis = True, outName = ""):
        self.unfold.drawCorrelation(var, steering, useAxis, outName)

    def getRawHist(self, var = "Mass", histName = "", steering = None, useAxis = True, sysName = "", outName = "", massBin = 0, binWidth = False):
        self.unfold.getRawHist(var, self.inHistDic['hist'], "Detector", histName, outName, steering, useAxis, binWidth)

    def getUnfInHist(self, var = "Mass", steering = None, useAxis = True, sysName = "", outName = "", massBin = 0, binWidth = False):
        return self.unfold.getUnfInput(var, steering, useAxis, massBin, binWidth)

    def getUnfHist(self, var = "Mass", outName = "" , steering = None, useAxis = True, binWidth = False):
        return self.unfold.getUnfoldedHists(var, outName, steering, useAxis, binWidth)

    def getGenMCHist(self, var = "Mass", steering = None, useAxis = True, sysName = "", outName = "", massBin = 0, binWidth = False):
        return self.unfold.getGenMCHist(var, steering, useAxis, massBin, binWidth)

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
        
    def getMeanVectors(self, whichLevel="Unfolded", var="Mass", sysName="", variationName="") :

        if whichLevel is "Unfolded" : 
            if var is "Mass" :
                return self.unfold.getUnfoldedMeanMassVectors(sysName, variationName)     
            if var is "Pt" :
                return self.unfold.getUnfoldedMeanPtVectors(sysName, variationName)     

        if whichLevel is "Acceptance" :
            if var is "Mass" :
                return self.unfold.getAccCorrectedMeanMassVectors(sysName, variationName)     
            if var is "Pt" :
                return self.unfold.getAccCorrectedMeanPtVectors(sysName, variationName)     

    def makeMeanDataFrames(self, whichLevel="Unfolded", systematics=None) :

        nominal_mass_vector=self.getMeanVectors(whichLevel, var="Mass")
        nominal_pt_vector=self.getMeanVectors(whichLevel, var="Pt")
        
        nBin=len(nominal_mass_vector)
        pd_series_binIndex=pd.Series(range(1,nBin+1), range(1, nBin+1), name="bin_index")
        
        pd_mass_mean=pd.DataFrame(index=pd_series_binIndex)
        pd_pt_mean=pd.DataFrame(index=pd_series_binIndex)

        pd_mass_mean["Nominal"]=nominal_mass_vector
        pd_pt_mean["Nominal"]=nominal_pt_vector
       
        for sysName, postfixs in systematics.items() :
            for postfix in postfixs :
                if postfix is "Nominal" : continue

                temp_mass_vector=self.getMeanVectors(whichLevel, "Mass", sysName, postfix)
                temp_pt_vector=self.getMeanVectors(whichLevel, "Pt", sysName, postfix)
                pd_mass_mean[postfix]=temp_mass_vector
                pd_pt_mean[postfix]=temp_pt_vector

        return pd_mass_mean, pd_pt_mean

    def getPtVsMassTGraph(self, grTitle = "", isData = True, whichLevel = "Detector", doSys = False):
        meanMass, meanPt = array('d'), array('d')
        meanMassStatErr, meanPtStatErr = array('d'), array('d')
        meanMassSysErr, meanPtSysErr = array('d'), array('d')
   
        fileName = "2016.csv" 
        import csv
        import os.path
        from os import path
        if not path.exists(fileName) :
            with open(fileName, 'w', newline='') as csvfile:
                writeMeasurement = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                writeMeasurement.writerow(["Mean mass", "Mean pT", "Mean mass uncertatinty", "Mean pT uncertainty", "Year", "Channel", "Level", "Sample"])

        for ibin in range(self.nMassBins):
            tempRow = list()
          
            if isData: 
                if whichLevel == "Unfolded":
                    meanMass.append(self.unfold.getUnfMeanMass(ibin))
                    meanPt.append(self.unfold.getUnfMeanPt(ibin))
                    meanMassStatErr.append(self.unfold.getUnfMeanMassError(ibin))
                    meanPtStatErr.append(self.unfold.getUnfMeanPtError(ibin))
                    if doSys :
                        meanMassSysErr.append(self.unfold.getUnfMeanMassSysError(ibin))
                        meanPtSysErr.append(self.unfold.getUnfMeanPtSysError(ibin))

                elif whichLevel == "Acceptance":
                    meanMass.append(self.unfold.getAccMeanMass(ibin))
                    meanPt.append(self.unfold.getAccMeanPt(ibin))
                    meanMassStatErr.append(self.unfold.getAccMeanMassError(ibin))
                    meanPtStatErr.append(self.unfold.getAccMeanPtError(ibin))
                    if doSys :
                        meanMassSysErr.append(self.unfold.getAccMeanMassTotError(ibin))
                        meanPtSysErr.append(self.unfold.getAccMeanPtTotError(ibin))

                elif whichLevel == "Detector":
                    meanMass.append(self.unfold.getDetMeanMass(ibin))
                    meanPt.append(self.unfold.getDetMeanPt(ibin))
                    meanMassStatErr.append(self.unfold.getDetMeanMassError(ibin))
                    meanPtStatErr.append(self.unfold.getDetMeanPtError(ibin))
                else :
                    print("Please check phase space name")

            if not isData:
                if whichLevel == "Detector":
                    meanMass.append(self.unfold.getMCDetMeanMass(ibin))
                    meanPt.append(self.unfold.getMCDetMeanPt(ibin))
                    meanMassStatErr.append(self.unfold.getMCDetMeanMassError(ibin))
                    meanPtStatErr.append(self.unfold.getMCDetMeanPtError(ibin))
                    #print("detector mc check pt", self.unfold.getMCDetMeanPt(ibin))
                    #print("detector mc check mass", self.unfold.getMCDetMeanMass(ibin))

            if doSys :
                tempRow.extend(list((meanMass[ibin], meanPt[ibin], meanMassSysErr[ibin], meanPtSysErr[ibin])))
                tempRow.extend(list((self.year, self.channel, whichLevel, "Data")))        

                with open(fileName, 'a', newline='') as csvfile: 
                    writeMeasurement = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                    writeMeasurement.writerow(tempRow)
 
        if doSys == False:
            gr = rt.TGraphErrors(self.nMassBins, meanMass, meanPt, meanMassStatErr, meanPtStatErr)
        else :
            gr = rt.TGraphErrors(self.nMassBins, meanMass, meanPt, meanMassSysErr, meanPtSysErr)
    
        gr.SetName(grTitle)
        return gr

    def getThoryPtVsMassTGraph(self, grTitle = "", sysName = "stat"):
        meanMass, meanPt       = array('d'), array('d')
        meanMassErr, meanPtErr = array('d'), array('d')

        for ibin in range(self.nMassBins):
            meanMass.append(self.unfold.getMCFullPhaseMeanMass(ibin, ""))
            meanMassErr.append(self.unfold.getMCFullPhaseMeanMass(ibin, sysName))

            meanPt.append(self.unfold.getMCFullPhaseMeanPt(ibin, ""))
            meanPtErr.append(self.unfold.getMCFullPhaseMeanPt(ibin, sysName))
 
        gr = rt.TGraphErrors(self.nMassBins, meanMass, meanPt, meanMassErr, meanPtErr)
        gr.SetName(grTitle)
        return gr
