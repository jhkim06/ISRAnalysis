import os
import sys
from array import array
import ROOT as rt
import math

import pyScripts.unfoldUtil as unfoldutil

class ISRAnalysis:
    
    def __init__(self, year_ = "2016", channel_= "electron", sys_ = False, matrix_filekey_ = "matrix",
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
        dataHistPostfix = "Double"
        if self.channel == "electron":
            dataHistPostfix = dataHistPostfix + "EG"
            if self.year == "2018":
                dataHistPostfix = "EGamma"
        if self.channel == "muon":
            dataHistPostfix = dataHistPostfix + "Muon"
        self.dataHistName = "histo_"+dataHistPostfix
            
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
        self.mode = 0 # 레귤라이제이션 모드
        
        # Create ISRUnfold object
        self.unfold = rt.ISRUnfold(self.channel, int(self.year), int(self.mode), sys_)
        self.unfold.setOutputBaseDir(self.outDirPath)
        self.unfold.setBias(self.bias)
        
        # Set response matrix
        self.unfold.setNomResMatrix("Pt", self.inHistDic[self.matrix_filekey], self.matrix_dirPath, self.matrix_histName, self.binDef)
        self.unfold.setNomResMatrix("Mass", self.inHistDic[self.matrix_filekey], self.matrix_dirPath, self.matrix_histName, self.binDef)
       
    def checkMatrixCond(self, var = "Mass"):
        return self.unfold.checkMatrixCond(var)

    def setInputHist(self, useMCInput = False, useUnfoldOut = False, unfoldObj = None, dirName = "Detector", isSys = False, sysName = "nominal", sysPostfix = "", isFSR = False, useMadgraph = False):
        
        inputHistName = self.dataHistName
        if useMCInput == True:
            #
            temp_filePath = self.inHistDic['hist'] 

            if isFSR == False:
                if self.channel == "electron":
                    inputHistName = "histo_DYJetsToEE"
                else :
                    inputHistName = "histo_DYJetsToMuMu"

                if useMadgraph == False:

                    self.unfold.setUnfInput("Pt",   self.binDef, temp_filePath, dirName, inputHistName, isSys, sysName, sysPostfix)
                    self.unfold.setUnfInput("Mass", self.binDef, temp_filePath, dirName, inputHistName, isSys, sysName, sysPostfix)
                else :
                    temp_filePath = self.inHistDic['hist_DYMG']
 
                    self.unfold.setUnfInput("Pt",   self.binDef, temp_filePath, dirName, inputHistName, isSys, sysName, sysPostfix)
                    self.unfold.setUnfInput("Mass", self.binDef, temp_filePath, dirName, inputHistName, isSys, sysName, sysPostfix)
            else :
                inputHistName = "histo_DYJets"

                temp_filePath = self.inHistDic['hist_accept_drp1']

                self.unfold.setUnfInput("Pt",   "Gen_FineCoarse", temp_filePath, dirName, inputHistName, isSys, sysName, sysPostfix, isFSR)
                self.unfold.setUnfInput("Mass", "Gen_FineCoarse", temp_filePath, dirName, inputHistName, isSys, sysName, sysPostfix, isFSR)

        else :
            if useUnfoldOut == False:
                if "LepMom" in sysName :
                    if sysPostfix == "Nominal" :
                        self.unfold.setUnfInput("Pt",   self.binDef, self.inHistDic['hist'], dirName, inputHistName, isSys, sysName, sysPostfix)
                        self.unfold.setUnfInput("Mass", self.binDef, self.inHistDic['hist'], dirName, inputHistName, isSys, sysName, sysPostfix)
                    else :
                        self.unfold.setUnfInput("Pt",   self.binDef, self.inHistDic['hist_lepScale'], dirName+"_"+sysPostfix, inputHistName, isSys, sysName, sysPostfix)
                        self.unfold.setUnfInput("Mass", self.binDef, self.inHistDic['hist_lepScale'], dirName+"_"+sysPostfix, inputHistName, isSys, sysName, sysPostfix)
                else :
                    self.unfold.setUnfInput("Pt",   self.binDef, self.inHistDic['hist'], dirName, inputHistName, isSys, sysName, sysPostfix)
                    self.unfold.setUnfInput("Mass", self.binDef, self.inHistDic['hist'], dirName, inputHistName, isSys, sysName, sysPostfix)
            else:
                # Set unfold input from previous unfold 
                self.unfold.setUnfInput(unfoldObj, "Pt", isSys, sysName, sysPostfix, True)
                self.unfold.setUnfInput(unfoldObj, "Mass", isSys, sysName, sysPostfix, True)

    def setFromPreviousUnfold(self, unfoldObj) :
        self.unfold.setFromPrevUnfResult(unfoldObj, True)
            
    def subFake(self, isSys = False, systName = "nominal", sysPostfix = "", isFSR = False, dirName = "detector_level_DY_Fake"):
            
        fakeList = {"DYJets": "DY", self.dy10to50HistName:"DY"}
       
        if not isFSR: 
            for fake in fakeList.items():
                if "LepMom" in systName :
                    if sysPostfix == "Nominal" : 
                        self.unfold.subBkgs(self.inHistDic['matrix'], fake, isSys, self.binDef, dirName, systName, sysPostfix, isFSR)
                    else :
                        self.unfold.subBkgs(self.inHistDic['matrix_lepScale'], fake, isSys, self.binDef, dirName+"_"+sysPostfix, systName, sysPostfix, isFSR)
                elif "ZptCorr" in systName : 
                    if sysPostfix == "Nominal" : 
                        self.unfold.subBkgs(self.inHistDic['matrix'], fake, isSys, self.binDef, dirName, systName, sysPostfix, isFSR)
                    else :
                        self.unfold.subBkgs(self.inHistDic['matrix_zptcorr'], fake, isSys, self.binDef, dirName, systName, sysPostfix, isFSR)
                else :
                    self.unfold.subBkgs(self.inHistDic['matrix'], fake, isSys, self.binDef, dirName, systName, sysPostfix, isFSR)
        else:
            for fake in fakeList.items():
                self.unfold.subBkgs(self.inHistDic['fsr_matrix'], fake, isSys, self.binDef, dirName, systName, sysPostfix, isFSR)

        
    def setUnfoldBkgs(self, doSystematic = False , dirName = "Detector",systName = "nominal", sysPostfix = ""):
   
        bkgList = {}
        # 2016 데이터만 single top 샘플을 갖고 있다 
        if self.year == "2016" or self.year == "2017" or self.year == "2018":
            #bkgList = {"QCD": "Fake", "WJet": "Fake",\
            bkgList = {
                        "WJets_MG": "WJets",
                        "DYJets10to50ToTauTau":"EWK", "DYJetsToTauTau":"EWK", 
                        "WW_pythia": "EWK", "WZ_pythia": "EWK", "ZZ_pythia": "EWK", 
                        "TTLL_powheg": "Top",
                        } 
        else :
            bkgList = {"WJets_MG": "WJets", 
                       "WW_pythia": "EWK", "WZ_pythia": "EWK", "ZZ_pythia": "EWK", 
                       "DYJets10to50ToTauTau":"EWK", "DYJetsToTauTau":"EWK", 
                       "TTLL_powheg": "Top"}
        
        for bkg in bkgList.items():
            # TODO If sysPostfix is "Nominal", then use self.inHistDic[self.matrix_filekey], self.matrix_dirPath
            if "LepMom" in systName :
                if sysPostfix == "Nominal" :
                    self.unfold.subBkgs(self.inHistDic['hist'], bkg, doSystematic, self.binDef, dirName, systName, sysPostfix)
                else : 
                    self.unfold.subBkgs(self.inHistDic['hist_lepScale'], bkg, doSystematic, self.binDef, dirName+"_"+sysPostfix, systName, sysPostfix)
            elif "Unfold_det" == systName :
                self.unfold.subBkgs(self.inHistDic['hist'], bkg, doSystematic, self.binDef, dirName, systName, sysPostfix)
            else :
                self.unfold.subBkgs(self.inHistDic['hist'], bkg, doSystematic, self.binDef, dirName, systName, sysPostfix)
            
    def setSystematics(self, sysName, sysHistName):
        self.unfold.setSystematics(sysName, sysHistName)

        if "Closure" in sysName:
            if sysHistName == "Nominal" :
                self.unfold.setSysTUnfoldDensity("Pt",   self.inHistDic["matrix"], self.matrix_dirPath, self.matrix_histName, sysName, sysHistName, self.binDef)
                self.unfold.setSysTUnfoldDensity("Mass", self.inHistDic["matrix"], self.matrix_dirPath, self.matrix_histName, sysName, sysHistName, self.binDef)
        
        elif "LepMom" in sysName :
            if sysHistName == "Nominal" :
                self.unfold.setSysTUnfoldDensity("Pt",   self.inHistDic[self.matrix_filekey], self.matrix_dirPath, self.matrix_histName, sysName, sysHistName, self.binDef)
                self.unfold.setSysTUnfoldDensity("Mass", self.inHistDic[self.matrix_filekey], self.matrix_dirPath, self.matrix_histName, sysName, sysHistName, self.binDef)
            else :
                self.unfold.setSysTUnfoldDensity("Pt",   self.inHistDic["matrix_lepScale"], self.matrix_dirPath+"_"+sysHistName, self.matrix_histName, sysName, sysHistName, self.binDef)
                self.unfold.setSysTUnfoldDensity("Mass", self.inHistDic["matrix_lepScale"], self.matrix_dirPath+"_"+sysHistName, self.matrix_histName, sysName, sysHistName, self.binDef)
        elif "Unfold" in sysName :
            if sysHistName == "Nominal" :
                self.unfold.setSysTUnfoldDensity("Pt",   self.inHistDic[self.matrix_filekey], self.matrix_dirPath, self.matrix_histName, sysName, sysHistName, self.binDef)
                self.unfold.setSysTUnfoldDensity("Mass", self.inHistDic[self.matrix_filekey], self.matrix_dirPath, self.matrix_histName, sysName, sysHistName, self.binDef)
            else :
                self.unfold.setSysTUnfoldDensity("Pt",   self.inHistDic["fsr_matrix_powheg_pythia"], self.matrix_dirPath, self.matrix_histName, sysName, sysHistName, self.binDef)
                self.unfold.setSysTUnfoldDensity("Mass", self.inHistDic["fsr_matrix_powheg_pythia"], self.matrix_dirPath, self.matrix_histName, sysName, sysHistName, self.binDef)

        elif "ZptCorr" in sysName :
            if sysHistName == "Nominal" :
                self.unfold.setSysTUnfoldDensity("Pt",   self.inHistDic[self.matrix_filekey], self.matrix_dirPath, self.matrix_histName, sysName, sysHistName, self.binDef)
                self.unfold.setSysTUnfoldDensity("Mass", self.inHistDic[self.matrix_filekey], self.matrix_dirPath, self.matrix_histName, sysName, sysHistName, self.binDef)
            else :
                self.unfold.setSysTUnfoldDensity("Pt",   self.inHistDic["matrix_zptcorr"], self.matrix_dirPath, self.matrix_histName, sysName, sysHistName, self.binDef)
                self.unfold.setSysTUnfoldDensity("Mass", self.inHistDic["matrix_zptcorr"], self.matrix_dirPath, self.matrix_histName, sysName, sysHistName, self.binDef)

        elif "FSR" in sysName :
            if sysHistName == "PHOTOS" :
                self.unfold.setSysTUnfoldDensity("Pt",   self.inHistDic["fsr_matrix_powheg_photos"], self.matrix_dirPath, self.matrix_histName, sysName, sysHistName, self.binDef)
                self.unfold.setSysTUnfoldDensity("Mass", self.inHistDic["fsr_matrix_powheg_photos"], self.matrix_dirPath, self.matrix_histName, sysName, sysHistName, self.binDef)
            else :
                self.unfold.setSysTUnfoldDensity("Pt",   self.inHistDic["fsr_matrix_powheg_pythia"], self.matrix_dirPath, self.matrix_histName, sysName, sysHistName, self.binDef)
                self.unfold.setSysTUnfoldDensity("Mass", self.inHistDic["fsr_matrix_powheg_pythia"], self.matrix_dirPath, self.matrix_histName, sysName, sysHistName, self.binDef)

        elif "PDF" in sysName:
            self.unfold.setSysTUnfoldDensity("Pt",   self.inHistDic["fsr_matrix_pdf"], self.matrix_dirPath, self.matrix_histName, sysName, sysHistName, self.binDef)
            self.unfold.setSysTUnfoldDensity("Mass", self.inHistDic["fsr_matrix_pdf"], self.matrix_dirPath, self.matrix_histName, sysName, sysHistName, self.binDef)

        else :
            self.unfold.setSysTUnfoldDensity("Pt",   self.inHistDic[self.matrix_filekey], self.matrix_dirPath, self.matrix_histName, sysName, sysHistName, self.binDef)
            self.unfold.setSysTUnfoldDensity("Mass", self.inHistDic[self.matrix_filekey], self.matrix_dirPath, self.matrix_histName, sysName, sysHistName, self.binDef)

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

    def drawDetPlot(self, var = "Mass", dirName = "Detector", steering = None, useAxis = True, sysName = "", outName = "", massBin = 0, binWidth = False):
        if "LepMom" in sysName :
            self.unfold.drawFoldedHists(var, self.inHistDic['hist'], dirName, steering, useAxis, sysName, outName, massBin, binWidth, self.inHistDic['hist_lepScale'])
        elif "ZptCorr" in sysName :
            self.unfold.drawFoldedHists(var, self.inHistDic['hist'], dirName, steering, useAxis, sysName, outName, massBin, binWidth, self.inHistDic['hist_zptcorr'])
        else : 
            self.unfold.drawFoldedHists(var, self.inHistDic['hist'], dirName, steering, useAxis, sysName, outName, massBin, binWidth)

    def drawUnfPlot(self, var = "Mass", steering = None, useAxis = True, sysName = "", outName = "", massBin = 0, binWidth = False, isType3Closure = False):
        self.unfold.drawUnfoldedHists(var, steering, useAxis, sysName, outName, massBin, binWidth, isType3Closure)

    def drawSystematics(self, var = "Pt", isAccept = False) :
        if isAccept:
            self.unfold.drawSystematics_Acceptance(var)
        else :
            self.unfold.drawSystematics(var) 
    # Do unfold! 
    def doUnfold(self, doSystematic = False, doRegularize = False):
        self.unfold.doISRUnfold(doSystematic, doRegularize)

    def doStatUnfold(self):
        print ("doStatUnfold")
        self.unfold.doStatUnfold()

    def setMeanValues(self, setDetector = False):
        # Set mean mass and pt values
        if setDetector :
            self.nMassBins = self.unfold.setMeanMass(self.inHistDic['hist'], "Detector")
            self.unfold.setMeanPt(self.inHistDic['hist'], "Detector")
        else :
            self.nMassBins = self.unfold.setMeanMass()
            self.unfold.setMeanPt()

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

    def drawSysVar(self, sysName, var = "Pt"):

        for ibin in range(self.nMassBins): 
            self.unfold.drawSysVariation(sysName, var, ibin)

    def doAcceptance(self, doSys = False, isFSR = False, outName = "") :
        if isFSR :
            self.unfold.doAcceptCorr(self.inHistDic['hist_accept_fullPhase'], "_FineCoarse", doSys, outName)
        else :
            self.unfold.doAcceptCorr(self.inHistDic['hist_accept_drp1'], "_FineCoarse", doSys, outName)

    def drawAcceptPlot(self, var = "Mass", steering = None, useAxis = True, sysName = "", outName = "", massBin = 0, binWidth = False, isFSR = False): 
        if isFSR :
            self.unfold.drawAcceptCorrHists(var, self.inHistDic['hist_accept_fullPhase'], self.binDef, steering, useAxis, sysName, outName, massBin, binWidth)  
        else : 
            self.unfold.drawAcceptCorrHists(var, self.inHistDic['hist_accept_drp1'], self.binDef, steering, useAxis, sysName, outName, massBin, binWidth)  

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
             

    def getPtVsMassTGraph(self, grTitle = "", isData = True, whichLevel = "Detector", doSys = False):
        meanMass, meanPt = array('d'), array('d')
        meanMassStatErr, meanPtStatErr = array('d'), array('d')
        meanMassSysErr, meanPtSysErr = array('d'), array('d')

        for ibin in range(self.nMassBins):
          
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
            if isData == False:
                if whichLevel == "Detector":
                    meanMass.append(self.unfold.getMCDetMeanMass(ibin))
                    meanPt.append(self.unfold.getMCDetMeanPt(ibin))
                    meanMassStatErr.append(self.unfold.getMCDetMeanMassError(ibin))
                    meanPtStatErr.append(self.unfold.getMCDetMeanPtError(ibin))
                    #print("detector mc check pt", self.unfold.getMCDetMeanPt(ibin))
                    #print("detector mc check mass", self.unfold.getMCDetMeanMass(ibin))
 
        if doSys == False:
            gr = rt.TGraphErrors(self.nMassBins, meanMass, meanPt, meanMassStatErr, meanPtStatErr)
        else :
            gr = rt.TGraphErrors(self.nMassBins, meanMass, meanPt, meanMassSysErr, meanPtSysErr)
    
        gr.SetName(grTitle)
        return gr
