import os
import sys
from array import array
import ROOT as rt

import pyScripts.unfoldUtil as unfoldutil

class ISRAnalysis:
    
    def __init__(self, year_ = "2016", channel_= "electron", matrix_filekey_ = "matrix",
                 matrix_dirPath_ = "Detector_Dressed_DRp1_Fiducial", matrix_histName_ = "Detector_Dressed_DRp1"):
        
        # 디텍터 언폴딩 디렉토리 경로 & 매트릭스 이름 
        self.matrix_filekey = matrix_filekey_
        self.matrix_dirPath  = matrix_dirPath_
        self.matrix_histName = matrix_histName_
        
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
        self.dataHistName = "histo_"+dataHistPostfix+"nominal"
            
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
        self.unfold = rt.ISRUnfold(self.channel, int(self.year), int(self.mode))
        self.unfold.setOutputBaseDir(self.outDirPath)
        self.unfold.setBias(self.bias)
        
        # Set response matrix
        self.unfold.setNomResMatrix("Pt", self.inHistDic[self.matrix_filekey], self.matrix_dirPath, self.matrix_histName)
        self.unfold.setNomResMatrix("Mass", self.inHistDic[self.matrix_filekey], self.matrix_dirPath, self.matrix_histName)
        
    def setInputHist(self, useMCInput = False, useUnfoldOut = False, unfoldObj = None):
        
        inputHistName = self.dataHistName
        if useMCInput == True:
            if self.channel == "electron":
                inputHistName = "histo_DYJetsToEEnominal"
            else :
                inputHistName = "histo_DYJetsToMuMunominal"
        
        if useUnfoldOut == False:
            self.unfold.setUnfInput("Pt",   self.inHistDic['hist'], "Detector", inputHistName, False, "nominal", 0)
            self.unfold.setUnfInput("Mass", self.inHistDic['hist'], "Detector", inputHistName, False, "nominal", 0)
        else:
            self.unfold.setUnfInput(unfoldObj, "Pt", False, "", 0)
            self.unfold.setUnfInput(unfoldObj, "Mass", False, "", 0)
            
    def subFake(self):
        
        fakeList = {"DYJets": "DY", self.dy10to50HistName:"DY"}
        
        for fake in fakeList.items():
            self.unfold.subBkgs(self.inHistDic['matrix'], fake, False, "", 0, -1, "detector_level_DY_Fake")
        
    def setUnfoldBkgs(self, doSystematic = False , systName = "nominal", nthSys = 0, nTotSys = -1):
   
        bkgList = {}
        # 2016 데이터만 single top 샘플을 갖고 있다 
        if self.year == "2016" :
            bkgList = {"WJets_MG": "WJets", \
                       "WW_pythia": "EWK", "WZ_pythia": "EWK", "ZZ_pythia": "EWK", \
                       "DYJets10to50ToTauTau":"EWK", "DYJetsToTauTau":"EWK", \
                       "TTLL_powheg": "Top", "SingleTop_tW_top_Incl": "Top", "SingleTop_tW_antitop_Incl": "Top"}
        else :
            bkgList = {"WJets_MG": "WJets", \
                       "WW_pythia": "EWK", "WZ_pythia": "EWK", "ZZ_pythia": "EWK", \
                       "DYJets10to50ToTauTau":"EWK", "DYJetsToTauTau":"EWK", \
                       "TTLL_powheg": "Top"}
        
        for bkg in bkgList.items():
            self.unfold.subBkgs(self.inHistDic['hist'], bkg, doSystematic, systName, nTotSys, nthSys, "Detector")
            
    # Do unfold! 
    def doUnfold(self, doSystematic = False):
        self.unfold.doISRUnfold(doSystematic)

    def doStatUnfold(self):
        self.unfold.doStatUnfold()

    def setMeanValues(self):
        # Set mean mass and pt values
        self.nMassBins = self.unfold.setMeanMass()
        self.unfold.setMeanPt()

    def setStatError(self):
        self.unfold.setStatError()
    
    def getISRUnfold(self):
        
        return self.unfold
    
    def drawStatVar(self, isPt = True):

        for ibin in range(self.nMassBins): 
            self.unfold.drawStatVariation(isPt, ibin)

    # Get histograms
    def getPtVsMassTGraph(self, grTitle = "", isUnfolded = True):
        meanMass, meanPt = array('d'), array('d')
        meanMassStatErr, meanPtStatErr = array('d'), array('d')

        for ibin in range(self.nMassBins):
           
            if isUnfolded:
                meanMass.append(self.unfold.getUnfMeanMass(ibin))
                meanPt.append(self.unfold.getUnfMeanPt(ibin))
                meanMassStatErr.append(self.unfold.getUnfMeanMassError(ibin))
                meanPtStatErr.append(self.unfold.getUnfMeanPtError(ibin))
            else:
                meanMass.append(self.unfold.getDetMeanMass(ibin))
                meanPt.append(self.unfold.getDetMeanPt(ibin))
                meanMassStatErr.append(self.unfold.getDetMeanMassError(ibin))
                meanPtStatErr.append(self.unfold.getDetMeanPtError(ibin))
            
        gr = rt.TGraphErrors(self.nMassBins, meanMass, meanPt, meanMassStatErr, meanPtStatErr)
        gr.SetName(grTitle)
        return gr

detUnfold_muon_2016 = ISRAnalysis("2016", "muon")
detUnfold_muon_2017 = ISRAnalysis("2017", "muon")
detUnfold_muon_2018 = ISRAnalysis("2018", "muon")

fsrUnfold_muon_2016 = ISRAnalysis("2016", "muon", "fsr_matrix", 
                                      "Dressed_DRp1_Dressed_DR4PI_FullPhase", "Dressed_DRp1_Dressed_DR4PI")

fsrUnfold_muon_2017 = ISRAnalysis("2017", "muon", "fsr_matrix", 
                                      "Dressed_DRp1_Dressed_DR4PI_FullPhase", "Dressed_DRp1_Dressed_DR4PI")

fsrUnfold_muon_2018 = ISRAnalysis("2018", "muon", "fsr_matrix", 
                                      "Dressed_DRp1_Dressed_DR4PI_FullPhase", "Dressed_DRp1_Dressed_DR4PI")

detUnfold_muon_2016.setInputHist()
detUnfold_muon_2016.setUnfoldBkgs()

detUnfold_muon_2017.setInputHist()
detUnfold_muon_2017.setUnfoldBkgs()

detUnfold_muon_2018.setInputHist()
detUnfold_muon_2018.setUnfoldBkgs()

detUnfold_muon_2016.subFake()
detUnfold_muon_2017.subFake()
detUnfold_muon_2018.subFake()

detUnfold_muon_2016.doUnfold()
detUnfold_muon_2017.doUnfold()
detUnfold_muon_2018.doUnfold()

detUnfold_muon_2016.setMeanValues()
detUnfold_muon_2017.setMeanValues()
detUnfold_muon_2018.setMeanValues()

detUnfold_muon_2016.doStatUnfold()
detUnfold_muon_2017.doStatUnfold()
detUnfold_muon_2018.doStatUnfold()

detUnfold_muon_2016.setStatError()
detUnfold_muon_2017.setStatError()
detUnfold_muon_2018.setStatError()

gr_2016_muon_det_unfold =  detUnfold_muon_2016.getPtVsMassTGraph("2016MuonDetUnf")
gr_2017_muon_det_unfold =  detUnfold_muon_2017.getPtVsMassTGraph("2017MuonDefUnf")
gr_2018_muon_det_unfold =  detUnfold_muon_2018.getPtVsMassTGraph("2018MuonDefUnf")

gr_2016_muon_det =  detUnfold_muon_2016.getPtVsMassTGraph("2016MuonDet", False)
gr_2017_muon_det =  detUnfold_muon_2017.getPtVsMassTGraph("2017MuonDet", False)
gr_2018_muon_det =  detUnfold_muon_2018.getPtVsMassTGraph("2018MuonDet", False)

fsrUnfold_muon_2016.setInputHist(False, True, detUnfold_muon_2016.getISRUnfold())
fsrUnfold_muon_2017.setInputHist(False, True, detUnfold_muon_2017.getISRUnfold())
fsrUnfold_muon_2018.setInputHist(False, True, detUnfold_muon_2018.getISRUnfold())

fsrUnfold_muon_2016.doUnfold()
fsrUnfold_muon_2017.doUnfold()
fsrUnfold_muon_2018.doUnfold()

fsrUnfold_muon_2016.setMeanValues()
fsrUnfold_muon_2017.setMeanValues()
fsrUnfold_muon_2018.setMeanValues()

fsrUnfold_muon_2016.doStatUnfold()
fsrUnfold_muon_2017.doStatUnfold()
fsrUnfold_muon_2018.doStatUnfold()

fsrUnfold_muon_2016.setStatError()
fsrUnfold_muon_2017.setStatError()
fsrUnfold_muon_2018.setStatError()

gr_2016_muon_fsr_unfold = fsrUnfold_muon_2016.getPtVsMassTGraph("2016MuonFSRUnf")
gr_2017_muon_fsr_unfold = fsrUnfold_muon_2017.getPtVsMassTGraph("2017MuonFSRUnf")
gr_2018_muon_fsr_unfold = fsrUnfold_muon_2018.getPtVsMassTGraph("2018MuonFSRUnf")

c_PtVsMass = rt.TCanvas("PtVsMass","PtVsMass", 1000, 600)
c_PtVsMass.SetGridx()
c_PtVsMass.SetGridy()
c_PtVsMass.SetLogx()

gr_2016_muon_det_unfold.SetTitle("2016, 2017, 2018 muon")
gr_2016_muon_det_unfold.Draw("AP")
gr_2017_muon_det_unfold.Draw("P same")
gr_2018_muon_det_unfold.Draw("P same")
gr_2016_muon_det_unfold.SetMarkerStyle(20)
gr_2016_muon_det_unfold.SetMarkerSize(1)
gr_2017_muon_det_unfold.SetMarkerStyle(25)
gr_2017_muon_det_unfold.SetMarkerSize(1)
gr_2018_muon_det_unfold.SetMarkerStyle(26)
gr_2018_muon_det_unfold.SetMarkerSize(1)

gr_2016_muon_det.Draw("P SAME")
gr_2016_muon_det.SetMarkerStyle(20)
gr_2016_muon_det.SetMarkerSize(1)
gr_2016_muon_det.SetMarkerColor(rt.kBlue)
gr_2016_muon_det.SetLineColor(rt.kBlue)
gr_2017_muon_det.Draw("P SAME")
gr_2017_muon_det.SetMarkerStyle(25)
gr_2017_muon_det.SetMarkerSize(1)
gr_2017_muon_det.SetMarkerColor(rt.kBlue)
gr_2017_muon_det.SetLineColor(rt.kBlue)
gr_2018_muon_det.Draw("P SAME")
gr_2018_muon_det.SetMarkerStyle(26)
gr_2018_muon_det.SetMarkerSize(1)
gr_2018_muon_det.SetMarkerColor(rt.kBlue)
gr_2018_muon_det.SetLineColor(rt.kBlue)

gr_2016_muon_fsr_unfold.Draw("P same")
gr_2016_muon_fsr_unfold.SetMarkerStyle(20)
gr_2016_muon_fsr_unfold.SetMarkerColor(rt.kRed)
gr_2016_muon_fsr_unfold.SetLineColor(rt.kRed)
gr_2016_muon_fsr_unfold.SetMarkerSize(1)
gr_2017_muon_fsr_unfold.Draw("P same")
gr_2017_muon_fsr_unfold.SetMarkerStyle(25)
gr_2017_muon_fsr_unfold.SetMarkerColor(rt.kRed)
gr_2017_muon_fsr_unfold.SetLineColor(rt.kRed)
gr_2017_muon_fsr_unfold.SetMarkerSize(1)
gr_2018_muon_fsr_unfold.Draw("P same")
gr_2018_muon_fsr_unfold.SetMarkerStyle(26)
gr_2018_muon_fsr_unfold.SetMarkerColor(rt.kRed)
gr_2018_muon_fsr_unfold.SetLineColor(rt.kRed)
gr_2018_muon_fsr_unfold.SetMarkerSize(1)

gr_2016_muon_det_unfold.GetYaxis().SetRangeUser(5., 40.)
gr_2016_muon_det_unfold.GetXaxis().SetLimits(10., 500.)

# Linear fit
fitLinear = rt.TF1("f_lepton", "[0]+[1]*log(x)", 10., 400.);
fitLinear.SetLineStyle(2)
fitLinear.SetLineColor(rt.kBlue)
fitLinear.SetLineWidth(1)
gr_2016_muon_fsr_unfold.Fit(fitLinear, "R0")

c_PtVsMass.Draw()

c_PtVsMass.SaveAs("Muon.pdf")

fsrUnfold_muon_2016.drawStatVar()
fsrUnfold_muon_2017.drawStatVar()
fsrUnfold_muon_2018.drawStatVar()

fsrUnfold_muon_2016.drawStatVar(False)
fsrUnfold_muon_2017.drawStatVar(False)
fsrUnfold_muon_2018.drawStatVar(False)
