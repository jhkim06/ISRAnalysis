import os
import sys
import ROOT as rt

rt.gSystem.Load("/home/jhkim/ISR_Run2/unfolding/TUnfoldISR2016/lib/libisrunfold.so")

#import etc.histDef as fHistDef

import gc
gc.collect()

from enum import Enum

class ptOrMassEnum(Enum):
	PT   = 1
	MASS = 2

# input file paths, binning, output directory, for data and bkg reco histograms only
def makeRecoHists(sample, outputDirectory, channel):

    rt.gROOT.SetBatch()
    print ("####################### makeRecoHists #################################")

    outDic = {}

    outfile = rt.TFile(outputDirectory + sample.name+".root",'recreate')
    if not sample.isMC: outDic[sample.name] = fHistDef.inputfHists(sample.name, outputDirectory + sample.name+".root", "data")
    else : outDic[sample.name] = fHistDef.inputfHists(sample.name, outputDirectory + sample.name+".root", "bkg")

    # need to create histogram before the file loop to save one histogram from the mutiple input files
    recoHists = rt.histTUnfold(rt.std.map("TString, TH1*")())

    ######################################################################################
    test = rt.histTUnfold(rt.std.map("TString, TH1*")())
    # https://root.cern.ch/faq/how-generate-dictionary
    vrecoHists = rt.std.vector("histTUnfold")()
    vrecoHists.push_back(test)
    vrecoHists.at(0).SetPtBinningRec(channel)
    vrecoHists.at(0).SetMassBinningRec(channel)
    vrecoHists.at(0).CreateHistMap(ptOrMassEnum.PT.value, "test") # 1: pt
    vrecoHists.at(0).CreateHistMap(ptOrMassEnum.MASS.value, "test") # 2: mass
    ######################################################################################

    #recoHists.SetPtBinningRec(channel)
    #recoHists.SetMassBinningRec(channel)

    recoHists.SetPtBinningRec("electron")
    recoHists.SetMassBinningRec("electron")

    recoHists.CreateHistMap(ptOrMassEnum.PT.value, "nominal")
    recoHists.CreateHistMap(ptOrMassEnum.MASS.value, "nominal")

    sysDict = {"PU": 2, "trgSF": 2, "recoSF": 2, "IdSF": 2, "IsoSF": 2, "unfoldsys": 1, "AlphaS": 2, "Scale": 9, "PDFerror": 100, "L1Prefire": 2, "LepScale": 2, "LepRes": 2, "FSRDR": 30}

    for sysName, nSys in sysDict.items():
    	recoHists.SetsysMap(sysName, nSys);

    	for nthSys in range(0,nSys):
    		recoHists.CreateHistMap(ptOrMassEnum.PT.value,   sysName+"_"+str(nthSys))
    		recoHists.CreateHistMap(ptOrMassEnum.MASS.value, sysName+"_"+str(nthSys))


    for filepath in sample.path:	
        
        infile = rt.TFile(filepath,'r')
        print (filepath)
        recoHists.saveRecoHists(infile, outfile, channel)
        del infile

    outfile.Write()
    outfile.Delete()

    del recoHists
    del vrecoHists
    return outDic

def makeSigHists(sample, outputDirectory, channel):
    rt.gROOT.SetBatch()
    print ("####################### makeSigHists #################################")
    
    outDic = {}
    
    outfile = rt.TFile(outputDirectory + sample.name+".root",'recreate')
    outfile_ = None
    
    if not sample.isAlt:  outDic[sample.name] = fHistDef.inputfHists(sample.name, outputDirectory + sample.name+".root", "sig")
    if sample.isAlt:      outDic[sample.name] = fHistDef.inputfHists(sample.name, outputDirectory + sample.name+".root", "sigAlt")
    # need to create histogram before the file loop to save one histogram from the multile input files
    sigHists  = rt.histTUnfold(rt.std.map("TString, TH1*")(), rt.std.map("TString, TH2*")(), sample.isInc, sample.isAlt)
    
    sigHists.SetPtBinningRec(channel)
    sigHists.SetMassBinningRec(channel)
    sigHists.SetPtBinningGen(channel)
    sigHists.SetMassBinningGen(channel)
    
    #sigHists.SetPtBinningRec("electron")
    #sigHists.SetMassBinningRec("electron")
    #sigHists.SetPtBinningGen("electron")
    #sigHists.SetMassBinningGen("electron")
    
    sigHists.CreateHist2DMap(ptOrMassEnum.PT.value, "nominal")
    sigHists.CreateHist2DMap(ptOrMassEnum.MASS.value, "nominal")
    
    sigHists.CreateHistMap(ptOrMassEnum.PT.value, "nominal")
    sigHists.CreateHistMap(ptOrMassEnum.MASS.value, "nominal")
    
    sigHists.CreateHistMap(ptOrMassEnum.PT.value, "nominal", "", False)
    sigHists.CreateHistMap(ptOrMassEnum.MASS.value, "nominal", "", False)
    
    if not sample.isAlt:
    	sysDict = {"PU": 2, "trgSF": 2, "recoSF": 2, "IdSF": 2, "IsoSF": 2, "unfoldsys": 1, "AlphaS": 2, "Scale": 9, "PDFerror": 100, "L1Prefire": 2, "LepScale": 2, "LepRes": 2, "FSRDR": 30}
    
    	for sysName, nSys in sysDict.items():
    	    sigHists.SetsysMap(sysName, nSys);
    
    	    for nthSys in range(0,nSys):
    	        sigHists.CreateHist2DMap(ptOrMassEnum.PT.value,   sysName+"_"+str(nthSys))
    	        sigHists.CreateHist2DMap(ptOrMassEnum.MASS.value, sysName+"_"+str(nthSys))
    
    	        sigHists.CreateHistMap(ptOrMassEnum.PT.value,     sysName+"_"+str(nthSys))
    	        sigHists.CreateHistMap(ptOrMassEnum.MASS.value,   sysName+"_"+str(nthSys))
    
    
    	if sample.isInc: # for DY to tautau, make one more histogram
    	    outfile_ = rt.TFile(outputDirectory + sample.name+"tau.root",'recreate')
    	    outDic[sample.name+"tau"] = fHistDef.inputfHists(sample.name+"tau", outputDirectory + sample.name+"tau.root", "bkg")
    
    	    sigHists.SetPtBinningRec(channel)
    	    sigHists.SetMassBinningRec(channel)
    
    	    #sigHists.SetPtBinningRec("electron")
    	    #sigHists.SetMassBinningRec("electron")
    
    	    sigHists.CreateHistMap(ptOrMassEnum.PT.value,  "nominal", "_tau")
    	    sigHists.CreateHistMap(ptOrMassEnum.MASS.value, "nominal", "_tau")
    
    
    	    for sysName, nSys in sysDict.items():
    	        sigHists.SetsysMap(sysName, nSys);
    
    	        for nthSys in range(0,nSys):
    	            sigHists.CreateHistMap(ptOrMassEnum.PT.value,     sysName+"_"+str(nthSys), "_tau")
    	            sigHists.CreateHistMap(ptOrMassEnum.MASS.value,   sysName+"_"+str(nthSys), "_tau")
 
    for filepath in sample.path:	
        
        infile = rt.TFile(filepath,'r')
        print (filepath)
        temp_kfactor = 1.;
        print ("temp_kfactor: " + str(temp_kfactor))
        sigHists.saveSigHists(infile, outfile, outfile_, channel, temp_kfactor)

        del infile
    
    outfile.Write()
    outfile.Delete()
    
    if sample.isInc and not sample.isAlt:
    	outfile_.Write()
    	outfile_.Delete()
    
    del sigHists
    
    return outDic
def makeMigrationM():
	pass
