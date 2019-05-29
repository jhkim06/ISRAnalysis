import os
import sys
import ROOT as rt

rt.gSystem.Load("/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/lib/libisrunfold.so")

import etc.histDef as fHistDef

import gc
gc.collect()

from enum import Enum

class ptOrMassEnum(Enum):
	PT   = 1
	MASS = 2

# input file paths, binning, output directory, for data and bkg reco histograms only
def makeRecoHists(sample, outputDirectory, channel):


        rt.gROOT.SetBatch()
        print "####################### makeRecoHists #################################"

        outDic = {}

	outfile = rt.TFile(outputDirectory + sample.name+".root",'recreate')
        if not sample.isMC: outDic[sample.name] = fHistDef.inputfHists(sample.name, outputDirectory + sample.name+".root", "data")
        else : outDic[sample.name] = fHistDef.inputfHists(sample.name, outputDirectory + sample.name+".root", "bkg")

        # need to create histogram before the file loop to save one histogram from the mutiple input files
        recoHists = rt.histTUnfold(rt.std.map("TString, TH1*")())

        test = rt.histTUnfold(rt.std.map("TString, TH1*")())
	# https://root.cern.ch/faq/how-generate-dictionary
	vrecoHists = rt.std.vector("histTUnfold")()
	vrecoHists.push_back(test)
	vrecoHists.at(0).SetPtBinningRec()
	vrecoHists.at(0).SetMassBinningRec()
	vrecoHists.at(0).CreateHistMap(ptOrMassEnum.PT.value, "test") # 1: pt
	vrecoHists.at(0).CreateHistMap(ptOrMassEnum.MASS.value, "test") # 2: mass

	recoHists.SetPtBinningRec()
	recoHists.SetMassBinningRec()

	recoHists.CreateHistMap(ptOrMassEnum.PT.value, "nominal")
	recoHists.CreateHistMap(ptOrMassEnum.MASS.value, "nominal")

        for filepath in sample.path:	

		infile = rt.TFile(filepath,'r')
		print filepath
		recoHists.saveRecoHists(infile, outfile, channel)

                del infile

        outfile.Write()
        outfile.Delete()

        del recoHists
	del vrecoHists
        return outDic

def makeSigHists(sample, outputDirectory, channel):
        rt.gROOT.SetBatch()
        print "####################### makeSigHists #################################"

        outDic = {}

	outfile = rt.TFile(outputDirectory + sample.name+".root",'recreate')
        outfile_ = None

        outDic[sample.name] = fHistDef.inputfHists(sample.name, outputDirectory + sample.name+".root", "sig")
        # need to create histogram before the file loop to save one histogram from the multile input files
        sigHists  = rt.histTUnfold(rt.std.map("TString, TH1*")(), rt.std.map("TString, TH2*")(), sample.isInc) 

        sigHists.SetPtBinningRec()
        sigHists.SetMassBinningRec()
        sigHists.SetPtBinningGen()
        sigHists.SetMassBinningGen()
      
	sigHists.CreateHist2DMap(ptOrMassEnum.PT.value, "nominal") 
	sigHists.CreateHist2DMap(ptOrMassEnum.MASS.value, "nominal") 

	sigHists.CreateHistMap(ptOrMassEnum.PT.value, "nominal") 
	sigHists.CreateHistMap(ptOrMassEnum.MASS.value, "nominal") 
	
	# PU
	for i in range(2):
		sigHists.CreateHist2DMap(ptOrMassEnum.PT.value, "PU_"+str(i)) 
		sigHists.CreateHist2DMap(ptOrMassEnum.MASS.value, "PU_"+str(i)) 

		sigHists.CreateHist2DMap(ptOrMassEnum.PT.value, "trgSF_"+str(i)) 
		sigHists.CreateHist2DMap(ptOrMassEnum.MASS.value, "trgSF_"+str(i)) 

		sigHists.CreateHist2DMap(ptOrMassEnum.PT.value, "recoSF_"+str(i)) 
		sigHists.CreateHist2DMap(ptOrMassEnum.MASS.value, "recoSF_"+str(i)) 

		sigHists.CreateHist2DMap(ptOrMassEnum.PT.value, "IdSF_"+str(i)) 
		sigHists.CreateHist2DMap(ptOrMassEnum.MASS.value, "IdSF_"+str(i)) 

		sigHists.CreateHist2DMap(ptOrMassEnum.PT.value, "IsoSF_"+str(i)) 
		sigHists.CreateHist2DMap(ptOrMassEnum.MASS.value, "IsoSF_"+str(i)) 

	# unfolding systematic using Z pt correction (need to be checked)
	sigHists.CreateHist2DMap(ptOrMassEnum.PT.value, "unfoldsys_0") 
	sigHists.CreateHist2DMap(ptOrMassEnum.MASS.value, "unfoldsys_0") 

        for i in range(19):
		if i < 9:
			sigHists.CreateHist2DMap(ptOrMassEnum.PT.value, "FSRDR0p"+str(i+1)) 
			sigHists.CreateHist2DMap(ptOrMassEnum.MASS.value, "FSRDR0p"+str(i+1)) 
		else:
			sigHists.CreateHist2DMap(ptOrMassEnum.PT.value, "FSRDR1p"+str((i+1)%10)) 
			sigHists.CreateHist2DMap(ptOrMassEnum.MASS.value, "FSRDR1p"+str((i+1)%10)) 


	# AlphaS systematic
	for i in range(2):
		sigHists.CreateHist2DMap(ptOrMassEnum.PT.value, "AlphaS_"+str(i)) 
		sigHists.CreateHist2DMap(ptOrMassEnum.MASS.value, "AlphaS_"+str(i)) 

	# Scale systematic
	for i in range(9):
		sigHists.CreateHist2DMap(ptOrMassEnum.PT.value, "Scale_"+str(i)) 
		sigHists.CreateHist2DMap(ptOrMassEnum.MASS.value, "Scale_"+str(i)) 

	# PDF error
        for i in range(100):
                sigHists.CreateHist2DMap(ptOrMassEnum.PT.value, "PDFerror_"+str(i))
                sigHists.CreateHist2DMap(ptOrMassEnum.MASS.value, "PDFerror_"+str(i))


        if sample.isInc: # for DY to tautau, make one more histogram
	        outfile_ = rt.TFile(outputDirectory + sample.name+"tau.root",'recreate')
        	outDic[sample.name+"tau"] = fHistDef.inputfHists(sample.name+"tau", outputDirectory + sample.name+"tau.root", "bkg")

        	sigHists.SetPtBinningRec()
        	sigHists.SetMassBinningRec()

	        sigHists.CreateHistMap(ptOrMassEnum.PT.value, "nominal", "_tau") 
	        sigHists.CreateHistMap(ptOrMassEnum.MASS.value, "nominal", "_tau") 

        for filepath in sample.path:	

		infile = rt.TFile(filepath,'r')
		print filepath
                temp_kfactor = 1.;
		#if filepath.find("M-10to50") > 0:
		#   print "this is M-10to50 DY sample"
                #   temp_kfactor = 6225.42/5765.4
		print "temp_kfactor: " + str(temp_kfactor)
		sigHists.saveSigHists(infile, outfile, outfile_, channel, temp_kfactor)

                del infile

        outfile.Write()
        outfile.Delete()

        if sample.isInc:
        	outfile_.Write()
        	outfile_.Delete()

        del sigHists

	return outDic
def makeMigrationM():
	pass
