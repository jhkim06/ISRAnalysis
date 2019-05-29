import os
import sys
import ROOT as rt

rt.gSystem.Load("/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/lib/libisrunfold.so")

import etc.histDef as fHistDef

import gc
gc.collect()

# input file paths, binning, output directory, for data and bkg reco histograms only
def makeRecoHists(sample, outputDirectory, channel):

        rt.gROOT.SetBatch()
        print "####################### makeRecoHists #################################"

        outDic = {}

        # don't need to check isMC for reco histogram 
	outfile = rt.TFile(outputDirectory + sample.name+".root",'recreate')
        if not sample.isMC: outDic[sample.name] = fHistDef.inputfHists(sample.name, outputDirectory + sample.name+".root", "data")
        else : outDic[sample.name] = fHistDef.inputfHists(sample.name, outputDirectory + sample.name+".root", "bkg")

        # need to create histogram before the file loop to save one histogram from the input files
        recoHists = rt.histTUnfold(rt.std.map("TString, TH1*")())

        test = rt.histTUnfold(rt.std.map("TString, TH1*")())
	# https://root.cern.ch/faq/how-generate-dictionary
	vrecoHists = rt.std.vector("histTUnfold")()
	vrecoHists.push_back(test)
	vrecoHists.at(0).SetPtBinningRec()
	vrecoHists.at(0).SetMassBinningRec()
	vrecoHists.at(0).CreateHistMap(1, "test") # 1: pt
	vrecoHists.at(0).CreateHistMap(2, "test") # 2: mass

	recoHists.SetPtBinningRec()
	recoHists.SetMassBinningRec()

	recoHists.CreateHistMap(1, "nominal")
	recoHists.CreateHistMap(2, "nominal")

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

        # don't need to check isMC for reco histogram 
	outfile = rt.TFile(outputDirectory + sample.name+".root",'recreate')
        outfile_ = None

        outDic[sample.name] = fHistDef.inputfHists(sample.name, outputDirectory + sample.name+".root", "sig")
        # need to create histogram before the file loop to save one histogram from the input files
        sigHists  = rt.histTUnfold(rt.std.map("TString, TH1*")(), rt.std.map("TString, TH2*")(), sample.isInc) 

        sigHists.SetPtBinningRec()
        sigHists.SetMassBinningRec()
        sigHists.SetPtBinningGen()
        sigHists.SetMassBinningGen()
      
	sigHists.CreateHist2DMap(1, "nominal") 
	sigHists.CreateHist2DMap(2, "nominal") 
	
	# PU
	for i in range(2):
		sigHists.CreateHist2DMap(1, "PU_"+str(i)) 
		sigHists.CreateHist2DMap(2, "PU_"+str(i)) 

		sigHists.CreateHist2DMap(1, "trgSF_"+str(i)) 
		sigHists.CreateHist2DMap(2, "trgSF_"+str(i)) 

		sigHists.CreateHist2DMap(1, "recoSF_"+str(i)) 
		sigHists.CreateHist2DMap(2, "recoSF_"+str(i)) 

		sigHists.CreateHist2DMap(1, "IdSF_"+str(i)) 
		sigHists.CreateHist2DMap(2, "IdSF_"+str(i)) 

		sigHists.CreateHist2DMap(1, "IsoSF_"+str(i)) 
		sigHists.CreateHist2DMap(2, "IsoSF_"+str(i)) 

	# unfolding systematic using Z pt correction (need to be checked)
	sigHists.CreateHist2DMap(1, "unfoldsys_0") 
	sigHists.CreateHist2DMap(2, "unfoldsys_0") 

        #for i in range(19):
	#	if i < 9:
	#		sigHists.CreateHist2DMap(1, "FSRDR0p"+str(i+1)) 
	#		sigHists.CreateHist2DMap(2, "FSRDR0p"+str(i+1)) 
	#	else:
	#		sigHists.CreateHist2DMap(1, "FSRDR1p"+str((i+1)%10)) 
	#		sigHists.CreateHist2DMap(2, "FSRDR1p"+str((i+1)%10)) 

	sigHists.CreateHistMap(1, "nominal") 
	sigHists.CreateHistMap(2, "nominal") 

	# AlphaS systematic
	for i in range(2):
		sigHists.CreateHist2DMap(1, "AlphaS_"+str(i)) 
		sigHists.CreateHist2DMap(2, "AlphaS_"+str(i)) 

	# Scale systematic
	for i in range(9):
		sigHists.CreateHist2DMap(1, "Scale_"+str(i)) 
		sigHists.CreateHist2DMap(2, "Scale_"+str(i)) 

        #for i in range(100):
        #        sigHists.CreateHist2DMap(1, "PDFerror_"+str(i))
        #        sigHists.CreateHist2DMap(2, "PDFerror_"+str(i))


        if sample.isInc: # for DY to tautau, make one more histogram
	        outfile_ = rt.TFile(outputDirectory + sample.name+"tau.root",'recreate')
        	outDic[sample.name+"tau"] = fHistDef.inputfHists(sample.name+"tau", outputDirectory + sample.name+"tau.root", "bkg")

        	sigHists.SetPtBinningRec()
        	sigHists.SetMassBinningRec()

	        sigHists.CreateHistMap(1, "nominal", "_tau") 
	        sigHists.CreateHistMap(2, "nominal", "_tau") 

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
