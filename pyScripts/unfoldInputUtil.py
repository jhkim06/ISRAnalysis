import os
import sys
import ROOT as rt

rt.gSystem.Load("/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/rootScripts/libunfold.so")
rt.gROOT.ProcessLine(".L /home/jhkim/ISR2016/unfolding/TUnfoldISR2016/rootScripts/makeRecoPlots.C++")

import gc
gc.collect()

import etc.histDef as fHistDef

# input file paths, binning, output directory, for data and bkg reco histograms only
def makeRecoHists(sample, outputDirectory):

        rt.gROOT.SetBatch()
        print "####################### makeRecoHists #################################"

        outDic = {}

        # don't need to check isMC for reco histogram 
	outfile = rt.TFile(outputDirectory + sample.name+".root",'recreate')
        if not sample.isMC: outDic[sample.name] = fHistDef.inputfHists(sample.name, outputDirectory + sample.name+".root", "data")
        else : outDic[sample.name] = fHistDef.inputfHists(sample.name, outputDirectory + sample.name+".root", "bkg")

        # need to create histogram before the file loop to save one histogram from the input files
        recoHists = rt.recoHistsinfo(rt.vector('TH1*')(), rt.vector('TH1*')(), rt.vector('TString')()) 
        
        recoHists.ptHists.push_back(rt.ptHistogram("norminal"))
        recoHists.massHists.push_back(rt.massHistogram("norminal"))
        recoHists.sysNames.push_back("norminal")

        recoHists.ptHists.push_back(rt.ptHistogram("fiducialPreFSR"))
        recoHists.massHists.push_back(rt.massHistogram("fiducialPreFSR"))
        recoHists.sysNames.push_back("fiducialPreFSR")

        recoHists.ptHists.push_back(rt.ptHistogram("ZptReweight"))
        recoHists.massHists.push_back(rt.massHistogram("ZptReweight"))
        recoHists.sysNames.push_back("ZptReweight")

        recoHists.ptHists.push_back(rt.ptHistogram("noBveto"))
        recoHists.massHists.push_back(rt.massHistogram("noBveto"))
        recoHists.sysNames.push_back("noBveto")

        for filepath in sample.path:	

		infile = rt.TFile(filepath,'r')
		print filepath
		rt.recoHists(infile, outfile, recoHists)

                del infile

        outfile.Write()
        outfile.Delete()

        del recoHists
        return outDic

def makeSigHists(sample, outputDirectory):
        rt.gROOT.SetBatch()
        print "####################### makeRecoHists #################################"

        outDic = {}

        # don't need to check isMC for reco histogram 
	outfile = rt.TFile(outputDirectory + sample.name+".root",'recreate')
        outfile_ = None

        outDic[sample.name] = fHistDef.inputfHists(sample.name, outputDirectory + sample.name+".root", "sig")
        # need to create histogram before the file loop to save one histogram from the input files
        sigHists  = rt.sigHistsinfo(rt.vector('TH1*')(), rt.vector('TH1*')(), rt.vector('TH2*')(), rt.vector('TH2*')(), rt.vector('TString')(), sample.isInc) 
        recoHists = rt.recoHistsinfo(rt.vector('TH1*')(), rt.vector('TH1*')(), rt.vector('TString')())
        
        sigHists.ptHists.push_back(rt.ptHistogram("norminal"))
        sigHists.ptMatrixs.push_back(rt.ptMatrix("norminal"))
        sigHists.massHists.push_back(rt.massHistogram("norminal"))
        sigHists.massMatrixs.push_back(rt.massMatrix("norminal"))
        sigHists.sysNames.push_back("norminal")

        sigHists.ptHists.push_back(rt.ptHistogram("fiducialPreFSR"))
        sigHists.ptMatrixs.push_back(rt.ptMatrix("fiducialPreFSR"))
        sigHists.massHists.push_back(rt.massHistogram("fiducialPreFSR"))
        sigHists.massMatrixs.push_back(rt.massMatrix("fiducialPreFSR"))
        sigHists.sysNames.push_back("fiducialPreFSR")

        sigHists.ptHists.push_back(rt.ptHistogram("ZptReweight"))
        sigHists.ptMatrixs.push_back(rt.ptMatrix("ZptReweight"))
        sigHists.massHists.push_back(rt.massHistogram("ZptReweight"))
        sigHists.massMatrixs.push_back(rt.massMatrix("ZptReweight"))
        sigHists.sysNames.push_back("ZptReweight")

        sigHists.ptHists.push_back(rt.ptHistogram("noBveto"))
        sigHists.ptMatrixs.push_back(rt.ptMatrix("noBveto"))
        sigHists.massHists.push_back(rt.massHistogram("noBveto"))
        sigHists.massMatrixs.push_back(rt.massMatrix("noBveto"))
        sigHists.sysNames.push_back("noBveto")

        if sample.isInc: # for DY to tautau, make one more histogram
	        outfile_ = rt.TFile(outputDirectory + sample.name+"tau.root",'recreate')
        	outDic[sample.name+"tau"] = fHistDef.inputfHists(sample.name+"tau", outputDirectory + sample.name+"tau.root", "bkg")

        	recoHists.ptHists.push_back(rt.ptHistogram("norminal"))
        	recoHists.massHists.push_back(rt.massHistogram("norminal"))
	        recoHists.sysNames.push_back("norminal")

                recoHists.ptHists.push_back(rt.ptHistogram("fiducialPreFSR"))
                recoHists.massHists.push_back(rt.massHistogram("fiducialPreFSR"))
                recoHists.sysNames.push_back("fiducialPreFSR")

                recoHists.ptHists.push_back(rt.ptHistogram("ZptReweight"))
                recoHists.massHists.push_back(rt.massHistogram("ZptReweight"))
                recoHists.sysNames.push_back("ZptReweight")

                recoHists.ptHists.push_back(rt.ptHistogram("noBveto"))
                recoHists.massHists.push_back(rt.massHistogram("noBveto"))
                recoHists.sysNames.push_back("noBveto")

        for filepath in sample.path:	

		infile = rt.TFile(filepath,'r')
		print filepath
		rt.sigHists(infile, outfile, outfile_, sigHists, recoHists)

                del infile

        outfile.Write()
        outfile.Delete()

        if sample.isInc:
        	outfile_.Write()
        	outfile_.Delete()

        del sigHists
        del recoHists

	return outDic
def makeMigrationM():
	pass
