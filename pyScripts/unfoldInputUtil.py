import os
import sys
import ROOT as rt

rt.gSystem.Load("/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/rootScripts/libunfold.so")
rt.gROOT.ProcessLine(".L /home/jhkim/ISR2016/unfolding/TUnfoldISR2016/rootScripts/makeRecoPlots.C++")

import gc
gc.collect()

# input file paths, binning, output directory, for data and bkg reco histograms only
def makeRecoHists(sample):

        rt.gROOT.SetBatch()
        print "####################### makeRecoHists #################################"

        # don't need to check isMC for reco histogram 
	outfile = rt.TFile(sample.name+".root",'recreate')

        # need to create histogram before the file loop to save one histogram from the input files
        recoHists = rt.recoHistsinfo(rt.vector('TH1*')(), rt.vector('TString')()) 
        
        recoHists.hists.push_back(rt.histogram("norminal"))
        recoHists.sysNames.push_back("norminal")

        for filepath in sample.path:	

		infile = rt.TFile(filepath,'r')
		print filepath
		rt.recoHists(infile, outfile, recoHists)

                del infile


        outfile.Write()
        outfile.Delete()

        del recoHists

def makeSigHists(sample):
        rt.gROOT.SetBatch()
        print "####################### makeRecoHists #################################"

        # don't need to check isMC for reco histogram 
	outfile = rt.TFile(sample.name+".root",'recreate')
        outfile_ = None

        # need to create histogram before the file loop to save one histogram from the input files
        sigHists  = rt.sigHistsinfo(rt.vector('TH1*')(), rt.vector('TH2*')(), rt.vector('TString')(), sample.isInc) 
        recoHists = rt.recoHistsinfo(rt.vector('TH1*')(), rt.vector('TString')())

        
        sigHists.hists.push_back(rt.histogram("norminal"))
        sigHists.matrixs.push_back(rt.matrix("norminal"))
        sigHists.sysNames.push_back("norminal")

        if sample.isInc: # for DY to tautau, make one more histogram
	        outfile_ = rt.TFile(sample.name+"tau.root",'recreate')

        	recoHists.hists.push_back(rt.histogram("norminal"))
	        recoHists.sysNames.push_back("norminal")

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

def makeMigrationM():
	pass
