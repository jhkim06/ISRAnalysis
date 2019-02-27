import os
import sys
import ROOT as rt

rt.gSystem.Load("/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/rootScripts/libunfold.so")
rt.gROOT.ProcessLine(".L /home/jhkim/ISR2016/unfolding/TUnfoldISR2016/rootScripts/makeRecoPlots.C++")

import gc
gc.collect()

# input file paths, binning, output directory
def makeRecoHists(sample):

        rt.gROOT.SetBatch()
        print "####################### makeRecoHists #################################"

        # don't need to check isMC for reco histogram 
	outfile = rt.TFile(sample.name+".root",'recreate')
        outfile_ = None

        # need to create histogram before the file loop to save one histogram from the input files
        recoHists = rt.recoTH1info(rt.vector('TH1*')(), rt.vector('TString')(), sample.isSig) 
        recoHists_ = rt.recoTH1info(rt.vector('TH1*')(), rt.vector('TString')(), sample.isSig) # just for DY to tautau
        
        recoHists.hists.push_back(rt.histogram(sample.name))
        recoHists.sysNames.push_back("norminal")

        if sample.isMC and sample.isSig: # for DY to tautau, make one more histogram
	        outfile_ = rt.TFile(sample.name+"tau.root",'recreate')

        	recoHists_.hists.push_back(rt.histogram("DYtotautau"))
	        recoHists_.sysNames.push_back("norminal")

        for filepath in sample.path:	

		infile = rt.TFile(filepath,'r')
		print filepath
		rt.recoHists(infile, outfile, outfile_, recoHists, recoHists_)

                del infile


        outfile.Write()
        outfile.Delete()

        if outfile_ is not None:
        	outfile_.Write()
        	outfile_.Delete()

        del recoHists
        del recoHists_

def makeMigrationM():
	pass
