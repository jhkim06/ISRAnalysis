import os
import sys
import ROOT as rt

rt.gSystem.Load("/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/rootScripts/libunfold.so")
rt.gROOT.ProcessLine(".L /home/jhkim/ISR2016/unfolding/TUnfoldISR2016/rootScripts/makeRecoPlots.C++")

# input file paths, binning, output directory
def makeRecoHists(sample):
        print "####################### makeRecoHists #################################"

        # don't need to check isMC for reco histogram 
	outfile = rt.TFile(sample.name+".root",'recreate')

        # need to create histogram before the file loop to save one histogram from the input files
        recoHists = rt.recoTH1info(rt.vector('TH1*')(), rt.vector('TString')()) 
        recoHists.hists.push_back(rt.histogram())
        recoHists.sysNames.push_back("norminal")

        for filepath in sample.path:	

		infile = rt.TFile(filepath,'r')
		print filepath
		rt.recoHists(infile, outfile, recoHists)

        outfile.Write()
        outfile.Delete()

def makeMigrationM():
	pass
