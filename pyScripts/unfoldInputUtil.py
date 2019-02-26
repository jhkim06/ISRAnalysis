import os
import sys
import ROOT as rt

rt.gSystem.Load("/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/rootScripts/libunfold.so")
rt.gROOT.ProcessLine(".L /home/jhkim/ISR2016/unfolding/TUnfoldISR2016/rootScripts/makeRecoPlots.C++")


#input file paths, binning, output directory
def makeRecoHists(sample):
        print "####################### makeRecoHists #################################"

        for filepath in sample.path:
        	if not sample.isMC:
			infile = rt.TFile(filepath,'r')
			print filepath
                        rt.recoHists(infile)
                        # call function to save 
        	if sample.isMC:
			#infile = rt.TFile(sample.path,'r')
			print filepath

def makeMigrationM():
	pass
