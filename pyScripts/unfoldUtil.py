import os
import sys
import ROOT as rt

rt.gSystem.Load("/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/rootScripts/libunfold.so")
rt.gROOT.ProcessLine(".L /home/jhkim/ISR2016/unfolding/TUnfoldISR2016/rootScripts/unfoldUtils.C++")

import gc
gc.collect()


def setUnfold(filepath):

        print "Open " + filepath
	infile = rt.TFile(filepath,'r')
        
        unfold = rt.setTUnfoldDensity(infile) # set the migration matrix TODO allow to select different unfolding option

	return unfold
