import os
import sys
import ROOT as rt

rt.gSystem.Load("/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/rootScripts/libunfold.so")
rt.gROOT.ProcessLine(".L /home/jhkim/ISR2016/unfolding/TUnfoldISR2016/rootScripts/drawUtils.C++")

import gc
gc.collect()



def basicRatio(outpdf, unfold, filepath):

        print "############################ draw ###############################################"
        print "Open " + filepath
        infile = rt.TFile(filepath,'r')

	rt.drawRatio(outpdf, unfold, infile)	
        del infile
