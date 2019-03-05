import os
import sys
import ROOT as rt

rt.gSystem.Load("/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/rootScripts/libunfold.so")
rt.gROOT.ProcessLine(".L /home/jhkim/ISR2016/unfolding/TUnfoldISR2016/rootScripts/drawUtils.C++")

import gc
gc.collect()

def basicRatio(outpdf, unfold_pt, unfold_mass, filepath):

        print "############################ draw ###############################################"
        print "Open " + filepath
        infile = rt.TFile(filepath,'r')

	rt.drawRatio(outpdf, unfold_pt, unfold_mass, infile)	
        del infile

def basicRatioMass(outpdf, unfold, filepath):

        print "############################ draw ###############################################"
        print "Open " + filepath
        infile = rt.TFile(filepath,'r')

        rt.drawMassRatio(outpdf, unfold, infile)
        del infile

def responseM(outpdf, unfold):

	rt.responseM(outpdf, unfold)


def efficiency(outpdf, unfold):

        rt.efficiency(outpdf, unfold)

