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

        del infile
	return unfold

def setUnfoldInput(unfold, filepath):

	print "Set input to unfold"
        infile = rt.TFile(filepath,'r')

        unfold = rt.setInput(unfold, infile)

        del infile


def subtractBkgs(unfold, filepath, name):

        print "Subtract background from data"
        infile = rt.TFile(filepath,'r')

        unfold = rt.subBkgs(unfold, infile, name)

        del infile

def doUnfold(unfold):

	print "do unfold!"

	unfold = rt.doUnfold(unfold)

