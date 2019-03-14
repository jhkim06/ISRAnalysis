import os
import sys
import ROOT as rt

rt.gSystem.Load("/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/rootScripts/libunfold.so")
rt.gROOT.ProcessLine(".L /home/jhkim/ISR2016/unfolding/TUnfoldISR2016/rootScripts/unfoldUtils.C++")

import gc
gc.collect()


def setUnfold(filepath, var, matrixName):

        print "Open " + filepath
	infile = rt.TFile(filepath,'r')
        
        unfold = rt.setTUnfoldDensity(infile, var, matrixName) # set the migration matrix TODO allow to select different unfolding option

        del infile
	return unfold

def setUnfoldInput(unfold, var, postfix, filepath):

	print "Set input to unfold"
        infile = rt.TFile(filepath,'r')

        unfold = rt.setInput(unfold, var, postfix, infile)

        del infile


def subtractBkgs(unfold, var, hname, filepath, name):

        print "Subtract background from data"
        infile = rt.TFile(filepath,'r')

        unfold = rt.subBkgs(unfold, var, hname, infile, name)

        del infile

def doUnfold(unfold):

	print "do unfold!"

	unfold = rt.doUnfold(unfold)

