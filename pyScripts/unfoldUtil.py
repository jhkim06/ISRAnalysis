import os
import sys
import ROOT as rt

rt.gSystem.Load("/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/lib/libunfoldUtils_C.so")

import gc
gc.collect()

def setUnfold(filepath, var, matrixName, isfsr):

        print "Open " + filepath
	infile = rt.TFile(filepath,'r')
        
        unfold = rt.setTUnfoldDensity(infile, var, matrixName, isfsr) # set the migration matrix 
							       #TODO allow to select different unfolding option
        del infile
	return unfold

def setVectorSys(filepath, unfold, var, sysMatrixName, size):
	print "Set vertor type systematic" 

	infile = rt.TFile(filepath, 'r')

        rt.setVetorSystematic(infile, unfold, var, sysMatrixName, size)

def setUnfoldInput(unfold, var, postfix, filepath):

	print "Set input to unfold"
        infile = rt.TFile(filepath,'r')

        unfold = rt.setInput(unfold, var, postfix, infile)

        del infile

def setUnfoldInputHist(unfold, inputhist):

	print "Set unfold input from histogram"
	rt.setInputHist(unfold, inputhist);


def subtractBkgs(unfold, var, hname, filepath, name):

        print "Subtract background from data"
        infile = rt.TFile(filepath,'r')

        unfold = rt.subBkgs(unfold, var, hname, infile, name)

        del infile

def doUnfold(unfold):

	print "do unfold!"

	unfold = rt.doUnfold(unfold)

def getUnfoldedHist(unfold, postfix):
	print "Get unfolded histogram"

	return rt.getHist(unfold, postfix)

