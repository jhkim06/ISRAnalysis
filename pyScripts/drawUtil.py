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

def drawTest(outpdf, unfold_pt):

        rt.drawTest(outpdf, unfold_pt)

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

def recoPt(outpdf, postfix, fdatapath, fdysigpath, fdybkgpath, fttbarpath, fvvpath, fwjetspath):

	fdata = rt.TFile(fdatapath,'r')
	fdysig = rt.TFile(fdysigpath,'r')
	fdybkg = rt.TFile(fdybkgpath,'r')
	fttbar = rt.TFile(fttbarpath,'r')
	fvv = rt.TFile(fvvpath,'r')
	fwjets = rt.TFile(fwjetspath,'r')

        rt.drawPtReco(outpdf, postfix, fdata, fdysig, fdybkg, fttbar, fvv, fwjets)
        del fdata
        del fdysig
        del fdybkg
        del fttbar
        del fvv
        del fwjets

def recoMass(outpdf, postfix, fdatapath, fdysigpath, fdybkgpath, fttbarpath, fvvpath, fwjetspath):

        fdata = rt.TFile(fdatapath,'r')
        fdysig = rt.TFile(fdysigpath,'r')
        fdybkg = rt.TFile(fdybkgpath,'r')
        fttbar = rt.TFile(fttbarpath,'r')
        fvv = rt.TFile(fvvpath,'r')
        fwjets = rt.TFile(fwjetspath,'r')

        rt.drawMassReco(outpdf, postfix, fdata, fdysig, fdybkg, fttbar, fvv, fwjets)
        del fdata
        del fdysig
        del fdybkg
        del fttbar
        del fvv
        del fwjets

