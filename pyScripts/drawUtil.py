import os
import sys
import ROOT as rt

rt.gSystem.Load("/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/lib/libisrunfold.so")

import gc
gc.collect()

def basicRatio(outpdf, unfold_pt, unfold_mass, filepath, channel):

        print "############################ draw ###############################################"
        print "Open " + filepath
        infile = rt.TFile(filepath,'r')

	rt.drawRatio(outpdf, unfold_pt, unfold_mass, infile, channel)	
        del infile

def drawUnfoldedPt(outpdf, unfold_pt):

	rt.drawUnfoldedPtDistWithSys(outpdf, unfold_pt)

def drawCombinedISR(outpdf, unfold_pt2016, unfold_mass2016, unfold_pt2017, unfold_mass2017 ):

	rt.drawCombinedISR(outpdf, unfold_pt2016, unfold_mass2016, unfold_pt2017, unfold_mass2017 )

def drawEMuCombinedISR(outpdf, unfold_ptElectron, unfold_massElectron, unfold_ptMuon, unfold_massMuon ):

        rt.drawEMuCombinedISR(outpdf, unfold_ptElectron, unfold_massElectron, unfold_ptMuon, unfold_massMuon )

def isrFit(outpdf, unfold_pt, unfold_mass, filepath):

        print "############################ draw fit ###############################################"
        print "Open " + filepath
        infile = rt.TFile(filepath,'r')

        rt.drawISRfit(outpdf, unfold_pt, unfold_mass, infile)
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

def recoPt(outpdf, postfix, fdatapath, fdysigpath, fdybkgpath, fttbarpath, fvvpath, fwjetspath, fqcdpath, channel):

	fdata = rt.TFile(fdatapath,'r')
	fdysig = rt.TFile(fdysigpath,'r')
	fdybkg = rt.TFile(fdybkgpath,'r')
	fttbar = rt.TFile(fttbarpath,'r')
	fvv = rt.TFile(fvvpath,'r')
	fwjets = rt.TFile(fwjetspath,'r')
	if fqcdpath != None:
		fqcd = rt.TFile(fqcdpath,'r')
	else :
		fqcd = None
		

        rt.drawPtReco(outpdf, postfix, fdata, fdysig, fdybkg, fttbar, fvv, fwjets, fqcd, channel)
        del fdata
        del fdysig
        del fdybkg
        del fttbar
        del fvv
        del fwjets

def recoMass(outpdf, postfix, fdatapath, fdysigpath, fdybkgpath, fttbarpath, fvvpath, fwjetspath, fqcdpath, channel):

        fdata = rt.TFile(fdatapath,'r')
        fdysig = rt.TFile(fdysigpath,'r')
        fdybkg = rt.TFile(fdybkgpath,'r')
        fttbar = rt.TFile(fttbarpath,'r')
        fvv = rt.TFile(fvvpath,'r')
        fwjets = rt.TFile(fwjetspath,'r')
        if fqcdpath != None:
                fqcd = rt.TFile(fqcdpath,'r')
        else :
                fqcd = None


        rt.drawMassReco(outpdf, postfix, fdata, fdysig, fdybkg, fttbar, fvv, fwjets, fqcd, channel)
        del fdata
        del fdysig
        del fdybkg
        del fttbar
        del fvv
        del fwjets

