import os
import sys
import ROOT as rt

rt.gSystem.Load("/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/lib/libisrunfold.so")

import gc
gc.collect()

def basicRatio(outpdf, unfold_pt, unfold_mass, channel):

        print "############################ draw ###############################################"

	rt.drawRatio(outpdf, unfold_pt, unfold_mass, channel)	

def drawUnfoldedPt(outpdf, unfold_pt, sysName, sysSize):

	rt.drawUnfoldedPtDistWithSys(outpdf, unfold_pt, sysName, sysSize)

def drawCombinedISR(outpdf, unfold_pt2016, unfold_mass2016, unfold_pt2017, unfold_mass2017 ):

	rt.drawCombinedISR(outpdf, unfold_pt2016, unfold_mass2016, unfold_pt2017, unfold_mass2017 )

def drawEMuCombinedISR(outpdf, unfold_ptElectron, unfold_massElectron, unfold_ptMuon, unfold_massMuon, unfold_ptElectron_MG, unfold_massElectron_MG, unfold_ptMuon_MG, unfold_massMuon_MG ):

        rt.drawEMuCombinedISR(outpdf, unfold_ptElectron, unfold_massElectron, unfold_ptMuon, unfold_massMuon, unfold_ptElectron_MG, unfold_massElectron_MG, unfold_ptMuon_MG, unfold_massMuon_MG )

def isrFit(outpdf, unfold_pt, unfold_mass):

        print "############################ draw fit ###############################################"

        rt.drawISRfit(outpdf, unfold_pt, unfold_mass)

def drawSys(outpdf, unfold_pt, unfold_mass, channel):

	rt.drawSystematicISR(outpdf, unfold_pt, unfold_mass, channel)

def drawMCSys(outpdf, unfold_pt, unfold_mass, channel):

        rt.drawMCSystematicISR(outpdf, unfold_pt, unfold_mass, channel)


def drawTest(outpdf, unfold_pt):

        rt.drawTest(outpdf, unfold_pt)

def basicRatioMass(outpdf, unfold):

        print "############################ draw ###############################################"

        rt.drawMassRatio(outpdf, unfold)

def responseM(outpdf, unfold, channel, var):

	rt.responseM(outpdf, unfold, channel, var)


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

