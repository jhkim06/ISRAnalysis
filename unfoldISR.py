import argparse
import os
import sys
import ROOT


parser = argparse.ArgumentParser(description='Unfolding for ISR analysis')
parser.add_argument('--createInputHists'  , action='store_true'  , help = 'create input histograms')
parser.add_argument('--setResMatrix'  , action='store_true'  , help = 'set response matrix')

args = parser.parse_args()

# get input root files information using sampleDef.py

import etc.sampleDef as isrSamples
import pyScripts.unfoldInputUtil as histUtil

outputDirectory = 'output/' 
inputfhisttxtName = outputDirectory + "fhist.txt"
if not os.path.exists( outputDirectory ):
	os.makedirs( outputDirectory )

# TODO allow to use several binning definitions

if args.createInputHists:
	import etc.histDef as inputfHist
	
	inputfHistDic = {}
	
	for sampleType in isrSamples.samplesDef_electron.keys():
	    sample =  isrSamples.samplesDef_electron[sampleType]
	    if sample is None : continue
	    print 'creating histogram for sample '
	    sample.dump()
	
	    if not sample.isMC:
	    	inputfHistDic.update(histUtil.makeRecoHists(sample, outputDirectory)) # TODO systematic list to consider
	
	    if sample.isMC and not sample.isSig:
	    	inputfHistDic.update(histUtil.makeRecoHists(sample, outputDirectory)) # TODO systematic list to consider
	
	    if sample.isMC and sample.isSig:
		inputfHistDic.update(histUtil.makeSigHists(sample, outputDirectory))
	
	fOutTxt = open( inputfhisttxtName,'w')
	
	for fhist in inputfHistDic.keys():
	        astr = inputfHistDic[fhist].name + ' ' + inputfHistDic[fhist].htype + ' ' + (inputfHistDic[fhist].path)[0]
	        fOutTxt.write( astr + '\n' )
	 	#print fhist
		#print inputfHistDic[fhist].name + ' ' + inputfHistDic[fhist].htype 
	        #print inputfHistDic[fhist].path
	
	fOutTxt.close()

# read text file
fOutTxtCheck = open( inputfhisttxtName,'r')
unfoldInputList = {}

for line in fOutTxtCheck:
	modifiedLine = line.lstrip(' ').rstrip(' ').rstrip('\n')
        if modifiedLine.split(' ')[1] == "data":
		unfoldInputList['data'] = modifiedLine.split(' ')[2]
        if modifiedLine.split(' ')[1] == "sig":
		unfoldInputList['sig'] = modifiedLine.split(' ')[2]
        if modifiedLine.split(' ')[1] == "bkg": # use the sample name as keyword for background
		unfoldInputList[modifiedLine.split(' ')[0]] = modifiedLine.split(' ')[2]

print unfoldInputList

import pyScripts.unfoldUtil as unfoldutil
import pyScripts.drawUtil as drawutil

if args.setResMatrix:

        # set unfolding class 
	#unfold_pt = unfoldutil.setUnfold(unfoldInputList['sig'], "Pt", "norminal")
	#unfold_pt = unfoldutil.setUnfold(unfoldInputList['sig'], "Pt", "fiducialPreFSR")
	unfold_pt = unfoldutil.setUnfold(unfoldInputList['sig'], "Pt", "ZptReweight")

	#unfold_mass = unfoldutil.setUnfold(unfoldInputList['sig'], "Mass", "norminal")
	#unfold_mass = unfoldutil.setUnfold(unfoldInputList['sig'], "Mass", "fiducialPreFSR")
	unfold_mass = unfoldutil.setUnfold(unfoldInputList['sig'], "Mass", "ZptReweight")

	# print out response matrix and efficiency plot

        outpdf = outputDirectory + "response.png"
        drawutil.responseM(outpdf, unfold_pt)
        outpdf = outputDirectory + "efficiency.png"
        drawutil.efficiency(outpdf, unfold_pt)

        outpdf_mass = outputDirectory + "response_mass.png"
        drawutil.responseM(outpdf_mass, unfold_mass)
        outpdf_mass = outputDirectory + "efficiency_mass.png"
        drawutil.efficiency(outpdf_mass, unfold_mass)

	# get bkg subtracted data
        unfoldutil.setUnfoldInput(unfold_pt, "Pt", unfoldInputList['data'])

        unfoldutil.subtractBkgs(unfold_pt, "Pt", unfoldInputList['DYtoEEtau'], "DYtoEEtau")
        unfoldutil.subtractBkgs(unfold_pt, "Pt", unfoldInputList['TTbar'], "TTbar")
        unfoldutil.subtractBkgs(unfold_pt, "Pt", unfoldInputList['VV'], "VV")
        unfoldutil.subtractBkgs(unfold_pt, "Pt", unfoldInputList['Wjets'], "Wjets")

	# for mass unfolding
        unfoldutil.setUnfoldInput(unfold_mass, "Mass", unfoldInputList['data'])

        unfoldutil.subtractBkgs(unfold_mass, "Mass", unfoldInputList['DYtoEEtau'], "DYtoEEtau")
        unfoldutil.subtractBkgs(unfold_mass, "Mass", unfoldInputList['TTbar'], "TTbar")
        unfoldutil.subtractBkgs(unfold_mass, "Mass", unfoldInputList['VV'], "VV")
        unfoldutil.subtractBkgs(unfold_mass, "Mass", unfoldInputList['Wjets'], "Wjets")

	# do unfolding with bkg subtracted data and the migration matrix
        unfoldutil.doUnfold(unfold_pt)
        unfoldutil.doUnfold(unfold_mass)

        outpdf = outputDirectory + "ratio.png"
        outpdf_mass = outputDirectory + "ratio_mass.png"
        # check unfolded distribution
        drawutil.basicRatio(outpdf, unfold_pt, unfold_mass, unfoldInputList['sig'])
        drawutil.basicRatioMass(outpdf_mass, unfold_mass, unfoldInputList['sig'])

        del unfold_pt
        del unfold_mass

def makeRecoPlots():
        # load TUnfold library 
	#ROOT.gSystem.Load("/home/jhkim/ISR2016/unfolding/TUnfold/libunfold.so")	
	#ROOT.gROOT.ProcessLine(".L /home/jhkim/ISR2016/unfolding/TUnfoldISR2016/rootScripts/makeRecoPlots.C++")
        pass

# test 
if __name__ == "__main__": 
	makeRecoPlots()
