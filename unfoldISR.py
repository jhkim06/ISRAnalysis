import os
import sys
import ROOT

# get input root files information using sampleDef.py

import etc.sampleDef as isrSamples
import pyScripts.unfoldInputUtil as histUtil

outputDirectory = 'output/' 

if not os.path.exists( outputDirectory ):
	os.makedirs( outputDirectory )

# use histDef to save info about input histograms 

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

for fhist in inputfHistDic.keys():
 	#print fhist
	print inputfHistDic[fhist].name + ' ' + inputfHistDic[fhist].htype 
        print inputfHistDic[fhist].path

# import pyScript.unfoldUtil as doUnfold
# get bkg subtracted data
# do unfolding with bkg subtracted data and the migration matrix

def makeRecoPlots():
        # load TUnfold library 
	ROOT.gSystem.Load("/home/jhkim/ISR2016/unfolding/TUnfold/libunfold.so")	
	#ROOT.gROOT.ProcessLine(".L /home/jhkim/ISR2016/unfolding/TUnfoldISR2016/rootScripts/makeRecoPlots.C++")

# test 
if __name__ == "__main__": 
	makeRecoPlots()
