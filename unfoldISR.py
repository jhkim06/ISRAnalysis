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

# use histDef to save info about input histograms 
#TODO save the histogram info as txt

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
	print unfoldInputList['sig'] 
	print unfoldInputList['data'] 
	unfold = unfoldutil.setUnfold(unfoldInputList['sig'])
	# get bkg subtracted data
        unfoldutil.setUnfoldInput(unfold, unfoldInputList['data'])
	# do unfolding with bkg subtracted data and the migration matrix
        unfoldutil.doUnfold(unfold)

        outpdf = outputDirectory + "ratio.png"
        # check unfolded distribution
        drawutil.basicRatio(outpdf, unfold, unfoldInputList['sig'])
        del unfold

def makeRecoPlots():
        # load TUnfold library 
	#ROOT.gSystem.Load("/home/jhkim/ISR2016/unfolding/TUnfold/libunfold.so")	
	#ROOT.gROOT.ProcessLine(".L /home/jhkim/ISR2016/unfolding/TUnfoldISR2016/rootScripts/makeRecoPlots.C++")
        pass

# test 
if __name__ == "__main__": 
	makeRecoPlots()
