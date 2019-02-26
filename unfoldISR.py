import os
import sys
import ROOT

# Class for unfolding
# Reco histograms (unfolding inputs): one histogram for data, one histogram for signal MC, dynamic size for background MCs
# Migration matrix

#class UnfoldingData:

#class UnfoldingBackground:

#class UnfoldingSigal:
#	def __init__(self, isSystematic):
#		self.systematic = isSystematic
#                 

# get input root files information using sampleDef.py

import etc.sampleDef as isrSamples
import pyScripts.unfoldInputUtil as histUtil

outputDirectory = 'output/' 

if not os.path.exists( outputDirectory ):
	os.makedirs( outputDirectory )

for sampleType in isrSamples.samplesDef_electron.keys():
    sample =  isrSamples.samplesDef_electron[sampleType]
    if sample is None : continue
    #if sampleType == args.sample or args.sample == 'all' :
    print 'creating histogram for sample '
    sample.dump()

    if not sample.isMC: 
	print "data, so only reco histograms"
        histUtil.makeRecoHists(sample) # FIXME systematic list to consider
    if sample.isMC and sample.isSig: 
        print "make migration matrix also"
        histUtil.makeRecoHists(sample)
	

    # path list, check if it is data or sig mc, bkg mc, binning info
    #tnpHist.makePassFailHistograms( sample, tnpConf.flags[args.flag], tnpBins, var )

# create histograms from the root files
# make test histogram(response matrix) using DY MC
# make response matrix: from detector to post FSR, from post to pre FSR, from detector to pre FSR
# make directories for each migration matrix

def makeRecoPlots():
        # load TUnfold library 
	ROOT.gSystem.Load("/home/jhkim/ISR2016/unfolding/TUnfold/libunfold.so")	
	#ROOT.gROOT.ProcessLine(".L /home/jhkim/ISR2016/unfolding/TUnfoldISR2016/rootScripts/makeRecoPlots.C++")

        #ROOT.dataHist()
# create migration matrix

# unfolding: inputs: histogram to unfold, migration matrix option to subtract background

# test 
if __name__ == "__main__": 
	makeRecoPlots()
