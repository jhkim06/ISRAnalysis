import os
import sys
import ROOT

# How to use TUnfold: try to load library (" gSystem->Load("./libunfold.so")
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

# create histograms from the root files
# try to use libunfold.so and make test histogram(response matrix) using DY MC
# make response matrix: from detector to post FSR, from post to pre FSR, from detector to pre FSR
# make directories for each migration matrix

def makeRecoPlots():
        # load TUnfold library 
	ROOT.gSystem.Load("/home/jhkim/ISR2016/unfolding/TUnfold/libunfold.so")	
	ROOT.gROOT.ProcessLine(".L /home/jhkim/ISR2016/unfolding/TUnfoldISR2016/rootScripts/makeRecoPlots.C++")

        ROOT.test()
# create migration matrix

# unfolding: inputs: histogram to unfold, migration matrix option to subtract background

# test 
if __name__ == "__main__": 
	makeRecoPlots()
