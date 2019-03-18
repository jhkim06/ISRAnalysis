import argparse
import os
import sys
import ROOT


parser = argparse.ArgumentParser(description='Unfolding for ISR analysis')

parser.add_argument('--channel' , dest = 'channel', default = 'electron', help = 'select channel electron or muon')
parser.add_argument('--createInputHists'  , action='store_true'  , help = 'create input histograms')
parser.add_argument('--createMatrixOnly'  , action='store_true'  , default = False, help = 'create histograms only for signal sample')
parser.add_argument('--setResMatrix'  , action='store_true'  , help = 'set response matrix')

args = parser.parse_args()

# get input root files information using sampleDef.py
import etc.sampleDef as isrSamples
import pyScripts.unfoldInputUtil as histUtil

outputDirectory = 'output/' + args.channel + "/"  
inputfhisttxtName = outputDirectory + "fhist.txt"
if not os.path.exists( outputDirectory ):
	os.makedirs( outputDirectory )

print "channel to run: " + args.channel

selectedSample = {} 
if args.channel == "electron":
	selectedSample = isrSamples.samplesDef_electronLegacy
if args.channel == "muon":
	selectedSample = isrSamples.samplesDef_muonLegacy 

if args.createInputHists:
	import etc.histDef as inputfHist
	
	inputfHistDic = {}
	
	for sampleType in selectedSample.keys():
	    sample =  selectedSample[sampleType]
	    if sample is None : continue
	
 	    if not args.createMatrixOnly:
	    	 print 'creating histogram for sample '
	    	 sample.dump()
	   	 if not sample.isMC:
	    		inputfHistDic.update(histUtil.makeRecoHists(sample, outputDirectory, args.channel)) # TODO systematic list to consider
	
	    	 if sample.isMC and not sample.isSig:
	    		inputfHistDic.update(histUtil.makeRecoHists(sample, outputDirectory, args.channel)) # TODO systematic list to consider
	
	    if sample.isMC and sample.isSig:
	    	print 'creating histogram for sample '
	    	sample.dump()
		inputfHistDic.update(histUtil.makeSigHists(sample, outputDirectory, args.channel))
	
	fOutTxt = open( inputfhisttxtName,'w')
	
	for fhist in inputfHistDic.keys():
	        astr = inputfHistDic[fhist].name + ' ' + inputfHistDic[fhist].htype + ' ' + (inputfHistDic[fhist].path)[0]
	        fOutTxt.write( astr + '\n' )
	
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

        postfix = "norminal"

        # set unfolding class 
	unfold_pt = unfoldutil.setUnfold(unfoldInputList['sig'], "Pt", postfix)

	unfold_mass = unfoldutil.setUnfold(unfoldInputList['sig'], "Mass", postfix)

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
        unfoldutil.setUnfoldInput(unfold_pt, "Pt", postfix, unfoldInputList['data'])
        #unfoldutil.setUnfoldInput(unfold_pt, "Pt", unfoldInputList['sig'])

        unfoldutil.subtractBkgs(unfold_pt, "Pt", postfix, unfoldInputList['DYtoEEtau'], "DYtoEEtau")
        unfoldutil.subtractBkgs(unfold_pt, "Pt", postfix, unfoldInputList['TTbar'], "TTbar")
        unfoldutil.subtractBkgs(unfold_pt, "Pt", postfix, unfoldInputList['VV'], "VV")
        unfoldutil.subtractBkgs(unfold_pt, "Pt", postfix, unfoldInputList['Wjets'], "Wjets")

	# for mass unfolding
        unfoldutil.setUnfoldInput(unfold_mass, "Mass", postfix, unfoldInputList['data'])
        #unfoldutil.setUnfoldInput(unfold_mass, "Mass", unfoldInputList['sig'])

        unfoldutil.subtractBkgs(unfold_mass, "Mass", postfix, unfoldInputList['DYtoEEtau'], "DYtoEEtau")
        unfoldutil.subtractBkgs(unfold_mass, "Mass", postfix, unfoldInputList['TTbar'], "TTbar")
        unfoldutil.subtractBkgs(unfold_mass, "Mass", postfix, unfoldInputList['VV'], "VV")
        unfoldutil.subtractBkgs(unfold_mass, "Mass", postfix, unfoldInputList['Wjets'], "Wjets")

	# do unfolding with bkg subtracted data and the migration matrix
        unfoldutil.doUnfold(unfold_pt)
        unfoldutil.doUnfold(unfold_mass)

        outpdf = outputDirectory + "ratio.png"
        outpdf_mass = outputDirectory + "ratio_mass.png"
        # check unfolded distribution
        drawutil.basicRatio(outpdf, unfold_pt, unfold_mass, unfoldInputList['sig'])
        drawutil.basicRatioMass(outpdf_mass, unfold_mass, unfoldInputList['sig'])

        outpdf = outputDirectory + "recoPt.png"
        drawutil.recoPt(outpdf, postfix, unfoldInputList['data'], unfoldInputList['sig'], unfoldInputList['DYtoEEtau'], unfoldInputList['TTbar'], unfoldInputList['VV'], unfoldInputList['Wjets']);

        outpdf = outputDirectory + "recoMass.png"
        drawutil.recoMass(outpdf, postfix, unfoldInputList['data'], unfoldInputList['sig'], unfoldInputList['DYtoEEtau'], unfoldInputList['TTbar'], unfoldInputList['VV'], unfoldInputList['Wjets']);

	outpdf = outputDirectory + "test.png" 
	drawutil.drawTest(outpdf, unfold_pt)

        del unfold_pt
        del unfold_mass

def makeRecoPlots():
        pass

# test 
if __name__ == "__main__": 
	makeRecoPlots()
