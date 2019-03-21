import argparse
import os
import sys
import ROOT


parser = argparse.ArgumentParser(description='Unfolding for ISR analysis')

parser.add_argument('--channel' , dest = 'channel', default = 'electron', help = 'select channel electron or muon')
parser.add_argument('--postfix' , dest = 'postfix', default = 'norminal', help = 'select histogram name')
parser.add_argument('--createInputHists'  , action='store_true'  , help = 'create input histograms')
parser.add_argument('--createMatrixOnly'  , action='store_true'  , default = False, help = 'create histograms only for signal sample')
parser.add_argument('--getUnfoldResults'  , action='store_true'  , help = 'Get unfolding resutls')
parser.add_argument('--doSeperateUnfold'  , action='store_true'  , default = False, help = 'Seperate unfolding steps')

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

if args.getUnfoldResults:

	import pyScripts.unfoldUtil as unfoldutil
	import pyScripts.drawUtil as drawutil
	
        postfix = args.postfix
	postfix_matrix = args.postfix
	if args.doSeperateUnfold:
		postfix = "norminal"
		postfix_matrix = "detector"

        # set unfolding class 
	unfold_pt =   unfoldutil.setUnfold(unfoldInputList['sig'], "Pt",   postfix_matrix, False)
	unfold_mass = unfoldutil.setUnfold(unfoldInputList['sig'], "Mass", postfix_matrix, False)

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

        # set systematic sources
        # two types: 
        # 1) vetor type
        #       input: sys name, vector size 
        # 2) up/down type
        #       input: sys name

        unfoldutil.setVectorSys(unfoldInputList['sig'], unfold_pt,   "Pt",   "AlphaS", 2)
        unfoldutil.setVectorSys(unfoldInputList['sig'], unfold_mass, "Mass", "AlphaS", 2)

        unfoldutil.setVectorSys(unfoldInputList['sig'], unfold_pt,   "Pt",   "Scale", 9)
        unfoldutil.setVectorSys(unfoldInputList['sig'], unfold_mass, "Mass", "Scale", 9)

        outpdf = outputDirectory + "ratio.png"
        outpdf_mass = outputDirectory + "ratio_mass.png"
        # check unfolded distribution
        drawutil.basicRatio(outpdf, unfold_pt, unfold_mass, unfoldInputList['sig'])
        drawutil.basicRatioMass(outpdf_mass, unfold_mass, unfoldInputList['sig'])

	# check reco level distribution
	# FIXME use TUnfoldDensityV17::GetBackground()
        outpdf = outputDirectory + "recoPt.png"
        drawutil.recoPt(outpdf, postfix, unfoldInputList['data'], unfoldInputList['sig'], unfoldInputList['DYtoEEtau'], unfoldInputList['TTbar'], unfoldInputList['VV'], unfoldInputList['Wjets'], args.channel);

        outpdf = outputDirectory + "recoMass.png"
        drawutil.recoMass(outpdf, postfix, unfoldInputList['data'], unfoldInputList['sig'], unfoldInputList['DYtoEEtau'], unfoldInputList['TTbar'], unfoldInputList['VV'], unfoldInputList['Wjets'], args.channel);

	#outpdf = outputDirectory + "test.png" 
	#drawutil.drawTest(outpdf, unfold_pt)

	if args.doSeperateUnfold:
		# FSR unfolding
		# get unfolded histogram from detector unfolding results
		hunfolded_pt = unfoldutil.getUnfoldedHist(unfold_pt, "detector")
		hunfolded_mass = unfoldutil.getUnfoldedHist(unfold_mass, "detector")

		# set migration matrix for FSR unfolding
        	unfoldFSR_pt = unfoldutil.setUnfold(unfoldInputList['sig'],   "Pt",   "FSR", True)
        	unfoldFSR_mass = unfoldutil.setUnfold(unfoldInputList['sig'], "Mass", "FSR", True)

        	outpdf = outputDirectory + "response_FSR.png"
        	drawutil.responseM(outpdf, unfoldFSR_pt)
        	outpdf = outputDirectory + "efficiency_FSR.png"
        	drawutil.efficiency(outpdf, unfoldFSR_pt)

        	outpdf_mass = outputDirectory + "response_mass_FSR.png"
        	drawutil.responseM(outpdf_mass, unfoldFSR_mass)
        	outpdf_mass = outputDirectory + "efficiency_mass_FSR.png"
        	drawutil.efficiency(outpdf_mass, unfoldFSR_mass)

		unfoldutil.setUnfoldInputHist(unfoldFSR_pt, hunfolded_pt)
		unfoldutil.setUnfoldInputHist(unfoldFSR_mass, hunfolded_mass)

        	unfoldutil.doUnfold(unfoldFSR_pt)
        	unfoldutil.doUnfold(unfoldFSR_mass)

        	outpdf = outputDirectory + "ratio_FSR.png"
        	outpdf_mass = outputDirectory + "ratio_mass_FSR.png"
        	# check unfolded distribution
        	drawutil.basicRatio(outpdf, unfoldFSR_pt, unfoldFSR_mass, unfoldInputList['sig'])
        	drawutil.basicRatioMass(outpdf_mass, unfoldFSR_mass, unfoldInputList['sig'])
		
		del unfoldFSR_pt
		del unfoldFSR_mass


        del unfold_pt
        del unfold_mass

def makeRecoPlots():
        pass

# test 
if __name__ == "__main__": 
	makeRecoPlots()
