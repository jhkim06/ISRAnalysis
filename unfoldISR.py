import argparse
import os
import sys
import ROOT as rt

parser = argparse.ArgumentParser(description='Unfolding for ISR analysis')

parser.add_argument('--channel' , dest = 'channel', default = 'electron', help = 'select channel electron or muon')
parser.add_argument('--year' , dest = 'year', default = '2016', help = 'select year')
parser.add_argument('--postfix' , dest = 'postfix', default = 'nominal', help = 'select histogram name')
parser.add_argument('--createInputHists'  , action='store_true'  , help = 'create input histograms')
parser.add_argument('--createMatrixOnly'  , action='store_true'  , default = False, help = 'create histograms only for signal sample')
parser.add_argument('--getUnfoldResults'  , action='store_true'  , help = 'Get unfolding resutls')
parser.add_argument('--getUnfoldResultsV2'  , action='store_true'  , help = 'Get unfolding resutls ver. 2')
parser.add_argument('--closure'  , action='store_true'  , help = 'Clousre test with MC')
parser.add_argument('--getCombinedResults'  , action='store_true'  , help = 'Combine 2016 and 2017')
parser.add_argument('--getEMuCombinedResults'  , action='store_true'  , help = 'Combine electron and muon')
parser.add_argument('--doSeperateUnfold'  , action='store_true'  , default = False, help = 'Seperate unfolding steps')

args = parser.parse_args()

# get input root files information using sampleDef.py
import etc.sampleDef as isrSamples
import pyScripts.unfoldInputUtil as histUtil

outputDirectory = 'output/'+args.year+'/' + args.channel + "/"  
inputfhisttxtName = outputDirectory + "fhist.txt"
if not os.path.exists( outputDirectory ):
	os.makedirs( outputDirectory )

print "channel to run: " + args.channel

selectedSample = {} 
if args.channel == "electron":
	if args.year == "2016":
		selectedSample = isrSamples.samplesDef_electron2016Legacy
	if args.year == "2017":
		selectedSample = isrSamples.samplesDef_electron2017Legacy
if args.channel == "muon":
	if args.year == "2016":
		selectedSample = isrSamples.samplesDef_muon2016Legacy 
	if args.year == "2017":
		selectedSample = isrSamples.samplesDef_muon2017Legacy 

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


if args.getUnfoldResultsV2:
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

        postfix = args.postfix
        postfix_matrix = args.postfix

        unfoldClass = rt.ISRUnfold()
        unfoldClass.setNomTUnfoldDensity(unfoldInputList['sig'], "Pt",     postfix_matrix, False)
        unfoldClass.setNomTUnfoldDensity(unfoldInputList['sig'], "Mass",   postfix_matrix, False)

	unfoldClass.setInput("Pt",   postfix, unfoldInputList['data'])
	unfoldClass.setInput("Mass", postfix, unfoldInputList['data'])

	# 
        unfoldClass.subBkgs("Pt", postfix, unfoldInputList['DYtoEEtau'], "DYtoEEtau")
        unfoldClass.subBkgs("Pt", postfix, unfoldInputList['TTbar'], "TTbar")
        unfoldClass.subBkgs("Pt", postfix, unfoldInputList['VV'], "VV")
        unfoldClass.subBkgs("Pt", postfix, unfoldInputList['Wjets'], "Wjets")

        unfoldClass.subBkgs("Mass", postfix, unfoldInputList['DYtoEEtau'], "DYtoEEtau")
        unfoldClass.subBkgs("Mass", postfix, unfoldInputList['TTbar'], "TTbar")
        unfoldClass.subBkgs("Mass", postfix, unfoldInputList['VV'], "VV")
        unfoldClass.subBkgs("Mass", postfix, unfoldInputList['Wjets'], "Wjets")


	sysDict = {"PU": 2, "trgSF": 2, "recoSF": 2, "IdSF": 2, "IsoSF": 2, "unfoldsys": 1, "AlphaS": 2, "Scale": 9, "PDFerror": 100}
	#sysDict = {"recoSF": 2}

	for sysName, nSys in sysDict.items():
		for nthSys in range(0,nSys):

			print "sysName: " + sysName + " nthSys: " + str(nthSys) + " #####################################################"

			# systematic test
			postfix = sysName
        		unfoldClass.setSysTUnfoldDensity(unfoldInputList['sig'], "Pt",     postfix, nthSys, False)
        		unfoldClass.setSysTUnfoldDensity(unfoldInputList['sig'], "Mass",   postfix, nthSys, False)

			unfoldClass.setInput("Pt",   postfix, unfoldInputList['data'], nthSys, True)
			unfoldClass.setInput("Mass", postfix, unfoldInputList['data'], nthSys, True)

        		unfoldClass.subBkgs("Pt", postfix, unfoldInputList['DYtoEEtau'], "DYtoEEtau", nthSys, True)
        		unfoldClass.subBkgs("Pt", postfix, unfoldInputList['TTbar'], "TTbar", nthSys, True)
        		unfoldClass.subBkgs("Pt", postfix, unfoldInputList['VV'], "VV", nthSys, True)
        		unfoldClass.subBkgs("Pt", postfix, unfoldInputList['Wjets'], "Wjets", nthSys, True)

        		unfoldClass.subBkgs("Mass", postfix, unfoldInputList['DYtoEEtau'], "DYtoEEtau", nthSys, True)
        		unfoldClass.subBkgs("Mass", postfix, unfoldInputList['TTbar'], "TTbar", nthSys, True)
        		unfoldClass.subBkgs("Mass", postfix, unfoldInputList['VV'], "VV", nthSys, True)
        		unfoldClass.subBkgs("Mass", postfix, unfoldInputList['Wjets'], "Wjets", nthSys, True)

	unfoldClass.doISRUnfold()

	unfoldClass.setMeanPt()
	unfoldClass.setMeanMass()

        for sysName, nSys in sysDict.items():

                        unfoldClass.drawNominalPlots(outputDirectory + "UnfoldedPt_"+args.channel+sysName, "Pt", 0, sysName);
                        unfoldClass.drawNominalPlots(outputDirectory + "UnfoldedPt_"+args.channel+sysName, "Pt", 1, sysName);
                        unfoldClass.drawNominalPlots(outputDirectory + "UnfoldedPt_"+args.channel+sysName, "Pt", 2, sysName);
                        unfoldClass.drawNominalPlots(outputDirectory + "UnfoldedPt_"+args.channel+sysName, "Pt", 3, sysName);
                        unfoldClass.drawNominalPlots(outputDirectory + "UnfoldedPt_"+args.channel+sysName, "Pt", 4, sysName);

                        unfoldClass.drawInputPlots(outputDirectory + "UnfoldedPt_"+args.channel+sysName, "Pt", 0, sysName);
                        unfoldClass.drawInputPlots(outputDirectory + "UnfoldedPt_"+args.channel+sysName, "Pt", 1, sysName);
                        unfoldClass.drawInputPlots(outputDirectory + "UnfoldedPt_"+args.channel+sysName, "Pt", 2, sysName);
                        unfoldClass.drawInputPlots(outputDirectory + "UnfoldedPt_"+args.channel+sysName, "Pt", 3, sysName);
                        unfoldClass.drawInputPlots(outputDirectory + "UnfoldedPt_"+args.channel+sysName, "Pt", 4, sysName);

	#unfoldClass.drawNominalPlots("Mass");
	del unfoldClass

def makeRecoPlots():
        pass
# test 
if __name__ == "__main__": 
	makeRecoPlots()
