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
parser.add_argument('--altResponse'  , action='store_true'  , default = False, help = 'create response matrix with alternative signal MC')
parser.add_argument('--getUnfoldResults'  , action='store_true'  , help = 'Get unfolding resutls')
parser.add_argument('--getUnfoldResultsV2'  , action='store_true'  , help = 'Get unfolding resutls ver. 2')
parser.add_argument('--doSys'  , action='store_true'  , default = False, help = 'Calculate systematics')
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
		if args.altResponse and not sample.isAlt: continue
	    	print 'creating histogram for sample '
	    	sample.dump()
		inputfHistDic.update(histUtil.makeSigHists(sample, outputDirectory, args.channel))
	
	fOutTxt = open( inputfhisttxtName,'w')
	
	for fhist in inputfHistDic.keys():
	        astr = inputfHistDic[fhist].name + ' ' + inputfHistDic[fhist].htype + ' ' + (inputfHistDic[fhist].path)[0]
	        fOutTxt.write( astr + '\n' )
	
	fOutTxt.close()


if args.getUnfoldResultsV2:
        # read text file including input histograms for unfolding
        fOutTxtCheck = open( inputfhisttxtName,'r')
        unfoldInputList = {}

        for line in fOutTxtCheck:
                modifiedLine = line.lstrip(' ').rstrip(' ').rstrip('\n')
                if modifiedLine.split()[1] == "data":
                        unfoldInputList['data'] = modifiedLine.split()[2]
                if modifiedLine.split()[1] == "sig":
                        unfoldInputList['sig'] = modifiedLine.split()[2]
                if modifiedLine.split()[1] == "sigAlt":
                        unfoldInputList['sigAlt'] = modifiedLine.split()[2]
                if modifiedLine.split()[1] == "matrix":
                        unfoldInputList['matrix'] = modifiedLine.split()[2]
                if modifiedLine.split()[1] == "bkg": # use the sample name as keyword for background
                        unfoldInputList[modifiedLine.split()[0]] = modifiedLine.split()[2]

        print unfoldInputList

        import pyScripts.unfoldUtil as unfoldutil
        import pyScripts.drawUtil as drawutil

        postfix = args.postfix

        unfoldClass = rt.ISRUnfold()

        # set response matrix
        unfoldClass.setNomTUnfoldDensity(unfoldInputList['matrix'], "Pt",     postfix, True)
        unfoldClass.setNomTUnfoldDensity(unfoldInputList['matrix'], "Mass",   postfix, True)

        #unfoldClass.setNomTUnfoldDensity(unfoldInputList['sig'], "Pt",     postfix, False)
        #unfoldClass.setNomTUnfoldDensity(unfoldInputList['sig'], "Mass",   postfix, False)

        # set unfolding input histogram
	unfoldClass.setInput("Pt",   postfix, unfoldInputList['data'])
	unfoldClass.setInput("Mass", postfix, unfoldInputList['data'])

	# set background histograms 
        unfoldClass.subBkgs("Pt", postfix, unfoldInputList['DYtoEEtau'], "DYtoEEtau")
        unfoldClass.subBkgs("Pt", postfix, unfoldInputList['TTbar'],     "TTbar")
        unfoldClass.subBkgs("Pt", postfix, unfoldInputList['VV'],        "VV")
        unfoldClass.subBkgs("Pt", postfix, unfoldInputList['Wjets'],     "Wjets")

        unfoldClass.subBkgs("Mass", postfix, unfoldInputList['DYtoEEtau'], "DYtoEEtau")
        unfoldClass.subBkgs("Mass", postfix, unfoldInputList['TTbar'],     "TTbar")
        unfoldClass.subBkgs("Mass", postfix, unfoldInputList['VV'],        "VV")
        unfoldClass.subBkgs("Mass", postfix, unfoldInputList['Wjets'],     "Wjets")


        if args.doSys == True:
	    sysDict = {"PU": 2, "trgSF": 2, "recoSF": 2, "IdSF": 2, "IsoSF": 2, "unfoldsys": 1, "AlphaS": 2, "Scale": 9, "PDFerror": 100, "Alt": 1, "L1Prefire": 2, "LepScale": 2, "LepRes": 2, "FSRDR": 30, "unfoldBias": 1, "unfoldScan": 1, "Closure": 1}
	    #sysDict = {"PU": 2, "trgSF": 2, "recoSF": 2, "IdSF": 2, "IsoSF": 2, "L1Prefire": 2, "Closure": 1,"PDFerror": 100}

	    for sysName, nSys in sysDict.items():
	    	for nthSys in range(0,nSys):

	    		print "sysName: " + sysName + " nthSys: " + str(nthSys) + " #####################################################"

	    		# set systematic response matrix
	    		postfix = sysName
            		if sysName != "Alt": 
                                # use systematic response matrix saved in the signal root file 
	    			unfoldClass.setSysTUnfoldDensity(unfoldInputList['sig'], "Pt",     postfix, nthSys, False)
            			unfoldClass.setSysTUnfoldDensity(unfoldInputList['sig'], "Mass",   postfix, nthSys, False)
	    		if sysName == "Alt":
	    			unfoldClass.setSysTUnfoldDensity(unfoldInputList['sigAlt'], "Pt",     postfix, nthSys, False)
            			unfoldClass.setSysTUnfoldDensity(unfoldInputList['sigAlt'], "Mass",   postfix, nthSys, False)

	    		bias = 1.;
	    		if sysName == "unfoldBias": bias = 0.95 
	    	
	    		if sysName != "Closure":	
                                # set systematic input histograms
	    			unfoldClass.setInput("Pt",   postfix, unfoldInputList['data'], nthSys, True, bias)
	    			unfoldClass.setInput("Mass", postfix, unfoldInputList['data'], nthSys, True, bias)

                                # set systematic background histograms
            			unfoldClass.subBkgs("Pt", postfix, unfoldInputList['DYtoEEtau'], "DYtoEEtau", nthSys, True)
            			unfoldClass.subBkgs("Pt", postfix, unfoldInputList['TTbar'],     "TTbar",     nthSys, True)
            			unfoldClass.subBkgs("Pt", postfix, unfoldInputList['VV'],        "VV",        nthSys, True)
            			unfoldClass.subBkgs("Pt", postfix, unfoldInputList['Wjets'],     "Wjets",     nthSys, True)

            			unfoldClass.subBkgs("Mass", postfix, unfoldInputList['DYtoEEtau'], "DYtoEEtau", nthSys, True)
            			unfoldClass.subBkgs("Mass", postfix, unfoldInputList['TTbar'],     "TTbar",     nthSys, True)
            			unfoldClass.subBkgs("Mass", postfix, unfoldInputList['VV'],        "VV",        nthSys, True)
            			unfoldClass.subBkgs("Mass", postfix, unfoldInputList['Wjets'],     "Wjets",     nthSys, True)

	    		if sysName == "Closure":
                                # for closure test, use signal histogram as unfolding input
	    			unfoldClass.setInput("Pt",   postfix, unfoldInputList['sig'], nthSys, True, bias)
	    			unfoldClass.setInput("Mass", postfix, unfoldInputList['sig'], nthSys, True, bias)

	# unfold up to pre FSR level 
	unfoldClass.doISRUnfold(args.doSys)

	# set nominal value and also systematic values
	unfoldClass.setMeanPt(args.doSys, False)
	unfoldClass.setMeanMass(args.doSys, False)

        for massBin in range(0,5):
            
                unfoldClass.drawClosurePlots(outputDirectory + "Closure_"+args.channel, "Pt", massBin); 
                unfoldClass.drawClosurePlots(outputDirectory + "Closure_"+args.channel, "Mass", massBin); 

                unfoldClass.drawNominalPlots(outputDirectory + "Unfolded_"+args.channel, "Pt", massBin);
                unfoldClass.drawNominalPlots(outputDirectory + "Unfolded_"+args.channel, "Mass", massBin);



        # draw plots including systematic 
        for sysName, nSys in sysDict.items():

            if sysName == "Closure" : continue

	    for massBin in range(0,5):

                unfoldClass.drawNominalPlots(outputDirectory + "Unfolded_"+args.channel+sysName, "Pt", massBin, sysName, args.doSys);
                unfoldClass.drawNominalPlots(outputDirectory + "Unfolded_"+args.channel+sysName, "Mass", massBin, sysName, args.doSys);
                unfoldClass.drawInputPlots(outputDirectory + args.channel + sysName, "Pt", massBin, sysName);
		unfoldClass.drawSysPlots(outputDirectory + "Sys_" + args.channel , massBin, sysName)


	#unfoldClass.drawNominalRecoPlots(outputDirectory + "RecoPt_" + args.channel + ".pdf", unfoldInputList['sig'], "Pt", 0)
	#unfoldClass.drawNominalRecoPlots(outputDirectory + "RecoPt_" + args.channel + ".pdf", unfoldInputList['sig'], "Pt", 1)
	#unfoldClass.drawNominalRecoPlots(outputDirectory + "RecoPt_" + args.channel + ".pdf", unfoldInputList['sig'], "Pt", 2)
	#unfoldClass.drawNominalRecoPlots(outputDirectory + "RecoPt_" + args.channel + ".pdf", unfoldInputList['sig'], "Pt", 3)
	#unfoldClass.drawNominalRecoPlots(outputDirectory + "RecoPt_" + args.channel + ".pdf", unfoldInputList['sig'], "Pt", 4)

	#unfoldClass.studyFSRDRPlots(outputDirectory + "FSRPt_" + args.channel, "Pt", 0)
	#unfoldClass.studyFSRDRPlots(outputDirectory + "FSRPt_" + args.channel, "Pt", 1)
	#unfoldClass.studyFSRDRPlots(outputDirectory + "FSRPt_" + args.channel, "Pt", 2)
	#unfoldClass.studyFSRDRPlots(outputDirectory + "FSRPt_" + args.channel, "Pt", 3)
	#unfoldClass.studyFSRDRPlots(outputDirectory + "FSRPt_" + args.channel, "Pt", 4)
	#unfoldClass.studyFSRDRPlots(outputDirectory + "FSRPt_" + args.channel, "Mass", -1)

        unfoldClass.drawLCurve(outputDirectory + "LCurve_" + args.channel + ".pdf", "Pt")
        unfoldClass.drawLCurve(outputDirectory + "LCurveMass_" + args.channel + ".pdf", "Mass")

        unfoldClass.drawRhoLog(outputDirectory + "RhoLog_" + args.channel + ".pdf", "Pt")
        unfoldClass.drawRhoLog(outputDirectory + "RhoLogMass_" + args.channel + ".pdf", "Mass")


	unfoldClass.drawISRresult(outputDirectory + "ISRfit_" + args.channel + ".pdf", False)

	del unfoldClass

def makeRecoPlots():
        pass
# test 
if __name__ == "__main__": 
	makeRecoPlots()
