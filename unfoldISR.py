import argparse
import os
import sys
import ROOT

parser = argparse.ArgumentParser(description='Unfolding for ISR analysis')

parser.add_argument('--channel' , dest = 'channel', default = 'electron', help = 'select channel electron or muon')
parser.add_argument('--year' , dest = 'year', default = '2016', help = 'select year')
parser.add_argument('--postfix' , dest = 'postfix', default = 'nominal', help = 'select histogram name')
parser.add_argument('--createInputHists'  , action='store_true'  , help = 'create input histograms')
parser.add_argument('--createMatrixOnly'  , action='store_true'  , default = False, help = 'create histograms only for signal sample')
parser.add_argument('--getUnfoldResults'  , action='store_true'  , help = 'Get unfolding resutls')
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


if args.getCombinedResults:

	outputDirectory2016 = 'output/2016/' + args.channel + "/"
	inputfhisttxtName2016 = outputDirectory2016 + "fhist.txt"

	# read text file
	fOutTxtCheck2016 = open(inputfhisttxtName2016,'r')
	unfoldInputList2016 = {}
	
	for line in fOutTxtCheck2016:
	        modifiedLine = line.lstrip(' ').rstrip(' ').rstrip('\n')
	        if modifiedLine.split(' ')[1] == "data":
	                unfoldInputList2016['data'] = modifiedLine.split(' ')[2]
	        if modifiedLine.split(' ')[1] == "sig":
	                unfoldInputList2016['sig'] = modifiedLine.split(' ')[2]
	        if modifiedLine.split(' ')[1] == "bkg": # use the sample name as keyword for background
	                unfoldInputList2016[modifiedLine.split(' ')[0]] = modifiedLine.split(' ')[2]
	
	print unfoldInputList2016

        outputDirectory2017 = 'output/2017/' + args.channel + "/"
        inputfhisttxtName2017 = outputDirectory2017 + "fhist.txt"

        # read text file
        fOutTxtCheck2017 = open(inputfhisttxtName2017,'r')
        unfoldInputList2017 = {}
        
        for line in fOutTxtCheck2017:
                modifiedLine = line.lstrip(' ').rstrip(' ').rstrip('\n')
                if modifiedLine.split(' ')[1] == "data":
                        unfoldInputList2017['data'] = modifiedLine.split(' ')[2]
                if modifiedLine.split(' ')[1] == "sig":
                        unfoldInputList2017['sig'] = modifiedLine.split(' ')[2]
                if modifiedLine.split(' ')[1] == "bkg": # use the sample name as keyword for background
                        unfoldInputList2017[modifiedLine.split(' ')[0]] = modifiedLine.split(' ')[2]
        
        print unfoldInputList2017

	import pyScripts.unfoldUtil as unfoldutil
	import pyScripts.drawUtil as drawutil

        postfix = args.postfix
        postfix_matrix = args.postfix

	## 2016
        # set unfolding class 
        unfold_pt2016 =   unfoldutil.setUnfold(unfoldInputList2016['sig'], "Pt",   postfix_matrix, False)
        unfold_mass2016 = unfoldutil.setUnfold(unfoldInputList2016['sig'], "Mass", postfix_matrix, False)

        # get bkg subtracted data
        unfoldutil.setUnfoldInput(unfold_pt2016, "Pt", postfix, unfoldInputList2016['data'])
        #unfoldutil.setUnfoldInput(unfold_pt2016, "Pt", unfoldInputList2016['sig'])

        unfoldutil.subtractBkgs(unfold_pt2016, "Pt", postfix, unfoldInputList2016['DYtoEEtau'], "DYtoEEtau")
        unfoldutil.subtractBkgs(unfold_pt2016, "Pt", postfix, unfoldInputList2016['TTbar'], "TTbar")
        unfoldutil.subtractBkgs(unfold_pt2016, "Pt", postfix, unfoldInputList2016['VV'], "VV")
        unfoldutil.subtractBkgs(unfold_pt2016, "Pt", postfix, unfoldInputList2016['Wjets'], "Wjets")

        # for mass unfolding
        unfoldutil.setUnfoldInput(unfold_mass2016, "Mass", postfix, unfoldInputList2016['data'])
        #unfoldutil.setUnfoldInput(unfold_mass2016, "Mass", unfoldInputList2016['sig'])

        unfoldutil.subtractBkgs(unfold_mass2016, "Mass", postfix, unfoldInputList2016['DYtoEEtau'], "DYtoEEtau")
        unfoldutil.subtractBkgs(unfold_mass2016, "Mass", postfix, unfoldInputList2016['TTbar'], "TTbar")
        unfoldutil.subtractBkgs(unfold_mass2016, "Mass", postfix, unfoldInputList2016['VV'], "VV")
        unfoldutil.subtractBkgs(unfold_mass2016, "Mass", postfix, unfoldInputList2016['Wjets'], "Wjets")

        # do unfolding with bkg subtracted data and the migration matrix
        unfoldutil.doUnfold(unfold_pt2016)
        unfoldutil.doUnfold(unfold_mass2016)

        unfoldutil.setVectorSys(unfoldInputList2016['sig'], unfold_pt2016,   "Pt",   "AlphaS", 2)
        unfoldutil.setVectorSys(unfoldInputList2016['sig'], unfold_mass2016, "Mass", "AlphaS", 2)

        unfoldutil.setVectorSys(unfoldInputList2016['sig'], unfold_pt2016,   "Pt",   "Scale", 9)
        unfoldutil.setVectorSys(unfoldInputList2016['sig'], unfold_mass2016, "Mass", "Scale", 9)

	## 2017
        # set unfolding class 
        unfold_pt2017 =   unfoldutil.setUnfold(unfoldInputList2017['sig'], "Pt",   postfix_matrix, False)
        unfold_mass2017 = unfoldutil.setUnfold(unfoldInputList2017['sig'], "Mass", postfix_matrix, False)

        # get bkg subtracted data
        unfoldutil.setUnfoldInput(unfold_pt2017, "Pt", postfix, unfoldInputList2017['data'])
        #unfoldutil.setUnfoldInput(unfold_pt2017, "Pt", unfoldInputList2017['sig'])

        unfoldutil.subtractBkgs(unfold_pt2017, "Pt", postfix, unfoldInputList2017['DYtoEEtau'], "DYtoEEtau")
        unfoldutil.subtractBkgs(unfold_pt2017, "Pt", postfix, unfoldInputList2017['TTbar'], "TTbar")
        unfoldutil.subtractBkgs(unfold_pt2017, "Pt", postfix, unfoldInputList2017['VV'], "VV")
        unfoldutil.subtractBkgs(unfold_pt2017, "Pt", postfix, unfoldInputList2017['Wjets'], "Wjets")

        # for mass unfolding
        unfoldutil.setUnfoldInput(unfold_mass2017, "Mass", postfix, unfoldInputList2017['data'])
        #unfoldutil.setUnfoldInput(unfold_mass2017, "Mass", unfoldInputList2017['sig'])

        unfoldutil.subtractBkgs(unfold_mass2017, "Mass", postfix, unfoldInputList2017['DYtoEEtau'], "DYtoEEtau")
        unfoldutil.subtractBkgs(unfold_mass2017, "Mass", postfix, unfoldInputList2017['TTbar'], "TTbar")
        unfoldutil.subtractBkgs(unfold_mass2017, "Mass", postfix, unfoldInputList2017['VV'], "VV")
        unfoldutil.subtractBkgs(unfold_mass2017, "Mass", postfix, unfoldInputList2017['Wjets'], "Wjets")

        # do unfolding with bkg subtracted data and the migration matrix
        unfoldutil.doUnfold(unfold_pt2017)
        unfoldutil.doUnfold(unfold_mass2017)

        unfoldutil.setVectorSys(unfoldInputList2017['sig'], unfold_pt2017,   "Pt",   "AlphaS", 2)
        unfoldutil.setVectorSys(unfoldInputList2017['sig'], unfold_mass2017, "Mass", "AlphaS", 2)

        unfoldutil.setVectorSys(unfoldInputList2017['sig'], unfold_pt2017,   "Pt",   "Scale", 9)
        unfoldutil.setVectorSys(unfoldInputList2017['sig'], unfold_mass2017, "Mass", "Scale", 9)

	drawutil.drawCombinedISR("./test_"+args.channel+".pdf", unfold_pt2016, unfold_mass2016, unfold_pt2017, unfold_mass2017 )


if args.getEMuCombinedResults:

	outputDirectoryElectron = 'output/2016/' + "electron" + "/"
	inputfhisttxtNameElectron = outputDirectoryElectron + "fhist.txt"

	# read text file
	fOutTxtCheckElectron = open(inputfhisttxtNameElectron,'r')
	unfoldInputListElectron = {}
	
	for line in fOutTxtCheckElectron:
	        modifiedLine = line.lstrip(' ').rstrip(' ').rstrip('\n')
	        if modifiedLine.split(' ')[1] == "data":
	                unfoldInputListElectron['data'] = modifiedLine.split(' ')[2]
	        if modifiedLine.split(' ')[1] == "sig":
	                unfoldInputListElectron['sig'] = modifiedLine.split(' ')[2]
	        if modifiedLine.split(' ')[1] == "bkg": # use the sample name as keyword for background
	                unfoldInputListElectron[modifiedLine.split(' ')[0]] = modifiedLine.split(' ')[2]
	
	print unfoldInputListElectron

        outputDirectoryMuon = 'output/2016/' + "muon" + "/"
        inputfhisttxtNameMuon = outputDirectoryMuon + "fhist.txt"

        # read text file
        fOutTxtCheckMuon = open(inputfhisttxtNameMuon,'r')
        unfoldInputListMuon = {}
        
        for line in fOutTxtCheckMuon:
                modifiedLine = line.lstrip(' ').rstrip(' ').rstrip('\n')
                if modifiedLine.split(' ')[1] == "data":
                        unfoldInputListMuon['data'] = modifiedLine.split(' ')[2]
                if modifiedLine.split(' ')[1] == "sig":
                        unfoldInputListMuon['sig'] = modifiedLine.split(' ')[2]
                if modifiedLine.split(' ')[1] == "bkg": # use the sample name as keyword for background
                        unfoldInputListMuon[modifiedLine.split(' ')[0]] = modifiedLine.split(' ')[2]
        
        print unfoldInputListMuon

	import pyScripts.unfoldUtil as unfoldutil
	import pyScripts.drawUtil as drawutil

        postfix = args.postfix
        postfix_matrix = args.postfix

	## Electron
        # set unfolding class 
        unfold_ptElectron =   unfoldutil.setUnfold(unfoldInputListElectron['sig'], "Pt",   postfix_matrix, False)
        unfold_massElectron = unfoldutil.setUnfold(unfoldInputListElectron['sig'], "Mass", postfix_matrix, False)

        # get bkg subtracted data
        unfoldutil.setUnfoldInput(unfold_ptElectron, "Pt", postfix, unfoldInputListElectron['data'])
        #unfoldutil.setUnfoldInput(unfold_ptElectron, "Pt", unfoldInputListElectron['sig'])

        unfoldutil.subtractBkgs(unfold_ptElectron, "Pt", postfix, unfoldInputListElectron['DYtoEEtau'], "DYtoEEtau")
        unfoldutil.subtractBkgs(unfold_ptElectron, "Pt", postfix, unfoldInputListElectron['TTbar'], "TTbar")
        unfoldutil.subtractBkgs(unfold_ptElectron, "Pt", postfix, unfoldInputListElectron['VV'], "VV")
        unfoldutil.subtractBkgs(unfold_ptElectron, "Pt", postfix, unfoldInputListElectron['Wjets'], "Wjets")

        # for mass unfolding
        unfoldutil.setUnfoldInput(unfold_massElectron, "Mass", postfix, unfoldInputListElectron['data'])
        #unfoldutil.setUnfoldInput(unfold_massElectron, "Mass", unfoldInputListElectron['sig'])

        unfoldutil.subtractBkgs(unfold_massElectron, "Mass", postfix, unfoldInputListElectron['DYtoEEtau'], "DYtoEEtau")
        unfoldutil.subtractBkgs(unfold_massElectron, "Mass", postfix, unfoldInputListElectron['TTbar'], "TTbar")
        unfoldutil.subtractBkgs(unfold_massElectron, "Mass", postfix, unfoldInputListElectron['VV'], "VV")
        unfoldutil.subtractBkgs(unfold_massElectron, "Mass", postfix, unfoldInputListElectron['Wjets'], "Wjets")

        # do unfolding with bkg subtracted data and the migration matrix
        unfoldutil.doUnfold(unfold_ptElectron)
        unfoldutil.doUnfold(unfold_massElectron)

        unfoldutil.setVectorSys(unfoldInputListElectron['sig'], unfold_ptElectron,   "Pt",   "AlphaS", 2)
        unfoldutil.setVectorSys(unfoldInputListElectron['sig'], unfold_massElectron, "Mass", "AlphaS", 2)

        unfoldutil.setVectorSys(unfoldInputListElectron['sig'], unfold_ptElectron,   "Pt",   "Scale", 9)
        unfoldutil.setVectorSys(unfoldInputListElectron['sig'], unfold_massElectron, "Mass", "Scale", 9)

	## Muon
        # set unfolding class 
        unfold_ptMuon =   unfoldutil.setUnfold(unfoldInputListMuon['sig'], "Pt",   postfix_matrix, False)
        unfold_massMuon = unfoldutil.setUnfold(unfoldInputListMuon['sig'], "Mass", postfix_matrix, False)

        # get bkg subtracted data
        unfoldutil.setUnfoldInput(unfold_ptMuon, "Pt", postfix, unfoldInputListMuon['data'])
        #unfoldutil.setUnfoldInput(unfold_ptMuon, "Pt", unfoldInputListMuon['sig'])

        unfoldutil.subtractBkgs(unfold_ptMuon, "Pt", postfix, unfoldInputListMuon['DYtoEEtau'], "DYtoEEtau")
        unfoldutil.subtractBkgs(unfold_ptMuon, "Pt", postfix, unfoldInputListMuon['TTbar'], "TTbar")
        unfoldutil.subtractBkgs(unfold_ptMuon, "Pt", postfix, unfoldInputListMuon['VV'], "VV")
        unfoldutil.subtractBkgs(unfold_ptMuon, "Pt", postfix, unfoldInputListMuon['Wjets'], "Wjets")

        # for mass unfolding
        unfoldutil.setUnfoldInput(unfold_massMuon, "Mass", postfix, unfoldInputListMuon['data'])
        #unfoldutil.setUnfoldInput(unfold_massMuon, "Mass", unfoldInputListMuon['sig'])

        unfoldutil.subtractBkgs(unfold_massMuon, "Mass", postfix, unfoldInputListMuon['DYtoEEtau'], "DYtoEEtau")
        unfoldutil.subtractBkgs(unfold_massMuon, "Mass", postfix, unfoldInputListMuon['TTbar'], "TTbar")
        unfoldutil.subtractBkgs(unfold_massMuon, "Mass", postfix, unfoldInputListMuon['VV'], "VV")
        unfoldutil.subtractBkgs(unfold_massMuon, "Mass", postfix, unfoldInputListMuon['Wjets'], "Wjets")

        # do unfolding with bkg subtracted data and the migration matrix
        unfoldutil.doUnfold(unfold_ptMuon)
        unfoldutil.doUnfold(unfold_massMuon)

        unfoldutil.setVectorSys(unfoldInputListMuon['sig'], unfold_ptMuon,   "Pt",   "AlphaS", 2)
        unfoldutil.setVectorSys(unfoldInputListMuon['sig'], unfold_massMuon, "Mass", "AlphaS", 2)

        unfoldutil.setVectorSys(unfoldInputListMuon['sig'], unfold_ptMuon,   "Pt",   "Scale", 9)
        unfoldutil.setVectorSys(unfoldInputListMuon['sig'], unfold_massMuon, "Mass", "Scale", 9)

	drawutil.drawEMuCombinedISR("./EMu.pdf", unfold_ptElectron, unfold_massElectron, unfold_ptMuon, unfold_massMuon )




if args.getUnfoldResults:

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
	#postfix_matrix = "unfoldsys_0"
	if args.doSeperateUnfold:
		postfix = "nominal"
		postfix_matrix = "detector"

        # set unfolding class 
	unfold_pt =   unfoldutil.setUnfold(unfoldInputList['sig'], "Pt",   postfix_matrix, False)
	unfold_mass = unfoldutil.setUnfold(unfoldInputList['sig'], "Mass", postfix_matrix, False)

	# print out response matrix and efficiency plot
        outpdf = outputDirectory + "response_"+args.channel+".pdf"
        drawutil.responseM(outpdf, unfold_pt)
        outpdf = outputDirectory + "efficiency_"+args.channel+".pdf"
        drawutil.efficiency(outpdf, unfold_pt)

        outpdf_mass = outputDirectory + "response_mass_"+args.channel+".pdf"
        drawutil.responseM(outpdf_mass, unfold_mass)
        outpdf_mass = outputDirectory + "efficiency_mass_"+args.channel+".pdf"
        drawutil.efficiency(outpdf_mass, unfold_mass)

	# get bkg subtracted data
        if not args.closure: 
		unfoldutil.setUnfoldInput(unfold_pt, "Pt", postfix, unfoldInputList['data'])
        	unfoldutil.subtractBkgs(unfold_pt, "Pt", postfix, unfoldInputList['DYtoEEtau'], "DYtoEEtau")
        	unfoldutil.subtractBkgs(unfold_pt, "Pt", postfix, unfoldInputList['TTbar'], "TTbar")
        	unfoldutil.subtractBkgs(unfold_pt, "Pt", postfix, unfoldInputList['VV'], "VV")
        	unfoldutil.subtractBkgs(unfold_pt, "Pt", postfix, unfoldInputList['Wjets'], "Wjets")
        	if args.channel == "muon": unfoldutil.subtractBkgs(unfold_pt, "Pt", postfix, unfoldInputList['QCD'], "QCD")
        else :
		unfoldutil.setUnfoldInput(unfold_pt, "Pt", postfix, unfoldInputList['sig'])


	# for mass unfolding
        if not args.closure:
		unfoldutil.setUnfoldInput(unfold_mass, "Mass", postfix, unfoldInputList['data'])
        	unfoldutil.subtractBkgs(unfold_mass, "Mass", postfix, unfoldInputList['DYtoEEtau'], "DYtoEEtau")
        	unfoldutil.subtractBkgs(unfold_mass, "Mass", postfix, unfoldInputList['TTbar'], "TTbar")
        	unfoldutil.subtractBkgs(unfold_mass, "Mass", postfix, unfoldInputList['VV'], "VV")
        	unfoldutil.subtractBkgs(unfold_mass, "Mass", postfix, unfoldInputList['Wjets'], "Wjets")
        	if args.channel == "muon": unfoldutil.subtractBkgs(unfold_pt, "Mass", postfix, unfoldInputList['QCD'], "QCD")
        else :
		unfoldutil.setUnfoldInput(unfold_mass, "Mass", postfix, unfoldInputList['sig'])


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

        unfoldutil.setVectorSys(unfoldInputList['sig'], unfold_pt, "Pt", "unfoldsys", 1)
        unfoldutil.setVectorSys(unfoldInputList['sig'], unfold_mass, "Mass", "unfoldsys", 1)

        # check unfolded distribution
        outpdf = outputDirectory + "ratio_"+args.channel+".pdf"
        outpdf_mass = outputDirectory + "ratio_mass_"+args.channel+".pdf"
        drawutil.basicRatio(outpdf, unfold_pt, unfold_mass, unfoldInputList['sig'])
        drawutil.basicRatioMass(outpdf_mass, unfold_mass, unfoldInputList['sig'])

        outpdf = outputDirectory + "isrFit_"+args.channel+".pdf"
        drawutil.isrFit(outpdf, unfold_pt, unfold_mass, unfoldInputList['sig'])

	# check reco level distribution
	# FIXME use TUnfoldDensityV17::GetBackground()


	if not args.closure:
        	outpdf = outputDirectory + "recoPt_"+args.channel+".pdf"
        	if args.channel == "muon": drawutil.recoPt(outpdf, postfix, unfoldInputList['data'], unfoldInputList['sig'], unfoldInputList['DYtoEEtau'], unfoldInputList['TTbar'], unfoldInputList['VV'], unfoldInputList['Wjets'], unfoldInputList['QCD'], args.channel) #None for electron channel
		else : drawutil.recoPt(outpdf, postfix, unfoldInputList['data'], unfoldInputList['sig'], unfoldInputList['DYtoEEtau'], unfoldInputList['TTbar'], unfoldInputList['VV'], unfoldInputList['Wjets'], None, args.channel)

        	outpdf = outputDirectory + "recoMass_"+args.channel+".pdf"
        	if args.channel == "muon": drawutil.recoMass(outpdf, postfix, unfoldInputList['data'], unfoldInputList['sig'], unfoldInputList['DYtoEEtau'], unfoldInputList['TTbar'], unfoldInputList['VV'], unfoldInputList['Wjets'], unfoldInputList['QCD'], args.channel);
		else : drawutil.recoMass(outpdf, postfix, unfoldInputList['data'], unfoldInputList['sig'], unfoldInputList['DYtoEEtau'], unfoldInputList['TTbar'], unfoldInputList['VV'], unfoldInputList['Wjets'], None, args.channel)

        	outpdf = outputDirectory + "recoPtdist_"+args.channel+".pdf"
		drawutil.drawUnfoldedPt(outpdf, unfold_pt)

	if args.doSeperateUnfold:
		# FSR unfolding
		# get unfolded histogram from detector unfolding results
		hunfolded_pt = unfoldutil.getUnfoldedHist(unfold_pt, "detector")
		hunfolded_mass = unfoldutil.getUnfoldedHist(unfold_mass, "detector")

		# set migration matrix for FSR unfolding
        	unfoldFSR_pt = unfoldutil.setUnfold(unfoldInputList['sig'],   "Pt",   "FSR", True)
        	unfoldFSR_mass = unfoldutil.setUnfold(unfoldInputList['sig'], "Mass", "FSR", True)

        	outpdf = outputDirectory + "response_FSR_"+args.channel+".pdf"
        	drawutil.responseM(outpdf, unfoldFSR_pt)
        	outpdf = outputDirectory + "efficiency_FSR_"+args.channel+".pdf"
        	drawutil.efficiency(outpdf, unfoldFSR_pt)

        	outpdf_mass = outputDirectory + "response_mass_FSR_"+args.channel+".pdf"
        	drawutil.responseM(outpdf_mass, unfoldFSR_mass)
        	outpdf_mass = outputDirectory + "efficiency_mass_FSR_"+args.channel+".pdf"
        	drawutil.efficiency(outpdf_mass, unfoldFSR_mass)

		unfoldutil.setUnfoldInputHist(unfoldFSR_pt, hunfolded_pt)
		unfoldutil.setUnfoldInputHist(unfoldFSR_mass, hunfolded_mass)

        	unfoldutil.doUnfold(unfoldFSR_pt)
        	unfoldutil.doUnfold(unfoldFSR_mass)

        	outpdf = outputDirectory + "ratio_FSR_"+args.channel+".pdf"
        	outpdf_mass = outputDirectory + "ratio_mass_FSR_"+args.channel+".pdf"
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
