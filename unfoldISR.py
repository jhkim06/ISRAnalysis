import argparse
import os
import sys
import ROOT as rt

# TODO make seperate helper module python file
def setUnfoldBkgs(unfold_class, channel, hfile_path, postfix_for_histname, isSys, nthSys):

        unfold_class.subBkgs("Pt", postfix_for_histname, hfile_path, "DYJetsToTauTau",       nthSys, isSys, "detector_level")
        unfold_class.subBkgs("Pt", postfix_for_histname, hfile_path, "DYJets10to50ToTauTau", nthSys, isSys, "detector_level")
        unfold_class.subBkgs("Pt", postfix_for_histname, hfile_path, "TTLL_powheg",          nthSys, isSys, "detector_level")
        unfold_class.subBkgs("Pt", postfix_for_histname, hfile_path, "WW_pythia",            nthSys, isSys, "detector_level")
        unfold_class.subBkgs("Pt", postfix_for_histname, hfile_path, "WZ_pythia",            nthSys, isSys, "detector_level")
        unfold_class.subBkgs("Pt", postfix_for_histname, hfile_path, "ZZ_pythia",            nthSys, isSys, "detector_level")
        unfold_class.subBkgs("Pt", postfix_for_histname, hfile_path, "WJets_MG",             nthSys, isSys, "detector_level")

        unfold_class.subBkgs("Mass", postfix_for_histname, hfile_path, "DYJetsToTauTau",       nthSys, isSys, "detector_level_dipt100_cut")
        unfold_class.subBkgs("Mass", postfix_for_histname, hfile_path, "DYJets10to50ToTauTau", nthSys, isSys, "detector_level_dipt100_cut")
        unfold_class.subBkgs("Mass", postfix_for_histname, hfile_path, "TTLL_powheg",          nthSys, isSys, "detector_level_dipt100_cut")
        unfold_class.subBkgs("Mass", postfix_for_histname, hfile_path, "WW_pythia",            nthSys, isSys, "detector_level_dipt100_cut")
        unfold_class.subBkgs("Mass", postfix_for_histname, hfile_path, "WZ_pythia",            nthSys, isSys, "detector_level_dipt100_cut")
        unfold_class.subBkgs("Mass", postfix_for_histname, hfile_path, "ZZ_pythia",            nthSys, isSys, "detector_level_dipt100_cut")
        unfold_class.subBkgs("Mass", postfix_for_histname, hfile_path, "WJets_MG",             nthSys, isSys, "detector_level_dipt100_cut")

parser = argparse.ArgumentParser(description='Unfolding for ISR analysis')

parser.add_argument('--channel' , dest = 'channel', default = 'electron', help = 'select channel electron or muon')
parser.add_argument('--year' , dest = 'year', default = '2016', help = 'select year')
parser.add_argument('--postfix' , dest = 'postfix', default = 'nominal', help = 'select histogram name')
parser.add_argument('--createInputHists'  , action='store_true'  , help = 'create input histograms')
parser.add_argument('--altResponse'  , action='store_true'  , default = False, help = 'create response matrix with alternative signal MC')
parser.add_argument('--getUnfoldResults'  , action='store_true'  , help = 'Get unfolding resutls')
parser.add_argument('--doSys'  , action='store_true'  , default = False, help = 'Calculate systematics')
parser.add_argument('--closure'  , action='store_true'  , help = 'Clousre test with MC')
parser.add_argument('--getCombinedResults'  , action='store_true'  , help = 'Combine 2016 and 2017')
parser.add_argument('--getEMuCombinedResults'  , action='store_true'  , help = 'Combine electron and muon')
parser.add_argument('--doSeperateUnfold'  , action='store_true'  , default = False, help = 'Seperate unfolding steps')

args = parser.parse_args()

# get input root files information using sampleDef.py
import etc.sampleDef as isrSamples
import pyScripts.unfoldInputUtil as histUtil

# set output directory
outputDirectory = 'output/'+args.year+'/' + args.channel + "/"  
inputfhisttxtName = outputDirectory + "fhist.txt"

# make output directory
if not os.path.exists( outputDirectory ):
	os.makedirs( outputDirectory )

print "channel to run: " + args.channel

if args.getUnfoldResults:
        # read text file including root file path for unfolding
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

                if modifiedLine.split()[1] == "hist":
                        unfoldInputList['hist'] = modifiedLine.split()[2]

                if modifiedLine.split()[1] == "bkg": # use the sample name as keyword for background
                        unfoldInputList[modifiedLine.split()[0]] = modifiedLine.split()[2]

        print unfoldInputList

        import pyScripts.unfoldUtil as unfoldutil
        import pyScripts.drawUtil as drawutil

        postfix = args.postfix

        # create unfold class
        unfoldClass = rt.ISRUnfold()

        # set response matrix
        unfoldClass.setNomTUnfoldDensity("Pt",  unfoldInputList['matrix'], "full_phase", "dressed_dR10")
        unfoldClass.setNomTUnfoldDensity("Mass",unfoldInputList['matrix'], "full_phase", "dressed_dR10")

        # for closure test
        unfoldClass.setSysTUnfoldDensity("Pt",   unfoldInputList['matrix'],  "Closure", 0, "full_phase", "dressed_dR10")
        unfoldClass.setSysTUnfoldDensity("Mass", unfoldInputList['matrix'],  "Closure", 0, "full_phase", "dressed_dR10")

        unfoldClass.setInput(args.channel, "Pt",   "Closure", unfoldInputList['matrix'], 0, True, 1., "full_phase")
        unfoldClass.setInput(args.channel, "Mass", "Closure", unfoldInputList['matrix'], 0, True, 1., "full_phase")

        # set unfolding input histogram
	unfoldClass.setInput(args.channel, "Pt",   postfix, unfoldInputList['hist'], 0, False, 1., "detector_level")
	unfoldClass.setInput(args.channel, "Mass", postfix, unfoldInputList['hist'], 0, False, 1., "detector_level_dipt100_cut")
        setUnfoldBkgs(unfoldClass, args.channel, unfoldInputList['hist'], postfix, False, 0) 

        # temp way to compare post/pre FSR unfolded results
        unfoldClass.setSysTUnfoldDensity("Pt",   unfoldInputList['matrix'],  "detector_temp", 0, "fiducial_phase_post_FSR", "post_fsr")
        unfoldClass.setSysTUnfoldDensity("Mass", unfoldInputList['matrix'],  "detector_temp", 0, "fiducial_phase_post_FSR", "post_fsr")
        unfoldClass.setInput(args.channel, "Pt",   "detector_temp", unfoldInputList['hist'], 0, True, 1., "detector_level")
        unfoldClass.setInput(args.channel, "Mass", "detector_temp", unfoldInputList['hist'], 0, True, 1., "detector_level_dipt100_cut")
        setUnfoldBkgs(unfoldClass, args.channel, unfoldInputList['hist'], "detector_temp", True, 0) 


        # set systematic response matrix and input histograms
        if args.doSys == True:
	    sysDict = {"PU": 2, "trgSF": 2, "recoSF": 2, "IdSF": 2, "IsoSF": 2, "unfoldsys": 1, "AlphaS": 2, "Scale": 9, "PDFerror": 100, "Alt": 1, "L1Prefire": 2, "LepScale": 2, "LepRes": 2, "FSRDR": 30, "unfoldBias": 1, "unfoldScan": 1}

	    for sysName, nSys in sysDict.items():
	    	for nthSys in range(0,nSys):

	    		print "sysName: " + sysName + " nthSys: " + str(nthSys) + " #####################################################"

	    		# set systematic response matrix
	    		postfix = sysName
            		if sysName != "Alt": 
                                # use systematic response matrix saved in the signal root file 
	    			unfoldClass.setSysTUnfoldDensity("Pt",  unfoldInputList['sig'],    postfix, nthSys)
            			unfoldClass.setSysTUnfoldDensity("Mass",unfoldInputList['sig'],    postfix, nthSys)
	    		if sysName == "Alt":
	    			unfoldClass.setSysTUnfoldDensity("Pt",   unfoldInputList['sigAlt'], postfix, nthSys)
            			unfoldClass.setSysTUnfoldDensity("Mass", unfoldInputList['sigAlt'], postfix, nthSys)

	    		bias = 1.;
	    		if sysName == "unfoldBias": bias = 0.95 
	    	
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

	# unfold 
	unfoldClass.doISRUnfold(True)

	# set nominal value and also systematic values
	unfoldClass.setMeanPt(args.doSys, False, True)
	unfoldClass.setMeanMass(args.doSys, False, True)

        for massBin in range(0,5):
            
                unfoldClass.drawClosurePlots(outputDirectory + "Closure_"+args.channel, "Pt", massBin) 
                unfoldClass.drawClosurePlots(outputDirectory + "Closure_"+args.channel, "Mass", massBin)

                unfoldClass.drawNominalPlots(outputDirectory + "Unfolded_"+args.channel, "Pt", massBin)
                unfoldClass.drawNominalPlots(outputDirectory + "Unfolded_"+args.channel, "Mass", massBin)



        # draw plots including systematic 
        if args.doSys == True:
            for sysName, nSys in sysDict.items():

                if sysName == "Closure" : continue

	        for massBin in range(0,5):

                    unfoldClass.drawNominalPlots(outputDirectory + "Unfolded_"+args.channel+sysName, "Pt", massBin, sysName, args.doSys)
                    unfoldClass.drawNominalPlots(outputDirectory + "Unfolded_"+args.channel+sysName, "Mass", massBin, sysName, args.doSys)
                    unfoldClass.drawInputPlots(outputDirectory + args.channel + sysName, "Pt", massBin, sysName)
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

        #unfoldClass.drawLCurve(outputDirectory + "LCurve_" + args.channel + ".pdf", "Pt")
        #unfoldClass.drawLCurve(outputDirectory + "LCurveMass_" + args.channel + ".pdf", "Mass")

        #unfoldClass.drawRhoLog(outputDirectory + "RhoLog_" + args.channel + ".pdf", "Pt")
        #unfoldClass.drawRhoLog(outputDirectory + "RhoLogMass_" + args.channel + ".pdf", "Mass")


	unfoldClass.drawISRresult(outputDirectory + "ISRfit_", args.channel, False)

	del unfoldClass

def setUnfolInput(unfold_class, channel, hfile_path, postfix_for_histname, isSys, nthSys, hist_dir):

        unfold_class.setInput(channel, "Pt",   postfix_for_histname, hfile_path, nthSys, isSys, 1., "detector_level")
        unfold_class.setInput(channel, "Mass", postfix_for_histname, hfile_path, nthSys, isSys, 1., "detector_level_dipt100_cut")


def makeRecoPlots():
        pass
# test 
if __name__ == "__main__": 
	makeRecoPlots()
