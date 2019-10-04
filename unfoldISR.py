import argparse
import os
import sys
import ROOT as rt

# TODO make seperate helper module python file
def setUnfoldBkgs(unfold_class, hfile_path, syst_name, isSys, nthSys, nTotSys):
    
    unfold_class.subBkgs("Pt", hfile_path, "DYJetsToTauTau",       isSys, syst_name, nTotSys, nthSys, "detector_level")
    unfold_class.subBkgs("Pt", hfile_path, "DYJets10to50ToTauTau", isSys, syst_name, nTotSys, nthSys, "detector_level")
    unfold_class.subBkgs("Pt", hfile_path, "TTLL_powheg",          isSys, syst_name, nTotSys, nthSys, "detector_level")
    unfold_class.subBkgs("Pt", hfile_path, "WW_pythia",            isSys, syst_name, nTotSys, nthSys, "detector_level")
    unfold_class.subBkgs("Pt", hfile_path, "WZ_pythia",            isSys, syst_name, nTotSys, nthSys, "detector_level")
    unfold_class.subBkgs("Pt", hfile_path, "ZZ_pythia",            isSys, syst_name, nTotSys, nthSys, "detector_level")
    unfold_class.subBkgs("Pt", hfile_path, "WJets_MG",             isSys, syst_name, nTotSys, nthSys, "detector_level")
    
    unfold_class.subBkgs("Mass", hfile_path, "DYJetsToTauTau",       isSys, syst_name, nTotSys, nthSys, "detector_level")
    unfold_class.subBkgs("Mass", hfile_path, "DYJets10to50ToTauTau", isSys, syst_name, nTotSys, nthSys, "detector_level")
    unfold_class.subBkgs("Mass", hfile_path, "TTLL_powheg",          isSys, syst_name, nTotSys, nthSys, "detector_level")
    unfold_class.subBkgs("Mass", hfile_path, "WW_pythia",            isSys, syst_name, nTotSys, nthSys, "detector_level")
    unfold_class.subBkgs("Mass", hfile_path, "WZ_pythia",            isSys, syst_name, nTotSys, nthSys, "detector_level")
    unfold_class.subBkgs("Mass", hfile_path, "ZZ_pythia",            isSys, syst_name, nTotSys, nthSys, "detector_level")
    unfold_class.subBkgs("Mass", hfile_path, "WJets_MG",             isSys, syst_name, nTotSys, nthSys, "detector_level")

parser = argparse.ArgumentParser(description='Unfolding for ISR analysis')

parser.add_argument('--channel' , dest = 'channel', default = 'electron', help = 'select channel electron or muon')
parser.add_argument('--year' , dest = 'year', default = '2016', help = 'select year')
parser.add_argument('--phase_space' , dest = 'phase_space', default = 'fiducial_phase_pre_FSR_dRp1', help = 'select unfolded phase space')
parser.add_argument('--FSR_dR' , dest = 'FSR_dR', default = 'dressed_dRp1', help = 'select size of dR for dressed lepton')
parser.add_argument('--createInputHists'  , action='store_true'  , help = 'create input histograms')
parser.add_argument('--getUnfoldResults'  , action='store_true'  , help = 'Get unfolding resutls')
parser.add_argument('--doSys'  , action='store_true'  , default = False, help = 'Calculate systematics')

# TODO options
parser.add_argument('--getCombinedResults'  , action='store_true'  , help = 'Combine 2016 and 2017')
parser.add_argument('--getEMuCombinedResults'  , action='store_true'  , help = 'Combine electron and muon')

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
                if modifiedLine.split()[1] == "matrix":
                        unfoldInputList['matrix'] = modifiedLine.split()[2]

                if modifiedLine.split()[1] == "fsr_matrix":
                        unfoldInputList['fsr_matrix'] = modifiedLine.split()[2]

                if modifiedLine.split()[1] == "hist":
                        unfoldInputList['hist'] = modifiedLine.split()[2]

        print unfoldInputList

        import pyScripts.unfoldUtil as unfoldutil
        import pyScripts.drawUtil as drawutil

        # create unfold class
        unfoldClass = rt.ISRUnfold(args.channel, unfoldInputList['hist'])

        # set response matrix
        unfoldClass.setNomTUnfoldDensity("Pt",  unfoldInputList['matrix'], args.phase_space, args.FSR_dR)
        unfoldClass.setNomTUnfoldDensity("Mass",unfoldInputList['matrix'], args.phase_space, args.FSR_dR)

        # for closure test
        unfoldClass.setSysTUnfoldDensity("Pt",   unfoldInputList['matrix'],  "Closure", -1, -1, args.phase_space, args.FSR_dR)
        unfoldClass.setSysTUnfoldDensity("Mass", unfoldInputList['matrix'],  "Closure", -1, -1, args.phase_space, args.FSR_dR)

        unfoldClass.setInput("Pt",   unfoldInputList['matrix'], True, "Closure", 0, 1., args.phase_space)
        unfoldClass.setInput("Mass", unfoldInputList['matrix'], True, "Closure", 0, 1., args.phase_space)

        # set unfolding input histogram
	unfoldClass.setInput("Pt",   unfoldInputList['hist'], False, "nominal", 0, 1., "detector_level")
	unfoldClass.setInput("Mass", unfoldInputList['hist'], False, "nominal", 0, 1., "detector_level")
        setUnfoldBkgs(unfoldClass, unfoldInputList['hist'], "nominal", False, 0, -1) 

        # set systematic response matrix and input histograms
	sysDict = {"PU": 2, "trgSF": 2, "recoSF": 2, "IdSF": 2, "L1Prefire": 2, "unfoldBias": 1, "unfoldScan": 1}
        if args.doSys == True:
	    #sysDict = {"PU": 2, "trgSF": 2, "recoSF": 2, "IdSF": 2, "IsoSF": 2, "unfoldsys": 1, "AlphaS": 2, "Scale": 9, "PDFerror": 100, "Alt": 1, "L1Prefire": 2, "LepScale": 2, "LepRes": 2, "FSRDR": 30, "unfoldBias": 1, "unfoldScan": 1}

	    for sysName, nSys in sysDict.items():
	    	for nthSys in range(0,nSys):
                    
                    print "sysName: " + sysName + " nthSys: " + str(nthSys) + " #####################################################"
                    
                    # set systematic response matrix
                    unfoldClass.setSysTUnfoldDensity("Pt",  unfoldInputList['matrix'],  sysName, nSys, nthSys, args.phase_space, args.FSR_dR)
                    unfoldClass.setSysTUnfoldDensity("Mass",unfoldInputList['matrix'],  sysName, nSys, nthSys, args.phase_space, args.FSR_dR)
                    
                    bias = 1.;
                    if sysName == "unfoldBias": bias = 0.95 
                    
                    # set systematic input histograms
                    unfoldClass.setInput("Pt",   unfoldInputList['hist'], True, sysName, nthSys, bias, "detector_level")
                    unfoldClass.setInput("Mass", unfoldInputList['hist'], True, sysName, nthSys, bias, "detector_level")
                    
                    # set systematic background histograms
                    setUnfoldBkgs(unfoldClass, unfoldInputList['hist'], sysName, True, nthSys, nSys)

	# unfold 
	unfoldClass.doISRUnfold(args.doSys)

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

        unfoldClass.setNomFSRTUnfoldDensity("Pt",    unfoldInputList['fsr_matrix'], "full_phase", "pre_fsr_dRp1")
        unfoldClass.setNomFSRTUnfoldDensity("Mass",  unfoldInputList['fsr_matrix'], "full_phase", "pre_fsr_dRp1")
        unfoldClass.setFSRUnfoldInput(unfoldInputList['fsr_matrix'], False, "", -1, "full_phase")

        if args.doSys == True:

            for sysName, nSys in sysDict.items():
                for nthSys in range(0,nSys):
                    unfoldClass.setSysFSRTUnfoldDensity("Pt",   unfoldInputList['fsr_matrix'], sysName, nSys, nthSys, "full_phase", "pre_fsr_dRp1") 
                    unfoldClass.setSysFSRTUnfoldDensity("Mass", unfoldInputList['fsr_matrix'], sysName, nSys, nthSys, "full_phase", "pre_fsr_dRp1") 

                    unfoldClass.setFSRUnfoldInput(unfoldInputList['fsr_matrix'], True, sysName, nthSys)
                    unfoldClass.setFSRUnfoldInput(unfoldInputList['fsr_matrix'], True, sysName, nthSys)


        unfoldClass.doISRQEDFSRUnfold(args.doSys)

        # set nominal value and also systematic values
        unfoldClass.setMeanPt(args.doSys, False, args.doSys)
        unfoldClass.setMeanMass(args.doSys, False, args.doSys)

        for massBin in range(0,5):

                if args.doSys:
                    unfoldClass.drawClosurePlots(outputDirectory + "Closure_"+args.channel, "Pt", massBin)
                    unfoldClass.drawClosurePlots(outputDirectory + "Closure_"+args.channel, "Mass", massBin)

                unfoldClass.drawNominalPlots(outputDirectory + "Unfolded_"+args.channel, "Pt", massBin)
                unfoldClass.drawNominalPlots(outputDirectory + "Unfolded_"+args.channel, "Mass", massBin)

        for massBin in range(0,5):

                unfoldClass.drawNominalPlots(outputDirectory + "Unfolded_"+args.channel, "Pt",   massBin, "", False, True)
                unfoldClass.drawNominalPlots(outputDirectory + "Unfolded_"+args.channel, "Mass", massBin, "", False, True)

        # draw plots including systematic 
        if args.doSys == True:
            for sysName, nSys in sysDict.items():

                if sysName == "Closure" : continue

                for massBin in range(0,5):

                    unfoldClass.drawNominalPlots(outputDirectory + "Unfolded_"+args.channel+sysName, "Pt", massBin, sysName, args.doSys)
                    unfoldClass.drawNominalPlots(outputDirectory + "Unfolded_"+args.channel+sysName, "Mass", massBin, sysName, args.doSys)
                    unfoldClass.drawInputPlots(outputDirectory + args.channel + sysName, "Pt", massBin, sysName)
                    unfoldClass.drawSysPlots(outputDirectory + "Sys_" + args.channel , massBin, sysName)


	unfoldClass.drawISRresult(outputDirectory + "ISRfit_", False, False)
        unfoldClass.drawISRMatrixInfo(outputDirectory, True)
        unfoldClass.drawISRMatrixInfo(outputDirectory, False)

	del unfoldClass

def makeRecoPlots():
        pass
# test 
if __name__ == "__main__": 
	makeRecoPlots()
