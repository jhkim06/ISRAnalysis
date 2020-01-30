import argparse
import os
import sys
import ROOT as rt

# TODO make seperate helper module python file
def setUnfoldBkgs(unfold_class, hfile_path, syst_name, isSys, nthSys, nTotSys, DYFake = False):

    if DYFake is False:
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
    else:
        unfold_class.subBkgs("Pt", hfile_path, "DYJets",       isSys, syst_name, nTotSys, nthSys, "detector_level_DY_Fake")
        unfold_class.subBkgs("Pt", hfile_path, "DYJets10to50", isSys, syst_name, nTotSys, nthSys, "detector_level_DY_Fake")
        unfold_class.subBkgs("Mass", hfile_path, "DYJets",       isSys, syst_name, nTotSys, nthSys, "detector_level_DY_Fake")
        unfold_class.subBkgs("Mass", hfile_path, "DYJets10to50", isSys, syst_name, nTotSys, nthSys, "detector_level_DY_Fake")

def doISRAnalysis(args, year, channel, doSys):

    # set output directory
    outputDirectory = 'output/'+year+'/' + channel + "/"
    inputfhisttxtName = outputDirectory + "fhist.txt"

    # make output directory
    if not os.path.exists( outputDirectory ):
        os.makedirs( outputDirectory )

    # read text file including root file path for unfolding
    fOutTxtCheck = open( inputfhisttxtName,'r')
    unfoldInputList = {}

    for line in fOutTxtCheck:

        modifiedLine = line.lstrip(' ').rstrip(' ').rstrip('\n')

        # root file for response matrix for detector unfolding
        if modifiedLine.split()[1] == "matrix":
            unfoldInputList['matrix'] = modifiedLine.split()[2]

        if modifiedLine.split()[1] == "closure_matrix":
            unfoldInputList['closure_matrix'] = modifiedLine.split()[2]

        if modifiedLine.split()[1] == "fsr_matrix":
            unfoldInputList['fsr_matrix'] = modifiedLine.split()[2]

        if modifiedLine.split()[1] == "fsr_photos_matrix":
            unfoldInputList['fsr_photos_matrix'] = modifiedLine.split()[2]

        if modifiedLine.split()[1] == "fsr_pythia_matrix":
            unfoldInputList['fsr_pythia_matrix'] = modifiedLine.split()[2]

        # root file for detector level histograms
        if modifiedLine.split()[1] == "hist":
            unfoldInputList['hist'] = modifiedLine.split()[2]

    print unfoldInputList

    import pyScripts.unfoldUtil as unfoldutil
    import pyScripts.drawUtil as drawutil

    # create unfold class
    unfoldClass = rt.ISRUnfold(channel, unfoldInputList['hist'], False, int(year))

    # set response matrix for detector unfolding
    unfoldClass.SetNomTUnfoldDensity("Pt",  unfoldInputList['matrix'], args.phase_space_detector, args.FSR_dR_detector)
    unfoldClass.SetNomTUnfoldDensity("Mass",unfoldInputList['matrix'], args.phase_space_detector, args.FSR_dR_detector)

    # for closure test
    if channel == "electron" and year == "2019":
        # closure test using othorgonal resonse matrix and detector histogram
        # currently available only for 2016
        unfoldClass.setSysTUnfoldDensity("Pt",   unfoldInputList['closure_matrix'],  "Closure", -1, -1, args.phase_space_detector+"_odd", args.FSR_dR_detector)
        unfoldClass.setSysTUnfoldDensity("Mass", unfoldInputList['closure_matrix'],  "Closure", -1, -1, args.phase_space_detector+"_odd", args.FSR_dR_detector)

        unfoldClass.setInput("Pt",   unfoldInputList['closure_matrix'], True, "Closure", 0, 1., args.phase_space_detector+"_even")
        unfoldClass.setInput("Mass", unfoldInputList['closure_matrix'], True, "Closure", 0, 1., args.phase_space_detector+"_even")
    else :
        unfoldClass.setSysTUnfoldDensity("Pt",   unfoldInputList['matrix'],  "Closure", -1, -1, args.phase_space_detector, args.FSR_dR_detector)
        unfoldClass.setSysTUnfoldDensity("Mass", unfoldInputList['matrix'],  "Closure", -1, -1, args.phase_space_detector, args.FSR_dR_detector)

        unfoldClass.setInput("Pt",   unfoldInputList['matrix'], True, "Closure", 0, 1., args.phase_space_detector)
        unfoldClass.setInput("Mass", unfoldInputList['matrix'], True, "Closure", 0, 1., args.phase_space_detector)

    # set unfolding input histogram
    unfoldClass.setInput("Pt",   unfoldInputList['hist'], False, "nominal", 0, 1., "detector_level")
    unfoldClass.setInput("Mass", unfoldInputList['hist'], False, "nominal", 0, 1., "detector_level")
    setUnfoldBkgs(unfoldClass, unfoldInputList['hist'], "nominal", False, 0, -1)

    # set systematic response matrix and input histograms
    if year == "2016":
        if channel == "electron" : sysDict = {"PU": 2, "trgSF": 2, "recoSF": 2, "IdSF": 2, "L1Prefire": 2, "AlphaS": 2, "Scale": 6, "PDFerror": 100, "unfoldBias": 1, "unfoldScan": 1}
        if channel == "muon" :     sysDict = {"PU": 2, "trgSF": 2, "IsoSF": 2, "IdSF": 2, "L1Prefire": 2, "AlphaS": 2, "Scale": 6, "PDFerror": 100, "unfoldBias": 1, "unfoldScan": 1}

    if year == "2017":
        if channel == "electron" : sysDict = {"PU": 2, "trgSF": 2, "recoSF": 2, "IdSF": 2, "L1Prefire": 2, "AlphaS": 2, "Scale": 6, "unfoldBias": 1, "unfoldScan": 1}
        if channel == "muon" :     sysDict = {"PU": 2, "trgSF": 2, "IsoSF": 2, "IdSF": 2, "L1Prefire": 2, "AlphaS": 2, "Scale": 6, "unfoldBias": 1, "unfoldScan": 1}

    if year == "2018":
        if channel == "electron" : sysDict = {"PU": 2, "recoSF": 2, "IdSF": 2, "L1Prefire": 2, "Scale": 6, "unfoldBias": 1, "unfoldScan": 1}
        if channel == "muon" :     sysDict = {"PU": 2, "IsoSF": 2, "IdSF": 2, "L1Prefire": 2, "Scale": 6, "unfoldBias": 1, "unfoldScan": 1}

    if doSys == True:

        for sysName, nSys in sysDict.items():
            for nthSys in range(0,nSys):

                #print "sysName: " + sysName + " nthSys: " + str(nthSys) + " #####################################################"
                # set systematic response matrix
                unfoldClass.setSysTUnfoldDensity("Pt",  unfoldInputList['matrix'],  sysName, nSys, nthSys, args.phase_space_detector, args.FSR_dR_detector)
                unfoldClass.setSysTUnfoldDensity("Mass",unfoldInputList['matrix'],  sysName, nSys, nthSys, args.phase_space_detector, args.FSR_dR_detector)

                bias = 1.
                if sysName == "unfoldBias": bias = 0.95

                # set systematic input histograms
                unfoldClass.setInput("Pt",   unfoldInputList['hist'], True, sysName, nthSys, bias, "detector_level")
                unfoldClass.setInput("Mass", unfoldInputList['hist'], True, sysName, nthSys, bias, "detector_level")

                # set systematic background histograms
                setUnfoldBkgs(unfoldClass, unfoldInputList['hist'], sysName, True, nthSys, nSys)

    # unfold for detector
    unfoldClass.doISRUnfold(doSys)

    if doSys == True:
        unfoldClass.drawLCurve(outputDirectory + "LCurve_" + channel + ".pdf", "Pt")
        unfoldClass.drawLCurve(outputDirectory + "LCurveMass_" + channel + ".pdf", "Mass")

        unfoldClass.drawRhoLog(outputDirectory + "RhoLog_" + channel + ".pdf", "Pt")
        unfoldClass.drawRhoLog(outputDirectory + "RhoLogMass_" + channel + ".pdf", "Mass")

    # set QED FSR unfolding response matrix and input
    unfoldClass.setNomFSRTUnfoldDensity("Pt",    unfoldInputList['fsr_matrix'], args.phase_space_fsr, args.FSR_dR_fsr)
    unfoldClass.setNomFSRTUnfoldDensity("Mass",  unfoldInputList['fsr_matrix'], args.phase_space_fsr, args.FSR_dR_fsr)
    unfoldClass.setFSRUnfoldInput(unfoldInputList['fsr_matrix'], False, "", -1, args.phase_space_fsr)

    if doSys == True:
        # systematic from detector unfolding
        for sysName, nSys in sysDict.items():
            for nthSys in range(0,nSys):
                unfoldClass.setSysFSRTUnfoldDensity("Pt",   unfoldInputList['fsr_matrix'], sysName, nSys, nthSys, args.phase_space_fsr, args.FSR_dR_fsr)
                unfoldClass.setSysFSRTUnfoldDensity("Mass", unfoldInputList['fsr_matrix'], sysName, nSys, nthSys, args.phase_space_fsr, args.FSR_dR_fsr)

                unfoldClass.setFSRUnfoldInput(unfoldInputList['fsr_matrix'], True, sysName, nthSys)

        # QED FSR ststematic
        sysDict["QED_FSR"] = 2

        unfoldClass.setSysFSRTUnfoldDensity("Pt",   unfoldInputList['fsr_photos_matrix'], "QED_FSR", 2, 0, args.phase_space_fsr, args.FSR_dR_fsr)
        unfoldClass.setSysFSRTUnfoldDensity("Mass", unfoldInputList['fsr_photos_matrix'], "QED_FSR", 2, 0, args.phase_space_fsr, args.FSR_dR_fsr)
        unfoldClass.setFSRUnfoldInput(unfoldInputList['fsr_matrix'], True, "QED_FSR", 0)

        unfoldClass.setSysFSRTUnfoldDensity("Pt",   unfoldInputList['fsr_pythia_matrix'], "QED_FSR", 2, 1, args.phase_space_fsr, args.FSR_dR_fsr)
        unfoldClass.setSysFSRTUnfoldDensity("Mass", unfoldInputList['fsr_pythia_matrix'], "QED_FSR", 2, 1, args.phase_space_fsr, args.FSR_dR_fsr)
        unfoldClass.setFSRUnfoldInput(unfoldInputList['fsr_matrix'], True, "QED_FSR", 1)

    # unfolding for QED FSR
    unfoldClass.doISRQEDFSRUnfold(doSys)

    # set nominal value and also systematic values
    unfoldClass.setMeanPt(doSys, False, doSys)
    unfoldClass.setMeanMass(doSys, False, doSys)
    #unfoldClass.setMCPreFSRMeanValues(unfoldInputList['fsr_matrix'])

    # now, draw plots
    # make output directory for closure test
    dirClosurePlots = "ClosurePlots/"

    for massBin in range(0,5):

        if doSys:
            if not os.path.exists( outputDirectory + dirClosurePlots ):
                os.makedirs( outputDirectory + dirClosurePlots )
            unfoldClass.drawClosurePlots(outputDirectory + dirClosurePlots + "Closure_"+channel, "Pt", massBin)
            unfoldClass.drawClosurePlots(outputDirectory + dirClosurePlots + "Closure_"+channel, "Mass", massBin)

        # detector unfolding
        unfoldClass.drawUnfoldedHists(outputDirectory + "Unfolded_"+channel, "Pt",   massBin)
        unfoldClass.drawUnfoldedHists(outputDirectory + "Unfolded_"+channel, "Mass", massBin)

        # FSR unfolding
        unfoldClass.drawUnfoldedHists(outputDirectory + "Unfolded_"+channel, "Pt",   massBin, "", False, True)
        unfoldClass.drawUnfoldedHists(outputDirectory + "Unfolded_"+channel, "Mass", massBin, "", False, True)

    # draw plots including systematic
    if doSys == True:
        # make output directory for systematic plots
        dirSysPlots = "SysPlots/"
        if not os.path.exists( outputDirectory + dirSysPlots ):
            os.makedirs( outputDirectory + dirSysPlots )

        for sysName, nSys in sysDict.items():

            if sysName is "Closure" : continue

            for massBin in range(0,5):
                unfoldClass.drawInputPlots(outputDirectory + channel + sysName, "Pt", massBin, sysName)

                unfoldClass.drawUnfoldedHists(outputDirectory + "Unfolded_"+channel+sysName, "Pt", massBin, sysName, doSys, True)
                unfoldClass.drawUnfoldedHists(outputDirectory + "Unfolded_"+channel+sysName, "Mass", massBin, sysName, doSys, True)

                if sysName is not "QED_FSR":

                    unfoldClass.drawSysPlots(outputDirectory + dirSysPlots + "Sys_" + channel , massBin, sysName, True)
                    unfoldClass.drawSysComparionPlots(outputDirectory + "Sys_" + channel , "Pt", massBin, sysName, True)
                    unfoldClass.drawSysComparionPlots(outputDirectory + "Sys_" + channel , "Mass", massBin, sysName, True)

                unfoldClass.drawSysPlots(outputDirectory + dirSysPlots + "Sys_" + channel , massBin, sysName, False)
                unfoldClass.drawSysComparionPlots(outputDirectory + "Sys_" + channel , "Pt", massBin, sysName, False)
                unfoldClass.drawSysComparionPlots(outputDirectory + "Sys_" + channel , "Mass", massBin, sysName, False)

        # draw with full systematic error
        for massBin in range(0,5):

            # draw all systematic points in one plot
            unfoldClass.drawSysPlots(outputDirectory + dirSysPlots + "Sys_" + channel , massBin, "full", True)
            unfoldClass.drawSysPlots(outputDirectory + dirSysPlots + "Sys_" + channel , massBin, "full", False)

            unfoldClass.drawUnfoldedHists(outputDirectory + "Unfolded_"+channel+"fullSys", "Pt", massBin, "full", doSys, True, True)
            unfoldClass.drawUnfoldedHists(outputDirectory + "Unfolded_"+channel+"fullSys", "Mass", massBin, "full", doSys, True, True)

    #unfoldClass.drawISRresult(outputDirectory + "ISRfit_", False, False)

    unfoldClass.drawISRMatrixInfo("Pt", outputDirectory, True)
    unfoldClass.drawISRMatrixInfo("Pt", outputDirectory, False)

    unfoldClass.drawISRMatrixInfo("Mass", outputDirectory, True)
    unfoldClass.drawISRMatrixInfo("Mass", outputDirectory, False)

    # FSR unfolding response matrix from Photos
    unfoldClass.drawISRMatrixInfo("Pt", outputDirectory, False, True)
    unfoldClass.drawISRMatrixInfo("Mass", outputDirectory, False, True)

    return unfoldClass

###
parser = argparse.ArgumentParser(description='Unfolding for ISR analysis')

parser.add_argument('--channel' , dest = 'channel', default = 'electron', help = 'select channel electron or muon')
parser.add_argument('--year' , dest = 'year', default = '2016', help = 'select year')
parser.add_argument('--phase_space_detector' , dest = 'phase_space_detector', default = 'fiducial_phase_pre_FSR_dRp1', help = 'select unfolded phase space')
parser.add_argument('--FSR_dR_detector' , dest = 'FSR_dR_detector', default = 'dressed_dRp1', help = 'select size of dR for dressed lepton')
parser.add_argument('--FSR_dR_fsr' , dest = 'FSR_dR_fsr', default = 'dRp1_pre_fsr', help = 'select size of dR for dressed lepton')
parser.add_argument('--phase_space_fsr' , dest = 'phase_space_fsr', default = 'full_phase_dRp1', help = 'select unfolded phase space')
parser.add_argument('--createInputHists'  , action='store_true'  , help = 'create input histograms')
parser.add_argument('--getUnfoldResults'  , action='store_true'  , help = 'Get unfolding resutls')
parser.add_argument('--doSys'  , action='store_true'  , default = False, help = 'Calculate systematics')

# TODO options
parser.add_argument('--doISRAnalysis'  , action='store_true'  , default = False, help = 'Use ISR analysis function ')
parser.add_argument('--getCombinedResults'  , action='store_true'  , help = 'Combine 2016 and 2017')
parser.add_argument('--getEMuCombinedResults'  , action='store_true'  , help = 'Combine electron and muon')

args = parser.parse_args()

# get input root files information using sampleDef.py
import etc.sampleDef as isrSamples
import pyScripts.unfoldInputUtil as histUtil

print "channel to run: " + args.channel

# for full Run2 result
if args.doISRAnalysis:

    # set output directory
    outputDirectory = 'output/Run2/'

    # make output directory
    if not os.path.exists( outputDirectory ):
        os.makedirs( outputDirectory )

    print "use ISR function"
    result_electron_2016 = doISRAnalysis(args, "2016", "electron", True)
    result_electron_2017 = doISRAnalysis(args, "2017", "electron", True)
    result_electron_2018 = doISRAnalysis(args, "2018", "electron", True)

    result_electron_2018.SavePtMassHists();
    result_electron_2017.SavePtMassHists();
    result_electron_2016.SavePtMassHists();

    result_muon_2016 = doISRAnalysis(args, "2016", "muon", True)
    result_muon_2017 = doISRAnalysis(args, "2017", "muon", True)
    result_muon_2018 = doISRAnalysis(args, "2018", "muon", True)

    result_electron_2016.drawISRresult(outputDirectory + "ISRfit_", False, False)
    canvas_electron_2017 = result_electron_2017.drawISRresult(outputDirectory + "ISRfit_", False, False)
    canvas_electron_2018 = result_electron_2018.drawISRresult(outputDirectory + "ISRfit_", False, False)

    canvas_muon_2016 = result_muon_2016.drawISRresult(outputDirectory + "ISRfit_", False, False)
    canvas_muon_2017 = result_muon_2017.drawISRresult(outputDirectory + "ISRfit_", False, False)
    canvas_muon_2018 = result_muon_2018.drawISRresult(outputDirectory + "ISRfit_", False, False)

    result_electron_2016.drawISRRun2results(outputDirectory + "ISRfit_", canvas_electron_2017, canvas_electron_2018, canvas_muon_2016, canvas_muon_2017, canvas_muon_2018)

if args.getUnfoldResults and args.doISRAnalysis == False:

    # set output directory
    outputDirectory = 'output/'+args.year+'/' + args.channel + "/"
    inputfhisttxtName = outputDirectory + "fhist.txt"

    # make output directory
    if not os.path.exists( outputDirectory ):
        os.makedirs( outputDirectory )

    # read text file including root file path for unfolding
    fOutTxtCheck = open( inputfhisttxtName,'r')
    unfoldInputList = {}

    for line in fOutTxtCheck:
            modifiedLine = line.lstrip(' ').rstrip(' ').rstrip('\n')
            if modifiedLine.split()[1] == "matrix":
                    unfoldInputList['matrix'] = modifiedLine.split()[2]

            if modifiedLine.split()[1] == "matrix_alt":
                    unfoldInputList['matrix_alt'] = modifiedLine.split()[2]

            if modifiedLine.split()[1] == "fsr_matrix":
                    unfoldInputList['fsr_matrix'] = modifiedLine.split()[2]

            if modifiedLine.split()[1] == "fsr_matrix_alt":
                    unfoldInputList['fsr_matrix_alt'] = modifiedLine.split()[2]

            if modifiedLine.split()[1] == "fsr_photos_matrix":
                    unfoldInputList['fsr_photos_matrix'] = modifiedLine.split()[2]

            if modifiedLine.split()[1] == "fsr_pythia_matrix":
                    unfoldInputList['fsr_pythia_matrix'] = modifiedLine.split()[2]

            if modifiedLine.split()[1] == "gen_hist":
                    unfoldInputList['gen_hist'] = modifiedLine.split()[2]

            if modifiedLine.split()[1] == "hist":
                    unfoldInputList['hist'] = modifiedLine.split()[2]

    print unfoldInputList

    import pyScripts.unfoldUtil as unfoldutil
    import pyScripts.drawUtil as drawutil

    DetectorUnfold = 0;
    FSRUnfold = 1;

    bias = 1.0

    # create unfold class                                                                    # regularization mode
    unfoldClass = rt.ISRUnfold(args.channel, unfoldInputList['hist'], unfoldInputList['matrix'], unfoldInputList['fsr_matrix'], False, int(args.year), int(0))
    unfoldClass.setOutputBaseDir(outputDirectory)
    unfoldClass.setBias(bias)
    # set response matrix
    unfoldClass.SetNomTUnfoldDensity("Pt",  unfoldInputList['matrix'], args.phase_space_detector, args.FSR_dR_detector)
    unfoldClass.SetNomTUnfoldDensity("Mass",unfoldInputList['matrix'], args.phase_space_detector, args.FSR_dR_detector)

    # FIXME
    unfoldClass.doClosureTest(DetectorUnfold, unfoldInputList['matrix'], "detector_level")

    # set unfolding input histogram
    unfoldClass.setInput("Pt",   unfoldInputList['hist'], False, "nominal", 0, bias, "detector_level")
    unfoldClass.setInput("Mass", unfoldInputList['hist'], False, "nominal", 0, bias, "detector_level")
    setUnfoldBkgs(unfoldClass, unfoldInputList['hist'], "nominal", False, 0, -1)
    setUnfoldBkgs(unfoldClass, unfoldInputList['matrix'], "nominal", False, 0, -1, True)

    # set systematic response matrix and input histograms
    #if args.channel == "electron" : sysDict = {"PU": 2, "trgSF": 2, "recoSF": 2, "IdSF": 2, "L1Prefire": 2, "AlphaS": 2, "Scale": 6, "PDFerror": 100, "unfoldBias": 1, "unfoldScan": 1, "Alt": 1}
    #if args.channel == "electron" : sysDict = {"PU": 2, "trgSF": 2, "recoSF": 2, "IdSF": 2, "L1Prefire": 2, "AlphaS": 2, "Scale": 6, "lepMom": 2, "unfoldBias": 1, "unfoldScan": 1, "Alt": 1}
    if args.channel == "electron" : sysDict = {"PU": 2, "trgSF": 2, "recoSF": 2, "IdSF": 2, "L1Prefire": 2, "Scale": 6, "unfoldBias": 1, "unfoldScan": 1, "Alt": 1}
    if args.channel == "muon" :     sysDict = {"PU": 2, "trgSF": 2, "IsoSF": 2, "IdSF": 2, "L1Prefire": 2, "AlphaS": 2, "Scale": 6, "lepMom": 2, "unfoldBias": 1, "unfoldScan": 1, "Alt": 1}
    if args.doSys == True:

        for sysName, nSys in sysDict.items():
            for nthSys in range(0,nSys):

                print "sysName: " + sysName + " nthSys: " + str(nthSys) + " #####################################################"

                if sysName == "Alt":
                    unfoldClass.setSysTUnfoldDensity("Pt",  unfoldInputList['matrix_alt'],  sysName, nSys, nthSys, args.phase_space_detector, args.FSR_dR_detector)
                    unfoldClass.setSysTUnfoldDensity("Mass",unfoldInputList['matrix_alt'],  sysName, nSys, nthSys, args.phase_space_detector, args.FSR_dR_detector)
                else :
                    # set systematic response matrix
                    unfoldClass.setSysTUnfoldDensity("Pt",  unfoldInputList['matrix'],  sysName, nSys, nthSys, args.phase_space_detector, args.FSR_dR_detector)
                    unfoldClass.setSysTUnfoldDensity("Mass",unfoldInputList['matrix'],  sysName, nSys, nthSys, args.phase_space_detector, args.FSR_dR_detector)

                if sysName == "unfoldBias": bias = 0.95

                # set systematic input histograms
                unfoldClass.setInput("Pt",   unfoldInputList['hist'], True, sysName, nthSys, bias, "detector_level")
                unfoldClass.setInput("Mass", unfoldInputList['hist'], True, sysName, nthSys, bias, "detector_level")

                # set systematic background histograms
                setUnfoldBkgs(unfoldClass, unfoldInputList['hist'], sysName, True, nthSys, nSys)
                setUnfoldBkgs(unfoldClass, unfoldInputList['matrix'], sysName, True, nthSys, nSys, True)


    # unfold
    unfoldClass.doISRUnfold(DetectorUnfold, args.doSys)

    if args.doSys == True:
        unfoldClass.drawLCurve(outputDirectory + "LCurve_" + args.channel + ".pdf", "Pt")
        unfoldClass.drawLCurve(outputDirectory + "LCurveMass_" + args.channel + ".pdf", "Mass")

        unfoldClass.drawRhoLog(outputDirectory + "RhoLog_" + args.channel + ".pdf", "Pt")
        unfoldClass.drawRhoLog(outputDirectory + "RhoLogMass_" + args.channel + ".pdf", "Mass")

    # set QED FSR unfolding response matrix and input
    unfoldClass.setNomFSRTUnfoldDensity("Pt",    unfoldInputList['fsr_matrix'], args.phase_space_fsr, args.FSR_dR_fsr)
    unfoldClass.setNomFSRTUnfoldDensity("Mass",  unfoldInputList['fsr_matrix'], args.phase_space_fsr, args.FSR_dR_fsr)
    unfoldClass.setFSRUnfoldInput(unfoldInputList['fsr_matrix'], False, "", -1, args.phase_space_fsr)

    unfoldClass.doClosureTest(FSRUnfold, unfoldInputList['fsr_matrix'], "fiducial_phase_dRp1")

    if args.doSys == True:
        # systematic from detector unfolding
        for sysName, nSys in sysDict.items():
            for nthSys in range(0,nSys):

                if sysName == "Alt":
                    unfoldClass.setSysFSRTUnfoldDensity("Pt",   unfoldInputList['fsr_matrix_alt'], sysName, nSys, nthSys, args.phase_space_fsr, args.FSR_dR_fsr)
                    unfoldClass.setSysFSRTUnfoldDensity("Mass", unfoldInputList['fsr_matrix_alt'], sysName, nSys, nthSys, args.phase_space_fsr, args.FSR_dR_fsr)

                else :
                    unfoldClass.setSysFSRTUnfoldDensity("Pt",   unfoldInputList['fsr_matrix'], sysName, nSys, nthSys, args.phase_space_fsr, args.FSR_dR_fsr)
                    unfoldClass.setSysFSRTUnfoldDensity("Mass", unfoldInputList['fsr_matrix'], sysName, nSys, nthSys, args.phase_space_fsr, args.FSR_dR_fsr)

                unfoldClass.setFSRUnfoldInput(unfoldInputList['fsr_matrix'], True, sysName, nthSys)

        # QED FSR ststematic
        sysDict["QED_FSR"] = 2

        # QED FSR ststematic
        unfoldClass.setSysFSRTUnfoldDensity("Pt",   unfoldInputList['fsr_photos_matrix'], "QED_FSR", 2, 0, args.phase_space_fsr, args.FSR_dR_fsr)
        unfoldClass.setSysFSRTUnfoldDensity("Mass", unfoldInputList['fsr_photos_matrix'], "QED_FSR", 2, 0, args.phase_space_fsr, args.FSR_dR_fsr)

        unfoldClass.setFSRUnfoldInput(unfoldInputList['fsr_matrix'], True, "QED_FSR", 0)

        unfoldClass.setSysFSRTUnfoldDensity("Pt",   unfoldInputList['fsr_pythia_matrix'], "QED_FSR", 2, 1, args.phase_space_fsr, args.FSR_dR_fsr)
        unfoldClass.setSysFSRTUnfoldDensity("Mass", unfoldInputList['fsr_pythia_matrix'], "QED_FSR", 2, 1, args.phase_space_fsr, args.FSR_dR_fsr)

        unfoldClass.setFSRUnfoldInput(unfoldInputList['fsr_matrix'], True, "QED_FSR", 1)

    # unfolding for QED FSR
    unfoldClass.doISRUnfold(FSRUnfold, args.doSys)

    # set nominal value and also systematic values
    unfoldClass.setMeanPt(args.doSys, False, args.doSys)
    unfoldClass.setMeanMass(args.doSys, False, args.doSys)

    # NOW DRAW PLOTS
    dirClosurePlots = "ClosurePlots/"
    if not os.path.exists( outputDirectory + dirClosurePlots ):
        os.makedirs( outputDirectory + dirClosurePlots )

    dirUnfoldHists = "UnfoldHist/"
    if not os.path.exists( outputDirectory + dirUnfoldHists ):
        os.makedirs( outputDirectory + dirUnfoldHists )

    for massBin in range(0,5):
        # closure test
        unfoldClass.drawClosurePlots(DetectorUnfold, unfoldInputList['gen_hist'], outputDirectory + dirClosurePlots + "DetClosure_" + args.channel, "Pt", massBin)
        unfoldClass.drawClosurePlots(DetectorUnfold, unfoldInputList['gen_hist'], outputDirectory + dirClosurePlots + "DetClosure_" + args.channel, "Mass", massBin)
        unfoldClass.drawClosurePlots(FSRUnfold, unfoldInputList['gen_hist'], outputDirectory + dirClosurePlots + "FSRClosure_" + args.channel, "Pt", massBin)
        unfoldClass.drawClosurePlots(FSRUnfold, unfoldInputList['gen_hist'], outputDirectory + dirClosurePlots + "FSRClosure_" + args.channel, "Mass", massBin)

        if not args.doSys:
            # detector unfold
            unfoldClass.drawUnfoldedHists(outputDirectory + dirUnfoldHists + "Unfolded_"+args.channel, "Pt", massBin)
            unfoldClass.drawUnfoldedHists(outputDirectory + dirUnfoldHists + "Unfolded_"+args.channel, "Mass", massBin)

            # FSR unfold
            unfoldClass.drawUnfoldedHists(outputDirectory + dirUnfoldHists + "Unfolded_"+args.channel, "Pt",   massBin, "", False, True)
            unfoldClass.drawUnfoldedHists(outputDirectory + dirUnfoldHists + "Unfolded_"+args.channel, "Mass", massBin, "", False, True)

    # draw plots including systematic
    if args.doSys == True:
        for sysName, nSys in sysDict.items():
            
            for massBin in range(0,5):
                
                unfoldClass.drawUnfoldedHists(outputDirectory + dirUnfoldHists + "Unfolded_"+args.channel+sysName, "Pt", massBin, sysName, args.doSys, True)
                unfoldClass.drawUnfoldedHists(outputDirectory + dirUnfoldHists + "Unfolded_"+args.channel+sysName, "Mass", massBin, sysName, args.doSys, True)
                #unfoldClass.drawInputPlots(outputDirectory + args.channel + sysName, "Pt", massBin, sysName)
                unfoldClass.drawSysPlots(outputDirectory + "Sys_" + args.channel , massBin, sysName, True)
                unfoldClass.drawSysPlots(outputDirectory + "Sys_" + args.channel , massBin, sysName, False)

                unfoldClass.drawSysPlots(outputDirectory + "Sys_" + args.channel , massBin, "full", True)
                unfoldClass.drawSysPlots(outputDirectory + "Sys_" + args.channel , massBin, "full", False)

                unfoldClass.drawUnfoldedHists(outputDirectory + "Unfolded_"+ args.channel+"fullSys", "Pt", massBin, "full", args.doSys, True, True)
                unfoldClass.drawUnfoldedHists(outputDirectory + "Unfolded_"+ args.channel+"fullSys", "Mass", massBin, "full", args.doSys, True, True)


    unfoldClass.drawISRresult(outputDirectory + "ISRfit_", False, False)

    unfoldClass.drawISRMatrixInfo("Pt", outputDirectory, True)
    unfoldClass.drawISRMatrixInfo("Pt", outputDirectory, False)

    unfoldClass.drawISRMatrixInfo("Mass", outputDirectory, True)
    unfoldClass.drawISRMatrixInfo("Mass", outputDirectory, False)

    #unfoldClass.drawISRMatrixInfo("Pt", outputDirectory, False, True)
    #unfoldClass.drawISRMatrixInfo("Mass", outputDirectory, False, True)

    del unfoldClass

def makeRecoPlots():
        pass
# test
if __name__ == "__main__":
	makeRecoPlots()
