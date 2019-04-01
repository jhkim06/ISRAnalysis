import os
import sys
import ROOT as rt

rt.gSystem.Load("/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/rootScripts/libunfold.so")
rt.gROOT.ProcessLine(".L /home/jhkim/ISR2016/unfolding/TUnfoldISR2016/rootScripts/makeRecoPlots.C++")

import gc
gc.collect()

import etc.histDef as fHistDef

# input file paths, binning, output directory, for data and bkg reco histograms only
def makeRecoHists(sample, outputDirectory, channel):

        rt.gROOT.SetBatch()
        print "####################### makeRecoHists #################################"

        outDic = {}

        # don't need to check isMC for reco histogram 
	outfile = rt.TFile(outputDirectory + sample.name+".root",'recreate')
        if not sample.isMC: outDic[sample.name] = fHistDef.inputfHists(sample.name, outputDirectory + sample.name+".root", "data")
        else : outDic[sample.name] = fHistDef.inputfHists(sample.name, outputDirectory + sample.name+".root", "bkg")

        # need to create histogram before the file loop to save one histogram from the input files
        recoHists = rt.recoHistsinfo(rt.std.map("TString, TH1*")())
        
        recoHists.histMaps.insert(rt.std.pair("const TString, TH1*")("pt_norminal", rt.ptHistogram("norminal")))
        recoHists.histMaps.insert(rt.std.pair("const TString, TH1*")("mass_norminal", rt.massHistogram("norminal")))

        recoHists.histMaps.insert(rt.std.pair("const TString, TH1*")("pt_noBveto", rt.ptHistogram("noBveto")))
        recoHists.histMaps.insert(rt.std.pair("const TString, TH1*")("mass_noBveto", rt.massHistogram("noBveto")))

        recoHists.histMaps.insert(rt.std.pair("const TString, TH1*")("pt_ZptWeight", rt.ptHistogram("ZptWeight")))
        recoHists.histMaps.insert(rt.std.pair("const TString, TH1*")("mass_ZptWeight", rt.massHistogram("ZptWeight")))

        for filepath in sample.path:	

		infile = rt.TFile(filepath,'r')
		print filepath
		rt.recoHists(infile, outfile, recoHists, channel)

                del infile

        outfile.Write()
        outfile.Delete()

        del recoHists
        return outDic

def makeSigHists(sample, outputDirectory, channel):
        rt.gROOT.SetBatch()
        print "####################### makeRecoHists #################################"

        outDic = {}

        # don't need to check isMC for reco histogram 
	outfile = rt.TFile(outputDirectory + sample.name+".root",'recreate')
        outfile_ = None

        outDic[sample.name] = fHistDef.inputfHists(sample.name, outputDirectory + sample.name+".root", "sig")
        # need to create histogram before the file loop to save one histogram from the input files
        sigHists  = rt.sigHistsinfo(rt.std.map("TString, TH1*")(), rt.std.map("TString, TH2*")(), sample.isInc) 
        recoHists = rt.recoHistsinfo(rt.std.map("TString, TH1*")())
       
	########################################### 
	rt.gInterpreter.GenerateDictionary("std::pair<TString, TH1*>", "map;TString.h;TH1.h")
	rt.gInterpreter.GenerateDictionary("std::map<TString, TH1*>", "map;TString.h;TH1.h")
	rt.gInterpreter.GenerateDictionary("std::pair<std::map<TString, TH1*>::iterator,bool>", "map;TString.h;TH1.h")

        rt.gInterpreter.GenerateDictionary("std::pair<TString, TH2*>", "map;TString.h;TH2.h")
        rt.gInterpreter.GenerateDictionary("std::map<TString, TH2*>", "map;TString.h;TH2.h")
        rt.gInterpreter.GenerateDictionary("std::pair<std::map<TString, TH2*>::iterator,bool>", "map;TString.h;TH2.h")
	#############################################

        sigHists.hist2DMaps.insert(rt.std.pair("const TString, TH2*")("pt_norminal", rt.ptMatrix("norminal")))
        sigHists.hist2DMaps.insert(rt.std.pair("const TString, TH2*")("mass_norminal", rt.massMatrix("norminal")))

        sigHists.hist2DMaps.insert(rt.std.pair("const TString, TH2*")("pt_noBveto", rt.ptMatrix("noBveto")))
        sigHists.hist2DMaps.insert(rt.std.pair("const TString, TH2*")("mass_noBveto", rt.massMatrix("noBveto")))

        sigHists.hist2DMaps.insert(rt.std.pair("const TString, TH2*")("pt_ZptWeight", rt.ptMatrix("ZptWeight")))
        sigHists.hist2DMaps.insert(rt.std.pair("const TString, TH2*")("mass_ZptWeight", rt.massMatrix("ZptWeight")))

	# AlphaS systematic
	for i in range(2):
        	sigHists.hist2DMaps.insert(rt.std.pair("const TString, TH2*")("pt_AlphaS_"+str(i), rt.ptMatrix("AlphaS_"+str(i))))
        	sigHists.hist2DMaps.insert(rt.std.pair("const TString, TH2*")("mass_AlphaS_"+str(i), rt.massMatrix("AlphaS_"+str(i))))

	# Scale systematic
	for i in range(9):
        	sigHists.hist2DMaps.insert(rt.std.pair("const TString, TH2*")("pt_Scale_"+str(i), rt.ptMatrix("Scale_"+str(i))))
        	sigHists.hist2DMaps.insert(rt.std.pair("const TString, TH2*")("mass_Scale_"+str(i), rt.massMatrix("Scale_"+str(i))))

        sigHists.histMaps.insert(rt.std.pair("const TString, TH1*")("pt_norminal", rt.ptHistogram("norminal")))
        sigHists.histMaps.insert(rt.std.pair("const TString, TH1*")("mass_norminal", rt.massHistogram("norminal")))

        sigHists.histMaps.insert(rt.std.pair("const TString, TH1*")("pt_noBveto", rt.ptHistogram("noBveto")))
        sigHists.histMaps.insert(rt.std.pair("const TString, TH1*")("mass_noBveto", rt.massHistogram("noBveto")))

        sigHists.histMaps.insert(rt.std.pair("const TString, TH1*")("pt_ZptWeight", rt.ptHistogram("ZptWeight")))
        sigHists.histMaps.insert(rt.std.pair("const TString, TH1*")("mass_ZptWeight", rt.massHistogram("ZptWeight")))

	# Migration matrix for detector and FSR seperately
        sigHists.hist2DMaps.insert(rt.std.pair("const TString, TH2*")("pt_detector", rt.ptMatrix("detector")))
        sigHists.hist2DMaps.insert(rt.std.pair("const TString, TH2*")("mass_detector", rt.massMatrix("detector")))

        sigHists.hist2DMaps.insert(rt.std.pair("const TString, TH2*")("pt_FSR", rt.ptFSRMatrix("FSR")))
        sigHists.hist2DMaps.insert(rt.std.pair("const TString, TH2*")("mass_FSR", rt.massFSRMatrix("FSR")))


        if sample.isInc: # for DY to tautau, make one more histogram
	        outfile_ = rt.TFile(outputDirectory + sample.name+"tau.root",'recreate')
        	outDic[sample.name+"tau"] = fHistDef.inputfHists(sample.name+"tau", outputDirectory + sample.name+"tau.root", "bkg")

		recoHists.histMaps.insert(rt.std.pair("const TString, TH1*")("pt_norminal", rt.ptHistogram("norminal")))
		recoHists.histMaps.insert(rt.std.pair("const TString, TH1*")("mass_norminal", rt.massHistogram("norminal")))

		recoHists.histMaps.insert(rt.std.pair("const TString, TH1*")("pt_noBveto", rt.ptHistogram("noBveto")))
		recoHists.histMaps.insert(rt.std.pair("const TString, TH1*")("mass_noBveto", rt.massHistogram("noBveto")))

		recoHists.histMaps.insert(rt.std.pair("const TString, TH1*")("pt_ZptWeight", rt.ptHistogram("ZptWeight")))
		recoHists.histMaps.insert(rt.std.pair("const TString, TH1*")("mass_ZptWeight", rt.massHistogram("ZptWeight")))

        for filepath in sample.path:	

		infile = rt.TFile(filepath,'r')
		print filepath
                temp_kfactor = 1.;
		#if filepath.find("M-10to50") > 0:
		#   print "this is M-10to50 DY sample"
                #   temp_kfactor = 6225.42/5765.4
		print "temp_kfactor: " + str(temp_kfactor)
		rt.sigHists(infile, outfile, outfile_, sigHists, recoHists, channel, temp_kfactor)

                del infile

        outfile.Write()
        outfile.Delete()

        if sample.isInc:
        	outfile_.Write()
        	outfile_.Delete()

        del sigHists
        del recoHists

	return outDic
def makeMigrationM():
	pass
