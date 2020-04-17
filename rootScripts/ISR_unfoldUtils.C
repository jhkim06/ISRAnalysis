#include "ISR_unfoldUtils.h"
#include "tdrstyle.C"
#include "CMS_lumi.C"

void ISRUnfold::setBias(double bias)
{
   nominal_bias = bias;
}

void ISRUnfold::setOutputBaseDir(TString outPath)
{
   output_baseDir = outPath;
}

void ISRUnfold::setNomResMatrix(TString var, TString filepath, TString matrixName)
{
    TFile* filein = new TFile(filepath);

    TString Rec_binName = "Rec_"+var;
    TString Gen_binName = "Gen_"+var;
    Rec_binName = matrixName + "/" + var + "_ResMatrix_" + matrixName + "/" + Rec_binName;
    Gen_binName = matrixName + "/" + var + "_ResMatrix_" + matrixName + "/" + Gen_binName;

    if(var == "Pt")
    {
        pt_binning_Rec = (TUnfoldBinning*)filein->Get(Rec_binName);
        pt_binning_Gen = (TUnfoldBinning*)filein->Get(Gen_binName);
    }
    else if(var == "Mass")
    {
        mass_binning_Rec = (TUnfoldBinning*)filein->Get(Rec_binName);
        mass_binning_Gen = (TUnfoldBinning*)filein->Get(Gen_binName);
    }
    else
    {
        cout << "ISRUnfold::setNomResMatrix, only Pt and Mass available for var" << endl;
        exit (EXIT_FAILURE);
    }


    // Set response matrix
    TH2* hmcGenRec;
    hmcGenRec = (TH2*)filein->Get(matrixName + "/" + var + "_ResMatrix_" + matrixName +"/hmc" + var + "GenRecnominal");
    
    
    if( var == "Pt" )
    {
    	nomPtUnfold = new TUnfoldDensityV17(hmcGenRec,
    	                               TUnfold::kHistMapOutputHoriz,
    	                               regMode_detector,
    	                               TUnfold::kEConstraintArea,
    	                               TUnfoldDensityV17::kDensityModeBinWidth,
    	                               pt_binning_Gen,pt_binning_Rec);
    
    }
    else 
    {
        nomMassUnfold = new TUnfoldDensityV17(hmcGenRec,
                                        TUnfold::kHistMapOutputHoriz,
                                        regMode_detector,
                                        TUnfold::kEConstraintArea,
                                        TUnfoldDensityV17::kDensityModeBinWidth,
                                        mass_binning_Gen,mass_binning_Rec);
    }

    filein->Close();
    delete filein;

}

void ISRUnfold::setNomFSRResMatrix(TString var, TString filepath, TString migrationName, TString phaseSpace)
{
    
    TFile* filein = new TFile(filepath);

    // Set response matrix
    TH2* hmcGenGen;
    hmcGenGen = (TH2*)filein->Get(migrationName + "_" + phaseSpace + "/" + var + "_ResMatrix_" + migrationName + "/hmc" + var + "GenRecnominal");

    if( var == "Pt" )
    {
            nomPtFSRUnfold = new TUnfoldDensityV17(hmcGenGen,
                                           TUnfold::kHistMapOutputHoriz,
                                           regMode_FSR, 
                                           TUnfold::kEConstraintArea,
                                           TUnfoldDensityV17::kDensityModeBinWidth,
                                           pt_binning_Gen,pt_binning_Gen);
    }
    else if( var == "Mass" )
    {
            nomMassFSRUnfold = new TUnfoldDensityV17(hmcGenGen,
                                           TUnfold::kHistMapOutputHoriz,
                                           regMode_FSR,
                                           TUnfold::kEConstraintArea,
                                           TUnfoldDensityV17::kDensityModeBinWidth,
                                           mass_binning_Gen,mass_binning_Gen);
    }
    else
    {
        cout << "ISRUnfold::setNomResMatrix, only Pt and Mass available for var" << endl;
        exit (EXIT_FAILURE);
    }

    filein->Close();
    delete filein;
}

void ISRUnfold::setSysFSRTUnfoldDensity(TString var, TString filepath, TString sysName, int totSysN, int nth, TString phase_name, TString fsr_correction_name)
{

    TFile* filein = new TFile(filepath);

    TString systematic_postfix = sysName;

    // set proper postfix to read systematic histograms
    if(totSysN == 2)
    {
        if(nth == 0 ) systematic_postfix+="Up";
        if(nth == 1 ) systematic_postfix+="Down";
    }

    if(totSysN == 6 && (sysName == "Scale" || sysName == "pdfScale"))
    {
        if(nth == 0 ) systematic_postfix+="AUp";
        if(nth == 1 ) systematic_postfix+="ADown";
        if(nth == 2 ) systematic_postfix+="BUp";
        if(nth == 3 ) systematic_postfix+="BDown";
        if(nth == 4 ) systematic_postfix+="ABUp";
        if(nth == 5 ) systematic_postfix+="ABDown";
    }

    if(totSysN == 100 && sysName == "PDFerror")
    {
        TString nth_;
        nth_.Form ("%03d", nth);
        systematic_postfix+=nth_;
    }

    // set response matrix
    TH2* hmcGenGen;

    if( (sysName == "AlphaS" || sysName == "alphaS") || (sysName == "Scale" || sysName == "pdfScale") || sysName == "PDFerror")
    {

        if(var == "Pt")
            hmcGenGen = (TH2*)filein->Get(phase_name + "/ptll_gen_post_fsr_" + fsr_correction_name + "_response_matrix/hmc" + var + "GenRec_" + systematic_postfix);
        else if(var == "Mass")
            hmcGenGen = (TH2*)filein->Get(phase_name + "/mll_gen_post_fsr_" + fsr_correction_name + "_response_matrix/hmc" + var + "GenRec_" + systematic_postfix);
        else
        {
            cout << "ISRUnfold::setNomResMatrix, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
        }

    }
    else
    {
        // use nominal response matrix
        if(var == "Pt")
            hmcGenGen = (TH2*)filein->Get(phase_name + "/ptll_gen_post_fsr_" + fsr_correction_name + "_response_matrix/hmc" + var + "GenRecnominal");
        else if(var == "Mass")
            hmcGenGen = (TH2*)filein->Get(phase_name + "/mll_gen_post_fsr_" + fsr_correction_name + "_response_matrix/hmc" + var + "GenRecnominal");
        else
        {
            cout << "ISRUnfold::setNomResMatrix, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
        }
    }

    TUnfold::ERegMode mode = regMode_FSR;
    if( sysName =="unfoldScan" || sysName=="unfoldBias")
    {
        mode = TUnfold::kRegModeCurvature;
        //mode = TUnfold::kRegModeMixed;
    }


    if( var == "Pt" )
    {
            sysPtFSRUnfold[sysName].push_back( new TUnfoldDensityV17(hmcGenGen,
                                           TUnfold::kHistMapOutputHoriz,
                                           mode,
                                           TUnfold::kEConstraintArea,
                                           TUnfoldDensityV17::kDensityModeBinWidth,
                                           pt_binning_Gen,pt_binning_Gen));
    }
    else if( var == "Mass" )
    {
            sysMassFSRUnfold[sysName].push_back( new TUnfoldDensityV17(hmcGenGen,
                                           TUnfold::kHistMapOutputHoriz,
                                           mode,
                                           TUnfold::kEConstraintArea,
                                           TUnfoldDensityV17::kDensityModeBinWidth,
                                           mass_binning_Gen,mass_binning_Gen));
    }
    else
    {
        cout << "ISRUnfold::setNomResMatrix, only Pt and Mass available for var" << endl;
        exit (EXIT_FAILURE);
    }

    filein->Close();
    delete filein;
}

void ISRUnfold::setSysTUnfoldDensity(TString var, TString filepath, TString sysName, int totSysN, int nth, TString phase_name, TString fsr_correction_name)
{

    TFile* filein = new TFile(filepath);

    TString systematic_postfix = sysName;

    if(totSysN == 2)
    {
        if(nth == 0 ) systematic_postfix = "_" + systematic_postfix + "Up";
        if(nth == 1 ) systematic_postfix = "_" + systematic_postfix + "Down";
    }

    if(totSysN == 6 && (sysName =="Scale" || sysName =="pdfScale"))
    {
        if(nth == 0 ) systematic_postfix = "_" + systematic_postfix + "AUp";
        if(nth == 1 ) systematic_postfix = "_" + systematic_postfix + "ADown";
        if(nth == 2 ) systematic_postfix = "_" + systematic_postfix + "BUp";
        if(nth == 3 ) systematic_postfix = "_" + systematic_postfix + "BDown";
        if(nth == 4 ) systematic_postfix = "_" + systematic_postfix + "ABUp";
        if(nth == 5 ) systematic_postfix = "_" + systematic_postfix + "ABDown";
    }

    if(totSysN == 100 && sysName == "PDFerror")
    {
        TString nth_;
        nth_.Form ("%03d", nth);
        systematic_postfix = "_" + systematic_postfix + nth_;
    }

    // Set migration matrix
    TH2* hmcGenRec = NULL;
    // systematic using the same response matrix
    // may be better to use flag to use nominal response matrix or not
    if(sysName=="Alt" || sysName=="unfoldBias" || sysName=="unfoldScan" || sysName == "Stat")
    {
        if(var == "Pt")   hmcGenRec = (TH2*)filein->Get(phase_name + "/ptll_rec_gen_" + fsr_correction_name + "_response_matrix/hmc" + var + "GenRecnominal");
        else if(var == "Mass") hmcGenRec = (TH2*)filein->Get(phase_name + "/mll_rec_gen_" + fsr_correction_name + "_response_matrix/hmc" + var + "GenRecnominal");
        else
        {
            cout << "ISRUnfold::setSysTUnfoldDensity, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
        }
    }
    else
    {

        TString histDirPostfix = "";
        if(sysName == "lepMom"){

            if(nth == 0)
            {
                phase_name += "_lepMomUp";
                histDirPostfix = "_lepMomUp";
                systematic_postfix = "nominal";
            }
            else if(nth == 1)
            {
                phase_name += "_lepMomDown";
                histDirPostfix = "_lepMomDown";
                systematic_postfix = "nominal";
            }
            else
            {
                exit(EXIT_FAILURE);
            }
        }

        if(var == "Pt")
        {
            hmcGenRec = (TH2*)filein->Get(phase_name + "/ptll_rec_gen_" + fsr_correction_name + "_response_matrix" + histDirPostfix + "/hmc" + var + "GenRec" + systematic_postfix);
        }
        else if(var == "Mass")
        {
            hmcGenRec = (TH2*)filein->Get(phase_name + "/mll_rec_gen_" + fsr_correction_name + "_response_matrix" + histDirPostfix + "/hmc" + var + "GenRec" + systematic_postfix);
        }
        else
        {
            cout << "ISRUnfold::setSysTUnfoldDensity, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
        }
    }

    TUnfold::ERegMode mode = regMode_detector;
    if( sysName =="unfoldScan" || sysName=="unfoldBias")
    {
        mode = TUnfold::kRegModeCurvature;
        //mode = TUnfold::kRegModeMixed;
    }

    if( var == "Pt" )
    {
            sysPtUnfold[sysName].push_back( new TUnfoldDensityV17(hmcGenRec,
                                           TUnfold::kHistMapOutputHoriz,
                                           mode,
                                           TUnfold::kEConstraintArea,
                                           TUnfoldDensityV17::kDensityModeBinWidth,
                                           pt_binning_Gen,pt_binning_Rec));
    }

    else if( var == "Mass" )
    {
            sysMassUnfold[sysName].push_back( new TUnfoldDensityV17(hmcGenRec,
                                           TUnfold::kHistMapOutputHoriz,
                                           mode,
                                           TUnfold::kEConstraintArea,
                                           TUnfoldDensityV17::kDensityModeBinWidth,
                                           mass_binning_Gen,mass_binning_Rec));
    }
    else
    {
        cout << "ISRUnfold::setSysTUnfoldDensity, only Pt and Mass available for var" << endl;
        exit (EXIT_FAILURE);
    }

    filein->Close();
    delete filein;
}

void ISRUnfold::setFSRUnfInput(bool isSys, TString sysName, int nth)
{

    if(!isSys)
    {
        nomPtFSRUnfold->SetInput(nomPtUnfold->GetOutput("hpreFSR_pt", 0, 0, 0, false), 1.);
        nomMassFSRUnfold->SetInput(nomMassUnfold->GetOutput("hpreFSR_mass", 0, 0, 0, false), 1.);
    }
    else
    {
        if(sysName=="QED_FSR")
        {
            sysPtFSRUnfold[sysName].at(nth)  ->SetInput(nomPtUnfold->GetOutput("hpreFSR_pt", 0, 0, 0, false),   1.);
            sysMassFSRUnfold[sysName].at(nth)->SetInput(nomMassUnfold->GetOutput("hpreFSR_mass", 0, 0, 0, false),   1.);
        }
        else{
            sysPtFSRUnfold[sysName].at(nth)  ->SetInput(sysPtUnfold[sysName].at(nth)->GetOutput("hpreFSR_pt", 0, 0, 0, false),   1.);
            sysMassFSRUnfold[sysName].at(nth)->SetInput(sysMassUnfold[sysName].at(nth)->GetOutput("hpreFSR_mass", 0, 0, 0, false),   1.);
        }
    }
    
}


void ISRUnfold::setUnfInput(TString var, TString filepath, bool isSys, TString sysName, int nth, double bias, TString phase_name)
{

    TFile* filein = new TFile(filepath);
    TString nth_;
    nth_.Form("%d", nth);
    TH1* hRec = NULL;

    // nominal histograms
    // TODO use input covariance matrix
    if(!isSys)
    {
        if(var == "Pt")
        {
            if(channel_name == "muon")     hRec = (TH1*)filein->Get(phase_name + "/"+var+"/histo_DoubleMuonnominal");
            if(channel_name == "electron" && year != 2018) hRec = (TH1*)filein->Get(phase_name + "/hist_ptll/histo_DoubleEGnominal");
            if(channel_name == "electron" && year == 2018) hRec = (TH1*)filein->Get(phase_name + "/hist_ptll/histo_EGammanominal");

            nomPtUnfold->SetInput(hRec,   bias);
        }
        else if(var == "Mass")
        {
            if(channel_name == "muon")     hRec = (TH1*)filein->Get(phase_name + "/"+var+"/histo_DoubleMuonnominal");
            if(channel_name == "electron" && year != 2018) hRec = (TH1*)filein->Get(phase_name + "/hist_mll/histo_DoubleEGnominal");
            if(channel_name == "electron" && year == 2018) hRec = (TH1*)filein->Get(phase_name + "/hist_mll/histo_EGammanominal");
            nomMassUnfold->SetInput(hRec, bias);
        }
        else{
            cout << "ISRUnfold::setUnfInput, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
        }
    }
    // systematic histograms
    else
    {
        // for systematic, using the same input histogram as nominal, unless data changed in a systematic change
        if(sysName != "lepMom")
        {

            if(sysName != "Stat")
            {
                if(var == "Pt")
                {
                    if(channel_name == "muon")     hRec = (TH1*)filein->Get(phase_name + "/hist_ptll/histo_DoubleMuonnominal");
                    if(channel_name == "electron" && year != 2018) hRec = (TH1*)filein->Get(phase_name + "/hist_ptll/histo_DoubleEGnominal");
                    if(channel_name == "electron" && year == 2018) hRec = (TH1*)filein->Get(phase_name + "/hist_ptll/histo_EGammanominal");
                    sysPtUnfold[sysName].at(nth)  ->SetInput(hRec,   bias);
                }
                else if(var == "Mass")
                {
                    if(channel_name == "muon")     hRec = (TH1*)filein->Get(phase_name + "/hist_mll/histo_DoubleMuonnominal");
                    if(channel_name == "electron" && year != 2018) hRec = (TH1*)filein->Get(phase_name + "/hist_mll/histo_DoubleEGnominal");
                    if(channel_name == "electron" && year == 2018) hRec = (TH1*)filein->Get(phase_name + "/hist_mll/histo_EGammanominal");
                    sysMassUnfold[sysName].at(nth)->SetInput(hRec,   bias);
                }
                else
                {
                    cout << "ISRUnfold::setUnfInput, only Pt and Mass available for var" << endl;
                    exit (EXIT_FAILURE);
                }
            }
            else
            {
                if(var == "Pt")
                {
                    TH1* temp_ptHist;

                    TString nth_;
                    nth_.Form("%d", nth);
                    temp_ptHist = nomPtUnfold->GetInput("ptToy_" + nth_, 0, 0, 0, false);

                    // randomize histogram bins
                    for(int ibin = 1; ibin<temp_ptHist->GetNbinsX()+1;ibin++)
                    {
                        double err = temp_ptHist->GetBinError(ibin);
                        if(err > 0.0)
                        {
                            temp_ptHist->SetBinContent(ibin, temp_ptHist->GetBinContent(ibin) + gRandom->Gaus(0,err));
                        }
                    }
                    sysPtUnfold[sysName].at(nth)->SetInput(temp_ptHist, bias);
                }
                else if(var == "Mass")
                {
                    TH1* temp_massHist;

                    TString nth_;
                    nth_.Form("%d", nth);
                    temp_massHist = nomMassUnfold->GetInput("ptToy_" + nth_, 0, 0, 0, false);

                    // randomize histogram bins
                    for(int ibin = 1; ibin<temp_massHist->GetNbinsX()+1;ibin++)
                    {
                        double err = temp_massHist->GetBinError(ibin);
                        if(err > 0.0)
                        {
                            temp_massHist->SetBinContent(ibin, temp_massHist->GetBinContent(ibin) + gRandom->Gaus(0,err));
                        }
                    }
                    sysMassUnfold[sysName].at(nth)->SetInput(temp_massHist, bias);
                }
                else
                {
                    cout << "ISRUnfold::setUnfInput, only Pt and Mass available for var" << endl;
                    exit (EXIT_FAILURE);
                }

            }
        }
        else
        {
            //FIXME currently only consider lepton momentum systematic here
            TString histDirPostfix = "";
            if(nth == 0)
            {
                phase_name += "_lepMomUp";
                histDirPostfix = "_lepMomUp";
            }
            else if(nth == 1)
            {
                phase_name += "_lepMomDown";
                histDirPostfix = "_lepMomDown";
            }
            else
            {
                exit(EXIT_FAILURE);
            }

            if(var == "Pt")
            {
                if(channel_name == "muon")     hRec = (TH1*)filein->Get(phase_name + "/hist_ptll" + histDirPostfix + "/histo_DoubleMuonnominal");
                if(channel_name == "electron" && year != 2018) hRec = (TH1*)filein->Get(phase_name + "/hist_ptll" + histDirPostfix + "/histo_DoubleEGnominal");
                if(channel_name == "electron" && year == 2018) hRec = (TH1*)filein->Get(phase_name + "/hist_ptll" + histDirPostfix + "/histo_EGammanominal");
                sysPtUnfold[sysName].at(nth)  ->SetInput(hRec,   bias);
            }
            else if(var == "Mass")
            {
                if(channel_name == "muon")     hRec = (TH1*)filein->Get(phase_name + "/hist_mll" + histDirPostfix + "/histo_DoubleMuonnominal");
                if(channel_name == "electron" && year != 2018) hRec = (TH1*)filein->Get(phase_name + "/hist_mll" + histDirPostfix + "/histo_DoubleEGnominal");
                if(channel_name == "electron" && year == 2018) hRec = (TH1*)filein->Get(phase_name + "/hist_mll" + histDirPostfix + "/histo_EGammanominal");
                sysMassUnfold[sysName].at(nth)->SetInput(hRec,   bias);
            }
            else
            {
                cout << "ISRUnfold::setUnfInput, only Pt and Mass available for var" << endl;
                exit (EXIT_FAILURE);
            }
        }
    }
    filein->Close();
    delete filein;
}

void ISRUnfold::subBkgs(TString var, TString filepath, TString bkgName, bool isSys, TString sysName, int totSysN, int nth, TString phase_name)
{

	TFile* filein = new TFile(filepath);
        TH1* hRec = NULL;

        double bkg_scale = 1.;

        TString systematic_postfix = sysName;

        if(totSysN == 2){
            if(nth == 0 ) systematic_postfix="_"+systematic_postfix+"Up";
            if(nth == 1 ) systematic_postfix="_"+systematic_postfix+"Down";
        }

        if(totSysN == 6 && (sysName == "Scale" || sysName =="pdfScale")){
            if(nth == 0 ) systematic_postfix="_"+systematic_postfix+"AUp";
            if(nth == 1 ) systematic_postfix="_"+systematic_postfix+"ADown";
            if(nth == 2 ) systematic_postfix="_"+systematic_postfix+"BUp";
            if(nth == 3 ) systematic_postfix="_"+systematic_postfix+"BDown";
            if(nth == 4 ) systematic_postfix="_"+systematic_postfix+"ABUp";
            if(nth == 5 ) systematic_postfix="_"+systematic_postfix+"ABDown";
        }

        if(totSysN == 100 && sysName == "PDFerror"){
            TString nth_;
            nth_.Form ("%03d", nth);
            systematic_postfix="_"+systematic_postfix+nth_;
        }

        // nominal histograms
	if(!isSys)
        {
                if(var == "Pt")
                {
                    hRec = (TH1*)filein->Get(phase_name + "/"+var+"/histo_" + bkgName + "nominal");
                }
                if(var == "Mass")
                {
                    hRec = (TH1*)filein->Get(phase_name + "/"+var+"/histo_" + bkgName + "nominal");
                }

        	if( var == "Pt" )   nomPtUnfold->  SubtractBackground(hRec, bkgName, bkg_scale);
        	if( var == "Mass" ) nomMassUnfold->SubtractBackground(hRec, bkgName, bkg_scale);
	}
        // systematic histograms
	else
        {	

            // systematics using nominal detector distributions
            if(sysName=="Alt" || sysName=="unfoldBias" || sysName=="unfoldScan" || sysName == "Stat")
            {
                if(var == "Pt")
                {
                    hRec = (TH1*)filein->Get(phase_name + "/hist_ptll/histo_" + bkgName + "nominal");
                }
                if(var == "Mass")
                {
                    hRec = (TH1*)filein->Get(phase_name + "/hist_mll/histo_" + bkgName + "nominal");
                }
            }
            else
            {
                // WW WZ ZZ have no PDF systematics
                if( (bkgName=="WW_pythia" || bkgName=="WZ_pythia" || bkgName=="ZZ_pythia") && (sysName=="pdfScale"||sysName=="Scale" || sysName=="AlphaS"||sysName=="alphaS"))
                {
                    if(var == "Pt")
                    {
                        hRec = (TH1*)filein->Get(phase_name + "/hist_ptll/histo_" + bkgName + "nominal");
                    }
                    if(var == "Mass")
                    {
                        hRec = (TH1*)filein->Get(phase_name + "/hist_mll/histo_" + bkgName + "nominal");
                    }
                }
                else
                {

                    TString histDirPostfix = "";
                    if(sysName == "lepMom"){
                        if(nth == 0)
                        {
                            phase_name += "_lepMomUp";
                            histDirPostfix = "_lepMomUp";
                            systematic_postfix = "nominal";
                        }
                        else if(nth == 1)
                        {
                            phase_name += "_lepMomDown";
                            histDirPostfix = "_lepMomDown";
                            systematic_postfix = "nominal";
                        }
                        else
                        {
                            exit(EXIT_FAILURE);
                        }
                    }

                    if(var == "Pt")
                    {
                        hRec = (TH1*)filein->Get(phase_name + "/hist_ptll" + histDirPostfix + "/histo_" + bkgName + systematic_postfix);
                        cout << "path: " << phase_name + "/hist_ptll" + histDirPostfix + "/histo_" + bkgName + systematic_postfix << endl;
                    }
                    if(var == "Mass")
                    {
                        hRec = (TH1*)filein->Get(phase_name + "/hist_mll" + histDirPostfix + "/histo_" + bkgName + systematic_postfix);
                    }
                }
            }

            if( var == "Pt" )   sysPtUnfold[sysName].at(nth)  ->SubtractBackground(hRec, bkgName, bkg_scale);
            if( var == "Mass" ) sysMassUnfold[sysName].at(nth)->SubtractBackground(hRec, bkgName, bkg_scale);
	}
	
	filein->Close();
	delete filein;
}

void ISRUnfold::doISRUnfold(int detOrFSR_unfold, bool doSys){

    double tauMin=1.e-4;
    double tauMax=1.e-3;

    nScan=50;
    rhoLogTau=0;
    lCurve=0;
    iBest = 0;

    nScan_mass=50;
    rhoLogTau_mass=0;
    lCurve_mass=0;
    iBest_mass = 0;


    if(detOrFSR_unfold == 0)
    {
        if(regMode_detector == TUnfold::kRegModeNone)
        {
            // No regularisation, set tau as zero
            //const TUnfoldBinningV17* bin=nomPtUnfold->GetOutputBinning();
            //int istart=bin->GetGlobalBinNumber(0.1, 200.1);
            //int iend=bin->GetGlobalBinNumber(100-0.01,200.1);
            //nomPtUnfold->RegularizeBins(istart,1,iend-istart+1,TUnfoldV17::kRegModeCurvature);
            //iBest=nomPtUnfold->ScanTau(nScan,0.,0.,&rhoLogTau,
            //                           TUnfoldDensity::kEScanTauRhoAvgSys,
            //                           0,0,
            //                           &lCurve);
            nomPtUnfold->DoUnfold(0);
            nomMassUnfold->DoUnfold(0);
        }
        else
        {
            // regularization, use ScanTau() as a default method to find tau
            iBest=nomPtUnfold->ScanTau(nScan,0.,0.,&rhoLogTau,
                                       TUnfoldDensity::kEScanTauRhoAvgSys,
                                       0,0,
                                       &lCurve);

            iBest_mass=nomMassUnfold->ScanTau(nScan_mass,0.,0.,&rhoLogTau_mass,
                                       TUnfoldDensity::kEScanTauRhoAvgSys,
                                       0,0,
                                       &lCurve_mass);
        }
    }
    else if(detOrFSR_unfold == 1)
    // FSR unfolding
    {
        if(regMode_FSR == TUnfold::kRegModeNone)
        {
            // no regularisation, set tau as zero
            nomPtFSRUnfold->DoUnfold(0);
            nomMassFSRUnfold->DoUnfold(0);
        }
        else
        {
            // regularization option not ready yet
        }
    }
    else{

        cout << " invalid unfolding step..." << endl;
        exit(EXIT_FAILURE);
    }

    // for systematic
    if(doSys)
    {
        std::map<TString, std::vector<TUnfoldDensityV17*>>::iterator it;
        std::map<TString, std::vector<TUnfoldDensityV17*>>::iterator it_end;

        // unfolding for pt distribution
        if(detOrFSR_unfold == 0)
        {
            it = sysPtUnfold.begin();
            it_end = sysPtUnfold.end();
        }
        else if(detOrFSR_unfold == 1)
        {
            it = sysPtFSRUnfold.begin();
            it_end = sysPtFSRUnfold.end();
        }
        else
        {
            cout << " invalid unfolding step..." << endl;
            exit(EXIT_FAILURE);
        }

        while(it != it_end){
        	int nSys = it->second.size();
        	for(int i = 0; i < nSys; i++){

        		if((it->first)=="unfoldScan" || (it->first)=="unfoldBias" )
                    {
        			//it->second.at(i)->ScanLcurve(50,tauMin,tauMax,0);
        			iBest=it->second.at(i)->ScanTau(nScan,0.,0.,&rhoLogTau,
        			                           TUnfoldDensity::kEScanTauRhoAvgSys,
        			                           0,0,
        			                           &lCurve);
        		}
        		else
                    {
                        it->second.at(i)->DoUnfold(0);
                    }

        	}
        	it++;
        }

        // unfolding for mass distribution
        if(detOrFSR_unfold == 0)
        {
            it = sysMassUnfold.begin();
            it_end = sysMassUnfold.end();
        }
        else if(detOrFSR_unfold == 1)
        {
            it = sysMassFSRUnfold.begin();
            it_end = sysMassFSRUnfold.end();
        }
        else
        {
            cout << " invalid unfolding step..." << endl;
            exit(EXIT_FAILURE);
        }

        while(it != it_end){
                int nSys = it->second.size();
                for(int i = 0; i < nSys; i++){
                        if((it->first)=="unfoldScan" || (it->first)=="unfoldBias"){
        			//it->second.at(i)->ScanLcurve(50,tauMin,tauMax,0);

                                iBest_mass=it->second.at(i)->ScanTau(nScan_mass,0.,0.,&rhoLogTau_mass,
                                                           TUnfoldDensity::kEScanTauRhoAvgSys,
                                                           0,0,
                                                           &lCurve_mass);
        		}
        		else it->second.at(i)->DoUnfold(0);

                }
                it++;
        }
    }

}

double ISRUnfold::Chi2Test(TH1 *data, TH1 *mc){ // check if this is right way to get chi square value

 gROOT->SetBatch();

 double chi2 = 0.;
 double ndf  = 0.;

 for(int i=1;i<=mc->GetNbinsX();i++) {
    ndf += 1.;
    if(data->GetBinError(i) == 0){
      std::cout << "unfolded " << i << " bin: " << data->GetBinContent(i) << std::endl;
      std::cout << "error is zero in the " << i << " bin..." << std::endl;
      std::cout << "so skip this bin" << std::endl;
      continue;
    }
    double pull=(data->GetBinContent(i)-mc->GetBinContent(i))/data->GetBinError(i);
    chi2+= pull*pull;
 }
 return chi2/(ndf-1);

}

// need unfolded hist, rho matrix (GetRhoIJtotal), MC truth
double ISRUnfold::DoFit(TString var, int nthMassBin, bool isFSRUnfold)
{

    TH1 *g_fcnHist=0;
    TMatrixD *g_fcnMatrix=0;

    TH1* hpt_temp_data;
    TH1* hpt_temp_mc;
    TH2 *hrho;

    TString ibinMass;
    ibinMass.Form("%d", nthMassBin);

    if(var=="Pt"){
        if(!isFSRUnfold){
            hpt_temp_data = nomPtUnfold->GetOutput("hunfolded_pt_temp_chi",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
            hpt_temp_mc   = nomPtUnfold->GetBias("histMC_pt_temp_chi",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
            hrho          = nomPtUnfold->GetRhoIJtotal("histRho_chi", 0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
        }
        else{
            hpt_temp_data = nomPtFSRUnfold->GetOutput("hunfolded_pt_temp_chi",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
            hpt_temp_mc   = nomPtFSRUnfold->GetBias("histMC_pt_temp_chi",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
            hrho          = nomPtFSRUnfold->GetRhoIJtotal("histRho_chi", 0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
        }
    }
    else{
        if(!isFSRUnfold){
            hpt_temp_data = nomMassUnfold->GetOutput("hunfolded_pt_temp_chi",0,0,"mass[UO];pt[UOC0]",kTRUE);
            hpt_temp_mc   = nomMassUnfold->GetBias("histMC_pt_temp_chi",0,0,"mass[UO];pt[UOC0]",kTRUE);
            hrho          = nomMassUnfold->GetRhoIJtotal("histRho_chi", 0,0,"mass[UO];pt[UOC0]",kTRUE);
        }
        else{
            hpt_temp_data = nomMassFSRUnfold->GetOutput("hunfolded_pt_temp_chi",0,0,"mass[UO];pt[UOC0]",kTRUE);
            hpt_temp_mc   = nomMassFSRUnfold->GetBias("histMC_pt_temp_chi",0,0,"mass[UO];pt[UOC0]",kTRUE);
            hrho          = nomMassFSRUnfold->GetRhoIJtotal("histRho_chi", 0,0,"mass[UO];pt[UOC0]",kTRUE);
        }
    }

    //cout << "rho matrix dimention: nbinx " << hrho->GetXaxis()->GetNbins() << " nbiny: " << hrho->GetXaxis()->GetNbins() << endl;

    int n=hpt_temp_data->GetNbinsX();

    TMatrixDSym v(n);
    //g_fcnMatrix=new TMatrixD(n,n);
    for(int i=0;i<n;i++) {
       for(int j=0;j<n;j++) {
          v(i,j)=hrho->GetBinContent(i+1,j+1)*(hpt_temp_data->GetBinError(i+1)*hpt_temp_data->GetBinError(j+1));
       }
    }

    TMatrixDSymEigen ev(v);
    TMatrixD d(n,n);
    TVectorD di(ev.GetEigenValues());
    for(int i=0;i<n;i++) {
       if(di(i)>0.0) {
          d(i,i)=1./di(i);
       } else {
          cout<<"bad eigenvalue i="<<i<<" di="<<di(i)<<"\n";
          exit(0);
       }
    }

    TMatrixD O(ev.GetEigenVectors());
    TMatrixD DOT(d,TMatrixD::kMultTranspose,O);
    g_fcnMatrix=new TMatrixD(O,TMatrixD::kMult,DOT);
    TMatrixD test(*g_fcnMatrix,TMatrixD::kMult,v);
    int error=0;

    for(int i=0;i<n;i++) {
       if(TMath::Abs(test(i,i)-1.0)>1.E-7) {
          error++;
       }
       for(int j=0;j<n;j++) {
          if(i==j) continue;
          if(TMath::Abs(test(i,j)>1.E-7)) error++;
       }
    }

    // calculate chi2
    double chi2 = 0.;
    double ndf = 0.;
    for(int i=0;i<hpt_temp_data->GetNbinsX();i++) {
       double di_=hpt_temp_data->GetBinContent(i+1)-hpt_temp_mc->GetBinContent(i+1);
       if(g_fcnMatrix) {
          for(int j=0;j<hpt_temp_data->GetNbinsX();j++) {
             double dj=hpt_temp_data->GetBinContent(j+1)-hpt_temp_mc->GetBinContent(j+1);
             chi2+=di_*dj*(*g_fcnMatrix)(i,j);
          }
       } else {
          double pull=di_/hpt_temp_data->GetBinError(i+1);
          chi2+=pull*pull;
       }
       ndf+=1.0;
    }

    delete g_fcnHist;
    delete g_fcnMatrix;

    delete hpt_temp_data;
    delete hpt_temp_mc;
    delete hrho;

    return chi2/ndf;
}

int ISRUnfold::setMeanMass(bool isDet)
{
    cout << "ISRUnfold::setMeanMass()   Save mean of dilepton..." << endl;
    
    // Get mass bin definition for post FSR level from TUnfoldDensity object
    // Nothe currently post and pre FSR level use the same bin definition
    const TUnfoldBinningV17* temp_binning_gen_pt = pt_binning_Gen;
    const TVectorD* temp_tvecd = temp_binning_gen_pt->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;
    const Double_t* massBins = temp_tvecd->GetMatrixArray();

    TUnfoldDensityV17* p_unfold = NULL;
    if(isDet)
        p_unfold = nomMassUnfold;
    else
        p_unfold = nomMassFSRUnfold;

    TH1 * hdetector_mass = p_unfold->GetInput("hdetector_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);
    TH1* hunfolded_mass =  p_unfold->GetOutput("hunfolded_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);
    TH1* hMC_mass =  p_unfold->GetBias("hunfolded_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);

    // Loop over mass bins
    for(int ibin = 0; ibin < nMassBin; ibin++)
    {
        // set x-axis range
        hdetector_mass->GetXaxis()->SetRange(hdetector_mass->GetXaxis()->FindBin(massBins[ibin]+0.01), hdetector_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
        hunfolded_mass->GetXaxis()->SetRange(hunfolded_mass->GetXaxis()->FindBin(massBins[ibin]+0.01), hunfolded_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
        hMC_mass->GetXaxis()->SetRange(hunfolded_mass->GetXaxis()->FindBin(massBins[ibin]+0.01), hunfolded_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));

        if(isDet)
        {
            // Get mean values
            cout << "Detector, " << ibin << " th mass bin, mean: " << hdetector_mass->GetMean() << " +/- " << hdetector_mass->GetMeanError() << endl;
            meanMass_data_detector. push_back(hdetector_mass->GetMean());
            meanMassStatErr_data_detector.push_back(hdetector_mass->GetMeanError());

            cout << "Unfolded, " << ibin << " th mass bin, mean: " << hunfolded_mass->GetMean() << " +/- " << hunfolded_mass->GetMeanError() << endl;
            meanMass_data_det_unf.   push_back(hunfolded_mass->GetMean());
            meanMassStatErr_data_det_unf.push_back(hunfolded_mass->GetMeanError());

            cout << "MC, " << ibin << " th mass bin, mean: " << hMC_mass->GetMean() << " +/- " << hMC_mass->GetMeanError() << endl;
            meanMass_mc_det_unf.   push_back(hMC_mass->GetMean());
            meanMassStatErr_mc_det_unf.push_back(hMC_mass->GetMeanError());
        }
        else
        {
            meanMass_data_pre_fsr.   push_back(hunfolded_mass->GetMean());
            meanMassStatErr_data_pre_fsr.push_back(hunfolded_mass->GetMeanError());
        }

    }// end of mass bin loop

    delete hdetector_mass;
    delete hunfolded_mass;
    delete hMC_mass;

    return nMassBin;
}

double ISRUnfold::getDetMeanMass(int ibin)
{

    int size = meanMass_data_detector.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);     
    }
    else
    {
        return meanMass_data_detector.at(ibin);
    }
}

double ISRUnfold::getUnfMeanMass(int ibin)
{

    int size = meanMass_data_det_unf.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);     
    }
    else
    {
        return meanMass_data_det_unf.at(ibin);
    }
}

double ISRUnfold::getFSRUnfMeanMass(int ibin)
{

    int size = meanMass_data_pre_fsr.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);     
    }
    else
    {
        return meanMass_data_pre_fsr.at(ibin);
    }
}

double ISRUnfold::getUnfMeanMassError(int ibin)
{

    int size = meanMassStatErr_data_det_unf.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);     
    }
    else
    {
        return meanMassStatErr_data_det_unf.at(ibin);
    }
}

double ISRUnfold::getFSRUnfMeanMassError(int ibin)
{

    int size = meanMassStatErr_data_pre_fsr.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);     
    }
    else
    {
        return meanMassStatErr_data_pre_fsr.at(ibin);
    }
}

double ISRUnfold::getMCGenMeanMass(int ibin)
{

    int size = meanMass_mc_det_unf.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);     
    }
    else
    {
        return meanMass_mc_det_unf.at(ibin);
    }
}

// set mean pt from mass and DY mc
int ISRUnfold::setMeanPt(bool isDet)
{
    cout << "ISRUnfold::setMeanPt()   Save mean of dilepton momentum..." << endl;

    // Find number of mass bins
    const TUnfoldBinningV17* temp_binning_gen_pt = pt_binning_Gen;
    const TVectorD* temp_tvecd = temp_binning_gen_pt->GetDistributionBinning(1);
    int nMassBin = nMassBin = temp_tvecd->GetNrows() - 1;

    TUnfoldDensityV17* p_unfold = NULL;
    if(isDet)
        p_unfold = nomPtUnfold;
    else
        p_unfold = nomPtFSRUnfold;

    // Save mean pt 
    for(int i = 0; i < nMassBin; i++)
    {
        TString ibinMass;
        ibinMass.Form("%d", i);

        TH1* hpt_temp_data;
        TH1* hpt_temp_mc;

        // get histograms to set mean values
        TH1* hdetector_data = p_unfold->GetInput("h_detector_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
        hpt_temp_data = p_unfold->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
        hpt_temp_mc   = p_unfold->GetBias("histMC_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE); 

        if(isDet)
        {
            cout << "Detector, " << i << " th mass bin, mean: " << hdetector_data->GetMean() << " +/- " << hdetector_data->GetMeanError() << endl;
            meanPt_data_detector.push_back(hdetector_data->GetMean());
            meanPtStatErr_data_detector.push_back(hdetector_data->GetMeanError());

            cout << "Unfolded, " << i << " th mass bin, mean: " << hpt_temp_data->GetMean() << " +/- " << hpt_temp_data->GetMeanError() << endl;
            meanPt_data_det_unf.push_back(hpt_temp_data->GetMean());
            meanPtStatErr_data_det_unf.push_back(hpt_temp_data->GetMeanError());

            meanPt_mc_det_unf.push_back(hpt_temp_mc->GetMean());
            meanPtErr_mc_det_unf.push_back(hpt_temp_mc->GetMeanError());
        }
        else
        {
            meanPt_data_pre_fsr.push_back(hpt_temp_data->GetMean());
            meanPtStatErr_data_pre_fsr.push_back(hpt_temp_data->GetMeanError());
        }

        delete hdetector_data;
        delete hpt_temp_data;
        delete hpt_temp_mc;
    }

    return nMassBin;
}

double ISRUnfold::getDetMeanPt(int ibin)
{

    int size = meanPt_data_detector.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);     
    }
    else
    {
        return meanPt_data_detector.at(ibin);
    }
}

double ISRUnfold::getUnfMeanPt(int ibin)
{

    int size = meanPt_data_det_unf.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);     
    }
    else
    {
        return meanPt_data_det_unf.at(ibin);
    }
}

double ISRUnfold::getFSRUnfMeanPt(int ibin)
{

    int size = meanPt_data_pre_fsr.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);     
    }
    else
    {
        return meanPt_data_pre_fsr.at(ibin);
    }
}

double ISRUnfold::getUnfMeanPtError(int ibin)
{

    int size = meanPtStatErr_data_det_unf.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);     
    }
    else
    {
        return meanPtStatErr_data_det_unf.at(ibin);
    }
}

double ISRUnfold::getFSRUnfMeanPtError(int ibin)
{

    int size = meanPtStatErr_data_pre_fsr.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);     
    }
    else
    {
        return meanPtStatErr_data_pre_fsr.at(ibin);
    }
}

double ISRUnfold::getMCGenMeanPt(int ibin)
{

    int size = meanPt_mc_det_unf.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);     
    }
    else
    {
        return meanPt_mc_det_unf.at(ibin);
    }
}

// https://root.cern.ch/root/html/tutorials/graphs/graphtext.C.html
void ISRUnfold::drawtext(TGraph *g)
{
    Int_t i,n;
    Double_t x,y;
    TLatex *l;

    n = g->GetN();
    for (i=0; i<n; i++)
    {
       g->GetPoint(i,x,y);
       l = new TLatex(x-0.1,y+0.02,Form("%d",i+1));
       l->SetTextSize(0.015);
       l->SetTextFont(42);
       l->SetTextAlign(21);
       l->Draw();
   }
}

void ISRUnfold::doNorm(TH1* hist, bool norm)
{
    for(int ibin = 1; ibin < hist->GetXaxis()->GetNbins()+1; ibin++)
    {
        double binWidth = hist->GetBinWidth(ibin);
        hist->SetBinContent(ibin, hist->GetBinContent(ibin)/ binWidth);
        hist->SetBinError(ibin, hist->GetBinError(ibin)/ binWidth);
    }
    if(norm)
        hist->Scale(1./ hist->Integral());
}

TH1* ISRUnfold::getDetUnfoldedHists(TString var, TString outHistName, TString steering, bool useAxis)
{
    if(var == "Mass")
        return nomMassUnfold->GetOutput(outHistName,0,0,steering,useAxis);
    else
        return nomPtUnfold->GetOutput(outHistName,0,0,steering,useAxis);
}

TH1* ISRUnfold::getFSRUnfoldedHists(TString var, TString outHistName, TString steering, bool useAxis)
{
    if(var == "Mass")
        return nomMassFSRUnfold->GetOutput(outHistName,0,0,steering,useAxis);
    else
        return nomPtFSRUnfold->GetOutput(outHistName,0,0,steering,useAxis);
}

TH1* ISRUnfold::getMCHists(TString var, TString outHistName, TString steering, bool useAxis)
{
    if(var == "Mass")
        return nomMassUnfold->GetBias(outHistName,0,0,steering,useAxis);
    else 
        return nomPtUnfold->GetBias(outHistName,0,0,steering,useAxis);
}

TH1* ISRUnfold::getDetHists(TString var, TString outHistName, TString steering, bool useAxis)
{
    
    if(var == "Mass")
        return nomMassUnfold->GetInput(outHistName,0,0,steering,useAxis);
    else
        return nomPtUnfold->GetInput(outHistName,0,0,steering,useAxis);
}

TH1* ISRUnfold::getRawHist(TString filePath, TString histName, TString outHistName, TString steering)
{
    //const TUnfoldBinning* recM_binning = nomMassUnfold->GetInputBinning("Rec_Mass");

    TH1::AddDirectory(kFALSE);
    TFile* filein = new TFile(filePath);
    //const TUnfoldBinning* recM_binning = (TUnfoldBinning*)filein->Get("detector_level/hist_mll/Rec_Mass");

    TH1* raw_hist = (TH1*)filein->Get(histName);
    raw_hist->SetName(outHistName);
    //TH1* hist = recM_binning->ExtractHistogram(outHistName, raw_hist, 0, true, steering);
    delete filein;
    
    return raw_hist; 
}
