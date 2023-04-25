#include "ISR_unfoldUtils.h"

void ISRUnfold::setBias(double bias)
{
   nominal_bias = bias;
}

TMatrixD ISRUnfold::makeMatrixFromHist(TH2F*hist)
{
    int nBinsX = hist->GetNbinsX();
    int nBinsY = hist->GetNbinsY();
    TMatrixD matrix(nBinsY,nBinsX);
    for(int i=1;i<=nBinsX;i++)
    {
        for(int j=1;j<=nBinsY;j++)
        {
            matrix(j-1,i-1) = hist->GetBinContent(i,j);
        }
    }
    return matrix;
}//end makeMatrixFromHist

// Set the nominal response matrix
void ISRUnfold::setNominalRM(TString file_path, TString channel_name)
{
    TH1::AddDirectory(kFALSE);
    TFile* filein = new TFile(file_path, "READ");

    // bin definition
    TString Rec_binName = var + "_smeared_bin";
    TString Gen_binName = var + "_truth_bin";

    // Set bin definition
    binningFine = (TUnfoldBinning*)filein->Get(channel_name + "/" + Rec_binName);
    binningCoarse = (TUnfoldBinning*)filein->Get(channel_name + "/" + Gen_binName);

    // Set response matrix
    // First, get the response matrix
    TH2* hmcGenRec = (TH2*)filein->Get(channel_name + "/" + var + "_responseM");

    nominalTUnfold = new TUnfoldDensity(hmcGenRec,
                                        TUnfold::kHistMapOutputHoriz,
                                        regMode,
                                        TUnfold::kEConstraintArea,
                                        densityMode,
                                        binningCoarse,
                                        binningFine);
    hResponseM = (TH2*) hmcGenRec->Clone("hResponseM");
    save_hists_from_responseM();

    cout << "TUnfold version: " << nominalTUnfold->GetTUnfoldVersion() << endl;

    filein->Close();
    delete filein;
}

// Option for unfold options
void ISRUnfold::setSystematicRM(TString file_path, TString channel_name, TString sys_name, TString hist_postfix)
{
    TH1::AddDirectory(kFALSE);
    TFile* filein = new TFile(file_path, "READ");

    TH2* hmcGenRec = (TH2*)filein->Get(channel_name + "/" + var + "_responseM" + hist_postfix);

    if(sys_name.Contains("IterEM"))
    {
        iterEMTUnfold = new TUnfoldIterativeEM(hmcGenRec,TUnfoldDensity::kHistMapOutputHoriz,binningCoarse,binningFine);
    }
    else
    {
        //cout << sys_name << " crate TUnfoldDensity..." << endl;
        if(hmcGenRec == NULL) cout << "check input file" << endl;
        systematicTUnfold[sys_name] = new TUnfoldDensity(hmcGenRec,
                                                        TUnfold::kHistMapOutputHoriz,
                                                        regMode,
                                                        TUnfold::kEConstraintArea,
                                                        densityMode,
                                                        binningCoarse,
                                                        binningFine);
    }
    
    TDirectory* varDir;
    varDir=fUnfoldOut->GetDirectory("unfolded/"+var);
    varDir->cd();
    
    TH1D* hProjectedTruth = (TH1D*) hmcGenRec->ProjectionX("histo_DY_"+sys_name, 0, -1, "e");
    hProjectedTruth->Write();

    filein->Close();
    
    delete filein;
    delete hProjectedTruth;
}

void ISRUnfold::save_hists_from_responseM()
{
    TDirectory* varDir;
    TDirectory* varDirForReco;

    varDirForReco=fUnfoldOut->GetDirectory("folded/"+var);
    varDir=fUnfoldOut->GetDirectory("matrix/"+var);
    varDir->cd();
    binningCoarse->Write();
    binningFine->Write();

    TH2F* hMigrationM = (TH2F*) nominalTUnfold->GetProbabilityMatrix("hMigrationM");
    hMigrationM->Write();
    hResponseM->SetName("hResponseM");
    hResponseM->Write();

    // Save projection of the migration matrix
    varDir=fUnfoldOut->GetDirectory("unfolded/"+var);
    varDir->cd();

    TH1D* hProjectedTruth = (TH1D*) hResponseM->ProjectionX("histo_DY", 0, -1, "e"); 
    TH1D* hProjectedBinZero = (TH1D*) hResponseM->ProjectionX("histo_ProjectedBinZero", 0, 0, "e");  //
    TH1D* hProjectedReco = (TH1D*) hResponseM->ProjectionY("histo_ProjectedReco", 1, -1, "e");  //
    hProjectedTruth->Write();
    hProjectedBinZero->Write();

    varDirForReco->cd();
    binningFine->Write();
    hProjectedReco->Write();

    varDir->cd();
}

// Set input histogram from root file
void ISRUnfold::setUnfInput(TString file_path, TString channel_name, TString sys_type, TString sys_name, TString sys_hist_postfix, bool isFSR)
{
    TH1::AddDirectory(kFALSE);

    TFile* filein = new TFile(file_path);

    TString full_hist_path = channel_name + "/" + var + "_smeared" + sys_hist_postfix;
    TH1* hRec = (TH1*)filein->Get(full_hist_path);

    // Nominal
    if(sys_type == "Type_0")
    {
        nominalTUnfold->SetInput(hRec, nominal_bias);
    }
    else
    {
        // Systematic histograms
        if(sys_name.Contains("IterEM"))
        {
            iterEMTUnfold->SetInput(hRec, 1.);
        }
        else
        {
            systematicTUnfold[sys_name]->SetInput(hRec, nominal_bias);
        }
    }

    filein->Close();
    delete filein;
}

void ISRUnfold::subBkgs(TString file_path, TString channel_name, TString bkg_name, TString sys_type, TString sys_name, TString hist_postfix)
{
    TFile* filein = new TFile(file_path);
    TH1* hRec = (TH1*)filein->Get(channel_name + "/" + var + "_smeared" + hist_postfix);

    // Nominal histograms
    if(sys_type=="Type_0")
    {
        nominalTUnfold->SubtractBackground(hRec, bkg_name);

        // Save Unfolding fake histogram in folded directory
        if(channel_name.Contains("Fake"))
        {
            TDirectory* varDirForReco;

            varDirForReco=fUnfoldOut->GetDirectory("folded/"+var);
            varDirForReco->cd();

            hRec->SetName("histo_Fake"+bkg_name);
            hRec->Write();
        }
    }
    else
    // Systematic
    {
        if(sys_name.Contains("IterEM"))
        {
            iterEMTUnfold->SubtractBackground(hRec, bkg_name);
        }
        else
        {
            // FIXME temporary method for background systematic
            if(sys_name.Contains("Background"))
            {
                if(sys_name.Contains("Up"))
                {
                    systematicTUnfold[sys_name]->SubtractBackground(hRec, bkg_name, 1.05);
                }
                if(sys_name.Contains("Down"))
                {
                    systematicTUnfold[sys_name]->SubtractBackground(hRec, bkg_name, 0.95);
                }
            }
            else
            {
                cout << "bkg name " << bkg_name << " sys type " << sys_type << " sys name " << sys_name << endl;
                if (hRec == NULL) cout << "check histogram!" << endl;
                systematicTUnfold[sys_name]->SubtractBackground(hRec, bkg_name);
            }
        }
    }

    filein->Close();
    delete filein;
}

void ISRUnfold::set_systematics(TString sysHistName)
{
    sysVector.push_back(sysHistName);
}

void ISRUnfold::setPartialRegularize2D(TUnfold::ERegMode partialRegMode, double startMass, double startPt, double endMass, double endPt)
{

    int istart=binningCoarse->GetGlobalBinNumber(startPt, startMass);
    int iend=binningCoarse->GetGlobalBinNumber(endPt, endMass);

    //int iend=binningCoarse->GetGlobalBinNumber(50-0.01,320);
    //unfold->RegularizeBins(istart, 1, iend-istart+1, TUnfoldV17::kRegModeSize);
    //unfold->RegularizeBins(istart, 1, iend-istart+1, TUnfoldV17::kRegModeDerivative);

    nominalTUnfold->RegularizeBins(istart, 1, iend-istart+1, partialRegMode);

    std::vector<TString>::iterator it = sysVector.begin();
    while(it != sysVector.end())
    {
        systematicTUnfold[*it]->RegularizeBins(istart, 1, iend-istart+1, partialRegMode);
        it++;
    }
}
//void ISRUnfold::setPartialRegularize1D()

// Unfolding
void ISRUnfold::doISRUnfold(bool partialReg)
{

    TDirectory* topDir;
    TDirectory* varDir;

    topDir=fUnfoldOut->GetDirectory("unfolded");
    varDir=fUnfoldOut->GetDirectory("unfolded/"+var);
    varDir->cd();
    binningCoarse->Write();

    topDir->cd();

    // No regularisation
    if(regMode == TUnfold::kRegModeNone)
    {
        if (partialReg) {

            //setPartialRegularize2D(TUnfold::kRegModeCurvature, 200., 0., 1000., 100.);
            setPartialRegularize2D(TUnfold::kRegModeCurvature, 320., 0., 1000., 100.);
            Int_t nScan=1000;
            Double_t tauMin = 1e-5; //If tauMin=tauMax, TUnfold automatically chooses a range
            Double_t tauMax = 1e-1; //Not certain how they choose the range

            // L-curve scan
            // much senstive to tauMin and tauMax when using ScanLcurve
            iBest_nominal=nominalTUnfold->ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);

            // Tau scan
            //TSpline *rhoScan=0;
            //nominalTUnfold->ScanTau(nScan, tauMin, tauMax, &rhoScan, TUnfoldDensity::kEScanTauRhoAvg);
            //cout << "tau: " << nominalTUnfold->GetTau() << endl;

            tau = nominalTUnfold->GetTau();
        }
        else {
            nominalTUnfold->DoUnfold(tau);
        }
    }
    else // apply regularization for all bins
    {
        Int_t nScan=100;
        Double_t tauMin = 1e-5; //If tauMin=tauMax, TUnfold automatically chooses a range
        Double_t tauMax = 1.; //Not certain how they choose the range

        // ScanLcurve
        //iBest_nominal=nominalTUnfold->ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);
        //tau=nominalTUnfold->GetTau(); // Use this tau for systematic unfolding
        //cout<< "tau=" << nominalTUnfold->GetTau() << endl;
        //cout<<"chi**2="<<nominalTUnfold->GetChi2A()<<"+"<<nominalTUnfold->GetChi2L()<<" / "<<nominalTUnfold->GetNdf()<<"\n";

        TSpline *rhoScan=0;
        nominalTUnfold->ScanTau(nScan, tauMin, tauMax, &rhoScan, TUnfoldDensity::kEScanTauRhoAvg);
        cout << "tau: " << nominalTUnfold->GetTau() << endl;
        tau = nominalTUnfold->GetTau();
    }

    varDir->cd();

    bool useAxisBinning = false;
    if(var == "1D_dimass")
    {
        useAxisBinning = true;
    }

    nominalTUnfold->GetRhoIJtotal("hCorrelation", 0, 0, 0, useAxisBinning)->Write();
    nominalTUnfold->GetEmatrixTotal("hCovariance", 0, 0, 0, useAxisBinning)->Write();
    nominalTUnfold->GetOutput("histo_Data",0,0, "*[*]", useAxisBinning)->Write();

    // TODO make a function for saving
    // Save unfolding input histogram
    TDirectory* varDirForReco;

    varDirForReco=fUnfoldOut->GetDirectory("folded/"+var);
    varDirForReco->cd();

    nominalTUnfold->GetInput("histo_UnfoldInput", 0, 0, 0, false)->Write();

    varDir->cd();

    // For systematic
    std::vector<TString>::iterator it = sysVector.begin();
    while(it != sysVector.end())
    {
        if( (*it).Contains("IterEM"))
        {
            iBest=iterEMTUnfold->ScanSURE(NITER_Iterative, &graph_SURE_IterativeSURE, &graph_DFdeviance_IterativeSURE);
            //cout << "iBest: " << iBest << endl;

            varDir->cd();
            iterEMTUnfold->GetOutput("histo_Data_"+(*it),0,0, "*[*]", false)->Write();
        }
        else{
            systematicTUnfold[*it]->DoUnfold(tau);

            varDir->cd();
            systematicTUnfold[*it]->GetOutput("histo_Data_"+(*it),0,0, "*[*]", false)->Write();
        }
        it++;
    }

    topDir->Write();
}


// In this function, acceptance factors change only for PDF, Scale, AlphaS systematics
void ISRUnfold::doAcceptCorr(TString filePath, TString filePath_for_accept)
{
    TDirectory* topDir;
    TDirectory* varDir;

    TH1* hAcceptance_raw = NULL;
    TFile* filein = new TFile(filePath);

    //if(!gSystem->AccessPathName(fullPath, kFileExists))
    topDir=fUnfoldOut->GetDirectory("acceptance");
    varDir=fUnfoldOut->GetDirectory("acceptance/"+var);
    varDir->cd();
    binningCoarse->Write();
    topDir->cd();

    TH1* hFiducialPhaseMC = NULL;

    hFullPhaseMC = (TH1*) filein->Get("Acceptance/"+var+ "_" + "/histo_DYJets");
    if(year.Contains("2016"))
        hFullPhaseMC->Add((TH1*) filein->Get("Acceptance/"+var+ "_" + "/histo_DYJets10to50"));
    else
        hFullPhaseMC->Add((TH1*) filein->Get("Acceptance/"+var+ "_" + "/histo_DYJets10to50_MG"));

    hFiducialPhaseMC = (TH1*)fUnfoldOut->Get("unfolded/"+var+"/"+"histo_DY");
    TString mass_binned_sample_prefix[9] = {"M-100to200", "M-200to400", "M-400to500", "M-500to700", "M-700to800",
                                            "M-800to1000", "M-1000to1500", "M-1500to2000", "M-2000to3000"};

    TH1* hFullPhaseMC_massBinned = NULL;
    if(filePath_for_accept != "")
    {
        int istart = binningCoarse->GetGlobalBinNumber(0, 200);

        // full phase histogram from mass binned sample
        for(int i = 0; i < 9; i++)
        {
            if(i==0)
            {
                hFullPhaseMC_massBinned = (TH1*) filein->Get("Acceptance/"+var+ "_" + "/histo_DYJets_" + mass_binned_sample_prefix[i]);
            }
            else
            {
                hFullPhaseMC_massBinned->Add((TH1*) filein->Get("Acceptance/"+var+ "_" + "/histo_DYJets_" + mass_binned_sample_prefix[i]));
            }
        }

        //
        for(int j = istart; j < hFullPhaseMC_massBinned->GetNbinsX()+1; j++)
        {
            hFullPhaseMC->SetBinContent(j, hFullPhaseMC_massBinned->GetBinContent(j));
            hFullPhaseMC->SetBinError(j,   hFullPhaseMC_massBinned->GetBinError(j));
        }

        // fiducial phase histogram from mass binned sample
        TH1* hFiducialPhaseMC_massBinned = NULL;
        for(int i = 0; i < 9; i++)
        {
            if(i==0)
            {
                hFiducialPhaseMC_massBinned = (TH1*) filein->Get("Acceptance_Efficiency/"+var+ "_" + "/histo_DYJets_" + mass_binned_sample_prefix[i]);
            }
            else
            {
                hFiducialPhaseMC_massBinned->Add((TH1*) filein->Get("Acceptance_Efficiency/"+var+ "_" + "/histo_DYJets_" + mass_binned_sample_prefix[i]));
            }
        }

        for(int j = istart; j < hFullPhaseMC_massBinned->GetNbinsX()+1; j++)
        {
            hFiducialPhaseMC->SetBinContent(j, hFiducialPhaseMC_massBinned->GetBinContent(j));
            hFiducialPhaseMC->SetBinError(j,   hFiducialPhaseMC_massBinned->GetBinError(j));
        }

        TFile* filein_temp = new TFile(filePath_for_accept);

        hAcceptance = (TH1*) filein_temp->Get("updated_accept_hist");
        filein_temp->Close();

        hAcceptance_raw = (TH1*) hFullPhaseMC->Clone("hAcceptance_raw");
        hAcceptance_raw->Divide(hFiducialPhaseMC);
    }
    else
    {
        hAcceptance = (TH1*) hFullPhaseMC->Clone("hAcceptance"+var);
        hAcceptance->Divide(hFiducialPhaseMC); // Nominal acceptance factor

        hAcceptance_raw = (TH1*) hAcceptance->Clone("hAcceptance_raw");
    }

    hFullPhaseData = nominalTUnfold->GetOutput("histo_Data",0,0, "*[*]", false);
    hFullPhaseData->Multiply(hAcceptance);

    varDir->cd();
    hFullPhaseData->Write();
    hFullPhaseMC->SetName("histo_DY");
    hFullPhaseMC->Write();
    hAcceptance->Write();

    std::vector<TString>::iterator it = sysVector.begin();
    while(it != sysVector.end())
    {
        //cout << (*it) << endl;
        TH1* hFullPhaseMC_raw_sys = NULL;
        TH1* hFiducialPhaseMC_sys = NULL;

        if((*it).Contains("IterEM"))
        {
            hSysFullPhaseData[*it]   = iterEMTUnfold->GetOutput("histo_Data_"+(*it),0,0, "*[*]", false);
            hFiducialPhaseMC_sys=hFiducialPhaseMC;
            hFiducialPhaseMC_sys->SetName("histo_DY_"+(*it));
        }
        else
        {
            hSysFullPhaseData[*it]   = systematicTUnfold[*it]->GetOutput("histo_Data_"+(*it),0,0, "*[*]", false);
            hFiducialPhaseMC_sys = (TH1*)fUnfoldOut->Get("unfolded/"+var+"/"+"histo_DY_"+(*it));
            //hFiducialPhaseMC_sys = hFiducialPhaseMC;
            hFiducialPhaseMC_sys->SetName("histo_DY_"+(*it));

            if(filePath_for_accept != "")
            {
                // fiducial phase histogram from mass binned sample

                int istart = binningCoarse->GetGlobalBinNumber(0, 200);
                TH1* hFiducialPhaseMC_sys_massBinned = NULL;
                for(int i = 0; i < 9; i++)
                {
                    if(i==0)
                    {
                        if( (((*it).Contains("Scale") && !(*it).Contains("Lep")) || (*it).Contains("PDF") || (*it).Contains("AlphaS")) )
                            hFiducialPhaseMC_sys_massBinned = (TH1*) filein->Get("Acceptance_Efficiency/"+var+ "_" +  "/histo_DYJets_" + mass_binned_sample_prefix[i] + "_" + (*it));
                        else
                            hFiducialPhaseMC_sys_massBinned = (TH1*) filein->Get("Acceptance_Efficiency/"+var+ "_" + "/histo_DYJets_" + mass_binned_sample_prefix[i]);
                    }
                    else
                    {
                        if( (((*it).Contains("Scale") && !(*it).Contains("Lep")) || (*it).Contains("PDF") || (*it).Contains("AlphaS")) )
                            hFiducialPhaseMC_sys_massBinned->Add((TH1*) filein->Get("Acceptance_Efficiency/"+var+ "_" + "/histo_DYJets_" + mass_binned_sample_prefix[i] + "_" + (*it)));
                        else
                            hFiducialPhaseMC_sys_massBinned->Add((TH1*) filein->Get("Acceptance_Efficiency/"+var+ "_" +  "/histo_DYJets_" + mass_binned_sample_prefix[i]));
                    }

                }

                for(int j = istart; j < hFullPhaseMC_massBinned->GetNbinsX()+1; j++)
                {
                    hFiducialPhaseMC_sys->SetBinContent(j, hFiducialPhaseMC_sys_massBinned->GetBinContent(j));
                    hFiducialPhaseMC_sys->SetBinError(j,   hFiducialPhaseMC_sys_massBinned->GetBinError(j));
                }
            }
        }

        // For PDF, AlphaS, Scale etc, nominator (of acceptance) also changes
        if( (((*it).Contains("Scale") && !(*it).Contains("Lep")) || (*it).Contains("PDF") || (*it).Contains("AlphaS")) )
        {
            hFullPhaseMC_raw_sys = (TH1*) filein->Get("Acceptance/"+var+ "_" +  "/histo_DYJets_"+(*it));
            if(year.Contains(2016))
                hFullPhaseMC_raw_sys->Add((TH1*) filein->Get("Acceptance/"+var+ "_" +  "/histo_DYJets10to50_"+(*it)));
            else
                hFullPhaseMC_raw_sys->Add((TH1*) filein->Get("Acceptance/"+var+ "_" + "/histo_DYJets10to50_MG_"+(*it)));

            if(filePath_for_accept != "")
            {
                int istart = binningCoarse->GetGlobalBinNumber(0, 200);
                TH1* hFullPhaseMC_sys_massBinned = NULL;
                for(int i = 0; i < 9; i++)
                {
                    if(i==0)
                    {
                        hFullPhaseMC_sys_massBinned = (TH1*) filein->Get("Acceptance/"+var+ "_" +  "/histo_DYJets_" + mass_binned_sample_prefix[i] + "_" + (*it));
                    }
                    else
                    {
                        hFullPhaseMC_sys_massBinned->Add((TH1*) filein->Get("Acceptance/"+var+ "_" +  "/histo_DYJets_" + mass_binned_sample_prefix[i] + "_" + (*it)));
                    }
                }
                for(int j = istart; j < hFullPhaseMC_massBinned->GetNbinsX()+1; j++)
                {
                    hFullPhaseMC_raw_sys->SetBinContent(j, hFullPhaseMC_sys_massBinned->GetBinContent(j));
                    hFullPhaseMC_raw_sys->SetBinError(j,   hFullPhaseMC_sys_massBinned->GetBinError(j));
                }
            }
        }
        else
        {
            hFullPhaseMC_raw_sys=hFullPhaseMC;
        }

        TH1* hAcceptance_sys = (TH1*) hFullPhaseMC_raw_sys->Clone("hAcceptance_sys");
        hAcceptance_sys->Divide(hFiducialPhaseMC_sys); // systematic acceptance

        // update
        hAcceptance_sys->Add(hAcceptance_raw, -1); // systematic - raw
        hAcceptance_sys->Add(hAcceptance, 1); // nominal + detla


        hSysFullPhaseData[*it]->Multiply(hAcceptance_sys);
        hSysFullPhaseMC[*it] = hFullPhaseMC_raw_sys;

        delete hAcceptance_sys;

        varDir->cd();

        hSysFullPhaseData[*it]->Write();
        hFullPhaseMC_raw_sys->SetName("histo_DY_"+(*it));
        hFullPhaseMC_raw_sys->Write();
        it++;
    }

    topDir->Write();;
    delete hAcceptance_raw;
    delete hFiducialPhaseMC;
}

TH1* ISRUnfold::getUnfoldedHists(TString outHistName, TString steering, bool useAxis)
{
    TH1* outHist = NULL;
    outHist = nominalTUnfold->GetOutput(outHistName,0,0,steering,useAxis);
    return outHist;
}
