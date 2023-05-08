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
void ISRUnfold::setNominalRM(TString file_path, TString top_dir)
{
    TH1::AddDirectory(kFALSE);
    TFile* filein = new TFile(file_path, "READ");

    // bin definition
    TString Rec_binName = "[tunfold:bin]_"+var+"_"+folded_bin_name;
    TString Gen_binName = "[tunfold:bin]_"+var+"_"+unfolded_bin_name;

    // Set bin definition
    binningFine = (TUnfoldBinning*)filein->Get(top_dir + "/" + Rec_binName);
    binningCoarse = (TUnfoldBinning*)filein->Get(top_dir + "/" + Gen_binName);

    // Set response matrix
    // First, get the response matrix
    TH2* hmcGenRec = (TH2*)filein->Get(top_dir + "/[tunfold:matrix]_"+var+"_"+folded_bin_name+"_"+unfolded_bin_name); // hist_suffix

    nominalTUnfold = new TUnfoldDensity(hmcGenRec,
                                        TUnfold::kHistMapOutputHoriz,
                                        regMode,
                                        TUnfold::kEConstraintArea,
                                        densityMode,
                                        binningCoarse,
                                        binningFine);
    hResponseM = (TH2*) hmcGenRec->Clone("[tunfold:matrix]_"+var+"_"+folded_bin_name+"_"+unfolded_bin_name);
    save_hists_from_responseM(filein);

    cout << "TUnfold version: " << nominalTUnfold->GetTUnfoldVersion() << endl;

    filein->Close();
    delete filein;
}

// Option for unfold options
void ISRUnfold::setSystematicRM(TString file_path, TString top_dir, TString sys_name, TString hist_postfix)
{
    TH1::AddDirectory(kFALSE);
    TFile* filein = new TFile(file_path, "READ");

    TH2* hmcGenRec = (TH2*)filein->Get(top_dir + "/" + var + "_responseM" + hist_postfix);

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

void ISRUnfold::save_hists_from_responseM(TFile* file)
{
    fUnfoldOut->cd(channel+year);
    binningCoarse->Write();
    binningFine->Write();

    TH2* truth_unfolded = (TH2*)file->Get(channel+year+"/[tunfold:hist]_"+var+"_"+unfolded_bin_name);
    fUnfoldOut->cd(channel+year+"/DY");
    truth_unfolded->Write("[tunfold:unfolded_hist]_"+var+"_"+unfolded_bin_name);


    TH2F* hMigrationM = (TH2F*) nominalTUnfold->GetProbabilityMatrix("[tunfold:hprob_matrix]_"+var+"_"+folded_bin_name+"_"+unfolded_bin_name);
    hMigrationM->Write();
    hResponseM->Write();

    // Save projection of the migration matrix
    TH1D* hProjectedTruth   = (TH1D*) hResponseM->ProjectionX("[tunfold:projX_hist]_"+var+"_"+unfolded_bin_name, 0, -1, "e"); 
    TH1D* hProjectedBinZero = (TH1D*) hResponseM->ProjectionX("[tunfold:projX_bin_zero_hist]_"+var+"_"+unfolded_bin_name, 0, 0, "e");  
    TH1D* hProjectedReco    = (TH1D*) hResponseM->ProjectionY("[tunfold:projY_hist]_"+var+"_"+folded_bin_name, 1, -1, "e");
    hProjectedTruth->Write();
    hProjectedBinZero->Write();
    hProjectedReco->Write();
}

// Set input histogram from root file
void ISRUnfold::setUnfInput(TString file_path, TString top_dir, TString sys_type, TString sys_name, TString sys_hist_postfix)
{
    TH1::AddDirectory(kFALSE);

    TFile* filein = new TFile(file_path);

    TString full_hist_path = top_dir+"/[tunfold:hist]_"+var+"_"+folded_bin_name+sys_hist_postfix;
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

void ISRUnfold::subBkgs(TString file_path, TString top_dir, TString bkg_name, TString sys_type, TString sys_name, TString hist_postfix)
{
    TFile* filein = new TFile(file_path);
    TH1* hRec = (TH1*)filein->Get(top_dir+"/[tunfold:hist]_"+var+"_"+folded_bin_name+hist_postfix);

    // Nominal histograms
    if(sys_type=="Type_0")
    {
        nominalTUnfold->SubtractBackground(hRec, bkg_name);

        // Save Unfolding fake histogram in folded directory
        if(top_dir.Contains("Fake"))
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

    // No regularisation
    if(regMode == TUnfold::kRegModeNone)
    {
        if (partialReg) {

            //setPartialRegularize2D(TUnfold::kRegModeCurvature, 200., 0., 1000., 100.);
            setPartialRegularize2D(TUnfold::kRegModeCurvature, 320., 0., 1000., 1000.);
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

    bool useAxisBinning = false;
    if(var.Contains("1D")) // FIXME check bin structure to decide whether to use axis_binning or not
    {
        useAxisBinning = true;
    }

    fUnfoldOut->cd(channel+year+"/DY"); 
    nominalTUnfold->GetRhoIJtotal("[tunfold:corr_hist]_"+var+"_"+folded_bin_name+"_"+unfolded_bin_name, 0, 0, 0, useAxisBinning)->Write();
    nominalTUnfold->GetEmatrixTotal("[tunfold:cov_hist]_"+var+"_"+folded_bin_name+"_"+unfolded_bin_name, 0, 0, 0, useAxisBinning)->Write();

    fUnfoldOut->cd(channel+year+"/data"); 
    nominalTUnfold->GetOutput("[tunfold:unfolded_hist]_"+var+"_"+unfolded_bin_name,0,0, "*[*]", useAxisBinning)->Write();
    nominalTUnfold->GetInput("[tunfold:input_hist]_"+var+"_"+folded_bin_name, 0, 0, 0, false)->Write();

    // For systematic
    std::vector<TString>::iterator it = sysVector.begin();
    while(it != sysVector.end())
    {
        if( (*it).Contains("IterEM"))
        {
            iBest=iterEMTUnfold->ScanSURE(NITER_Iterative, &graph_SURE_IterativeSURE, &graph_DFdeviance_IterativeSURE);
            //cout << "iBest: " << iBest << endl;

            iterEMTUnfold->GetOutput("histo_Data_"+(*it),0,0, "*[*]", false)->Write();
        }
        else{
            systematicTUnfold[*it]->DoUnfold(tau);

            systematicTUnfold[*it]->GetOutput("histo_Data_"+(*it),0,0, "*[*]", false)->Write();
        }
        it++;
    }
}

// In this function, acceptance factors change only for PDF, Scale, AlphaS systematics
void ISRUnfold::doAcceptCorr(TString filePath, TString top_dir)
{
    TFile* filein = new TFile(filePath);


    TH1* hFiducialPhaseMC = NULL;

    hFullPhaseMC     = (TH1*) filein->Get(top_dir+"/[tunfold:fullphase_hist]_"+var+"_"+unfolded_bin_name);
    hFiducialPhaseMC = (TH1*) fUnfoldOut->Get(channel+year+"/DY/[tunfold:projX_hist]_"+var+"_"+unfolded_bin_name);

    hAcceptance = (TH1*) hFullPhaseMC->Clone("[tunfold:acc_hist]_"+var+"_"+unfolded_bin_name);
    hAcceptance->Divide(hFiducialPhaseMC); // Nominal acceptance correction factor

    TH1* hAcceptance_raw = (TH1*) hAcceptance->Clone("hAcceptance_raw");

    fUnfoldOut->cd(channel+year+"/data"); 
    hFullPhaseData = nominalTUnfold->GetOutput("[tunfold:unfolded_fullphase_hist]_"+var+"_"+unfolded_bin_name,0,0, "*[*]", false);
    hFullPhaseData->Multiply(hAcceptance); // acceptance corrected data
    hFullPhaseData->Write();

    fUnfoldOut->cd(channel+year+"/DY"); 
    hFullPhaseMC->SetName("[tunfold:unfolded_fullphase_hist]_"+var+"_"+unfolded_bin_name);
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
        }

        // For PDF, AlphaS, Scale etc, nominator (of acceptance) also changes
        if( (((*it).Contains("Scale") && !(*it).Contains("Lep")) || (*it).Contains("PDF") || (*it).Contains("AlphaS")) )
        {
            hFullPhaseMC_raw_sys = (TH1*) filein->Get("Acceptance/"+var+ "_" +  "/histo_DYJets_"+(*it));
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

        hSysFullPhaseData[*it]->Write();
        hFullPhaseMC_raw_sys->SetName("histo_DY_"+(*it));
        hFullPhaseMC_raw_sys->Write();
        it++;
    }

    delete hAcceptance_raw;
    delete hFiducialPhaseMC;
}

TH1* ISRUnfold::getUnfoldedHists(TString outHistName, TString steering, bool useAxis)
{
    TH1* outHist = NULL;
    outHist = nominalTUnfold->GetOutput(outHistName,0,0,steering,useAxis);
    return outHist;
}
