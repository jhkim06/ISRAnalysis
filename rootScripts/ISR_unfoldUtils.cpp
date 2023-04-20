#include "ISR_unfoldUtils.h"

void ISRUnfold::setBias(double bias)
{
   nominalBias = bias;
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

void ISRUnfold::checkMatrixCond()
{
    TH2F* hProb = NULL;
    hProb = (TH2F*) nominalTUnfold->GetProbabilityMatrix("hProb_"+var);

    //int nBinsX = hProb->GetNbinsX();
    //int nBinsY = hProb->GetNbinsY();
    TMatrixD matrix = makeMatrixFromHist(hProb);

    TDecompSVD decomp(matrix);
    double condition = decomp.Condition();
    cout << "The condition number for " << ": " << condition << endl;

    double determinant;
    TMatrixD mInverse = matrix.Invert(&determinant);
    cout << "The determinant of " << " is " << determinant << endl;
    //return condition;

/*
    // TODO check hProb->GetNbinsY()+2,hProb->GetNbinsX()+2 without "F"
    TMatrixD* mProb=new TMatrixD(hProb->GetNbinsX()+2,hProb->GetNbinsY()+2,hProb->GetArray(),"F"); // +2 for under/overflow bins
    TMatrixD* mtProb=new TMatrixD(hProb->GetNbinsY()+2,hProb->GetNbinsX()+2);
    mtProb->Transpose(*mProb);
    TDecompSVD* svdProb=new TDecompSVD(*mtProb);
    //mtProb->Print();
    cout << "Decompose(), successed? " << svdProb->Decompose() << endl;
    const Int_t colLwb = svdProb->GetColLwb();
    const Int_t nCols  = svdProb->GetNcols();
    const TVectorD& singularValues = svdProb->GetSig();

    // Find minimum (>0)
    int i = colLwb+nCols-1;
    while(i >= colLwb)
    {
        if(singularValues[i] > 0)
        {
            double min = singularValues[i];
            cout << var << ", Cond(): " << singularValues[colLwb]/min << endl;
            break;
        }
        else
        {
           i--;
        }
    }
    return svdProb->GetSig();
*/
}

void ISRUnfold::checkIterEMUnfold(void) 
{
    double yMin=1.;
    double yLine=10.;
    double yMax=graph_SURE_IterativeSURE->GetMaximum();

    gStyle->SetPadBottomMargin(0.2);
    TCanvas *canvas1=new TCanvas("compare","",3600,1200);
    canvas1->Divide(2,1);

    canvas1->cd(1);
    gPad->SetLogy();
    graph_SURE_IterativeSURE->GetYaxis()->SetRangeUser(yMin,yMax);
    graph_SURE_IterativeSURE->GetXaxis()->SetRangeUser(-1.5,100.5);
    graph_SURE_IterativeSURE->GetXaxis()->SetTitle("iteration");
    graph_SURE_IterativeSURE->GetXaxis()->SetTitleOffset(1.2);
    graph_SURE_IterativeSURE->GetXaxis()->SetTitleFont(43);
    graph_SURE_IterativeSURE->GetXaxis()->SetTitleSize(100);
    graph_SURE_IterativeSURE->SetMarkerColor(kBlue);
    graph_SURE_IterativeSURE->SetMarkerStyle(20);
    graph_SURE_IterativeSURE->SetMarkerSize(2);
    graph_SURE_IterativeSURE->DrawClone("APW");
    int n_scanSURE_iterative=graph_SURE_IterativeSURE->GetN();
    double const *nIter_scanSURE_iterative=graph_SURE_IterativeSURE->GetX();
    double const *DF_scanSURE_iterative=graph_DFdeviance_IterativeSURE->GetX();
    double const *deviance_scanSURE=graph_DFdeviance_IterativeSURE->GetY();
    TGraph *DF_iterative=new TGraph
       (n_scanSURE_iterative,nIter_scanSURE_iterative,DF_scanSURE_iterative);
    TGraph *deviance_iterative=new TGraph
       (n_scanSURE_iterative,nIter_scanSURE_iterative,deviance_scanSURE);
    DF_iterative->SetMarkerColor(kRed);
    DF_iterative->SetMarkerStyle(24);
    DF_iterative->SetMarkerSize(2);
    DF_iterative->DrawClone("P");
    deviance_iterative->SetMarkerColor(kMagenta);
    deviance_iterative->SetMarkerStyle(22);
    deviance_iterative->SetMarkerSize(2);
    deviance_iterative->DrawClone("P");
    TLine *line2=new TLine(iBest,yLine,iBest,1e4);
    line2->SetLineStyle(1);
    line2->Draw();
    TLegend *legend3=new TLegend(0.25,0.2,0.9,0.45,"Iterative EM, minimize SURE");
    legend3->SetBorderSize(0);
    legend3->SetFillStyle(0);
    legend3->SetTextSize(0.045);
    legend3->AddEntry(graph_SURE_IterativeSURE,"SURE","p");
    legend3->AddEntry(DF_iterative,"D.F.","p");
    legend3->AddEntry(deviance_iterative,"deviance","p");

    legend3->AddEntry(line2,TString::Format
                      ("min(SURE) at iteration=%d",iBest),"l");
    legend3->AddEntry((TObject *)0,TString::Format
                      ("D.F.=%3g",DF_scanSURE_iterative[iBest]),"");
    legend3->Draw();
    canvas1->SaveAs("ISR_scan.pdf");
}

// Set the nominal response matrix
void ISRUnfold::setNominalRM(TString filepath, TString dirName, TString binDef)
{
    TH1::AddDirectory(kFALSE);
    TFile* filein = new TFile(filepath, "READ");

    TString fullDirPath = dirName + "/";

    // bin definition
    TString Rec_binName = var + "_smeared_bin"; 
    TString Gen_binName = var + "_truth_bin";
    Rec_binName = fullDirPath + Rec_binName;
    Gen_binName = fullDirPath + Gen_binName;

    // Set bin definition
    binningFine = (TUnfoldBinning*)filein->Get(Rec_binName);
    binningCoarse = (TUnfoldBinning*)filein->Get(Gen_binName);

    // Set response matrix
    // First, get the response matrix
    TH2* hmcGenRec = (TH2*)filein->Get(fullDirPath + var + "_responseM");

    nominalTUnfold = new TUnfoldDensity(hmcGenRec, 
                                        TUnfold::kHistMapOutputHoriz, 
                                        regMode, 
                                        TUnfold::kEConstraintArea, 
                                        densityMode, 
                                        binningCoarse, 
                                        binningFine);
    hResponseM = (TH2*) hmcGenRec->Clone("hResponseM");

    cout << "TUnfold version: " << nominalTUnfold->GetTUnfoldVersion() << endl;

    // For statistical uncertainty
    if(doInputStatUnc)
    {
        // cout << "Create response matrix for statistical uncertainty..." << endl;
        for(int i = 0; i < statSize; i++)
        {
            unfInputStatTUnfold.push_back(new TUnfoldDensity(hmcGenRec, TUnfold::kHistMapOutputHoriz, regMode, TUnfold::kEConstraintArea, densityMode, binningCoarse, binningFine));
        }
    }
    if(doRMStatUnc)
    {
        for(int i = 0; i < statSize; i++)
        {
            TString nth_;
            nth_.Form("%d", i);
            TH2* tempRM = (TH2*) hmcGenRec->Clone("hRM_stat_" + nth_);

            for(int xbin=1; xbin <= tempRM->GetXaxis()->GetNbins(); xbin++)
            {
                for(int ybin=0; ybin <= tempRM->GetYaxis()->GetNbins(); ybin++)
                {
                    double err = tempRM->GetBinError(xbin, ybin);
                    if(err >= 0.0)
                        tempRM->SetBinContent(xbin, ybin, tempRM->GetBinContent(xbin, ybin) + gRandom->Gaus(0, err));
                }
            }

            unfMatrixStatTUnfold.push_back(new TUnfoldDensity(tempRM, TUnfold::kHistMapOutputHoriz, regMode, TUnfold::kEConstraintArea, densityMode, binningCoarse, binningFine));
        }
    }

    TH2* hmcGenRec_wo_binZero = NULL;
    if(ignoreBinZero)
    {
        hmcGenRec_wo_binZero = (TH2*) hmcGenRec->Clone("matrix_wo_binzero");
        for(int xbin=0; xbin <= hmcGenRec_wo_binZero->GetXaxis()->GetNbins(); xbin++)
        {
            hmcGenRec_wo_binZero->SetBinContent(xbin, 0, 0);
            hmcGenRec_wo_binZero->SetBinError(xbin, 0, 0);
        }
        ignoreBinZeroTUnfold = new TUnfoldDensity(hmcGenRec_wo_binZero, TUnfold::kHistMapOutputHoriz, regMode, TUnfold::kEConstraintArea, densityMode, binningCoarse, binningFine);
    }

    if(doModelUnc)
    {
        modelUncertaintyTUnfold = new TUnfoldDensity(hReweightSF, TUnfold::kHistMapOutputHoriz, regMode, TUnfold::kEConstraintArea, densityMode, binningCoarse, binningFine);
    }

    // TODO make function for saving histograms
    // Save migration and response matrix
    TDirectory* varDir;
    TDirectory* varDirForReco;

    varDirForReco=fUnfoldOut->GetDirectory("folded/"+var);
    varDir=fUnfoldOut->GetDirectory("matrix/"+var);
    varDir->cd();
    binningCoarse->Write();
    binningFine->Write();

    TH2F* hMigrationM = (TH2F*) nominalTUnfold->GetProbabilityMatrix("hMigrationM");
    hMigrationM->Write();
    hmcGenRec->SetName("hMigrationM");
    hmcGenRec->Write();

    // Save projection of the migration matrix
    varDir=fUnfoldOut->GetDirectory("unfolded/"+var);
    varDir->cd();

    TH1D* hProjectedTruth = (TH1D*) hmcGenRec->ProjectionX("histo_DY", 0, -1, "e");  //
    TH1D* hProjectedBinZero = (TH1D*) hmcGenRec->ProjectionX("histo_ProjectedBinZero", 0, 0, "e");  //

    TH1D* hProjectedReco = (TH1D*) hmcGenRec->ProjectionY("histo_ProjectedReco", 1, -1, "e");  //
    TH1D* hProjectedReweightedReco = NULL;
    if(doModelUnc)
    {
        hProjectedReweightedReco = hReweightSF->ProjectionY("histo_ProjectedReweightedReco", 1, -1, "e"); 
    }
    
    if(ignoreBinZero)
    {
        TH1D* hProjectedBinZeroWoBinZero = (TH1D*)hmcGenRec_wo_binZero->ProjectionX("histo_ProjectedBinZeroWoBinZero", 0, 0, "e");
        TH1D* hProjectedTruthWoBinZero = (TH1D*)hmcGenRec_wo_binZero->ProjectionX("histo_ProjectedTruthWoBinZero", 0, -1, "e");
        TH1D* hProjectedRecoWoBinZero = (TH1D*)hmcGenRec_wo_binZero->ProjectionY("histo_ProjectedRecoWoBinZero", 1, -1, "e");

        hProjectedBinZeroWoBinZero->Write();
        hProjectedTruthWoBinZero->Write();

        varDirForReco->cd();
        hProjectedRecoWoBinZero->Write();
        varDir->cd();
    }

    hProjectedTruth->Write();
    hProjectedBinZero->Write();

    varDirForReco->cd();
    binningFine->Write();
    hProjectedReco->Write();
    if(doModelUnc)
    {
        hProjectedReweightedReco->Write();
    }
    varDir->cd();

    filein->Close();
    delete filein;
}

void ISRUnfold::setFromPrevUnfResult(ISRUnfold* unfold, bool useAccept)
{
    TDirectory* varDir;

    varDir=fUnfoldOut->GetDirectory("unfolded/"+var);
    varDir->cd();
    //cout << "setFromPrevUnfResult(), useAccept? " << useAccept << endl;
    // Loop over sytematics considered in the previous unfold class
    // So first get sysVector map object
    std::vector<TString> sysVector_previous = unfold->getSystematicVector();
    std::vector<TString>::iterator it = sysVector_previous.begin();
    while(it != sysVector_previous.end())
    {
        //cout << "Systematic name: " << it->first << endl;
        std::vector<TString>::iterator found = find(this->sysVector.begin(), this->sysVector.end(), *it);
        if(found == this->sysVector.end())
        //if(this->sysVector.find(it->first) == this->sysVector.end())
        {
            // Not found in this ISRUnfold class, but exits in the previous one
            // Create TUnfoldDensity using the DEFAULT response matrix
            //cout << "Systematic variation, " << sysVector_previous[it->first][ith] << endl;
            if((*it).Contains("IterEM"))
            {
                this->iterEMTUnfold   = new TUnfoldIterativeEM(hResponseM,TUnfoldDensity::kHistMapOutputHoriz,binningCoarse,binningFine);

                if(!useAccept)
                {
                    this->iterEMTUnfold->SetInput(unfold->iterEMTUnfold->GetOutput("hUnfolded" + var + "_"+ *it + "_" + *it,0,0,"*[*]",false), nominalBias);
                }
                else
                {
                    //cout << "use acceptance corrected output!" << endl;
                    this->iterEMTUnfold->SetInput(unfold->hSysFullPhaseData[*it], nominalBias);
                }
            }
            else
            {
                this->systematicTUnfold[*it] = new TUnfoldDensity(hResponseM,TUnfold::kHistMapOutputHoriz,regMode, TUnfold::kEConstraintArea, densityMode, binningCoarse,binningFine);
                TH1D* hProjectedTruth = (TH1D*) hResponseM->ProjectionX("histo_DY_"+(*it), 0, -1, "e");
                hProjectedTruth->Write();

                if(!useAccept)
                {
                     this->systematicTUnfold[*it]->SetInput(unfold->systematicTUnfold[*it]->GetOutput("hUnfolded" + var + "_"+ *it + "_" + *it,0,0,"*[*]",false), nominalBias);
                }
                else
                {
                    //cout << "use acceptance corrected output!" << endl;
                    this->systematicTUnfold[*it]->SetInput(unfold->hSysFullPhaseData[*it], nominalBias);
                }
            }
            this->sysVector.push_back(*it);
        }
        it++;
    }
    // Loop over variations of event selection efficiency correction
}

// Option for unfold options
void ISRUnfold::setSystematicRM(TString filepath, TString dirName, TString binDef, TString sysName, TString histPostfix)
{
    TFile* filein = new TFile(filepath, "READ");
    TH2* hmcGenRec = NULL;

    TString histNameWithSystematic = "hmc" + var + "GenRec" + histPostfix;
    hmcGenRec = (TH2*)filein->Get(dirName + "/" + var + "_ResMatrix_" + binDef + "/" + histNameWithSystematic);
    //cout << "ISRUnfold::setSystematicRM " << filepath << " " << dirName + "/" + var + "_ResMatrix_" + binDef + "/" + histNameWithSystematic << endl;

    TDirectory* varDir;

    varDir=fUnfoldOut->GetDirectory("unfolded/"+var);
    varDir->cd();

    if(sysName.Contains("IterEM"))
    {
        iterEMTUnfold = new TUnfoldIterativeEM(hmcGenRec,TUnfoldDensity::kHistMapOutputHoriz,binningCoarse,binningFine);
    }
    else
    {
        //cout << sysName << " crate TUnfoldDensity..." << endl;
        if(hmcGenRec == NULL) cout << "check input file" << endl; 
        systematicTUnfold[sysName] = new TUnfoldDensity(hmcGenRec, TUnfold::kHistMapOutputHoriz, regMode, TUnfold::kEConstraintArea, densityMode, binningCoarse, binningFine);
    }
    TH1D* hProjectedTruth = (TH1D*) hmcGenRec->ProjectionX("histo_DY_"+sysName, 0, -1, "e");  //
    hProjectedTruth->Write();

    filein->Close();
    delete filein;
    delete hProjectedTruth;
}

// Set input histogram using the nominal output of the previous unfolding
void ISRUnfold::setUnfInput(ISRUnfold* unfold, TString thisSysType, TString sysName, bool useAccept)
{
    TH1::AddDirectory(kFALSE);

    if(!useAccept)
    {
        if(thisSysType=="Type_0")
        {
            nominalTUnfold->SetInput(unfold->getUnfoldedHists("UnfoldOut_"+var, "*[*]", false), nominalBias);
        }
        else
        {
            // FIXME
            systematicTUnfold[sysName]->SetInput(unfold->getUnfoldedHists("UnfoldOut_"+var+thisSysType+sysName, "*[*]", false), nominalBias);
        }
    }
    else
    {
        
        //cout << "set from previous unfold class, isSys " << isSys << endl;
        if(thisSysType=="Type_0")
        {
            nominalTUnfold->SetInput(unfold->hFullPhaseData, nominalBias);
        }
        else
        {
            //
            if(sysName.Contains("IterEM"))
            {
                iterEMTUnfold->SetInput(unfold->hFullPhaseData, 1.);
            }
            else
            {
                if(sysName.Contains("fsr"))
                    systematicTUnfold[sysName]->SetInput(unfold->hFullPhaseData, nominalBias);
                else
                {
                    TH1* htemp = unfold->hSysFullPhaseData[sysName];
                    systematicTUnfold[sysName]->SetInput(htemp, nominalBias);
                }
            }
        }
    }
}

void ISRUnfold::setUnfInputUnfSys()
{
    //
    // input statistical
    if(doInputStatUnc)
    {
        for(int istat = 0; istat < statSize; istat++)
        {
            TH1* tempInput;

            TString nth_;
            nth_.Form("%d", istat);
            tempInput = nominalTUnfold->GetInput("temp" + var + "Hist_" + nth_, 0, 0, 0, false);

            // randomize input histogram bin content
            for(int ibin = 1; ibin<tempInput->GetNbinsX()+1;ibin++)
            {
                double err = tempInput->GetBinError(ibin);
                if(err > 0.0)
                {
                    tempInput->SetBinContent(ibin, tempInput->GetBinContent(ibin) + gRandom->Gaus(0,err));
                }
            }
            unfInputStatTUnfold.at(istat)->SetInput(tempInput, nominalBias);
            delete tempInput;
        }
    }
    
    // response matrix statistical
    TH1* tempInput = nominalTUnfold->GetInput("BkgSubtractedInput", 0, 0, 0, false);

    if(doRMStatUnc)
    {
        for(int istat = 0; istat < statSize; istat++)
        {
            TString nth_;
            nth_.Form("%d", istat);

            unfMatrixStatTUnfold.at(istat)->SetInput(tempInput, nominalBias);
        }
    }
    
    if(ignoreBinZero)
    {
        ignoreBinZeroTUnfold->SetInput(tempInput, nominalBias);
    }

    if(doModelUnc)
    {
        modelUncertaintyTUnfold->SetInput(tempInput, nominalBias);
    }
    
    delete tempInput;
    
}

// Set input histogram from root file
void ISRUnfold::setUnfInput(TString filepath, TString dirName, TString binDef, TString sysType, TString sysName, TString histPostfix, bool isFSR)
{
    TH1::AddDirectory(kFALSE);

    TFile* filein = new TFile(filepath);
    TH1* hRec = (TH1*)filein->Get(dirName + "/" + var + "_smeared"); 
    
    // Nominal
    if(sysType == "Type_0")
    {
        nominalTUnfold->SetInput(hRec, nominalBias, nominalBias);
    }
    else
    {
        // Systematic histograms
        if(sysName.Contains("IterEM"))
        {
            iterEMTUnfold->SetInput(hRec, 1.);
        }
        else
        {
            systematicTUnfold[sysName]->SetInput(hRec, nominalBias);
        }
    }

    filein->Close();
    delete filein;
}

void ISRUnfold::subBkgs(TString filepath, TString dirName, TString binDef, TString bkgName, TString sysType, TString sysName, TString histPostfix)
{
    TFile* filein = new TFile(filepath);
    TH1* hRec = NULL;

    hRec = (TH1*)filein->Get(dirName + "/" + var + "_smeared");

    // Nominal histograms
    if(sysType=="Type_0")
    {
        nominalTUnfold->SubtractBackground(hRec, bkgName);
        
        // Save Unfolding fake histogram in folded directory
        if(dirName.Contains("Fake"))
        {
            TDirectory* varDirForReco;

            varDirForReco=fUnfoldOut->GetDirectory("folded/"+var);
            varDirForReco->cd();

            hRec->SetName("histo_Fake"+bkgName);
            hRec->Write();
        }
    }
    else
    // Systematic
    {
        //cout << "file path: " << filepath << endl;
        if(sysName.Contains("IterEM"))
        {
            iterEMTUnfold->SubtractBackground(hRec, bkgName);
        }
        else
        {
            // FIXME temporary method for background systematic
            if(sysName.Contains("Background"))
            {
                if(sysName.Contains("Up"))
                {
                    systematicTUnfold[sysName]->SubtractBackground(hRec, bkgName, 1.05);
                }
                if(sysName.Contains("Down"))
                {
                    systematicTUnfold[sysName]->SubtractBackground(hRec, bkgName, 0.95);
                }
            }
            else
            {
                systematicTUnfold[sysName]->SubtractBackground(hRec, bkgName);
            }
        }
    }

    filein->Close();
    delete filein;
}

void ISRUnfold::setSystematics(TString sysHistName)
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

    if(doInputStatUnc)
    {
        for(int istat = 0; istat < statSize; istat++)
        {
            TString nth_;
            nth_.Form("%d", istat);
            
            unfInputStatTUnfold.at(istat)->RegularizeBins(istart, 1, iend-istart+1, partialRegMode);
        }
    }

    if(doRMStatUnc)
    {
        for(int istat = 0; istat < statSize; istat++)
        {
            TString nth_;
            nth_.Form("%d", istat);

            unfMatrixStatTUnfold.at(istat)->RegularizeBins(istart, 1, iend-istart+1, partialRegMode);
        }
    }

    if(ignoreBinZero)
    {
        ignoreBinZeroTUnfold->RegularizeBins(istart, 1, iend-istart+1, partialRegMode);;
    }

    if(doModelUnc)
    {
        modelUncertaintyTUnfold->RegularizeBins(istart, 1, iend-istart+1, partialRegMode);;
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

    if(doInputStatUnc)
    {
        for(int istat = 0; istat < statSize; istat++)
        {
            TString nth_;
            nth_.Form("%d", istat);
            
            unfInputStatTUnfold.at(istat)->DoUnfold(tau);

            varDir->cd();
            unfInputStatTUnfold.at(istat)->GetOutput("histo_Data_UnfoldingInputStat_" + nth_, 0, 0, "*[*]", false)->Write();

            delete unfInputStatTUnfold.at(istat);
        }
        unfInputStatTUnfold.clear();
    }

    if(doRMStatUnc)
    {
        for(int istat = 0; istat < statSize; istat++)
        {
            TString nth_;
            nth_.Form("%d", istat);

            unfMatrixStatTUnfold.at(istat)->DoUnfold(tau);
            varDir->cd();
            unfMatrixStatTUnfold.at(istat)->GetOutput("histo_Data_UnfoldingMatrixStat_" + nth_, 0, 0, "*[*]", false)->Write();

            delete unfMatrixStatTUnfold.at(istat);
        }
        unfMatrixStatTUnfold.clear();
    }

    if(ignoreBinZero)
    {
        ignoreBinZeroTUnfold->DoUnfold(tau);

        varDir->cd();
        ignoreBinZeroTUnfold->GetOutput("histo_DataNoBinZeroTUnfold", 0, 0, "*[*]", false)->Write();
    }

    if(doModelUnc)
    {
        TH1D* temp_projectedTruthReweighted = (TH1D*) hReweightSF->ProjectionX("histo_DY_reweighted", 0, -1, "e");  //
        
        modelUncertaintyTUnfold->DoUnfold(tau);

        varDir->cd();
        modelUncertaintyTUnfold->GetOutput("histo_Data_UnfoldModel", 0, 0, "*[*]", false)->Write();
        temp_projectedTruthReweighted->Write();

        delete temp_projectedTruthReweighted;
    }
    
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
void ISRUnfold::doAcceptCorr(TString filePath, TString binDef, TString filePath_for_accept)
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

    hFullPhaseMC = (TH1*) filein->Get("Acceptance/"+var+ "_" + binDef + "/histo_DYJets");
    if(year==2016)
        hFullPhaseMC->Add((TH1*) filein->Get("Acceptance/"+var+ "_" + binDef + "/histo_DYJets10to50"));
    else
        hFullPhaseMC->Add((TH1*) filein->Get("Acceptance/"+var+ "_" + binDef + "/histo_DYJets10to50_MG"));

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
                hFullPhaseMC_massBinned = (TH1*) filein->Get("Acceptance/"+var+ "_" + binDef + "/histo_DYJets_" + mass_binned_sample_prefix[i]);
            }
            else
            {
                hFullPhaseMC_massBinned->Add((TH1*) filein->Get("Acceptance/"+var+ "_" + binDef + "/histo_DYJets_" + mass_binned_sample_prefix[i]));
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
                hFiducialPhaseMC_massBinned = (TH1*) filein->Get("Acceptance_Efficiency/"+var+ "_" + binDef + "/histo_DYJets_" + mass_binned_sample_prefix[i]);
            }
            else
            {
                hFiducialPhaseMC_massBinned->Add((TH1*) filein->Get("Acceptance_Efficiency/"+var+ "_" + binDef + "/histo_DYJets_" + mass_binned_sample_prefix[i]));
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

    if(doInputStatUnc)
    {
        TH1::AddDirectory(kFALSE);
        TH1* hFullPhaseDataTemp = NULL;
        for(int istat = 0; istat < statSize; istat++)
        {
            TString nth_;
            nth_.Form("%d", istat);

            hFullPhaseDataTemp=(TH1*)fUnfoldOut->Get("unfolded/"+var+"/"+"histo_Data_UnfoldingInputStat_" + nth_);
            hFullPhaseDataTemp->Multiply(hAcceptance);

            varDir->cd();
            hFullPhaseDataTemp->Write();

            delete hFullPhaseDataTemp;
        }
    }

    if(doRMStatUnc)
    {
        TH1::AddDirectory(kFALSE);
        TH1* hFullPhaseDataTemp = NULL;
        for(int istat = 0; istat < statSize; istat++)
        {
            TString nth_;
            nth_.Form("%d", istat);

            //hFullPhaseDataTemp=unfMatrixStatTUnfold.at(istat)->GetOutput("histo_Data_UnfoldingMatrixStat_" + nth_, 0, 0, "*[*]", false);
            hFullPhaseDataTemp=(TH1*)fUnfoldOut->Get("unfolded/"+var+"/"+"histo_Data_UnfoldingMatrixStat_" + nth_);
            hFullPhaseDataTemp->Multiply(hAcceptance);

            varDir->cd();
            hFullPhaseDataTemp->Write();

            delete hFullPhaseDataTemp;
        }
    }

    if(doModelUnc)
    {
        TH1::AddDirectory(kFALSE);
        TH1* hFullPhaseDataTemp = NULL;

        hFullPhaseDataTemp=(TH1*)fUnfoldOut->Get("unfolded/"+var+"/"+"histo_Data_UnfoldModel");   
        hFullPhaseDataTemp->Multiply(hAcceptance);

        varDir->cd(); 
        hFullPhaseDataTemp->Write(); 
        delete hFullPhaseDataTemp; 
    }

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
                            hFiducialPhaseMC_sys_massBinned = (TH1*) filein->Get("Acceptance_Efficiency/"+var+ "_" + binDef + "/histo_DYJets_" + mass_binned_sample_prefix[i] + "_" + (*it));
                        else 
                            hFiducialPhaseMC_sys_massBinned = (TH1*) filein->Get("Acceptance_Efficiency/"+var+ "_" + binDef + "/histo_DYJets_" + mass_binned_sample_prefix[i]);
                    }
                    else
                    {
                        if( (((*it).Contains("Scale") && !(*it).Contains("Lep")) || (*it).Contains("PDF") || (*it).Contains("AlphaS")) ) 
                            hFiducialPhaseMC_sys_massBinned->Add((TH1*) filein->Get("Acceptance_Efficiency/"+var+ "_" + binDef + "/histo_DYJets_" + mass_binned_sample_prefix[i] + "_" + (*it)));
                        else
                            hFiducialPhaseMC_sys_massBinned->Add((TH1*) filein->Get("Acceptance_Efficiency/"+var+ "_" + binDef + "/histo_DYJets_" + mass_binned_sample_prefix[i]));
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
            hFullPhaseMC_raw_sys = (TH1*) filein->Get("Acceptance/"+var+ "_" + binDef + "/histo_DYJets_"+(*it));
            if(year==2016)
                hFullPhaseMC_raw_sys->Add((TH1*) filein->Get("Acceptance/"+var+ "_" + binDef + "/histo_DYJets10to50_"+(*it)));
            else
                hFullPhaseMC_raw_sys->Add((TH1*) filein->Get("Acceptance/"+var+ "_" + binDef + "/histo_DYJets10to50_MG_"+(*it)));

            if(filePath_for_accept != "")
            {
                int istart = binningCoarse->GetGlobalBinNumber(0, 200);
                TH1* hFullPhaseMC_sys_massBinned = NULL;
                for(int i = 0; i < 9; i++)
                {
                    if(i==0)
                    {
                        hFullPhaseMC_sys_massBinned = (TH1*) filein->Get("Acceptance/"+var+ "_" + binDef + "/histo_DYJets_" + mass_binned_sample_prefix[i] + "_" + (*it));
                    }
                    else
                    {
                        hFullPhaseMC_sys_massBinned->Add((TH1*) filein->Get("Acceptance/"+var+ "_" + binDef + "/histo_DYJets_" + mass_binned_sample_prefix[i] + "_" + (*it)));
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

void ISRUnfold::varyHistWithStatError(TH1* hist, int sys)
{
    for(int ibin = 1; ibin < hist->GetNbinsX()+1; ibin++)
    {
        hist->SetBinContent(ibin, hist->GetBinContent(ibin) + double(sys) * hist->GetBinError(ibin));
    }
}

TH1* ISRUnfold::getUnfoldedHists(TString outHistName, TString steering, bool useAxis)
{
    TH1* outHist = NULL;
    outHist = nominalTUnfold->GetOutput(outHistName,0,0,steering,useAxis);
    return outHist;
}

TH1* ISRUnfold::getRawHist(TString filePath, TString dirName, TString histName, TString outHistName, TString steering, bool useAxis)
{
    TH1::AddDirectory(kFALSE);
    TFile* filein = new TFile(filePath);
    TH1* hist = NULL;

    if(steering != "")
    {
        TH1* raw_hist = (TH1*)filein->Get(dirName+"/"+var+"/"+histName);
        if(histName.Contains("DYJetsTo") && !histName.Contains("Tau"))
        {
            histName.ReplaceAll("DYJetsTo", "DYJets10to50To");
            raw_hist->Add((TH1*)filein->Get(dirName+"/"+var+"/"+histName));
        }

        hist = binningFine->ExtractHistogram(outHistName, raw_hist, 0, useAxis, steering);

        delete raw_hist;
    }
    else
    {
        //cout << "Steering not provided, get raw histogram." << endl;
        //cout << dirName+"/"+var+"/"+histName << endl;
        hist = (TH1*)filein->Get(dirName+"/"+var+"/"+histName);
        if(histName.Contains("DYJetsTo") && !histName.Contains("Tau"))
        {
            histName.ReplaceAll("DYJetsTo", "DYJets10to50To");
            hist->Add((TH1*)filein->Get(dirName+"/"+var+"/"+histName));
        }
    }

    delete filein;
    return hist;
}

