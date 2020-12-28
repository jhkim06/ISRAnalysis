#include "ISR_unfoldUtils.h"

void ISRUnfold::setBias(double bias)
{
   nominal_bias = bias;
}

void ISRUnfold::setOutputBaseDir(TString outPath)
{
    output_baseDir = outPath;
}

const TVectorD& ISRUnfold::checkMatrixCond()
{
    TH2D* hProb = NULL;
    hProb = (TH2D*) nominalTUnfold->GetProbabilityMatrix("hProb_"+var);

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
}

double ISRUnfold::getSmearedChi2(TString filePath, TString dirName, TString steering, bool useAxis, bool divBinWidth)
{
    double chi2 = 0.;
    double ndf  = 0.;

    TH1* hData; // Data - Bkg
    TH1* hDY; // DY MC

    TString DYHistName_ = "histo_DYJetsToMuMu";
    if(channel_name == "electron") DYHistName_ = "histo_DYJetsToEE";
    if(var.Contains("Mass"))
    {
        hData = nomMassUnfold->GetInput("hData_"+var, 0, 0, steering, useAxis);
        hDY = getRawHist(filePath, dirName, DYHistName_, "Signal_"+var, steering, useAxis, divBinWidth); ;
    }
    else
    {
        hData = nominalTUnfold->GetInput("hData_"+var, 0, 0, steering, useAxis);
        hDY = getRawHist(filePath, dirName, DYHistName_, "Signal_"+var, steering, useAxis, divBinWidth); ;
    }

    for(int i=1;i<=hDY->GetNbinsX();i++)
    {
        ndf += 1.;
        if(hData->GetBinError(i) == 0)
        {
          std::cout << "unfolded " << i << " bin: " << hData->GetBinContent(i) << std::endl;
          std::cout << "error is zero in the " << i << " bin..." << std::endl;
          std::cout << "so skip this bin" << std::endl;
          continue;
        }
        // TODO add option to use input covariance matrix
        double pull=(hData->GetBinContent(i)-hDY->GetBinContent(i))/hData->GetBinError(i);
        //cout << "data: " << hData->GetBinContent(i) << " mc: " << hDY->GetBinContent(i) << " data error: " << hData->GetBinError(i) << endl;
        chi2+= pull*pull;
    }
    //cout << "chi^{2}, " << chi2 << endl;
    return chi2;
}

double ISRUnfold::getUnfoldedChi2(TString steering, bool useAxis, bool divBinWidth)
{
    divBinWidth = false;
    double chi2 = 0.;
    double ndf  = 0.;

    TH1* hData; // Data - Bkg
    TH1* hDY; // DY MC
    TH2* hRho;

    //TH1 *g_fcnHist=0;
    TMatrixD *g_fcnMatrix=0;

    if(var.Contains("Mass"))
    {
        hData = nomMassUnfold->GetOutput("hData_"+var, 0, 0, steering, useAxis);
        hDY = nomMassUnfold->GetBias("hData_"+var, 0, 0, steering, useAxis);
        hRho = nomMassUnfold->GetRhoIJtotal("histRho_chi_"+var, 0,0, steering, useAxis);
    }
    else
    {
        hData = nominalTUnfold->GetOutput("hData_"+var, 0, 0, steering, useAxis);
        hDY = nominalTUnfold->GetBias("hData_"+var, 0, 0, steering, useAxis);
        hRho = nominalTUnfold->GetRhoIJtotal("histRho_chi_"+var, 0,0, steering, useAxis);
    }

    // FIXME check if "n" or "n+1"
    int n = hData->GetNbinsX();
    TMatrixDSym v(n);
    for(int i=0;i<n;i++)
    {
       for(int j=0;j<n;j++)
        {
            v(i,j)=hRho->GetBinContent(i+1,j+1)*(hData->GetBinError(i+1)*hData->GetBinError(j+1));
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

    for(int i=0;i<n;i++)
    {
        if(TMath::Abs(test(i,i)-1.0)>1.E-7)
        {
            error++;
        }
        for(int j=0;j<n;j++)
        {
            if(i==j) continue;
            if(TMath::Abs(test(i,j)>1.E-7)) error++;
        }
    }


    // Calculate chi2
    //for(int i=0;i<hData->GetNbinsX();i++)
    //{
    //
    //    double di_=hData->GetBinContent(i+1)-hDY->GetBinContent(i+1);
    //    if(g_fcnMatrix)
    //    {
    //        for(int j=0;j<hData->GetNbinsX();j++)
    //        {
    //            double dj=hData->GetBinContent(j+1)-hDY->GetBinContent(j+1);
    //            chi2+=di_*dj*(*g_fcnMatrix)(i,j);
    //        }
    //    }
    //    else
    //    {
    //        double pull=di_/hData->GetBinError(i+1);
    //        chi2+=pull*pull;
    //    }
    //    ndf+=1.0;
    //}

    for(int i=1;i<=hDY->GetNbinsX();i++)
    {
        ndf += 1.;
        if(hData->GetBinError(i) == 0)
        {
          std::cout << "unfolded " << i << " bin: " << hData->GetBinContent(i) << std::endl;
          std::cout << "error is zero in the " << i << " bin..." << std::endl;
          std::cout << "so skip this bin" << std::endl;
          continue;
        }
        // TODO add option to use input covariance matrix
        double pull=(hData->GetBinContent(i)-hDY->GetBinContent(i))/hData->GetBinError(i);
        //cout << "data: " << hData->GetBinContent(i) << " mc: " << hDY->GetBinContent(i) << " data error: " << hData->GetBinError(i) << endl;
        chi2+= pull*pull;
    }
    //cout << "chi^{2}, " << chi2 << endl;
    return chi2;

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
void ISRUnfold::setNominalRM(TString filepath, TString dirName, TString histName, TString binDef)
{
    //cout << "ISRUnfold::setNominalRM set response matrix..." << endl;
    TH1::AddDirectory(kFALSE);
    TFile* filein = new TFile(filepath);

    TString fullDirPath = dirName + "/" + var + "_ResMatrix_" + histName + binDef + "/";

    TString Rec_binName = "Rec_"+var;
    TString Gen_binName = "Gen_"+var;
    Rec_binName = fullDirPath + Rec_binName;
    Gen_binName = fullDirPath + Gen_binName;

    // Set bin definition
    binning_Rec = (TUnfoldBinning*)filein->Get(Rec_binName);
    binning_Gen = (TUnfoldBinning*)filein->Get(Gen_binName);

    // Set response matrix
    // First, get the response matrix
    TH2* hmcGenRec;
    hmcGenRec = (TH2*)filein->Get(fullDirPath + "hmc" + var + "GenRec");

    nominalTUnfold = new TUnfoldDensity(hmcGenRec, TUnfold::kHistMapOutputHoriz, regMode, TUnfold::kEConstraintArea, TUnfoldDensity::kDensityModeNone, binning_Gen, binning_Rec);
    cout << "Check TUnfold version " << nominalTUnfold->GetTUnfoldVersion() << endl;
    hResponseM = (TH2*) hmcGenRec->Clone("hResponseM");

    // For statistical uncertainty
    if(makeStatUnfold)
    {
        // cout << "Create response matrix for statistical uncertainty..." << endl;
        for(int i = 0; i < statSize; i++)
        {
            statisticalTUnfold.push_back(new TUnfoldDensity(hmcGenRec, TUnfold::kHistMapOutputHoriz, regMode, TUnfold::kEConstraintArea, TUnfoldDensity::kDensityModeNone, binning_Gen, binning_Rec));
        }
    }

    filein->Close();
    delete filein;
}

void ISRUnfold::setFromPrevUnfResult(ISRUnfold* unfold, bool useAccept)
{
    //cout << "setFromPrevUnfResult(), useAccept? " << useAccept << endl;
    // Loop over sytematics considered in the previous unfold class
    // So first get sysMap map object
    std::map<TString, std::vector<TString>> sysMap_previous = unfold->getSystematicMap();
    std::map<TString, std::vector<TString>>::iterator it = sysMap_previous.begin();
    while(it != sysMap_previous.end())
    {
        //cout << "Systematic name: " << it->first << endl;
        if(this->sysMap.find(it->first) == this->sysMap.end())
        {
            // Not found in this ISRUnfold class, but exits in the previous one
            int size = sysMap_previous[it->first].size();
            for(int ith = 0; ith < size; ith++)
            {
                // Create TUnfoldDensity using the DEFAULT response matrix
                //cout << "Systematic variation, " << sysMap_previous[it->first][ith] << endl;
                if((it->first).Contains("Unfolding") && !(sysMap_previous[it->first][ith]).Contains("Nominal"))
                {
                    this->iterEMTUnfold   = new TUnfoldIterativeEM(hResponseM,TUnfoldDensity::kHistMapOutputHoriz,binning_Gen,binning_Rec);

                    if(!useAccept)
                    {
                        this->iterEMTUnfold->SetInput(unfold->iterEMTUnfold->GetOutput("hUnfolded" + var + "_"+ it->first + "_" + sysMap_previous[it->first][ith],0,0,"*[*]",false), nominal_bias);
                    }
                    else
                    {
                        //cout << "use acceptance corrected output!" << endl;
                        this->iterEMTUnfold->SetInput(unfold->hSysFullPhaseData[it->first][sysMap_previous[it->first][ith]], nominal_bias);
                    }
                }
                else
                {
                    this->systematicTUnfold[it->first][sysMap_previous[it->first][ith]]   = new TUnfoldDensity(hResponseM,TUnfold::kHistMapOutputHoriz,regMode, TUnfold::kEConstraintArea, TUnfoldDensity::kDensityModeNone, binning_Gen,binning_Rec);

                    if(!useAccept)
                    {
                         this->systematicTUnfold[it->first][sysMap_previous[it->first][ith]]->SetInput(unfold->systematicTUnfold[it->first][sysMap_previous[it->first][ith]]->GetOutput("hUnfolded" + var + "_"+ it->first + "_" + sysMap_previous[it->first][ith],0,0,"*[*]",false), nominal_bias);
                    }
                    else
                    {
                        //cout << "use acceptance corrected output!" << endl;
                        this->systematicTUnfold[it->first][sysMap_previous[it->first][ith]]->SetInput(unfold->hSysFullPhaseData[it->first][sysMap_previous[it->first][ith]], nominal_bias);
                    }
                }

                this->sysMap[it->first].push_back(sysMap_previous[it->first][ith]);
            }
        }
        else
        {
            // Found
            // Systematic both considered in this and previous unfolding
            // Loop over systematic varations
            int size = sysMap_previous[it->first].size();
            for(int ith = 0; ith < size; ith++)
            {
                if(!useAccept)
                {
                    this->systematicTUnfold[it->first][sysMap_previous[it->first][ith]]->SetInput(unfold->systematicTUnfold[it->first][sysMap_previous[it->first][ith]]->GetOutput("hUnfolded"+var+"_"+ it->first + "_" + sysMap_previous[it->first][ith],0,0,"*[*]",false), nominal_bias);
                }
                else
                {

                    if((it->first).Contains("Unfolding") && !(sysMap_previous[it->first][ith]).Contains("Nominal"))
                    {
                        this->iterEMTUnfold->SetInput(unfold->hSysFullPhaseData[it->first][sysMap_previous[it->first][ith]], nominal_bias);
                    }
                    else
                    {
                        this->systematicTUnfold[it->first][sysMap_previous[it->first][ith]]->SetInput(unfold->hSysFullPhaseData[it->first][sysMap_previous[it->first][ith]], nominal_bias);
                    }
                }
            }
        }
        it++;
    }
    // Loop over variations of event selection efficiency correction
}

void ISRUnfold::setSystematicRM(TString filepath, TString dirName, TString histName, TString sysName, TString sysPostfix, TString histPostfix, TString binDef)
{
    TFile* filein = new TFile(filepath);
    TH2* hmcGenRec = NULL;

    TString fullHistName = "hmc" + var + "GenRec";
    if(histPostfix != "")
        fullHistName = fullHistName + "_" + sysPostfix;

    hmcGenRec = (TH2*)filein->Get(dirName + "/" + var + "_ResMatrix_" + histName + binDef + "/" + fullHistName);

    
    if(sysName.Contains("Unfolding") && !sysPostfix.Contains("Nominal"))
    {
        iterEMTUnfold = new TUnfoldIterativeEM(hmcGenRec,TUnfoldDensity::kHistMapOutputHoriz,binning_Gen,binning_Rec);
    }
    else
    {
        systematicTUnfold[sysName][sysPostfix] = new TUnfoldDensity(hmcGenRec, TUnfold::kHistMapOutputHoriz, regMode, TUnfold::kEConstraintArea, TUnfoldDensity::kDensityModeNone, binning_Gen, binning_Rec);
    }

    filein->Close();
    delete filein;
}

// Set input histogram from unfolding output
void ISRUnfold::setUnfInput(ISRUnfold* unfold, bool isSys, TString sysName, TString sysPostfix, bool useAccept)
{
    TH1::AddDirectory(kFALSE);

    if(!useAccept)
    {
        if(!isSys)
        {
            nominalTUnfold->SetInput(unfold->getUnfoldedHists(var, "UnfoldOut_"+var, "*[*]", false), 1.);
            
        }
        else
        {
            
            systematicTUnfold[sysName][sysPostfix]->SetInput(unfold->getUnfoldedHists(var, "UnfoldOut_"+var+sysName+sysPostfix, "*[*]", false), 1.);
            
        }
    }
    else
    {
        //cout << "set from previous unfold class, isSys " << isSys << endl;
        if(!isSys)
        {
            nominalTUnfold->SetInput(unfold->hFullPhaseData, 1.);
        }
        else
        {
            if(sysName.Contains("Unfolding") && !sysPostfix.Contains("Nominal"))
            {
                iterEMTUnfold->SetInput(unfold->hFullPhaseData, 1.);
            }
            else
            {
                systematicTUnfold[sysName][sysPostfix]->SetInput(unfold->hFullPhaseData, 1.);
            }
        }
    }
}

// Set input histogram from root file
void ISRUnfold::setUnfInput(TString varPostfix, TString filepath, TString dirName, TString histName, bool isSys, TString sysName, TString sysPostfix, bool isFSR)
{
    TH1::AddDirectory(kFALSE);

    TFile* filein = new TFile(filepath);
    TH1* hRec = NULL;
    //cout << dirName+"/"+var+varPostfix+"/"+histName << endl;
    hRec = (TH1*)filein->Get(dirName+"/"+var+varPostfix+"/"+histName);

    // Use DY MC as unfolding input, i.e. simple closure test
    if(!isFSR)
    {
        if(histName.Contains("DYJetsTo"))
        {
            histName.ReplaceAll("DYJetsTo", "DYJets10to50To");
            hRec->Add((TH1*)filein->Get(dirName+"/"+var+varPostfix+"/"+histName));
        }
    }
    else
    {
        if(histName.Contains("DYJets"))
        {
            histName.ReplaceAll("DYJets", "DYJets10to50");
            hRec->Add((TH1*)filein->Get(dirName+"/"+var+varPostfix+"/"+histName));
        }
    }

    // Very preliminary test for input covariance using ID SF
    //TFile* fcov = new TFile("/home/jhkim/ISR_Run2/unfolding/TUnfoldISR2016/rootScripts/covariance.root");
    //TFile* fcov_pt = new TFile("/home/jhkim/ISR_Run2/unfolding/TUnfoldISR2016/rootScripts/covariance_pt.root");
    //TH2* hCov = (TH2*) fcov->Get("cov");
    //TH2* hCov_pt = (TH2*) fcov_pt->Get("cov");

    // Nominal
    if(!isSys)
    {
        nominalTUnfold->SetInput(hRec,   nominal_bias, 0);
    }
    else
    // Systematic histograms
    {
        if(sysName.Contains("Unfolding") && !sysPostfix.Contains("Nominal"))
        {
            iterEMTUnfold->SetInput(hRec, nominal_bias);
        }
        else
        {
            systematicTUnfold[sysName][sysPostfix]->SetInput(hRec, nominal_bias);
        }
    }

    filein->Close();
    delete filein;
}

void ISRUnfold::subBkgs(TString filepath, std::pair<TString, TString>& bkgInfo, bool isSys, TString binDef, TString dirName, TString sysName, TString sysPostfix, TString histPostfix)
{
    TFile* filein = new TFile(filepath);
    TH1* hRec = NULL;

    // Nominal histograms
    if(!isSys)
    {
        hRec = (TH1*)filein->Get(dirName + "/" + var + binDef+"/histo_" + bkgInfo.first);
        nominalTUnfold->  SubtractBackground(hRec, bkgInfo.first);
    }
    else
    // Systematic
    {
        TString fullHistName = bkgInfo.first + "_" + sysPostfix;
        if(histPostfix == "")
            fullHistName = bkgInfo.first;

        hRec = (TH1*)filein->Get(dirName + "/" + var + binDef+"/histo_" + fullHistName);

        //cout << "file path: " << filepath << endl;
        //cout << dirName + "/Pt"+binDef+"/histo_" + fullHistName << endl;

        if(sysName.Contains("Unfolding") && !sysPostfix.Contains("Nominal"))
        {
            iterEMTUnfold->SubtractBackground(hRec, bkgInfo.first);
        }
        else
        {
            // FIXME temporary method for background systematic
            if(sysName.Contains("Background"))
            {
                if(sysPostfix.Contains("Up"))
                {
                    systematicTUnfold[sysName][sysPostfix]->SubtractBackground(hRec, bkgInfo.first, 1.05);
                }
                if(sysPostfix.Contains("Down"))
                {
                    systematicTUnfold[sysName][sysPostfix]->SubtractBackground(hRec, bkgInfo.first, 0.95);
                }
            }
            else
            {
                systematicTUnfold[sysName][sysPostfix]->SubtractBackground(hRec, bkgInfo.first);
            }
        }
    }

    filein->Close();
    delete filein;
}

void ISRUnfold::setSystematics(TString sysName, TString sysHistName)
{
    sysMap[sysName].push_back(sysHistName);
}

void ISRUnfold::divideByBinWidth(TH1* hist, bool norm)
{
    for(int ibin = 1; ibin < hist->GetXaxis()->GetNbins()+1; ibin++)
    {
        double binWidth = hist->GetBinWidth(ibin);
        hist->SetBinContent(ibin, hist->GetBinContent(ibin)/ binWidth);
        hist->SetBinError(ibin, hist->GetBinError(ibin)/ binWidth);
    }
    if(norm)
    {
        hist->Scale(1./ hist->Integral());
    }
}

void ISRUnfold::doStatUnfold()
{
    //cout << "ISRUnfold::doStatUnfold() " << endl;
    for(int istat = 0; istat < statSize; istat++)
    {
        //cout << istat << " th stat.." << endl;
        TH1* tempInput;

        TString nth_;
        nth_.Form("%d", istat);
        tempInput = nominalTUnfold->GetInput("temp" + var + "Hist_" + nth_, 0, 0, 0, false);

        // randomize histogram bin content
        for(int ibin = 1; ibin<tempInput->GetNbinsX()+1;ibin++)
        {
            double err = tempInput->GetBinError(ibin);
            if(err > 0.0)
            {
                tempInput->SetBinContent(ibin, tempInput->GetBinContent(ibin) + gRandom->Gaus(0,err));
            }
        }

        statisticalTUnfold.at(istat)->SetInput(tempInput, nominal_bias);
        statisticalTUnfold.at(istat)->DoUnfold(0);

        //fillPtStatVariationHist(istat);
        //fillMassStatVariationHist(istat);

        delete tempInput;
        delete statisticalTUnfold.at(istat);
    }
    statisticalTUnfold.clear();
}

void ISRUnfold::doISRUnfold(bool doSys)
{

    TString yearStr;
    yearStr.Form("%d", (int)year);

    TFile* f;
    TDirectory* topDir;
    TDirectory* varDir;

    TString fullPath=output_baseDir+unfold_name+"_"+channel_name+"_"+yearStr+"_"+var+".root";
    cout << fullPath << endl;
    if(!gSystem->AccessPathName(fullPath, kFileExists))
    {
        f=TFile::Open(fullPath, "UPDATE"); 

        topDir=f->GetDirectory("unfolded");
        varDir=f->GetDirectory("unfolded/"+var);
    }
    else
    {
        f=new TFile(fullPath, "CREATE");

        // Create directory
        topDir=f->mkdir("unfolded");
        varDir=topDir->mkdir(var); 

        varDir->cd();
        binning_Gen->Write();
    }

    // cout << "ISRUnfold::doISRUnfold!!" << endl;
    // Nominal
    if(!doSys)
    {
        //cout << "Unfold without systematic" << endl;
        // No regularisation
        if(regMode == TUnfold::kRegModeNone)
        {
            // Nominal unfolding
            nominalTUnfold->DoUnfold(0);
        }

        /*
        else
        {
            nomMassUnfold->DoUnfold(0);

            if(var=="Pt") 
            {
            int istart = binning_Gen->GetGlobalBinNumber(0., 200.);
            int iend = binning_Gen->GetGlobalBinNumber(99., 200.);
            nominalTUnfold->RegularizeBins(istart, 1, iend-istart+1, regMode);
            
            double tauMin=1.e-4;
            double tauMax=1.e-1;
            nominalTUnfold->ScanLcurve(100, tauMin, tauMax, 0);
            
            TH2 *histL= nominalTUnfold->GetL("L");
            if(histL)
            {
                for(Int_t j=1;j<=histL->GetNbinsY();j++)
                {
                    cout<<"L["<<nominalTUnfold->GetLBinning()->GetBinName(j)<<"]";
                    for(Int_t i=1;i<=histL->GetNbinsX();i++) {
                        Double_t c=histL->GetBinContent(i,j);
                        if(c!=0.0) cout<<" ["<<i<<"]="<<c;
                    }
                    cout<<"\n";
                }
            }
            }
        }
        */

        varDir->cd(); 
        nominalTUnfold->GetOutput("histo_Data",0,0, "*[*]", false)->Write(); 
        nominalTUnfold->GetBias("histo_DY", 0, 0, "*[*]", false)->Write(); 
       
    }
    // For systematic
    else
    {
        std::map<TString, std::vector<TString>>::iterator it = sysMap.begin();
        while(it != sysMap.end())
        {
            int size = (it->second).size();
            for(int i = 0; i < size; i++)
            {
                if( (it->first).Contains("Unfolding") && !((it->second).at(i)).Contains("Nominal"))
                {
                    iBest=iterEMTUnfold->ScanSURE(NITER_Iterative, &graph_SURE_IterativeSURE, &graph_DFdeviance_IterativeSURE);
                    //cout << "iBest pt, Mass: " << iBest_pt << " " << iBest_mass << endl;

                    varDir->cd();
                    iterEMTUnfold->GetOutput("histo_Data_"+(it->second).at(i),0,0, "*[*]", false)->Write(); 
                    nominalTUnfold->GetBias("histo_DY_"+(it->second).at(i), 0, 0, "*[*]", false)->Write(); 
                }
                else
                {

                    if(regMode == TUnfold::kRegModeNone)
                    {
                        systematicTUnfold[it->first][(it->second).at(i)]->DoUnfold(0);
                    }
                    else
                    {
                        double tauMin=1.e-4;
                        double tauMax=1.e-1;
                        systematicTUnfold[it->first][(it->second).at(i)]->ScanLcurve(100, tauMin, tauMax, 0);
                    }

                    varDir->cd();
                    systematicTUnfold[it->first][(it->second).at(i)]->GetOutput("histo_Data_"+(it->second).at(i),0,0, "*[*]", false)->Write(); 
                    systematicTUnfold[it->first][(it->second).at(i)]->GetBias("histo_DY_"+(it->second).at(i), 0, 0, "*[*]", false)->Write(); 
                }
            }
            it++;
        }
    }// Unfold for systematic
    topDir->Write();
    f->Close();
}

void ISRUnfold::doAcceptCorr(TString filePath, TString binDef, bool doSys, TString outName, bool isAccept)
{
    TString yearStr;
    yearStr.Form("%d", (int)year);

    TFile* f; 
    TDirectory* topDir;
    TDirectory* varDir;

    TString fullPath=output_baseDir+unfold_name+"_Acceptance_"+channel_name+"_"+yearStr+"_"+var+".root";
    if(!gSystem->AccessPathName(fullPath, kFileExists))
    {

        f=TFile::Open(fullPath, "UPDATE"); 

        topDir=f->GetDirectory("acceptance");
        varDir=f->GetDirectory("acceptance/"+var);
    }
    else
    {
        f=new TFile(fullPath, "CREATE");

        // Create directory
        topDir=f->mkdir("acceptance");
        varDir=topDir->mkdir(var);

        varDir->cd();
        binning_Gen->Write();
    }

    TFile* filein = new TFile(filePath);

    TString accepCorrOrEffCorr;
    if(isAccept)
        accepCorrOrEffCorr = "Acceptance";
    else
        accepCorrOrEffCorr = "Efficiency";

    TH1* hFiducialPhaseMC = NULL;

    // Pt
    hFullPhaseMC = (TH1*) filein->Get("Acceptance/"+var+"Gen" + binDef + "/histo_DYJets");
    if(year==2016)
        hFullPhaseMC->Add((TH1*) filein->Get("Acceptance/"+var+"Gen" + binDef + "/histo_DYJets10to50"));
    else
        hFullPhaseMC->Add((TH1*) filein->Get("Acceptance/"+var+"Gen" + binDef + "/histo_DYJets10to50_MG"));

    hFiducialPhaseMC = nominalTUnfold->GetBias("hFiducial"+var, 0, 0, "*[*]", false);
    hAcceptance = (TH1*) hFullPhaseMC->Clone("hAcceptance"+var);

    hAcceptanceFraction = (TH1*) hFiducialPhaseMC->Clone("hAcceptanceFraction"+var);
    hAcceptanceFraction->Divide(hFullPhaseMC);

    hAcceptance->Divide(hFiducialPhaseMC);
    hFullPhaseData = nominalTUnfold->GetOutput("histo_Data",0,0, "*[*]", false);
    hFullPhaseData->Multiply(hAcceptance);

    varDir->cd();
    hFullPhaseData->Write();
    hFullPhaseMC->SetName("histo_DY");
    hFullPhaseMC->Write();
    hAcceptance->Write();

    if(doSys)
    {
        std::map<TString, vector<TString>> sysMapForAcceptance;

        std::map<TString, std::vector<TString>>::iterator it = sysMap.begin();
        while(it != sysMap.end())
        {
            int size = (it->second).size();
            for(int i = 0; i < size; i++)
            {
                TH1* hFullPhaseMassMC_raw_sys = NULL;
                TH1* hFullPhaseMC_raw_sys = NULL;

                TH1* hFiducialPhaseMassMC_sys = NULL;
                TH1* hFiducialPhaseMC_sys = NULL;
                
                if((it->first).Contains("Unfolding") && !((it->second).at(i)).Contains("Nominal"))
                {
                    hSysFullPhaseData[it->first][(it->second).at(i)]   = iterEMTUnfold->GetOutput("histo_Data_"+(it->second).at(i),0,0, "*[*]", false);
                    hFiducialPhaseMC_sys=hFiducialPhaseMC;
                    hFiducialPhaseMC_sys->SetName("histo_DY_"+(it->second).at(i));
                }
                else
                {
                    hSysFullPhaseData[it->first][(it->second).at(i)]   = systematicTUnfold[it->first][(it->second).at(i)]->GetOutput("histo_Data_"+(it->second).at(i),0,0, "*[*]", false);
                    hFiducialPhaseMC_sys = systematicTUnfold[it->first][(it->second).at(i)]->GetBias("hFiducial"+var+"_sys"+(it->second).at(i), 0, 0, "*[*]", false);
                }

                // For PDF, AlphaS, Scale etc, denominator changed
                if( (((it->first).Contains("Scale") && !(it->first).Contains("Lep")) || (it->first).Contains("PDF") || (it->first).Contains("AlphaS")) && !(it->first).Contains("_") )
                {
                    hFullPhaseMC_raw_sys = (TH1*) filein->Get("Acceptance/"+var+ "Gen" + binDef + "/histo_DYJets_"+(it->second).at(i));
                    if(year==2016)
                        hFullPhaseMC_raw_sys->Add((TH1*) filein->Get("Acceptance/"+var+"Gen" + binDef + "/histo_DYJets10to50_"+(it->second).at(i)));
                    else
                        hFullPhaseMC_raw_sys->Add((TH1*) filein->Get("Acceptance/"+var+"Gen" + binDef + "/histo_DYJets10to50_MG_"+(it->second).at(i)));
                }
                else
                {
                    hFullPhaseMC_raw_sys=hFullPhaseMC;
                }

                // For pt
                TH1* hAcceptance_sys = (TH1*) hFullPhaseMC_raw_sys->Clone("hAcceptance_sys");
                hAcceptance_sys->Divide(hFiducialPhaseMC_sys);

                TH1* hAcceptanceFraction_sys = (TH1*) hFiducialPhaseMC_sys->Clone("hAcceptanceFraction_sys");
                hAcceptanceFraction_sys->Divide(hFullPhaseMC_raw_sys);

                //hSysFullPhaseData[accepCorrOrEffCorr + "_" + it->first][(it->second).at(i)] = nominalTUnfold->GetOutput("hAcceptPtData" +it->first+(it->second).at(i),0,0, "*[*]", false);
                hSysFullPhaseData[it->first][(it->second).at(i)]->Multiply(hAcceptance_sys);
                hSysFullPhaseMC[it->first][(it->second).at(i)] = hFullPhaseMC_raw_sys;
                sysMapForAcceptance[it->first].push_back((it->second).at(i)); // Update sysMapForAcceptance 

                //hSysAcceptance[it->first][(it->second).at(i)] = (TH1*) hAcceptance_sys->Clone("Pt_" + it->first + "_" + (it->second).at(i));
                hSysAcceptanceFraction[it->first][(it->second).at(i)] = (TH1*) hAcceptanceFraction_sys->Clone("Fraction"+var+"_" + it->first + "_" + (it->second).at(i));
                delete hAcceptance_sys;
                delete hAcceptanceFraction_sys;

                //hSysFullPhaseMassData[it->first][(it->second).at(i)]->Multiply(hAcceptanceMass);
                //hSysFullPhaseData[it->first][(it->second).at(i)]->Multiply(hAcceptance);
                hSysFullPhaseMC[it->first][(it->second).at(i)]   = hFullPhaseMC_raw_sys;

                varDir->cd();

                hSysFullPhaseData[it->first][(it->second).at(i)]->Write();
                hFullPhaseMC_raw_sys->SetName("histo_DY_"+(it->second).at(i));
                hFullPhaseMC_raw_sys->Write();
                hSysAcceptance[it->first][(it->second).at(i)] = (TH1*) hFullPhaseMC->Clone(var + "_" + it->first + "_" + (it->second).at(i));
            }
            it++;
        }

        // Update sys map
        // Sys. unc. of acceptance correction
        // PDF, alphaS, Scale
        //it = sysMapForAcceptance.begin();
        //while(it != sysMapForAcceptance.end())
        //{
        //    int size = (it->second).size();
        //    for(int i = 0; i < size; i++)
        //    {
        //       sysMap[it->first].push_back((it->second).at(i)); 
        //    }
        //    it++;
        //}
  
        // Stat. unc. of acceptance correction

        hSysFullPhaseData[accepCorrOrEffCorr + "_Stat"]["Up"] = nominalTUnfold->GetOutput("hFullPhaseData"+accepCorrOrEffCorr+"StatUp",0,0, "*[*]", false); 
        hSysFullPhaseData[accepCorrOrEffCorr + "_Stat"]["Down"] = nominalTUnfold->GetOutput("hFullPhaseData"+accepCorrOrEffCorr+"StatDown",0,0, "*[*]", false); 

        TH1* hAcceptance_statUp = (TH1*) hAcceptance->Clone("h"+accepCorrOrEffCorr+"StatUp"+var); 
        TH1* hAcceptance_statDown = (TH1*) hAcceptance->Clone("h"+accepCorrOrEffCorr+"StatDown"+var); 
        varyHistWithStatError(hAcceptance_statUp, 1);
        varyHistWithStatError(hAcceptance_statDown, -1);

        hSysFullPhaseData[accepCorrOrEffCorr + "_Stat"]["Up"]->Multiply(hAcceptance_statUp);
        hSysFullPhaseData[accepCorrOrEffCorr + "_Stat"]["Down"]->Multiply(hAcceptance_statDown);
        hSysFullPhaseMC[accepCorrOrEffCorr + "_Stat"]["Up"] = hFullPhaseMC;
        hSysFullPhaseMC[accepCorrOrEffCorr + "_Stat"]["Down"] = hFullPhaseMC;
        //setSystematics(accepCorrOrEffCorr, "StatUp");
        //setSystematics(accepCorrOrEffCorr, "StatDown");
    }

    topDir->Write();;
    f->Close();

    delete hFiducialPhaseMC;
}

void ISRUnfold::varyHistWithStatError(TH1* hist, int sys)
{
    for(int ibin = 1; ibin < hist->GetNbinsX()+1; ibin++)  
    {
        hist->SetBinContent(ibin, hist->GetBinContent(ibin) + double(sys) * hist->GetBinError(ibin));
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

TH1* ISRUnfold::getUnfoldedHists(TString outHistName, TString steering, bool useAxis, bool binWidth)
{
    TH1* outHist = NULL;
    if(var == "Mass")
    {
        outHist = nomMassUnfold->GetOutput(outHistName,0,0,steering,useAxis);
        if(binWidth)
        {
            divideByBinWidth(outHist, false);
        }
        return outHist;
    }
    else
    {
        outHist = nominalTUnfold->GetOutput(outHistName,0,0,steering,useAxis);
        if(binWidth)
        {
            divideByBinWidth(outHist, false);
        }
        else
            return outHist;
    }
    return outHist;
}

TH1* ISRUnfold::getRawHist(TString filePath, TString dirName, TString histName, TString outHistName, TString steering, bool useAxis, bool divBinWidth)
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

        hist = binning_Rec->ExtractHistogram(outHistName, raw_hist, 0, useAxis, steering);

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

    if(divBinWidth)
        divideByBinWidth(hist, false);
    return hist;
}

