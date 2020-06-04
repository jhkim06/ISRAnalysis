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

const TVectorD& ISRUnfold::checkMatrixCond(TString var)
{
    TH2D* hProb;
    if(var == "Pt") hProb = (TH2D*) nomPtUnfold->GetProbabilityMatrix("hProb_"+var);
    if(var == "Mass") hProb = (TH2D*) nomMassUnfold->GetProbabilityMatrix("hProb_"+var);

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
            cout << var << ", Cond(): " << singularValues[colLwb]/min;
            break;
        }
        else
        {
           i--;
        }
    }
    return svdProb->GetSig();
}

double ISRUnfold::getSmearedChi2(TString var, TString filePath, TString dirName, TString steering, bool useAxis, bool divBinWidth)
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
        hDY = getRawHist(var, filePath, dirName, DYHistName_, "Signal_"+var, steering, useAxis, divBinWidth); ;
    }
    else
    {
        hData = nomPtUnfold->GetInput("hData_"+var, 0, 0, steering, useAxis);
        hDY = getRawHist(var, filePath, dirName, DYHistName_, "Signal_"+var, steering, useAxis, divBinWidth); ;
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

double ISRUnfold::getUnfoldedChi2(TString var, TString steering, bool useAxis, bool divBinWidth)
{

    double chi2 = 0.;
    double ndf  = 0.;

    TH1* hData; // Data - Bkg
    TH1* hDY; // DY MC
    TH2* hRho;

    TH1 *g_fcnHist=0;
    TMatrixD *g_fcnMatrix=0;

    if(var.Contains("Mass"))
    {
        hData = nomMassUnfold->GetOutput("hData_"+var, 0, 0, steering, useAxis);
        hDY = nomMassUnfold->GetBias("hData_"+var, 0, 0, steering, useAxis);
        hRho = nomMassUnfold->GetRhoIJtotal("histRho_chi_"+var, 0,0, steering, useAxis);
    }
    else
    {
        hData = nomPtUnfold->GetOutput("hData_"+var, 0, 0, steering, useAxis);
        hDY = nomPtUnfold->GetBias("hData_"+var, 0, 0, steering, useAxis);
        hRho = nomPtUnfold->GetRhoIJtotal("histRho_chi_"+var, 0,0, steering, useAxis);
    }

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
    cout << "chi^{2}, " << chi2 << endl;
    return chi2;

}

void ISRUnfold::setNomResMatrix(TString var, TString filepath, TString dirName, TString histName, TString binDef)
{
    cout << "ISRUnfold::setNomResMatrix set response matrix..." << endl;
    TFile* filein = new TFile(filepath);

    TString Rec_binName = "Rec_"+var;
    TString Gen_binName = "Gen_"+var;
    Rec_binName = dirName + "/" + var + "_ResMatrix_" + histName + binDef + "/" + Rec_binName;
    Gen_binName = dirName + "/" + var + "_ResMatrix_" + histName + binDef + "/" + Gen_binName;

    cout << "Rec_binName: " << Rec_binName << endl;
    cout << "Gen_binName: " << Gen_binName << endl;
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

    // Set mass bin edges
    // setBassBinEdges
    setMassBindEdges();

    // Set response matrix
    TH2* hmcGenRec;
    hmcGenRec = (TH2*)filein->Get(dirName + "/" + var + "_ResMatrix_" + histName + binDef + "/hmc" + var + "GenRec");

    if( var == "Pt" )
    {
    	nomPtUnfold = new TUnfoldDensityV17(hmcGenRec,
    	                                    TUnfold::kHistMapOutputHoriz,
    	                                    regMode,
    	                                    TUnfold::kEConstraintArea,
    	                                    TUnfoldDensityV17::kDensityModeBinWidth,
    	                                    pt_binning_Gen,pt_binning_Rec);

        // For statistical uncertainty
        if(makeStatUnfold)
        {
            cout << "Create response matrix for statistical uncertainty..." << endl;
            for(int i = 0; i < statSize; i++)
            {
                statPtUnfold.push_back(new TUnfoldDensityV17(hmcGenRec,
    	                                                    TUnfold::kHistMapOutputHoriz,
    	                                                    regMode,
    	                                                    TUnfold::kEConstraintArea,
    	                                                    TUnfoldDensityV17::kDensityModeBinWidth,
    	                                                    pt_binning_Gen,pt_binning_Rec));
            }
        }
    }
    else
    {
        nomMassUnfold = new TUnfoldDensityV17(hmcGenRec,
                                              TUnfold::kHistMapOutputHoriz,
                                              regMode,
                                              TUnfold::kEConstraintArea,
                                              TUnfoldDensityV17::kDensityModeBinWidth,
                                              mass_binning_Gen,mass_binning_Rec);

        // For statistical uncertainty
        if(makeStatUnfold)
        {
            for(int i = 0; i < statSize; i++)
            {
                statMassUnfold.push_back(new TUnfoldDensityV17(hmcGenRec,
                                                               TUnfold::kHistMapOutputHoriz,
                                                               regMode,
                                                               TUnfold::kEConstraintArea,
                                                               TUnfoldDensityV17::kDensityModeBinWidth,
                                                               mass_binning_Gen,mass_binning_Rec));
            }
        }
    }

    filein->Close();
    delete filein;
}

void ISRUnfold::setMassBindEdges()
{
    cout << "ISRUnfold::setMassBindEdges " << endl;
    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    const Double_t* massBins = temp_tvecd->GetMatrixArray();
    int nMassBinEdges = temp_tvecd->GetNrows();

    if(massBinEdges.size() == 0)
    {
        for(int i = 0 ; i < nMassBinEdges; i++)
        {
            massBinEdges.push_back(massBins[i]);
            cout << i << " th mass bin edge: " << massBins[i] << endl;
        }
    }
    else
    {
        cout << "ISRUnfold::setMassBindEdges massBinEdges already set." << endl;
    }

}

void ISRUnfold::setFromPrevUnfResult(ISRUnfold* unfold)
{
    cout << "setFromPrevUnfResult()" << endl;
    // Loop over sytematics considered in the previous unfold class
    // So first get sysMap map object
    std::map<TString, std::vector<TString>> sysMap_ = unfold->getSystematicMap();

    std::map<TString, std::vector<TString>>::iterator it = sysMap_.begin();
    while(it != sysMap_.end())
    {
        cout << "Systematic name: " << it->first << endl;

        if(sysMap.find(it->first) == sysMap.end())
        {
          // not found
        } else
        {
            // found
            sysPtUnfold[it->first][sysMap_[it->first][0]]->SetInput(unfold->sysPtUnfold[it->first][sysMap_[it->first][0]]->GetOutput("hUnfoldedPt_"+ it->first + "_" + sysMap_[it->first][0],0,0,"*[*]",false), 1.);
            //sysMassUnfold[it->first][sysPostfix]->SetInput(unfold->getDetUnfoldedHists("Mass", "UnfoldOut_Mass_"+it->first, "*[*]", false), 1.);
            //sysPtUnfold[sysName][sysMap[sysName][0]]->GetOutput("hUnfoldedPt_up",0,0,steering,useAxis);
        }
        it++;
    }

    // Check if systematic also exist this unfold class
    // If so, just set input
    // If not, create TUnfoldDensity and set input
    // Finally, update this sysMap!

}

void ISRUnfold::setSysTUnfoldDensity(TString var, TString filepath, TString dirName, TString histName, TString sysName, TString sysPostfix, TString binDef)
{
    TFile* filein = new TFile(filepath);
    TH2* hmcGenRec = NULL;
    if(sysName != "Fake")
        hmcGenRec = (TH2*)filein->Get(dirName + "/" + var + "_ResMatrix_" + histName + binDef +"/hmc" + var + "GenRec_" + sysPostfix);
    else
        hmcGenRec = (TH2*)filein->Get(dirName + "/" + var + "_ResMatrix_" + histName + binDef +"/hmc" + var + "GenRec");

    if( var == "Pt" )
    {
        //cout << "sys: " << sysName << " postfix: " << sysPostfix << endl;
        sysPtUnfold[sysName][sysPostfix] = new TUnfoldDensityV17(hmcGenRec,
                                                                 TUnfold::kHistMapOutputHoriz,
                                                                 regMode,
                                                                 TUnfold::kEConstraintArea,
                                                                 TUnfoldDensityV17::kDensityModeBinWidth,
                                                                 pt_binning_Gen,pt_binning_Rec);
    }
    else if( var == "Mass" )
    {
        sysMassUnfold[sysName][sysPostfix] = new TUnfoldDensityV17(hmcGenRec,
                                                                   TUnfold::kHistMapOutputHoriz,
                                                                   regMode,
                                                                   TUnfold::kEConstraintArea,
                                                                   TUnfoldDensityV17::kDensityModeBinWidth,
                                                                   mass_binning_Gen,mass_binning_Rec);
    }
    else
    {
        cout << "ISRUnfold::setSysTUnfoldDensity, only Pt and Mass available for var" << endl;
        exit (EXIT_FAILURE);
    }

    filein->Close();
    delete filein;
}

// Set input histogram from unfolding output
void ISRUnfold::setUnfInput(ISRUnfold* unfold, TString var, bool isSys, TString sysName, int nth)
{
    TH1::AddDirectory(kFALSE);

    if(!isSys)
    {
        if(var=="Pt")
        {
            nomPtUnfold->SetInput(unfold->getDetUnfoldedHists("Pt", "UnfoldOut_Pt", "*[*]", false), 1.);
        }
        else
        {
            nomMassUnfold->SetInput(unfold->getDetUnfoldedHists("Mass", "UnfoldOut_Mass", "*[*]", false), 1.);
        }
    }
    else
    {
        cout << "ISRUnfold::setUnfInput not ready for systematic..." << endl;
        exit(EXIT_FAILURE);
    }
}

// Set input histogram from root file
void ISRUnfold::setUnfInput(TString var, TString varPostfix, TString filepath, TString dirName, TString histName, bool isSys, TString sysName, TString sysPostfix)
{
    TH1::AddDirectory(kFALSE);

    TFile* filein = new TFile(filepath);
    TH1* hRec = NULL;
    hRec = (TH1*)filein->Get(dirName+"/"+var+varPostfix+"/"+histName);
    // Use DY MC as unfolding input, i.e. simple closure test
    if(histName.Contains("DYJetsTo"))
    {
        histName.ReplaceAll("DYJetsTo", "DYJets10to50To");
        hRec->Add((TH1*)filein->Get(dirName+"/"+var+"/"+histName));
    }

    // Nominal histograms
    // First input covariance test (ID)
    TFile* fcov = new TFile("/home/jhkim/ISR_Run2/unfolding/TUnfoldISR2016/rootScripts/covariance.root");
    TH2* hCov = (TH2*) fcov->Get("cov");
    if(!isSys)
    {
        if(var == "Pt")
        {
            nomPtUnfold->SetInput(hRec,   nominal_bias);
        }
        else if(var == "Mass")
        {
            nomMassUnfold->SetInput(hRec, nominal_bias, 0, hCov);
        }
        else{
            cout << "ISRUnfold::setUnfInput, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
        }
    }
    else
    // Systematic histograms
    {
        if(var == "Pt")
        {
            sysPtUnfold[sysName][sysPostfix]->SetInput(hRec, nominal_bias);
        }
        else if(var == "Mass")
        {
            //cout << "sysName: " << sysName << " postfix: " << sysPostfix << endl;
            sysMassUnfold[sysName][sysPostfix]->SetInput(hRec, nominal_bias);
        }
        else
        {
            cout << "ISRUnfold::setUnfInput, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
        }
    }
    filein->Close();
    delete filein;
}

void ISRUnfold::subBkgs(TString filepath, std::pair<TString, TString>& bkgInfo, bool isSys, TString binDef, TString dirName, TString sysName, TString sysPostfix)
{
    TFile* filein = new TFile(filepath);
    TH1* hPtRec = NULL;
    TH1* hMassRec = NULL;

    // Nominal histograms
    if(!isSys)
    {

        bkgNames.push_back(bkgInfo.first);
        bkgTypes.push_back(bkgInfo.second);

        hPtRec = (TH1*)filein->Get(dirName + "/Pt"+binDef+"/histo_" + bkgInfo.first);
        nomPtUnfold->  SubtractBackground(hPtRec, bkgInfo.first);

        hMassRec = (TH1*)filein->Get(dirName + "/Mass"+binDef+"/histo_" + bkgInfo.first);
        nomMassUnfold->SubtractBackground(hMassRec, bkgInfo.first);
    }
    else
    // Systematic
    {
        TString histPostfix = bkgInfo.first + "_" + sysPostfix;
        if(bkgInfo.second != "Fake" && sysName == "Fake") histPostfix = bkgInfo.first;
        if(bkgInfo.second == "Fake" && sysName != "Fake") histPostfix = bkgInfo.first;

        hPtRec = (TH1*)filein->Get(dirName + "/Pt"+binDef+"/histo_" + histPostfix);
        //cout << dirName + "/Pt"+binDef+"/histo_" + histPostfix << endl;
        sysPtUnfold[sysName][sysPostfix]->SubtractBackground(hPtRec, bkgInfo.first);

        hMassRec = (TH1*)filein->Get(dirName + "/Mass"+binDef+"/histo_" + histPostfix);
        sysMassUnfold[sysName][sysPostfix]->SubtractBackground(hMassRec, bkgInfo.first);
    }

    filein->Close();
    delete filein;
}

void ISRUnfold::setSystematics(TString sysName, TString sysHistName)
{
    sysMap[sysName].push_back(sysHistName);
}

// Draw detector distributions using input root file
TCanvas* ISRUnfold::drawFoldedHists(TString var, TString filePath, TString dirName, TString steering, bool useAxis, TString sysName, TString outName, int nthMassBin, bool divBinWidth)
{
    // If steering == "", then usual TH1 histogram
    // If seering != "", TH1 from TUnfold

    double meanDipt = 0.;
    double meanDipt_bkgsub = 0.;

    setTDRStyle();
    writeExtraText = true;
    extraText  = "Work in progress";
    gStyle->SetLineWidth(3);
    gStyle->SetFrameLineWidth(3);
    gROOT->ForceStyle();

    TH1::AddDirectory(kFALSE);
    cout << "ISRUnfold::drawFoldedHists, Draw plot!" << endl;

    TFile* filein = new TFile(filePath);

    // For nominal histogram
    TH1* hData = NULL;
    TH1* hDY = NULL;
    TH1* hMCtotal = NULL;
    TH1* hRatio = NULL;

    // For systematic
    // TODO consider PDF uncertainty etc.
    TH1* hDY_up = NULL;
    TH1* hMCtotal_up = NULL;
    TH1* hRatio_up = NULL;

    TH1* hDY_down = NULL;
    TH1* hMCtotal_down = NULL;
    TH1* hRatio_down = NULL;

    // FIXME consider electron channel name
    TString dataHistName_ = "histo_DoubleMuon";
    TString DYHistName_ = "histo_DYJetsToMuMu";
    if(channel_name == "electron")
    {
        dataHistName_ = "histo_DoubleEG";
        DYHistName_  = "histo_DYJetsToEE";
    }

    hData = getRawHist(var, filePath, dirName, dataHistName_, "Data", steering, useAxis, divBinWidth);
    hDY = getRawHist(var, filePath, dirName, DYHistName_, "Signal", steering, useAxis, divBinWidth);

    hMCtotal = (TH1*) hDY->Clone("hMCtotal");
    hRatio = (TH1*) hData->Clone("hRatio");

    if(sysName != "")
    {
        hDY_up = getRawHist(var, filePath, dirName, "histo_DYJetsToMuMu_" + sysMap[sysName][0], "Signal"+sysMap[sysName][0], steering, useAxis, divBinWidth);
        hMCtotal_up = (TH1*) hDY_up->Clone("hMCtotal_up");
        hRatio_up = (TH1*) hData->Clone("hRatio_up");

        hDY_down = getRawHist(var, filePath, dirName, "histo_DYJetsToMuMu_" + sysMap[sysName][1], "Signal"+sysMap[sysName][1], steering, useAxis, divBinWidth);
        hMCtotal_down = (TH1*) hDY_down->Clone("hMCtotal_down");
        hRatio_down = (TH1*) hData->Clone("hRatio_down");

    }

    // Create canvas
    TCanvas* c_out = new TCanvas("detector_level_"+var, "detector_level_"+var, 3200, 2800);
    c_out->Draw();
    c_out->cd();

    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetBottomMargin(0.01);
    pad1->SetTopMargin(0.1);
    pad1->SetTicks(1);
    pad1->SetLogy();
    if(var.Contains("Mass"))
        pad1->SetLogx();
    pad1->Draw();
    pad1->cd();

    hData->SetTitle("");
    hData->SetStats(false);
    hData->GetXaxis()->SetMoreLogLabels(true);
    hData->Draw("p9histe");
    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(3.2);
    hData->SetLineColor(kBlack);
    hData->GetYaxis()->SetTitle("Events/Bin");
    hData->SetMaximum(5e9);
    hData->SetMinimum(1e-1);

    hDY->SetFillColor(kYellow);

    TLegend* leg = new TLegend(0.7, 0.65, 0.95, 0.85,"","brNDC");
    //leg->SetNColumns(2);
    leg->SetTextFont(43);
    leg->SetTextSize(100);
    leg->SetFillStyle(0); // transparent
    leg->SetBorderSize(0);
    leg->AddEntry(hData, "Data", "pe");
    leg->AddEntry(hDY, "Drell-Yan", "F");

    THStack* hsMC = new THStack("hsMC", "hsMC");
    setTHStack(var, filePath, dirName, *hsMC, *hMCtotal, *leg, steering, useAxis, "", divBinWidth);
    hsMC->Add(hDY);

    THStack* hsMC_up;
    THStack* hsMC_down;
    if(sysName != "")
    {
        hsMC_up = new THStack("hsMC_up", "hsMC_up");
        hsMC_down = new THStack("hsMC_down", "hsMC_down");
        setTHStack(var, filePath, dirName, *hsMC_up, *hMCtotal_up, *leg, steering, useAxis, sysMap[sysName][0], divBinWidth);
        setTHStack(var, filePath, dirName, *hsMC_down, *hMCtotal_down, *leg, steering, useAxis, sysMap[sysName][1], divBinWidth);
    }

    hsMC->Draw("hist same");
    hData->Draw("p9histe same");
    pad1->RedrawAxis();

    // Get average transeverse momentum values
    meanDipt = hData->GetMean();
    TH1* hData_ = (TH1*) hData->Clone("DataBKGsubtracted");
    TH1* hMCtotal_ = (TH1*) hMCtotal->Clone("MCtotalDYsubtracted");
    hMCtotal_->Add(hDY, -1);
    hData_->Add(hMCtotal_, -1);
    meanDipt_bkgsub = hData_->GetMean();
    double binnedMean = getBinnedMean(hData); // To check how average value changes with binned histogram

    //TLatex meanDipt_;
    //TLatex meanDipt_bkgsub_;
    //meanDipt_.SetTextFont(43);
    //meanDipt_.SetTextSize(40);
    //meanDipt_bkgsub_.SetTextFont(43);
    //meanDipt_bkgsub_.SetTextSize(40);

    //if(var.Contains("Pt"))
    //{
    //    TString itos, itos_;
    //    itos.Form ("%.2f", meanDipt);
    //    itos_.Form ("%.2f", binnedMean);
    //    meanDipt_.DrawLatexNDC(0.2, 0.6, "avg. p_{T}^{Data}: "+itos+"("+itos_+") GeV");
    //    itos.Form ("%.2f", meanDipt_bkgsub);
    //    meanDipt_bkgsub_.DrawLatexNDC(0.2, 0.6-0.07, "avg. p_{T}^{Data-Bkg}: "+itos);
    //}

    //getSmearedChi2(TString var, TString filePath, TString dirName, TString steering, bool useAxis, bool divBinWidth)
    TLatex smearedChi2;
    smearedChi2.SetTextFont(43);
    smearedChi2.SetTextSize(100);
    TString chi2_;
    chi2_.Form("%.2f", getSmearedChi2(var, filePath, dirName, steering,  useAxis, false));
    smearedChi2.DrawLatexNDC(0.2, 0.6, "#chi^{2}: " + chi2_);

    leg->Draw();

    int iPeriod_ = 4;
    if(year == 2017)
        iPeriod_ = 5;
    if(year == 2018)
        iPeriod_ = 6;
    CMS_lumi(pad1, iPeriod_, 11);
    // writeCutInfo(pad, var, nthMassBin);
    writeCutInfo(pad1, var, nthMassBin);

    TLine massEdgeLine;
    if(var=="Mass")
    {
        for(unsigned int i = 0; i < massBinEdges.size(); i++)
        {
            if(i==0) continue;
            if(i==massBinEdges.size()-1) continue;
            massEdgeLine.SetLineStyle(2);
            massEdgeLine.DrawLine(massBinEdges[i], hData->GetMinimum(), massBinEdges[i], hData->GetMaximum());
        }
    }
    c_out->cd();

    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.35);
    if(var.Contains("Mass"))
        pad2->SetLogx();
    pad2->SetTicks(1);
    pad2->SetGridy(1);
    pad2->Draw();
    pad2->cd();

    hRatio->SetStats(false);
    hRatio->Divide(hMCtotal);
    hRatio->Draw("p9histe");
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(3.2);
    hRatio->SetLineColor(kBlack);
    hRatio->GetYaxis()->SetTitle("Data/MC");

    hRatio->SetMinimum(0.7);
    hRatio->SetMaximum(1.3);

    setXaxisTitle(hRatio, var, useAxis);

    // TODO Save systematic histograms
    TH1* sysBand_ratio = NULL;
    if(sysName != "")
    {
        hRatio_up->Divide(hMCtotal_up);
        hRatio_down->Divide(hMCtotal_down);
        sysBand_ratio = (TH1*)hRatio_up->Clone("sysBand_ratio");
        for(int ibin = 1; ibin < sysBand_ratio->GetNbinsX()+1; ibin++)
        {
            double delta = fabs(hRatio_up->GetBinContent(ibin) - hRatio_down->GetBinContent(ibin));
            sysBand_ratio->SetBinError(ibin, delta);
            sysBand_ratio->SetBinContent(ibin, 1.);
        }
        sysBand_ratio->SetFillColorAlpha(kRed,0.8);
        sysBand_ratio->SetFillStyle(3004);
        sysBand_ratio->SetMarkerSize(0.);
        sysBand_ratio->Draw("E2 same");

        sysRelHist_detector[sysName] = sysBand_ratio; 
    }
    hRatio->Draw("p9histe same");

    TLine* l_ = new TLine(hRatio->GetXaxis()->GetXmin(),1,hRatio->GetXaxis()->GetXmax(),1);
    l_->SetLineColor(kRed);
    l_->Draw("same");
    l_->SetLineStyle(3);

    // Save canvas
    c_out->cd();
    c_out->SaveAs(outName!=""?outName+".png":"detector_"+var+".png");

    delete filein;

    return c_out;
}

TCanvas* ISRUnfold::drawUnfoldedHists(TString var, TString steering, bool useAxis, TString sysName, TString outName, int nthMassBin, bool divBinWidth)
{
    // If steering == "", then usual TH1 histogram
    // If seering != "", TH1 from TUnfold

    setTDRStyle();
    writeExtraText = true;
    extraText  = "Work in progress";
    gStyle->SetLineWidth(3.);
    gStyle->SetFrameLineWidth(3.);
    gROOT->ForceStyle();

    TH1::AddDirectory(kFALSE);
    cout << "ISRUnfold::drawFoldedHists, Draw plot!" << endl;

    // For nominal histogram
    TH1* hData = NULL;
    TH1* hDY = NULL;
    TH1* hRatio = NULL;

    // For systematic
    // TODO consider PDF uncertainty etc.
    // Make a function returning systematic band for systematic source
    TH1* hDY_up = NULL;
    TH1* hRatio_up = NULL;

    TH1* hDY_down = NULL;
    TH1* hRatio_down = NULL;

    TH1* hData_up = NULL;
    TH1* hData_down = NULL;

    if(var.Contains("Pt"))
    {
        hData = nomPtUnfold->GetOutput("hUnfoldedPt",0,0,steering,useAxis);
        hDY = nomPtUnfold->GetBias("hDYMCPt",0,0,steering,useAxis);;
    }
    else
    {
        hData = nomMassUnfold->GetOutput("hUnfoldedMass",0,0,steering,useAxis);
        hDY = nomMassUnfold->GetBias("hDYMCMass",0,0,steering,useAxis);;
    }
    if(divBinWidth)
    {
        divideByBinWidth(hData, false);
        divideByBinWidth(hDY, false);
    }
    hRatio = (TH1*) hData->Clone("hRatio");

    if(sysName != "")
    {
        if(var.Contains("Pt"))
        {
            hDY_up = sysPtUnfold[sysName][sysMap[sysName][0]]->GetBias("hDYMCPt_up",0,0,steering,useAxis);;
            hDY_down = sysPtUnfold[sysName][sysMap[sysName][1]]->GetBias("hDYMCPt_down",0,0,steering,useAxis);;

            hData_up = sysPtUnfold[sysName][sysMap[sysName][0]]->GetOutput("hUnfoldedPt_up",0,0,steering,useAxis);
            hData_down = sysPtUnfold[sysName][sysMap[sysName][1]]->GetOutput("hUnfoldedPt_down",0,0,steering,useAxis);
        }
        else
        {
            hDY_up = sysMassUnfold[sysName][sysMap[sysName][0]]->GetBias("hDYMCMass_up",0,0,steering,useAxis);;
            hDY_down = sysMassUnfold[sysName][sysMap[sysName][1]]->GetBias("hDYMCMass_down",0,0,steering,useAxis);;

            hData_up = sysMassUnfold[sysName][sysMap[sysName][0]]->GetOutput("hUnfoldedMass_up",0,0,steering,useAxis);
            hData_down = sysMassUnfold[sysName][sysMap[sysName][1]]->GetOutput("hUnfoldedMass_down",0,0,steering,useAxis);
        }
        if(divBinWidth)
        {
            divideByBinWidth(hDY_up, false);
            divideByBinWidth(hDY_down, false);
            divideByBinWidth(hData_up, false);
            divideByBinWidth(hData_down, false);
        }
        hRatio_up = (TH1*) hData_up->Clone("hRatio_up");
        hRatio_down = (TH1*) hData_down->Clone("hRatio_down");
    }

    // Create canvas
    TCanvas* c_out = new TCanvas("unfolded_level_"+var, "unfoled_level_"+var, 3200, 2800);
    c_out->Draw();
    c_out->cd();

    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetBottomMargin(0.01);
    pad1->SetTopMargin(0.1);
    pad1->SetTicks(1);
    pad1->SetLogy();
    if(var.Contains("Mass"))
        pad1->SetLogx();
    pad1->Draw();
    pad1->cd();

    hData->SetTitle("");
    hData->SetStats(false);
    hData->GetXaxis()->SetMoreLogLabels(true);
    hData->Draw("p9histe");
    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(3.2);
    hData->SetLineColor(kBlack);
    hData->GetYaxis()->SetTitle("Events/Bin");
    hData->SetMaximum(5e9);
    hData->SetMinimum(2e-1);

    if(sysName != "")
    {
        hData_up->SetLineStyle(2);
        hData_down->SetLineStyle(2);
        hData_up->Draw("HIST SAME");
        hData_down->Draw("HIST SAME");
    }

    hDY->SetFillColor(kYellow);
    hDY->Draw("hist same");

    TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9,"","brNDC");
    //leg->SetNColumns(2);
    leg->SetTextFont(43);
    leg->SetTextSize(80);
    leg->SetFillStyle(0); // transparent
    leg->SetBorderSize(0);
    leg->AddEntry(hData, "Unfolded data", "pe");
    leg->AddEntry(hDY, "Drell-Yan (MG5_aMC@NLO)", "F");

    hData->Draw("p9histe same");
    pad1->RedrawAxis();

    leg->Draw();

    int iPeriod_ = 4;
    if(year == 2017)
        iPeriod_ = 5;
    if(year == 2018)
        iPeriod_ = 6;
    CMS_lumi(pad1, iPeriod_, 11);
    // writeCutInfo(pad, var, nthMassBin);
    writeCutInfo(pad1, var, nthMassBin);

    // getUnfoldedChi2(TString var, TString steering, bool useAxis, bool divBinWidth)
    TLatex unfoldedChi2;
    unfoldedChi2.SetTextFont(43);
    unfoldedChi2.SetTextSize(50);
    TString chi2_;
    chi2_.Form("%.2f", getUnfoldedChi2(var, steering, useAxis, divBinWidth));
    unfoldedChi2.DrawLatexNDC(0.2, 0.6, "#chi^{2}: " + chi2_);

    TLine massEdgeLine;
    if(var=="Mass")
    {
        for(unsigned int i = 0; i < massBinEdges.size(); i++)
        {
            if(i==0) continue;
            if(i==massBinEdges.size()-1) continue;
            massEdgeLine.SetLineStyle(2);
            massEdgeLine.DrawLine(massBinEdges[i], hData->GetMinimum(), massBinEdges[i], hData->GetMaximum());
        }
    }
    c_out->cd();

    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.35);
    if(var=="Mass")
    {
        pad2->SetLogx();
        hRatio->GetXaxis()->SetMoreLogLabels(true);
    }
    pad2->SetTicks(1);
    pad2->SetGridy(1);
    pad2->Draw();
    pad2->cd();

    hRatio->SetStats(false);
    hRatio->Divide(hDY);
    hRatio->Draw("p9histe");
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(3.2);
    hRatio->SetLineColor(kBlack);
    hRatio->GetYaxis()->SetTitle("Unfolded/MC");

    hRatio->SetMinimum(0.7);
    hRatio->SetMaximum(1.3);

    setXaxisTitle(hRatio, var, useAxis);

    // TODO Save systematic histograms
    TH1* sysBand_ratio = NULL;
    if(sysName != "")
    {
        hRatio_up->Divide(hDY);
        hRatio_down->Divide(hDY);
        sysBand_ratio = (TH1*)hRatio_up->Clone("sysBand_ratio");
        for(int ibin = 1; ibin < sysBand_ratio->GetNbinsX()+1; ibin++)
        {
            double delta = fabs(hRatio_up->GetBinContent(ibin) - hRatio_down->GetBinContent(ibin));
            sysBand_ratio->SetBinError(ibin, delta);
            sysBand_ratio->SetBinContent(ibin, hRatio->GetBinContent(ibin));
        }
        sysBand_ratio->SetFillColorAlpha(kBlack,0.8);
        sysBand_ratio->SetFillStyle(3004);
        sysBand_ratio->SetMarkerSize(0.);
        sysBand_ratio->Draw("E2 same");
    }
    hRatio->Draw("p9histe same");

    TLine* l_ = new TLine(hRatio->GetXaxis()->GetXmin(),1,hRatio->GetXaxis()->GetXmax(),1);
    l_->SetLineColor(kRed);
    l_->Draw("same");
    l_->SetLineStyle(3);

    // Save canvas
    c_out->cd();
    c_out->SaveAs(outName!=""?outName+".png":"unfolded_"+var+".png");

    return c_out;
}

TH1* ISRUnfold::getDetectorSystematicBand(TString var, TString steering, bool useAxis, bool isRatio, bool isMC)
{
    //
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

double ISRUnfold::getBinnedMean(TH1* hist)
{
    TH1* hBinned = (TH1*) hist->Clone("Binned");
    TH1* hDummy = (TH1*) hist->Clone("Dummy");
    hDummy->Reset("ICES");

    hBinned->Add(hDummy, -1); // After Add() Sumwx information removed, so can get binned mean.
    return hBinned->GetMean();
}

void ISRUnfold::setXaxisTitle(TH1* hist, TString var, bool useAxis, TString title)
{
    TString channel_name_;
    if(channel_name=="electron") channel_name_ = "ee";
    else channel_name_ = "#mu#mu";

    if(useAxis)
    {
        TString title_;
        if(var.Contains("Mass"))
        {
            title_ = "Mass^{" + channel_name_  + "} [GeV]";
            if(title != "") title_ = title;
            hist->GetXaxis()->SetTitle(title_);
        }
        if(var.Contains("Pt"))
        {
            title_ = "p_{T}^{" + channel_name_  + "} [GeV]";
            if(title != "") title_ = title;
            hist->GetXaxis()->SetTitle(title_);
        }
    }
    else
    {
        if(var.Contains("Mass"))
        {
            hist->GetXaxis()->SetTitle("Mass bin index");
        }
        if(var.Contains("Pt"))
        {
            hist->GetXaxis()->SetTitle("p_{T} bin index");
        }
    }
}

void ISRUnfold::writeCutInfo(TPad* pad, TString var, int nthMassBin)
{
    pad->cd();

    double x_ = 0.2;
    double y_ = 0.7;
    TString mass_cut_info;
    TString lepton_cut_info;
    TString lepton_type;
    if(channel_name=="electron"){
        lepton_type = "ee";
        lepton_cut_info = "p_{T}>25(20) GeV, |#eta|<2.5";
    }
    else
    {
        lepton_type = "#mu#mu";
        lepton_cut_info = "p_{T}>20(10) GeV, |#eta|<2.4";
    }

    TString low_bound_, upper_bound_;
    low_bound_.Form("%d", (int)massBinEdges[nthMassBin]);
    upper_bound_.Form("%d", (int)massBinEdges[nthMassBin+1]);
    mass_cut_info = low_bound_ + " < M(" + lepton_type + ") < " + upper_bound_ + " (GeV)";

    TLatex dimass_cut;
    TLatex lepton_cut;
    dimass_cut.SetTextFont(43);
    dimass_cut.SetTextSize(40 * 2.5);
    lepton_cut.SetTextFont(43);
    lepton_cut.SetTextSize(40 * 2.5);
    if(var.Contains("Pt"))
    {
        dimass_cut.DrawLatexNDC(x_, y_, mass_cut_info);
        lepton_cut.DrawLatexNDC(x_, y_-0.05, lepton_cut_info);
    }
    if(var.Contains("Mass"))
    {
        //dimass_cut.DrawLatexNDC(x_, y_, mass_cut_info);
        lepton_cut.DrawLatexNDC(x_, y_, lepton_cut_info);
    }
}

void ISRUnfold::setTHStack(TString var, TString filePath, TString dirName, THStack& hs, TH1& hMCtotal, TLegend& leg, TString steering, bool useAxis, TString sysName, bool divBinWidth)
{
    TH1::AddDirectory(kFALSE);
    int bkgSize = bkgNames.size();

    // Count total number of each background type N
    map<TString, int> bkgTypeN;
    //cout << "N bkg: " << bkgSize << endl;
    for(int i = 0; i < bkgSize; i++)
    {
        if(bkgTypes[i] == "DY") continue;
        map<TString, int>::iterator it = bkgTypeN.find(bkgTypes[i]);
        if(it != bkgTypeN.end())
        {
            bkgTypeN[bkgTypes[i]]++;
        }
        else
        {
            //cout << bkgTypes[i] << " first found" << endl;
            bkgTypeN[bkgTypes[i]] = 1;
        }
    }

    TH1* htemp = NULL;
    bool isFirstBkg = true;
    int nthBkg = 0;

    for(int i = 0; i < bkgSize; i++)
    {
        if(bkgTypes[i] == "DY") continue;
        TString histName_ = "histo_" + bkgNames[i];;

        if(isFirstBkg)
        {
            if(sysName == "")
            {
                htemp = getRawHist(var, filePath, dirName, histName_, "h"+bkgNames[i], steering, useAxis, divBinWidth);
            }
            else
            {
                TString histPostfix = "histo_"+bkgNames[i] + "_" + sysName;
                if(bkgTypes[i] != "Fake" && sysName == "Fake") histPostfix = "histo_"+bkgNames[i];
                if(bkgTypes[i] == "Fake" && sysName != "Fake") histPostfix = "histo_"+bkgNames[i];
                htemp = getRawHist(var, filePath, dirName, histPostfix, "h"+bkgNames[i], steering, useAxis, divBinWidth);
            }
            isFirstBkg = false;
            nthBkg++;
        }
        else
        {
            if(sysName == "")
            {
                htemp->Add(getRawHist(var, filePath, dirName, histName_, "h"+bkgNames[i], steering, useAxis, divBinWidth));
            }
            else
            {
                TString histPostfix = "histo_"+bkgNames[i] + "_" + sysName;
                if(bkgTypes[i] != "Fake" && sysName == "Fake") histPostfix = "histo_"+bkgNames[i];
                if(bkgTypes[i] == "Fake" && sysName != "Fake") histPostfix = "histo_"+bkgNames[i];
                htemp->Add(getRawHist(var, filePath, dirName, histPostfix, "h"+bkgNames[i], steering, useAxis, divBinWidth));
            }
            nthBkg++;
        }

        // This type of backgrounds all added, so add them to THStack
        if(nthBkg == bkgTypeN[bkgTypes[i]])
        {
            //cout << bkgTypes[i] << " " << bkgTypeN[bkgTypes[i]] << endl;
            htemp->SetFillColor(bkgColors[bkgTypes[i]]);
            hs.Add(htemp);
            hMCtotal.Add(htemp);

            if(sysName == "")
                leg.AddEntry(htemp, bkgTypes[i], "F");

            isFirstBkg = true;
            nthBkg = 0;
        }
    }
}
void ISRUnfold::doStatUnfold()
{
    //double tauMin=1.e-4;
    //double tauMax=1.e-3;

    nScan=50;
    rhoLogTau=0;
    lCurve=0;
    iBest = 0;

    nScan_mass=50;
    rhoLogTau_mass=0;
    lCurve_mass=0;
    iBest_mass = 0;

    if(regMode == TUnfold::kRegModeNone)
    {
        for(int istat = 0; istat < statSize; istat++)
        {
            TH1* tempMassInput;
            TH1* tempPtInput;

            TString nth_;
            nth_.Form("%d", istat);
            tempPtInput = nomPtUnfold->GetInput("tempPtHist_" + nth_, 0, 0, 0, false);
            tempMassInput = nomMassUnfold->GetInput("tempMassHist_" + nth_, 0, 0, 0, false);

            // randomize histogram bin content
            for(int ibin = 1; ibin<tempPtInput->GetNbinsX()+1;ibin++)
            {
                double err = tempPtInput->GetBinError(ibin);
                if(err > 0.0)
                {
                    tempPtInput->SetBinContent(ibin, tempPtInput->GetBinContent(ibin) + gRandom->Gaus(0,err));
                }
            }
            for(int ibin = 1; ibin<tempMassInput->GetNbinsX()+1;ibin++)
            {
                double err = tempMassInput->GetBinError(ibin);
                if(err > 0.0)
                {
                    tempMassInput->SetBinContent(ibin, tempMassInput->GetBinContent(ibin) + gRandom->Gaus(0,err));
                }
            }
            statPtUnfold.at(istat)->SetInput(tempPtInput, nominal_bias);
            statMassUnfold.at(istat)->SetInput(tempMassInput, nominal_bias);

            statPtUnfold.at(istat)->DoUnfold(0);
            statMassUnfold.at(istat)->DoUnfold(0);

            fillPtStatVariationHist(istat);
            fillMassStatVariationHist(istat);

            delete tempMassInput;
            delete tempPtInput;
            delete statPtUnfold.at(istat);
            delete statMassUnfold.at(istat);
        }
    }
}

void ISRUnfold::doISRUnfold(bool doSys){

    //double tauMin=1.e-4;
    //double tauMax=1.e-3;

    nScan=50;
    rhoLogTau=0;
    lCurve=0;
    iBest = 0;

    nScan_mass=50;
    rhoLogTau_mass=0;
    lCurve_mass=0;
    iBest_mass = 0;

    if(!doSys)
    {
        // Nominal unfolding
        if(regMode == TUnfold::kRegModeNone)
        {
            nomPtUnfold->DoUnfold(0);
            nomMassUnfold->DoUnfold(0);
        }
        else
        {
            // Regularization, use ScanTau() as a default method to find tau
            // ScanLcurve
            //iBest=nomPtUnfold->ScanTau(nScan,0.,0.,&rhoLogTau,
            //                           TUnfoldDensity::kEScanTauRhoAvgSys,
            //                           0,0,
            //                           &lCurve);

            nomPtUnfold->DoUnfold(0.); // Tau
            nomMassUnfold->DoUnfold(0);
            //iBest_mass=nomMassUnfold->ScanTau(nScan_mass,0.,0.,&rhoLogTau_mass,
            //                           TUnfoldDensity::kEScanTauRhoAvgSys,
            //                           0,0,
            //                           &lCurve_mass);
        }
    }
    else
    {
        cout << "Do systematic unfold!" << endl;
        // For systematic
        std::map<TString, std::vector<TString>>::iterator it = sysMap.begin();
        while(it != sysMap.end())
        {
            cout << "Unfold for " << it->first << " systematic." << endl;
            int size = (it->second).size();
            cout << size << " systematic variation exist." << endl;

            for(int i = 0; i < size; i++)
            {
                cout << "posfix: " << (it->second).at(i) << endl;
                sysPtUnfold[it->first][(it->second).at(i)]->DoUnfold(0);
                sysMassUnfold[it->first][(it->second).at(i)]->DoUnfold(0);
            }
            it++;
        }
    }// Unfold for systematic

}

int ISRUnfold::setMeanMass()
{
    cout << "ISRUnfold::setMeanMass()   Save mean of dilepton..." << endl;

    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;
    const Double_t* massBins = temp_tvecd->GetMatrixArray();

    TUnfoldDensityV17* p_unfold = NULL;
    p_unfold = nomMassUnfold;

    TH1 * hdetector_mass = p_unfold->GetInput("hdetector_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);
    TH1* hunfolded_mass =  p_unfold->GetOutput("hunfolded_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);
    TH1* hMC_mass =  p_unfold->GetBias("hunfolded_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);

    // Loop over mass bins
    for(int ibin = 0; ibin < nMassBin; ibin++)
    {
        // Set x-axis range
        hdetector_mass->GetXaxis()->SetRange(hdetector_mass->GetXaxis()->FindBin(massBins[ibin]+0.01), hdetector_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
        hunfolded_mass->GetXaxis()->SetRange(hunfolded_mass->GetXaxis()->FindBin(massBins[ibin]+0.01), hunfolded_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
        hMC_mass->GetXaxis()->SetRange(hunfolded_mass->GetXaxis()->FindBin(massBins[ibin]+0.01), hunfolded_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));

        {
            // Get mean values
            //cout << "Detector, " << ibin << " th mass bin, mean: " << hdetector_mass->GetMean() << " +/- " << hdetector_mass->GetMeanError() << endl;
            meanMass_data_folded. push_back(hdetector_mass->GetMean());
            meanMassStatErr_data_folded.push_back(hdetector_mass->GetMeanError());

            //cout << "Unfolded, " << ibin << " th mass bin, mean: " << hunfolded_mass->GetMean() << " +/- " << hunfolded_mass->GetMeanError() << endl;
            meanMass_data_unfoled.   push_back(hunfolded_mass->GetMean());
            meanMassStatErr_data_unfoled.push_back(hunfolded_mass->GetMeanError());

            //cout << "MC, " << ibin << " th mass bin, mean: " << hMC_mass->GetMean() << " +/- " << hMC_mass->GetMeanError() << endl;
            meanMass_mc_unfoled.   push_back(hMC_mass->GetMean());
            meanMassStatErr_mc_unfoled.push_back(hMC_mass->GetMeanError());
        }

    }// End of mass bin loop

    delete hdetector_mass;
    delete hunfolded_mass;
    delete hMC_mass;

    return nMassBin;
}

int ISRUnfold::setSysMeanMass()
{
    cout << "ISRUnfold::setSysMeanMass()   Save mean of dilepton..." << endl;

    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;
    const Double_t* massBins = temp_tvecd->GetMatrixArray();

    // For systematic
    std::map<TString, std::vector<TString>>::iterator it = sysMap.begin();
    while(it != sysMap.end())
    {
        TUnfoldDensityV17* p_unfold = NULL;

        cout << "Unfold for " << it->first << " systematic." << endl;
        int size = (it->second).size();
        cout << size << " systematic variation exist." << endl;

        for(int i = 0; i < size; i++)
        {
            p_unfold = sysMassUnfold[it->first][(it->second).at(i)];
            TH1* hunfolded_mass =  p_unfold->GetOutput("hunfolded_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);

            // Loop over mass bins
            for(int ibin = 0; ibin < nMassBin; ibin++)
            {
                // Set x-axis range
                hunfolded_mass->GetXaxis()->SetRange(hunfolded_mass->GetXaxis()->FindBin(massBins[ibin]+0.01), hunfolded_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
                meanMass_data_unfolded_sysVariation[it->first][(it->second).at(i)].push_back(hunfolded_mass->GetMean());
                cout << it->first << " " << (it->second).at(i) << " " << hunfolded_mass->GetMean() << endl;
            }// End of mass bin loop
            delete hunfolded_mass;
        }
        it++;
    }
    return nMassBin;
}

void ISRUnfold::setSysError()
{
    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;

    // For systematic
    std::map<TString, std::vector<TString>>::iterator it = sysMap.begin();
    while(it != sysMap.end())
    {

        for(int ibin = 0; ibin < nMassBin; ibin++)
        {
            double error_mass = 0.;
            double error_pt = 0.;
            int size = (it->second).size();
            for(int i = 0; i < size; i++)
            {
                double temp_error_mass = fabs(meanMass_data_unfolded_sysVariation[it->first][(it->second).at(i)].at(ibin) - meanMass_data_unfoled.at(ibin));
                double temp_error_pt = fabs(meanPt_data_unfolded_sysVariation[it->first][(it->second).at(i)].at(ibin) - meanPt_data_unfoled.at(ibin));
                error_mass = error_mass < temp_error_mass?temp_error_mass:error_mass;
                error_pt = error_pt < temp_error_pt?temp_error_pt:error_pt;

            }
            meanMass_data_folded_systematic[it->first].push_back(error_mass);
            meanPt_data_folded_systematic[it->first].push_back(error_pt);
            cout << "Systematic: " << it->first << " " << ibin << " mass bin " << endl;
            cout << "mass: " << error_mass << " pt: " << error_pt << endl;
        }
        it++;
    }
}

void ISRUnfold::setTotSysError()
{

    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;

    // For systematic
    for(int ibin = 0; ibin < nMassBin; ibin++)
    {
        double totSys_mass = 0.;
        double totSys_pt = 0.;
        std::map<TString, std::vector<TString>>::iterator it = sysMap.begin();
        while(it != sysMap.end())
        {
            double sysErr_mass = meanMass_data_folded_systematic[it->first].at(ibin);
            double sysErr_pt = meanPt_data_folded_systematic[it->first].at(ibin);

            totSys_mass += pow(sysErr_mass, 2);
            totSys_pt += pow(sysErr_pt, 2);
            it++;
        }
        meanMassSysErr_data_unfoled.push_back(sqrt(totSys_mass));
        meanPtSysErr_data_unfoled.push_back(sqrt(totSys_pt));
    }
}

void ISRUnfold::setStatError()
{
    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;

    // Loop over mass bins
    for(int ibin = 0; ibin < nMassBin; ibin++)
    {
        meanMassStatErr_data_unfoled.push_back(meanPtStatVariation.at(ibin)->GetRMS());
        meanPtStatErr_data_unfoled.push_back(meanPtStatVariation.at(ibin)->GetRMS());
    }
}

double ISRUnfold::getDetMeanMass(int ibin)
{

    int size = meanMass_data_folded.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanMass_data_folded.at(ibin);
    }
}

double ISRUnfold::getDetMeanMassError(int ibin)
{

    int size = meanMass_data_folded.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanMassStatErr_data_folded.at(ibin);
    }
}

double ISRUnfold::getUnfMeanMass(int ibin)
{

    int size = meanMass_data_unfoled.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanMass_data_unfoled.at(ibin);
    }
}


double ISRUnfold::getUnfMeanMassError(int ibin)
{

    int size = meanMassStatErr_data_unfoled.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanMassStatErr_data_unfoled.at(ibin);
    }
}

double ISRUnfold::getUnfMeanMassSysError(int ibin)
{

    int size = meanMassStatErr_data_unfoled.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanMassSysErr_data_unfoled.at(ibin);
    }
}

double ISRUnfold::getMCGenMeanMass(int ibin)
{

    int size = meanMass_mc_unfoled.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanMass_mc_unfoled.at(ibin);
    }
}

void ISRUnfold::fillPtStatVariationHist(int istat)
{
    cout << "ISRUnfold::fillPtStatVariationHist()  " << endl;

    // Find number of mass bins
    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;

    // Save mean pt
    for(int i = 0; i < nMassBin; i++)
    {
        TString ibinMass;
        ibinMass.Form("%d", i);

        if(istat == 0)
        {
            meanPtStatVariation.push_back(new TH1F("MeanPtStat_bin"+ibinMass, "MeanPtStat_bin"+ibinMass, 40, meanPt_data_unfoled.at(i)-1., meanPt_data_unfoled.at(i)+1.));
        }

        TH1* hpt_temp_data;

        // Get histograms to set mean values
        hpt_temp_data = statPtUnfold.at(istat)->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
        meanPtStatVariation.at(i)->Fill(hpt_temp_data->GetMean());

        delete hpt_temp_data;
    }
}

void ISRUnfold::fillMassStatVariationHist(int istat)
{

    cout << "ISRUnfold::setMeanMass()   Save mean of dilepton..." << endl;

    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;
    const Double_t* massBins = temp_tvecd->GetMatrixArray();

    TH1* hunfolded_mass = statMassUnfold.at(istat)->GetOutput("hunfolded_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);

    // Loop over mass bins
    for(int ibin = 0; ibin < nMassBin; ibin++)
    {
        TString ibinMass;
        ibinMass.Form("%d", ibin);

        // set x-axis range
        hunfolded_mass->GetXaxis()->SetRange(hunfolded_mass->GetXaxis()->FindBin(massBins[ibin]+0.01), hunfolded_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));

        if(istat == 0)
        {
            meanMassStatVariation.push_back(new TH1F("MeanMassStat_bin"+ibinMass, "MeanMassStat_bin"+ibinMass, 80, meanMass_data_unfoled.at(ibin)-2., meanMass_data_unfoled.at(ibin)+2.));
        }
        meanMassStatVariation.at(ibin)->Fill(hunfolded_mass->GetMean());
    }// end of mass bin loop

    delete hunfolded_mass;
}

int ISRUnfold::setMeanPt()
{
    cout << "ISRUnfold::setMeanPt()   Save mean of dilepton momentum..." << endl;

    // Find number of mass bins
    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;

    TUnfoldDensityV17* p_unfold = NULL;
    p_unfold = nomPtUnfold;

    // Save mean pt
    for(int i = 0; i < nMassBin; i++)
    {
        TString ibinMass;
        ibinMass.Form("%d", i);

        TH1* hpt_temp_data;
        TH1* hpt_temp_mc;

        // get histograms to set mean values
        TH1* hdetector_data = p_unfold->GetInput("h_folded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
        hpt_temp_data = p_unfold->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
        hpt_temp_mc   = p_unfold->GetBias("histMC_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);

        //cout << "Detector, " << i << " th mass bin, mean: " << hdetector_data->GetMean() << " +/- " << hdetector_data->GetMeanError() << endl;
        meanPt_data_folded.push_back(hdetector_data->GetMean());
        meanPtStatErr_data_folded.push_back(hdetector_data->GetMeanError());

        //cout << "Unfolded, " << i << " th mass bin, mean: " << hpt_temp_data->GetMean() << " +/- " << hpt_temp_data->GetMeanError() << endl;
        meanPt_data_unfoled.push_back(hpt_temp_data->GetMean());
        meanPtStatErr_data_unfoled.push_back(hpt_temp_data->GetMeanError());

        meanPt_mc_unfoled.push_back(hpt_temp_mc->GetMean());
        meanPtStatErr_mc_unfoled.push_back(hpt_temp_mc->GetMeanError());

        delete hdetector_data;
        delete hpt_temp_data;
        delete hpt_temp_mc;
    }

    return nMassBin;
}

int ISRUnfold::setSysMeanPt()
{
    cout << "ISRUnfold::setSysMeanPt()   Save mean of dilepton momentum..." << endl;

    // Find number of mass bins
    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;

    std::map<TString, std::vector<TString>>::iterator it = sysMap.begin();
    while(it != sysMap.end())
    {
        TUnfoldDensityV17* p_unfold = NULL;

        cout << "Unfold for " << it->first << " systematic." << endl;
        int size = (it->second).size();
        cout << size << " systematic variation exist." << endl;

        for(int i = 0; i < size; i++)
        {
            p_unfold = sysPtUnfold[it->first][(it->second).at(i)];
            // Save mean pt
            for(int j = 0; j < nMassBin; j++)
            {
                TString ibinMass;
                ibinMass.Form("%d", j);

                TH1* hunfolded_pt;

                // Get histograms to set mean values
                hunfolded_pt = p_unfold->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
                meanPt_data_unfolded_sysVariation[it->first][(it->second).at(i)].push_back(hunfolded_pt->GetMean());
                cout << it->first << " " << (it->second).at(i) << " " << hunfolded_pt->GetMean() << endl;
                delete hunfolded_pt;
            }
        }
        it++;
    }

    return nMassBin;
}

double ISRUnfold::getDetMeanPt(int ibin)
{

    int size = meanPt_data_folded.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanPt_data_folded.at(ibin);
    }
}

double ISRUnfold::getDetMeanPtError(int ibin)
{

    int size = meanPt_data_folded.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanPtStatErr_data_folded.at(ibin);
    }
}

double ISRUnfold::getUnfMeanPt(int ibin)
{

    int size = meanPt_data_unfoled.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanPt_data_unfoled.at(ibin);
    }
}

double ISRUnfold::getUnfMeanPtError(int ibin)
{

    int size = meanPtStatErr_data_unfoled.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanPtStatErr_data_unfoled.at(ibin);
    }
}

double ISRUnfold::getUnfMeanPtSysError(int ibin)
{

    int size = meanPtStatErr_data_unfoled.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanPtSysErr_data_unfoled.at(ibin);
    }
}

double ISRUnfold::getMCGenMeanPt(int ibin)
{

    int size = meanPt_mc_unfoled.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanPt_mc_unfoled.at(ibin);
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

TH1* ISRUnfold::getMCHists(TString var, TString outHistName, TString steering, bool useAxis)
{
    if(var == "Mass")
        return nomMassUnfold->GetBias(outHistName,0,0,steering,useAxis);
    else
        return nomPtUnfold->GetBias(outHistName,0,0,steering,useAxis);
}

TH1* ISRUnfold::getDetHists(TString var, TString outHistName, TString steering, bool useAxis)
{
    // Return background subtracted data

    if(var == "Mass")
        return nomMassUnfold->GetInput(outHistName,0,0,steering,useAxis);
    else
        return nomPtUnfold->GetInput(outHistName,0,0,steering,useAxis);
}

TH1* ISRUnfold::getRawHist(TString var, TString filePath, TString dirName, TString histName, TString outHistName, TString steering, bool useAxis, bool divBinWidth)
{
    TH1::AddDirectory(kFALSE);
    TFile* filein = new TFile(filePath);
    TH1* hist = NULL;

    if(steering != "")
    {
        //cout << dirName+"/"+var+"/"+histName << endl;
        TH1* raw_hist = (TH1*)filein->Get(dirName+"/"+var+"/"+histName);
        //cout << "# of bins: " << raw_hist->GetNbinsX() << endl;
        if(histName.Contains("DYJetsTo") && !histName.Contains("Tau"))
        {
            histName.ReplaceAll("DYJetsTo", "DYJets10to50To");
            raw_hist->Add((TH1*)filein->Get(dirName+"/"+var+"/"+histName));
        }

        if(var.Contains("Pt"))
        {
            hist = pt_binning_Rec->ExtractHistogram(outHistName, raw_hist, 0, useAxis, steering);
            //cout << "steering: " << steering << endl;
            //cout << "# of bins: " << hist->GetNbinsX() << endl;
        }
        else
        {
            hist = mass_binning_Rec->ExtractHistogram(outHistName, raw_hist, 0, useAxis, steering);
        }

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

void ISRUnfold::drawStatVariation(bool isPt, int massBin)
{
    TCanvas* c = new TCanvas("c","c", 800, 800);
    c->SetLogy();
    c->cd();

    TString nth;
    nth.Form("%d", massBin);
    TString year_;
    year_.Form("%d", year);

    if(isPt)
    {
        meanPtStatVariation.at(massBin)->Draw("pe");
        c->SaveAs("MeanPtStat_" + nth + year_ + ".pdf");
    }
    else
    {
        meanMassStatVariation.at(massBin)->Draw("pe");
        c->SaveAs("MeanMassStat_" + nth + year_ + ".pdf");
    }
    delete c;
}
