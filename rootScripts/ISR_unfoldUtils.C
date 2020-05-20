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

void ISRUnfold::setNomResMatrix(TString var, TString filepath, TString dirName, TString histName, bool isSquareMatrix)
{
    cout << "ISRUnfold::setNomResMatrix set response matrix..." << endl;
    TFile* filein = new TFile(filepath);

    TString Rec_binName = "Rec_"+var;
    TString Gen_binName = "Gen_"+var;
    Rec_binName = dirName + "/" + var + "_ResMatrix_" + histName + "/" + Rec_binName;

    // Use the same bin definition for both rec and gen histogram
    if(isSquareMatrix)
    {
        Rec_binName = dirName + "/" + var + "_ResMatrix_" + histName + "/" + Gen_binName;
    }
    Gen_binName = dirName + "/" + var + "_ResMatrix_" + histName + "/" + Gen_binName;

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
    hmcGenRec = (TH2*)filein->Get(dirName + "/" + var + "_ResMatrix_" + histName +"/hmc" + var + "GenRecnominal");
    
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

void ISRUnfold::setSysTUnfoldDensity(TString var, TString filepath, TString dirName, TString histName, TString sysName, TString sysPostfix, int nth)
{
    TFile* filein = new TFile(filepath);
    TH2* hmcGenRec;
    hmcGenRec = (TH2*)filein->Get(dirName + "/" + var + "_ResMatrix_" + histName +"/hmc" + var + "GenRec_" + sysPostfix);

    if( var == "Pt" )
    {
        sysPtUnfold[sysName].push_back( new TUnfoldDensityV17(hmcGenRec,
                                                              TUnfold::kHistMapOutputHoriz,
                                                              regMode,
                                                              TUnfold::kEConstraintArea,
                                                              TUnfoldDensityV17::kDensityModeBinWidth,
                                                              pt_binning_Gen,pt_binning_Rec));
    }
    else if( var == "Mass" )
    {
        sysMassUnfold[sysName].push_back( new TUnfoldDensityV17(hmcGenRec,
                                                                TUnfold::kHistMapOutputHoriz,
                                                                regMode,
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
void ISRUnfold::setUnfInput(TString var, TString filepath, TString dirName, TString histName, bool isSys, TString sysName, int nth)
{
    TH1::AddDirectory(kFALSE);

    TFile* filein = new TFile(filepath);
    TString nth_;
    nth_.Form("%d", nth);
    TH1* hRec = NULL;

    // Nominal histograms
    // TODO Use input covariance matrix
    if(!isSys)
    {
        hRec = (TH1*)filein->Get(dirName+"/"+var+"/"+histName);

        // Use DY MC as unfolding input, i.e. simple closure test
        if(histName.Contains("DYJetsTo"))
        {
            histName.ReplaceAll("DYJetsTo", "DYJets10to50To");
            hRec->Add((TH1*)filein->Get(dirName+"/"+var+"/"+histName));
        }

        if(var == "Pt")
        {
            nomPtUnfold->SetInput(hRec,   nominal_bias);
        }
        else if(var == "Mass")
        {
            nomMassUnfold->SetInput(hRec, nominal_bias);
        }
        else{
            cout << "ISRUnfold::setUnfInput, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
        }
    }
    // Systematic histograms
    else
    {
        hRec = (TH1*)filein->Get(dirName+"/"+var+"/"+histName);
        // For systematic, using the same input histogram as nominal, unless data changed in a systematic change
        if(var == "Pt")
        {
            sysPtUnfold[sysName].at(nth)->SetInput(hRec,   nominal_bias);
        }
        else if(var == "Mass")
        {
            sysMassUnfold[sysName].at(nth)->SetInput(hRec,   nominal_bias);
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

void ISRUnfold::subBkgs(TString filepath, std::pair<TString, TString>& bkgInfo, bool isSys, TString sysName, int totSysN, int nth, TString dirName, double bkgScale)
{

    TFile* filein = new TFile(filepath);
    TH1* hPtRec = NULL;
    TH1* hMassRec = NULL;
    
    double bkg_scale = bkgScale;
    TString systematic_postfix = sysName;
    
    if(totSysN == 2)
    {
        if(nth == 0 ) systematic_postfix="_"+systematic_postfix+"Up";
        if(nth == 1 ) systematic_postfix="_"+systematic_postfix+"Down";
    }
    
    if(totSysN == 6 && (sysName == "Scale" || sysName =="pdfScale"))
    {
        if(nth == 0 ) systematic_postfix="_"+systematic_postfix+"AUp";
        if(nth == 1 ) systematic_postfix="_"+systematic_postfix+"ADown";
        if(nth == 2 ) systematic_postfix="_"+systematic_postfix+"BUp";
        if(nth == 3 ) systematic_postfix="_"+systematic_postfix+"BDown";
        if(nth == 4 ) systematic_postfix="_"+systematic_postfix+"ABUp";
        if(nth == 5 ) systematic_postfix="_"+systematic_postfix+"ABDown";
    }
    
    if(totSysN == 100 && sysName == "PDFerror")
    {
        TString nth_;
        nth_.Form ("%03d", nth);
        systematic_postfix="_"+systematic_postfix+nth_;
    }
    
    // Nominal histograms
    if(!isSys)
    {

        bkgNames.push_back(bkgInfo.first);
        bkgTypes.push_back(bkgInfo.second);

        hPtRec = (TH1*)filein->Get(dirName + "/Pt/histo_" + bkgInfo.first + "_nominal");
        nomPtUnfold->  SubtractBackground(hPtRec, bkgInfo.first, bkg_scale);

        hMassRec = (TH1*)filein->Get(dirName + "/Mass/histo_" + bkgInfo.first + "_nominal");
        nomMassUnfold->SubtractBackground(hMassRec, bkgInfo.first, bkg_scale);
    }

/*
    // Systematic histograms
    else
    {	
    
        // Systematics using nominal detector distributions
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
 */   
    filein->Close();
    delete filein;
}

void ISRUnfold::setSystematics(TString sysName, TString sysHistName)
{
    sysMap[sysName].push_back(sysHistName);
}

// Draw detector distributions using input root file
TCanvas* ISRUnfold::drawFoldedHists(TString var, TString filePath, TString steering, bool useAxis, TString sysName)
{
    TH1::AddDirectory(kFALSE);
    cout << "ISRUnfold::drawFoldedHists, Draw plot!" << endl; 

    TFile* filein = new TFile(filePath);

    TH1* hData = NULL;
    TH1* hDY = NULL;
    TH1* hMCtotal = NULL;
    TH1* hRatio = NULL;
    // For systematic
    TH1* hDY_up = NULL;
    TH1* hMCtotal_up = NULL;
    TH1* hRatio_up = NULL;
    TH1* hDY_down = NULL;
    TH1* hMCtotal_down = NULL;
    TH1* hRatio_down = NULL;

    hData = getRawHist(var, filePath, "Detector", "histo_DoubleMuon_nominal", "Data", steering, useAxis);
    hDY = getRawHist(var, filePath, "Detector", "histo_DYJetsToMuMu_nominal", "Signal", steering, useAxis);
    hMCtotal = (TH1*) hDY->Clone("hMCtotal");
    hRatio = (TH1*) hData->Clone("hRatio");

    if(sysName != "")
    {
        hDY_up = getRawHist(var, filePath, "Detector", "histo_DYJetsToMuMu_" + sysMap[sysName][0], "Signal"+sysMap[sysName][0], steering, useAxis);
        hMCtotal_up = (TH1*) hDY_up->Clone("hMCtotal_up");
        hRatio_up = (TH1*) hData->Clone("hRatio_up");

        hDY_down = getRawHist(var, filePath, "Detector", "histo_DYJetsToMuMu_" + sysMap[sysName][1], "Signal"+sysMap[sysName][1], steering, useAxis);
        hMCtotal_down = (TH1*) hDY_down->Clone("hMCtotal_down");
        hRatio_down = (TH1*) hData->Clone("hRatio_down");
    }

    setTDRStyle();
    writeExtraText = true;
    extraText  = "work in progress";

    TCanvas* c_out = new TCanvas("detector_level_"+var, "detector_level_"+var, 50, 50, 1600, 1400);
    c_out->Draw();
    c_out->cd();

    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetBottomMargin(0.01);
    pad1->SetTopMargin(0.1);
    pad1->SetTicks(1);
    pad1->SetLogy();
    if(var=="Mass")
        pad1->SetLogx();
    pad1->Draw();
    pad1->cd();

    hData->SetTitle("");
    hData->SetStats(false);
    hData->Draw("p9histe");
    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(.7);
    hData->SetLineColor(kBlack);
    hData->GetYaxis()->SetTitle("Events/bin");
    hData->SetMaximum(1e8);
    hData->SetMinimum(1e-1);

    hDY->SetFillColor(kYellow);

    THStack* hsMC = new THStack("hsMC", "hsMC");
    setTHStack(var, filePath, *hsMC, *hMCtotal, steering, useAxis);
    hsMC->Add(hDY);

    THStack* hsMC_up;
    THStack* hsMC_down;
    if(sysName != "")
    {
        hsMC_up = new THStack("hsMC_up", "hsMC_up");
        hsMC_down = new THStack("hsMC_down", "hsMC_down");
        setTHStack(var, filePath, *hsMC_up, *hMCtotal_up, steering, useAxis);
        setTHStack(var, filePath, *hsMC_down, *hMCtotal_down, steering, useAxis);
    }
    
    hsMC->Draw("hist same");
    hData->Draw("p9histe same");
    pad1->RedrawAxis();

    TLine massEdgeLine;
    if(var=="Mass")
    {
        for(int i = 0; i < massBinEdges.size(); i++)
        {
            if(i==0) continue;
            if(i==massBinEdges.size()-1) continue;
            massEdgeLine.DrawLine(massBinEdges[i], hData->GetMinimum(), massBinEdges[i], hData->GetMaximum());
        }
    }
    c_out->cd();

    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.2);
    if(var=="Mass")
        pad2->SetLogx();
    pad2->SetTicks(1);
    pad2->SetGridy(1);
    pad2->Draw();
    pad2->cd();
   
    hRatio->SetStats(false);  
    hRatio->Divide(hMCtotal);
    hRatio->Draw("p9histe");
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(0.7);
    hRatio->SetLineColor(kBlack);
    hRatio->GetYaxis()->SetTitle("Data/MC");

    hRatio->SetMinimum(0.7);
    hRatio->SetMaximum(1.3);

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
        sysBand_ratio->SetFillColorAlpha(kBlack,0.8);
        sysBand_ratio->SetFillStyle(3004);
        sysBand_ratio->SetMarkerSize(0.);
        sysBand_ratio->Draw("E2 same");
    }

    c_out->cd();
    CMS_lumi(c_out, 7, 11);
    c_out->SaveAs("detector_"+var+".png");

    delete filein;

    return c_out;
}

void ISRUnfold::setTHStack(TString var, TString filePath, THStack& hs, TH1& hMCtotal, TString steering, bool useAxis, TString sysName)
{
    TH1::AddDirectory(kFALSE); 
    int bkgSize = bkgNames.size();

    // Count total number of each background type N
    map<TString, int> bkgTypeN;
    cout << "N bkg: " << bkgSize << endl;
    for(int i = 0; i < bkgSize; i++)
    {
        map<TString, int>::iterator it = bkgTypeN.find(bkgTypes[i]);
        if(it != bkgTypeN.end())
        {
            bkgTypeN[bkgTypes[i]]++;
        }
        else
        {
            cout << bkgTypes[i] << " first found" << endl;
            bkgTypeN[bkgTypes[i]] = 1;
        }        
    }

    TH1* htemp = NULL; 
    bool isFirstBkg = true;
    int nthBkg = 0;

    for(int i = 0; i < bkgSize; i++)
    {
        if(isFirstBkg)
        {
            if(sysName == "")
                htemp = getRawHist(var, filePath, "Detector", "histo_"+bkgNames[i]+"_nominal", "h"+bkgNames[i], steering, useAxis);  
            else
                htemp = getRawHist(var, filePath, "Detector", "histo_"+bkgNames[i]+"_"+sysName, "h"+bkgNames[i], steering, useAxis);  
            isFirstBkg = false;
            nthBkg++;
        }
        else
        {
            if(sysName == "")
                htemp->Add(getRawHist(var, filePath, "Detector", "histo_"+bkgNames[i]+"_nominal", "h"+bkgNames[i], steering, useAxis));  
            else
                htemp->Add(getRawHist(var, filePath, "Detector", "histo_"+bkgNames[i]+"_"+sysName, "h"+bkgNames[i], steering, useAxis));  
            nthBkg++;
        }

        // This type of backgrounds all added, so add them to THStack
        if(nthBkg == bkgTypeN[bkgTypes[i]])
        {
            cout << bkgTypes[i] << " " << bkgTypeN[bkgTypes[i]] << endl;
            htemp->SetFillColor(bkgColors[bkgTypes[i]]); 
            hs.Add(htemp);
            hMCtotal.Add(htemp);

            isFirstBkg = true;
            nthBkg = 0;
        }
    }
}
void ISRUnfold::doStatUnfold()
{
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

    if(regMode == TUnfold::kRegModeNone)
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

    // For systematic
    if(doSys)
    {
        std::map<TString, std::vector<TUnfoldDensityV17*>>::iterator it;
        std::map<TString, std::vector<TUnfoldDensityV17*>>::iterator it_end;

        // Unfold for pt distribution
        it = sysPtUnfold.begin();
        it_end = sysPtUnfold.end();

        while(it != it_end)
        {
            int nSys = it->second.size();
            for(int i = 0; i < nSys; i++)
            {
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

        // Unfolding for mass distribution
        it = sysMassUnfold.begin();
        it_end = sysMassUnfold.end();

        while(it != it_end)
        {
            int nSys = it->second.size();
            for(int i = 0; i < nSys; i++)
            {
                if((it->first)=="unfoldScan" || (it->first)=="unfoldBias")
                {
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
    }// Unfold for systematic
}

int ISRUnfold::setMeanMass()
{
    cout << "ISRUnfold::setMeanMass()   Save mean of dilepton..." << endl;
    
    const TUnfoldBinningV17* temp_binning_gen_pt = pt_binning_Gen;
    const TVectorD* temp_tvecd = temp_binning_gen_pt->GetDistributionBinning(1);
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
        // set x-axis range
        hdetector_mass->GetXaxis()->SetRange(hdetector_mass->GetXaxis()->FindBin(massBins[ibin]+0.01), hdetector_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
        hunfolded_mass->GetXaxis()->SetRange(hunfolded_mass->GetXaxis()->FindBin(massBins[ibin]+0.01), hunfolded_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
        hMC_mass->GetXaxis()->SetRange(hunfolded_mass->GetXaxis()->FindBin(massBins[ibin]+0.01), hunfolded_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));

        {
            // Get mean values
            //cout << "Detector, " << ibin << " th mass bin, mean: " << hdetector_mass->GetMean() << " +/- " << hdetector_mass->GetMeanError() << endl;
            meanMass_data_detector. push_back(hdetector_mass->GetMean());
            meanMassStatErr_data_detector.push_back(hdetector_mass->GetMeanError());

            //cout << "Unfolded, " << ibin << " th mass bin, mean: " << hunfolded_mass->GetMean() << " +/- " << hunfolded_mass->GetMeanError() << endl;
            meanMass_data_unfoled.   push_back(hunfolded_mass->GetMean());
            meanMassStatErr_data_unfoled.push_back(hunfolded_mass->GetMeanError());

            //cout << "MC, " << ibin << " th mass bin, mean: " << hMC_mass->GetMean() << " +/- " << hMC_mass->GetMeanError() << endl;
            meanMass_mc_unfoled.   push_back(hMC_mass->GetMean());
            meanMassStatErr_mc_unfoled.push_back(hMC_mass->GetMeanError());
        }

    }// end of mass bin loop

    delete hdetector_mass;
    delete hunfolded_mass;
    delete hMC_mass;

    return nMassBin;
}

void ISRUnfold::setStatError()
{
    const TUnfoldBinningV17* temp_binning_gen_pt = pt_binning_Gen;
    const TVectorD* temp_tvecd = temp_binning_gen_pt->GetDistributionBinning(1);
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

double ISRUnfold::getDetMeanMassError(int ibin)
{

    int size = meanMass_data_detector.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);     
    }
    else
    {
        return meanMassStatErr_data_detector.at(ibin);
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
    const TUnfoldBinningV17* temp_binning_gen_pt = pt_binning_Gen;
    const TVectorD* temp_tvecd = temp_binning_gen_pt->GetDistributionBinning(1);
    int nMassBin = nMassBin = temp_tvecd->GetNrows() - 1;

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
        TH1* hpt_temp_mc;

        // Get histograms to set mean values
        hpt_temp_data = statPtUnfold.at(istat)->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
        meanPtStatVariation.at(i)->Fill(hpt_temp_data->GetMean());

        delete hpt_temp_data;
    }
}

void ISRUnfold::fillMassStatVariationHist(int istat)
{

    cout << "ISRUnfold::setMeanMass()   Save mean of dilepton..." << endl;
    
    const TUnfoldBinningV17* temp_binning_gen_pt = pt_binning_Gen;
    const TVectorD* temp_tvecd = temp_binning_gen_pt->GetDistributionBinning(1);
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
// Set mean pt from mass and DY mc
int ISRUnfold::setMeanPt()
{
    cout << "ISRUnfold::setMeanPt()   Save mean of dilepton momentum..." << endl;

    // Find number of mass bins
    const TUnfoldBinningV17* temp_binning_gen_pt = pt_binning_Gen;
    const TVectorD* temp_tvecd = temp_binning_gen_pt->GetDistributionBinning(1);
    int nMassBin = nMassBin = temp_tvecd->GetNrows() - 1;

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
        TH1* hdetector_data = p_unfold->GetInput("h_detector_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
        hpt_temp_data = p_unfold->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
        hpt_temp_mc   = p_unfold->GetBias("histMC_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE); 

        //cout << "Detector, " << i << " th mass bin, mean: " << hdetector_data->GetMean() << " +/- " << hdetector_data->GetMeanError() << endl;
        meanPt_data_detector.push_back(hdetector_data->GetMean());
        meanPtStatErr_data_detector.push_back(hdetector_data->GetMeanError());

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

double ISRUnfold::getDetMeanPtError(int ibin)
{

    int size = meanPt_data_detector.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);     
    }
    else
    {
        return meanPtStatErr_data_detector.at(ibin);
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

TH1* ISRUnfold::getRawHist(TString var, TString filePath, TString dirName, TString histName, TString outHistName, TString steering, bool useAxis)
{
    TH1::AddDirectory(kFALSE);
    TFile* filein = new TFile(filePath);

    TH1* raw_hist = (TH1*)filein->Get(dirName+"/"+var+"/"+histName);
    if(histName.Contains("DYJetsTo") && !histName.Contains("Tau"))
    {
        histName.ReplaceAll("DYJetsTo", "DYJets10to50To");
        raw_hist->Add((TH1*)filein->Get(dirName+"/"+var+"/"+histName));
    }

    TH1* hist;
    if(var=="Pt")
        hist = pt_binning_Rec->ExtractHistogram(outHistName, raw_hist, 0, useAxis, steering);
    else
        hist = mass_binning_Rec->ExtractHistogram(outHistName, raw_hist, 0, useAxis, steering);

    delete filein;
    delete raw_hist;
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
