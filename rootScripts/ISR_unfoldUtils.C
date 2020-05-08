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

// FIXME Update!!
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
        // For systematic, using the same input histogram as nominal, unless data changed in a systematic change
        if(sysName != "lepMom")
        {

            if(sysName != "Stat")
            {
                if(var == "Pt")
                {
                    if(channel_name == "muon")     hRec = (TH1*)filein->Get(dirName + "/hist_ptll/histo_DoubleMuonnominal");
                    if(channel_name == "electron" && year != 2018) hRec = (TH1*)filein->Get(dirName + "/hist_ptll/histo_DoubleEGnominal");
                    if(channel_name == "electron" && year == 2018) hRec = (TH1*)filein->Get(dirName + "/hist_ptll/histo_EGammanominal");
                    sysPtUnfold[sysName].at(nth)  ->SetInput(hRec,   nominal_bias);
                }
                else if(var == "Mass")
                {
                    if(channel_name == "muon")     hRec = (TH1*)filein->Get(dirName + "/hist_mll/histo_DoubleMuonnominal");
                    if(channel_name == "electron" && year != 2018) hRec = (TH1*)filein->Get(dirName + "/hist_mll/histo_DoubleEGnominal");
                    if(channel_name == "electron" && year == 2018) hRec = (TH1*)filein->Get(dirName + "/hist_mll/histo_EGammanominal");
                    sysMassUnfold[sysName].at(nth)->SetInput(hRec,   nominal_bias);
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

                    // Randomize histogram bins
                    for(int ibin = 1; ibin<temp_ptHist->GetNbinsX()+1;ibin++)
                    {
                        double err = temp_ptHist->GetBinError(ibin);
                        if(err > 0.0)
                        {
                            temp_ptHist->SetBinContent(ibin, temp_ptHist->GetBinContent(ibin) + gRandom->Gaus(0,err));
                        }
                    }
                    sysPtUnfold[sysName].at(nth)->SetInput(temp_ptHist, nominal_bias);
                }
                else if(var == "Mass")
                {
                    TH1* temp_massHist;

                    TString nth_;
                    nth_.Form("%d", nth);
                    temp_massHist = nomMassUnfold->GetInput("ptToy_" + nth_, 0, 0, 0, false);

                    // Randomize histogram bins
                    for(int ibin = 1; ibin<temp_massHist->GetNbinsX()+1;ibin++)
                    {
                        double err = temp_massHist->GetBinError(ibin);
                        if(err > 0.0)
                        {
                            temp_massHist->SetBinContent(ibin, temp_massHist->GetBinContent(ibin) + gRandom->Gaus(0,err));
                        }
                    }
                    sysMassUnfold[sysName].at(nth)->SetInput(temp_massHist, nominal_bias);
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
                dirName += "_lepMomUp";
                histDirPostfix = "_lepMomUp";
            }
            else if(nth == 1)
            {
                dirName += "_lepMomDown";
                histDirPostfix = "_lepMomDown";
            }
            else
            {
                exit(EXIT_FAILURE);
            }

            if(var == "Pt")
            {
                if(channel_name == "muon")     hRec = (TH1*)filein->Get(dirName + "/hist_ptll" + histDirPostfix + "/histo_DoubleMuonnominal");
                if(channel_name == "electron" && year != 2018) hRec = (TH1*)filein->Get(dirName + "/hist_ptll" + histDirPostfix + "/histo_DoubleEGnominal");
                if(channel_name == "electron" && year == 2018) hRec = (TH1*)filein->Get(dirName + "/hist_ptll" + histDirPostfix + "/histo_EGammanominal");
                sysPtUnfold[sysName].at(nth)  ->SetInput(hRec,   nominal_bias);
            }
            else if(var == "Mass")
            {
                if(channel_name == "muon")     hRec = (TH1*)filein->Get(dirName + "/hist_mll" + histDirPostfix + "/histo_DoubleMuonnominal");
                if(channel_name == "electron" && year != 2018) hRec = (TH1*)filein->Get(dirName + "/hist_mll" + histDirPostfix + "/histo_DoubleEGnominal");
                if(channel_name == "electron" && year == 2018) hRec = (TH1*)filein->Get(dirName + "/hist_mll" + histDirPostfix + "/histo_EGammanominal");
                sysMassUnfold[sysName].at(nth)->SetInput(hRec,   nominal_bias);
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

void ISRUnfold::subBkgs(TString filepath, std::pair<TString, TString>& bkgInfo, bool isSys, TString sysName, int totSysN, int nth, TString dirName)
{

    TFile* filein = new TFile(filepath);
    TH1* hPtRec = NULL;
    TH1* hMassRec = NULL;
    
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
    
    // Nominal histograms
    if(!isSys)
    {

        bkgNames.push_back(bkgInfo.first);
        bkgTypes.push_back(bkgInfo.second);

        hPtRec = (TH1*)filein->Get(dirName + "/Pt/histo_" + bkgInfo.first + "nominal");
        nomPtUnfold->  SubtractBackground(hPtRec, bkgInfo.first, bkg_scale);

        hMassRec = (TH1*)filein->Get(dirName + "/Mass/histo_" + bkgInfo.first + "nominal");
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

// Draw detector distributions using input root file
TCanvas* ISRUnfold::drawFoldedHists(TString var, TString filePath)
{
    TH1::AddDirectory(kFALSE);
    cout << "ISRUnfold::drawFoldedHists, Draw plot!" << endl; 

    TFile* filein = new TFile(filePath);

    TH1* hData = NULL;
    TH1* hDY = NULL;
    TH1* hMCtotal = NULL;
    TH1* hRatio = NULL;

    TString steering = "mass[UO];pt[UOC0]";
    bool useAxis = true; // 
    if(var == "Pt")
    {
        steering = "pt[UO];mass[UO]";
        useAxis = false;
    }
    hData = getRawHist(var, filePath, "Detector", "histo_DoubleMuonnominal", "Data", steering, useAxis);
    hDY = getRawHist(var, filePath, "Detector", "histo_DYJetsToMuMunominal", "Signal", steering, useAxis);
    hMCtotal = (TH1*) hDY->Clone("hMCtotal");
    hRatio = (TH1*) hData->Clone("hRatio");

    TCanvas* c_out = new TCanvas("detector_level", "detector_level", 50, 50, 800, 700);
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
    // TODO stack function
    setTHStack(var, filePath, *hsMC, *hMCtotal);
    hsMC->Add(hDY);
    
    hsMC->Draw("hist same");
    hData->Draw("p9histe same");

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

    c_out->cd();
    c_out->SaveAs("detector_"+var+".pdf");

    delete filein;

    return c_out;
}

void ISRUnfold::setTHStack(TString var, TString filePath, THStack& hs, TH1& hMCtotal)
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

    TString steering = "mass[UO];pt[UOC0]";
    bool useAxis = true; // 
    if(var == "Pt")
    {
        steering = "pt[UO];mass[UO]";
        useAxis = false;
    }

    for(int i = 0; i < bkgSize; i++)
    {
        if(isFirstBkg)
        {
            htemp = getRawHist(var, filePath, "Detector", "histo_"+bkgNames[i]+"nominal", "h"+bkgNames[i], steering, useAxis);  
            isFirstBkg = false;
            nthBkg++;
        }
        else
        {
            htemp->Add(getRawHist(var, filePath, "Detector", "histo_"+bkgNames[i]+"nominal", "h"+bkgNames[i], steering, useAxis));  
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

    // For systematic
    if(doSys)
    {
        std::map<TString, std::vector<TUnfoldDensityV17*>>::iterator it;
        std::map<TString, std::vector<TUnfoldDensityV17*>>::iterator it_end;

        // Unfolding for pt distribution
        {
            it = sysPtUnfold.begin();
            it_end = sysPtUnfold.end();
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

        // Unfolding for mass distribution
        {
            it = sysMassUnfold.begin();
            it_end = sysMassUnfold.end();
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
            meanMass_data_det_unf.   push_back(hunfolded_mass->GetMean());
            meanMassStatErr_data_det_unf.push_back(hunfolded_mass->GetMeanError());

            //cout << "MC, " << ibin << " th mass bin, mean: " << hMC_mass->GetMean() << " +/- " << hMC_mass->GetMeanError() << endl;
            meanMass_mc_det_unf.   push_back(hMC_mass->GetMean());
            meanMassStatErr_mc_det_unf.push_back(hMC_mass->GetMeanError());
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

        {
            //cout << "Detector, " << i << " th mass bin, mean: " << hdetector_data->GetMean() << " +/- " << hdetector_data->GetMeanError() << endl;
            meanPt_data_detector.push_back(hdetector_data->GetMean());
            meanPtStatErr_data_detector.push_back(hdetector_data->GetMeanError());

            //cout << "Unfolded, " << i << " th mass bin, mean: " << hpt_temp_data->GetMean() << " +/- " << hpt_temp_data->GetMeanError() << endl;
            meanPt_data_det_unf.push_back(hpt_temp_data->GetMean());
            meanPtStatErr_data_det_unf.push_back(hpt_temp_data->GetMeanError());

            meanPt_mc_det_unf.push_back(hpt_temp_mc->GetMean());
            meanPtErr_mc_det_unf.push_back(hpt_temp_mc->GetMeanError());
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
