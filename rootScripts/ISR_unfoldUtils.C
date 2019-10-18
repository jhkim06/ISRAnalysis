#include "ISR_unfoldUtils.h"

#include "tdrstyle.C"
#include "CMS_lumi.C"

// TODO use one function for setting response matrix
void ISRUnfold::setNomTUnfoldDensity(TString var, TString filepath, TString phase_name, TString fsr_correction_name){

	TFile* filein = new TFile(filepath);

        // set response matrix
        TH2* hmcGenRec;
        
        if(var == "Pt")
            hmcGenRec = (TH2*)filein->Get(phase_name + "/ptll_rec_gen_" + fsr_correction_name + "_response_matrix/hmc" + var + "GenRecnominal");
        else if(var == "Mass")
            hmcGenRec = (TH2*)filein->Get(phase_name + "/mll_rec_gen_" + fsr_correction_name + "_response_matrix/hmc" + var + "GenRecnominal");
        else{
            cout << "ISRUnfold::setNomTUnfoldDensity, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
        }

        if(do_normalization){
            if(channel_name =="electron") hmcGenRec->Scale(0.959939);
            else hmcGenRec->Scale(0.953026);
        }

        // set binning definition
        TUnfoldBinning* binning_Rec = NULL;
        TUnfoldBinning* binning_Gen = NULL;

        if( var == "Pt" ){

          TString Rec_Pt = "Rec_Pt";
          TString Gen_Pt = "Gen_Pt";

          Rec_Pt = phase_name + "/ptll_rec_gen_" + fsr_correction_name + "_response_matrix/" + Rec_Pt;
          Gen_Pt = phase_name + "/ptll_rec_gen_" + fsr_correction_name + "_response_matrix/" + Gen_Pt;

          binning_Rec = (TUnfoldBinning*)filein->Get(Rec_Pt);
          binning_Gen = (TUnfoldBinning*)filein->Get(Gen_Pt);
        }
        else if( var == "Mass" ){

          TString Rec_Mass = "Rec_Mass";
          TString Gen_Mass = "Gen_Mass";

          Rec_Mass = phase_name + "/mll_rec_gen_" + fsr_correction_name + "_response_matrix/" + Rec_Mass;
          Gen_Mass = phase_name + "/mll_rec_gen_" + fsr_correction_name + "_response_matrix/" + Gen_Mass;

          binning_Rec = (TUnfoldBinning*)filein->Get(Rec_Mass);
          binning_Gen = (TUnfoldBinning*)filein->Get(Gen_Mass);

        }
        else{
            cout << "ISRUnfold::setNomTUnfoldDensity, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
        }

	if( var == "Pt" ){ 
        	nomPtUnfold = new TUnfoldDensityV17(hmcGenRec,
        	                               TUnfold::kHistMapOutputHoriz,
        	                               TUnfold::kRegModeNone, // fixed to use no regularisation temporary
        	                               TUnfold::kEConstraintArea,
        	                               TUnfoldDensityV17::kDensityModeBinWidth,
        	                               binning_Gen,binning_Rec);
	}
        else if( var == "Mass" ){
                nomMassUnfold = new TUnfoldDensityV17(hmcGenRec,
                                               TUnfold::kHistMapOutputHoriz,
                                               TUnfold::kRegModeNone, // fixed to use no regularisation temporary
                                               TUnfold::kEConstraintArea,
                                               TUnfoldDensityV17::kDensityModeBinWidth,
                                               binning_Gen,binning_Rec);
        }
        else{
            cout << "ISRUnfold::setNomTUnfoldDensity, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
        }
}

void ISRUnfold::setNomFSRTUnfoldDensity(TString var, TString filepath, TString phase_name, TString fsr_correction_name){

        TFile* filein = new TFile(filepath);

        // set response matrix
        TH2* hmcGenGen;

        if(var == "Pt"){
            hmcGenGen = (TH2*)filein->Get(phase_name + "/ptll_gen_post_fsr_" + fsr_correction_name + "_response_matrix/hmc" + var + "GenRecnominal");
            cout << phase_name + "/ptll_gen_post_fsr_" + fsr_correction_name + "_response_matrix/hmc" + var + "GenRecnominal" << endl;
        }
        else if(var == "Mass")
            hmcGenGen = (TH2*)filein->Get(phase_name + "/mll_gen_post_fsr_" + fsr_correction_name + "_response_matrix/hmc" + var + "GenRecnominal");
        else{
            cout << "ISRUnfold::setNomTUnfoldDensity, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
        }

        if(do_normalization){
            if(channel_name =="electron") hmcGenGen->Scale(0.959939);
            else hmcGenGen->Scale(0.953026);
        }

        // set binning definition
        TUnfoldBinning* binning_Gen = NULL;

        if( var == "Pt" ){

          TString Gen_Pt = "Gen_Pt";

          Gen_Pt = phase_name + "/ptll_gen_post_fsr_" + fsr_correction_name + "_response_matrix/" + Gen_Pt;

          binning_Gen = (TUnfoldBinning*)filein->Get(Gen_Pt);
        }
        else if( var == "Mass" ){

          TString Gen_Mass = "Gen_Mass";

          Gen_Mass = phase_name + "/mll_gen_post_fsr_" + fsr_correction_name + "_response_matrix/" + Gen_Mass;

          binning_Gen = (TUnfoldBinning*)filein->Get(Gen_Mass);

        }
        else{
            cout << "ISRUnfold::setNomTUnfoldDensity, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
        }

        if( var == "Pt" ){
                nomPtFSRUnfold = new TUnfoldDensityV17(hmcGenGen,
                                               TUnfold::kHistMapOutputHoriz,
                                               TUnfold::kRegModeNone, // fixed to use no regularisation temporary
                                               TUnfold::kEConstraintArea,
                                               TUnfoldDensityV17::kDensityModeBinWidth,
                                               binning_Gen,binning_Gen);
        }
        else if( var == "Mass" ){
                nomMassFSRUnfold = new TUnfoldDensityV17(hmcGenGen,
                                               TUnfold::kHistMapOutputHoriz,
                                               TUnfold::kRegModeNone, // fixed to use no regularisation temporary
                                               TUnfold::kEConstraintArea,
                                               TUnfoldDensityV17::kDensityModeBinWidth,
                                               binning_Gen,binning_Gen);
        }
        else{
            cout << "ISRUnfold::setNomTUnfoldDensity, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
        }
}

void ISRUnfold::setSysFSRTUnfoldDensity(TString var, TString filepath, TString sysName, int totSysN, int nth, TString phase_name, TString fsr_correction_name){

        TFile* filein = new TFile(filepath);

        TString systematic_postfix = sysName;

        if(totSysN == 2){
            if(nth == 0 ) systematic_postfix+="Up";
            if(nth == 1 ) systematic_postfix+="Down";
        }

        if(totSysN == 6 && sysName == "Scale"){
            if(nth == 0 ) systematic_postfix+="AUp";
            if(nth == 1 ) systematic_postfix+="ADown";
            if(nth == 2 ) systematic_postfix+="BUp";
            if(nth == 3 ) systematic_postfix+="BDown";
            if(nth == 4 ) systematic_postfix+="ABUp";
            if(nth == 5 ) systematic_postfix+="ABDown";
        }

        // set response matrix
        TH2* hmcGenGen;

        if( sysName == "AlphaS"){

            if(var == "Pt")
                hmcGenGen = (TH2*)filein->Get(phase_name + "/ptll_gen_post_fsr_" + fsr_correction_name + "_response_matrix/hmc" + var + "GenRec_" + systematic_postfix);
            else if(var == "Mass")
                hmcGenGen = (TH2*)filein->Get(phase_name + "/mll_gen_post_fsr_" + fsr_correction_name + "_response_matrix/hmc" + var + "GenRec_" + systematic_postfix);
            else{
                cout << "ISRUnfold::setNomTUnfoldDensity, only Pt and Mass available for var" << endl;
                exit (EXIT_FAILURE);
            }

        }
        else{
            if(var == "Pt")
                hmcGenGen = (TH2*)filein->Get(phase_name + "/ptll_gen_post_fsr_" + fsr_correction_name + "_response_matrix/hmc" + var + "GenRecnominal");
            else if(var == "Mass")
                hmcGenGen = (TH2*)filein->Get(phase_name + "/mll_gen_post_fsr_" + fsr_correction_name + "_response_matrix/hmc" + var + "GenRecnominal");
            else{
                cout << "ISRUnfold::setNomTUnfoldDensity, only Pt and Mass available for var" << endl;
                exit (EXIT_FAILURE);
            }
        }

        if(do_normalization){
            if(channel_name =="electron") hmcGenGen->Scale(0.959939);
            else hmcGenGen->Scale(0.953026);
        }

        // set binning definition
        TUnfoldBinning* binning_Gen = NULL;

        if( var == "Pt" ){

          TString Gen_Pt = "Gen_Pt";

          Gen_Pt = phase_name + "/ptll_gen_post_fsr_" + fsr_correction_name + "_response_matrix/" + Gen_Pt;

          binning_Gen = (TUnfoldBinning*)filein->Get(Gen_Pt);
        }
        else if( var == "Mass" ){

          TString Gen_Mass = "Gen_Mass";

          Gen_Mass = phase_name + "/mll_gen_post_fsr_" + fsr_correction_name + "_response_matrix/" + Gen_Mass;

          binning_Gen = (TUnfoldBinning*)filein->Get(Gen_Mass);

        }
        else{
            cout << "ISRUnfold::setNomTUnfoldDensity, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
        }

        if( var == "Pt" ){
                sysPtFSRUnfold[sysName].push_back( new TUnfoldDensityV17(hmcGenGen,
                                               TUnfold::kHistMapOutputHoriz,
                                               TUnfold::kRegModeNone, // fixed to use no regularisation temporary
                                               TUnfold::kEConstraintArea,
                                               TUnfoldDensityV17::kDensityModeBinWidth,
                                               binning_Gen,binning_Gen));
        }
        else if( var == "Mass" ){
                sysMassFSRUnfold[sysName].push_back( new TUnfoldDensityV17(hmcGenGen,
                                               TUnfold::kHistMapOutputHoriz,
                                               TUnfold::kRegModeNone, // fixed to use no regularisation temporary
                                               TUnfold::kEConstraintArea,
                                               TUnfoldDensityV17::kDensityModeBinWidth,
                                               binning_Gen,binning_Gen));
        }
        else{
            cout << "ISRUnfold::setNomTUnfoldDensity, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
        }
}

void ISRUnfold::setSysTUnfoldDensity(TString var, TString filepath, TString sysName, int totSysN, int nth, TString phase_name, TString fsr_correction_name){

        TFile* filein = new TFile(filepath);

        TString systematic_postfix = sysName;

        if(totSysN == 2){
            if(nth == 0 ) systematic_postfix+="Up";
            if(nth == 1 ) systematic_postfix+="Down";
        }

        if(totSysN == 6 && sysName =="Scale"){
            if(nth == 0 ) systematic_postfix+="AUp";
            if(nth == 1 ) systematic_postfix+="ADown";
            if(nth == 2 ) systematic_postfix+="BUp";
            if(nth == 3 ) systematic_postfix+="BDown";
            if(nth == 4 ) systematic_postfix+="ABUp";
            if(nth == 5 ) systematic_postfix+="ABDown";
        }

        // set migration matrix
        TH2* hmcGenRec = NULL;
	if(sysName=="Alt" || sysName=="unfoldBias" || sysName=="unfoldScan" || sysName=="Closure"){

            if(var == "Pt")   hmcGenRec = (TH2*)filein->Get(phase_name + "/ptll_rec_gen_" + fsr_correction_name + "_response_matrix/hmc" + var + "GenRecnominal");
            else if(var == "Mass") hmcGenRec = (TH2*)filein->Get(phase_name + "/mll_rec_gen_" + fsr_correction_name + "_response_matrix/hmc" + var + "GenRecnominal");
            else{
                cout << "ISRUnfold::setSysTUnfoldDensity, only Pt and Mass available for var" << endl;
                exit (EXIT_FAILURE);
            }
        }
        else{
            if(var == "Pt"){
                hmcGenRec = (TH2*)filein->Get(phase_name + "/ptll_rec_gen_" + fsr_correction_name + "_response_matrix/hmc" + var + "GenRec_" + systematic_postfix);
            }
            else if(var == "Mass"){
                hmcGenRec = (TH2*)filein->Get(phase_name + "/mll_rec_gen_" + fsr_correction_name + "_response_matrix/hmc" + var + "GenRec_" + systematic_postfix);
            }
            else{
                cout << "ISRUnfold::setSysTUnfoldDensity, only Pt and Mass available for var" << endl;
                exit (EXIT_FAILURE);
            }
        }
        
        if(do_normalization){
            if(channel_name =="electron") hmcGenRec->Scale(0.959939);
            else hmcGenRec->Scale(0.953026);
        }

        // set bin definition
        TUnfoldBinning* binning_Rec = NULL;
        TUnfoldBinning* binning_Gen = NULL;

        if( var == "Pt" ){

          TString Rec_Pt = "Rec_Pt";
          TString Gen_Pt = "Gen_Pt";

          Rec_Pt = phase_name + "/ptll_rec_gen_" + fsr_correction_name + "_response_matrix/" + Rec_Pt;
          Gen_Pt = phase_name + "/ptll_rec_gen_" + fsr_correction_name + "_response_matrix/" + Gen_Pt;

          binning_Rec = (TUnfoldBinning*)filein->Get(Rec_Pt);
          binning_Gen = (TUnfoldBinning*)filein->Get(Gen_Pt);
        }
        else if( var == "Mass" ){

          TString Rec_Mass = "Rec_Mass";
          TString Gen_Mass = "Gen_Mass";

          Rec_Mass = phase_name + "/mll_rec_gen_" + fsr_correction_name + "_response_matrix/" + Rec_Mass;
          Gen_Mass = phase_name + "/mll_rec_gen_" + fsr_correction_name + "_response_matrix/" + Gen_Mass;

          binning_Rec = (TUnfoldBinning*)filein->Get(Rec_Mass);
          binning_Gen = (TUnfoldBinning*)filein->Get(Gen_Mass);

        }
        else{
            cout << "ISRUnfold::setSysTUnfoldDensity, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
        }

	TUnfold::ERegMode mode=TUnfold::kRegModeNone;
	if( sysName =="unfoldScan" || sysName=="unfoldBias") mode = TUnfold::kRegModeCurvature;

        if( var == "Pt" ){
                sysPtUnfold[sysName].push_back( new TUnfoldDensityV17(hmcGenRec,
                                               TUnfold::kHistMapOutputHoriz,
                                               mode, // fixed to use no regularisation temporary
                                               TUnfold::kEConstraintArea,
                                               TUnfoldDensityV17::kDensityModeBinWidth,
                                               binning_Gen,binning_Rec));
        }

        else if( var == "Mass" ){
                sysMassUnfold[sysName].push_back( new TUnfoldDensityV17(hmcGenRec,
                                               TUnfold::kHistMapOutputHoriz,
                                               mode, // fixed to use no regularisation temporary
                                               TUnfold::kEConstraintArea,
                                               TUnfoldDensityV17::kDensityModeBinWidth,
                                               binning_Gen,binning_Rec));
        }
        else{
            cout << "ISRUnfold::setSysTUnfoldDensity, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
        }
}

void ISRUnfold::drawISRMatrixInfo(TString var, TString outpdf, bool detector_unfold, bool fsr_systematic){

    c1 = new TCanvas("c1","c1", 50, 50, 800, 700);
    gStyle->SetOptFit(0);
    gStyle->SetPalette(55);
    c1->cd();

    c1->SetBottomMargin(0.2);
    c1->SetRightMargin(0.2);
    c1->SetTopMargin(0.08);
    c1->SetTicks(1);
    c1->SetLogz();

    TH2 *histProb = NULL; 
    TH1 *histEfficiency = NULL;

    if(var=="Pt"){
        if(detector_unfold){
            histProb = nomPtUnfold->GetProbabilityMatrix("Migration prob. for pt mass bin",";p_T(gen);p_T(Rec)");
        }
        else{
            if(!fsr_systematic) histProb = nomPtFSRUnfold->GetProbabilityMatrix("Migration prob. for pt mass bin",";p_T(gen);p_T(Rec)");
            else histProb = sysPtFSRUnfold["QED_FSR"].at(0)->GetProbabilityMatrix("Migration prob. for pt mass bin",";p_T(gen);p_T(Rec)");
        }
    }
    else if(var=="Mass"){
        if(detector_unfold){
            histProb = nomMassUnfold->GetProbabilityMatrix("Migration prob. for mass bin",";mass(gen);mass(Rec)");
        }    
        else{
            if(!fsr_systematic) histProb = nomMassFSRUnfold->GetProbabilityMatrix("Migration prob. for mass bin",";mass(gen);mass(Rec)");
            else histProb = sysMassFSRUnfold["QED_FSR"].at(0)->GetProbabilityMatrix("Migration prob. for mass bin",";mass(gen);mass(Rec)");
        }    
    }
    else{
            cout << "ISRUnfold::drawISRMatrixInfo, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
    }

    histEfficiency=histProb->ProjectionX("histEfficiency");
    histProb->Draw("COLZ");

    if(var=="Pt"){
        if(detector_unfold) histProb->SetTitle("migration probabilities;p_{T} mass bin index (post FSR) ;p_{T} mass bin index (Rec)");
        else histProb->SetTitle("migration probabilities;p_{T} mass bin index (pre FSR) ;p_{T} mass bin index (post FSR)");
    }
    else if(var=="Mass"){
        if(detector_unfold) histProb->SetTitle("migration probabilities;mass bin index (post FSR) ;mass bin index (Rec)");
        else histProb->SetTitle("migration probabilities;mass bin index (pre FSR) ;mass bin index (post FSR)");
    }
    else{
            cout << "ISRUnfold::drawISRMatrixInfo, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
    }

    const TUnfoldBinningV17* xaxis_binning;
    const TVectorD* xaxis1_tvecd;
    const TVectorD* xaxis2_tvecd;
    int xaxis1_nbin;
    int xaxis2_nbin;

    const TUnfoldBinningV17* yaxis_binning = NULL;
    const TVectorD* yaxis1_tvecd = NULL;
    const TVectorD* yaxis2_tvecd = NULL;
    int yaxis1_nbin;
    int yaxis2_nbin;

    if(var=="Pt"){
        if(detector_unfold){
            xaxis_binning = nomPtUnfold->GetOutputBinning("Gen_Pt");
            xaxis1_tvecd = xaxis_binning->GetDistributionBinning(0);
            xaxis2_tvecd = xaxis_binning->GetDistributionBinning(1);
            xaxis1_nbin = xaxis1_tvecd->GetNrows() - 1;
            xaxis2_nbin = xaxis2_tvecd->GetNrows() - 1;

            yaxis_binning = nomPtUnfold->GetInputBinning("Rec_Pt");            
            yaxis1_tvecd = yaxis_binning->GetDistributionBinning(0);
            yaxis2_tvecd = yaxis_binning->GetDistributionBinning(1);
            yaxis1_nbin = yaxis1_tvecd->GetNrows() - 1;
            yaxis2_nbin = yaxis2_tvecd->GetNrows() - 1;
        }
        else{
            xaxis_binning = nomPtFSRUnfold->GetOutputBinning("Gen_Pt");
            xaxis1_tvecd = xaxis_binning->GetDistributionBinning(0);
            xaxis2_tvecd = xaxis_binning->GetDistributionBinning(1);
            xaxis1_nbin = xaxis1_tvecd->GetNrows() - 1;
            xaxis2_nbin = xaxis2_tvecd->GetNrows() - 1;
            
            yaxis_binning = nomPtFSRUnfold->GetInputBinning("Gen_Pt"); 
            yaxis1_tvecd = yaxis_binning->GetDistributionBinning(0);
            yaxis2_tvecd = yaxis_binning->GetDistributionBinning(1);
            yaxis1_nbin = yaxis1_tvecd->GetNrows() - 1;
            yaxis2_nbin = yaxis2_tvecd->GetNrows() - 1;
        }
    }
    else if(var=="Mass"){
        if(detector_unfold){
            xaxis_binning = nomMassUnfold->GetOutputBinning("Gen_Mass");
            xaxis1_tvecd = xaxis_binning->GetDistributionBinning(0);
            xaxis2_tvecd = xaxis_binning->GetDistributionBinning(1);
            xaxis1_nbin = xaxis1_tvecd->GetNrows() - 1; 
            xaxis2_nbin = xaxis2_tvecd->GetNrows() - 1; 

            yaxis_binning = nomMassUnfold->GetInputBinning("Rec_Mass");     
            yaxis1_tvecd = yaxis_binning->GetDistributionBinning(0);
            yaxis2_tvecd = yaxis_binning->GetDistributionBinning(1);
            yaxis1_nbin = yaxis1_tvecd->GetNrows() - 1; 
            yaxis2_nbin = yaxis2_tvecd->GetNrows() - 1; 
        }    
        else{
            xaxis_binning = nomMassFSRUnfold->GetOutputBinning("Gen_Mass");
            xaxis1_tvecd = xaxis_binning->GetDistributionBinning(0);
            xaxis2_tvecd = xaxis_binning->GetDistributionBinning(1);
            xaxis1_nbin = xaxis1_tvecd->GetNrows() - 1; 
            xaxis2_nbin = xaxis2_tvecd->GetNrows() - 1; 
     
            yaxis_binning = nomMassFSRUnfold->GetInputBinning("Gen_Mass"); 
            yaxis1_tvecd = yaxis_binning->GetDistributionBinning(0);
            yaxis2_tvecd = yaxis_binning->GetDistributionBinning(1);
            yaxis1_nbin = yaxis1_tvecd->GetNrows() - 1; 
            yaxis2_nbin = yaxis2_tvecd->GetNrows() - 1; 
        }    
    }    
    else{
            cout << "ISRUnfold::drawISRMatrixInfo, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
    }

    // find points to draw lines
    // count how many axis 1 bins -> then find n points
    // check axis 1 have under flow 

    TLine grid_;
    TLine grid_bin_boundary;
    //grid_.SetLineColor(kGray+2);
    grid_.SetLineColorAlpha(kGray+2, 0.35);;
    grid_.SetLineStyle(kSolid);
    grid_bin_boundary.SetLineColorAlpha(kRed, 1.);;
    grid_bin_boundary.SetLineStyle(kSolid);

    int boundarybin_x = xaxis1_nbin; 
    if(xaxis_binning->HasUnderflow(0)) boundarybin_x++;
    int boundarybin_y = yaxis1_nbin; 
    if(yaxis_binning->HasUnderflow(0)) boundarybin_y++;
    bool add_UO = true;
    int count_drawn_boundary = 0;

    for( int ii=0; ii<histProb->GetXaxis()->GetNbins(); ii++ )
    {
        Int_t i_bin = ii+1;
        Double_t binEdge = histProb->GetXaxis()->GetBinUpEdge(i_bin);
        grid_.DrawLine(binEdge, histProb->GetYaxis()->GetBinUpEdge(0), binEdge, histProb->GetYaxis()->GetBinUpEdge(histProb->GetYaxis()->GetNbins()) );

        if(boundarybin_x == i_bin && i_bin != histProb->GetXaxis()->GetNbins()){ 
            if(count_drawn_boundary == 0){
                grid_bin_boundary.DrawLine(histProb->GetXaxis()->GetBinUpEdge(i_bin-xaxis1_nbin), histProb->GetYaxis()->GetBinUpEdge(boundarybin_y - yaxis1_nbin), histProb->GetXaxis()->GetBinUpEdge(i_bin-xaxis1_nbin), histProb->GetYaxis()->GetBinUpEdge(boundarybin_y) );
            }

            if(add_UO){
                if(xaxis_binning->HasUnderflow(0)) boundarybin_x++;
                if(xaxis_binning->HasOverflow(0)) boundarybin_x++;
                add_UO = false;
            }
            else{
                boundarybin_x += xaxis1_nbin;
                add_UO = true;
            }
            grid_bin_boundary.DrawLine(binEdge, histProb->GetYaxis()->GetBinUpEdge(boundarybin_y - yaxis1_nbin), binEdge, histProb->GetYaxis()->GetBinUpEdge(boundarybin_y) );
            if(count_drawn_boundary%2 == 0){
                if(yaxis_binning->HasUnderflow(0)) boundarybin_y++;
                if(yaxis_binning->HasOverflow(0)) boundarybin_y++;
                boundarybin_y += yaxis1_nbin;
            }
            count_drawn_boundary++; 
        }
    }

    boundarybin_x = xaxis1_nbin; 
    boundarybin_y = yaxis1_nbin; 
    if(xaxis_binning->HasUnderflow(0)) boundarybin_x++;
    if(yaxis_binning->HasUnderflow(0)) boundarybin_y++;
    add_UO = true;
    count_drawn_boundary = 0;

    for( int ii=0; ii<histProb->GetYaxis()->GetNbins(); ii++ )
    {
        Int_t i_bin = ii+1;
        Double_t binEdge = histProb->GetYaxis()->GetBinUpEdge(i_bin);
        grid_.DrawLine(histProb->GetXaxis()->GetBinUpEdge(0), binEdge, histProb->GetXaxis()->GetBinUpEdge(histProb->GetXaxis()->GetNbins()),binEdge );

        if(boundarybin_y == i_bin && i_bin != histProb->GetYaxis()->GetNbins()){ 
            if(count_drawn_boundary == 0){
                grid_bin_boundary.DrawLine(histProb->GetXaxis()->GetBinUpEdge(boundarybin_x - xaxis1_nbin), histProb->GetYaxis()->GetBinUpEdge(i_bin-yaxis1_nbin), histProb->GetXaxis()->GetBinUpEdge(boundarybin_x), histProb->GetYaxis()->GetBinUpEdge(i_bin-yaxis1_nbin));
            }
            if(add_UO){
                if(yaxis_binning->HasUnderflow(0)) boundarybin_y++;
                if(yaxis_binning->HasOverflow(0)) boundarybin_y++;
                add_UO = false;
            }    
            else{
                boundarybin_y += yaxis1_nbin;
                add_UO = true;
            }    
            grid_bin_boundary.DrawLine(histProb->GetXaxis()->GetBinUpEdge(boundarybin_x - xaxis1_nbin), binEdge, histProb->GetXaxis()->GetBinUpEdge(boundarybin_x),binEdge);
            if(count_drawn_boundary%2 == 0){
                if(xaxis_binning->HasUnderflow(0)) boundarybin_x++;
                if(xaxis_binning->HasOverflow(0)) boundarybin_x++;
                boundarybin_x += xaxis1_nbin;
            }    
            count_drawn_boundary++; 
        } 
    }

    TLatex mcName;
    mcName.SetTextSize(0.035);
    if(!fsr_systematic){
        mcName.DrawLatexNDC(c1->GetLeftMargin(), 1-(c1->GetTopMargin()) + 0.2 * c1->GetTopMargin(), "Madgraph5 aMC@NLO + PYTHIA");
    }
    else mcName.DrawLatexNDC(c1->GetLeftMargin(), 1-(c1->GetTopMargin()) + 0.2 * c1->GetTopMargin(), "Powheg + PHOTOS");

    if(var=="Pt"){
        if(detector_unfold) c1->SaveAs(outpdf + "/detector_pt_matrix.pdf");
        else{
            if(!fsr_systematic) c1->SaveAs(outpdf + "/fsr_pt_matrix.pdf");
            else c1->SaveAs(outpdf + "/fsr_pt_photos_matrix.pdf");
        }
    }
    if(var=="Mass"){
        if(detector_unfold) c1->SaveAs(outpdf + "/detector_mass_matrix.pdf");
        else{
            if(!fsr_systematic) c1->SaveAs(outpdf + "/fsr_mass_matrix.pdf");
            else c1->SaveAs(outpdf + "/fsr_mass_photos_matrix.pdf");
        }
    }

    delete c1;

    c1 = new TCanvas("c1","c1", 50, 50, 800, 700);
    gStyle->SetOptFit(0);
    gStyle->SetPalette(55);
    c1->cd();

    c1->SetBottomMargin(0.2);
    c1->SetRightMargin(0.2);
    c1->SetTopMargin(0.08);
    c1->SetTicks(1);

    c1->cd();
    c1->cd()->SetBottomMargin(0.2);
    c1->cd()->SetRightMargin(0.2);
    c1->cd()->SetTopMargin(0.08);
    c1->cd()->SetTicks(1);
    c1->cd()->SetGridy();

    if(var=="Pt") histEfficiency->SetTitle("efficiency;p_{T} mass bin (gen) ;A #times #epsilon");
    if(var=="Mass") histEfficiency->SetTitle("efficiency;mass bin (gen) ;A #times #epsilon");
    histEfficiency->GetYaxis()->SetTitleOffset(1.5);
    histEfficiency->Draw();

    if(var=="Pt"){
        if(detector_unfold) c1->SaveAs(outpdf + "/detector_pt_accxeff.pdf");
        else{ 
            if(!fsr_systematic) c1->SaveAs(outpdf + "/fsr_pt_accxeff.pdf");
            else c1->SaveAs(outpdf + "/fsr_pt_photos_accxeff.pdf");
        }
    }
    if(var=="Mass"){
        if(detector_unfold) c1->SaveAs(outpdf + "/detector_mass_accxeff.pdf");
        else{
            if(!fsr_systematic) c1->SaveAs(outpdf + "/fsr_mass_accxeff.pdf");
            else c1->SaveAs(outpdf + "/fsr_mass_photos_accxeff.pdf");
        }
    }
    delete histProb;
    delete histEfficiency;
    delete c1;
}

void ISRUnfold::setFSRUnfoldInput(TString filepath, bool isSys, TString sysName, int nth, TString phase_name){

    TFile* filein = new TFile(filepath);
    TH1* hmass_input = NULL;
    TH1* hpt_input = NULL;

    //hpt_input= (TH1*)filein->Get(hist_dir + "/hist_ptll_post_fsr/histo_DYJetsnominal");
    //hpt_input->Add((TH1*)filein->Get(hist_dir + "/hist_ptll_post_fsr/histo_DYJets10to50nominal"));

    //hmass_input = (TH1*)filein->Get(hist_dir + "/hist_mll_post_fsr/histo_DYJetsnominal");
    //hmass_input->Add((TH1*)filein->Get(hist_dir + "/hist_mll_post_fsr/histo_DYJets10to50nominal"));

    //nomPtFSRUnfold->SetInput(hpt_input, 1.);
    //nomMassFSRUnfold->SetInput(hmass_input, 1.);

    if(!isSys){
        nomPtFSRUnfold->SetInput(nomPtUnfold->GetOutput("hpreFSR_pt", 0, 0, 0, false), 1.);
        nomMassFSRUnfold->SetInput(nomMassUnfold->GetOutput("hpreFSR_mass", 0, 0, 0, false), 1.);
    }
    else{
        if(sysName=="QED_FSR" || sysName=="QED_FSR_dRp1"){
            sysPtFSRUnfold[sysName].at(nth)  ->SetInput(nomPtUnfold->GetOutput("hpreFSR_pt", 0, 0, 0, false),   1.);
            sysMassFSRUnfold[sysName].at(nth)->SetInput(nomMassUnfold->GetOutput("hpreFSR_mass", 0, 0, 0, false),   1.);
        }
        else{
            sysPtFSRUnfold[sysName].at(nth)  ->SetInput(sysPtUnfold[sysName].at(nth)->GetOutput("hpreFSR_pt", 0, 0, 0, false),   1.);
            sysMassFSRUnfold[sysName].at(nth)->SetInput(sysMassUnfold[sysName].at(nth)->GetOutput("hpreFSR_mass", 0, 0, 0, false),   1.);
        }
    }

}

void ISRUnfold::doISRQEDFSRUnfold(bool doSys){

    nomPtFSRUnfold->DoUnfold(0);
    nomMassFSRUnfold->DoUnfold(0);

    if(doSys){
        std::map<TString, std::vector<TUnfoldDensityV17*>>::iterator it = sysPtFSRUnfold.begin();

        while(it != sysPtFSRUnfold.end()){
            int nSys = it->second.size();
            for(int i = 0; i < nSys; i++){
               it->second.at(i)->DoUnfold(0);
            }
            it++;
        }

        it = sysMassFSRUnfold.begin();
        while(it != sysMassFSRUnfold.end()){
                int nSys = it->second.size();
                for(int i = 0; i < nSys; i++){
                    it->second.at(i)->DoUnfold(0);
                }
                it++;
        }
    }
}

void ISRUnfold::setInput(TString var, TString filepath, bool isSys, TString sysName, int nth, double bias, TString phase_name){
	// No effects on the unfolded results respect to bias factor 
        
	TFile* filein = new TFile(filepath);
        TString nth_;
        nth_.Form("%d", nth);
        TH1* hRec = NULL;

        // nominal histograms
	if(!isSys){
            if(var == "Pt"){
                if(channel_name == "muon")     hRec = (TH1*)filein->Get(phase_name + "/hist_ptll/histo_DoubleMuonnominal");
                if(channel_name == "electron") hRec = (TH1*)filein->Get(phase_name + "/hist_ptll/histo_DoubleEGnominal");
                nomPtUnfold->SetInput(hRec,   bias);
            }
            else if(var == "Mass"){
                if(channel_name == "muon")     hRec = (TH1*)filein->Get(phase_name + "/hist_mll/histo_DoubleMuonnominal");
                if(channel_name == "electron") hRec = (TH1*)filein->Get(phase_name + "/hist_mll/histo_DoubleEGnominal");
                nomMassUnfold->SetInput(hRec, bias);
            }
            else{
                cout << "ISRUnfold::setInput, only Pt and Mass available for var" << endl;
                exit (EXIT_FAILURE);
            }
	}
        // systematic histograms
	else{
            // use DY MC histograms as unfolding input       
	    if(sysName=="Closure"){

                if(var == "Pt"){
                    hRec = (TH1*)filein->Get(phase_name + "/hist_ptll/histo_DYJetsnominal");
                    hRec->Add((TH1*)filein->Get(phase_name + "/hist_ptll/histo_DYJets10to50nominal"));
                    if(do_normalization){
                        if(channel_name=="electron") hRec->Scale(0.959939); 
                        if(channel_name=="muon") hRec->Scale(0.953026); 
                    }
                    sysPtUnfold[sysName].at(nth)  ->SetInput(hRec,   bias);
                }
                else if(var == "Mass"){
                    hRec = (TH1*)filein->Get(phase_name + "/hist_mll/histo_DYJetsnominal");
                    hRec->Add((TH1*)filein->Get(phase_name + "/hist_mll/histo_DYJets10to50nominal"));
                    if(do_normalization){
                        if(channel_name=="electron") hRec->Scale(0.959939);
                        if(channel_name=="muon") hRec->Scale(0.953026);
                    }
                    sysMassUnfold[sysName].at(nth)->SetInput(hRec,   bias);
                }
                else{
                    cout << "ISRUnfold::setInput, only Pt and Mass available for var" << endl;
                    exit (EXIT_FAILURE);

                }
            }
            // input histogram(data) for each systematic source
            else{
                // for systematic, using the same input histogram as nominal 
                if(var == "Pt"){
                    if(channel_name == "muon")     hRec = (TH1*)filein->Get(phase_name + "/hist_ptll/histo_DoubleMuonnominal");
                    if(channel_name == "electron") hRec = (TH1*)filein->Get(phase_name + "/hist_ptll/histo_DoubleEGnominal");
                    sysPtUnfold[sysName].at(nth)  ->SetInput(hRec,   bias);
                }
                else if(var == "Mass"){
                    if(channel_name == "muon")     hRec = (TH1*)filein->Get(phase_name + "/hist_mll/histo_DoubleMuonnominal");
                    if(channel_name == "electron") hRec = (TH1*)filein->Get(phase_name + "/hist_mll/histo_DoubleEGnominal");
                    sysMassUnfold[sysName].at(nth)->SetInput(hRec,   bias); 
                }
                else{
                    cout << "ISRUnfold::setInput, only Pt and Mass available for var" << endl;
                    exit (EXIT_FAILURE);
                }
            }
	}
	filein->Close();
	delete filein;
}

void ISRUnfold::subBkgs(TString var, TString filepath, TString bkgName, bool isSys, TString sysName, int totSysN, int nth, TString phase_name){

	TFile* filein = new TFile(filepath);
        TH1* hRec = NULL;

        double bkg_scale = 1.;
        if(do_normalization){
            if(channel_name=="electron") bkg_scale = 0.959939;
            if(channel_name=="muon") bkg_scale = 0.953026;
        }

        TString systematic_postfix = sysName;

        if(totSysN == 2){
            if(nth == 0 ) systematic_postfix+="Up";
            if(nth == 1 ) systematic_postfix+="Down";
        }

        if(totSysN == 6 && sysName == "Scale"){
            if(nth == 0 ) systematic_postfix+="AUp";
            if(nth == 1 ) systematic_postfix+="ADown";
            if(nth == 2 ) systematic_postfix+="BUp";
            if(nth == 3 ) systematic_postfix+="BDown";
            if(nth == 4 ) systematic_postfix+="ABUp";
            if(nth == 5 ) systematic_postfix+="ABDown";
        }

        // nominal histograms
	if(!isSys){
                if(var == "Pt"){
                    hRec = (TH1*)filein->Get(phase_name + "/hist_ptll/histo_" + bkgName + "nominal");
                }
                if(var == "Mass"){
                    hRec = (TH1*)filein->Get(phase_name + "/hist_mll/histo_" + bkgName + "nominal");
                }

        	if( var == "Pt" )   nomPtUnfold->  SubtractBackground(hRec, bkgName, bkg_scale);
        	if( var == "Mass" ) nomMassUnfold->SubtractBackground(hRec, bkgName, bkg_scale);
	}
        // systematic histograms
	else{	
            
            // systematics using nominal detector distributions
            if(sysName=="Alt" || sysName=="unfoldBias" || sysName=="unfoldScan"){
                if(var == "Pt"){
                    hRec = (TH1*)filein->Get(phase_name + "/hist_ptll/histo_" + bkgName + "nominal");
                }
                if(var == "Mass"){
                    hRec = (TH1*)filein->Get(phase_name + "/hist_mll/histo_" + bkgName + "nominal");
                }
            }
            else{
                if(var == "Pt"){
                    hRec = (TH1*)filein->Get(phase_name + "/hist_ptll/histo_" + bkgName + "_" + systematic_postfix);
                }
                if(var == "Mass"){
                    hRec = (TH1*)filein->Get(phase_name + "/hist_mll/histo_" + bkgName + "_" + systematic_postfix);
                }
            }

            if( var == "Pt" )   sysPtUnfold[sysName].at(nth)  ->SubtractBackground(hRec, bkgName, bkg_scale);
            if( var == "Mass" ) sysMassUnfold[sysName].at(nth)->SubtractBackground(hRec, bkgName, bkg_scale);
	}
	
	filein->Close();
	delete filein;
}

void ISRUnfold::drawLCurve(TString outpdf, TString var){

        gROOT->SetBatch();

	TGraph *lCurve_temp;
	Int_t iBest_temp = 0;
	if( var == "Pt"){
		lCurve_temp = lCurve;
		iBest_temp = iBest;
	}

	else if( var == "Mass"){
                lCurve_temp = lCurve_mass;
                iBest_temp = iBest_mass;
        }

	Double_t x[1],y[1];
	lCurve_temp->GetPoint(iBest_temp,x[0],y[0]);
	TGraph *bestLCurve=new TGraph(1,x,y);

        c1 = new TCanvas("c1","c1", 50, 50, 850, 700);
        c1->cd();
        gStyle->SetOptFit(0);

        c1->SetBottomMargin(0.2);
        c1->SetTopMargin(0.08);
        c1->SetTicks(1);

	lCurve_temp->Draw("AL");
  	bestLCurve->SetMarkerColor(kRed);
  	bestLCurve->Draw("*");
	c1->SaveAs(outpdf);	
	delete c1;
	

}

void ISRUnfold::drawRhoLog(TString outpdf, TString var){

        gROOT->SetBatch();

	TSpline *rhoLogTau_temp;
        Int_t iBest_temp = 0;
	Int_t nScan_temp = 0;
        if( var == "Pt"){
                rhoLogTau_temp = rhoLogTau;
                iBest_temp = iBest;
		nScan_temp = nScan;
        }

        else if( var == "Mass"){
                rhoLogTau_temp = rhoLogTau_mass;
                iBest_temp = iBest_mass;
		nScan_temp = nScan_mass;
        }


        Double_t t[1],rho[1];
        rhoLogTau_temp->GetKnot(iBest_temp,t[0],rho[0]);
	TGraph *bestRhoLogTau=new TGraph(1,t,rho);

  	Double_t *tAll=new Double_t[nScan_temp],*rhoAll=new Double_t[nScan_temp];
  	for(Int_t i=0;i<nScan_temp;i++) {
  	   rhoLogTau_temp->GetKnot(i,tAll[i],rhoAll[i]);
	   cout << i << " tAll[i]: " << tAll[i] << " rhoAll[i]: " << rhoAll[i] << endl;
  	}
  	TGraph *knots=new TGraph(nScan_temp,tAll,rhoAll);

        c1 = new TCanvas("c1","c1", 50, 50, 850, 700);
        c1->cd();
        gStyle->SetOptFit(0);

        c1->SetBottomMargin(0.2);
        c1->SetTopMargin(0.08);
        c1->SetTicks(1);

  	//rhoLogTau_temp->Draw();
  	knots->Draw("AP");
  	bestRhoLogTau->SetMarkerColor(kRed);
  	bestRhoLogTau->Draw("*");

        c1->SaveAs(outpdf);
        delete c1;


}

void ISRUnfold::doISRUnfold(bool doSys){

	double tauMin=1.e-4;
	double tauMax=1.e-3;

  	nScan=30;
  	rhoLogTau=0;
  	lCurve=0;
	iBest = 0;

        nScan_mass=30;
        rhoLogTau_mass=0;
        lCurve_mass=0;
        iBest_mass = 0;

	nomPtUnfold->DoUnfold(0);
	nomMassUnfold->DoUnfold(0);


        if(doSys){
	    std::map<TString, std::vector<TUnfoldDensityV17*>>::iterator it = sysPtUnfold.begin();

	    while(it != sysPtUnfold.end()){
	    	int nSys = it->second.size();
	    	for(int i = 0; i < nSys; i++){
	    		if((it->first)=="unfoldScan" || (it->first)=="unfoldBias" ){
	    			//it->second.at(i)->ScanLcurve(50,tauMin,tauMax,0);
  	    			iBest=it->second.at(i)->ScanTau(nScan,0.,0.,&rhoLogTau,
  	    			                           TUnfoldDensity::kEScanTauRhoAvgSys,
  	    			                           0,0,
  	    			                           &lCurve);
	    		}
	    		else it->second.at(i)->DoUnfold(0);

	    	}
	    	it++;
	    }

	    it = sysMassUnfold.begin();
            while(it != sysMassUnfold.end()){
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
double ISRUnfold::DoFit(TString var, int nthMassBin, bool isFSRUnfold){

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
            hpt_temp_mc   = nomPtUnfold->GetBias("histMCTruth_pt_temp_chi",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
            hrho          = nomPtUnfold->GetRhoIJtotal("histRho_chi", 0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE); 
        }
        else{
            hpt_temp_data = nomPtFSRUnfold->GetOutput("hunfolded_pt_temp_chi",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
            hpt_temp_mc   = nomPtFSRUnfold->GetBias("histMCTruth_pt_temp_chi",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
            hrho          = nomPtFSRUnfold->GetRhoIJtotal("histRho_chi", 0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE); 
        } 
    }
    else{
        if(!isFSRUnfold){
            hpt_temp_data = nomMassUnfold->GetOutput("hunfolded_pt_temp_chi",0,0,"mass[UO];pt[UOC0]",kTRUE);
            hpt_temp_mc   = nomMassUnfold->GetBias("histMCTruth_pt_temp_chi",0,0,"mass[UO];pt[UOC0]",kTRUE);
            hrho          = nomMassUnfold->GetRhoIJtotal("histRho_chi", 0,0,"mass[UO];pt[UOC0]",kTRUE);
        }
        else{
            hpt_temp_data = nomMassFSRUnfold->GetOutput("hunfolded_pt_temp_chi",0,0,"mass[UO];pt[UOC0]",kTRUE);
            hpt_temp_mc   = nomMassFSRUnfold->GetBias("histMCTruth_pt_temp_chi",0,0,"mass[UO];pt[UOC0]",kTRUE);
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

void ISRUnfold::setMCPreFSRMeanValues(TString filepath){

    TFile* filein = new TFile(filepath);
    TH2* hmcPtGenGen;
    hmcPtGenGen = (TH2*)filein->Get("full_phase/ptll_gen_post_fsr_pre_MC_fsr_response_matrix/hmcPtGenRecnominal");
    TH2* hmcMassGenGen;
    hmcMassGenGen = (TH2*)filein->Get("full_phase/mll_gen_post_fsr_pre_fsr_response_matrix/hmcMassGenRecnominal");

    TUnfoldBinning* binningPt_Gen = (TUnfoldBinning*)filein->Get("full_phase/ptll_gen_post_fsr_pre_MC_fsr_response_matrix/Gen_Pt");
    TUnfoldBinning* binningMass_Gen = (TUnfoldBinning*)filein->Get("full_phase/mll_gen_post_fsr_pre_fsr_response_matrix/Gen_Mass");

    TUnfoldDensityV17* tempPtFSRUnfold = new TUnfoldDensityV17(hmcPtGenGen,
                                               TUnfold::kHistMapOutputHoriz,
                                               TUnfold::kRegModeNone, // fixed to use no regularisation temporary
                                               TUnfold::kEConstraintArea,
                                               TUnfoldDensityV17::kDensityModeBinWidth,
                                               binningPt_Gen,binningPt_Gen);
    
    TUnfoldDensityV17* tempMassFSRUnfold = new TUnfoldDensityV17(hmcMassGenGen,
                                               TUnfold::kHistMapOutputHoriz,
                                               TUnfold::kRegModeNone, // fixed to use no regularisation temporary
                                               TUnfold::kEConstraintArea,
                                               TUnfoldDensityV17::kDensityModeBinWidth,
                                               binningMass_Gen,binningMass_Gen);

    TH1* hpt_input = NULL;

    hpt_input= (TH1*)filein->Get("full_phase/hist_ptll_post_fsr_MC/histo_DYJetsnominal");
    hpt_input->Add((TH1*)filein->Get("full_phase/hist_ptll_post_fsr_MC/histo_DYJets10to50nominal"));

    // just set input to get bias distribution
    tempPtFSRUnfold->SetInput(hpt_input, 1.);
    tempMassFSRUnfold->SetInput(nomMassUnfold->GetOutput("hpreFSR_mass", 0, 0, 0, false), 1.);

    int nMassBin = -1;

    const TUnfoldBinningV17* temp_binning = tempPtFSRUnfold->GetOutputBinning("Gen_Pt");
    const TVectorD* temp_tvecd = temp_binning->GetDistributionBinning(1); // get mass distribution
    nMassBin = temp_tvecd->GetNrows() - 1;
    const Double_t* massBins = temp_tvecd->GetMatrixArray();

    // set mean mass and error
    TH1 *histMCTruth_mass= tempMassFSRUnfold->GetBias("histMCTruth_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);
    for(int ibin = 0; ibin < nMassBin; ibin++){

        TString ibinMass;
        ibinMass.Form("%d", ibin);

        histMCTruth_mass->GetXaxis()->SetRange(histMCTruth_mass->GetXaxis()->FindBin(massBins[ibin]+0.01),histMCTruth_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
        meanMass_mc_pre_fsr_.   push_back(histMCTruth_mass->GetMean());
        meanMassErr_mc_pre_fsr_.push_back(histMCTruth_mass->GetMeanError());

        TH1* hpt_temp_mc = tempPtFSRUnfold->GetBias("histMCTruth_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);

        meanPt_mc_pre_fsr_.   push_back(hpt_temp_mc->GetMean());
        meanPtErr_mc_pre_fsr_.push_back(hpt_temp_mc->GetMeanError());

        cout << ibin << " pt: " << hpt_temp_mc->GetMean() << " mass: " << histMCTruth_mass->GetMean() << endl;
        delete hpt_temp_mc;
    }

    delete tempPtFSRUnfold, tempMassFSRUnfold;
}

void ISRUnfold::setMeanMass(bool doSys, bool altMC, bool detector_unfold){

        int nMassBin = -1;

        const TUnfoldBinningV17* temp_binning = nomPtUnfold->GetOutputBinning("Gen_Pt");
	const TVectorD* temp_tvecd = temp_binning->GetDistributionBinning(1);
        nMassBin = temp_tvecd->GetNrows() - 1; 
	const Double_t* massBins = temp_tvecd->GetMatrixArray();

        // set nominal detector & QED FSR unfolded results
        TH1* hunfolded_mass =  nomMassUnfold->GetOutput("hunfolded_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);
        TH1 *histMCTruth_mass= nomMassUnfold->GetBias("histMCTruth_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);

        TH1* h_fsr_unfolded_mass =  nomMassFSRUnfold->GetOutput("h_fsr_unfolded_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);
        TH1 *histMCTruth_pre_fsr_mass= nomMassFSRUnfold->GetBias("histMCTruth_pre_fsr_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);

        TH1* h_fsr_dRp1_unfolded_mass =  sysMassFSRUnfold["QED_FSR_dRp1"].at(0)->GetOutput("h_fsr_dRp1_unfolded_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);

        // loop over mass bins
        for(int ibin = 0; ibin < nMassBin; ibin++){
                hunfolded_mass->GetXaxis()->  SetRange(hunfolded_mass->GetXaxis()->  FindBin(massBins[ibin]+0.01),hunfolded_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
                histMCTruth_mass->GetXaxis()->SetRange(histMCTruth_mass->GetXaxis()->FindBin(massBins[ibin]+0.01),histMCTruth_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));

                h_fsr_unfolded_mass->GetXaxis()->  SetRange(hunfolded_mass->GetXaxis()->  FindBin(massBins[ibin]+0.01),hunfolded_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
                histMCTruth_pre_fsr_mass->GetXaxis()->SetRange(histMCTruth_mass->GetXaxis()->FindBin(massBins[ibin]+0.01),histMCTruth_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));

                h_fsr_dRp1_unfolded_mass->GetXaxis()->  SetRange(hunfolded_mass->GetXaxis()->  FindBin(massBins[ibin]+0.01),hunfolded_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));

                meanMass_data.   push_back(hunfolded_mass->GetMean());
                meanMassStatErr_data.push_back(hunfolded_mass->GetMeanError());

                meanMass_data_pre_fsr.   push_back(h_fsr_unfolded_mass->GetMean());
                meanMassStatErr_data_pre_fsr.push_back(h_fsr_unfolded_mass->GetMeanError());

                meanMass_data_pre_fsr_dRp1.   push_back(h_fsr_dRp1_unfolded_mass->GetMean());
                meanMassStatErr_data_pre_fsr_dRp1.push_back(h_fsr_dRp1_unfolded_mass->GetMeanError());

                meanMass_mc.   push_back(histMCTruth_mass->GetMean());
                meanMassErr_mc.push_back(histMCTruth_mass->GetMeanError());

                if(altMC){
                    TH1 *histMCTruth_massAlt= sysMassUnfold["Alt"].at(0)->GetBias("histMCTruth_massAlt",0,0,"mass[UO];pt[UOC0]",kTRUE);
                    histMCTruth_massAlt->GetXaxis()->SetRange(histMCTruth_mass->GetXaxis()->FindBin(massBins[ibin]+0.01),histMCTruth_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
                    meanMass_mcAlt.   push_back(histMCTruth_massAlt->GetMean());
                    meanMassErr_mcAlt.push_back(histMCTruth_massAlt->GetMeanError());

	            delete histMCTruth_massAlt;
                }

                // save systematic mean mass values 
                if(doSys){
                    // detector unfolding systematics 
	            std::map<TString, std::vector<TUnfoldDensityV17*>>::iterator it = sysMassUnfold.begin();
                    std::map<TString, std::vector<Double_t>> temp_map; // temp map to save systematic results for a mass bin

           	    while(it != sysMassUnfold.end()){
           	            int nSys = it->second.size();
           	            TH1* hdatasys_temp;
           	            for(int i = 0; i < nSys; i++){
           	                 hdatasys_temp = sysMassUnfold[it->first].at(i)->GetOutput("hunfolded_mass_systemp",0,0,"mass[UO];pt[UOC0]",kTRUE);
		                 hdatasys_temp->GetXaxis()->SetRange(hdatasys_temp->GetXaxis()->FindBin(massBins[ibin]+0.01), hdatasys_temp->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
           	                 temp_map[it->first].push_back(hdatasys_temp->GetMean());

           	                 delete hdatasys_temp;
           	            }
           	            it++;
           	    }
           	    meanMass_sysdata.push_back(temp_map);
           	    temp_map.clear();

                    // QED FSR unfolding systematics
                    it = sysMassFSRUnfold.begin();
                    while(it != sysMassFSRUnfold.end()){
                            int nSys = it->second.size();
                            TH1* hdatasys_temp;
                            for(int i = 0; i < nSys; i++){
                                 hdatasys_temp = sysMassFSRUnfold[it->first].at(i)->GetOutput("hunfolded_mass_systemp",0,0,"mass[UO];pt[UOC0]",kTRUE);
                                 hdatasys_temp->GetXaxis()->SetRange(hdatasys_temp->GetXaxis()->FindBin(massBins[ibin]+0.01), hdatasys_temp->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
                                 temp_map[it->first].push_back(hdatasys_temp->GetMean());

                                 delete hdatasys_temp;
                            }
                            it++;
                    }
                    meanMass_sysdata_pre_fsr.push_back(temp_map);
                    temp_map.clear();
                }
        }

	delete hunfolded_mass;
	delete histMCTruth_mass;
        delete h_fsr_unfolded_mass;
        delete histMCTruth_pre_fsr_mass; 

        // calculate systematic uncertainty
        if(doSys){
            for(int i = 0; i < nMassBin; i++){
            
                // detector unfolding
               std::map<TString, Double_t> temp_map_; 
               std::map<TString, Int_t> temp_map_sysname_index;

               std::map<TString, std::vector<Double_t>>::iterator it = meanMass_sysdata.at(i).begin();
               while(it != meanMass_sysdata.at(i).end()){
                    int size_ = it->second.size();
                    if( (it->first) == "Closure" ){ 
                        it++; continue;
                    }

                    TH1F *hpdfsys = NULL;
                    if((it->first)=="PDFerror") hpdfsys = new TH1F("pdfsys", "pdfsys", 100, meanMass_data.at(i)-0.2, meanMass_data.at(i)+0.2); // temp histogram to contain PDF variations

                    Double_t err = -999.; // 
                    // use maximum variation as systematic
                    for(int j = 0; j < size_; j++){
                            //if( (i==5 || i==7) && (it->first)=="Scale") continue;

                            if((it->first)=="PDFerror"){
                                    hpdfsys->Fill(it->second.at(j));
                            }

                            Double_t temp_err =  fabs(meanMass_data.at(i) - it->second.at(j));
                            if( temp_err > err){
                                   err = temp_err;
                            }
                            //cout << i << " th mass bin, " << it->first << j << " th sys value: " << it->second.at(j) << endl;
                    }
                    if((it->first)=="PDFerror"){
                            err = hpdfsys->GetRMS();
                            delete hpdfsys;
                    }

                    temp_map_[it->first] = err;
                    it++;
               }
               meanMassErr_sysdata.push_back(temp_map_);
               temp_map_.clear();

               it = meanMass_sysdata_pre_fsr.at(i).begin();
               while(it != meanMass_sysdata_pre_fsr.at(i).end()){
                    int size_ = it->second.size();
                    if( (it->first) == "Closure" || (it->first) == "QED_FSR_dRp1"){
                        it++; continue;
                    }

                    TH1F *hpdfsys = NULL;
                    if((it->first)=="PDFerror") 
                        hpdfsys = new TH1F("pdfsys", "pdfsys", 100, meanMass_data_pre_fsr.at(i)-0.2, meanMass_data_pre_fsr.at(i)+0.2); // temp histogram to contain PDF variations

                    Double_t err = -999.; // 
                    Int_t index_to_save = -1;
                    // use maximum variation as systematic
                    for(int j = 0; j < size_; j++){
                            //if( (i==5 || i==7) && (it->first)=="Scale") continue;

                            if((it->first)=="PDFerror"){
                                    hpdfsys->Fill(it->second.at(j));
                            }

                            Double_t temp_err =  0.;
                            if((it->first) == "QED_FSR"){
                                temp_err =  fabs(it->second.at(0) - it->second.at(1));
                            }
                            else{
                                temp_err =  fabs(meanMass_data_pre_fsr.at(i) - it->second.at(j));
                            }

                            if( temp_err > err){
                                err = temp_err;
                                index_to_save = j;
                            }
                            //cout << i << " th mass bin, " << it->first << j << " th sys value: " << it->second.at(j) << endl;
                    }
                    if((it->first)=="PDFerror"){
                            err = hpdfsys->GetRMS();
                            delete hpdfsys;
                    }

                    temp_map_[it->first] = err;
                    if((it->first)!="PDFerror" || (it->first)!="QED_FSR"){
                        temp_map_sysname_index[it->first] = index_to_save;
                    }
                    it++;
               }
               meanMassErr_sysdata_pre_fsr.push_back(temp_map_);
               meanMassErrIdx_sysdata_pre_fsr.push_back(temp_map_sysname_index);
               temp_map_.clear();
               temp_map_sysname_index.clear();

            }// loop for mass bins
        }
        
        for(int i = 0; i < nMassBin; i++){
            Double_t totalSys = 0.;
            Double_t totalSys_pre_fsr = 0.;
           
            cout << i << " th mass bin... " << endl; 
            if(doSys){
                std::map<TString, Double_t>::iterator it = meanMassErr_sysdata.at(i).begin();
                while(it != meanMassErr_sysdata.at(i).end()){
                    
                     //cout << i << " th mass bin, mass" << it->first << " " << it->second << endl;
                     totalSys += pow(it->second, 2);
                     it++;
                }
                //cout << i << " th mass bin, total mass systematic uncertainty: " << sqrt(totalSys) << " statistical error: " << meanMassStatErr_data.at(i) << endl;

                it = meanMassErr_sysdata_pre_fsr.at(i).begin();
                while(it != meanMassErr_sysdata_pre_fsr.at(i).end()){
               
                     cout << "mass systematic: " << it->first << " " << it->second/ meanMass_data_pre_fsr.at(i) * 100.<< endl; 
                     //cout << i << " th mass bin, mass" << it->first << " " << it->second << endl;
                     totalSys_pre_fsr += pow(it->second, 2);
                     it++;
                }

            }
            meanMassSysErr_data.push_back(sqrt(totalSys));
	    meanMassTotErr_data.push_back(sqrt(totalSys + pow(meanMassStatErr_data.at(i),2)));

            meanMassSysErr_data_pre_fsr.push_back(sqrt(totalSys_pre_fsr));
	    meanMassTotErr_data_pre_fsr.push_back(sqrt(totalSys_pre_fsr + pow(meanMassStatErr_data_pre_fsr.at(i),2)));
        }
}


// set mean pt from mass and DY mc
void ISRUnfold::setMeanPt(bool doSys, bool altMC, bool detector_unfold){
      
    int nMassBin = -1;
    
    const TUnfoldBinningV17* temp_binning = nomPtUnfold->GetOutputBinning("Gen_Pt");
    const TVectorD* temp_tvecd = temp_binning->GetDistributionBinning(1);
    nMassBin = temp_tvecd->GetNrows() - 1;

    // save mean pt for each systematic variation
    for(int i = 0; i < nMassBin; i++){
         
        TString ibinMass;
        ibinMass.Form("%d", i);

        TH1* hpt_temp_data;
        TH1* hpt_temp_mc;
        TH1* hpt_temp_mcAlt;

        hpt_temp_data = nomPtUnfold->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
        hpt_temp_mc   = nomPtUnfold->GetBias("histMCTruth_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);

        TH1* h_pre_fsr_pt_temp_data = nomPtFSRUnfold->GetOutput("h_fsr_unfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
        TH1* h_pre_fsr_pt_temp_mc   = nomPtFSRUnfold->GetBias("histMCTruth_pt_fsr_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);

        TH1* h_pre_fsr_dRp1_pt_temp_data = sysPtFSRUnfold["QED_FSR_dRp1"].at(0)->GetOutput("h_fsr_dRp1_unfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);

        meanPt_data.   push_back(hpt_temp_data->GetMean());
        meanPtStatErr_data.push_back(hpt_temp_data->GetMeanError());

        meanPt_data_pre_fsr.   push_back(h_pre_fsr_pt_temp_data->GetMean());
        meanPtStatErr_data_pre_fsr.push_back(h_pre_fsr_pt_temp_data->GetMeanError());

        meanPt_data_pre_fsr_dRp1.   push_back(h_pre_fsr_dRp1_pt_temp_data->GetMean());
        meanPtStatErr_data_pre_fsr_dRp1.push_back(h_pre_fsr_dRp1_pt_temp_data->GetMeanError());

        meanPt_mc.   push_back(hpt_temp_mc->GetMean());
        meanPtErr_mc.push_back(hpt_temp_mc->GetMeanError());

        if(altMC){
            hpt_temp_mcAlt   = sysPtUnfold["Alt"].at(0)->GetBias("histMCTruth_pt_tempAlt",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
            meanPt_mcAlt.   push_back(hpt_temp_mcAlt->GetMean());
            meanPtErr_mcAlt.push_back(hpt_temp_mcAlt->GetMeanError());
            delete hpt_temp_mcAlt;
        }
         
        if(doSys){
            // detector unfolding
            std::map<TString, std::vector<TUnfoldDensityV17*>>::iterator it = sysPtUnfold.begin();
            std::map<TString, std::vector<Double_t>> temp_map; // temp map to save systematic results for a mass bin
                                                               // format: systematic name, mean pt
            while(it != sysPtUnfold.end()){
                    int nSys = it->second.size();
                    TH1* hdatasys_temp;
                    for(int isys = 0; isys < nSys; isys++){
                 	hdatasys_temp = sysPtUnfold[it->first].at(isys)->GetOutput("hunfolded_pt_systemp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
                 	temp_map[it->first].push_back(hdatasys_temp->GetMean());
    
                 	delete hdatasys_temp;
                    }
                    it++;
            }
            meanPt_sysdata.push_back(temp_map);
            temp_map.clear();

            // QED FSR unfolding
            it = sysPtFSRUnfold.begin();

            while(it != sysPtFSRUnfold.end()){
                    int nSys = it->second.size();
                    TH1* hdatasys_temp;
                    for(int isys = 0; isys < nSys; isys++){
                        hdatasys_temp = sysPtFSRUnfold[it->first].at(isys)->GetOutput("hunfolded_pt_systemp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
                        temp_map[it->first].push_back(hdatasys_temp->GetMean());
                        delete hdatasys_temp;
                    }
                    it++;
            }   
            meanPt_sysdata_pre_fsr.push_back(temp_map);
            temp_map.clear();
        }
        
        delete hpt_temp_data;
        delete hpt_temp_mc;
        delete h_pre_fsr_pt_temp_data;
        delete h_pre_fsr_pt_temp_mc;
        delete h_pre_fsr_dRp1_pt_temp_data;
    }

    // calculate systematic uncertainty for each systematic variation
    if(doSys){
	for(int i = 0; i < nMassBin; i++){
           
            // detector unfolding 
            std::map<TString, Double_t> temp_map_;
            std::map<TString, Int_t> temp_map_sysname_index;

            std::map<TString, std::vector<Double_t>>::iterator it = meanPt_sysdata.at(i).begin();
            while(it != meanPt_sysdata.at(i).end()){
	    	int size_ = it->second.size(); // size of systematic variations
	    	if((it->first)=="Closure"){
                    it++; continue;
                }
                	
	    	TH1F *hpdfsys = NULL;
	    	if((it->first)=="PDFerror")
                    hpdfsys = new TH1F("pdfsys", "pdfsys", 100, meanPt_data.at(i)-0.2, meanPt_data.at(i)+0.2); // temp histogram to contain PDF variations

	    	Double_t err = -999.; // 
	    	for(int j = 0; j < size_; j++){
	    	    //if( (i==5 || i==7) && (it->first)=="Scale") continue;

	    	    if((it->first)=="PDFerror"){
	     		hpdfsys->Fill(it->second.at(j));
	    	    }

                    Double_t temp_err =  fabs(meanPt_data.at(i) - it->second.at(j));
                    if( temp_err > err){
                        err = temp_err;
                    }
	    	    //cout << i << " th mass bin, " << it->first << j << " th sys value: " << it->second.at(j) << endl; 
	    	}// loop for systematic variations
	    	if((it->first)=="PDFerror"){
	    	    err = hpdfsys->GetRMS();
	    	    delete hpdfsys;
	    	}

	    	temp_map_[it->first] = err;
	    	it++;
            }// loop for systematic sources
	    meanPtErr_sysdata.push_back(temp_map_);
	    temp_map_.clear();

            // QED FSR unfolding
            it = meanPt_sysdata_pre_fsr.at(i).begin();
            while(it != meanPt_sysdata_pre_fsr.at(i).end()){
                int size_ = it->second.size(); // size of systematic variations
                if((it->first)=="Closure" || (it->first)=="QED_FSR_dRp1"){
                    it++; continue;
                }   
                
                TH1F *hpdfsys = NULL;
                if((it->first)=="PDFerror")
                    hpdfsys = new TH1F("pdfsys", "pdfsys", 100, meanPt_data_pre_fsr.at(i)-0.2, meanPt_data_pre_fsr.at(i)+0.2); // temp histogram to contain PDF variations

                Double_t err = -999.; // 
                Int_t index_to_save = -1;
                for(int j = 0; j < size_; j++){

                    if((it->first)=="PDFerror"){
                        hpdfsys->Fill(it->second.at(j));
                    }   

                    Double_t temp_err =  0.;
                    if((it->first) == "QED_FSR"){
                        temp_err =  fabs(it->second.at(0) - it->second.at(1));
                    }
                    else{
                        temp_err =  fabs(meanPt_data_pre_fsr.at(i) - it->second.at(j));
                    }

                    if( temp_err > err){
                        err = temp_err;
                        index_to_save = j;
                    }   
                    //cout << i << " th mass bin, " << it->first << j << " th sys value: " << it->second.at(j) << endl; 
                }// loop for systematic variations

                if((it->first)=="PDFerror"){
                    err = hpdfsys->GetRMS(); // only for 2016 DY MC
                    delete hpdfsys;
                }
                //if((it->first) == "QED_FSR") cout << "QED FSR error (pt): " << err << endl;
                temp_map_[it->first] = err;

                if((it->first)!="PDFerror" || (it->first)!="QED_FSR"){
                    temp_map_sysname_index[it->first] = index_to_save;
                }
                it++;
            }// loop for systematic sources
            meanPtErr_sysdata_pre_fsr.push_back(temp_map_);
            meanPtErrIdx_sysdata_pre_fsr.push_back(temp_map_sysname_index);
            temp_map_.clear();
            temp_map_sysname_index.clear();

        }// loop for mass binss
    }

    for(int i = 0; i < nMassBin; i++){
        Double_t totalSys = 0.;
        Double_t totalSys_pre_fsr = 0.;

        cout << i << " th mass bin... " << endl;
        if(doSys){
            std::map<TString, Double_t>::iterator it = meanPtErr_sysdata.at(i).begin();
            while(it != meanPtErr_sysdata.at(i).end()){
                //cout << i << " th mass bin, pt" << it->first << " " << it->second << endl;
                totalSys += pow(it->second, 2);	
                it++;
            }
     
            //cout << i << " th mass bin, total pt systematic uncertainty: " << sqrt(totalSys) << " statistical error: " << meanPtStatErr_data.at(i) << endl;
            it = meanPtErr_sysdata_pre_fsr.at(i).begin();
            while(it != meanPtErr_sysdata_pre_fsr.at(i).end()){
                cout << "systematic: " << it->first << " " << it->second/ meanPt_data_pre_fsr.at(i) * 100.<< endl;
                totalSys_pre_fsr += pow(it->second, 2);	
                it++;
            }
        }
        meanPtSysErr_data.push_back(sqrt(totalSys));
        meanPtTotErr_data.push_back(sqrt(totalSys + pow(meanPtStatErr_data.at(i),2)));

        meanPtSysErr_data_pre_fsr.push_back(sqrt(totalSys_pre_fsr));
        meanPtTotErr_data_pre_fsr.push_back(sqrt(totalSys_pre_fsr + pow(meanPtStatErr_data_pre_fsr.at(i),2)));
    }// loop for mass bins
}

void ISRUnfold::drawISRresult(TString outpdf, bool altMC, bool doFit){

        gROOT->SetBatch();

        int marker_ = 20;
        if(channel_name=="muon") marker_ = 20;

        setTDRStyle();
        writeExtraText = true;       // if extra text
        extraText  = "work in progress";

        c1 = new TCanvas("c1","c1", 50, 50, 900, 700);
        c1->cd();
        gStyle->SetOptFit(0);

        c1->SetBottomMargin(0.2);
        c1->SetTopMargin(0.08);
        //c1->SetTicks();
        c1->SetLogx();
        //c1->SetGridy();
        //c1->SetGridx();

        TGraphErrors *grUnfolded = new TGraphErrors(5, &meanMass_data[0], &meanPt_data[0], &meanMassTotErr_data[0], &meanPtTotErr_data[0]);
        grUnfolded->SetLineColor(kGray+1);
        grUnfolded->SetMarkerColor(kGray+1);
        grUnfolded->SetMarkerStyle(marker_);
        grUnfolded->SetMarkerSize(0.9);
        grUnfolded->SetLineStyle(1);
        grUnfolded->Draw("apZ");
        grUnfolded->GetYaxis()->SetRangeUser(10.,30.);
        grUnfolded->GetXaxis()->SetLimits(30.,500.);
        grUnfolded->GetXaxis()->SetMoreLogLabels(true);
        grUnfolded->GetYaxis()->SetTitle("Average p_{T} (GeV)");
        grUnfolded->GetXaxis()->SetTitle("Average Mass (GeV)");

        TGraphErrors *grMC = new TGraphErrors(5, &meanMass_mc[0], &meanPt_mc[0], &meanMassErr_mc[0], &meanPtErr_mc[0]);
        grMC->SetLineColor(kRed);
        grMC->SetMarkerColor(kRed);
        grMC->SetMarkerStyle(marker_);
        grMC->SetMarkerSize(.9);
        grMC->SetLineStyle(1);
        grMC->SetLineColor(kRed);
        grMC->Draw("pZ same");

        TGraphErrors *grMC_preFSR = new TGraphErrors(9, &meanMass_mc_pre_fsr_[0], &meanPt_mc_pre_fsr_[0], &meanMassErr_mc_pre_fsr_[0], &meanPtErr_mc_pre_fsr_[0]);
        grMC_preFSR->SetLineColor(kRed);
        grMC_preFSR->SetMarkerColor(kRed);
        grMC_preFSR->SetMarkerStyle(marker_);
        grMC_preFSR->SetMarkerSize(.2);
        grMC_preFSR->SetLineStyle(1);
        grMC_preFSR->SetLineColor(kRed);
        grMC_preFSR->Draw("pC same");

        TGraphErrors *grUnfolded_pre_fsr_dRp1 = new TGraphErrors(5, &meanMass_data_pre_fsr_dRp1[0], &meanPt_data_pre_fsr_dRp1[0], &meanMassStatErr_data_pre_fsr_dRp1[0], &meanPtStatErr_data_pre_fsr_dRp1[0]);
        grUnfolded_pre_fsr_dRp1->SetLineColor(kGray+2);
        grUnfolded_pre_fsr_dRp1->SetMarkerColor(kGray+2);
        grUnfolded_pre_fsr_dRp1->SetMarkerStyle(marker_);
        grUnfolded_pre_fsr_dRp1->SetMarkerSize(0.9);
        grUnfolded_pre_fsr_dRp1->SetLineStyle(2);
        grUnfolded_pre_fsr_dRp1->Draw("pCZ same");

        TGraphErrors *grUnfolded_pre_fsr = new TGraphErrors(5, &meanMass_data_pre_fsr[0], &meanPt_data_pre_fsr[0], &meanMassTotErr_data_pre_fsr[0], &meanPtTotErr_data_pre_fsr[0]);
        grUnfolded_pre_fsr->SetLineColor(kBlack);
        grUnfolded_pre_fsr->SetMarkerColor(kBlack);
        grUnfolded_pre_fsr->SetMarkerStyle(marker_);
        grUnfolded_pre_fsr->SetMarkerSize(0.9);
        grUnfolded_pre_fsr->SetLineStyle(1);
        grUnfolded_pre_fsr->Draw("pZ same");

        TGraphErrors *grMCAlt;
        if(altMC){
            grMCAlt = new TGraphErrors(5, &meanMass_mcAlt[0], &meanPt_mcAlt[0], &meanMassErr_mcAlt[0], &meanPtErr_mcAlt[0]);
            grMCAlt->SetLineColor(kRed);
            grMCAlt->SetMarkerColor(kRed);
            grMCAlt->SetMarkerStyle(marker_);
            grMCAlt->SetMarkerSize(1.);
            grMCAlt->SetLineStyle(1);
            grMCAlt->SetLineColor(kRed);
            grMCAlt->Draw("pe same");
        }

        //grUnfolded->Draw("pe same");

        TLegend* leg_ = new TLegend(0.53, 0.25, 0.9, 0.5,"","brNDC");
        leg_->SetTextSize(0.027);
        leg_->SetFillStyle(0); // transparent 
        leg_->SetBorderSize(0);
        leg_->AddEntry(grUnfolded, "Detector unfolded data ", "pe");
        leg_->AddEntry(grUnfolded_pre_fsr_dRp1, "QED FSR (dR<0.1)unfolded data ", "pe");
        leg_->AddEntry(grUnfolded_pre_fsr, "QED FSR unfolded data ", "pe");
        leg_->AddEntry(grMC, "DY MC at post FSR (aMC@NLO)", "pe");
        if(altMC) leg_->AddEntry(grMCAlt, "DY MC at pre FSR (Madgraph)", "p");
        leg_->Draw();

        TF1 *f1 = NULL;
        if(doFit){
            f1 = new TF1("f1", "[0]+[1]*log(x)", 42., 300.);
            f1->GetXaxis()->SetRangeUser(42., 300.);
            f1->SetLineColor(kBlack);
            f1->SetLineWidth(1);
            //grUnfolded->Fit(f1, "R0"); // R: fitting sub range
            grUnfolded_pre_fsr->Fit(f1, "R0"); // R: fitting sub range
            f1->Draw("same");
        }

        CMS_lumi( c1, 4, 11 );
        c1->SaveAs(outpdf + channel_name + ".pdf");
	delete grUnfolded;
	delete grMC;
	delete f1;
	delete c1;
}

void ISRUnfold::drawInputPlots(TString outpdf, TString var, int nthMassBin, TString sysName){

   gROOT->SetBatch();

   setTDRStyle();
   writeExtraText = true;
   extraText  = "work in progress";

   TString ibinMass;
   ibinMass.Form("%d", nthMassBin);

   TH1* hpt_temp_data;
   TH1F *ratio = NULL;
   
   hpt_temp_data   = nomPtUnfold->GetInput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);

   c1=new TCanvas("c1", "c1", 50, 50, 800, 800);
   c1->cd();
   gStyle->SetOptStat(0);

   TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
   pad1->SetBottomMargin(0.01);
   pad1->SetTopMargin(0.1);
   pad1->SetTicks(1);
   pad1->SetLogy();
   pad1->Draw();
   pad1->cd();

   hpt_temp_data->SetTitle("");
   hpt_temp_data->Draw("p9histe");
   hpt_temp_data->SetMarkerStyle(20);
   hpt_temp_data->SetMarkerSize(.7);
   hpt_temp_data->SetLineColor(kBlack);
   hpt_temp_data->GetYaxis()->SetTitle("Events/bin");

   TH1* hpt_sys_temp;
   int sysSize = sysPtUnfold[sysName].size();
   for(int i = 0; i < sysSize; i++){
           //if((i==5 || i==7) && sysName=="Scale") continue;

           TString isys;
           isys.Form("%d", i);

           TH1 * hsyspt_temp = sysPtUnfold[sysName].at(i)->GetInput("hunfolded_pt_systemp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
           hpt_sys_temp = ((TH1F*)hsyspt_temp->Clone("pt_temp"));
           hpt_sys_temp->Draw("histsame");
           hpt_sys_temp->SetLineColor(kBlack);
           hpt_sys_temp->SetLineStyle(2);

           delete hsyspt_temp;
   }

   TString mean_nom;
   mean_nom.Form("%.5f", hpt_temp_data->GetMean());

   TLegend* leg_nom = new TLegend(0.45, 0.70, 0.75, 0.9,"","brNDC");
   leg_nom->SetNColumns(2);
   leg_nom->SetTextSize(0.055);
   leg_nom->SetFillStyle(0);
   leg_nom->SetBorderSize(0);

   leg_nom->AddEntry(hpt_temp_data, "Bkg subtracted data (mean: " + mean_nom + ")", "pl");
   leg_nom->Draw();
   c1->cd();

   TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.4);
   pad2->SetTopMargin(0.05);
   pad2->SetBottomMargin(0.2);
   pad2->SetTicks(1);
   pad2->SetGridy(1);
   pad2->Draw();
   pad2->cd();

   for(int i = 0; i < sysSize; i++){
        //if((i==5 || i==7) && sysName=="Scale") continue;

        TString isys;
        isys.Form("%d", i); 

        TH1 * hsyspt_temp = sysPtUnfold[sysName].at(i)->GetInput("hunfolded_pt_systemp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
        ratio = ((TH1F*)hsyspt_temp->Clone("pt_temp"));
   	ratio->Divide(hpt_temp_data);
        if(i==0 ){
   	    ratio->Draw("hist");
            ratio->GetYaxis()->SetTitle("Systematic/ Nominal input");
   	    ratio->GetXaxis()->SetTitle("p_{T} at pre FSR(GeV)");
   	    ratio->SetMinimum(0.9);
   	    ratio->SetMaximum(1.1);
   	    ratio->SetTitle("");
   	    ratio->GetXaxis()->SetTitleOffset(1.5);
   	    ratio->GetYaxis()->SetNdivisions(515);
            ratio->SetLineColor(kBlack);
            ratio->SetLineStyle(2);
   	}
        else{
   	    ratio->Draw("histsame");
            ratio->SetLineStyle(2);
   	}

        delete hsyspt_temp;
   }   

   CMS_lumi( c1, 4, 0 );
   c1->cd();
   c1->SaveAs(outpdf+"_input_"+ibinMass+".pdf");

   delete hpt_temp_data;
   delete ratio;
   delete pad1;
   delete pad2;
   delete c1;
}

// https://root.cern.ch/root/html/tutorials/graphs/graphtext.C.html
void ISRUnfold::drawtext(TGraph *g)
{
   Int_t i,n;
   Double_t x,y;
   TLatex *l;

   n = g->GetN();
   for (i=0; i<n; i++) {
       g->GetPoint(i,x,y);
       l = new TLatex(x-0.1,y+0.02,Form("%d",i+1));
       l->SetTextSize(0.015);
       l->SetTextFont(42);
       l->SetTextAlign(21);
       l->Draw();
   }
}

void ISRUnfold::drawSysPlots(TString outpdf, int nthMassBin, TString sysName, bool detector_unfold){

        gROOT->SetBatch();

        TString ibinMass;
        ibinMass.Form("%d", nthMassBin);

        c1 = new TCanvas("c1","c1", 50, 50, 850, 700);
        c1->cd();
        gStyle->SetOptFit(0);

        c1->SetBottomMargin(0.2);
        c1->SetTopMargin(0.08);
        c1->SetTicks(1);

        int sysSize = 0;
        double amaxMass = 0.;
        double aminMass = 0.;
        double amaxPt = 0.;
        double aminPt = 0.;

        TGraphErrors *grNominal;
        TGraphErrors *grStatErr;
        TGraph *grSys;

        if(detector_unfold){
            // set axis min and max range
	    sysSize = sysPtUnfold[sysName].size();
            amaxMass=meanMass_data[nthMassBin]+ 2. * meanMassTotErr_data[nthMassBin];
            aminMass=meanMass_data[nthMassBin]- 2. * meanMassTotErr_data[nthMassBin];
            amaxPt=meanPt_data[nthMassBin]+ 2. * meanPtTotErr_data[nthMassBin];
            aminPt=meanPt_data[nthMassBin]- 2. * meanPtTotErr_data[nthMassBin];

            for(int i=0;i<sysSize;i++) {
               amaxMass=TMath::Max(amaxMass, (meanMass_sysdata.at(nthMassBin)[sysName])[i]);
               aminMass=TMath::Min(aminMass, (meanMass_sysdata.at(nthMassBin)[sysName])[i]);

               amaxPt=TMath::Max(amaxPt, (meanPt_sysdata.at(nthMassBin)[sysName])[i]);
               aminPt=TMath::Min(aminPt, (meanPt_sysdata.at(nthMassBin)[sysName])[i]);
            }

            // nominal point with total uncertainty
            TGaxis::SetMaxDigits(4);   
            grNominal = new TGraphErrors(1, &meanMass_data[nthMassBin], &meanPt_data[nthMassBin], &meanMassTotErr_data[nthMassBin], &meanPtTotErr_data[nthMassBin]);
            grNominal->SetLineColor(kBlack);
            grNominal->SetMarkerColor(kBlack);
            grNominal->SetMarkerStyle(20);
            grNominal->SetMarkerSize(1.);
            grNominal->SetLineStyle(1);
            grNominal->Draw("ape");
            grNominal->GetYaxis()->SetRangeUser(aminPt*0.995,amaxPt*1.005);
            grNominal->GetXaxis()->SetLimits(aminMass*0.995,amaxMass*1.005);
            grNominal->GetXaxis()->SetNdivisions(505, false);
            grNominal->GetYaxis()->SetTitle("Average p_{T} (GeV)");
            grNominal->GetXaxis()->SetTitle("Average Mass (GeV)");

	    // meanPtStatErr_data
            grStatErr = new TGraphErrors(1, &meanMass_data[nthMassBin], &meanPt_data[nthMassBin], &meanMassStatErr_data[nthMassBin], &meanPtStatErr_data[nthMassBin]);
            grStatErr->SetFillColor(kBlack);;
	    grStatErr->SetFillColorAlpha(kBlack,0.3);
            grStatErr->Draw("sameE2");

            grSys = new TGraph(sysSize, &(meanMass_sysdata.at(nthMassBin)[sysName])[0], &(meanPt_sysdata.at(nthMassBin)[sysName])[0] );
            grSys->SetLineColor(kRed);
            grSys->SetMarkerColor(kRed);
            grSys->SetMarkerStyle(20);
            grSys->SetMarkerSize(0.8);
            grSys->SetLineStyle(1);
            grSys->Draw("pe same ");
	    drawtext(grSys);

            grNominal->Draw("pe same");

            TLegend* leg = new TLegend(0.6, 0.70, 0.85, 0.9,"","brNDC");
            leg->SetTextSize(0.02);
            leg->SetFillStyle(0);
            leg->SetBorderSize(0);

            leg->AddEntry(grNominal, "Measurment with total uncertainty", "lep");
            leg->AddEntry(grSys,      "Systematic source : " + sysName, "pl");
            leg->AddEntry(grStatErr,  "Statistical uncertainty" , "f");
            leg->Draw();
        }
        else{

            // set axis min and max range
            sysSize = sysPtFSRUnfold[sysName].size();
            amaxMass=meanMass_data_pre_fsr[nthMassBin] + 2. * meanMassTotErr_data_pre_fsr[nthMassBin];
            aminMass=meanMass_data_pre_fsr[nthMassBin] - 2. * meanMassTotErr_data_pre_fsr[nthMassBin];
            amaxPt=meanPt_data_pre_fsr[nthMassBin] + 2. * meanPtTotErr_data_pre_fsr[nthMassBin];
            aminPt=meanPt_data_pre_fsr[nthMassBin] - 2. * meanPtTotErr_data_pre_fsr[nthMassBin];

            for(int i=0;i<sysSize;i++) {
               amaxMass=TMath::Max(amaxMass, (meanMass_sysdata_pre_fsr.at(nthMassBin)[sysName])[i]);
               aminMass=TMath::Min(aminMass, (meanMass_sysdata_pre_fsr.at(nthMassBin)[sysName])[i]);

               amaxPt=TMath::Max(amaxPt, (meanPt_sysdata_pre_fsr.at(nthMassBin)[sysName])[i]);
               aminPt=TMath::Min(aminPt, (meanPt_sysdata_pre_fsr.at(nthMassBin)[sysName])[i]);
            }

            // nominal point with total uncertainty
            TGaxis::SetMaxDigits(4);
            grNominal = new TGraphErrors(1, &meanMass_data_pre_fsr[nthMassBin], &meanPt_data_pre_fsr[nthMassBin], &meanMassTotErr_data_pre_fsr[nthMassBin], &meanPtTotErr_data_pre_fsr[nthMassBin]);
            grNominal->SetLineColor(kBlack);
            grNominal->SetMarkerColor(kBlack);
            grNominal->SetMarkerStyle(20);
            grNominal->SetMarkerSize(1.);
            grNominal->SetLineStyle(1);
            grNominal->Draw("ape");
            grNominal->GetYaxis()->SetRangeUser(aminPt*0.995,amaxPt*1.005);
            grNominal->GetXaxis()->SetLimits(aminMass*0.995,amaxMass*1.005);
            grNominal->GetXaxis()->SetNdivisions(505, false);
            grNominal->GetYaxis()->SetTitle("Average p_{T} (GeV)");
            grNominal->GetXaxis()->SetTitle("Average Mass (GeV)");

            // meanPtStatErr_data
            grStatErr = new TGraphErrors(1, &meanMass_data_pre_fsr[nthMassBin], &meanPt_data_pre_fsr[nthMassBin], &meanMassStatErr_data_pre_fsr[nthMassBin], &meanPtStatErr_data_pre_fsr[nthMassBin]);
            grStatErr->SetFillColor(kBlack);;
            grStatErr->SetFillColorAlpha(kBlack,0.3);
            grStatErr->Draw("sameE2");

            grSys = new TGraph(sysSize, &(meanMass_sysdata_pre_fsr.at(nthMassBin)[sysName])[0], &(meanPt_sysdata_pre_fsr.at(nthMassBin)[sysName])[0] );
            grSys->SetLineColor(kRed);
            grSys->SetMarkerColor(kRed);
            grSys->SetMarkerStyle(20);
            grSys->SetMarkerSize(0.8);
            grSys->SetLineStyle(1);
            grSys->Draw("pe same ");
            drawtext(grSys);

            grNominal->Draw("pe same");

            TLegend* leg = new TLegend(0.6, 0.70, 0.85, 0.9,"","brNDC");
            leg->SetTextSize(0.02);
            leg->SetFillStyle(0);
            leg->SetBorderSize(0);

            leg->AddEntry(grNominal, "Measurment with total uncertainty", "lep");
            leg->AddEntry(grSys,      "Systematic source : " + sysName, "pl");
            leg->AddEntry(grStatErr,  "Statistical uncertainty" , "f");
            leg->Draw();

        }

        CMS_lumi( c1, 4, 11 );
        if (detector_unfold) c1->SaveAs(outpdf+"_"+ibinMass+"_"+sysName+".pdf");
        else c1->SaveAs(outpdf+"_"+ibinMass+"_pre_fsr_"+sysName+".pdf");

        delete grNominal, grStatErr, grSys;
	delete c1;
}

void ISRUnfold::drawNominalRecoPlots(TString outpdf, TString filepath, TString var, int nthMassBin, TString sysName){

        TString ibinMass;
        ibinMass.Form("%d", nthMassBin);

        gROOT->SetBatch();

        setTDRStyle();
        writeExtraText = true;
        extraText  = "work in progress";

	// get binning definition
        const TUnfoldBinning* temp_binning = nomPtUnfold->GetInputBinning("Rec_Pt");

	TFile* filein = new TFile(filepath);
        TH1* hdysig = (TH1*)filein->Get("hPtRecnominal"); // get DY 
        TH1* hdysigNoUO = temp_binning->ExtractHistogram("hdysig", hdysig, 0, kTRUE, "pt[UO];mass[UOC"+ibinMass+"]");

        TH1* hpt_temp_data;
        TH1F *ratio;

	// bkg subracted data
        hpt_temp_data   = nomPtUnfold->GetInput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);

        c1=new TCanvas("c1", "c1", 50, 50, 800, 800);
        c1->cd();
        gStyle->SetOptStat(0);

        TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
        pad1->SetBottomMargin(0.01);
        pad1->SetTopMargin(0.1);
        pad1->SetTicks(1);
        pad1->SetLogy();
        pad1->Draw();
        pad1->cd();

        hpt_temp_data->SetTitle("");
        hpt_temp_data->Draw("p9histe");
        hpt_temp_data->SetMarkerStyle(20);
        hpt_temp_data->SetMarkerSize(.7);
        hpt_temp_data->SetLineColor(kBlack);
        hpt_temp_data->GetYaxis()->SetTitle("Events/bin");
        hdysigNoUO->Draw("hist same");
        hdysigNoUO->SetLineColor(kRed);


        TString mean_nom;
        mean_nom.Form("%.5f", hpt_temp_data->GetMean());

        TLegend* leg_nom = new TLegend(0.45, 0.70, 0.75, 0.9,"","brNDC");
        leg_nom->SetNColumns(2);
        leg_nom->SetTextSize(0.055);
        leg_nom->SetFillStyle(0);
        leg_nom->SetBorderSize(0);

        leg_nom->AddEntry(hpt_temp_data, "Bkg subtracted data (mean: " + mean_nom + ")", "pl");
        leg_nom->Draw();

        TLatex chi2_norm;
        TString chi2_;

        // TODO add Mass option 
        if(var == "Pt" ){
                chi2_.Form("%f", Chi2Test(hpt_temp_data, hdysigNoUO));
                chi2_norm.DrawLatexNDC(0.2, 0.3, "#chi^{2}/NDOF= " + chi2_);
        }

        c1->cd();

        TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.4);
        pad2->SetTopMargin(0.05);
        pad2->SetBottomMargin(0.2);
        pad2->SetTicks(1);
        pad2->SetGridy(1);
        pad2->Draw();

        pad2->cd();

	ratio = ((TH1F*)hpt_temp_data->Clone("pt_temp"));
	ratio->Divide(hdysigNoUO);
        ratio->Draw("hist");
        ratio->GetYaxis()->SetTitle("Data (bkg sub)/ DY MC");
        ratio->GetXaxis()->SetTitle("p_{T} (GeV)");
        ratio->SetMinimum(0.8);
        ratio->SetMaximum(1.2);
        ratio->SetTitle("");
        ratio->GetXaxis()->SetTitleOffset(1.5);
        ratio->GetYaxis()->SetNdivisions(515);
        ratio->SetLineColor(kBlack);
        //ratio->SetLineStyle(2);

        CMS_lumi( c1, 4, 0 );
        c1->cd();
        c1->SaveAs(outpdf+"_input_"+ibinMass+".pdf");

        delete hpt_temp_data;
	delete hdysigNoUO;
        delete ratio;
        delete pad1;
        delete pad2;
        delete c1;

}

void ISRUnfold::studyFSRDRPlots(TString outpdf, TString var, int nthMassBin){

        gROOT->SetBatch();

        setTDRStyle();
        writeExtraText = true;
        extraText  = "work in progress";

        TString ibinMass;
        ibinMass.Form("%d", nthMassBin);

        TH1* hunfolded_data = NULL;
        TH1* hpreFSR_dressed0p1 = NULL;
        TH1F *ratio = NULL;

        // get nominal unfoled result
        if(var == "Pt" ){
                hunfolded_data  = nomPtUnfold->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
                hpreFSR_dressed0p1   =    sysPtUnfold["FSRDR"].at(0)->GetBias("histMCTruth_pt_tempAlt",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
        }
        if(var == "Mass" ){
                hunfolded_data  = nomMassUnfold->GetOutput("hunfolded_mass_temp",0,0,"*[UO]",kTRUE);
                hpreFSR_dressed0p1   =    sysMassUnfold["FSRDR"].at(0)->GetOutput("hunfolded_mass_systemp",0,0,"*[UO]",kTRUE);
        }

        ratio= ((TH1F*)hunfolded_data->Clone("ratio"));
        ratio->Divide(hpreFSR_dressed0p1);

        c1=new TCanvas("c1", "c1", 50, 50, 800, 800);
        c1->cd();
        gStyle->SetOptStat(0);

        TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
        pad1->SetBottomMargin(0.01);
        pad1->SetTopMargin(0.1);
        pad1->SetTicks(1);
        pad1->SetLogy();
        pad1->Draw();
        pad1->cd();

        hunfolded_data->SetTitle("");
        hunfolded_data->Draw("p9histe");
        hpreFSR_dressed0p1->Draw("histsame");
        hunfolded_data->SetMarkerStyle(20);
        hunfolded_data->SetMarkerSize(.7);
        hunfolded_data->SetLineColor(kBlack);
        hpreFSR_dressed0p1->SetLineColor(kRed);
        hunfolded_data->GetYaxis()->SetTitle("Events/bin");
        hunfolded_data->SetMinimum(10.);

        c1->cd();

        TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.4);
        pad2->SetTopMargin(0.05);
        pad2->SetBottomMargin(0.2);
        pad2->SetTicks(1);
        pad2->SetGridy(1);
        pad2->Draw();

        pad2->cd();
        ratio->Draw("hist");
        ratio->GetYaxis()->SetTitle("#DeltaR=X/ #DeltaR=0.1");
        ratio->GetXaxis()->SetTitle("p_{T} at pre FSR(GeV)");
        ratio->SetMinimum(0.5);
        ratio->SetMaximum(1.5);
        ratio->SetTitle("");
        ratio->GetXaxis()->SetTitleOffset(1.5);
        ratio->GetYaxis()->SetNdivisions(515);
        ratio->SetLineColor(kRed);

        TH1* hratio_dr_temp;
        int sysSize = sysPtUnfold["FSRDR"].size();
        for(int i = 1; i < sysSize; i++){

		TH1* htemp = NULL;
        	if(var == "Pt" ){
        	        htemp   =    sysPtUnfold["FSRDR"].at(i)->GetBias("histMCTruth_pt_tempAlt_",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
        	}
        	if(var == "Mass" ){
        	        htemp   =    sysMassUnfold["FSRDR"].at(i)->GetOutput("hunfolded_mass_systemp_",0,0,"*[UO]",kTRUE);
        	}
		hratio_dr_temp = ((TH1F*)htemp->Clone("htemp"));
		hratio_dr_temp->Divide(hpreFSR_dressed0p1);
                hratio_dr_temp->Draw("histsame");
                hratio_dr_temp->SetLineColor(i);
                hratio_dr_temp->SetLineStyle(2);
                delete htemp;
	}

        CMS_lumi( c1, 4, 0 );
        c1->cd();
        c1->SaveAs(outpdf+"_"+ibinMass+"_"+var+".pdf");

	//delete hratio_dr_temp;
        delete hunfolded_data;
        delete hpreFSR_dressed0p1;
        delete pad1;
        delete pad2;
        delete c1;
}

void ISRUnfold::drawClosurePlots(TString outpdf, TString var, int nthMassBin){

    const TUnfoldBinningV17* temp_binning = nomPtUnfold->GetOutputBinning("Gen_Pt");
    // get mass bin definition from (pt, mass) bin definition
    const TVectorD* temp_tvecd = temp_binning->GetDistributionBinning(1);
    const Double_t* massBins = temp_tvecd->GetMatrixArray();

    gROOT->SetBatch();

    setTDRStyle();
    writeExtraText = true;
    extraText  = "work in progress";

    TString ibinMass;
    ibinMass.Form("%d", nthMassBin);

    TH1* hunfolded_data = NULL;
    TH1* hpreFSR_mc = NULL;
    TH1F *ratio = NULL;

    // get nominal unfoled result
    if(var == "Pt" ){
        hunfolded_data  = sysPtUnfold["Closure"].at(0)->GetOutput("hunfolded_pt_closure",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
        hpreFSR_mc   =    sysPtUnfold["Closure"].at(0)->GetBias("histMCTruth_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
    }
    if(var == "Mass" ){
        hunfolded_data  = sysMassUnfold["Closure"].at(0)->GetOutput("hunfolded_mass_closure",0,0,"mass[UO];pt[UOC0]",kTRUE);
        hunfolded_data->GetXaxis()->SetRange(hunfolded_data->GetXaxis()->FindBin(massBins[nthMassBin]+0.01),hunfolded_data->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));

        hpreFSR_mc   = sysMassUnfold["Closure"].at(0)->GetBias("histMCTruth_mass_temp",0,0,"mass[UO];pt[UOC0]",kTRUE);
        hpreFSR_mc->GetXaxis()->SetRange(hunfolded_data->GetXaxis()->FindBin(massBins[nthMassBin]+0.01),hunfolded_data->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));
    }

    ratio= ((TH1F*)hunfolded_data->Clone("ratio"));
    ratio->Divide(hpreFSR_mc);

    c1=new TCanvas("c1", "c1", 50, 50, 800, 800);
    c1->cd();
    gStyle->SetOptStat(0);

    TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
    pad1->SetBottomMargin(0.01);
    pad1->SetTopMargin(0.1);
    pad1->SetTicks(1);
    pad1->SetLogy();
    pad1->Draw();
    pad1->cd();

    hunfolded_data->SetTitle("");
    hunfolded_data->Draw("p9histe");
    hpreFSR_mc->Draw("histsame");

    hunfolded_data->SetMarkerStyle(20);
    hunfolded_data->SetMarkerSize(.7);
    hunfolded_data->SetLineColor(kBlack);
    hunfolded_data->GetYaxis()->SetTitle("Events/bin");
    hunfolded_data->SetMinimum(10.);

    hpreFSR_mc->SetLineColor(kRed);

    TString mean_nom;
    TString meanMC_nom;
    mean_nom.Form("%.5f", hunfolded_data->GetMean());
    meanMC_nom.Form("%.5f", hpreFSR_mc->GetMean());

    TLegend* leg_nom = new TLegend(0.45, 0.4, 0.75, 0.6,"","brNDC");
    //leg_nom->SetNColumns(2);
    leg_nom->SetTextSize(0.055);
    leg_nom->SetFillStyle(0);
    leg_nom->SetBorderSize(0);
    leg_nom->AddEntry(hunfolded_data, "Unfolded MC (mean: " + mean_nom + ")", "pl");
    leg_nom->AddEntry(hpreFSR_mc, "Gen MC (mean: " + meanMC_nom + ")", "l");

    TLatex chi2_norm;
    TString chi2_;

    // TODO add Mass option 
    if(var == "Pt" ){
            chi2_.Form("%f",DoFit("Pt", nthMassBin));
            chi2_norm.DrawLatexNDC(0.2, 0.3, "#chi^{2}/NDOF= " + chi2_);
    }
    leg_nom->Draw();

    c1->cd();

    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.4);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.2);
    pad2->SetTicks(1);
    pad2->SetGridy(1);
    pad2->Draw();

    pad2->cd();
    ratio->Draw("histe");
    ratio->GetYaxis()->SetTitle("Unfold/ MC Gen");
    ratio->GetXaxis()->SetTitle("p_{T} at pre FSR(GeV)");
    ratio->SetMinimum(0.8);
    ratio->SetMaximum(1.2);
    ratio->SetTitle("");
    ratio->GetXaxis()->SetTitleOffset(1.5);
    ratio->GetYaxis()->SetNdivisions(515);

    CMS_lumi( c1, 4, 0 );
    c1->cd();
    c1->SaveAs(outpdf+"_"+ibinMass+"_"+var+".pdf");

    delete hunfolded_data;
    delete hpreFSR_mc;
    delete pad1;
    delete pad2;
    delete c1;

}

void ISRUnfold::drawNominalPlots(TString outpdf, TString var, int nthMassBin, TString sysName, bool systematic, bool isFSRUnfold){

        const TUnfoldBinningV17* temp_binning = nomPtUnfold->GetOutputBinning("Gen_Pt");
        const TUnfoldBinningV17* temp_binning_rec = nomPtUnfold->GetInputBinning("Rec_Pt");
        const TUnfoldBinningV17* temp_binning_rec_mass = nomMassUnfold->GetInputBinning("Rec_Mass");
        const TUnfoldBinningV17* temp_binning_gen_mass = nomMassUnfold->GetOutputBinning("Gen_Mass");
        // get mass bin definition from (pt, mass) bin definition
        const TVectorD* temp_tvecd = temp_binning->GetDistributionBinning(1);
        const Double_t* massBins = temp_tvecd->GetMatrixArray();

        gROOT->SetBatch();

	setTDRStyle();
	writeExtraText = true;
	extraText  = "work in progress";

        TString ibinMass;
        ibinMass.Form("%d", nthMassBin);

        TH1* hunfolded_data = NULL;
        TH1* h_data_detector = NULL;
        TH1* h_mc_detector = NULL;
	TH1* hunfolded_sys_err = NULL;
        TH1* hpreFSR_mc = NULL;

        TH1F *ratio = NULL;
        TH1F *ratio_detector = NULL;
        TH1F *ratio_sys_err = NULL;

        // get DY MC histograms at detector
        TFile* filein = new TFile(hist_file_path);
        TH1* hdysig_pt = NULL; // get DY 
        if(channel_name == "electron"){
            hdysig_pt = (TH1*)filein->Get("detector_level/hist_ptll/histo_DYJetsToEEnominal"); // get DY 
            hdysig_pt->Add((TH1*)filein->Get("detector_level/hist_ptll/histo_DYJets10to50ToEEnominal"));
            if(do_normalization) hdysig_pt->Scale(0.959939);
        }
        if(channel_name == "muon"){
            hdysig_pt = (TH1*)filein->Get("detector_level/hist_ptll/histo_DYJetsToMuMunominal"); // get DY 
            hdysig_pt->Add((TH1*)filein->Get("detector_level/hist_ptll/histo_DYJets10to50ToMuMunominal"));
            if(do_normalization) hdysig_pt->Scale(0.953026);
        }

        TH1* hdysig_mass = NULL;
        if(channel_name == "electron"){ 
            hdysig_mass = (TH1*)filein->Get("detector_level/hist_mll/histo_DYJetsToEEnominal");
            hdysig_mass->Add((TH1*)filein->Get("detector_level/hist_mll/histo_DYJets10to50ToEEnominal"));
            if(do_normalization) hdysig_mass->Scale(0.959939);
        }
        if(channel_name == "muon"){
            hdysig_mass = (TH1*)filein->Get("detector_level/hist_mll/histo_DYJetsToMuMunominal");
            hdysig_mass->Add((TH1*)filein->Get("detector_level/hist_mll/histo_DYJets10to50ToMuMunominal"));
            if(do_normalization) hdysig_mass->Scale(0.953026);
        }

        //TFile* filein = new TFile("/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/output/2016/electron/DY_FSR_v3.root");
        //TH1* hdysig_pt = NULL; // get DY 
        //if(channel_name == "electron"){
        //    hdysig_pt = (TH1*)filein->Get("fiducial_phase_post_FSR/hist_ptll_post_fsr/histo_DYJetsnominal"); // get DY 
        //    hdysig_pt->Add((TH1*)filein->Get("fiducial_phase_post_FSR/hist_ptll_post_fsr/histo_DYJets10to50nominal"));
        //}
        //if(channel_name == "muon"){
        //    hdysig_pt = (TH1*)filein->Get("fiducial_phase_post_FSR/hist_ptll_post_fsr/histo_DYJetsToMuMunominal"); // get DY 
        //    hdysig_pt->Add((TH1*)filein->Get("fiducial_phase_post_FSR/hist_ptll_post_fsr/histo_DYJets10to50ToMuMunominal"));
        //}

        //TH1* hdysig_mass = NULL;
        //if(channel_name == "electron"){
        //    hdysig_mass = (TH1*)filein->Get("fiducial_phase_post_FSR/hist_mll_post_fsr/histo_DYJetsnominal");
        //    hdysig_mass->Add((TH1*)filein->Get("fiducial_phase_post_FSR/hist_mll_post_fsr/histo_DYJets10to50nominal"));
        //}
        //if(channel_name == "muon"){
        //    hdysig_mass = (TH1*)filein->Get("fiducial_phase_post_FSR/hist_mll_post_fsr/histo_DYJetsToMuMunominal");
        //    hdysig_mass->Add((TH1*)filein->Get("fiducial_phase_post_FSR/hist_mll_post_fsr/histo_DYJets10to50ToMuMunominal"));
        //}

	// get nominal unfoled result
	if(var == "Pt" ){
            if(!isFSRUnfold){
            hunfolded_data  = nomPtUnfold->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
            h_data_detector  = nomPtUnfold->GetInput("hdata_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
	    hunfolded_sys_err= ((TH1F*)hunfolded_data->Clone("sysErr")); 
            // get gen histogram from response matrix
	    hpreFSR_mc = nomPtUnfold->GetBias("histMCTruth_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
            h_mc_detector = temp_binning_rec->ExtractHistogram("hdysig", hdysig_pt, 0, kTRUE, "pt[UO];mass[UOC"+ibinMass+"]");
            }
            else{
                hunfolded_data  = nomPtFSRUnfold->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
                hunfolded_sys_err= ((TH1F*)hunfolded_data->Clone("sysErr"));
                // get gen histogram from response matrix
                hpreFSR_mc = nomPtFSRUnfold->GetBias("histMCTruth_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
                
                h_data_detector  = nomPtFSRUnfold->GetInput("hdata_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
                h_mc_detector  = nomPtUnfold->GetBias("hunfolded_pt_temp_",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
            }
	}
	if(var == "Mass" ){
                if(!isFSRUnfold){
                    hunfolded_data  = nomMassUnfold->GetOutput("hunfolded_mass_temp",0,0,"mass[UO];pt[UOC0]",kTRUE);
                    h_data_detector = nomMassUnfold->GetInput("hdata_mass_temp",0,0,"mass[UO];pt[UOC0]",kTRUE);
		    hunfolded_sys_err= (TH1F*)hunfolded_data->Clone("sysErr"); 

		    hunfolded_data->GetXaxis()->SetRange(hunfolded_data->GetXaxis()->FindBin(massBins[nthMassBin]+0.01),hunfolded_data->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));
		    h_data_detector->GetXaxis()->SetRange(h_data_detector->GetXaxis()->FindBin(massBins[nthMassBin]+0.01),h_data_detector->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));
                    // get gen histogram from response matrix
		    hpreFSR_mc   = nomMassUnfold->GetBias("histMCTruth_mass_temp",0,0,"mass[UO];pt[UOC0]",kTRUE);
                    h_mc_detector = temp_binning_rec_mass->ExtractHistogram("hdysig_mass", hdysig_mass, 0, kTRUE, "mass[UO];pt[UOC0]");

                    hpreFSR_mc->GetXaxis()->SetRange(hunfolded_data->GetXaxis()->FindBin(massBins[nthMassBin]+0.01), hunfolded_data->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));
                    h_mc_detector->GetXaxis()->SetRange(h_mc_detector->GetXaxis()->FindBin(massBins[nthMassBin]+0.01), h_mc_detector->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));

                    if(nthMassBin == 2) cout << "normalisation factor: " << h_data_detector->Integral()/h_mc_detector->Integral() << endl;
                }
                else{
                    hunfolded_data  = nomMassFSRUnfold->GetOutput("hunfolded_mass_temp",0,0,"mass[UO];pt[UOC0]",kTRUE);
                    hunfolded_data->GetXaxis()->SetRange(hunfolded_data->GetXaxis()->FindBin(massBins[nthMassBin]+0.01),hunfolded_data->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));
                    hunfolded_sys_err= ((TH1F*)hunfolded_data->Clone("sysErr"));
               
                    h_data_detector = nomMassFSRUnfold->GetInput("hdata_mass_temp",0,0,"mass[UO];pt[UOC0]",kTRUE);
                    h_data_detector->GetXaxis()->SetRange(h_data_detector->GetXaxis()->FindBin(massBins[nthMassBin]+0.01),h_data_detector->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));
                    // get gen histogram from response matrix
                    hpreFSR_mc   = nomMassFSRUnfold->GetBias("histMCTruth_mass_temp",0,0,"mass[UO];pt[UOC0]",kTRUE);
                    h_mc_detector  = nomMassUnfold->GetBias("hunfolded_mass_temp_",0,0,"mass[UO];pt[UOC0]",kTRUE);

                    hpreFSR_mc->GetXaxis()->SetRange(hunfolded_data->GetXaxis()->FindBin(massBins[nthMassBin]+0.01),hunfolded_data->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));
                    h_mc_detector->GetXaxis()->SetRange(hunfolded_data->GetXaxis()->FindBin(massBins[nthMassBin]+0.01),hunfolded_data->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));
                }
	}

        bool data_over_mc = false;
        if(data_over_mc){
            ratio = (TH1F*)hunfolded_data->Clone("ratio");
            ratio_detector = (TH1F*)h_data_detector->Clone("ratio_detector");
            ratio_sys_err= (TH1F*)hunfolded_data->Clone("ratio_sys");
            ratio->Divide(hpreFSR_mc);
            ratio_detector->Divide(h_mc_detector);
            ratio_sys_err->Divide(hpreFSR_mc);
        }
        else{
            ratio = (TH1F*)hpreFSR_mc->Clone("ratio");
            ratio_detector = (TH1F*)h_mc_detector->Clone("ratio_detector");
            ratio_sys_err= (TH1F*)hpreFSR_mc->Clone("ratio_sys");
            ratio->Divide(hunfolded_data);
            ratio_detector->Divide(h_data_detector);
            ratio_sys_err->Divide(hunfolded_data);
        }

        c1=new TCanvas("c1", "c1", 50, 50, 800, 800);
        c1->cd();
        gStyle->SetOptStat(0);

        TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
        pad1->SetBottomMargin(0.01);
        pad1->SetTopMargin(0.1);
        pad1->SetTicks(1);
        pad1->SetLogy();
        pad1->Draw();
        pad1->cd();

        hunfolded_data->SetTitle("");
        hunfolded_data->Draw("p9histe");
        h_data_detector->Draw("p9histesame");
        hpreFSR_mc->Draw("histsame");
        h_mc_detector->Draw("histsame");

        hunfolded_data->SetMarkerStyle(20);
        hunfolded_data->SetMarkerSize(.7);
        hunfolded_data->SetLineColor(kBlack);
        h_data_detector->SetMarkerColor(kGray+2);
        h_data_detector->SetLineColor(kGray+2);
        h_data_detector->SetMarkerSize(.7);
        hunfolded_data->GetYaxis()->SetTitle("Events/bin");
        hunfolded_data->SetMinimum(10.);
        hunfolded_data->SetMaximum(9e9);

        hpreFSR_mc->SetLineColor(kRed);
        h_mc_detector->SetLineColor(kMagenta);

        TString mean_nom;
        mean_nom.Form("%.5f", hunfolded_data->GetMean());

        //TLegend* leg_nom = new TLegend(0.45, 0.7, 0.9, 0.85,"","brNDC");
        TLegend* leg_nom = new TLegend(0.6, 0.65, 0.9, 0.85,"","brNDC");
        leg_nom->SetTextFont(63);
        leg_nom->SetTextSize(18);
        //leg_nom->SetTextSize(0.04);
        leg_nom->SetFillStyle(0);
        leg_nom->SetBorderSize(0);
        //leg_nom->SetNColumns(2);
        //leg_nom->AddEntry(hunfolded_data, "#splitline{Unfolded data}{(mean: " + mean_nom + ")}", "pl");
        if(!isFSRUnfold){
            leg_nom->AddEntry(hunfolded_data, "Detector unfolded data", "pl");
            leg_nom->AddEntry(hpreFSR_mc, "MC at post FSR", "l");

            leg_nom->AddEntry(h_data_detector, "Detector level data", "pl");
            leg_nom->AddEntry(h_mc_detector, "MC at detector", "l");
        }
        else{
            leg_nom->AddEntry(hunfolded_data, "QED FSR unfolded data", "pl");
            leg_nom->AddEntry(hpreFSR_mc, "MC at pre FSR", "l");
            
            leg_nom->AddEntry(h_data_detector, "Detector unfodled data", "pl");
            leg_nom->AddEntry(h_mc_detector, "MC at post FSR", "l");
        }
        leg_nom->Draw();

        TString cut_info;
        TString lepton_type;
        if(channel_name=="electron") lepton_type = "ee"; 
        else lepton_type = "#mu#mu"; 

        if(var=="Pt"){
            TString low_bound_, upper_bound_;
            low_bound_.Form("%d", (int)massBins[nthMassBin]);
            upper_bound_.Form("%d", (int)massBins[nthMassBin+1]);
            cut_info = low_bound_ + " < M(" + lepton_type + ") < " + upper_bound_ + " (GeV)";
        }
        else{
            cut_info = "p_{T}(" + lepton_type + ") < 100 (GeV) ";
        }

        TLatex cut_info_;
        cut_info_.SetTextFont(63);
        cut_info_.SetTextSize(20);
        cut_info_.DrawLatexNDC(0.2, 0.8, cut_info);
        

        TLatex chi2_;
        chi2_.SetTextFont(63);
        chi2_.SetTextSize(20);
        TString after_unfolding_chi2, before_unfolding_chi2;

        // this chi2 calculation don't consider correlation between bins
        double chi2_without_correlation = Chi2Test(h_data_detector, h_mc_detector);

        if(!isFSRUnfold){
            if(var == "Pt" ){
                after_unfolding_chi2.Form("%.2f", DoFit("Pt", nthMassBin));
                before_unfolding_chi2.Form("%.2f", chi2_without_correlation);
                chi2_.DrawLatexNDC(0.2, 0.15, "#splitline{Norm. #chi_{Detector unfolded}^{2}: " + after_unfolding_chi2 + "}{Norm. #chi^{2}: " + before_unfolding_chi2 + "}");
            }
            else{
                after_unfolding_chi2.Form("%.2f", DoFit("Mass", nthMassBin));
                before_unfolding_chi2.Form("%.2f", chi2_without_correlation);
                chi2_.DrawLatexNDC(0.2, 0.15, "#splitline{Norm. #chi_{Detector unfolded}^{2}: " + after_unfolding_chi2 + "}{Norm. #chi^{2}: " + before_unfolding_chi2 + "}");
            }
        }
        else{
            if(var == "Pt" ){
                after_unfolding_chi2.Form("%.2f", DoFit("Pt", nthMassBin, isFSRUnfold));
                before_unfolding_chi2.Form("%.2f", DoFit("Pt", nthMassBin));
                chi2_.DrawLatexNDC(0.2, 0.15, "#splitline{Norm. #chi_{FSR unfolded}^{2}: " + after_unfolding_chi2 + "}{Norm. #chi_{detector unfolded}^{2}: " + before_unfolding_chi2 + "}");
            }
            else{
                after_unfolding_chi2.Form("%.2f", DoFit("Mass", nthMassBin, isFSRUnfold));
                before_unfolding_chi2.Form("%.2f", DoFit("Mass", nthMassBin));
                chi2_.DrawLatexNDC(0.2, 0.15, "#splitline{Norm. #chi_{FSR unfolded}^{2}: " + after_unfolding_chi2 + "}{Norm. #chi_{detector unfolded}^{2}: " + before_unfolding_chi2 + "}");
            }
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////// systematics ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if(isFSRUnfold && systematic && (sysName != "PDFerror")){
             
            TH1 * hdatasys0_temp = NULL;
            TH1 * hdatasys1_temp = NULL;
            int systematic_variation_index = 0;
      
            // read histogram giving maximum variation on mean value 
            if(var == "Pt"){
                if(sysName != "QED_FSR"){
                    systematic_variation_index = meanPtErrIdx_sysdata_pre_fsr.at(nthMassBin)[sysName];
                    hdatasys0_temp = sysPtFSRUnfold[sysName].at(systematic_variation_index)->GetOutput("hunfolded_pt_systemp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
                }
                else{
                    hdatasys0_temp = sysPtFSRUnfold[sysName].at(0)->GetOutput("hunfolded_pt_systemp_fsr0",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
                    hdatasys1_temp = sysPtFSRUnfold[sysName].at(1)->GetOutput("hunfolded_pt_systemp_fsr1",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
                }
            }

            if(var == "Mass"){
                if(sysName != "QED_FSR"){
                    systematic_variation_index = meanMassErrIdx_sysdata_pre_fsr.at(nthMassBin)[sysName];
                    hdatasys0_temp = sysMassFSRUnfold[sysName].at(systematic_variation_index)->GetOutput("hunfolded_mass_systemp",0,0,"mass[UO];pt[UOC0]",kTRUE);
                    hdatasys0_temp->GetXaxis()->SetRange(hdatasys0_temp->GetXaxis()->FindBin(massBins[nthMassBin]+0.01), hdatasys0_temp->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));
                }
                else{
                    hdatasys0_temp = sysMassFSRUnfold[sysName].at(0)->GetOutput("hunfolded_mass_systemp_fsr0",0,0,"mass[UO];pt[UOC0]",kTRUE);
                    hdatasys1_temp = sysMassFSRUnfold[sysName].at(1)->GetOutput("hunfolded_mass_systemp_fsr1",0,0,"mass[UO];pt[UOC0]",kTRUE);
                    hdatasys0_temp->GetXaxis()->SetRange(hdatasys0_temp->GetXaxis()->FindBin(massBins[nthMassBin]+0.01), hdatasys0_temp->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));
                    hdatasys1_temp->GetXaxis()->SetRange(hdatasys1_temp->GetXaxis()->FindBin(massBins[nthMassBin]+0.01), hdatasys1_temp->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));
                }
            }

            TH1F * ratio0_temp = NULL;
            TH1F * ratio1_temp = NULL;

            if(sysName != "QED_FSR"){
                if(data_over_mc){
                    ratio0_temp = (TH1F*)hdatasys0_temp->Clone("ratio");
                    ratio0_temp->Divide(hpreFSR_mc);
                }
                else{
                    ratio0_temp = (TH1F*)hpreFSR_mc->Clone("ratio");
                    ratio0_temp->Divide(hdatasys0_temp);
                }
            }
            else{
                if(data_over_mc){
                    ratio0_temp = (TH1F*)hdatasys0_temp->Clone("ratio");
                    ratio0_temp->Divide(hpreFSR_mc);

                    ratio1_temp = (TH1F*)hdatasys1_temp->Clone("ratio");
                    ratio1_temp->Divide(hpreFSR_mc);
                }
                else{
                    ratio0_temp = (TH1F*)hpreFSR_mc->Clone("ratio");
                    ratio0_temp->Divide(hdatasys0_temp);

                    ratio1_temp = (TH1F*)hpreFSR_mc->Clone("ratio");
                    ratio1_temp->Divide(hdatasys1_temp);
                }
            }
           


            // loop over each bin of unfolded histogram
            for(int ibin = 1; ibin<hunfolded_sys_err->GetNbinsX()+1;ibin++){
                Double_t temp_err;
                Double_t temp_ratio_err;
                // get "envelope"
                // absolute difference between nominal unfolded histogram and systematic unfolded histogram
                if(sysName != "QED_FSR"){
                    temp_err =  fabs(hunfolded_data->GetBinContent(ibin) - hdatasys0_temp->GetBinContent(ibin));
                    temp_ratio_err = fabs(ratio->GetBinContent(ibin) - ratio0_temp->GetBinContent(ibin));
                }
                else{   
                    temp_err =  fabs(hdatasys0_temp->GetBinContent(ibin) - hdatasys1_temp->GetBinContent(ibin));
                    temp_ratio_err = fabs(ratio0_temp->GetBinContent(ibin) - ratio1_temp->GetBinContent(ibin));
                }

                hunfolded_sys_err->SetBinError(ibin, temp_err);
                ratio_sys_err->SetBinContent(ibin, 1.);
	        if(temp_ratio_err < 5.e-6) temp_ratio_err = 1.e-6;
                ratio_sys_err->SetBinError(ibin, temp_ratio_err);

            }// loop for bin contents

            delete ratio0_temp;
            delete ratio1_temp;
            delete hdatasys0_temp;
            delete hdatasys1_temp;

            // draw systematic envelope for systematic source 
            hunfolded_sys_err->SetLineColor(12);
            hunfolded_sys_err->SetFillColor(12);
            hunfolded_sys_err->SetLineWidth(5);
            hunfolded_sys_err->SetFillStyle(3005);
            hunfolded_sys_err->SetMarkerSize(0);
            hunfolded_sys_err->Draw("E2 same");

	    //delete hdata_sys_temp;
        }
        ////////////////////////////////////////////////////// systematic ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	c1->cd();

        TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.4);
        pad2->SetTopMargin(0.05);
        pad2->SetBottomMargin(0.2);
        pad2->SetTicks(1);
        pad2->SetGridy(1);
        pad2->Draw();

        pad2->cd();

        if(data_over_mc){
            ratio->Draw("pe");
            ratio_detector->Draw("pesame");
            //ratio_detector->SetMarkerColor(kGray+2);
            //ratio_detector->SetLineColor(kGray+2);
            ratio->GetYaxis()->SetTitle("Data/ MC");
        }
        else{
            ratio->Draw("hist");
            ratio_detector->Draw("histsame");
            ratio->SetLineColor(kRed);
            ratio_detector->SetLineColor(kMagenta);
            ratio->GetYaxis()->SetTitle("MC/ Data");

        }

        if(isFSRUnfold && systematic && (sysName != "PDFerror" || sysName != "QED_FSR")){
            ratio_sys_err->SetLineColor(12);
            ratio_sys_err->SetFillColor(12);
            ratio_sys_err->SetLineWidth(5);
            ratio_sys_err->SetFillStyle(3005);
            ratio_sys_err->SetMarkerSize(0);
            ratio_sys_err->Draw("E2 same");
            //ratio_sys_err->SetFillColorAlpha(kBlack,0.3);
        }
        if(var=="Pt") ratio->GetXaxis()->SetTitle("p_{T} (GeV)");
        if(var=="Mass") ratio->GetXaxis()->SetTitle("mass (GeV)");
        ratio->SetMinimum(0.5);
        ratio->SetMaximum(1.5);
        ratio->SetTitle("");
        //ratio->GetXaxis()->SetTitleOffset(1.5);
        ratio->GetYaxis()->SetNdivisions(505);

        TLine *l_;
        l_ = new TLine(ratio->GetXaxis()->GetXmin(),1,ratio->GetXaxis()->GetXmax(),1);
        if(var=="Mass") l_ = new TLine(massBins[nthMassBin],1,massBins[nthMassBin+1],1);
        l_->Draw("same");
        l_->SetLineStyle(1);
        if(data_over_mc) l_->SetLineColor(kRed);
        else l_->SetLineColor(kBlack);

	CMS_lumi( c1, 4, 0 );
        c1->cd();
        if(!isFSRUnfold) c1->SaveAs(outpdf+"_"+ibinMass+"_"+var+".pdf");
        else c1->SaveAs(outpdf+"_"+ibinMass+"_"+var+"_FSRUnfold.pdf");

        delete hunfolded_data;
	delete hunfolded_sys_err;
        delete hpreFSR_mc;
        delete pad1;
        delete pad2;
        delete c1;
}
