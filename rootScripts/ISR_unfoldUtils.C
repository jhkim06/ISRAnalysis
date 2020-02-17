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

void ISRUnfold::SetNomTUnfoldDensity(TString var, TString filepath, TString phase_name, TString fsr_correction_name, TString filepath_closure)
{

	TFile* filein = new TFile(filepath);
        // root file for closure test only
        // TODO
        // make an array with three TUnfoldDensity objects for three closure tests
	//TFile* filein_closure = new TFile("/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/output/2016/electron/ISR_response_closureTest.root");
        TFile* filein_closure = new TFile(filepath_closure);

        // set response matrix
        TH2* hmcGenRec;
        TH2* hmcGenRec_closure;

        if(var == "Pt")
        {
            hmcGenRec = (TH2*)filein->Get(phase_name + "/ptll_rec_gen_" + fsr_correction_name + "_response_matrix/hmc" + var + "GenRecnominal");
            hmcGenRec_closure = (TH2*)filein_closure->Get(phase_name + "_unfold_split_p5/ptll_rec_gen_" + fsr_correction_name + "_response_matrix/hmc" + var + "GenRecnominal");
        }
        else if(var == "Mass")
        {
            hmcGenRec = (TH2*)filein->Get(phase_name + "/mll_rec_gen_" + fsr_correction_name + "_response_matrix/hmc" + var + "GenRecnominal");
            hmcGenRec_closure = (TH2*)filein_closure->Get(phase_name + "_unfold_split_p5/mll_rec_gen_" + fsr_correction_name + "_response_matrix/hmc" + var + "GenRecnominal");
        }
        else{
            cout << "ISRUnfold::SetNomTUnfoldDensity, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
        }

        if(do_normalization)
        {
            if(channel_name =="electron") hmcGenRec->Scale(0.959939);
            else hmcGenRec->Scale(0.953026);
        }

        // set binning definition
        TUnfoldBinning* binning_Rec = NULL;
        TUnfoldBinning* binning_Gen = NULL;

        if( var == "Pt" )
        {

          TString Rec_Pt = "Rec_Pt";
          TString Gen_Pt = "Gen_Pt";

          Rec_Pt = phase_name + "/ptll_rec_gen_" + fsr_correction_name + "_response_matrix/" + Rec_Pt;
          Gen_Pt = phase_name + "/ptll_rec_gen_" + fsr_correction_name + "_response_matrix/" + Gen_Pt;

          binning_Rec = (TUnfoldBinning*)filein->Get(Rec_Pt);
          binning_Gen = (TUnfoldBinning*)filein->Get(Gen_Pt);
        }
        else if( var == "Mass" )
        {

          TString Rec_Mass = "Rec_Mass";
          TString Gen_Mass = "Gen_Mass";

          Rec_Mass = phase_name + "/mll_rec_gen_" + fsr_correction_name + "_response_matrix/" + Rec_Mass;
          Gen_Mass = phase_name + "/mll_rec_gen_" + fsr_correction_name + "_response_matrix/" + Gen_Mass;

          binning_Rec = (TUnfoldBinning*)filein->Get(Rec_Mass);
          binning_Gen = (TUnfoldBinning*)filein->Get(Gen_Mass);

        }
        else{
            cout << "ISRUnfold::SetNomTUnfoldDensity, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
        }

	if( var == "Pt" )
        {
        	nomPtUnfold = new TUnfoldDensityV17(hmcGenRec,
        	                               TUnfold::kHistMapOutputHoriz,
        	                               regMode_detector,
        	                               TUnfold::kEConstraintArea,
        	                               TUnfoldDensityV17::kDensityModeBinWidth,
        	                               binning_Gen,binning_Rec);

                // for closure test
                nomPtUnfold_closure = new TUnfoldDensityV17(hmcGenRec_closure,
                                               TUnfold::kHistMapOutputHoriz,
                                               regMode_detector,
                                               TUnfold::kEConstraintArea,
                                               TUnfoldDensityV17::kDensityModeBinWidth,
                                               binning_Gen,binning_Rec);
	}
        else if( var == "Mass" )
        {
                nomMassUnfold = new TUnfoldDensityV17(hmcGenRec,
                                               TUnfold::kHistMapOutputHoriz,
                                               regMode_detector,
                                               TUnfold::kEConstraintArea,
                                               TUnfoldDensityV17::kDensityModeBinWidth,
                                               binning_Gen,binning_Rec);

                // for closure
                nomMassUnfold_closure = new TUnfoldDensityV17(hmcGenRec_closure,
                                               TUnfold::kHistMapOutputHoriz,
                                               regMode_detector,
                                               TUnfold::kEConstraintArea,
                                               TUnfoldDensityV17::kDensityModeBinWidth,
                                               binning_Gen,binning_Rec);
        }
        else{
            cout << "ISRUnfold::SetNomTUnfoldDensity, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
        }

    filein->Close();
    delete filein;

}

void ISRUnfold::setNomFSRTUnfoldDensity(TString var, TString filepath, TString phase_name, TString fsr_correction_name,TString filepath_closure)
{

        TFile* filein = new TFile(filepath);
        TFile* filein_closure = new TFile(filepath_closure);

        // set response matrix
        TH2* hmcGenGen;
        TH2* hmcGenGen_closure;

        if(var == "Pt")
        {
            hmcGenGen = (TH2*)filein->Get(phase_name + "/ptll_gen_post_fsr_" + fsr_correction_name + "_response_matrix/hmc" + var + "GenRecnominal");
            hmcGenGen_closure = (TH2*)filein_closure->Get(phase_name + "_split_p5/ptll_gen_post_fsr_" + fsr_correction_name + "_response_matrix/hmc" + var + "GenRecnominal");
        }
        else if(var == "Mass")
        {
            hmcGenGen = (TH2*)filein->Get(phase_name + "/mll_gen_post_fsr_" + fsr_correction_name + "_response_matrix/hmc" + var + "GenRecnominal");
            hmcGenGen_closure = (TH2*)filein_closure->Get(phase_name + "_split_p5/mll_gen_post_fsr_" + fsr_correction_name + "_response_matrix/hmc" + var + "GenRecnominal");
        }
        else{
            cout << "ISRUnfold::SetNomTUnfoldDensity, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
        }

        if(do_normalization){
            if(channel_name =="electron") hmcGenGen->Scale(0.959939);
            else hmcGenGen->Scale(0.953026);
        }

        // set binning definition
        TUnfoldBinning* binning_Gen = NULL;

        if( var == "Pt" )
        {

          TString Gen_Pt = "Gen_Pt";

          Gen_Pt = phase_name + "/ptll_gen_post_fsr_" + fsr_correction_name + "_response_matrix/" + Gen_Pt;

          binning_Gen = (TUnfoldBinning*)filein->Get(Gen_Pt);
        }
        else if( var == "Mass" )
        {

          TString Gen_Mass = "Gen_Mass";

          Gen_Mass = phase_name + "/mll_gen_post_fsr_" + fsr_correction_name + "_response_matrix/" + Gen_Mass;

          binning_Gen = (TUnfoldBinning*)filein->Get(Gen_Mass);

        }
        else{
            cout << "ISRUnfold::SetNomTUnfoldDensity, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
        }

        if( var == "Pt" )
        {
                nomPtFSRUnfold = new TUnfoldDensityV17(hmcGenGen,
                                               TUnfold::kHistMapOutputHoriz,
                                               regMode_FSR, // fixed to use no regularisation temporary
                                               TUnfold::kEConstraintArea,
                                               TUnfoldDensityV17::kDensityModeBinWidth,
                                               binning_Gen,binning_Gen);

                // for closure test
                nomPtFSRUnfold_closure = new TUnfoldDensityV17(hmcGenGen_closure,
                                               TUnfold::kHistMapOutputHoriz,
                                               regMode_FSR,
                                               TUnfold::kEConstraintNone,
                                               TUnfoldDensityV17::kDensityModeBinWidth,
                                               binning_Gen,binning_Gen);
        }
        else if( var == "Mass" )
        {
                nomMassFSRUnfold = new TUnfoldDensityV17(hmcGenGen,
                                               TUnfold::kHistMapOutputHoriz,
                                               regMode_FSR,
                                               TUnfold::kEConstraintArea,
                                               TUnfoldDensityV17::kDensityModeBinWidth,
                                               binning_Gen,binning_Gen);

                // for closure test
                nomMassFSRUnfold_closure = new TUnfoldDensityV17(hmcGenGen_closure,
                                               TUnfold::kHistMapOutputHoriz,
                                               regMode_FSR,
                                               TUnfold::kEConstraintNone,
                                               TUnfoldDensityV17::kDensityModeBinWidth,
                                               binning_Gen,binning_Gen);
        }
        else{
            cout << "ISRUnfold::SetNomTUnfoldDensity, only Pt and Mass available for var" << endl;
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
            cout << "ISRUnfold::SetNomTUnfoldDensity, only Pt and Mass available for var" << endl;
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
            cout << "ISRUnfold::SetNomTUnfoldDensity, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
        }
    }

    if(do_normalization)
    {
        if(channel_name =="electron") hmcGenGen->Scale(0.959939);
        else hmcGenGen->Scale(0.953026);
    }

    // set binning definition
    TUnfoldBinning* binning_Gen = NULL;

    if( var == "Pt" )
    {

      TString Gen_Pt = "Gen_Pt";

      Gen_Pt = phase_name + "/ptll_gen_post_fsr_" + fsr_correction_name + "_response_matrix/" + Gen_Pt;

      binning_Gen = (TUnfoldBinning*)filein->Get(Gen_Pt);
    }
    else if( var == "Mass" )
    {

      TString Gen_Mass = "Gen_Mass";

      Gen_Mass = phase_name + "/mll_gen_post_fsr_" + fsr_correction_name + "_response_matrix/" + Gen_Mass;

      binning_Gen = (TUnfoldBinning*)filein->Get(Gen_Mass);

    }
    else{
        cout << "ISRUnfold::SetNomTUnfoldDensity, only Pt and Mass available for var" << endl;
        exit (EXIT_FAILURE);
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
                                           binning_Gen,binning_Gen));
    }
    else if( var == "Mass" )
    {
            sysMassFSRUnfold[sysName].push_back( new TUnfoldDensityV17(hmcGenGen,
                                           TUnfold::kHistMapOutputHoriz,
                                           mode,
                                           TUnfold::kEConstraintArea,
                                           TUnfoldDensityV17::kDensityModeBinWidth,
                                           binning_Gen,binning_Gen));
    }
    else
    {
        cout << "ISRUnfold::SetNomTUnfoldDensity, only Pt and Mass available for var" << endl;
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

    // set migration matrix
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

    if(do_normalization)
    {
        if(channel_name =="electron") hmcGenRec->Scale(0.959939);
        else hmcGenRec->Scale(0.953026);
    }

    // set bin definition
    TUnfoldBinning* binning_Rec = NULL;
    TUnfoldBinning* binning_Gen = NULL;

    if( var == "Pt" )
    {

      TString Rec_Pt = "Rec_Pt";
      TString Gen_Pt = "Gen_Pt";

      Rec_Pt = phase_name + "/ptll_rec_gen_" + fsr_correction_name + "_response_matrix/" + Rec_Pt;
      Gen_Pt = phase_name + "/ptll_rec_gen_" + fsr_correction_name + "_response_matrix/" + Gen_Pt;

      binning_Rec = (TUnfoldBinning*)filein->Get(Rec_Pt);
      binning_Gen = (TUnfoldBinning*)filein->Get(Gen_Pt);
    }
    else if( var == "Mass" )
    {

      TString Rec_Mass = "Rec_Mass";
      TString Gen_Mass = "Gen_Mass";

      Rec_Mass = phase_name + "/mll_rec_gen_" + fsr_correction_name + "_response_matrix/" + Rec_Mass;
      Gen_Mass = phase_name + "/mll_rec_gen_" + fsr_correction_name + "_response_matrix/" + Gen_Mass;

      binning_Rec = (TUnfoldBinning*)filein->Get(Rec_Mass);
      binning_Gen = (TUnfoldBinning*)filein->Get(Gen_Mass);

    }
    else
    {
        cout << "ISRUnfold::setSysTUnfoldDensity, only Pt and Mass available for var" << endl;
        exit (EXIT_FAILURE);
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
                                           binning_Gen,binning_Rec));
    }

    else if( var == "Mass" )
    {
            sysMassUnfold[sysName].push_back( new TUnfoldDensityV17(hmcGenRec,
                                           TUnfold::kHistMapOutputHoriz,
                                           mode,
                                           TUnfold::kEConstraintArea,
                                           TUnfoldDensityV17::kDensityModeBinWidth,
                                           binning_Gen,binning_Rec));
    }
    else
    {
        cout << "ISRUnfold::setSysTUnfoldDensity, only Pt and Mass available for var" << endl;
        exit (EXIT_FAILURE);
    }

    filein->Close();
    delete filein;
}

void ISRUnfold::drawISRMatrixInfo(TString var, TString outpdf, bool detector_unfold, bool fsr_systematic)
{
    const TUnfoldBinningV17* xaxis_binning;
    const TVectorD* xaxis1_tvecd;
    int xaxis1_nbin;

    const TUnfoldBinningV17* yaxis_binning = NULL;
    const TVectorD* yaxis1_tvecd = NULL;
    const TVectorD* yaxis2_tvecd = NULL;
    int yaxis1_nbin;


    if(var=="Pt"){
        if(detector_unfold){
            xaxis_binning = nomPtUnfold->GetOutputBinning("Gen_Pt");
            xaxis1_tvecd = xaxis_binning->GetDistributionBinning(0);
            xaxis1_nbin = xaxis1_tvecd->GetNrows() - 1;

            yaxis_binning = nomPtUnfold->GetInputBinning("Rec_Pt");
            yaxis1_tvecd = yaxis_binning->GetDistributionBinning(0);
            yaxis2_tvecd = yaxis_binning->GetDistributionBinning(1);
            yaxis1_nbin = yaxis1_tvecd->GetNrows() - 1;
        }
        else{
            xaxis_binning = nomPtFSRUnfold->GetOutputBinning("Gen_Pt");
            xaxis1_tvecd = xaxis_binning->GetDistributionBinning(0);
            xaxis1_nbin = xaxis1_tvecd->GetNrows() - 1;

            yaxis_binning = nomPtFSRUnfold->GetInputBinning("Gen_Pt");
            yaxis1_tvecd = yaxis_binning->GetDistributionBinning(0);
            yaxis2_tvecd = yaxis_binning->GetDistributionBinning(1);
            yaxis1_nbin = yaxis1_tvecd->GetNrows() - 1;
        }
    }
    else if(var=="Mass"){
        if(detector_unfold){
            xaxis_binning = nomMassUnfold->GetOutputBinning("Gen_Mass");
            xaxis1_tvecd = xaxis_binning->GetDistributionBinning(0);
            xaxis1_nbin = xaxis1_tvecd->GetNrows() - 1;

            yaxis_binning = nomMassUnfold->GetInputBinning("Rec_Mass");
            yaxis1_tvecd = yaxis_binning->GetDistributionBinning(0);
            yaxis2_tvecd = yaxis_binning->GetDistributionBinning(1);
            yaxis1_nbin = yaxis1_tvecd->GetNrows() - 1;
        }
        else{
            xaxis_binning = nomMassFSRUnfold->GetOutputBinning("Gen_Mass");
            xaxis1_tvecd = xaxis_binning->GetDistributionBinning(0);
            xaxis1_nbin = xaxis1_tvecd->GetNrows() - 1;

            yaxis_binning = nomMassFSRUnfold->GetInputBinning("Gen_Mass");
            yaxis1_tvecd = yaxis_binning->GetDistributionBinning(0);
            yaxis2_tvecd = yaxis_binning->GetDistributionBinning(1);
            yaxis1_nbin = yaxis1_tvecd->GetNrows() - 1;
        }
    }
    else{
            cout << "ISRUnfold::drawISRMatrixInfo, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
    }

    setTDRStyle();
    writeExtraText = true;
    extraText  = "work in progress";

    c1 = new TCanvas("c1","c1", 50, 50, 800, 700);
    gStyle->SetOptFit(0);
    gStyle->SetPalette(55);
    c1->cd();

    c1->SetBottomMargin(0.2);
    c1->SetRightMargin(0.15);
    c1->SetLeftMargin(0.2);
    c1->SetTopMargin(0.08);
    c1->SetTicks(1);
    c1->SetLogz();

    TH2 *histProb = NULL;
    TH1 *histEfficiency = NULL;

    if(var=="Pt")
    {
        if(detector_unfold)
        {
            histProb = nomPtUnfold->GetProbabilityMatrix("Migration prob. for pt mass bin",";p_T(gen);p_T(Rec)");
        }
        else
        {
            if(!fsr_systematic) histProb = nomPtFSRUnfold->GetProbabilityMatrix("Migration prob. for pt mass bin",";p_T(gen);p_T(Rec)");
            else histProb = sysPtFSRUnfold["QED_FSR"].at(0)->GetProbabilityMatrix("Migration prob. for pt mass bin",";p_T(gen);p_T(Rec)");
        }
    }
    else if(var=="Mass")
    {
        if(detector_unfold)
        {
            histProb = nomMassUnfold->GetProbabilityMatrix("Migration prob. for mass bin",";mass(gen);mass(Rec)");
        }
        else
        {
            if(!fsr_systematic) histProb = nomMassFSRUnfold->GetProbabilityMatrix("Migration prob. for mass bin",";mass(gen);mass(Rec)");
            else histProb = sysMassFSRUnfold["QED_FSR"].at(0)->GetProbabilityMatrix("Migration prob. for mass bin",";mass(gen);mass(Rec)");
        }
    }
    else{
            cout << "ISRUnfold::drawISRMatrixInfo, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
    }

    histEfficiency=histProb->ProjectionX("histEfficiency");

    TH2D* ticks_ = NULL;
    TString draw_option = "COLZ";

    //if(var == "Pt")
    {
        ticks_ = new TH2D("tick", "tick", histProb->GetNbinsX(), histProb->GetXaxis()->GetXmin(), histProb->GetXaxis()->GetXmax(),
                                          histProb->GetNbinsY(), histProb->GetYaxis()->GetXmin(), histProb->GetYaxis()->GetXmax());

        const TVectorD* temp_tvecd = xaxis_binning->GetDistributionBinning(1);
        const Double_t* massBins = temp_tvecd->GetMatrixArray();
        int nMassBin = temp_tvecd->GetNrows() - 1;

        int center = (xaxis1_nbin + xaxis_binning->HasUnderflow(0) + xaxis_binning->HasOverflow(0)) / 2;
        int totalBins = xaxis1_nbin + xaxis_binning->HasUnderflow(0) + xaxis_binning->HasOverflow(0);
        int center_y = (yaxis1_nbin + yaxis_binning->HasUnderflow(0) + yaxis_binning->HasOverflow(0)) / 2;
        int totalBins_y = yaxis1_nbin + yaxis_binning->HasUnderflow(0) + yaxis_binning->HasOverflow(0);

        if(xaxis_binning->HasUnderflow(1))
        {
            ticks_->GetXaxis()->SetBinLabel(histProb->GetXaxis()->FindBin(center), "Underflow");  // 5 = # of pt bins / 2  
            ticks_->GetYaxis()->SetBinLabel(histProb->GetYaxis()->FindBin(center_y), "Underflow");  // 5 = # of pt bins / 2  

            center += totalBins;
            center_y += totalBins_y;
        }

        for(int ibin = 0; ibin < nMassBin; ibin++)
        {
            TString lowMassEdge;
            TString highMassEdge;

            lowMassEdge.Form("%d", (int)massBins[ibin]);
            highMassEdge.Form("%d", (int)massBins[ibin+1]);
            TString var_name = "M";
            if(var == "Mass") var_name = "p_{T}";
            ticks_->GetXaxis()->SetBinLabel(histProb->GetXaxis()->FindBin(center), lowMassEdge+"<"+var_name+"<"+highMassEdge+" GeV");  // TODO option to set mass region 

            ticks_->GetYaxis()->SetBinLabel(histProb->GetYaxis()->FindBin(center_y), lowMassEdge+"<"+var_name+"<"+highMassEdge+" GeV");   

            center += totalBins;
            center_y += totalBins_y;
        }
        if(xaxis_binning->HasOverflow(1))
        {
            ticks_->GetXaxis()->SetBinLabel(histProb->GetXaxis()->FindBin(center), "Overflow");   
            ticks_->GetYaxis()->SetBinLabel(histProb->GetYaxis()->FindBin(center_y), "Overflow");   
        }

        if(detector_unfold) ticks_->GetYaxis()->SetTitle("Detector level bin");
        else ticks_->GetYaxis()->SetTitle("Dressed level bin");

        ticks_->GetYaxis()->SetTitleFont(43);
        ticks_->GetYaxis()->SetTitleSize(40);
        ticks_->GetYaxis()->SetTitleOffset(1.2);
        if(detector_unfold) ticks_->GetXaxis()->SetTitle("Dressed level bin");
        else ticks_->GetXaxis()->SetTitle("Pre FSR level bin");
        ticks_->GetXaxis()->SetTitleFont(43);
        ticks_->GetXaxis()->SetTitleSize(40);
        ticks_->GetXaxis()->SetTitleOffset(1.2);

        ticks_->GetXaxis()->SetLabelColor(kGray+2);
        ticks_->GetYaxis()->SetLabelColor(kGray+2);
        ticks_->GetXaxis()->SetLabelSize(15);
        ticks_->GetYaxis()->SetLabelSize(15);
        ticks_->GetXaxis()->LabelsOption("v");
        ticks_->GetXaxis()->SetTickSize(0);
        ticks_->GetYaxis()->SetTickSize(0);
        ticks_->GetZaxis()->SetRangeUser(1e-3, 0.9);
        ticks_->Draw("colz"); // 
        histProb->Draw("same axis COLZ"); // "colz same" not working
        draw_option = "same COLZ";
    }

    histProb->Draw(draw_option); 
    histProb->GetZaxis()->SetRangeUser(1e-3, 0.9);
    histProb->GetYaxis()->SetTitle("Detector bin");
    histProb->GetYaxis()->SetTitleFont(43);
    histProb->GetYaxis()->SetTitleSize(40);
    histProb->GetYaxis()->SetTitleOffset(1.2);
    histProb->GetXaxis()->SetTitle("Generator bin");
    histProb->GetXaxis()->SetTitleFont(43);
    histProb->GetXaxis()->SetTitleSize(40);
    histProb->GetXaxis()->SetTitleOffset(1.2);

    c1->RedrawAxis();
    
    if(var=="Pt")
    {
        if(detector_unfold) histProb->SetTitle("migration probabilities;p_{T} mass bin index (post FSR) ;p_{T} mass bin index (Rec)");
        else histProb->SetTitle("migration probabilities;p_{T} mass bin index (pre FSR) ;p_{T} mass bin index (post FSR)");
    }
    else if(var=="Mass")
    {
        if(detector_unfold) histProb->SetTitle("migration probabilities;mass bin index (post FSR) ;mass bin index (Rec)");
        else histProb->SetTitle("migration probabilities;mass bin index (pre FSR) ;mass bin index (post FSR)");
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
    grid_.SetLineColor(kGray+2);
    //grid_.SetLineColorAlpha(kGray+2, 0.35);;
    grid_.SetLineStyle(2);
    grid_bin_boundary.SetLineColorAlpha(kRed, 1.);;
    grid_bin_boundary.SetLineStyle(kSolid);


    //int boundarybin_x = xaxis1_nbin;
    int boundarybin_x = 1;
    int boundarybin_y = yaxis1_nbin + 1;
    if(xaxis_binning->HasUnderflow(0)) boundarybin_x++;
    if(yaxis_binning->HasUnderflow(0)) boundarybin_y++;
    int countDrawnRedBoundary = 0;

    for( int ii=0; ii<histProb->GetXaxis()->GetNbins(); ii++ )
    {
        Int_t i_bin = ii+1;
        Double_t binEdge = histProb->GetXaxis()->GetBinLowEdge(i_bin);

        if(boundarybin_x == i_bin)
        {
         
            grid_bin_boundary.DrawLine(binEdge, histProb->GetYaxis()->GetBinLowEdge(boundarybin_y - yaxis1_nbin), binEdge, histProb->GetYaxis()->GetBinLowEdge(boundarybin_y) );
            countDrawnRedBoundary++;
            if(countDrawnRedBoundary%2 == 0)
            {
                if(!xaxis_binning->HasOverflow(0))
                    grid_.DrawLine(binEdge, histProb->GetYaxis()->GetBinUpEdge(0), binEdge, histProb->GetYaxis()->GetBinUpEdge(histProb->GetYaxis()->GetNbins()) );
                else
                    grid_.DrawLine(binEdge + 1, histProb->GetYaxis()->GetBinUpEdge(0), binEdge + 1, histProb->GetYaxis()->GetBinUpEdge(histProb->GetYaxis()->GetNbins()) );

                boundarybin_y += yaxis1_nbin;
                boundarybin_y ++;

                if(yaxis_binning->HasUnderflow(0)) boundarybin_y++;
                //if(yaxis_binning->HasOverflow(0)) boundarybin_y++;

                if(xaxis_binning->HasUnderflow(0)) boundarybin_x++;
                if(xaxis_binning->HasOverflow(0)) boundarybin_x++;
            }
            else
            {
                boundarybin_x += xaxis1_nbin; // next edge to draw
            }
        }
    }

    boundarybin_x = xaxis1_nbin + 1;
    boundarybin_y = 1;
    countDrawnRedBoundary = 0;
    if(xaxis_binning->HasUnderflow(0)) boundarybin_x++;
    if(yaxis_binning->HasUnderflow(0)) boundarybin_y++;

    for( int ii=0; ii<histProb->GetYaxis()->GetNbins(); ii++ )
    {
        Int_t i_bin = ii+1;
        Double_t binEdge = histProb->GetYaxis()->GetBinLowEdge(i_bin);

        if(boundarybin_y == i_bin)
        {
            grid_bin_boundary.DrawLine(histProb->GetXaxis()->GetBinLowEdge(boundarybin_x - xaxis1_nbin), binEdge, histProb->GetXaxis()->GetBinLowEdge(boundarybin_x),binEdge);
            countDrawnRedBoundary++;
            if(countDrawnRedBoundary%2 == 0)
            {
                if(!yaxis_binning->HasOverflow(0))
                    grid_.DrawLine(histProb->GetXaxis()->GetBinUpEdge(0), binEdge, histProb->GetXaxis()->GetBinUpEdge(histProb->GetXaxis()->GetNbins()),binEdge);
                else
                    grid_.DrawLine(histProb->GetXaxis()->GetBinUpEdge(0), binEdge+1, histProb->GetXaxis()->GetBinUpEdge(histProb->GetXaxis()->GetNbins()),binEdge+1);

                boundarybin_x += xaxis1_nbin;
                boundarybin_x ++;
                if(xaxis_binning->HasUnderflow(0)) boundarybin_x++;
                //if(xaxis_binning->HasOverflow(0)) boundarybin_x++;

                if(yaxis_binning->HasUnderflow(0)) boundarybin_y++;
                if(yaxis_binning->HasOverflow(0)) boundarybin_y++;
            }
            else
            {
                boundarybin_y += yaxis1_nbin; // next edge to draw
            }
        }

    }


    CMS_lumi(c1, 6, 11);

    TLatex mcName;
    mcName.SetTextSize(0.04);
    if(!fsr_systematic)
    {
        mcName.DrawLatexNDC(c1->GetLeftMargin(), 1-(c1->GetTopMargin()) + 0.2 * c1->GetTopMargin(), "Response Matrix (aMC@NLO+PYTHIA)");
    }
    else mcName.DrawLatexNDC(c1->GetLeftMargin(), 1-(c1->GetTopMargin()) + 0.2 * c1->GetTopMargin(), "Response Matrix (Powheg+PHOTOS)");

    if(var=="Pt")
    {
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

    if(var=="Pt")
    {
        histEfficiency->SetTitle("efficiency;p_{T} mass bin (gen) ;A #times #epsilon");
    }
    if(var=="Mass")
    {
        histEfficiency->SetTitle("efficiency;mass bin (gen) ;A #times #epsilon");
    }
    histEfficiency->GetYaxis()->SetTitleOffset(1.5);
    histEfficiency->GetYaxis()->SetRangeUser(0., 1.2);
    histEfficiency->Draw();

    TLine grid_bin_boundary_;
    grid_bin_boundary_.SetLineColorAlpha(kRed, 1.);;
    grid_bin_boundary_.SetLineStyle(kSolid);

    boundarybin_x = xaxis1_nbin;
    boundarybin_y = yaxis1_nbin;
    if(xaxis_binning->HasUnderflow(0)) boundarybin_x++;
    if(yaxis_binning->HasUnderflow(0)) boundarybin_y++;
    int count_drawn_boundary = 0;
    bool add_UO = true;

    // draw y lines
    for( int ii=0; ii<histEfficiency->GetXaxis()->GetNbins(); ii++ )
    {
        Int_t i_bin = ii+1;
        Double_t binEdge = histEfficiency->GetXaxis()->GetBinUpEdge(i_bin);

        if(boundarybin_x == i_bin && i_bin != histEfficiency->GetXaxis()->GetNbins())
        {
            if(count_drawn_boundary == 0)
            {
                grid_bin_boundary_.DrawLine(histEfficiency->GetXaxis()->GetBinUpEdge(i_bin-xaxis1_nbin), 0., histEfficiency->GetXaxis()->GetBinUpEdge(i_bin-xaxis1_nbin), 1.2 );
            }

            if(add_UO)
            {
                if(xaxis_binning->HasUnderflow(0)) boundarybin_x++;
                if(xaxis_binning->HasOverflow(0)) boundarybin_x++;
                add_UO = false;
            }
            else
            {
                boundarybin_x += xaxis1_nbin;
                add_UO = true;
            }

            grid_bin_boundary_.DrawLine(binEdge, 0., binEdge, 1.2 );
            if(count_drawn_boundary%2 == 0){
                if(yaxis_binning->HasUnderflow(0)) boundarybin_y++;
                if(yaxis_binning->HasOverflow(0)) boundarybin_y++;
                boundarybin_y += yaxis1_nbin;
            }
            count_drawn_boundary++;
        }
    }

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
    if(ticks_ != NULL);
        delete ticks_;
    delete histProb;
    delete histEfficiency;
    delete c1;
}

void ISRUnfold::setFSRUnfoldInput(TString filepath, bool isSys, TString sysName, int nth, TString phase_name)
{

    TFile* filein = new TFile(filepath);

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

    filein->Close();
    delete filein;
}


void ISRUnfold::setInput(TString var, TString filepath, bool isSys, TString sysName, int nth, double bias, TString phase_name)
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
            if(channel_name == "muon")     hRec = (TH1*)filein->Get(phase_name + "/hist_ptll/histo_DoubleMuonnominal");
            if(channel_name == "electron" && year != 2018) hRec = (TH1*)filein->Get(phase_name + "/hist_ptll/histo_DoubleEGnominal");
            if(channel_name == "electron" && year == 2018) hRec = (TH1*)filein->Get(phase_name + "/hist_ptll/histo_EGammanominal");

            nomPtUnfold->SetInput(hRec,   bias);
        }
        else if(var == "Mass")
        {
            if(channel_name == "muon")     hRec = (TH1*)filein->Get(phase_name + "/hist_mll/histo_DoubleMuonnominal");
            if(channel_name == "electron" && year != 2018) hRec = (TH1*)filein->Get(phase_name + "/hist_mll/histo_DoubleEGnominal");
            if(channel_name == "electron" && year == 2018) hRec = (TH1*)filein->Get(phase_name + "/hist_mll/histo_EGammanominal");
            nomMassUnfold->SetInput(hRec, bias);
        }
        else{
            cout << "ISRUnfold::setInput, only Pt and Mass available for var" << endl;
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
                    cout << "ISRUnfold::setInput, only Pt and Mass available for var" << endl;
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
                    cout << "ISRUnfold::setInput, only Pt and Mass available for var" << endl;
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
                cout << "ISRUnfold::setInput, only Pt and Mass available for var" << endl;
                exit (EXIT_FAILURE);
            }
        }
    }
    filein->Close();
    delete filein;
}

void ISRUnfold::doClosureTest(int detOrFSR_unfold, TString filepath, TString dir_name)
{

    nScan = 50;
    nScan_mass = 50;
    TFile* filein = new TFile(filepath);
    TH1* hpt_DYMC = NULL;
    TH1* hmass_DYMC = NULL;

    // use DY MC as input for closure test
    // for FSR closure test, use dressed level distribution as unfolding input

    // detector unfold
    if(detOrFSR_unfold == 0)
    {
        hpt_DYMC = (TH1*)filein->Get(dir_name + "_unfold_split_p5/hist_ptll/histo_DYJetsnominal");
        hpt_DYMC->Add((TH1*)filein->Get(dir_name + "_unfold_split_p5/hist_ptll/histo_DYJets10to50nominal"));
        hpt_DYMC->Add((TH1*)filein->Get(dir_name + "_DY_Fake_unfold_split_p5/hist_ptll/histo_DYJetsnominal"), -1);
        hpt_DYMC->Add((TH1*)filein->Get(dir_name + "_DY_Fake_unfold_split_p5/hist_ptll/histo_DYJets10to50nominal"), -1);
        nomPtUnfold_closure->SetInput(hpt_DYMC,   nominal_bias);

        hmass_DYMC = (TH1*)filein->Get(dir_name + "_unfold_split_p5/hist_mll/histo_DYJetsnominal");
        hmass_DYMC->Add((TH1*)filein->Get(dir_name + "_unfold_split_p5/hist_mll/histo_DYJets10to50nominal"));
        hmass_DYMC->Add((TH1*)filein->Get(dir_name + "_DY_Fake_unfold_split_p5/hist_mll/histo_DYJetsnominal"), -1);
        hmass_DYMC->Add((TH1*)filein->Get(dir_name + "_DY_Fake_unfold_split_p5/hist_mll/histo_DYJets10to50nominal"), -1);
        nomMassUnfold_closure->SetInput(hmass_DYMC,   nominal_bias);

        if(regMode_detector == TUnfold::kRegModeNone)
        {
            // no regularisation, set tau as zero
            nomPtUnfold_closure->DoUnfold(0);
            nomMassUnfold_closure->DoUnfold(0);
        }
        else
        {
            // regularization, use ScanTau() as a default method to find tau
            iBest=nomPtUnfold_closure->ScanTau(nScan,0.,0.,&rhoLogTau,
                                       TUnfoldDensity::kEScanTauRhoAvgSys,
                                       0,0,
                                       &lCurve);

            iBest_mass=nomMassUnfold_closure->ScanTau(nScan_mass,0.,0.,&rhoLogTau_mass,
                                       TUnfoldDensity::kEScanTauRhoAvgSys,
                                       0,0,
                                       &lCurve_mass);
        }
    }
    else if(detOrFSR_unfold == 1)
    // FSR unfolding
    {
        hpt_DYMC = (TH1*)filein->Get(dir_name + "_split_p5/hist_ptll/histo_DYJetsnominal");
        hpt_DYMC->Add((TH1*)filein->Get(dir_name + "_split_p5/hist_ptll/histo_DYJets10to50nominal"));
        nomPtFSRUnfold_closure->SetInput(hpt_DYMC,   nominal_bias);

        hmass_DYMC = (TH1*)filein->Get(dir_name + "_split_p5/hist_mll/histo_DYJetsnominal");
        hmass_DYMC->Add((TH1*)filein->Get(dir_name + "_split_p5/hist_mll/histo_DYJets10to50nominal"));
        nomMassFSRUnfold_closure->SetInput(hmass_DYMC,   nominal_bias);

        if(regMode_FSR == TUnfold::kRegModeNone)
        {
            // no regularisation, set tau as zero
            nomPtFSRUnfold_closure->DoUnfold(0);
            nomMassFSRUnfold_closure->DoUnfold(0);
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
}

void ISRUnfold::subBkgs(TString var, TString filepath, TString bkgName, bool isSys, TString sysName, int totSysN, int nth, TString phase_name)
{

	TFile* filein = new TFile(filepath);
        TH1* hRec = NULL;

        double bkg_scale = 1.;
        if(do_normalization){
            if(channel_name=="electron") bkg_scale = 0.959939;
            if(channel_name=="muon") bkg_scale = 0.953026;
        }

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
                    cout << "histo: " << phase_name + "/hist_ptll/histo_" + bkgName + "nominal" << endl;
                    hRec = (TH1*)filein->Get(phase_name + "/hist_ptll/histo_" + bkgName + "nominal");
                }
                if(var == "Mass")
                {
                    hRec = (TH1*)filein->Get(phase_name + "/hist_mll/histo_" + bkgName + "nominal");
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

void ISRUnfold::drawLCurve(TString outpdf, TString var)
{
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

void ISRUnfold::drawRhoLog(TString outpdf, TString var)
{
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
            // no regularisation, set tau as zero
            const TUnfoldBinningV17* bin=nomPtUnfold->GetOutputBinning();
            int istart=bin->GetGlobalBinNumber(0.1, 200.1);
            int iend=bin->GetGlobalBinNumber(100-0.01,200.1);
            nomPtUnfold->RegularizeBins(istart,1,iend-istart+1,TUnfoldV17::kRegModeCurvature);
            //nomPtUnfold->DoUnfold(1e-3);
            iBest=nomPtUnfold->ScanTau(nScan,0.,0.,&rhoLogTau,
                                       TUnfoldDensity::kEScanTauRhoAvgSys,
                                       0,0,
                                       &lCurve);
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

void ISRUnfold::setMeanMass(bool doSys, bool altMC, bool detector_unfold)
{

    int nMassBin = -1;

    // get mass bin definition for post FSR level from TUnfoldDensity object
    // nothe currently post and pre FSR level use the same bin definition
    const TUnfoldBinningV17* temp_binning_gen_pt = nomPtUnfold->GetOutputBinning("Gen_Pt");
    const TVectorD* temp_tvecd = temp_binning_gen_pt->GetDistributionBinning(1);
    nMassBin = temp_tvecd->GetNrows() - 1;
    const Double_t* massBins = temp_tvecd->GetMatrixArray();

    // get detector level distribution
    TH1 * hdetector_mass = nomMassUnfold->GetInput("hdetector_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);

    TFile* filein = new TFile(hist_file_path_DYHists);
    TString genHistPath;

    TString massLowEdge[5] = {"50.", "65.", "80.", "100.", "200."};
    TString massHighEdge[5] = {"65.", "80.", "100.", "200.", "350."};

    if(channel_name=="muon")
    {
        massLowEdge[0] = "40.";
        massLowEdge[1] = "60.";

        massHighEdge[0] = "60.";
        massHighEdge[1] = "80.";
    }

    // get nominal detector & QED FSR unfolded results
    TH1* hunfolded_mass =  nomMassUnfold->GetOutput("hunfolded_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);
    //TH1 *histMC_mass= nomMassUnfold->GetBias("histMC_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);
    TH1 *histMC_mass;
    TH1* hFSRunfolded_mass =  nomMassFSRUnfold->GetOutput("h_fsr_unfolded_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);
    //TH1 *histMC_pre_fsr_mass= nomMassFSRUnfold->GetBias("histMC_pre_fsr_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);
    TH1 *histMC_pre_fsr_mass;

    // loop over mass bins
    for(int ibin = 0; ibin < nMassBin; ibin++)
    {
        genHistPath = channel_name + "_full_phase_m" + massLowEdge[ibin] + "to" + massHighEdge[ibin] + "/mll_preFSR_" + channel_name + "/";
        histMC_pre_fsr_mass = (TH1*)filein->Get(genHistPath + "histo_DYJets");
        histMC_pre_fsr_mass->SetDirectory(0);
        histMC_pre_fsr_mass->Add((TH1*)filein->Get(genHistPath + "histo_DYJets10to50"));

        genHistPath = channel_name + "_fiducial_phase_post_FSR_dRp1_m" + massLowEdge[ibin] + "to" + massHighEdge[ibin] + "/mll_" + channel_name + "/";
        histMC_mass = (TH1*)filein->Get(genHistPath + "histo_DYJets");
        histMC_mass->SetDirectory(0);
        histMC_mass->Add((TH1*)filein->Get(genHistPath + "histo_DYJets10to50"));

        // set x-axis range
        hdetector_mass->GetXaxis()->  SetRange(hdetector_mass->GetXaxis()->  FindBin(massBins[ibin]+0.01),hdetector_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
        hunfolded_mass->GetXaxis()->  SetRange(hunfolded_mass->GetXaxis()->  FindBin(massBins[ibin]+0.01),hunfolded_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
        histMC_mass->GetXaxis()->SetRange(histMC_mass->GetXaxis()->FindBin(massBins[ibin]+0.01),histMC_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
        hFSRunfolded_mass->GetXaxis()->  SetRange(hunfolded_mass->GetXaxis()->  FindBin(massBins[ibin]+0.01),hunfolded_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
        histMC_pre_fsr_mass->GetXaxis()->SetRange(histMC_mass->GetXaxis()->FindBin(massBins[ibin]+0.01),histMC_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));

        // get mean values
        meanMass_data_detector. push_back(hdetector_mass->GetMean());
        meanMassStatErr_data_detector.push_back(hdetector_mass->GetMeanError());

        meanMass_data_det_unf.   push_back(hunfolded_mass->GetMean());
        meanMassStatErr_data_det_unf.push_back(hunfolded_mass->GetMeanError());

        meanMass_data_pre_fsr.   push_back(hFSRunfolded_mass->GetMean());
        meanMassStatErr_data_pre_fsr.push_back(hFSRunfolded_mass->GetMeanError());
        meanMassStatRelErr_data_pre_fsr.push_back(hFSRunfolded_mass->GetMeanError()/ hFSRunfolded_mass->GetMean() * 100.);

        meanMass_mc_det_unf.   push_back(histMC_mass->GetMean());
        meanMassErr_mc_det_unf.push_back(histMC_mass->GetMeanError()); //FIXME change to meanMassStatErr_mc_det_unf

        meanMass_mc_pre_fsr.   push_back(histMC_pre_fsr_mass->GetMean());

        if(altMC)
        {
            TH1 *histMC_massAlt= sysMassUnfold["Alt"].at(0)->GetBias("histMC_massAlt",0,0,"mass[UO];pt[UOC0]",kTRUE);
            histMC_massAlt->GetXaxis()->SetRange(histMC_mass->GetXaxis()->FindBin(massBins[ibin]+0.01),histMC_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
            meanMass_mcAlt.   push_back(histMC_massAlt->GetMean());
            meanMassErr_mcAlt.push_back(histMC_massAlt->GetMeanError());

            delete histMC_massAlt;
        }

        // save systematic mean mass values
        if(doSys)
        {
            // detector unfolding systematics
            std::map<TString, std::vector<TUnfoldDensityV17*>>::iterator it = sysMassUnfold.begin();
            std::map<TString, std::vector<Double_t>> temp_map; // temp map to save systematic results for a mass bin
            std::map<TString, std::vector<Double_t>> temp_map_mc; // temp map to save systematic results for a mass bin

       	    while(it != sysMassUnfold.end())
            {

       	        int nSys = it->second.size();
       	        TH1* hdatasys_temp;
       	        for(int i = 0; i < nSys; i++)
                {

       	            hdatasys_temp = sysMassUnfold[it->first].at(i)->GetOutput("hunfolded_mass_systemp",0,0,"mass[UO];pt[UOC0]",kTRUE);
    	            hdatasys_temp->GetXaxis()->SetRange(hdatasys_temp->GetXaxis()->FindBin(massBins[ibin]+0.01), hdatasys_temp->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
       	            temp_map[it->first].push_back(hdatasys_temp->GetMean());

                    
       	            delete hdatasys_temp;

                    hdatasys_temp = sysMassUnfold[it->first].at(i)->GetBias("hunfolded_mass_systemp",0,0,"mass[UO];pt[UOC0]",kTRUE);
                    hdatasys_temp->GetXaxis()->SetRange(hdatasys_temp->GetXaxis()->FindBin(massBins[ibin]+0.01), hdatasys_temp->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
                    temp_map_mc[it->first].push_back(hdatasys_temp->GetMean());

                    delete hdatasys_temp;
       	        }
       	        it++;
       	    }
       	    meanMass_sysdata_det_unf.push_back(temp_map);
       	    meanMass_sysmc_det_unf.push_back(temp_map_mc);
       	    temp_map.clear();
       	    temp_map_mc.clear();

            // QED FSR unfolding systematics
            it = sysMassFSRUnfold.begin();
            while(it != sysMassFSRUnfold.end())
            {
                int nSys = it->second.size();
                TH1* hdatasys_temp;
                for(int i = 0; i < nSys; i++)
                {
                    hdatasys_temp = sysMassFSRUnfold[it->first].at(i)->GetOutput("hunfolded_mass_systemp",0,0,"mass[UO];pt[UOC0]",kTRUE);
                    hdatasys_temp->GetXaxis()->SetRange(hdatasys_temp->GetXaxis()->FindBin(massBins[ibin]+0.01), hdatasys_temp->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
                    temp_map[it->first].push_back(hdatasys_temp->GetMean());

                    delete hdatasys_temp;

                    hdatasys_temp = sysMassFSRUnfold[it->first].at(i)->GetBias("hunfolded_mass_systemp",0,0,"mass[UO];pt[UOC0]",kTRUE);
                    hdatasys_temp->GetXaxis()->SetRange(hdatasys_temp->GetXaxis()->FindBin(massBins[ibin]+0.01), hdatasys_temp->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
                    temp_map_mc[it->first].push_back(hdatasys_temp->GetMean());
                    //cout << it->first << endl;
                    //cout << "mean: " << hdatasys_temp->GetMean() << endl;

                    delete hdatasys_temp;
                }
                it++;
            }
            meanMass_sysdata_pre_fsr.push_back(temp_map);
            meanMass_sysmc_pre_fsr.push_back(temp_map_mc);
            temp_map.clear();
            temp_map_mc.clear();
        }
    }// end of mass bin loop

    delete hdetector_mass;
    delete hunfolded_mass;
    delete histMC_mass;
    delete hFSRunfolded_mass;
    delete histMC_pre_fsr_mass;

    // calculate systematic uncertainty
    if(doSys)
    {
        for(int i = 0; i < nMassBin; i++)
        {
            // detector unfolding
            std::map<TString, Double_t> temp_map_;
            std::map<TString, Double_t> temp_map_mc_;
            std::map<TString, Int_t> temp_map_sysname_index;
            std::map<TString, Int_t> temp_map_mc_sysname_index;

            std::map<TString, std::vector<Double_t>>::iterator it = meanMass_sysdata_det_unf.at(i).begin();
            while(it != meanMass_sysdata_det_unf.at(i).end())
            {

                int size_ = it->second.size();

                TH1F *hpdfsys = NULL;
                TH1F *hpdfsys_mc = NULL;
                if((it->first)=="PDFerror" || (it->first)=="Stat")
                {
                    hpdfsys = new TH1F("pdfsys", "pdfsys", 100, meanMass_data_det_unf.at(i)-0.2, meanMass_data_det_unf.at(i)+0.2); // temp histogram to contain PDF variations
                    hpdfsys_mc = new TH1F("pdfsys_mc", "pdfsys_mc", 100, meanMass_data_det_unf.at(i)-0.2, meanMass_data_det_unf.at(i)+0.2); // temp histogram to contain PDF variations
                }

                Double_t err = -999.; //
                Double_t err_mc = -999.; //
                Int_t index_to_save = -1;
                Int_t index_to_save_mc = -1;

                // use maximum variation as systematic
                for(int j = 0; j < size_; j++){

                    if( (i==5 || i==7) && (it->first)=="Scale") continue;

                    if((it->first)=="PDFerror" || (it->first)=="Stat")
                    {
                        hpdfsys->Fill(it->second.at(j));
                        hpdfsys_mc->Fill(meanMass_sysmc_det_unf.at(i)[it->first].at(j));
                    }

                    Double_t temp_err =  fabs(meanMass_data_det_unf.at(i) - it->second.at(j));
                    Double_t temp_err_mc = fabs(meanMass_mc_det_unf.at(i) - (meanMass_sysmc_det_unf.at(i)[it->first]).at(j));

                    if( temp_err > err)
                    {
                        // update error
                        err = temp_err;
                        index_to_save = j;
                    }
                    if( temp_err_mc > err_mc)
                    {
                        err_mc = temp_err_mc;
                        index_to_save_mc = j;
                    }
                }

                if((it->first)=="PDFerror" || (it->first)=="Stat")
                {
                    err = hpdfsys->GetRMS();
                    err_mc = hpdfsys_mc->GetRMS();
                    delete hpdfsys;
                    delete hpdfsys_mc;
                }

                temp_map_[it->first] = err;
                temp_map_mc_[it->first] = err_mc;

                temp_map_sysname_index[it->first] = index_to_save;
                temp_map_mc_sysname_index[it->first] = index_to_save_mc;
                it++;
            }
            meanMassErr_sysdata_det_unf.push_back(temp_map_);
            meanMassErrIdx_sysdata_det_unf.push_back(temp_map_sysname_index);

            meanMassErr_sysmc_det_unf.push_back(temp_map_mc_);
            meanMassErrIdx_sysmc_det_unf.push_back(temp_map_mc_sysname_index);

            temp_map_.clear();
            temp_map_mc_.clear();
            temp_map_sysname_index.clear();
            temp_map_mc_sysname_index.clear();

            // FSR unfolding
            it = meanMass_sysdata_pre_fsr.at(i).begin();
            while(it != meanMass_sysdata_pre_fsr.at(i).end())
            {

                int size_ = it->second.size();

                TH1F *hpdfsys = NULL;
                TH1F *hpdfsys_mc = NULL;

                if((it->first)=="PDFerror" || (it->first)=="Stat")
                {
                    hpdfsys = new TH1F("pdfsys", "pdfsys", 100, meanMass_data_pre_fsr.at(i)-0.2, meanMass_data_pre_fsr.at(i)+0.2); // temp histogram to contain PDF variations
                    hpdfsys_mc = new TH1F("pdfsys_mc", "pdfsys_mc", 100, meanMass_data_pre_fsr.at(i)-0.2, meanMass_data_pre_fsr.at(i)+0.2); // temp histogram to contain PDF variations
                }

                Double_t err = -999.; //
                Double_t err_mc = -999.; //
                Int_t index_to_save = -1;
                Int_t index_to_save_mc = -1;
                // use maximum variation as systematic
                for(int j = 0; j < size_; j++)
                {
                    if( (i==5 || i==7) && (it->first)=="Scale") continue;

                    if((it->first)=="PDFerror" || (it->first)=="Stat")
                    {
                            hpdfsys->Fill(it->second.at(j));
                            hpdfsys_mc->Fill(meanMass_sysmc_pre_fsr.at(i)[it->first].at(j));
                    }

                    Double_t temp_err =  0.;
                    Double_t temp_err_mc =  0.;
                    if((it->first) == "QED_FSR")
                    {
                        temp_err = fabs(it->second.at(0) - it->second.at(1));
                        temp_err_mc = fabs((meanMass_sysmc_pre_fsr.at(i)[it->first]).at(0) - (meanMass_sysmc_pre_fsr.at(i)[it->first]).at(1));
                    }
                    else
                    {
                        temp_err =  fabs(meanMass_data_pre_fsr.at(i) - it->second.at(j));
                        temp_err_mc = fabs(meanMass_mc_pre_fsr.at(i) - (meanMass_sysmc_pre_fsr.at(i)[it->first]).at(j));
                        //cout << it->first << " " << meanMass_mc_pre_fsr.at(i) << " " << (meanMass_sysmc_pre_fsr.at(i)[it->first]).at(j) << endl;
                    }

                    if( temp_err > err)
                    {
                        err = temp_err;
                        index_to_save = j;
                    }

                    if( temp_err_mc > err_mc)
                    {
                        err_mc = temp_err_mc;
                        index_to_save_mc = j;
                    }
                    //cout << i << " th mass bin, " << it->first << j << " th sys value: " << it->second.at(j) << endl;
                }

                if((it->first)=="PDFerror" || (it->first)=="Stat")
                {
                    err = hpdfsys->GetRMS();
                    err_mc = hpdfsys_mc->GetRMS();
                    delete hpdfsys;
                    delete hpdfsys_mc;
                }

                temp_map_[it->first] = err;
                temp_map_mc_[it->first] = err_mc;

                temp_map_sysname_index[it->first] = index_to_save;
                temp_map_mc_sysname_index[it->first] = index_to_save_mc;
                it++;
            }
            meanMassErr_sysdata_pre_fsr.push_back(temp_map_);
            meanMassErrIdx_sysdata_pre_fsr.push_back(temp_map_sysname_index);

            meanMassErr_sysmc_pre_fsr.push_back(temp_map_mc_);
            meanMassErrIdx_sysmc_pre_fsr.push_back(temp_map_mc_sysname_index);

            temp_map_.clear();
            temp_map_mc_.clear();
            temp_map_sysname_index.clear();
            temp_map_mc_sysname_index.clear();

        }// loop for mass bins
    }

    for(int i = 0; i < nMassBin; i++)
    {
        Double_t totalSys = 0.;
        Double_t totalSys_mc = 0.;
        Double_t totalSys_pre_fsr = 0.;
        Double_t totalSys_mc_pre_fsr = 0.;

        if(doSys){
            std::map<TString, Double_t>::iterator it = meanMassErr_sysdata_det_unf.at(i).begin();
            while(it != meanMassErr_sysdata_det_unf.at(i).end())
            {

                totalSys += pow(it->second, 2);
                totalSys_mc += pow(meanMassErr_sysmc_det_unf.at(i)[it->first], 2);
                it++;
            }

            it = meanMassErr_sysdata_pre_fsr.at(i).begin();
            while(it != meanMassErr_sysdata_pre_fsr.at(i).end())
            {

                totalSys_pre_fsr += pow(it->second, 2);
                totalSys_mc_pre_fsr += pow(meanMassErr_sysmc_pre_fsr.at(i)[it->first], 2);
                it++;
            }
        }
        meanMassSysErr_mc_det_unf.push_back(sqrt(totalSys_mc));
        meanMassSysErr_data_det_unf.push_back(sqrt(totalSys));
        meanMassTotErr_data_det_unf.push_back(sqrt(totalSys + pow(meanMassErr_sysdata_det_unf.at(i)["Stat"],2)));

        meanMassSysErr_mc_pre_fsr.push_back(sqrt(totalSys_mc_pre_fsr));
        meanMassSysErr_data_pre_fsr.push_back(sqrt(totalSys_pre_fsr));
        meanMassTotErr_data_pre_fsr.push_back(sqrt(totalSys_pre_fsr + pow(meanMassErr_sysdata_pre_fsr.at(i)["Stat"],2)));
        meanMassTotRelErr_data_pre_fsr.push_back(sqrt(totalSys_pre_fsr + pow(meanMassErr_sysdata_pre_fsr.at(i)["Stat"],2))/ meanMass_data_pre_fsr.at(i) * 100.);
    }
    filein->Close();
}


// set mean pt from mass and DY mc
void ISRUnfold::setMeanPt(bool doSys, bool altMC, bool detector_unfold)
{
    int nMassBin = -1;

    const TUnfoldBinningV17* temp_binning_gen_pt = nomPtUnfold->GetOutputBinning("Gen_Pt");
    const TVectorD* temp_tvecd = temp_binning_gen_pt->GetDistributionBinning(1);
    nMassBin = temp_tvecd->GetNrows() - 1;

    TFile* filein = new TFile(hist_file_path_DYHists);

    TString massLowEdge[5] = {"50.", "65.", "80.", "100.", "200."};
    TString massHighEdge[5] = {"65.", "80.", "100.", "200.", "350."};
    TString genHistPath;

    if(channel_name=="muon")
    {
        massLowEdge[0] = "40.";
        massLowEdge[1] = "60.";

        massHighEdge[0] = "60.";
        massHighEdge[1] = "80.";
    }

    // save mean pt for each systematic variation
    for(int i = 0; i < nMassBin; i++)
    {

        TString ibinMass;
        ibinMass.Form("%d", i);

        TH1* hpt_temp_data;
        TH1* hpt_temp_mc;
        TH1* hpt_temp_mcAlt;

        // get histograms to set mean values
        TH1* hdetector_data = nomPtUnfold->GetInput("h_detector_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
        hpt_temp_data = nomPtUnfold->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
        //hpt_temp_mc   = nomPtUnfold->GetBias("histMC_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);

        genHistPath = channel_name + "_fiducial_phase_post_FSR_dRp1_m" + massLowEdge[i] + "to" + massHighEdge[i] + "/ptll/";
        hpt_temp_mc = (TH1*)filein->Get(genHistPath + "histo_DYJets");
        hpt_temp_mc->SetDirectory(0);
        hpt_temp_mc->Add((TH1*)filein->Get(genHistPath + "histo_DYJets10to50"));

        TH1* h_pre_fsr_pt_temp_data = nomPtFSRUnfold->GetOutput("h_fsr_unfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
        //TH1* h_pre_fsr_pt_temp_mc   = nomPtFSRUnfold->GetBias("histMC_pt_fsr_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
        TH1* h_pre_fsr_pt_temp_mc = NULL;

        genHistPath = channel_name + "_full_phase_m" + massLowEdge[i] + "to" + massHighEdge[i] + "/ptll_preFSR/";
        h_pre_fsr_pt_temp_mc = (TH1*)filein->Get(genHistPath + "histo_DYJets");
        h_pre_fsr_pt_temp_mc->SetDirectory(0);
        h_pre_fsr_pt_temp_mc->Add((TH1*)filein->Get(genHistPath + "histo_DYJets10to50"));

        meanPt_data_detector.   push_back(hdetector_data->GetMean());
        meanPtStatErr_data_detector.push_back(hdetector_data->GetMeanError());

        meanPt_data_det_unf.   push_back(hpt_temp_data->GetMean());
        meanPtStatErr_data_det_unf.push_back(hpt_temp_data->GetMeanError());

        meanPt_data_pre_fsr.   push_back(h_pre_fsr_pt_temp_data->GetMean());
        meanPtStatErr_data_pre_fsr.push_back(h_pre_fsr_pt_temp_data->GetMeanError());

        //doNorm(h_pre_fsr_pt_temp_data);
        //TF1* landau_nom = new TF1("landau_nom","[0]*TMath::Landau(x,[1],[2]) + [3] * expo(x * [4])", 0., 100);
        //landau_nom->SetParameters(2., 6., 2., 0.05, -0.01);
        //landau_nom->SetParLimits(0, 0.0, 5.);
        //landau_nom->SetParLimits(1, 2.0, 9.);
        //landau_nom->SetParLimits(2, 1., 10.0);
        //landau_nom->SetParLimits(3, -0.05, 0.12);
        //landau_nom->SetParLimits(4, -0.2, -0.001);
        //h_pre_fsr_pt_temp_data->Fit(landau_nom);
        //TH1* h_ = landau_nom->GetHistogram();
        //meanPt_data_pre_fsr.push_back(h_->GetMean());
        //delete landau_nom;

        meanPtStatRelErr_data_pre_fsr.push_back(h_pre_fsr_pt_temp_data->GetMeanError()/ h_pre_fsr_pt_temp_data->GetMean() * 100.);

        meanPt_mc_det_unf.   push_back(hpt_temp_mc->GetMean());
        meanPtErr_mc_det_unf.push_back(hpt_temp_mc->GetMeanError());

        meanPt_mc_pre_fsr.   push_back(h_pre_fsr_pt_temp_mc->GetMean());

        if(altMC)
        {
            hpt_temp_mcAlt   = sysPtUnfold["Alt"].at(0)->GetBias("histMC_pt_tempAlt",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
            meanPt_mcAlt.   push_back(hpt_temp_mcAlt->GetMean());
            meanPtErr_mcAlt.push_back(hpt_temp_mcAlt->GetMeanError());
            delete hpt_temp_mcAlt;
        }

        if(doSys)
        {
            // detector unfolding
            // temp map to save systematic results for a mass bin
            // format: systematic name, mean pt
            std::map<TString, std::vector<Double_t>> temp_map;
            std::map<TString, std::vector<Double_t>> temp_map_mc;

            std::map<TString, std::vector<TUnfoldDensityV17*>>::iterator it = sysPtUnfold.begin();
            while(it != sysPtUnfold.end())
            {
                int nSys = it->second.size();
                TH1* hdatasys_temp;
                for(int isys = 0; isys < nSys; isys++)
                {
                    hdatasys_temp = sysPtUnfold[it->first].at(isys)->GetOutput("hunfolded_pt_systemp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
                    temp_map[it->first].push_back(hdatasys_temp->GetMean());

                    //doNorm(hdatasys_temp);
                    //TF1* landau = new TF1("landau","[0]*TMath::Landau(x,[1],[2]) + [3] * expo(x * [4])", 0., 100);
                    //landau->SetParameters(2., 6., 2., 0.05, -0.01);
                    //landau->SetParLimits(0, 0.0, 5.);
                    //landau->SetParLimits(1, 2.0, 9.);
                    //landau->SetParLimits(2, 1., 10.0);
                    //landau->SetParLimits(3, -0.05, 0.12);
                    //landau->SetParLimits(4, -0.2, -0.001);
                    //hdatasys_temp->Fit(landau);
                    //TH1* h = landau->GetHistogram();
                    //temp_map[it->first].push_back(h->GetMean());
                    

                    //delete landau;
                    delete hdatasys_temp;

                    hdatasys_temp = sysPtUnfold[it->first].at(isys)->GetBias("hunfolded_pt_systemp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
                    temp_map_mc[it->first].push_back(hdatasys_temp->GetMean());

                    delete hdatasys_temp;
                }
                it++;
            }
            meanPt_sysdata_det_unf.push_back(temp_map);
            meanPt_sysmc_det_unf.push_back(temp_map_mc);
            temp_map.clear();
            temp_map_mc.clear();

            // QED FSR unfolding
            it = sysPtFSRUnfold.begin();
            while(it != sysPtFSRUnfold.end())
            {

                int nSys = it->second.size();
                TH1* hdatasys_temp;
                for(int isys = 0; isys < nSys; isys++)
                {
                    hdatasys_temp = sysPtFSRUnfold[it->first].at(isys)->GetOutput("hunfolded_pt_systemp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
                    temp_map[it->first].push_back(hdatasys_temp->GetMean());

                    //doNorm(hdatasys_temp);
                    //TF1* landau = new TF1("landau","[0]*TMath::Landau(x,[1],[2]) + [3] * expo(x * [4])", 0., 100);
                    //landau->SetParameters(2., 6., 2., 0.05, -0.01);
                    //landau->SetParLimits(0, 0.0, 5.);
                    //landau->SetParLimits(1, 2.0, 9.);
                    //landau->SetParLimits(2, 1., 10.0);
                    //landau->SetParLimits(3, -0.05, 0.12);
                    //landau->SetParLimits(4, -0.2, -0.001);
                    //hdatasys_temp->Fit(landau);
                    //TH1* h = landau->GetHistogram();
                    //temp_map[it->first].push_back(h->GetMean());
                    
                    //delete landau;
                    delete hdatasys_temp;

                    hdatasys_temp = sysPtFSRUnfold[it->first].at(isys)->GetBias("hunfolded_pt_systemp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
                    temp_map_mc[it->first].push_back(hdatasys_temp->GetMean());
                    delete hdatasys_temp;

                }
                it++;
            }
            meanPt_sysdata_pre_fsr.push_back(temp_map);
            meanPt_sysmc_pre_fsr.push_back(temp_map_mc);
            temp_map.clear();
            temp_map_mc.clear();
        }

            
        delete hdetector_data;
        delete hpt_temp_data;
        delete hpt_temp_mc;
        delete h_pre_fsr_pt_temp_data;
        delete h_pre_fsr_pt_temp_mc;
    }

    // calculate systematic uncertainty for each systematic variation
    if(doSys)
    {
	for(int i = 0; i < nMassBin; i++)
        {
            // detector unfolding
            std::map<TString, Double_t> temp_map_;
            std::map<TString, Double_t> temp_map_mc_;
            std::map<TString, Int_t> temp_map_sysname_index;
            std::map<TString, Int_t> temp_map_mc_sysname_index;

            std::map<TString, std::vector<Double_t>>::iterator it = meanPt_sysdata_det_unf.at(i).begin();
            while(it != meanPt_sysdata_det_unf.at(i).end())
            {
	    	int size_ = it->second.size(); // size of systematic variations
                	
	    	TH1F *hpdfsys = NULL;
	    	TH1F *hpdfsys_mc = NULL;
	    	if((it->first)=="PDFerror" || (it->first)=="Stat")
                {
                    hpdfsys = new TH1F("pdfsys", "pdfsys", 100, meanPt_data_det_unf.at(i)-0.2, meanPt_data_det_unf.at(i)+0.2); // temp histogram to contain PDF variations
                    hpdfsys_mc = new TH1F("pdfsys_mc", "pdfsys_mc", 100, meanPt_data_det_unf.at(i)-0.2, meanPt_data_det_unf.at(i)+0.2); // temp histogram to contain PDF variations
                }

	    	Double_t err = -999.; //
	    	Double_t err_mc = -999.; //
                Int_t index_to_save = -1;
                Int_t index_to_save_mc = -1;
	    	for(int j = 0; j < size_; j++)
                {
	    	    if( (i==5 || i==7) && (it->first)=="Scale") continue;

	    	    if((it->first)=="PDFerror" || (it->first)=="Stat")
                    {
	     		hpdfsys->Fill(it->second.at(j));
	     		hpdfsys_mc->Fill(meanPt_sysdata_det_unf.at(i)[it->first].at(j));
	    	    }

                    Double_t temp_err =  fabs(meanPt_data_det_unf.at(i) - it->second.at(j));
                    Double_t temp_err_mc =  fabs(meanPt_mc_det_unf.at(i) - (meanPt_sysmc_det_unf.at(i)[it->first]).at(j));

                    if( temp_err > err){
                        err = temp_err;
                        index_to_save = j;
                    }
                    if( temp_err_mc > err_mc)
                    {
                        err_mc = temp_err_mc;
                        index_to_save_mc = j;
                    }
	    	    //cout << i << " th mass bin, " << it->first << j << " th sys value: " << it->second.at(j) << endl;
	    	}// loop for systematic variations
	    	if((it->first)=="PDFerror" || (it->first)=="Stat")
                {
	    	    err = hpdfsys->GetRMS();
	    	    err_mc = hpdfsys_mc->GetRMS();
	    	    delete hpdfsys;
	    	    delete hpdfsys_mc;
	    	}
	    	temp_map_[it->first] = err;
	    	temp_map_mc_[it->first] = err_mc;

                temp_map_sysname_index[it->first] = index_to_save;
                temp_map_mc_sysname_index[it->first] = index_to_save_mc;
	    	it++;
            }// loop for systematic sources
	    meanPtErr_sysdata_det_unf.push_back(temp_map_);
            meanPtErrIdx_sysdata_det_unf.push_back(temp_map_sysname_index);

            meanPtErr_sysmc_det_unf.push_back(temp_map_mc_);
            meanPtErrIdx_sysmc_det_unf.push_back(temp_map_mc_sysname_index);

            temp_map_.clear();
            temp_map_mc_.clear();
            temp_map_sysname_index.clear();
            temp_map_mc_sysname_index.clear();

            // QED FSR unfolding
            it = meanPt_sysdata_pre_fsr.at(i).begin();
            while(it != meanPt_sysdata_pre_fsr.at(i).end())
            {
                int size_ = it->second.size(); // size of systematic variations

                TH1F *hpdfsys = NULL;
                TH1F *hpdfsys_mc = NULL;
                if((it->first)=="PDFerror" || (it->first)=="Stat")
                {
                    hpdfsys = new TH1F("pdfsys", "pdfsys", 100, meanPt_data_pre_fsr.at(i)-0.2, meanPt_data_pre_fsr.at(i)+0.2); // temp histogram to contain PDF variations
                    hpdfsys_mc = new TH1F("pdfsys_mc", "pdfsys_mc", 100, meanPt_data_pre_fsr.at(i)-0.2, meanPt_data_pre_fsr.at(i)+0.2); // temp histogram to contain PDF variations
                }

                Double_t err = -999.; //
                Double_t err_mc = -999.; //
                Int_t index_to_save = -1;
                Int_t index_to_save_mc = -1;
                for(int j = 0; j < size_; j++)
                {
	    	    if( (i==5 || i==7) && (it->first)=="Scale") continue;

                    if((it->first)=="PDFerror" || (it->first)=="Stat")
                    {
                        hpdfsys->Fill(it->second.at(j));
                        hpdfsys_mc->Fill(meanPt_sysmc_pre_fsr.at(i)[it->first].at(j));
                    }

                    Double_t temp_err =  0.;
                    Double_t temp_err_mc =  0.;
                    if((it->first) == "QED_FSR")
                    {
                        temp_err =  fabs(it->second.at(0) - it->second.at(1));
                        temp_err_mc = fabs((meanPt_sysmc_pre_fsr.at(i)[it->first]).at(0) - (meanPt_sysmc_pre_fsr.at(i)[it->first]).at(1));
                        cout << "QED fsr check" << endl;
                        cout << "mean1: " << it->second.at(0) << " mean2: " << it->second.at(1) << endl;
                    }
                    else
                    {
                        temp_err =  fabs(meanPt_data_pre_fsr.at(i) - it->second.at(j));
                        temp_err_mc = fabs(meanPt_mc_pre_fsr.at(i) - (meanPt_sysmc_pre_fsr.at(i)[it->first]).at(j));
                    }

                    if( temp_err > err)
                    {
                        err = temp_err;
                        index_to_save = j;
                    }
                    if( temp_err_mc > err_mc)
                    {
                        err_mc = temp_err_mc;
                        index_to_save_mc = j;
                    }
                }// loop for systematic variations

                if((it->first)=="PDFerror" || (it->first)=="Stat")
                {
                    err = hpdfsys->GetRMS(); // only for 2016 DY MC
                    err_mc = hpdfsys_mc->GetRMS(); // only for 2016 DY MC

                    delete hpdfsys;
                    delete hpdfsys_mc;
                }

                temp_map_[it->first] = err;
                temp_map_mc_[it->first] = err_mc;

                // save the index of systematic distribution giving maximum variation on the mean value
                temp_map_sysname_index[it->first] = index_to_save;
                temp_map_mc_sysname_index[it->first] = index_to_save_mc;
                it++;
            }// loop for systematic sources
            meanPtErr_sysdata_pre_fsr.push_back(temp_map_);
            meanPtErrIdx_sysdata_pre_fsr.push_back(temp_map_sysname_index);

            meanPtErr_sysmc_pre_fsr.push_back(temp_map_mc_);
            meanPtErrIdx_sysmc_pre_fsr.push_back(temp_map_mc_sysname_index);

            temp_map_.clear();
            temp_map_sysname_index.clear();
            temp_map_mc_sysname_index.clear();

        }// loop for mass binss
    }

    for(int i = 0; i < nMassBin; i++)
    {
        Double_t totalSys = 0.;
        Double_t totalSys_mc = 0.;
        Double_t totalSys_pre_fsr = 0.;
        Double_t totalSys_mc_pre_fsr = 0.;

        if(doSys)
        {
            cout << i << " th mass bin... " << endl;
            std::map<TString, Double_t>::iterator it = meanPtErr_sysdata_det_unf.at(i).begin();
            while(it != meanPtErr_sysdata_det_unf.at(i).end())
            {
                cout << "systematic detector unfold: " << it->first << " " << it->second/ meanPt_data_det_unf.at(i) * 100.<< endl;
                totalSys += pow(it->second, 2);	
                totalSys_mc += pow(meanPtErr_sysmc_det_unf.at(i)[it->first], 2);	
                it++;
            }

            it = meanPtErr_sysdata_pre_fsr.at(i).begin();
            while(it != meanPtErr_sysdata_pre_fsr.at(i).end())
            {
                // TODO make plots(TGraph) to show relative and absolute  errors
                // make txt file to save systematic values
                cout << "systematic FSR unfold: " << it->first << " " << it->second/ meanPt_data_pre_fsr.at(i) * 100.<< endl;
                totalSys_pre_fsr += pow(it->second, 2);	
                totalSys_mc_pre_fsr += pow(meanPtErr_sysmc_pre_fsr.at(i)[it->first], 2);
                it++;
            }
        }
        meanPtSysErr_mc_det_unf.push_back(sqrt(totalSys_mc));
        meanPtSysErr_data_det_unf.push_back(sqrt(totalSys));
        meanPtTotErr_data_det_unf.push_back(sqrt(totalSys + pow(meanPtErr_sysdata_det_unf.at(i)["Stat"],2)));

        meanPtSysErr_mc_pre_fsr.push_back(sqrt(totalSys_mc_pre_fsr));
        meanPtSysErr_data_pre_fsr.push_back(sqrt(totalSys_pre_fsr));
        meanPtTotErr_data_pre_fsr.push_back(sqrt(totalSys_pre_fsr + pow(meanPtErr_sysdata_pre_fsr.at(i)["Stat"],2)));
        meanPtTotRelErr_data_pre_fsr.push_back(sqrt(totalSys_pre_fsr + pow(meanPtErr_sysdata_pre_fsr.at(i)["Stat"],2))/ meanPt_data_pre_fsr.at(i) * 100.);
    }// loop for mass bins
    filein->Close();
}

void ISRUnfold::drawISRRun2results(TString outpdf, TCanvas* c_2017, TCanvas* c_2018, TCanvas* c_muon_2016, TCanvas* c_muon_2017, TCanvas* c_muon_2018)
{

    TFile* filein_dy = new TFile("/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/rootScripts/gen_level_isr.root");
    TFile* filein_dy_mg = new TFile("/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/rootScripts/gen_level_isr_mg.root");
    TFile* filein_dy_8TeV = new TFile("/home/jhkim/ISR2016/8TeV/isr_8TeV.root");

    TH1* hmcPt;
    TH1* hmcMass;

    vector<double> mean_pt, mean_pt_stat_error;
    vector<double> mean_mass, mean_mass_stat_error;

    vector<double> mean_pt_100cut, mean_pt_100cut_stat_error, mean_pt_100cut_sys_error;
    vector<double> mean_mass_100cut, mean_mass_100cut_stat_error, mean_mass_100cut_sys_error;

    vector<double> mean_pt_100cut_MG, mean_pt_100cut_stat_error_MG, mean_pt_100cut_sys_error_MG;
    vector<double> mean_mass_100cut_MG, mean_mass_100cut_stat_error_MG, mean_mass_100cut_sys_error_MG;

    vector<double> mean_pt_100cut_8tev, mean_pt_100cut_stat_error_8tev, mean_pt_100cut_sys_error_8tev;
    vector<double> mean_mass_100cut_8tev, mean_mass_100cut_stat_error_8tev, mean_mass_100cut_sys_error_8tev;

    TString dir_name = "m80to100";

    TString mass_bin[19] = {"20","40","60","80","100","120","140","160","180","200","220","240","260","280","300","320","340","360","380"};
    std::map<TString, int> sys_map;
    sys_map["AlphaS"] = 2;
    sys_map["Scale"] = 6;

    for(int i = 0; i < 18; i++)
    {

        dir_name = "m" + mass_bin[i] + "to" + mass_bin[i+1];
        //cout << dir_name << endl;

        hmcPt = (TH1*)filein_dy->Get(dir_name + "/hist_ptll/histo_DYJets");
        hmcPt->Add((TH1*)filein_dy->Get(dir_name + "/hist_ptll/histo_DYJets10to50"));

        hmcMass = (TH1*)filein_dy->Get(dir_name + "/hist_mll/histo_DYJets");
        hmcMass->Add((TH1*)filein_dy->Get(dir_name + "/hist_mll/histo_DYJets10to50"));

        mean_pt.push_back(hmcPt->GetMean());
        mean_pt_stat_error.push_back(hmcPt->GetMeanError());

        mean_mass.push_back(hmcMass->GetMean());
        mean_mass_stat_error.push_back(hmcMass->GetMeanError());

        delete hmcPt, hmcMass;

        dir_name = "m" + mass_bin[i] + "to" + mass_bin[i+1] + "_pt100";
        hmcPt = (TH1*)filein_dy->Get(dir_name + "/hist_ptll/histo_DYJets");
        hmcPt->Add((TH1*)filein_dy->Get(dir_name + "/hist_ptll/histo_DYJets10to50"));

        hmcMass = (TH1*)filein_dy->Get(dir_name + "/hist_mll/histo_DYJets");
        hmcMass->Add((TH1*)filein_dy->Get(dir_name + "/hist_mll/histo_DYJets10to50"));

        mean_pt_100cut.push_back(hmcPt->GetMean());
        mean_pt_100cut_stat_error.push_back(hmcPt->GetMeanError());

        mean_mass_100cut.push_back(hmcMass->GetMean());
        mean_mass_100cut_stat_error.push_back(hmcMass->GetMeanError());

        delete hmcPt, hmcMass;

        dir_name = "m" + mass_bin[i] + "to" + mass_bin[i+1] + "_pt100";
        hmcPt = (TH1*)filein_dy_mg->Get(dir_name + "/hist_ptll/histo_DYJets_MG");
        hmcPt->Add((TH1*)filein_dy_mg->Get(dir_name + "/hist_ptll/histo_DYJets10to50_MG"));

        hmcMass = (TH1*)filein_dy_mg->Get(dir_name + "/hist_mll/histo_DYJets_MG");
        hmcMass->Add((TH1*)filein_dy_mg->Get(dir_name + "/hist_mll/histo_DYJets10to50_MG"));

        mean_pt_100cut_MG.push_back(hmcPt->GetMean());
        mean_pt_100cut_stat_error_MG.push_back(hmcPt->GetMeanError());

        mean_mass_100cut_MG.push_back(hmcMass->GetMean());
        mean_mass_100cut_stat_error_MG.push_back(hmcMass->GetMeanError());

        delete hmcPt, hmcMass;

        dir_name = "m" + mass_bin[i] + "to" + mass_bin[i+1];
        hmcPt = (TH1*)filein_dy_8TeV->Get("dimuPt_" + dir_name);

        hmcMass = (TH1*)filein_dy_8TeV->Get("dimuMass_" + dir_name);

        mean_pt_100cut_8tev.push_back(hmcPt->GetMean());
        mean_pt_100cut_stat_error_8tev.push_back(hmcPt->GetMeanError());

        mean_mass_100cut_8tev.push_back(hmcMass->GetMean());
        mean_mass_100cut_stat_error_8tev.push_back(hmcMass->GetMeanError());

        delete hmcPt, hmcMass;

        double tot_sys_pt_err = 0;
        double tot_sys_mass_err = 0;
        std::map<TString, int>::iterator it = sys_map.begin();
        while(it != sys_map.end())
        {
            
            double sys_pt_err = -999.;
            double sys_mass_err = -999.;

            for(int ith_sys = 0; ith_sys < it->second; ith_sys++)
            {
                TString systematic_postfix = it->first;

                double temp_sys_pt_err = -999.;
                double temp_sys_mass_err = -999.;

                if(it->first == "AlphaS"){
                    if(ith_sys == 0) systematic_postfix+="Up";
                    if(ith_sys == 1) systematic_postfix+="Down";
                }
                if(it->first == "Scale"){
                    if(ith_sys == 0 ) systematic_postfix+="AUp";
                    if(ith_sys == 1 ) systematic_postfix+="ADown";
                    if(ith_sys == 2 ) systematic_postfix+="BUp";
                    if(ith_sys == 3 ) systematic_postfix+="BDown";
                    if(ith_sys == 4 ) systematic_postfix+="ABUp";
                    if(ith_sys == 5 ) systematic_postfix+="ABDown";
                }
                dir_name = "m" + mass_bin[i] + "to" + mass_bin[i+1] + "_pt100";

                hmcPt = (TH1*)filein_dy->Get(dir_name + "/hist_ptll/histo_DYJets_"+systematic_postfix);
                //cout << dir_name + "/hist_ptll/histo_DYJets10to50_"+systematic_postfix << endl;
                hmcPt->Add((TH1*)filein_dy->Get(dir_name + "/hist_ptll/histo_DYJets10to50_"+systematic_postfix));

                hmcMass = (TH1*)filein_dy->Get(dir_name + "/hist_mll/histo_DYJets_"+systematic_postfix);
                hmcMass->Add((TH1*)filein_dy->Get(dir_name + "/hist_mll/histo_DYJets10to50_"+systematic_postfix));

                temp_sys_pt_err = fabs(hmcPt->GetMean()-mean_pt_100cut.at(i));
                temp_sys_mass_err = fabs(hmcMass->GetMean()-mean_mass_100cut.at(i));

                if(sys_pt_err < temp_sys_pt_err){
                    sys_pt_err = temp_sys_pt_err;
                }
                if(sys_mass_err < temp_sys_mass_err){
                    sys_mass_err = temp_sys_mass_err;
                }
            }
            tot_sys_pt_err += pow(sys_pt_err,2);
            tot_sys_mass_err += pow(sys_mass_err,2);
            it++;
        }
        mean_pt_100cut_sys_error.push_back(sqrt(tot_sys_pt_err));
        mean_mass_100cut_sys_error.push_back(sqrt(tot_sys_mass_err));

        tot_sys_pt_err = 0;
        tot_sys_mass_err = 0;
        it = sys_map.begin();
        while(it != sys_map.end()){

            double sys_pt_err = -999.;
            double sys_mass_err = -999.;

            for(int ith_sys = 0; ith_sys < it->second; ith_sys++){
                TString systematic_postfix = it->first;

                double temp_sys_pt_err = -999.;
                double temp_sys_mass_err = -999.;

                if(it->first == "AlphaS"){
                    if(ith_sys == 0) systematic_postfix+="Up";
                    if(ith_sys == 1) systematic_postfix+="Down";
                }
                if(it->first == "Scale"){
                    if(ith_sys == 0 ) systematic_postfix+="AUp";
                    if(ith_sys == 1 ) systematic_postfix+="ADown";
                    if(ith_sys == 2 ) systematic_postfix+="BUp";
                    if(ith_sys == 3 ) systematic_postfix+="BDown";
                    if(ith_sys == 4 ) systematic_postfix+="ABUp";
                    if(ith_sys == 5 ) systematic_postfix+="ABDown";
                }
                dir_name = "m" + mass_bin[i] + "to" + mass_bin[i+1] + "_pt100";

                hmcPt = (TH1*)filein_dy_mg->Get(dir_name + "/hist_ptll/histo_DYJets_MG_"+systematic_postfix);
                hmcPt->Add((TH1*)filein_dy_mg->Get(dir_name + "/hist_ptll/histo_DYJets10to50_MG_"+systematic_postfix));

                hmcMass = (TH1*)filein_dy_mg->Get(dir_name + "/hist_mll/histo_DYJets_MG_"+systematic_postfix);
                hmcMass->Add((TH1*)filein_dy_mg->Get(dir_name + "/hist_mll/histo_DYJets10to50_MG_"+systematic_postfix));

                temp_sys_pt_err = fabs(hmcPt->GetMean()-mean_pt_100cut.at(i));
                temp_sys_mass_err = fabs(hmcMass->GetMean()-mean_mass_100cut.at(i));

                if(sys_pt_err < temp_sys_pt_err){
                    sys_pt_err = temp_sys_pt_err;
                }
                if(sys_mass_err < temp_sys_mass_err){
                    sys_mass_err = temp_sys_mass_err;
                }
            }
            tot_sys_pt_err += pow(sys_pt_err,2);
            tot_sys_mass_err += pow(sys_mass_err,2);
            it++;
        }
        mean_pt_100cut_sys_error_MG.push_back(sqrt(tot_sys_pt_err));
        mean_mass_100cut_sys_error_MG.push_back(sqrt(tot_sys_mass_err));

    }

    bool doFit = false;
    gROOT->SetBatch();

    int marker_ = 20;

    setTDRStyle();
    writeExtraText = true;       // if extra text
    extraText  = "work in progress";

    TString year_string = "Run2";

    c1 = new TCanvas("c_Run2_preFSR", "c_Run2_preFSR", 50, 50, 900*1.5, 700*1.5);
    c1->cd();
    gStyle->SetOptFit(0);

    c1->SetBottomMargin(0.2);
    c1->SetTopMargin(0.08);
    //c1->SetTicks();
    c1->SetLogx();
    c1->SetGridy();
    c1->SetGridx();

    double marker_size = 1.5;

    // 2016 results
    TGraphErrors *grUnfolded_pre_fsr = new TGraphErrors(5, &meanMass_data_pre_fsr[0], &meanPt_data_pre_fsr[0], &meanMassTotErr_data_pre_fsr[0], &meanPtTotErr_data_pre_fsr[0]);
    grUnfolded_pre_fsr->SetLineColor(kBlack);
    grUnfolded_pre_fsr->SetMarkerColor(kBlack);
    grUnfolded_pre_fsr->SetMarkerStyle(marker_);
    grUnfolded_pre_fsr->SetMarkerSize(marker_size);
    grUnfolded_pre_fsr->SetLineStyle(1);
    grUnfolded_pre_fsr->Draw("apZ");
    grUnfolded_pre_fsr->GetYaxis()->SetRangeUser(11.,32.);
    grUnfolded_pre_fsr->GetXaxis()->SetLimits(30.,500.);
    grUnfolded_pre_fsr->GetXaxis()->SetMoreLogLabels(true);
    grUnfolded_pre_fsr->GetYaxis()->SetTitle("\\mbox{Average } p_{T}^{\\ell\\ell} \\mbox{ (GeV)}"); // \\ell not working with pdf
    grUnfolded_pre_fsr->GetYaxis()->SetLabelSize(23*2);
    grUnfolded_pre_fsr->GetYaxis()->SetTitleSize(23*2);
    grUnfolded_pre_fsr->GetYaxis()->SetTitleOffset(1.5);
    grUnfolded_pre_fsr->GetXaxis()->SetTitle("\\mbox{Average } Mass^{\\ell\\ell} \\mbox{ (GeV)}");
    grUnfolded_pre_fsr->GetXaxis()->SetLabelSize(23*2);
    grUnfolded_pre_fsr->GetXaxis()->SetTitleSize(23*2);
    grUnfolded_pre_fsr->GetXaxis()->SetTitleOffset(2.);

    TGraphErrors* grUnfolded_pre_fsr_2017 = (TGraphErrors*)c_2017->GetListOfPrimitives()->FindObject("preFSRUnfoldedData_electron_2017");
    TGraphErrors* grUnfolded_pre_fsr_2018 = (TGraphErrors*)c_2018->GetListOfPrimitives()->FindObject("preFSRUnfoldedData_electron_2018");

    TGraphErrors* grUnfolded_pre_fsr_muon_2016 = (TGraphErrors*)c_muon_2016->GetListOfPrimitives()->FindObject("preFSRUnfoldedData_muon_2016");
    TGraphErrors* grUnfolded_pre_fsr_muon_2017 = (TGraphErrors*)c_muon_2017->GetListOfPrimitives()->FindObject("preFSRUnfoldedData_muon_2017");
    TGraphErrors* grUnfolded_pre_fsr_muon_2018 = (TGraphErrors*)c_muon_2018->GetListOfPrimitives()->FindObject("preFSRUnfoldedData_muon_2018");

    // Get points
    double* p_x_electron_2016 = grUnfolded_pre_fsr->GetX();
    double* p_y_electron_2016 = grUnfolded_pre_fsr->GetY();
    double* p_ex_electron_2016 = grUnfolded_pre_fsr->GetEX();
    double* p_ey_electron_2016 = grUnfolded_pre_fsr->GetEY();

    double* p_x_electron_2017 = grUnfolded_pre_fsr_2017->GetX();
    double* p_y_electron_2017 = grUnfolded_pre_fsr_2017->GetY();
    double* p_ex_electron_2017 = grUnfolded_pre_fsr_2017->GetEX();
    double* p_ey_electron_2017 = grUnfolded_pre_fsr_2017->GetEY();

    double* p_x_electron_2018 = grUnfolded_pre_fsr_2018->GetX();
    double* p_y_electron_2018 = grUnfolded_pre_fsr_2018->GetY();
    double* p_ex_electron_2018 = grUnfolded_pre_fsr_2018->GetEX();
    double* p_ey_electron_2018 = grUnfolded_pre_fsr_2018->GetEY();

    double* p_x_muon_2016 = grUnfolded_pre_fsr_muon_2016->GetX();
    double* p_y_muon_2016 = grUnfolded_pre_fsr_muon_2016->GetY();
    double* p_ex_muon_2016 = grUnfolded_pre_fsr_muon_2016->GetEX();
    double* p_ey_muon_2016 = grUnfolded_pre_fsr_muon_2016->GetEY();

    double* p_x_muon_2017 = grUnfolded_pre_fsr_muon_2017->GetX();
    double* p_y_muon_2017 = grUnfolded_pre_fsr_muon_2017->GetY();
    double* p_ex_muon_2017 = grUnfolded_pre_fsr_muon_2017->GetEX();
    double* p_ey_muon_2017 = grUnfolded_pre_fsr_muon_2017->GetEY();

    double* p_x_muon_2018 = grUnfolded_pre_fsr_muon_2018->GetX();
    double* p_y_muon_2018 = grUnfolded_pre_fsr_muon_2018->GetY();
    double* p_ex_muon_2018 = grUnfolded_pre_fsr_muon_2018->GetEX();
    double* p_ey_muon_2018 = grUnfolded_pre_fsr_muon_2018->GetEY();

    // Save Run 2 electron points into one array
    double run2_electron_mass[15];
    double run2_electron_pt[15];
    double run2_electron_mass_error[15];
    double run2_electron_pt_error[15];

    double run2_muon_mass[15];
    double run2_muon_pt[15];
    double run2_muon_mass_error[15];
    double run2_muon_pt_error[15];

    double run2_lepton_mass[30];
    double run2_lepton_pt[30];
    double run2_lepton_mass_error[30];
    double run2_lepton_pt_error[30];

    for(int ip = 0; ip < 5; ip++)
    {

        run2_electron_mass[ip] = p_x_electron_2016[ip];
        run2_electron_pt[ip] = p_y_electron_2016[ip];
        run2_electron_mass_error[ip] = p_ex_electron_2016[ip];
        run2_electron_pt_error[ip] = p_ey_electron_2016[ip];

        run2_electron_mass[ip + 5] = p_x_electron_2017[ip];
        run2_electron_pt[ip + 5] = p_y_electron_2017[ip];
        run2_electron_mass_error[ip + 5] = p_ex_electron_2017[ip];
        run2_electron_pt_error[ip + 5] = p_ey_electron_2017[ip];

        run2_electron_mass[ip + 10] = p_x_electron_2018[ip];
        run2_electron_pt[ip + 10] = p_y_electron_2018[ip];
        run2_electron_mass_error[ip + 10] = p_ex_electron_2018[ip];
        run2_electron_pt_error[ip + 10] = p_ey_electron_2018[ip];

        run2_muon_mass[ip] = p_x_muon_2016[ip];
        run2_muon_pt[ip] = p_y_muon_2016[ip];
        run2_muon_mass_error[ip] = p_ex_muon_2016[ip];
        run2_muon_pt_error[ip] = p_ey_muon_2016[ip];

        run2_muon_mass[ip + 5] = p_x_muon_2017[ip];
        run2_muon_pt[ip + 5] = p_y_muon_2017[ip];
        run2_muon_mass_error[ip + 5] = p_ex_muon_2017[ip];
        run2_muon_pt_error[ip + 5] = p_ey_muon_2017[ip];

        run2_muon_mass[ip + 10] = p_x_muon_2018[ip];
        run2_muon_pt[ip + 10] = p_y_muon_2018[ip];
        run2_muon_mass_error[ip + 10] = p_ex_muon_2018[ip];
        run2_muon_pt_error[ip + 10] = p_ey_muon_2018[ip];

        // lepton
        run2_lepton_mass[ip] = p_x_electron_2016[ip];
        run2_lepton_pt[ip] = p_y_electron_2016[ip];
        run2_lepton_mass_error[ip] = p_ex_electron_2016[ip];
        run2_lepton_pt_error[ip] = p_ey_electron_2016[ip];

        run2_lepton_mass[ip + 5] = p_x_electron_2017[ip];
        run2_lepton_pt[ip + 5] = p_y_electron_2017[ip];
        run2_lepton_mass_error[ip + 5] = p_ex_electron_2017[ip];
        run2_lepton_pt_error[ip + 5] = p_ey_electron_2017[ip];

        run2_lepton_mass[ip + 10] = p_x_electron_2018[ip];
        run2_lepton_pt[ip + 10] = p_y_electron_2018[ip];
        run2_lepton_mass_error[ip + 10] = p_ex_electron_2018[ip];
        run2_lepton_pt_error[ip + 10] = p_ey_electron_2018[ip];

        run2_lepton_mass[ip + 15] = p_x_muon_2016[ip];
        run2_lepton_pt[ip + 15] = p_y_muon_2016[ip];
        run2_lepton_mass_error[ip + 15] = p_ex_muon_2016[ip];
        run2_lepton_pt_error[ip + 15] = p_ey_muon_2016[ip];

        run2_lepton_mass[ip + 20] = p_x_muon_2017[ip];
        run2_lepton_pt[ip + 20] = p_y_muon_2017[ip];
        run2_lepton_mass_error[ip + 20] = p_ex_muon_2017[ip];
        run2_lepton_pt_error[ip + 20] = p_ey_muon_2017[ip];

        run2_lepton_mass[ip + 25] = p_x_muon_2018[ip];
        run2_lepton_pt[ip + 25] = p_y_muon_2018[ip];
        run2_lepton_mass_error[ip + 25] = p_ex_muon_2018[ip];
        run2_lepton_pt_error[ip + 25] = p_ey_muon_2018[ip];

    }

    TGraphErrors *grUnfolded_pre_fsr_electron_run2 = new TGraphErrors(15, run2_electron_mass, run2_electron_pt, run2_electron_mass_error, run2_electron_pt_error);
    TGraphErrors *grUnfolded_pre_fsr_muon_run2 = new TGraphErrors(15, run2_muon_mass, run2_muon_pt, run2_muon_mass_error, run2_muon_pt_error);
    TGraphErrors *grUnfolded_pre_fsr_lepton_run2 = new TGraphErrors(30, run2_lepton_mass, run2_lepton_pt, run2_lepton_mass_error, run2_lepton_pt_error);

    //for(int ix = 0; ix < 5; ix++)
    //{
    //    cout << "p_x: " << p_x_electron_2017[ix] <<  " +/-" << p_ex_electron_2017[ix] << endl;
    //    cout << "p_y: " << p_y_electron_2017[ix] <<  " +/-" << p_ey_electron_2017[ix] << endl;
    //}

    TGraphErrors* grMC_pre_fsr_2017 = (TGraphErrors*)c_2017->GetListOfPrimitives()->FindObject("preFSRMC_electron_2017");
    TGraphErrors* grMC_pre_fsr_2018 = (TGraphErrors*)c_2018->GetListOfPrimitives()->FindObject("preFSRMC_electron_2018");

    grUnfolded_pre_fsr_2017->SetLineColor(kBlack);
    grUnfolded_pre_fsr_2017->SetMarkerColor(kBlack);
    grUnfolded_pre_fsr_2017->SetMarkerStyle(marker_+5);
    grUnfolded_pre_fsr_2017->SetMarkerSize(marker_size);
    grUnfolded_pre_fsr_2017->SetLineStyle(1);
    grUnfolded_pre_fsr_2017->Draw("pZ same");

    grUnfolded_pre_fsr_2018->SetLineColor(kBlack);
    grUnfolded_pre_fsr_2018->SetMarkerColor(kBlack);
    grUnfolded_pre_fsr_2018->SetMarkerStyle(marker_+6);
    grUnfolded_pre_fsr_2018->SetMarkerSize(marker_size);
    grUnfolded_pre_fsr_2018->SetLineStyle(1);
    grUnfolded_pre_fsr_2018->Draw("pZ same");

    grUnfolded_pre_fsr_muon_2016->SetLineColor(kRed);
    grUnfolded_pre_fsr_muon_2016->SetMarkerColor(kRed);
    grUnfolded_pre_fsr_muon_2016->SetMarkerStyle(marker_);
    grUnfolded_pre_fsr_muon_2016->SetMarkerSize(marker_size);
    grUnfolded_pre_fsr_muon_2016->SetLineStyle(1);
    grUnfolded_pre_fsr_muon_2016->Draw("pZ same");

    grUnfolded_pre_fsr_muon_2017->SetLineColor(kRed);
    grUnfolded_pre_fsr_muon_2017->SetMarkerColor(kRed);
    grUnfolded_pre_fsr_muon_2017->SetMarkerStyle(marker_+5);
    grUnfolded_pre_fsr_muon_2017->SetMarkerSize(marker_size);
    grUnfolded_pre_fsr_muon_2017->SetLineStyle(1);
    grUnfolded_pre_fsr_muon_2017->Draw("pZ same");

    grUnfolded_pre_fsr_muon_2018->SetLineColor(kRed);
    grUnfolded_pre_fsr_muon_2018->SetMarkerColor(kRed);
    grUnfolded_pre_fsr_muon_2018->SetMarkerStyle(marker_+6);
    grUnfolded_pre_fsr_muon_2018->SetMarkerSize(marker_size);
    grUnfolded_pre_fsr_muon_2018->SetLineStyle(1);
    grUnfolded_pre_fsr_muon_2018->Draw("pZ same");

    //grMC_pre_fsr_2017->SetFillColorAlpha(kBlack,0.3);
    //grMC_pre_fsr_2017->Draw("E3 same");

    grMC_pre_fsr_2017->SetMarkerColor(kRed);
    grMC_pre_fsr_2017->SetMarkerStyle(21);
    grMC_pre_fsr_2017->SetMarkerSize(marker_size);
    grMC_pre_fsr_2017->SetLineStyle(1);
    grMC_pre_fsr_2017->SetLineColor(kRed);
    //grMC_pre_fsr_2017->Draw("pZ same");

    //grMC_pre_fsr_2018->SetFillColorAlpha(kBlack,0.3);
    //grMC_pre_fsr_2018->Draw("E3 same");

    grMC_pre_fsr_2018->SetMarkerColor(kRed);
    grMC_pre_fsr_2018->SetMarkerStyle(marker_);
    grMC_pre_fsr_2018->SetMarkerSize(marker_size);
    grMC_pre_fsr_2018->SetLineStyle(1);
    grMC_pre_fsr_2018->SetLineColor(kRed);
    //grMC_pre_fsr_2018->Draw("pZ same");

    TGraphErrors *grMC_pre_fsr_2016 = new TGraphErrors(5, &meanMass_mc_pre_fsr[0], &meanPt_mc_pre_fsr[0], &meanMassSysErr_mc_pre_fsr[0], &meanPtSysErr_mc_pre_fsr[0]);
    //grMC_preFSR->SetFillColorAlpha(kBlack,0.3);
    //grMC_preFSR->Draw("E3 same");

    grMC_pre_fsr_2016->SetMarkerColor(kRed);
    grMC_pre_fsr_2016->SetMarkerStyle(20);
    grMC_pre_fsr_2016->SetMarkerSize(marker_size);
    grMC_pre_fsr_2016->SetLineStyle(1);
    grMC_pre_fsr_2016->SetLineColor(kRed);
    //grMC_pre_fsr_2016->Draw("pZ same");

    TGraphErrors *grMC = new TGraphErrors(18, &mean_mass[0], &mean_pt[0], &mean_mass_stat_error[0], &mean_pt_stat_error[0]);
    grMC->SetLineColor(kRed);
    grMC->SetMarkerColor(kRed);
    grMC->SetMarkerStyle(marker_+4);
    grMC->SetMarkerSize(marker_size);
    grMC->SetLineStyle(1);
    //grMC->Draw("pZ same");

    TGraphErrors *grMC_100cut = new TGraphErrors(18, &mean_mass_100cut[0], &mean_pt_100cut[0], &mean_mass_100cut_stat_error[0], &mean_pt_100cut_stat_error[0]);
    grMC_100cut->SetLineColor(kRed);
    grMC_100cut->SetMarkerColor(kRed);
    grMC_100cut->SetMarkerStyle(marker_+4);
    grMC_100cut->SetMarkerSize(marker_size);
    grMC_100cut->SetLineStyle(1);
    //grMC_100cut->Draw("pZ same");

    TGraphErrors *grMC_100cut_err = new TGraphErrors(17, &mean_mass_100cut[1], &mean_pt_100cut[1], &mean_mass_100cut_sys_error[1], &mean_pt_100cut_sys_error[1]);
    grMC_100cut_err->SetFillColorAlpha(kRed,0.3);
    //grMC_100cut_err->Draw("E4 same");

    TGraphErrors *grMC_100cut_MG = new TGraphErrors(18, &mean_mass_100cut_MG[0], &mean_pt_100cut_MG[0], &mean_mass_100cut_stat_error_MG[0], &mean_pt_100cut_stat_error_MG[0]);
    grMC_100cut_MG->SetLineColor(kMagenta);
    grMC_100cut_MG->SetMarkerColor(kMagenta);
    grMC_100cut_MG->SetMarkerStyle(marker_+4);
    grMC_100cut_MG->SetMarkerSize(marker_size);
    grMC_100cut_MG->SetLineStyle(1);
    //grMC_100cut_MG->Draw("pZ same");

    TGraphErrors *grMC_100cut_err_MG = new TGraphErrors(18, &mean_mass_100cut_MG[0], &mean_pt_100cut_MG[0], &mean_mass_100cut_sys_error_MG[0], &mean_pt_100cut_sys_error_MG[0]);
    grMC_100cut_err_MG->SetFillColorAlpha(kMagenta,0.3);
    //grMC_100cut_err_MG->Draw("E4 same");

    TGraphErrors *grMC_100cut_8tev = new TGraphErrors(18, &mean_mass_100cut_8tev[0], &mean_pt_100cut_8tev[0], &mean_mass_100cut_stat_error_8tev[0], &mean_pt_100cut_stat_error_8tev[0]);
    grMC_100cut_8tev->SetLineColor(kBlue);
    grMC_100cut_8tev->SetMarkerColor(kBlue);
    grMC_100cut_8tev->SetMarkerStyle(marker_);
    grMC_100cut_8tev->SetMarkerSize(marker_size);
    grMC_100cut_8tev->SetLineStyle(1);
    //grMC_100cut_8tev->Draw("plZ same");

    // grUnfolded_pre_fsr_electron_run2
    TF1 *f_electron = NULL;
    f_electron = new TF1("f_electron", "[0]+[1]*log(x)", 42., 300.);
    f_electron->GetXaxis()->SetRangeUser(42., 300.);
    f_electron->SetLineColor(kBlack);
    f_electron->SetLineWidth(2);
    grUnfolded_pre_fsr_electron_run2->Fit(f_electron, "R0"); // R: fitting sub range
    //f_electron->Draw("same");

    TF1 *f_muon = NULL;
    f_muon = new TF1("f_muon", "[0]+[1]*log(x)", 42., 300.);
    f_muon->GetXaxis()->SetRangeUser(42., 300.);
    f_muon->SetLineColor(kRed);
    f_muon->SetLineWidth(2);
    grUnfolded_pre_fsr_muon_run2->Fit(f_muon, "R0"); // R: fitting sub range
    //f_muon->Draw("same");
    //
    TF1 *f_lepton = NULL;
    f_lepton = new TF1("f_lepton", "[0]+[1]*log(x)", 42., 300.);
    f_lepton->GetXaxis()->SetRangeUser(42., 300.);
    f_lepton->SetLineColor(kBlue);
    f_lepton->SetLineWidth(2);
    grUnfolded_pre_fsr_lepton_run2->Fit(f_lepton, "R0"); // R: fitting sub range
    f_lepton->Draw("same");

    double chi2_nom = f_lepton->GetChisquare()/ f_lepton->GetNDF();
    cout << "NDF: " << f_lepton->GetNDF() << endl;
    TLegend* leg_ = new TLegend(0.2, 0.55, 0.5, 0.75,"","brNDC");
    leg_->SetNColumns(2);
    leg_->SetTextSize(0.027);
    leg_->SetFillStyle(0); // transparent
    leg_->SetBorderSize(0);
    //leg_->AddEntry(grMC_pre_fsr_2016, "DY MC aMC@NLO", "pe");
    //leg_->AddEntry(grMC_100cut_MG, "DY MC Madgraph (LO)", "pe");
    //leg_->AddEntry(grMC_100cut_8tev, "DY MC Madgraph (8TeV, LO)", "pe");
    //leg_->AddEntry(grMC_pre_fsr_2017, "2017,2018 DY MC", "pe");
    leg_->AddEntry(grUnfolded_pre_fsr, "2016 electron", "pe");
    leg_->AddEntry(grUnfolded_pre_fsr_muon_2016, "2016 muon", "pe");
    leg_->AddEntry(grUnfolded_pre_fsr_2017, "2017 electron ", "pe");
    leg_->AddEntry(grUnfolded_pre_fsr_muon_2017, "2017 muon ", "pe");
    leg_->AddEntry(grUnfolded_pre_fsr_2018, "2018 electron ", "pe");
    leg_->AddEntry(grUnfolded_pre_fsr_muon_2018, "2018 muon ", "pe");
    leg_->Draw();

    TLegend* leg_fit = new TLegend(0.65, 0.35, 0.9, 0.55,"","brNDC");
    leg_fit->SetNColumns(1);
    leg_fit->SetTextSize(0.027);
    leg_fit->SetFillStyle(0); // transparent
    leg_fit->SetBorderSize(0);
    leg_fit->AddEntry(f_lepton, Form("\\mbox{Linear fit } (\\chi^{2}/NDF=%.2f)", chi2_nom), "l");
    leg_fit->Draw("same");

    TF1* f1;
    TF1* f2;
    TF1* f3;
    TF1* f4;
    TF1* f5;
    TF1* f6;

    if(doFit){
        f1 = new TF1("f1", "[0]+[1]*log(x)", 42., 300.);
        f1->GetXaxis()->SetRangeUser(42., 300.);
        f1->SetLineColor(kBlack);
        f1->SetLineWidth(1);
        grUnfolded_pre_fsr->Fit(f1, "R0"); // R: fitting sub range
        f1->Draw("same");

        f2 = new TF1("f2", "[0]+[1]*log(x)", 42., 300.);
        f2->GetXaxis()->SetRangeUser(42., 300.);
        f2->SetLineColor(kRed);
        f2->SetLineWidth(1);
        grUnfolded_pre_fsr_2017->Fit(f2, "R0"); // R: fitting sub range
        f2->Draw("same");

        f3 = new TF1("f3", "[0]+[1]*log(x)", 42., 300.);
        f3->GetXaxis()->SetRangeUser(42., 300.);
        f3->SetLineColor(kBlue);
        f3->SetLineWidth(1);
        grUnfolded_pre_fsr_2018->Fit(f3, "R0"); // R: fitting sub range
        f3->Draw("same");

        f4 = new TF1("f4", "[0]+[1]*log(x)", 42., 300.);
        f4->GetXaxis()->SetRangeUser(42., 300.);
        f4->SetLineColor(kBlack);
        f4->SetLineWidth(1);
        grUnfolded_pre_fsr_muon_2016->Fit(f4, "R0"); // R: fitting sub range
        f4->Draw("same");

        f5 = new TF1("f5", "[0]+[1]*log(x)", 42., 300.);
        f5->GetXaxis()->SetRangeUser(42., 300.);
        f5->SetLineColor(kRed);
        f5->SetLineWidth(1);
        grUnfolded_pre_fsr_muon_2017->Fit(f5, "R0"); // R: fitting sub range
        f5->Draw("same");

        f6 = new TF1("f6", "[0]+[1]*log(x)", 42., 300.);
        f6->GetXaxis()->SetRangeUser(42., 300.);
        f6->SetLineColor(kBlue);
        f6->SetLineWidth(1);
        grUnfolded_pre_fsr_muon_2018->Fit(f6, "R0"); // R: fitting sub range
        f6->Draw("same");
    }

    filein_dy->Close();

    CMS_lumi( c1, 5, 10 );
    //c1->SaveAs(outpdf + "_" + year_string + "_nofit.png");
    c1->SaveAs(outpdf + "_" + year_string + "_fit.png");
    //delete grUnfolded;
    delete f1;
    delete c1;

    c1 = new TCanvas("c_Run2_detector", "c_Run2_detector", 50, 50, 900*1.5, 700*1.5);
    c1->cd();
    gStyle->SetOptFit(0);

    c1->SetBottomMargin(0.2);
    c1->SetTopMargin(0.08);
    //c1->SetTicks();
    c1->SetLogx();
    c1->SetGridy();
    c1->SetGridx();

    TGraphErrors* grDetector_2017 = (TGraphErrors*)c_2017->GetListOfPrimitives()->FindObject("detectorData_electron_2017");
    TGraphErrors* grDetector_2018 = (TGraphErrors*)c_2018->GetListOfPrimitives()->FindObject("detectorData_electron_2018");

    TGraphErrors* grDetector_muon_2016 = (TGraphErrors*)c_muon_2016->GetListOfPrimitives()->FindObject("detectorData_muon_2016");
    TGraphErrors* grDetector_muon_2017 = (TGraphErrors*)c_muon_2017->GetListOfPrimitives()->FindObject("detectorData_muon_2017");
    TGraphErrors* grDetector_muon_2018 = (TGraphErrors*)c_muon_2018->GetListOfPrimitives()->FindObject("detectorData_muon_2018");

    TGraphErrors *grDetector = new TGraphErrors(5, &meanMass_data_detector[0], &meanPt_data_detector[0], &meanMassStatErr_data_detector[0], &meanPtStatErr_data_detector[0]);
    grDetector->SetLineColor(kBlack);
    grDetector->SetMarkerColor(kBlack);
    grDetector->SetMarkerStyle(marker_);
    grDetector->SetMarkerSize(marker_size);
    grDetector->SetLineStyle(1);
    grDetector->Draw("apZ");
    grDetector->GetYaxis()->SetRangeUser(11.,32.);
    grDetector->GetXaxis()->SetLimits(30.,500.);
    grDetector->GetXaxis()->SetMoreLogLabels(true);
    grDetector->GetYaxis()->SetTitle("\\mbox{Average } p_{T}^{\\ell\\ell} \\mbox{ (GeV)}"); // \\ell not working with pdf
    grDetector->GetYaxis()->SetLabelSize(23*2);
    grDetector->GetYaxis()->SetTitleSize(23*2);
    grDetector->GetYaxis()->SetTitleOffset(1.5);
    grDetector->GetXaxis()->SetTitle("\\mbox{Average } Mass^{\\ell\\ell} \\mbox{ (GeV)}");
    grDetector->GetXaxis()->SetLabelSize(23*2);
    grDetector->GetXaxis()->SetTitleSize(23*2);
    grDetector->GetXaxis()->SetTitleOffset(2.);

    grDetector_2017->SetLineColor(kBlack);
    grDetector_2017->SetMarkerColor(kBlack);
    grDetector_2017->SetMarkerStyle(marker_+5);
    grDetector_2017->SetMarkerSize(marker_size);
    grDetector_2017->SetLineStyle(1);
    grDetector_2017->Draw("pZ same");

    grDetector_2018->SetLineColor(kBlack);
    grDetector_2018->SetMarkerColor(kBlack);
    grDetector_2018->SetMarkerStyle(marker_+6);
    grDetector_2018->SetMarkerSize(marker_size);
    grDetector_2018->SetLineStyle(1);
    grDetector_2018->Draw("pZ same");

    grDetector_muon_2016->SetLineColor(kRed);
    grDetector_muon_2016->SetMarkerColor(kRed);
    grDetector_muon_2016->SetMarkerStyle(marker_);
    grDetector_muon_2016->SetMarkerSize(marker_size);
    grDetector_muon_2016->SetLineStyle(1);
    grDetector_muon_2016->Draw("pZ same");

    grDetector_muon_2017->SetLineColor(kRed);
    grDetector_muon_2017->SetMarkerColor(kRed);
    grDetector_muon_2017->SetMarkerStyle(marker_+ 5);
    grDetector_muon_2017->SetMarkerSize(marker_size);
    grDetector_muon_2017->SetLineStyle(1);
    grDetector_muon_2017->Draw("pZ same");

    grDetector_muon_2018->SetLineColor(kBlue);
    grDetector_muon_2018->SetMarkerColor(kBlue);
    grDetector_muon_2018->SetMarkerStyle(marker_+ 6);
    grDetector_muon_2018->SetMarkerSize(marker_size);
    grDetector_muon_2018->SetLineStyle(1);
    grDetector_muon_2018->Draw("pZ same");

    TLegend* leg_detector = new TLegend(0.2, 0.55, 0.5, 0.75,"","brNDC");
    leg_detector->SetNColumns(2);
    leg_detector->SetTextSize(0.027);
    leg_detector->SetFillStyle(0); // transparent
    leg_detector->SetBorderSize(0);
    //leg_detector->AddEntry(grMC_pre_fsr_2016, "DY MC aMC@NLO", "pe");
    //leg_detector->AddEntry(grMC_100cut_MG, "DY MC Madgraph (LO)", "pe");
    //leg_detector->AddEntry(grMC_100cut_8tev, "DY MC Madgraph (8TeV, LO)", "pe");
    //leg_detector->AddEntry(grMC_pre_fsr_2017, "2017,2018 DY MC", "pe");
    leg_detector->AddEntry(grDetector, "2016 electron", "pe");
    leg_detector->AddEntry(grDetector_muon_2016, "2016 muon", "pe");
    leg_detector->AddEntry(grDetector_2017, "2017 electron ", "pe");
    leg_detector->AddEntry(grDetector_muon_2017, "2017 muon ", "pe");
    leg_detector->AddEntry(grDetector_2018, "2018 electron ", "pe");
    leg_detector->AddEntry(grDetector_muon_2018, "2018 muon ", "pe");
    leg_detector->Draw();

    CMS_lumi( c1, 5, 10 );
    c1->SaveAs(outpdf + "Detector_" + year_string + ".png");
    //delete grUnfolded;
    delete c1;

    c1 = new TCanvas("c_Run2_postFSR", "c_Run2_postFSR", 50, 50, 900*1.5, 700*1.5);
    c1->cd();
    gStyle->SetOptFit(0);

    c1->SetBottomMargin(0.2);
    c1->SetTopMargin(0.08);
    //c1->SetTicks();
    c1->SetLogx();
    c1->SetGridy();
    c1->SetGridx();

    TGraphErrors* grUnfolded_post_fsr_2017 = (TGraphErrors*)c_2017->GetListOfPrimitives()->FindObject("postFSRUnfoldedData_electron_2017");
    TGraphErrors* grUnfolded_post_fsr_2018 = (TGraphErrors*)c_2018->GetListOfPrimitives()->FindObject("postFSRUnfoldedData_electron_2018");

    TGraphErrors* grUnfolded_post_fsr_muon_2016 = (TGraphErrors*)c_muon_2016->GetListOfPrimitives()->FindObject("postFSRUnfoldedData_muon_2016");
    TGraphErrors* grUnfolded_post_fsr_muon_2017 = (TGraphErrors*)c_muon_2017->GetListOfPrimitives()->FindObject("postFSRUnfoldedData_muon_2017");
    TGraphErrors* grUnfolded_post_fsr_muon_2018 = (TGraphErrors*)c_muon_2018->GetListOfPrimitives()->FindObject("postFSRUnfoldedData_muon_2018");

    TGraphErrors *grUnfolded = new TGraphErrors(5, &meanMass_data_det_unf[0], &meanPt_data_det_unf[0], &meanMassTotErr_data_det_unf[0], &meanPtTotErr_data_det_unf[0]);
    grUnfolded->SetLineColor(kBlack);
    grUnfolded->SetMarkerColor(kBlack);
    grUnfolded->SetMarkerStyle(marker_);
    grUnfolded->SetMarkerSize(marker_size);
    grUnfolded->SetLineStyle(1);
    grUnfolded->Draw("apZ");
    grUnfolded->GetYaxis()->SetRangeUser(11.,32.);
    grUnfolded->GetXaxis()->SetLimits(30.,500.);
    grUnfolded->GetXaxis()->SetMoreLogLabels(true);
    grUnfolded->GetYaxis()->SetTitle("\\mbox{Average } p_{T}^{\\ell\\ell} \\mbox{ (GeV)}"); // \\ell not working with pdf
    grUnfolded->GetYaxis()->SetLabelSize(23*2);
    grUnfolded->GetYaxis()->SetTitleSize(23*2);
    grUnfolded->GetYaxis()->SetTitleOffset(1.5);
    grUnfolded->GetXaxis()->SetTitle("\\mbox{Average } Mass^{\\ell\\ell} \\mbox{ (GeV)}");
    grUnfolded->GetXaxis()->SetLabelSize(23*2);
    grUnfolded->GetXaxis()->SetTitleSize(23*2);
    grUnfolded->GetXaxis()->SetTitleOffset(2.);

    grUnfolded_post_fsr_2017->SetLineColor(kBlack);
    grUnfolded_post_fsr_2017->SetMarkerColor(kBlack);
    grUnfolded_post_fsr_2017->SetMarkerStyle(marker_+5);
    grUnfolded_post_fsr_2017->SetMarkerSize(marker_size);
    grUnfolded_post_fsr_2017->SetLineStyle(1);
    grUnfolded_post_fsr_2017->Draw("pZ same");

    grUnfolded_post_fsr_2018->SetLineColor(kBlack);
    grUnfolded_post_fsr_2018->SetMarkerColor(kBlack);
    grUnfolded_post_fsr_2018->SetMarkerStyle(marker_+6);
    grUnfolded_post_fsr_2018->SetMarkerSize(marker_size);
    grUnfolded_post_fsr_2018->SetLineStyle(1);
    grUnfolded_post_fsr_2018->Draw("pZ same");

    grUnfolded_post_fsr_muon_2016->SetLineColor(kRed);
    grUnfolded_post_fsr_muon_2016->SetMarkerColor(kRed);
    grUnfolded_post_fsr_muon_2016->SetMarkerStyle(marker_);
    grUnfolded_post_fsr_muon_2016->SetMarkerSize(marker_size);
    grUnfolded_post_fsr_muon_2016->SetLineStyle(1);
    grUnfolded_post_fsr_muon_2016->Draw("pZ same");

    grUnfolded_post_fsr_muon_2017->SetLineColor(kRed);
    grUnfolded_post_fsr_muon_2017->SetMarkerColor(kRed);
    grUnfolded_post_fsr_muon_2017->SetMarkerStyle(marker_+ 5);
    grUnfolded_post_fsr_muon_2017->SetMarkerSize(marker_size);
    grUnfolded_post_fsr_muon_2017->SetLineStyle(1);
    grUnfolded_post_fsr_muon_2017->Draw("pZ same");

    grUnfolded_post_fsr_muon_2018->SetLineColor(kRed);
    grUnfolded_post_fsr_muon_2018->SetMarkerColor(kRed);
    grUnfolded_post_fsr_muon_2018->SetMarkerStyle(marker_+ 6);
    grUnfolded_post_fsr_muon_2018->SetMarkerSize(marker_size);
    grUnfolded_post_fsr_muon_2018->SetLineStyle(1);
    grUnfolded_post_fsr_muon_2018->Draw("pZ same");

    TLegend* leg_unfolded = new TLegend(0.2, 0.55, 0.5, 0.75,"","brNDC");
    leg_unfolded->SetNColumns(2);
    leg_unfolded->SetTextSize(0.027);
    leg_unfolded->SetFillStyle(0); // transparent
    leg_unfolded->SetBorderSize(0);
    //leg_unfolded->AddEntry(grMC_pre_fsr_2016, "DY MC aMC@NLO", "pe");
    //leg_unfolded->AddEntry(grMC_100cut_MG, "DY MC Madgraph (LO)", "pe");
    //leg_unfolded->AddEntry(grMC_100cut_8tev, "DY MC Madgraph (8TeV, LO)", "pe");
    //leg_unfolded->AddEntry(grMC_pre_fsr_2017, "2017,2018 DY MC", "pe");
    leg_unfolded->AddEntry(grUnfolded, "2016 electron", "pe");
    leg_unfolded->AddEntry(grUnfolded_post_fsr_muon_2016, "2016 muon", "pe");
    leg_unfolded->AddEntry(grUnfolded_post_fsr_2017, "2017 electron ", "pe");
    leg_unfolded->AddEntry(grUnfolded_post_fsr_muon_2017, "2017 muon ", "pe");
    leg_unfolded->AddEntry(grUnfolded_post_fsr_2018, "2018 electron ", "pe");
    leg_unfolded->AddEntry(grUnfolded_post_fsr_muon_2018, "2018 muon ", "pe");
    leg_unfolded->Draw();

    CMS_lumi( c1, 5, 10 );
    c1->SaveAs(outpdf + "postFSR_" + year_string + ".png");
    //delete grUnfolded;
    delete c1;
}

TCanvas* ISRUnfold::drawISRresult(TString outpdf, bool doFit, TString genISRFile)
{

    gROOT->SetBatch();

    TFile* filein_dy = NULL;
    TFile* filein_dy_MG = NULL;
    TFile* filein_dy_Powheg = NULL;
    TH1* hmcPt;
    TH1* hmcMass;
    TH1* hmcPt_mg;
    TH1* hmcMass_mg;
    TH1* hmcPt_powheg;
    TH1* hmcMass_powheg;

    vector<double> mean_pt, mean_pt_stat_error, mean_pt_sys_error;
    vector<double> mean_mass, mean_mass_stat_error, mean_mass_sys_error;
    vector<double> mean_pt_mg, mean_pt_stat_error_mg;
    vector<double> mean_mass_mg, mean_mass_stat_error_mg;
    vector<double> mean_pt_powheg, mean_pt_stat_error_powheg;
    vector<double> mean_mass_powheg, mean_mass_stat_error_powheg;

    int points = 5;
    if(genISRFile != "")
    {
        filein_dy = new TFile(genISRFile + "DY_FSR.root");
        filein_dy_MG = new TFile(genISRFile + "DY_FSR_MG.root");
        filein_dy_Powheg = new TFile(genISRFile + "DY_FSR_powheg.root");

        std::map<TString, int> sys_map;
        sys_map["AlphaS"] = 2;
        sys_map["Scale"] = 6;
        sys_map["PDFerror"] = 100;

        int lowMassEdge = 35;
        int highMassEdge = 45;

        TString massLowEdge[5] = {"50.", "65.", "80.", "100.", "200."};
        TString massHighEdge[5] = {"65.", "80.", "100.", "200.", "350."};

        if(channel_name=="muon")
        {
            massLowEdge[0] = "40.";
            massLowEdge[1] = "60.";

            massHighEdge[0] = "60.";
            massHighEdge[1] = "80.";
        }

        for(int i = 0; i < points; i++)
        {
            int lowMassEdge_ = lowMassEdge + 5 * i;
            int highMassEdge_ = highMassEdge + 5 * i;
            TString s_lowMassEdge;;
            TString s_highMassEdge;;
            s_lowMassEdge.Form("%d", lowMassEdge_);
            s_highMassEdge.Form("%d", highMassEdge_);

            //TString dir_name = channel_name + "_full_phase_m" + s_lowMassEdge + "to" + s_highMassEdge;
            TString dir_name = "full_phase_m" + massLowEdge[i] + "to" + massHighEdge[i];
            cout << dir_name << endl;

            TString phase_spase = "";

            hmcPt = (TH1*)filein_dy->Get(dir_name + "/ptll" + phase_spase + "/histo_DYJets");
            hmcPt->Add((TH1*)filein_dy->Get(dir_name + "/ptll" + phase_spase + "/histo_DYJets10to50"));

            hmcMass = (TH1*)filein_dy->Get(dir_name + "/mll" + phase_spase + "/histo_DYJets");
            hmcMass->Add((TH1*)filein_dy->Get(dir_name + "/mll" + phase_spase + "/histo_DYJets10to50"));

            mean_pt.push_back(hmcPt->GetMean());
            mean_pt_stat_error.push_back(hmcPt->GetMeanError());

            mean_mass.push_back(hmcMass->GetMean());
            mean_mass_stat_error.push_back(hmcMass->GetMeanError());

            delete hmcPt, hmcMass;

            hmcPt_mg = (TH1*)filein_dy_MG->Get(dir_name + "/ptll" + phase_spase + "/histo_DYJets_MG");
            hmcPt_mg->Add((TH1*)filein_dy_MG->Get(dir_name + "/ptll" + phase_spase + "/histo_DYJets10to50_MG"));

            hmcMass_mg = (TH1*)filein_dy_MG->Get(dir_name + "/mll" + phase_spase + "/histo_DYJets_MG");
            hmcMass_mg->Add((TH1*)filein_dy_MG->Get(dir_name + "/mll" + phase_spase + "/histo_DYJets10to50_MG"));

            mean_pt_mg.push_back(hmcPt_mg->GetMean());
            mean_pt_stat_error_mg.push_back(hmcPt_mg->GetMeanError());

            mean_mass_mg.push_back(hmcMass_mg->GetMean());
            mean_mass_stat_error_mg.push_back(hmcMass_mg->GetMeanError());

            delete hmcPt_mg, hmcMass_mg;

            // calculate uncertainty
            double tot_sys_pt_err = 0;
            double tot_sys_mass_err = 0;
            double tot_sys_pt_err_mg = 0;
            double tot_sys_mass_err_mg = 0;
            double tot_sys_pt_err_powheg = 0;
            double tot_sys_mass_err_powheg = 0;
            std::map<TString, int>::iterator it = sys_map.begin();
            while(it != sys_map.end())
            {
                double sys_pt_err = -999.;
                double sys_mass_err = -999.;
                double sys_pt_err_mg = -999.;
                double sys_mass_err_mg = -999.;
                double sys_pt_err_powheg = -999.;
                double sys_mass_err_powheg = -999.;

                for(int ith_sys = 0; ith_sys < it->second; ith_sys++)
                {
                    TString systematic_postfix = it->first;

                    double temp_sys_pt_err = -999.;
                    double temp_sys_mass_err = -999.;
                    double temp_sys_pt_err_mg = -999.;
                    double temp_sys_mass_err_mg = -999.;
                    double temp_sys_pt_err_powheg = -999.;
                    double temp_sys_mass_err_powheg = -999.;
                    if(it->first == "AlphaS")
                    {
                        if(ith_sys == 0) systematic_postfix+="Up";
                        if(ith_sys == 1) systematic_postfix+="Down";
                    }
                    if(it->first == "Scale")
                    {
                        if(ith_sys == 0 ) systematic_postfix+="AUp";
                        if(ith_sys == 1 ) systematic_postfix+="ADown";
                        if(ith_sys == 2 ) systematic_postfix+="BUp";
                        if(ith_sys == 3 ) systematic_postfix+="BDown";
                        if(ith_sys == 4 ) systematic_postfix+="ABUp";
                        if(ith_sys == 5 ) systematic_postfix+="ABDown";
                    }
                    if(it->first  == "PDFerror")
                    {
                        TString nth_;
                        nth_.Form ("%03d", ith_sys);
                        systematic_postfix+=nth_;
                    }

                    hmcPt = (TH1*)filein_dy->Get(dir_name + "/ptll" + phase_spase + "/histo_DYJets_"+systematic_postfix);
                    hmcPt->Add((TH1*)filein_dy->Get(dir_name + "/ptll" + phase_spase + "/histo_DYJets10to50_"+systematic_postfix));

                    hmcMass = (TH1*)filein_dy->Get(dir_name + "/mll" + phase_spase + "/histo_DYJets_"+systematic_postfix);
                    hmcMass->Add((TH1*)filein_dy->Get(dir_name + "/mll" + phase_spase + "/histo_DYJets10to50_"+systematic_postfix));

                    temp_sys_pt_err = fabs(hmcPt->GetMean()-mean_pt.at(i));
                    temp_sys_mass_err = fabs(hmcMass->GetMean()-mean_mass.at(i));
                    delete hmcPt, hmcMass;

                    if(sys_pt_err < temp_sys_pt_err)
                    {
                        sys_pt_err = temp_sys_pt_err;
                    }
                    if(sys_mass_err < temp_sys_mass_err)
                    {
                        sys_mass_err = temp_sys_mass_err;
                    }

                }// loop over variations

                tot_sys_pt_err += pow(sys_pt_err,2);
                tot_sys_mass_err += pow(sys_mass_err,2);
                it++;
            }// loop over systematic sources
            mean_pt_sys_error.push_back(sqrt(tot_sys_pt_err));
            mean_mass_sys_error.push_back(sqrt(tot_sys_mass_err));

            TString channel_ = "EE";
            if(channel_name == "muon")
            {
                channel_ = "MuMu";
            }
            hmcPt_powheg = (TH1*)filein_dy_Powheg->Get(dir_name + "/ptll" + phase_spase + "/histo_ZTo" + channel_ + "_M_50_120_powheg");
            hmcPt_powheg->Add((TH1*)filein_dy_Powheg->Get(dir_name + "/ptll" + phase_spase + "/histo_ZTo" + channel_ + "_M_120_200_powheg"));
            hmcPt_powheg->Add((TH1*)filein_dy_Powheg->Get(dir_name + "/ptll" + phase_spase + "/histo_ZTo" + channel_ + "_M_200_400_powheg"));

            hmcMass_powheg = (TH1*)filein_dy_Powheg->Get(dir_name + "/mll" + phase_spase + "/histo_ZTo" + channel_ + "_M_50_120_powheg");
            hmcMass_powheg->Add((TH1*)filein_dy_Powheg->Get(dir_name + "/mll" + phase_spase + "/histo_ZTo" + channel_ + "_M_120_200_powheg"));
            hmcMass_powheg->Add((TH1*)filein_dy_Powheg->Get(dir_name + "/mll" + phase_spase + "/histo_ZTo" + channel_ + "_M_200_400_powheg"));

            mean_pt_powheg.push_back(hmcPt_powheg->GetMean());
            mean_pt_stat_error_powheg.push_back(hmcPt_powheg->GetMeanError());

            mean_mass_powheg.push_back(hmcMass_powheg->GetMean());
            mean_mass_stat_error_powheg.push_back(hmcMass_powheg->GetMeanError());

            delete hmcPt_powheg, hmcMass_powheg;
        }
    }

    int marker_ = 20;
    if(channel_name=="muon") marker_ = 20;

    setTDRStyle();
    writeExtraText = true;       // if extra text
    extraText  = "work in progress";

    TString year_string;
    year_string.Form ("%d", year);

    c1 = new TCanvas("c_"+channel_name+"_"+year_string, "c_"+channel_name+"_"+year_string, 50, 50, 900, 850);
    c1->cd();
    gStyle->SetOptFit(0);

    c1->SetBottomMargin(0.2);
    c1->SetTopMargin(0.08);
    //c1->SetTicks();
    c1->SetLogx();
    c1->SetGridy();
    c1->SetGridx();

    TString lepton_type;
    if(channel_name=="electron") lepton_type = "ee";
    else lepton_type = "#mu#mu";

    // unfolded data
    TGraphErrors *grUnfolded_pre_fsr = new TGraphErrors(5, &meanMass_data_pre_fsr[0], &meanPt_data_pre_fsr[0], &meanMassTotErr_data_pre_fsr[0], &meanPtTotErr_data_pre_fsr[0]);
    grUnfolded_pre_fsr->SetLineColor(kBlack);
    grUnfolded_pre_fsr->SetMarkerColor(kBlack);
    grUnfolded_pre_fsr->SetMarkerStyle(marker_);
    grUnfolded_pre_fsr->SetMarkerSize(1.2);
    grUnfolded_pre_fsr->SetLineStyle(1);
    grUnfolded_pre_fsr->GetYaxis()->SetRangeUser(12.,28.);
    grUnfolded_pre_fsr->GetXaxis()->SetLimits(40.,400.);
    grUnfolded_pre_fsr->GetXaxis()->SetMoreLogLabels(true);
    grUnfolded_pre_fsr->GetYaxis()->SetTitle("Average p_{T} (" + lepton_type + ") (GeV)");
    grUnfolded_pre_fsr->GetYaxis()->SetTitleSize(40);
    grUnfolded_pre_fsr->GetYaxis()->SetTitleOffset(1.2);
    grUnfolded_pre_fsr->GetYaxis()->SetLabelSize(40);
    grUnfolded_pre_fsr->GetXaxis()->SetLabelSize(40);
    grUnfolded_pre_fsr->GetXaxis()->SetTitle("Average Mass (" + lepton_type + ") (GeV)");
    grUnfolded_pre_fsr->GetXaxis()->SetTitleSize(40);
    grUnfolded_pre_fsr->GetXaxis()->SetTitleOffset(1.5);
    grUnfolded_pre_fsr->GetYaxis()->SetNdivisions(510);
    grUnfolded_pre_fsr->Draw("ap");
    grUnfolded_pre_fsr->SetName("preFSRUnfoldedData_"+channel_name+"_"+year_string);

    //TGraphErrors *grMC_preFSR = new TGraphErrors(points, &mean_mass[0], &mean_pt[0], &mean_mass_stat_error[0], &mean_pt_stat_error[0]);
    TGraphErrors *grMC_preFSR = new TGraphErrors(5, &meanMass_mc_pre_fsr[0], &meanPt_mc_pre_fsr[0], &mean_mass_stat_error[0], &mean_pt_stat_error[0]);
    grMC_preFSR->SetMarkerColor(kRed);
    grMC_preFSR->SetMarkerStyle(marker_+4);
    grMC_preFSR->SetMarkerSize(1.5);
    grMC_preFSR->SetLineStyle(1);
    grMC_preFSR->SetLineColor(kRed);
    grMC_preFSR->Draw("p same");
    //
    //grMC_preFSR->SetFillColor(kBlack);
    //grMC_preFSR->SetFillStyle(3001);
    //grMC_preFSR->SetFillColorAlpha(kBlack,0.5);
    //grMC_preFSR->Draw("E3 same");
    grMC_preFSR->SetName("preFSRMC_"+channel_name+"_"+year_string);

    TGraphErrors *grMC_preFSR_MG = new TGraphErrors(points, &mean_mass_mg[0], &mean_pt_mg[0], &mean_mass_stat_error_mg[0], &mean_pt_stat_error_mg[0]);
    grMC_preFSR_MG->SetMarkerColor(kMagenta);
    grMC_preFSR_MG->SetMarkerStyle(marker_+5);
    grMC_preFSR_MG->SetMarkerSize(1.5);
    grMC_preFSR_MG->SetLineStyle(1);
    grMC_preFSR_MG->SetLineColor(kMagenta);
    grMC_preFSR_MG->Draw("p same");

    TGraphErrors *grMC_preFSR_Powheg = new TGraphErrors(points, &mean_mass_powheg[0], &mean_pt_powheg[0], &mean_mass_stat_error_powheg[0], &mean_pt_stat_error_powheg[0]);
    grMC_preFSR_Powheg->SetMarkerColor(kBlue);
    grMC_preFSR_Powheg->SetMarkerStyle(marker_+6);
    grMC_preFSR_Powheg->SetMarkerSize(1.5);
    grMC_preFSR_Powheg->SetLineStyle(1);
    grMC_preFSR_Powheg->SetLineColor(kBlue);
    grMC_preFSR_Powheg->Draw("p same");

    // detector unfolded data
    TGraphErrors *grUnfolded = new TGraphErrors(5, &meanMass_data_det_unf[0], &meanPt_data_det_unf[0], &meanMassTotErr_data_det_unf[0], &meanPtTotErr_data_det_unf[0]);
    grUnfolded->SetLineColor(kGray+1);
    grUnfolded->SetMarkerColor(kGray+1);
    grUnfolded->SetMarkerStyle(marker_);
    grUnfolded->SetMarkerSize(0.9);
    grUnfolded->SetLineStyle(1);
    //grUnfolded->Draw("pZ same");
    grUnfolded->SetName("postFSRUnfoldedData_"+channel_name+"_"+year_string);

    TGraphErrors *grMC = new TGraphErrors(5, &meanMass_mc_det_unf[0], &meanPt_mc_det_unf[0], &meanMassErr_mc_det_unf[0], &meanPtErr_mc_det_unf[0]);
    grMC->SetLineColor(kRed);
    grMC->SetMarkerColor(kRed);
    grMC->SetMarkerStyle(marker_);
    grMC->SetMarkerSize(.9);
    grMC->SetLineStyle(1);
    grMC->SetLineColor(kRed);
    //grMC->Draw("pZ same");
    grMC->SetName("postFSRMC_"+channel_name+"_"+year_string);


    // detector level data
    TGraphErrors *grDetector = new TGraphErrors(5, &meanMass_data_detector[0], &meanPt_data_detector[0], &meanMassStatErr_data_detector[0], &meanPtStatErr_data_detector[0]);
    grDetector->SetLineColor(kGray+2);
    grDetector->SetMarkerColor(kGray+2);
    grDetector->SetMarkerStyle(marker_);
    grDetector->SetMarkerSize(0.9);
    grDetector->SetLineStyle(1);
    //grDetector->Draw("pZ same");
    grDetector->SetName("detectorData_"+channel_name+"_"+year_string);


    TF1 *f1 = NULL;
    TF1 *f2 = NULL;
    double chi2_nom_data = 0.;
    if(doFit)
    {

        f1 = new TF1("f1", "[0]+[1]*log(x)", 40., 400.);
        f1->GetXaxis()->SetRangeUser(30., 400.);
        f1->SetLineColor(kBlack);
        f1->SetLineWidth(1);
        //grUnfolded->Fit(f1, "R0"); // R: fitting sub range
        grUnfolded_pre_fsr->Fit(f1, "R0"); // R: fitting sub range
        f1->Draw("same");
        chi2_nom_data = f1->GetChisquare()/ f1->GetNDF();

        f2 = new TF1("f2", "[0]+[1]*log(x)", 20., 500.);
        f2->GetXaxis()->SetRangeUser(20., 500.);
        f2->SetLineColor(kRed);
        f2->SetLineWidth(1);
        f2->SetLineStyle(2);
        grMC_preFSR->Fit(f2, "R0"); // R: fitting sub range
        //f2->Draw("same");

    }
    grUnfolded_pre_fsr->Draw("p same");

    TLegend* leg_ = new TLegend(0.15, 0.5, 0.5, 0.9,"","brNDC");
    leg_->SetTextSize(0.035);
    leg_->SetFillStyle(0); // transparent
    leg_->SetBorderSize(0);
    leg_->AddEntry(grUnfolded_pre_fsr, "Unfolded data ", "pe");
    leg_->AddEntry(f1, Form("\\mbox{Linear fit } (\\chi^{2}/NDF=%.2f)", chi2_nom_data), "l");
    leg_->Draw();
    TLegend* leg_mc = new TLegend(0.6, 0.25, 0.8, 0.55,"","brNDC");
    leg_mc->SetTextSize(0.035);
    leg_mc->SetFillStyle(0); // transparent
    leg_mc->SetBorderSize(0);
    leg_mc->AddEntry(grMC_preFSR, "Drell-Yan(DY) aMC@NLO", "pe");
    leg_mc->AddEntry(grMC_preFSR_MG, "DY Madgraph", "pe");
    leg_mc->AddEntry(grMC_preFSR_Powheg, "DY Powheg", "pe");
    leg_mc->Draw();

    CMS_lumi( c1, 4, 0 );
    c1->SaveAs(outpdf + channel_name + "_" + year_string + ".pdf");
    //delete grUnfolded;
    delete grMC;
    delete f1;
    //delete c1;
    //
    mean_pt.clear();
    mean_pt_stat_error.clear();
    mean_mass.clear();
    mean_mass_stat_error.clear();
    return c1;
}

void ISRUnfold::SavePtMassHists()
{
    gROOT->SetBatch();

    TFile* p_fout = NULL;

    TString year_string;
    year_string.Form("%d", year);

    TH1* p_hist_detector;
    TH1* p_hist_detector_unfold;
    TH1* p_hist_fsr_unfold;

    p_fout = new TFile(year_string + "_" + channel_name + "_histograms.root","recreate");
    p_fout->cd();

    TString ibin_string;
    for(int ibin = 0; ibin < 5; ibin++)
    {
        ibin_string.Form("%d", ibin);

        p_hist_detector = nomPtUnfold->GetInput("h_detector_Pt_MassBin_"+ibin_string,0,0,"pt[UO];mass[UOC"+ibin_string+"]",kTRUE);
        p_hist_detector_unfold = nomPtUnfold->GetOutput("h_detector_unfold_Pt_MassBin_"+ibin_string,0,0,"pt[UO];mass[UOC"+ibin_string+"]",kTRUE);
        p_hist_fsr_unfold = nomPtFSRUnfold->GetOutput("h_fsr_unfold_Pt_MassBin_"+ibin_string,0,0,"pt[UO];mass[UOC"+ibin_string+"]",kTRUE);

        //p_hist_detector_unfold = sysPtUnfold["unfoldScan"].at(0)->GetOutput("h_detector_unfold_Pt_MassBin_"+ibin_string,0,0,"pt[UO];mass[UOC"+ibin_string+"]",kTRUE);
        //p_hist_fsr_unfold = sysPtFSRUnfold["unfoldScan"].at(0)->GetOutput("h_fsr_unfold_Pt_MassBin_"+ibin_string,0,0,"pt[UO];mass[UOC"+ibin_string+"]",kTRUE);

        p_hist_detector->Write();
        p_hist_detector_unfold->Write();
        p_hist_fsr_unfold->Write();
        delete p_hist_detector;
        delete p_hist_detector_unfold;
        delete p_hist_fsr_unfold;
    }


    p_fout->cd();
    p_fout->Write();
    p_fout->Close();
}

void ISRUnfold::drawInputPlots(TString outpdf, TString var, int nthMassBin, TString sysName)
{

    cout << " ISRUnfold::drawInputPlots " << endl;
    cout << "sysName: " << sysName << endl;
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
    for(int i = 0; i < sysSize; i++)
    {

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

    for(int i = 0; i < sysSize; i++)
    {
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
    delete leg_nom;
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

void ISRUnfold::drawSysSummaryPlots(TString outpdf, bool detector_unfold)
{
    TString var = "Pt";

    bool relative = true;
    typedef std::map<TString, Double_t> mapdata;
    const TUnfoldBinningV17* temp_binning_gen_pt = nomPtUnfold->GetOutputBinning("Gen_Pt");
    const TVectorD* temp_tvecd = temp_binning_gen_pt->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;
    int nPtBin = temp_binning_gen_pt->GetDistributionBinning(0)->GetNrows() - 1;

    gROOT->SetBatch();

    setTDRStyle();
    writeExtraText = true; 
    extraText  = "work in progress";

    TString year_string;
    year_string.Form ("%d", year);

    c1 = new TCanvas("c_SysSummary_"+channel_name+"_"+year_string, "c_"+channel_name+"_"+year_string, 50, 50, 950, 900);
    c1->cd();

    TPad* pad1 = new TPad("pad1","pad1",0,0.35,1,1);
    pad1->SetBottomMargin(0.05);
    pad1->SetTopMargin(0.1);
    pad1->SetTicks(1);
    pad1->Draw();
    pad1->cd();

    //c1->SetBottomMargin(0.2);
    //c1->SetTopMargin(0.08);
    //c1->SetTicks();
    //c1->SetLogx();
    //c1->SetGridy();
    //c1->SetGridx();

    double massBin[5] = {1., 2., 3., 4., 5.};
    int marker_style_sys = 24;
    int marker_color_sys = 900;

    std::map<TString, Double_t>::iterator it_for_sysName;
    std::map<TString, Double_t>::iterator end_for_sysName;

    std::vector<double*> systematicValues;
    std::vector<TH1*> pv_hists;

    std::vector<double> default_values;

    if(detector_unfold)
    {
        it_for_sysName = meanPtErr_sysdata_det_unf.at(0).begin();
        end_for_sysName = meanPtErr_sysdata_det_unf.at(0).end();
        default_values = meanPt_data_det_unf;
    }
    else
    {
        it_for_sysName = meanPtErr_sysdata_pre_fsr.at(0).begin();
        end_for_sysName = meanPtErr_sysdata_pre_fsr.at(0).end();
        default_values = meanPt_data_pre_fsr;
    }

    int nth_systematicValues = 0;

    TLegend* leg = new TLegend(0.45, 0.35, 0.95, 0.85,"","brNDC");
    leg->SetTextSize(0.04);
    leg->SetNColumns(2);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);

    TH1* hTotal_uncert = new TH1D("hTotalUncert","hTotalUncert", nMassBin, 0, nMassBin);
    for(int imass = 0; imass < nMassBin; imass++)
    {
        hTotal_uncert->Fill(imass, meanPtTotRelErr_data_pre_fsr[imass]);
    }
    hTotal_uncert->SetLineColor(kBlack);
    hTotal_uncert->SetMarkerColor(kBlack);
    hTotal_uncert->SetMarkerSize(0);
    hTotal_uncert->SetLineStyle(2);
    hTotal_uncert->SetLineWidth(2);
    hTotal_uncert->SetMinimum(0);
    hTotal_uncert->SetMaximum(3);
    hTotal_uncert->GetYaxis()->SetTitle("#splitline{Systematic uncertainty}{in average p_{T} (%)}");
    //hTotal_uncert->GetXaxis()->SetTitle("Mass bin");
    hTotal_uncert->GetXaxis()->SetTitleSize(30);
    hTotal_uncert->GetXaxis()->SetTitleOffset(1.2);
    hTotal_uncert->GetYaxis()->SetTitleSize(35);
    hTotal_uncert->GetYaxis()->SetTitleOffset(1.2);
    hTotal_uncert->Draw("p9histl");

    hTotal_uncert->GetXaxis()->SetBinLabel(hTotal_uncert->GetXaxis()->FindBin(0.), "");
    leg->AddEntry(hTotal_uncert, "Total uncertainty", "pl");
        
    double stat_uncert[5] = {0., 0., 0., 0., 0.};

    std::map<TString, int> marker_colors;
    std::map<TString, int> marker_styles;

    while(it_for_sysName != end_for_sysName)
    {
        if(it_for_sysName->first != "Stat")
            systematicValues.push_back(new double[5]);

        for(int i = 0; i < nMassBin; i++)
        {
            double default_value = 1.;

            if(relative)
            {
                default_value = default_value/ default_values.at(i) * 100.;
            }
            mapdata* p_map;
            if(detector_unfold)
            {
                p_map = &meanPtErr_sysdata_det_unf.at(i);
            }
            else
            {
                p_map = &meanPtErr_sysdata_pre_fsr.at(i);
            }
            if(it_for_sysName->first == "Stat")
            {
                stat_uncert[i] = (*p_map)[it_for_sysName->first] * default_value; 
            }
            else{
                systematicValues.at(nth_systematicValues)[i] = (*p_map)[it_for_sysName->first] * default_value;
            }
        }
        if(it_for_sysName->first == "Stat")
        {
            it_for_sysName++;
            //nth_systematicValues++;
            continue;
        }

        // draw
        TString draw_option = "p9hist same";

        double temp_marker_color = marker_color_sys;
        if(it_for_sysName->first == "PDFerror")
        {
            marker_color_sys = 870; //kAzure + 10
        }
        if(it_for_sysName->first == "Scale")
        {
            marker_color_sys = 601; //kBlue + 1
        }

        // marker style
        if(it_for_sysName->first == "AlphaS")
        {
            marker_style_sys = 24;
        }
        if(it_for_sysName->first == "Alt")
        {
            marker_style_sys = 25;
        }
        if(it_for_sysName->first == "IdSF")
        {
            marker_style_sys = 26;
        }
        if(it_for_sysName->first == "IsoSF")
        {
            marker_style_sys = 3;
        }
        if(it_for_sysName->first == "L1Prefire")
        {
            marker_style_sys = 27;
        }
        if(it_for_sysName->first == "PU")
        {
            marker_style_sys = 28;
        }
        if(it_for_sysName->first == "QED_FSR")
        {
            marker_style_sys = 29;
            marker_color_sys = kRed; 
        }
        if(it_for_sysName->first == "Scale")
        {
            marker_style_sys = 30;
        }
        if(it_for_sysName->first == "Stat")
        {
            marker_style_sys = 31;
        }
        if(it_for_sysName->first == "lepMom")
        {
            marker_style_sys = 32;
        }
        if(it_for_sysName->first == "recoSF")
        {
            marker_style_sys = 33;
        }
        if(it_for_sysName->first == "trgSF")
        {
            marker_style_sys = 34;
        }

        pv_hists.push_back(new TH1D(it_for_sysName->first, it_for_sysName->first, nMassBin, 0, nMassBin));
        for(int imass = 0; imass < nMassBin; imass++)
        {
            pv_hists.at(nth_systematicValues)->Fill(imass, systematicValues.at(nth_systematicValues)[imass]);
        }
        marker_colors[it_for_sysName->first] = marker_color_sys;
        marker_styles[it_for_sysName->first] = marker_style_sys;

        pv_hists.at(nth_systematicValues)->SetLineColor(marker_color_sys);
        pv_hists.at(nth_systematicValues)->SetMarkerColor(marker_color_sys);
        pv_hists.at(nth_systematicValues)->SetMarkerStyle(marker_style_sys);
        pv_hists.at(nth_systematicValues)->SetMarkerSize(1.2);
        pv_hists.at(nth_systematicValues)->SetLineStyle(1);
        pv_hists.at(nth_systematicValues)->Draw(draw_option);
        leg->AddEntry(pv_hists.at(nth_systematicValues), it_for_sysName->first=="Alt" ? "Unfold":it_for_sysName->first, "pl");

        it_for_sysName++;
        nth_systematicValues++;

        marker_color_sys = temp_marker_color;

        marker_style_sys++;
        marker_color_sys++;
    }

    TH1* hstat_uncert = new TH1D("stat", "stat", nMassBin, 0, nMassBin);
    for(int imass = 0; imass < nMassBin; imass++)
    {
        hstat_uncert->Fill(imass, stat_uncert[imass]);
    }
    hstat_uncert->Draw("p9histl same");
    hstat_uncert->SetLineColor(kGray+2);
    hstat_uncert->SetMarkerColor(kGray+2);
    hstat_uncert->SetMarkerSize(0);
    hstat_uncert->SetLineStyle(3);
    hstat_uncert->SetLineWidth(2);
    leg->AddEntry(hstat_uncert, "Stat. uncertainty", "pl");

    leg->Draw();

    CMS_lumi( pad1, 4, 0);
    c1->cd();

    TPad* pad2 = new TPad("pad2","pad2",0,0,1,0.35);
    pad2->SetBottomMargin(0.25);
    pad2->SetTopMargin(0.);
    pad2->SetTicks(1);
    //pad2->SetLogy();
    pad2->Draw();
    pad2->cd();

    if(detector_unfold)
    {
        it_for_sysName = meanPtErr_sysdata_det_unf.at(0).begin();
        end_for_sysName = meanPtErr_sysdata_det_unf.at(0).end();
    }
    else
    {
        it_for_sysName = meanPtErr_sysdata_pre_fsr.at(0).begin();
        end_for_sysName = meanPtErr_sysdata_pre_fsr.at(0).end();
    }

    std::vector<TH1*> pv_syshists;
    TH1* hunfolded_data = NULL;
    hunfolded_data  = nomPtFSRUnfold->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UO]",kFALSE);
    nth_systematicValues = 0;

    while(it_for_sysName != end_for_sysName)
    {


    //if(it_for_sysName->first == "Stat" || it_for_sysName->first == "PDEerror")
    //{ 
    //    it_for_sysName++;
    //    continue;
    //}

    // lets draw systematic in each bins...
    pv_syshists.push_back((TH1F*)hunfolded_data->Clone(it_for_sysName->first));

    bool firstIter = true;
    bool fullSys = false;
    for(int imass = 0; imass < nMassBin; imass++)
    {
        // initialize errors the first iteration
        if(firstIter)
        {
            for(int ibin = 1; ibin<pv_syshists[nth_systematicValues]->GetNbinsX()+1;ibin++)
            {

                pv_syshists[nth_systematicValues]->SetBinContent(ibin, 0.);
                pv_syshists[nth_systematicValues]->SetBinError(ibin, 0.);
            }
        }
        // update hunfolded_sys_err
        makeSystBand(var, imass, it_for_sysName->first, fullSys, hunfolded_data, pv_syshists[nth_systematicValues]);   
        firstIter = false;
    }

    TString option = "histe";
    if(nth_systematicValues != 0)
        option = "histsamee";

    pv_syshists[nth_systematicValues]->SetMaximum(18.);
    //pv_syshists[nth_systematicValues]->SetMinimum(1e-3);
    pv_syshists[nth_systematicValues]->SetLineColor(marker_colors[it_for_sysName->first]);
    pv_syshists[nth_systematicValues]->SetMarkerColor(marker_colors[it_for_sysName->first]);
    pv_syshists[nth_systematicValues]->SetMarkerStyle(marker_styles[it_for_sysName->first]);
    pv_syshists[nth_systematicValues]->SetMarkerSize(1.);
    pv_syshists[nth_systematicValues]->SetLineStyle(1);
    pv_syshists[nth_systematicValues]->SetLineWidth(1);
    pv_syshists[nth_systematicValues]->Draw(option);

    it_for_sysName++;
    nth_systematicValues++;
    }
    
    TLine binBoundary;
    binBoundary.SetLineColor(kGray+2); 
    binBoundary.SetLineStyle(2);
    for(int imass = 0; imass < nMassBin; imass++) 
    {
        binBoundary.DrawLine(pv_syshists[0]->GetXaxis()->GetBinUpEdge(nPtBin + nPtBin * imass), 0., pv_syshists[0]->GetXaxis()->GetBinUpEdge(nPtBin + nPtBin * imass), 18.);
    }
    pv_syshists[0]->GetXaxis()->SetTitle("p_{T} bin number");
    pv_syshists[0]->GetYaxis()->SetTitle("#splitline{Systematic uncertainty}{in #frac{N}{dp_{T}dM} (%)}");


    c1->cd();
    c1->SaveAs(outpdf + channel_name + "_" + year_string + ".pdf");

    int pv_size = systematicValues.size();
    for(int i = 0; i < pv_size; i++)
    {
        delete [] systematicValues[i];
        delete pv_hists[i];
    }
    int pv_sysSize = pv_syshists.size();
    for(int i = 0; i < pv_sysSize; i++)
    {
        delete pv_syshists[i];
    }

    delete hstat_uncert;
    delete pad1;
    delete c1;
}

void ISRUnfold::makeSystBand(const TString var, const int nthMassBin, const TString sysName, const bool fullSys, const TH1* hunfolded_data, TH1* hunfolded_sys_err)
{
    cout << "makeSystBanc called" << endl;

    bool isFoldedSys = false;
    const TUnfoldBinningV17* temp_binning_gen_pt = nomPtUnfold->GetOutputBinning("Gen_Pt");
    const TUnfoldBinningV17* temp_binning_rec_pt = nomPtUnfold->GetInputBinning("Rec_Pt");
    const TUnfoldBinningV17* temp_binning_gen_mass = nomMassUnfold->GetOutputBinning("Gen_Mass");
    const TUnfoldBinningV17* temp_binning_rec_mass = nomMassUnfold->GetInputBinning("Rec_Mass");

    // get mass bin definition from (pt, mass) bin definition
    const TVectorD* temp_tvecd = temp_binning_gen_pt->GetDistributionBinning(1);
    const Double_t* massBins = temp_tvecd->GetMatrixArray();

    // get number of pT bins
    const TVectorD* temp_tvecd_ = temp_binning_gen_pt->GetDistributionBinning(0);
    int nPtBin = temp_tvecd_->GetNrows() - 1;

    TString ibinMass;
    ibinMass.Form("%d", nthMassBin);

    std::map<TString, std::vector<TUnfoldDensityV17*>>::iterator it;
    std::map<TString, std::vector<TUnfoldDensityV17*>>::iterator it_end;

    if(!isFoldedSys)
    {
        it = sysPtFSRUnfold.begin();
        it_end = sysPtFSRUnfold.end();
    }
    else
    {
        it = sysPtUnfold.begin();
        it_end = sysPtUnfold.end();
    }

    while(it != it_end)
    {
        if(!fullSys && it->first != sysName)
        {
            it++;
            continue;
        }

        TH1 * hdatasys0_temp = NULL;
        TH1 * hdatasys1_temp = NULL;

        int systematic_variation_index = 0;
        int systematic_variation_index_mc = 0;

        std::vector<std::map<TString, Int_t>>* p_meanPtErrIdx_data;
        std::vector<std::map<TString, Int_t>>* p_meanPtErrIdx_mc;
        std::map<TString, std::vector<TUnfoldDensityV17*>> temp_sysTUnfDensity;

        // get distribution which gives the maximum variation on mean value
        if(var == "Pt")
        {
            if(!isFoldedSys)
            {
                p_meanPtErrIdx_data = &meanPtErrIdx_sysdata_pre_fsr;
                p_meanPtErrIdx_mc = &meanPtErrIdx_sysmc_pre_fsr;
                temp_sysTUnfDensity = sysPtFSRUnfold;
            }
            else
            {
                p_meanPtErrIdx_data = &meanPtErrIdx_sysdata_det_unf;
                p_meanPtErrIdx_mc = &meanPtErrIdx_sysmc_det_unf;
                temp_sysTUnfDensity = sysPtUnfold;
            }
            // note except QED_FSR, use the systematic histogram giving maximum variation on mean value
            // to calculate systematic error on each bin.
            // this is not correct for some systematic sources, for example PDF error.
            if(it->first != "QED_FSR")
            {
                systematic_variation_index = p_meanPtErrIdx_data->at(nthMassBin)[it->first];
                {
                    hdatasys0_temp = temp_sysTUnfDensity[it->first].at(systematic_variation_index)->GetOutput("hunfolded_pt_systemp",0,0,"pt[UO];mass[UO]",kFALSE);
                }

            }
            else
            {
                {
                    hdatasys0_temp = temp_sysTUnfDensity[it->first].at(0)->GetOutput("hunfolded_pt_systemp_fsr0",0,0,"pt[UO];mass[UO]",kFALSE);
                    hdatasys1_temp = temp_sysTUnfDensity[it->first].at(1)->GetOutput("hunfolded_pt_systemp_fsr1",0,0,"pt[UO];mass[UO]",kFALSE);
                }
            }
        }

        TH1F* ratio0_temp = NULL;
        TH1F* ratio1_temp = NULL;

        if(it->first != "QED_FSR")
        {
            {
                ratio0_temp = (TH1F*)hunfolded_data->Clone("ratio");
                ratio0_temp->Divide(hdatasys0_temp);
            }
        }
        else
        {
            {
                ratio0_temp = (TH1F*)hunfolded_data->Clone("ratio");
                ratio0_temp->Divide(hdatasys0_temp);

                ratio1_temp = (TH1F*)hunfolded_data->Clone("ratio");
                ratio1_temp->Divide(hdatasys1_temp);
            }
        }

        // loop over each bin of unfolded histogram and set error due to the systematic source
        // note errors of hunfolded_sys_err initialized as zero before looping over all systematic sources
        // for mass, update only the specific bins in the the current mass range

        int firstBin = 1;
        int lastBin = hunfolded_sys_err->GetNbinsX()+1;
        if(var == "Mass")
        {
            firstBin = hdatasys0_temp->GetXaxis()->FindBin(massBins[nthMassBin]+0.01);
            lastBin = hdatasys0_temp->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01)+1;
        }
        if(var == "Pt")
        {
            firstBin = firstBin + nPtBin * nthMassBin;
            lastBin = firstBin + nPtBin;
        }

        cout << "ibin: " << firstBin << " lastBin: " << lastBin << endl;
        for(int ibin = firstBin; ibin < lastBin; ibin++)
        {

            Double_t temp_ratio_err = 0.;
            Double_t previous_ratio_err = hunfolded_sys_err->GetBinContent(ibin);

            if(it->first != "QED_FSR")
            {
                temp_ratio_err = fabs(1. - ratio0_temp->GetBinContent(ibin)) * 100.;
            }
            else
            {
                temp_ratio_err = fabs(ratio0_temp->GetBinContent(ibin) - ratio1_temp->GetBinContent(ibin)) * 100;
            }

            cout << "uncertainty: " << sqrt(pow(previous_ratio_err, 2) + pow(temp_ratio_err, 2)) << endl;
            hunfolded_sys_err->SetBinContent(ibin, sqrt(pow(previous_ratio_err, 2) + pow(temp_ratio_err, 2)));
            hunfolded_sys_err->SetBinError(ibin, 1e-6);
            
        }// loop for bin contents

        delete ratio0_temp;
        delete ratio1_temp;
        delete hdatasys0_temp;
        delete hdatasys1_temp;

        it++;
    }// end of while

}

void ISRUnfold::drawSysPlots(TString outpdf, int nthMassBin, TString sysName, bool detector_unfold)
{

    gROOT->SetBatch();

    const TUnfoldBinningV17* temp_binning_gen_pt = nomPtUnfold->GetOutputBinning("Gen_Pt");
    const TUnfoldBinningV17* temp_binning_rec_pt = nomPtUnfold->GetInputBinning("Rec_Pt");
    const TUnfoldBinningV17* temp_binning_gen_mass = nomMassUnfold->GetOutputBinning("Gen_Mass");
    const TUnfoldBinningV17* temp_binning_rec_mass = nomMassUnfold->GetInputBinning("Rec_Mass");

    // get mass bin definition from (pt, mass) bin definition
    const TVectorD* temp_tvecd = temp_binning_gen_pt->GetDistributionBinning(1);
    const Double_t* massBins = temp_tvecd->GetMatrixArray();

    TString ibinMass;
    ibinMass.Form("%d", nthMassBin);

    c1 = new TCanvas("c1","c1", 50, 50, 800, 800);
    c1->cd();

    TPad *pad1 = new TPad("pad1","pad1",0,0.0,1,1);
    pad1->SetBottomMargin(0.12);
    pad1->SetTopMargin(0.1);
    pad1->SetTicks(1);
    pad1->Draw();
    pad1->cd();

    gStyle->SetOptFit(0);
    //c1->SetBottomMargin(0.12);
    //c1->SetTopMargin(0.1);
    //c1->SetTicks(1);

    int sysSize = 0;
    double amaxMass = 0.;
    double aminMass = 0.;
    double amaxPt = 0.;
    double aminPt = 0.;

    double amax = 0.;
    double amin = 0.;

    bool isFullSys = false;
    bool isRangeUpdated = false;

    int marker_color_sys = 632;
    int marker_style_sys = 24;

    if(sysName == "full")
    {
        isFullSys = true;
        marker_color_sys = 900;
    }

    TGraphErrors* p_grNominal;
    TGraphErrors* p_grStatErr;
    TGraph* grSys;
    vector<TGraph*> pv_grSys;

    if(detector_unfold)
    {
        // set axis min and max range
        amaxMass=meanMass_data_det_unf[nthMassBin] + (2. * meanMassTotErr_data_det_unf[nthMassBin]);
        aminMass=meanMass_data_det_unf[nthMassBin] - (2. * meanMassTotErr_data_det_unf[nthMassBin]);
        amaxPt=meanPt_data_det_unf[nthMassBin] + (2. * meanPtTotErr_data_det_unf[nthMassBin]);
        aminPt=meanPt_data_det_unf[nthMassBin] - (2. * meanPtTotErr_data_det_unf[nthMassBin]);

        // nominal point with total uncertainty
        TGaxis::SetMaxDigits(4);
        p_grNominal = new TGraphErrors(1, &meanMass_data_det_unf[nthMassBin], &meanPt_data_det_unf[nthMassBin], &meanMassTotErr_data_det_unf[nthMassBin], &meanPtTotErr_data_det_unf[nthMassBin]);
        p_grNominal->SetLineColor(1);
        p_grNominal->SetFillColor(12);
        p_grNominal->SetFillStyle(3003);
        p_grNominal->SetMarkerStyle(20);
        p_grNominal->SetMarkerSize(1.2);
        p_grNominal->Draw("a2");
        p_grNominal->GetYaxis()->SetRangeUser(aminPt*0.995,amaxPt*1.005);
        p_grNominal->GetXaxis()->SetLimits(meanMass_data_det_unf[nthMassBin] - (amaxPt*1.005-aminPt*0.995)/2., meanMass_data_det_unf[nthMassBin] + (amaxPt*1.005-aminPt*0.995)/2.);
        //p_grNominal->GetXaxis()->SetLimits(aminMass*0.995,amaxMass*1.005);
        //p_grNominal->GetXaxis()->SetNdivisions(505, false);
        p_grNominal->GetYaxis()->SetTitle("Average p_{T} (GeV)");
        p_grNominal->GetXaxis()->SetTitle("Average Mass (GeV)");
        p_grNominal->GetXaxis()->SetTitleOffset(1.5);

        // meanPtStatErr_data
        p_grStatErr = new TGraphErrors(1, &meanMass_data_det_unf[nthMassBin], &meanPt_data_det_unf[nthMassBin], &meanMassStatErr_data_det_unf[nthMassBin], &meanPtStatErr_data_det_unf[nthMassBin]);

        TLegend* leg = new TLegend(0.6, 0.4, 0.85, 0.9,"","brNDC");
        leg->SetTextSize(0.02);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);

        leg->AddEntry(p_grNominal, "Total uncertainty", "f");
        leg->AddEntry(p_grStatErr,  "Nominal with stat. uncertainty" , "pe");

        std::map<TString, std::vector<TUnfoldDensityV17*>>::iterator it;
        if(!isFullSys)
        {
            it = sysPtUnfold.find(sysName);
        }
        else
        {
            it = sysPtUnfold.begin();
        }

        int p_index = 0;
        while(it != sysPtUnfold.end())
        {
            sysSize = it->second.size();

            // check if ranges needed to be updated
            for(int i=0;i<sysSize;i++)
            {
                double temp_amaxMass = amaxMass;
                double temp_aminMass = aminMass;
                double temp_amaxPt = amaxPt;
                double temp_aminPt = aminPt;

                amaxMass=TMath::Max(amaxMass, (meanMass_sysdata_det_unf.at(nthMassBin)[it->first])[i]);
                aminMass=TMath::Min(aminMass, (meanMass_sysdata_det_unf.at(nthMassBin)[it->first])[i]);

                amaxPt=TMath::Max(amaxPt, (meanPt_sysdata_det_unf.at(nthMassBin)[it->first])[i]);
                aminPt=TMath::Min(aminPt, (meanPt_sysdata_det_unf.at(nthMassBin)[it->first])[i]);

                if( fabs(amaxMass-temp_amaxMass) > 1e-5 || fabs(aminMass-temp_aminMass) > 1e-5
                    || fabs(amaxPt-temp_amaxPt) > 1e-5 || fabs(aminPt-temp_aminPt) > 1e-5)
                {
                    isRangeUpdated = true;
                }
            }

            // TODO make option to draw all systematic points here
            pv_grSys.push_back(new TGraph(sysSize, &(meanMass_sysdata_det_unf.at(nthMassBin)[it->first])[0], &(meanPt_sysdata_det_unf.at(nthMassBin)[it->first])[0]));
            pv_grSys[p_index]->SetLineColor(marker_color_sys);
            pv_grSys[p_index]->SetMarkerColor(marker_color_sys);
            pv_grSys[p_index]->SetMarkerStyle(marker_style_sys++);
            pv_grSys[p_index]->SetMarkerSize(1.2);
            pv_grSys[p_index]->SetLineStyle(1);
            pv_grSys[p_index]->Draw("pe same ");
            drawtext(pv_grSys[p_index]);

            leg->AddEntry(pv_grSys[p_index],      "Systematic source : " + it->first, "p");

            if(isFullSys)
            {
                it++;
                p_index++;
                marker_color_sys++;
            }
            else
            {
                break;
            }
        }
        p_grStatErr->Draw("pesame");
        leg->Draw();
    }
    else
    {

        // set axis min and max range
        amaxMass=meanMass_data_pre_fsr[nthMassBin] + (2. * meanMassTotErr_data_pre_fsr[nthMassBin]);
        aminMass=meanMass_data_pre_fsr[nthMassBin] - (2. * meanMassTotErr_data_pre_fsr[nthMassBin]);
        amaxPt=meanPt_data_pre_fsr[nthMassBin] + (2. * meanPtTotErr_data_pre_fsr[nthMassBin]);
        aminPt=meanPt_data_pre_fsr[nthMassBin] - (2. * meanPtTotErr_data_pre_fsr[nthMassBin]);

        // nominal point with total uncertainty
        TGaxis::SetMaxDigits(4);
        p_grNominal = new TGraphErrors(1, &meanMass_data_pre_fsr[nthMassBin], &meanPt_data_pre_fsr[nthMassBin], &meanMassTotErr_data_pre_fsr[nthMassBin], &meanPtTotErr_data_pre_fsr[nthMassBin]);
        p_grNominal->SetLineColor(1);
        p_grNominal->SetFillColor(12);
        p_grNominal->SetFillStyle(3003);
        p_grNominal->SetMarkerStyle(20);
        p_grNominal->SetMarkerSize(1.2);
        p_grNominal->Draw("a2");
        p_grNominal->GetYaxis()->SetRangeUser(aminPt*0.995,amaxPt*1.005);
        p_grNominal->GetXaxis()->SetLimits(meanMass_data_pre_fsr[nthMassBin] - (amaxPt*1.005-aminPt*0.995)/2., meanMass_data_pre_fsr[nthMassBin] + (amaxPt*1.005-aminPt*0.995)/2.);
        //p_grNominal->GetXaxis()->SetLimits(aminMass*0.995,amaxMass*1.005);
        //p_grNominal->GetXaxis()->SetNdivisions(505, false);
        p_grNominal->GetYaxis()->SetTitle("Average p_{T} (GeV)");
        p_grNominal->GetXaxis()->SetTitle("Average Mass (GeV)");
        p_grNominal->GetXaxis()->SetTitleOffset(1.5);

        // meanPtStatErr_data_pre_fsr
        p_grStatErr = new TGraphErrors(1, &meanMass_data_pre_fsr[nthMassBin], &meanPt_data_pre_fsr[nthMassBin], &meanMassStatErr_data_pre_fsr[nthMassBin], &meanPtStatErr_data_pre_fsr[nthMassBin]);

        TLegend* leg = new TLegend(0.6, 0.4, 0.85, 0.9,"","brNDC");
        leg->SetTextSize(0.02);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);

        leg->AddEntry(p_grNominal, "Total uncertainty", "f");
        leg->AddEntry(p_grStatErr,  "Nominal with stat. uncertainty" , "pe");

        std::map<TString, std::vector<TUnfoldDensityV17*>>::iterator it;
        if(!isFullSys)
        {
            it = sysPtFSRUnfold.find(sysName);
        }
        else
        {
            it = sysPtFSRUnfold.begin();
        }

        int p_index = 0;
        while(it != sysPtFSRUnfold.end())
        {
            sysSize = it->second.size();

            // check if ranges needed to be updated
            for(int i=0;i<sysSize;i++)
            {
                double temp_amaxMass = amaxMass;
                double temp_aminMass = aminMass;
                double temp_amaxPt = amaxPt;
                double temp_aminPt = aminPt;

                amaxMass=TMath::Max(amaxMass, (meanMass_sysdata_pre_fsr.at(nthMassBin)[it->first])[i]);
                aminMass=TMath::Min(aminMass, (meanMass_sysdata_pre_fsr.at(nthMassBin)[it->first])[i]);

                amaxPt=TMath::Max(amaxPt, (meanPt_sysdata_pre_fsr.at(nthMassBin)[it->first])[i]);
                aminPt=TMath::Min(aminPt, (meanPt_sysdata_pre_fsr.at(nthMassBin)[it->first])[i]);

                if( fabs(amaxMass-temp_amaxMass) > 1e-5 || fabs(aminMass-temp_aminMass) > 1e-5
                    || fabs(amaxPt-temp_amaxPt) > 1e-5 || fabs(aminPt-temp_aminPt) > 1e-5)
                {
                    isRangeUpdated = true;
                }
            }

            double temp_marker_color = marker_color_sys;
            if(it->first == "PDFerror")
            {
                marker_color_sys = 870; //kAzure + 10
            }
            if(it->first == "Scale")
            {
                marker_color_sys = 601; //kBlue + 1
            }

            pv_grSys.push_back(new TGraph(sysSize, &(meanMass_sysdata_pre_fsr.at(nthMassBin)[it->first])[0], &(meanPt_sysdata_pre_fsr.at(nthMassBin)[it->first])[0]));
            pv_grSys[p_index]->SetLineColor(marker_color_sys);
            pv_grSys[p_index]->SetMarkerColor(marker_color_sys);
            pv_grSys[p_index]->SetMarkerStyle(marker_style_sys++);
            pv_grSys[p_index]->SetMarkerSize(1.2);
            pv_grSys[p_index]->SetLineStyle(1);
            pv_grSys[p_index]->Draw("pe same ");
            drawtext(pv_grSys[p_index]);

            marker_color_sys = temp_marker_color;

            leg->AddEntry(pv_grSys[p_index],      "Systematic source : " + it->first, "p");

            if(isFullSys)
            {
                it++;
                p_index++;
                marker_color_sys++;
            }
            else
            {
                break;
            }
        }
        p_grStatErr->Draw("pesame");
        leg->Draw();

    }

    TString mass_cut_info;
    TString pt_cut_info;
    TString lepton_type;
    if(channel_name=="electron") lepton_type = "ee";
    else lepton_type = "#mu#mu";

    TString low_bound_, upper_bound_;
    low_bound_.Form("%d", (int)massBins[nthMassBin]);
    upper_bound_.Form("%d", (int)massBins[nthMassBin+1]);
    mass_cut_info = low_bound_ + " < M(" + lepton_type + ") < " + upper_bound_ + " (GeV)";
    pt_cut_info = "p_{T}(" + lepton_type + ") < 100 (GeV) ";

    TLatex mass_cut_info_;
    mass_cut_info_.SetTextFont(63);
    mass_cut_info_.SetTextSize(23);
    mass_cut_info_.DrawLatexNDC(0.2, 0.8, mass_cut_info);

    TLatex pt_cut_info_;
    pt_cut_info_.SetTextFont(63);
    pt_cut_info_.SetTextSize(23);
    pt_cut_info_.DrawLatexNDC(0.2, 0.75, pt_cut_info);

    CMS_lumi( pad1, 4, 0 );
    c1->cd();
    if (detector_unfold) c1->SaveAs(outpdf+"_"+ibinMass+"_"+sysName+".pdf");
    else c1->SaveAs(outpdf+"_"+ibinMass+"_pre_fsr_"+sysName+".pdf");

    delete p_grNominal, p_grStatErr;
    delete c1;

    int pv_size = pv_grSys.size();
    for(int i = 0; i < pv_size; i++)
    {
        delete pv_grSys[i];
    }
    //if(!detector_unfold)
    //    delete grSys;
}

void ISRUnfold::drawNominalRecoPlots(TString outpdf, TString filepath, TString var, int nthMassBin, TString sysName)
{

    TString ibinMass;
    ibinMass.Form("%d", nthMassBin);

    const TUnfoldBinningV17* temp_binning_gen_pt = nomPtUnfold->GetOutputBinning("Gen_Pt");
    const TUnfoldBinningV17* temp_binning_rec_pt = nomPtUnfold->GetInputBinning("Rec_Pt");
    const TUnfoldBinningV17* temp_binning_gen_mass = nomMassUnfold->GetOutputBinning("Gen_Mass");
    const TUnfoldBinningV17* temp_binning_rec_mass = nomMassUnfold->GetInputBinning("Rec_Mass");

    // get mass bin definition from (pt, mass) bin definition
    const TVectorD* temp_tvecd = temp_binning_gen_pt->GetDistributionBinning(1);
    const Double_t* massBins = temp_tvecd->GetMatrixArray();
    int nMassBin = temp_tvecd->GetNrows() - 1;

    gROOT->SetBatch();

    setTDRStyle();
    writeExtraText = true;
    extraText  = "work in progress";

    Int_t linecolorZ   = kOrange-3;
    Int_t fillcolorZ   = kOrange-2;
    Int_t linecolorEWK = kOrange+10;
    Int_t fillcolorEWK = kOrange+7;
    Int_t linecolorTop = kGreen+2;
    Int_t fillcolorTop = kGreen-5;

    // get binning definition
    const TUnfoldBinning* rec_binning = nomPtUnfold->GetInputBinning("Rec_Pt");

    TString PostFixForHist = "EEnominal";
    if(channel_name == "muon")
    {
        PostFixForHist = "MuMunominal";
    }

    TFile* filein = new TFile(filepath);
    TH1* hdy_m50 = (TH1*)filein->Get("detector_level/hist_ptll/histo_DYJetsTo"+PostFixForHist); // get DY
    TH1* hdy_m10to50 = (TH1*)filein->Get("detector_level/hist_ptll/histo_DYJets10to50To"+PostFixForHist); // get DY
    TH1* hZZ = (TH1*)filein->Get("detector_level/hist_ptll/histo_ZZ_pythianominal"); // 
    TH1* hWZ = (TH1*)filein->Get("detector_level/hist_ptll/histo_WZ_pythianominal"); // 
    TH1* hWW = (TH1*)filein->Get("detector_level/hist_ptll/histo_WW_pythianominal"); // 
    TH1* hdy_m50_tau = (TH1*)filein->Get("detector_level/hist_ptll/histo_DYJetsToTauTaunominal"); // 
    TH1* hdy_m10to50_tau = (TH1*)filein->Get("detector_level/hist_ptll/histo_DYJets10to50ToTauTaunominal"); // 
    TH1* hTT = (TH1*)filein->Get("detector_level/hist_ptll/histo_TTLL_powhegnominal"); // 
    //
    TH1* hdy_m50_NoUO = rec_binning->ExtractHistogram("hdy_m50", hdy_m50, 0, kTRUE, "pt[UO];mass[UOC"+ibinMass+"]");
    TH1* hdy_m10to50_NoUO = rec_binning->ExtractHistogram("hdy_m10to50", hdy_m10to50, 0, kTRUE, "pt[UO];mass[UOC"+ibinMass+"]");
    hdy_m50_NoUO->Add(hdy_m10to50_NoUO);
    doNorm(hdy_m50_NoUO, false);
    TH1* hZZ_NoUO = rec_binning->ExtractHistogram("hZZ", hZZ, 0, kTRUE, "pt[UO];mass[UOC"+ibinMass+"]");
    TH1* hWZ_NoUO = rec_binning->ExtractHistogram("hWZ", hWZ, 0, kTRUE, "pt[UO];mass[UOC"+ibinMass+"]");
    TH1* hWW_NoUO = rec_binning->ExtractHistogram("hWW", hWW, 0, kTRUE, "pt[UO];mass[UOC"+ibinMass+"]");
    TH1* hdy_m50_tau_NoUO = rec_binning->ExtractHistogram("hdy_m50_tau", hdy_m50_tau, 0, kTRUE, "pt[UO];mass[UOC"+ibinMass+"]");
    TH1* hdy_m10to50_tau_NoUO = rec_binning->ExtractHistogram("hdy_m10to50_tau", hdy_m10to50_tau, 0, kTRUE, "pt[UO];mass[UOC"+ibinMass+"]");
    TH1* hTT_NoUO = rec_binning->ExtractHistogram("hTT", hTT, 0, kTRUE, "pt[UO];mass[UOC"+ibinMass+"]");
    hZZ_NoUO->Add(hWZ_NoUO);
    hZZ_NoUO->Add(hWW_NoUO);
    hZZ_NoUO->Add(hdy_m50_tau_NoUO);
    hZZ_NoUO->Add(hdy_m10to50_tau_NoUO);
    doNorm(hZZ_NoUO, false);
    doNorm(hTT_NoUO, false);

    TH1* hdata;
    TH1F *ratio;

    // bkg subracted data
    PostFixForHist = "EGnominal";
    if(channel_name == "muon")
    {
        PostFixForHist = "Muonnominal";
    }
    hdata = (TH1*)filein->Get("detector_level/hist_ptll/histo_Double"+PostFixForHist); 
    TH1* hdata_NoUO = rec_binning->ExtractHistogram("hdata_NoUO", hdata, 0, kTRUE, "pt[UO];mass[UOC"+ibinMass+"]"); 
    doNorm(hdata_NoUO, false);

    c1=new TCanvas("c1", "c1", 50, 50, 800, 800);
    c1->cd();
    gStyle->SetOptStat(0);

    TPad *pad1 = new TPad("pad1","pad1",0,0.35,1,1);
    pad1->SetBottomMargin(0.01);
    pad1->SetTopMargin(0.1);
    pad1->SetTicks(1);
    pad1->SetLogy();
    pad1->Draw();
    pad1->cd();

    hdata_NoUO->SetTitle("");
    hdata_NoUO->Draw("p9e");
    hdata_NoUO->SetMarkerStyle(20);
    hdata_NoUO->SetMarkerSize(1.2);
    hdata_NoUO->SetLineColor(kBlack);
    hdata_NoUO->SetMinimum(2e-1);
    hdata_NoUO->SetMaximum(1e10);  
    hdata_NoUO->GetYaxis()->SetTitle("Events/bin");

    THStack *hsMCs = new THStack("hsMCs","hsMCs");
    hdy_m50_NoUO->SetLineColor(linecolorZ);
    hdy_m50_NoUO->SetFillColor(fillcolorZ);
    hZZ_NoUO->SetLineColor(linecolorEWK);
    hZZ_NoUO->SetFillColor(fillcolorEWK);
    hTT_NoUO->SetLineColor(linecolorTop);
    hTT_NoUO->SetFillColor(fillcolorTop);
    hsMCs->Add(hTT_NoUO);
    hsMCs->Add(hZZ_NoUO);
    hsMCs->Add(hdy_m50_NoUO);
    hsMCs->Draw("histsame");
    hdata_NoUO->Draw("p9samee");
    pad1->RedrawAxis();


    TString lepton_type;
    TString cut_info;
    if(channel_name=="electron") lepton_type = "ee";
    else lepton_type = "#mu#mu";

    if(var=="Pt")
    {
        TString low_bound_, upper_bound_;
        low_bound_.Form("%d", (int)massBins[nthMassBin]);
        upper_bound_.Form("%d", (int)massBins[nthMassBin+1]);
        cut_info = low_bound_ + " < M(" + lepton_type + ") < " + upper_bound_ + " (GeV)";
    }
    else
    {
        cut_info = "p_{T}(" + lepton_type + ") < 100 (GeV) ";
    }

    TLatex cut_info_;
    cut_info_.SetTextFont(63);
    cut_info_.SetTextSize(20);
    cut_info_.DrawLatexNDC(0.2, 0.8, cut_info);

    TLegend* leg_nom = new TLegend(0.7, 0.45, 0.95, 0.85,"","brNDC");
    //leg_nom->SetNColumns(2);
    leg_nom->SetTextSize(0.055);
    leg_nom->SetFillStyle(0);
    leg_nom->SetBorderSize(0);

    leg_nom->AddEntry(hdata_NoUO, "Data", "pl");
    leg_nom->AddEntry(hdy_m50_NoUO, "Drell-Yan", "F");
    leg_nom->AddEntry(hZZ_NoUO, "EWK", "F");
    leg_nom->AddEntry(hTT_NoUO, "TT", "F");
    leg_nom->Draw();

    TLatex chi2_norm;
    TString chi2_;

    // TODO add Mass option
    //if(var == "Pt" )
    //{
    //    chi2_.Form("%f", Chi2Test(hdata_NoUO, hdysigNoUO));
    //    chi2_norm.DrawLatexNDC(0.2, 0.3, "#chi^{2}/NDOF= " + chi2_);
    //}

    c1->cd();

    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.35);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.25);
    pad2->SetTicks(1);
    pad2->SetGridy(1);
    pad2->Draw();

    pad2->cd();

    ratio = ((TH1F*)hdy_m50_NoUO->Clone("ratio"));
    ratio->Add(hZZ_NoUO);
    ratio->Add(hTT_NoUO);

    ratio->Divide(hdata_NoUO);
    ratio->Draw("pe");
    ratio->GetYaxis()->SetTitle("MC/ Data");
    ratio->GetXaxis()->SetTitle("p_{T} (GeV)");
    ratio->SetMinimum(0.6);
    ratio->SetMaximum(1.4);
    ratio->GetYaxis()->SetLabelFont(43);
    ratio->GetYaxis()->SetLabelSize(25);
    ratio->GetYaxis()->SetNdivisions(504);
    ratio->GetYaxis()->SetRangeUser(0.61, 1.39);
    ratio->GetYaxis()->SetTitleFont(43);
    ratio->GetYaxis()->SetTitleSize(25);
    ratio->GetYaxis()->SetTitleOffset(2.0);
    ratio->SetTitle("");
    ratio->SetLineColor(kBlack);
    //ratio->SetLineStyle(2);

    CMS_lumi( c1, 4, 0 );
    c1->cd();
    c1->SaveAs(outpdf+"_input_"+ibinMass+".pdf");

    delete hdata_NoUO;
    delete hdy_m50_NoUO;
    delete ratio;
    delete pad1;
    delete pad2;
    delete c1;

    filein->Close();
}

void ISRUnfold::drawClosurePlots(int detOrFSR_unfold, TString filepath, TString outpdf, TString var, int nthMassBin)
{
    gROOT->SetBatch();

    const TUnfoldBinningV17* temp_binning_gen_pt = nomPtUnfold->GetOutputBinning("Gen_Pt");
    // get mass bin definition from (pt, mass) bin definition
    const TVectorD* temp_tvecd = temp_binning_gen_pt->GetDistributionBinning(1);
    const Double_t* massBins = temp_tvecd->GetMatrixArray();
    int nMassBin = temp_tvecd->GetNrows() - 1;
    int nPtBin = temp_binning_gen_pt->GetDistributionBinning(0)->GetNrows() - 1;

    int firstMassBin = nthMassBin;
    int lastBin = nthMassBin + 1;
    bool allMassBins = false;
    if(nMassBin <= nthMassBin)
    {
        firstMassBin = 0;
        lastBin = nMassBin; 
        allMassBins = true;
    }

    TString low_mass_edge[5] = {"50.", "65.", "80.", "100.", "200."};
    TString high_mass_edge[5] = {"65.", "80.", "100.", "200.", "350."};

    setTDRStyle();
    writeExtraText = true;
    extraText  = "Simulation";
    TString ibinMass;
    ibinMass.Form("C%d", nthMassBin);
    if(allMassBins)
        ibinMass = "";

    TH1* hunfolded_data = NULL;
    TH1* hpreFSR_mc = NULL;
    TH1F *ratio = NULL;

    TUnfoldDensityV17* temp_tunfoldDensity = NULL;

    if(channel_name == "muon")
    {
        low_mass_edge[0] = "40.";
        low_mass_edge[1] = "60.";

        high_mass_edge[0] = "60.";
        high_mass_edge[1] = "80.";

    }

    // get nominal unfoled result
    if(var == "Pt" )
    {
        TFile* filein = new TFile(filepath);
        TString dir_path;

        bool useAxis = kTRUE;
        if(allMassBins)
            useAxis = kFALSE;

        if(detOrFSR_unfold == 0)
        {
            temp_tunfoldDensity = nomPtUnfold_closure;
        }
        else if(detOrFSR_unfold == 1)
        {
            temp_tunfoldDensity = nomPtFSRUnfold_closure;
        }
        else
        {
            exit(EXIT_FAILURE);
        }
            
        hunfolded_data  = temp_tunfoldDensity->GetOutput("hunfolded_pt_closure",0,0,"pt[UO];mass[UO"+ibinMass+"]", useAxis);
        hpreFSR_mc = (TH1*)hunfolded_data->Clone("copy_histogram");

        int ibin_global_mc = 1;
        for(int imass = firstMassBin; imass < lastBin; imass++)
        {
            if(detOrFSR_unfold == 0)
            {
                dir_path = channel_name + "_fiducial_phase_post_FSR_dRp1_split_p5_m"+ low_mass_edge[imass] + "to" + high_mass_edge[imass] + "/ptll/";
            }
            else if(detOrFSR_unfold == 1)
            {
                dir_path = channel_name + "_full_phase_split_p5_m"+ low_mass_edge[imass] + "to" + high_mass_edge[imass] + "/ptll_preFSR/";
            }
            else
            {
                exit(EXIT_FAILURE);
            }

            TH1* temp_hist = NULL;
            temp_hist = (TH1*)filein->Get(dir_path + "histo_DYJets");
            temp_hist->SetDirectory(0);
            temp_hist->Add((TH1*)filein->Get(dir_path + "histo_DYJets10to50"));

            for(int ibin = 1; ibin<temp_hist->GetNbinsX()+1;ibin++)
            {
                hpreFSR_mc->SetBinContent(ibin_global_mc, temp_hist->GetBinContent(ibin));
                hpreFSR_mc->SetBinError(ibin_global_mc, temp_hist->GetBinError(ibin));

                ibin_global_mc++;
            }
            delete temp_hist;
        }

        ratio= ((TH1F*)hunfolded_data->Clone("ratio"));
        ratio->Divide(hpreFSR_mc);
    }

    if(var == "Mass" )
    {

        TFile* filein = new TFile(filepath);
        TString dir_path;

        if(detOrFSR_unfold == 0)
        {
            temp_tunfoldDensity = nomMassUnfold_closure;
            dir_path = channel_name + "_fiducial_phase_post_FSR_dRp1_split_p5_m"+ low_mass_edge[nthMassBin] + "to" + high_mass_edge[nthMassBin] + "/mll_" + channel_name + "/";
        }
        else if(detOrFSR_unfold == 1)
        {
            temp_tunfoldDensity = nomMassFSRUnfold_closure;
            dir_path = channel_name + "_full_phase_split_p5_m"+ low_mass_edge[nthMassBin] + "to" + high_mass_edge[nthMassBin] + "/mll_preFSR_" + channel_name + "/";
        }
        else
        {
            exit(EXIT_FAILURE);
        }

        hunfolded_data  = temp_tunfoldDensity->GetOutput("hunfolded_mass_closure",0,0,"mass[UO];pt[UOC0]",kTRUE);
        hunfolded_data->GetXaxis()->SetRange(hunfolded_data->GetXaxis()->FindBin(massBins[nthMassBin]+0.01),hunfolded_data->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));

        hpreFSR_mc   = (TH1*)filein->Get(dir_path + "histo_DYJets");
        hpreFSR_mc->Add((TH1*)filein->Get(dir_path + "histo_DYJets10to50"));
        hpreFSR_mc->GetXaxis()->SetRange(hpreFSR_mc->GetXaxis()->FindBin(massBins[nthMassBin]+0.01), hpreFSR_mc->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));

        ratio= ((TH1F*)hunfolded_data->Clone("ratio"));
        ratio->Divide(hpreFSR_mc);
        ratio->GetXaxis()->SetRange(hpreFSR_mc->GetXaxis()->FindBin(massBins[nthMassBin]+0.01), hpreFSR_mc->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));
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
    hpreFSR_mc->Draw("histsamee");
    hpreFSR_mc->SetMarkerSize(0);

    hunfolded_data->SetMarkerStyle(20);
    hunfolded_data->SetMarkerSize(1.5);
    hunfolded_data->SetLineColor(kBlack);
    hunfolded_data->GetYaxis()->SetTitle("Events/bin");
    hunfolded_data->SetMinimum(10.);
    hunfolded_data->SetMaximum(9e9);

    if(allMassBins)
    {
        TLine binBoundary;
        binBoundary.SetLineColor(kGray+2); 
        binBoundary.SetLineStyle(2);

        for(int imass = firstMassBin; imass < lastBin; imass++) 
        {
            binBoundary.DrawLine(hunfolded_data->GetXaxis()->GetBinUpEdge(nPtBin + nPtBin * imass), 10., hunfolded_data->GetXaxis()->GetBinUpEdge(nPtBin + nPtBin * imass), 9e9);
        }
    }

    hpreFSR_mc->SetLineColor(kRed);

    TString mean_nom;
    TString meanMC_nom;
    mean_nom.Form("%.2f", hunfolded_data->GetMean());
    meanMC_nom.Form("%.2f", hpreFSR_mc->GetMean());

    TLegend* leg_nom = new TLegend(0.45, 0.7, 0.75, 0.9,"","brNDC");
    //leg_nom->SetNColumns(2);
    leg_nom->SetTextSize(0.05);
    leg_nom->SetFillStyle(0);
    leg_nom->SetBorderSize(0);
    if(!allMassBins)
    {
    leg_nom->AddEntry(hunfolded_data, "Unfolded MC (" + mean_nom + ")", "pl");
    leg_nom->AddEntry(hpreFSR_mc, "Gen MC (" + meanMC_nom + ")", "l");
    }
    else
    {
    leg_nom->AddEntry(hunfolded_data, "Unfolded MC", "pl");
    leg_nom->AddEntry(hpreFSR_mc, "MC Gen", "l");
    }
    leg_nom->Draw();

    c1->cd();

    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.4);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.3);
    pad2->SetTicks(1);
    //pad2->SetGridy(1);
    pad2->Draw();

    pad2->cd();
    ratio->Draw("histe");
    ratio->GetYaxis()->SetTitle("#splitline{Unfolded MC/}{MC Gen}");
    ratio->GetXaxis()->SetTitle("p_{T} at pre FSR(GeV)");
    if(allMassBins)
        ratio->GetXaxis()->SetTitle("p_{T} bin number");
    ratio->SetMinimum(0.6);
    ratio->SetMaximum(1.4);
    ratio->SetTitle("");
    ratio->GetXaxis()->SetTitleOffset(5.5);
    ratio->GetYaxis()->SetNdivisions(504);
    ratio->GetYaxis()->SetRangeUser(0.61, 1.39);
    ratio->SetMarkerSize(1.5);

    if(allMassBins)
    {
        TLine binBoundary;
        binBoundary.SetLineColor(kGray+2); 
        binBoundary.SetLineStyle(2);

        for(int imass = firstMassBin; imass < lastBin; imass++) 
        {
            binBoundary.DrawLine(hunfolded_data->GetXaxis()->GetBinUpEdge(nPtBin + nPtBin * imass), 0.61, hunfolded_data->GetXaxis()->GetBinUpEdge(nPtBin + nPtBin * imass), 1.39);
        }
    }

    CMS_lumi( c1, 4, 0 );
    c1->cd();
    c1->SaveAs(outpdf+"_"+ibinMass+"_"+var+".pdf");

    delete hunfolded_data;
    delete hpreFSR_mc;
    delete pad1;
    delete pad2;
    delete c1;
}

void ISRUnfold::drawSysComparionPlots(TString outpdf, TString var, int nthMassBin, TString sysName, bool isDetector)
{

    gROOT->SetBatch();

    const TUnfoldBinningV17* temp_binning_gen_pt = nomPtUnfold->GetOutputBinning("Gen_Pt");
    const TUnfoldBinningV17* temp_binning_rec_pt = nomPtUnfold->GetInputBinning("Rec_Pt");
    const TUnfoldBinningV17* temp_binning_rec_mass = nomMassUnfold->GetInputBinning("Rec_Mass");
    const TUnfoldBinningV17* temp_binning_gen_mass = nomMassUnfold->GetOutputBinning("Gen_Mass");
    // get mass bin definition from (pt, mass) bin definition
    const TVectorD* temp_tvecd = temp_binning_gen_pt->GetDistributionBinning(1);
    const Double_t* massBins = temp_tvecd->GetMatrixArray();


    setTDRStyle();
    writeExtraText = true;
    extraText  = "work in progress";

    TString ibinMass;
    ibinMass.Form("%d", nthMassBin);

    TH1*  h_nominal_unfolded_data = NULL;
    TH1F* ratio = NULL;
    TF1* landau_nom = new TF1("landau_nom","[0]*TMath::Landau(x,[1],[2]) + [3] * expo(x * [4])", 0., 100);
    TF1* landau_nom_ = new TF1("landau_nom_","[0]*TMath::Landau(x,[1],[2]) + [3] * expo(x * [4])", 0., 100);

    if(sysName!="QED_FSR")
    {
        if(isDetector)
        {
            if(var=="Pt")
            {
                h_nominal_unfolded_data  = nomPtUnfold->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
            }
            else
            {
                h_nominal_unfolded_data  = nomMassUnfold->GetOutput("hunfolded_mass_temp",0,0,"mass[UO];pt[UOC0]",kTRUE);
            }
        }
        else
        {
            if(var=="Pt")
            {
                h_nominal_unfolded_data  = nomPtFSRUnfold->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
                doNorm(h_nominal_unfolded_data);
            }
            else
            {
                h_nominal_unfolded_data  = nomMassFSRUnfold->GetOutput("hunfolded_mass_temp",0,0,"mass[UO];pt[UOC0]",kTRUE);
            }
        }
    }
    else
    {

        if(var=="Pt")
        {
            h_nominal_unfolded_data  = sysPtFSRUnfold[sysName].at(1)->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
            doNorm(h_nominal_unfolded_data);

            //TF1* landau_nom = new TF1("landau_nom","[0]*TMath::Landau(x,[1],[2]) + [3] * expo(x * [4])", 0., 100);
            //landau_nom->SetParameters(2., 6., 2., 0.05, -0.01);
            //landau_nom->SetParLimits(0, 0.0, 5.);
            //landau_nom->SetParLimits(1, 2.0, 9.);
            //landau_nom->SetParLimits(2, 1., 4.0);
            //landau_nom->SetParLimits(3, -0.05, 0.12);
            //landau_nom->SetParLimits(4, -0.2, -0.001);
            //h_nominal_unfolded_data->Fit(landau_nom);
            //TH1* h_ = landau_nom->GetHistogram();
            //delete landau_nom;
        }
        else
        {
            h_nominal_unfolded_data  = sysMassFSRUnfold[sysName].at(1)->GetOutput("hunfolded_mass_temp",0,0,"mass[UO];pt[UOC0]",kTRUE);
        }
    }

    c1=new TCanvas("c1", "c1", 50, 50, 800, 800);
    c1->cd();
    gStyle->SetOptStat(0);

    TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
    pad1->SetBottomMargin(0.);
    pad1->SetTopMargin(0.1);
    pad1->SetTicks(1);
    pad1->SetLogy();
    pad1->Draw();
    pad1->cd();

    h_nominal_unfolded_data->SetTitle("");
    h_nominal_unfolded_data->Draw("p9histe");
    

    h_nominal_unfolded_data->SetMarkerStyle(20);
    h_nominal_unfolded_data->SetMarkerSize(.7);
    h_nominal_unfolded_data->SetLineColor(kBlack);
    h_nominal_unfolded_data->SetMarkerSize(.7);
    h_nominal_unfolded_data->GetYaxis()->SetTitle("Events/bin");

    //h_nominal_unfolded_data->SetMinimum(10.);
    //h_nominal_unfolded_data->SetMaximum(9e9);
    //
    //if(sysName=="QED_FSR" && var=="Pt")
    //{
        //cout << "fit QED FSR systematic..." << endl;
        //landau_nom->SetParameters(2., 6., 2., 0.05, -0.01);
        //landau_nom->SetParLimits(0, 0.0, 5.);
        //landau_nom->SetParLimits(1, 2.0, 9.);
        //landau_nom->SetParLimits(2, 1., 10.0);
        //landau_nom->SetParLimits(3, -0.05, 0.12);
        //landau_nom->SetParLimits(4, -0.2, -0.001);
        //h_nominal_unfolded_data->Fit(landau_nom);
        //landau_nom->SetLineColor(kBlack);
        //landau_nom->SetLineStyle(1);
        //landau_nom->Draw("same");
        //TH1* h_ = landau_nom->GetHistogram();
        //cout << "mean: " << h_->GetMean() << endl;
    //}

    TH1* h_data_sys_temp;
    int sysSize;
    if(isDetector) sysSize = sysPtUnfold[sysName].size();
    else sysSize = sysPtFSRUnfold[sysName].size();

    for(int i = 0; i < sysSize; i++){

     if (sysName == "QED_FSR" && i == 1) break;
     TString isys;
     isys.Form("%d", i);

     TH1 * hsys_temp = NULL;

     if(isDetector){
         if(var=="Pt"){
             hsys_temp  = sysPtUnfold[sysName].at(i)->GetOutput("hunfolded_pt_temp_",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
         }
         else{
             hsys_temp  = sysMassUnfold[sysName].at(i)->GetOutput("hunfolded_mass_temp_",0,0,"mass[UO];pt[UOC0]",kTRUE);
         }
     }
     else{
         if(var=="Pt"){
            
            hsys_temp  = sysPtFSRUnfold[sysName].at(i)->GetOutput("hunfolded_pt_temp_",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
            doNorm(hsys_temp);

            //cout << "fit QED FSR systematic..." << endl;
            //TF1* landau_nom = new TF1("landau_nom","[0]*TMath::Landau(x,[1],[2]) + [3] * expo(x * [4])", 0., 100);
            //landau_nom->SetParameters(2., 6., 2., 0.05, -0.01);
            //landau_nom->SetParLimits(0, 0.0, 5.);
            //landau_nom->SetParLimits(1, 2.0, 9.);
            //landau_nom->SetParLimits(2, 1., 4.0);
            //landau_nom->SetParLimits(3, -0.05, 0.12);
            //landau_nom->SetParLimits(4, -0.2, -0.001);
            //hsys_temp->Fit(landau_nom);
            //TH1* h_ = landau_nom->GetHistogram();
            //delete landau_nom;
         }
         else{
             hsys_temp  = sysMassFSRUnfold[sysName].at(i)->GetOutput("hunfolded_mass_temp_",0,0,"mass[UO];pt[UOC0]",kTRUE);
         }
     }
     h_data_sys_temp = ((TH1F*)hsys_temp->Clone("sys_temp"));
     h_data_sys_temp->Draw("histsame");
     //h_data_sys_temp->SetLineColor(kBlack);
     h_data_sys_temp->SetLineColor(2+i);
     h_data_sys_temp->SetLineStyle(2);
        
    //if(sysName=="QED_FSR")
    //{
    //cout << "fit QED FSR systematic..." << endl;
    //landau_nom_->SetParameters(2., 6., 2., 0.05, -0.01);
    //landau_nom_->SetParLimits(0, 0.0, 5.);
    //landau_nom_->SetParLimits(1, 2.0, 9.);
    //landau_nom_->SetParLimits(2, 1., 10.0);
    //landau_nom_->SetParLimits(3, -0.05, 0.12);
    //landau_nom_->SetParLimits(4, -0.2, -0.001);
    //h_data_sys_temp->Fit(landau_nom_);
    //landau_nom->SetLineColor(kRed);
    //landau_nom_->SetLineStyle(2);
    //landau_nom_->Draw("same");
    //TH1* h_ = landau_nom_->GetHistogram();
    //cout << "mean: " << h_->GetMean() << endl;
    //}

     delete hsys_temp;
    }
    c1->cd();

    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.4);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.2);
    pad2->SetTicks(1);
    pad2->SetGridy(1);
    pad2->Draw();
    pad2->cd();

    for(int i = 0; i < sysSize; i++){

    if(sysName == "QED_FSR" && i == 1) break;

    TString isys;
    isys.Form("%d", i);

    TH1 * hsys_temp = NULL;

    if(isDetector){
        if(var=="Pt"){
            hsys_temp  = sysPtUnfold[sysName].at(i)->GetOutput("hunfolded_pt_temp_",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
        }
        else{
            hsys_temp  = sysMassUnfold[sysName].at(i)->GetOutput("hunfolded_mass_temp_",0,0,"mass[UO];pt[UOC0]",kTRUE);
        }
    }
    else{
        if(var=="Pt"){
            hsys_temp  = sysPtFSRUnfold[sysName].at(i)->GetOutput("hunfolded_pt_temp_",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
            doNorm(hsys_temp);
        }
        else{
            hsys_temp  = sysMassFSRUnfold[sysName].at(i)->GetOutput("hunfolded_mass_temp_",0,0,"mass[UO];pt[UOC0]",kTRUE);
        }
    }

    ratio = ((TH1F*)hsys_temp->Clone("ratio_temp"));
    ratio->Divide(h_nominal_unfolded_data);
    if(i==0 ){
        ratio->Draw("hist");
        if(sysName != "QED_FSR")ratio->GetYaxis()->SetTitle("Systematic/ Nominal input");
        else ratio->GetYaxis()->SetTitle("Photos/ Pythia");
        ratio->SetLineColor(2+i);
        ratio->SetMinimum(0.8);
        ratio->SetMaximum(1.2);
        ratio->SetTitle("");
        ratio->GetYaxis()->SetNdivisions(504);
        ratio->GetYaxis()->SetRangeUser(0.81, 1.19);
        ratio->SetLineStyle(2);
        ratio->GetXaxis()->SetTitleOffset(2.5);
        if(var=="Pt") ratio->GetXaxis()->SetTitle("p_{T} (GeV)");
        if(var=="Mass") ratio->GetXaxis()->SetTitle("mass (GeV)");
    }
    else{
        ratio->Draw("histsame");
        ratio->SetLineStyle(2);
    }
    delete hsys_temp;
    }

   CMS_lumi( c1, 4, 0 );
   c1->cd();
   if (isDetector) c1->SaveAs(outpdf+var+"_distribution_"+ibinMass+"_"+sysName+".pdf");
   else c1->SaveAs(outpdf+var+"distribution_pre_fsr_"+ibinMass+sysName+".pdf");

   delete h_nominal_unfolded_data;
   delete ratio;
   delete pad1;
   delete pad2;
   delete c1;
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

void ISRUnfold::drawUnfoldedHists(TString outpdf, TString var, int nthMassBin_, TString sysName, bool systematic, bool fullSys, bool doNormalisation)
{

    cout << "ISRUnfold::drawUnfoldedHists called" << endl;
    gROOT->SetBatch();
    gStyle->SetOptFit(0);

    setTDRStyle();
    writeExtraText = true;
    extraText  = "work in progress";

    bool combineMassBins = false;
    bool firstIter = true;

    const TUnfoldBinningV17* temp_binning_gen_pt = nomPtUnfold->GetOutputBinning("Gen_Pt");
    const TUnfoldBinningV17* temp_binning_rec_pt = nomPtUnfold->GetInputBinning("Rec_Pt");
    const TUnfoldBinningV17* temp_binning_gen_mass = nomMassUnfold->GetOutputBinning("Gen_Mass");
    const TUnfoldBinningV17* temp_binning_rec_mass = nomMassUnfold->GetInputBinning("Rec_Mass");

    // get mass bin definition from (pt, mass) bin definition
    const TVectorD* temp_tvecd = temp_binning_gen_pt->GetDistributionBinning(1);
    const Double_t* massBins = temp_tvecd->GetMatrixArray();
    int nMassBin = temp_tvecd->GetNrows() - 1;

    int lastBin = nthMassBin_ + 1;
    // if nthMassBin_ is larger than the number of mass bins, then draw one mass distribution including all mass bins
    if( nthMassBin_ >= nMassBin  && !doNormalisation)
    {
        lastBin = nMassBin;
        nthMassBin_ = 0;
        combineMassBins = true;
    }

    bool dilepPtShapeComparison = false;
    if(nthMassBin_ >= nMassBin  && doNormalisation)
    {
        lastBin = nMassBin;
        nthMassBin_ = 0;
        dilepPtShapeComparison = true;
    }

    TH1* hunfolded_data = NULL;
    TH1* hunfolded_mc = NULL;
    TH1* hfolded_data = NULL;
    TH1* hfolded_mc = NULL;

    TH1* hunfolded_sys_err = NULL;
    TH1* hunfolded_mc_sys_err = NULL;
    TH1* hfolded_sys_err = NULL;
    TH1* hfolded_mc_sys_err = NULL;

    TH1F *hunfolded_ratio = NULL;
    TH1F *hfolded_ratio = NULL;

    TH1F *hunfolded_ratio_sys_err = NULL;
    TH1F *hunfolded_ratio_sys_err_mc = NULL;
    TH1F *hfolded_ratio_sys_err = NULL;
    TH1F *hfolded_ratio_sys_err_mc = NULL;

    TFile* filein = new TFile(hist_file_path);

    TString PostFixForHist = "EEnominal";
    if(channel_name == "muon")
    {
        PostFixForHist = "MuMunominal";
    }

    TPad *pad1 = NULL;
    TPad *pad2 = NULL;
    TPad *pad3 = NULL;
    TLegend* leg_nom = NULL;
    TString cut_info;
    TString lepton_type;

    TLine *l_;
    double true_mc_mean = 1.;

    bool data_over_mc = false;
    for(int nthMassBin = nthMassBin_; nthMassBin < lastBin; nthMassBin++)
    {

        TString ibinMass;
        ibinMass.Form("%d", nthMassBin);
        TString collapse = "C";
        if(combineMassBins)
        {
            ibinMass = "full_range";
        }

        if(firstIter || dilepPtShapeComparison)
        {

            TString massLowEdge[5] = {"50.", "65.", "80.", "100.", "200."};
            TString massHighEdge[5] = {"65.", "80.", "100.", "200.", "350."};

            if(channel_name=="muon")
            {
                massLowEdge[0] = "40.";
                massLowEdge[1] = "60.";

                massHighEdge[0] = "60.";
                massHighEdge[1] = "80.";
            }
            TString genHistPath;

            // get nominal unfoled result
            if(var == "Pt" )
            {
                // FSR unfolding
                {
                    bool useAxis = kTRUE;
                    if(combineMassBins)
                    {
                        ibinMass = "";
                        collapse = "";
                        useAxis = kFALSE;

                    }

                    // get data distribution
                    hunfolded_data  = nomPtFSRUnfold->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UO"+collapse+ibinMass+"]",useAxis);
                    hfolded_data  = nomPtFSRUnfold->GetInput("hdata_pt_temp",0,0,"pt[UO];mass[UO"+collapse+ibinMass+"]",useAxis);
                    if(!combineMassBins && doNormalisation)
                    {
                        doNorm(hunfolded_data);
                        doNorm(hfolded_data);
                    }
                    hunfolded_sys_err= ((TH1F*)hunfolded_data->Clone("sysErr"));
                    hfolded_sys_err= ((TH1F*)hfolded_data->Clone("sysErr_folded"));

                    // get truth MC distribution
                    if(!combineMassBins)
                    {
                        genHistPath = channel_name + "_full_phase_m" + massLowEdge[nthMassBin] + "to" + massHighEdge[nthMassBin] + "/ptll_preFSR/";
                        filein = new TFile(hist_file_path_DYHists);
                        hunfolded_mc = (TH1*)filein->Get(genHistPath + "histo_DYJets");
                        hunfolded_mc->SetDirectory(0);
                        hunfolded_mc->Add((TH1*)filein->Get(genHistPath + "histo_DYJets10to50"));
                        true_mc_mean = hunfolded_mc->GetMean();
                        if(doNormalisation)
                            doNorm(hunfolded_mc);

                        filein->Close();

                        genHistPath = channel_name + "_fiducial_phase_post_FSR_dRp1_m" + massLowEdge[nthMassBin] + "to" + massHighEdge[nthMassBin] + "/ptll/";
                        filein = new TFile(hist_file_path_DYHists);
                        hfolded_mc = (TH1*)filein->Get(genHistPath + "histo_DYJets");
                        hfolded_mc->SetDirectory(0);
                        hfolded_mc->Add((TH1*)filein->Get(genHistPath + "histo_DYJets10to50"));
                        if(doNormalisation)
                            doNorm(hfolded_mc);

                        filein->Close();
                    }
                    else
                    {
                        // create the same histograms to save MC disbributions
                        hunfolded_mc = ((TH1F*)hunfolded_data->Clone("hunfolded_mc"));
                        hfolded_mc = ((TH1F*)hfolded_data->Clone("hfolded_mc"));

                        filein = new TFile(hist_file_path_DYHists);
                        int ibin_global_pre_fsr = 1;
                        int ibin_global_dressed = 1;
                        for(int imass = 0; imass < lastBin; imass++)
                        {

                            TH1* temp_hist = NULL;
                            genHistPath = channel_name + "_full_phase_m" + massLowEdge[imass] + "to" + massHighEdge[imass] + "/ptll_preFSR/";
                            temp_hist = (TH1*)filein->Get(genHistPath + "histo_DYJets");
                            temp_hist->SetDirectory(0);
                            temp_hist->Add((TH1*)filein->Get(genHistPath + "histo_DYJets10to50"));

                            for(int ibin = 1; ibin<temp_hist->GetNbinsX()+1;ibin++)
                            {
                                hunfolded_mc->SetBinContent(ibin_global_pre_fsr, temp_hist->GetBinContent(ibin));
                                hunfolded_mc->SetBinError(ibin_global_pre_fsr, temp_hist->GetBinError(ibin));

                                ibin_global_pre_fsr++;
                            }
                            delete temp_hist;

                            genHistPath = channel_name + "_fiducial_phase_post_FSR_dRp1_m" + massLowEdge[imass] + "to" + massHighEdge[imass] + "/ptll/";
                            temp_hist = (TH1*)filein->Get(genHistPath + "histo_DYJets");
                            temp_hist->SetDirectory(0);
                            temp_hist->Add((TH1*)filein->Get(genHistPath + "histo_DYJets10to50"));

                            for(int ibin = 1; ibin<temp_hist->GetNbinsX()+1;ibin++)
                            {
                                hfolded_mc->SetBinContent(ibin_global_dressed, temp_hist->GetBinContent(ibin));
                                hfolded_mc->SetBinError(ibin_global_dressed, temp_hist->GetBinError(ibin));

                                ibin_global_dressed++;
                            }

                        }
                    }

                    hunfolded_mc_sys_err = ((TH1F*)hunfolded_mc->Clone("mcsysErr"));
                    hfolded_mc_sys_err = ((TH1F*)hfolded_mc->Clone("mcsysErr_folded"));
                }
            }

            if(var == "Mass" )
            {
                {
                    // get data distribution
                    hunfolded_data  = nomMassFSRUnfold->GetOutput("hunfolded_mass_temp",0,0,"mass[UO];pt[UOC0]",kTRUE);
                    if(!combineMassBins)
                        hunfolded_data->GetXaxis()->SetRange(hunfolded_data->GetXaxis()->FindBin(massBins[nthMassBin]+0.01),hunfolded_data->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));
                    hunfolded_sys_err= ((TH1F*)hunfolded_data->Clone("sysErr"));

                    hfolded_data = nomMassFSRUnfold->GetInput("hdata_mass_temp",0,0,"mass[UO];pt[UOC0]",kTRUE);
                    if(!combineMassBins)
                        hfolded_data->GetXaxis()->SetRange(hfolded_data->GetXaxis()->FindBin(massBins[nthMassBin]+0.01),hfolded_data->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));
                    hfolded_sys_err= ((TH1F*)hfolded_data->Clone("sysErr_folded"));

                    // get truth MC distribution
                    if(!combineMassBins)
                    {
                        genHistPath = channel_name + "_full_phase_m" + massLowEdge[nthMassBin] + "to" + massHighEdge[nthMassBin] + "/mll_preFSR_" + channel_name + "/";
                        filein = new TFile(hist_file_path_DYHists);
                        hunfolded_mc = (TH1*)filein->Get(genHistPath + "histo_DYJets");
                        hunfolded_mc->SetDirectory(0);
                        hunfolded_mc->Add((TH1*)filein->Get(genHistPath + "histo_DYJets10to50"));

                        filein->Close();

                        genHistPath = channel_name + "_fiducial_phase_post_FSR_dRp1_m" + massLowEdge[nthMassBin] + "to" + massHighEdge[nthMassBin] + "/mll_" + channel_name + "/";
                        filein = new TFile(hist_file_path_DYHists);
                        hfolded_mc = (TH1*)filein->Get(genHistPath + "histo_DYJets");
                        hfolded_mc->SetDirectory(0);
                        hfolded_mc->Add((TH1*)filein->Get(genHistPath + "histo_DYJets10to50"));

                        filein->Close();
                    }
                    else
                    {
                        filein = new TFile(hist_file_path_DYHists);
                        for(int imass = 0; imass < lastBin; imass++)
                        {
                            if(imass == 0)
                            {
                                genHistPath = channel_name + "_full_phase_m" + massLowEdge[imass] + "to" + massHighEdge[imass] + "/mll_preFSR_" + channel_name + "/";
                                hunfolded_mc = (TH1*)filein->Get(genHistPath + "histo_DYJets");
                                hunfolded_mc->SetDirectory(0);
                                hunfolded_mc->Add((TH1*)filein->Get(genHistPath + "histo_DYJets10to50"));

                                genHistPath = channel_name + "_fiducial_phase_post_FSR_dRp1_m" + massLowEdge[imass] + "to" + massHighEdge[imass] + "/mll_" + channel_name + "/";
                                hfolded_mc = (TH1*)filein->Get(genHistPath + "histo_DYJets");
                                hfolded_mc->SetDirectory(0);
                                hfolded_mc->Add((TH1*)filein->Get(genHistPath + "histo_DYJets10to50"));
                            }
                            else
                            {
                                genHistPath = channel_name + "_full_phase_m" + massLowEdge[imass] + "to" + massHighEdge[imass] + "/mll_preFSR_" + channel_name + "/";
                                hunfolded_mc->Add((TH1*)filein->Get(genHistPath + "histo_DYJets"));
                                hunfolded_mc->Add((TH1*)filein->Get(genHistPath + "histo_DYJets10to50"));

                                genHistPath = channel_name + "_fiducial_phase_post_FSR_dRp1_m" + massLowEdge[imass] + "to" + massHighEdge[imass] + "/mll_" + channel_name + "/";
                                hfolded_mc->Add((TH1*)filein->Get(genHistPath + "histo_DYJets"));
                                hfolded_mc->Add((TH1*)filein->Get(genHistPath + "histo_DYJets10to50"));
                            }
                        }
                        filein->Close();
                    }

                    if(!combineMassBins)
                    {
                        hunfolded_mc->GetXaxis()->SetRange(hunfolded_data->GetXaxis()->FindBin(massBins[nthMassBin]+0.01),hunfolded_data->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));
                        hfolded_mc->GetXaxis()->SetRange(hunfolded_data->GetXaxis()->FindBin(massBins[nthMassBin]+0.01),hunfolded_data->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));
                    }
                    hunfolded_mc_sys_err = ((TH1F*)hunfolded_mc->Clone("mcsysErr"));
                    hfolded_mc_sys_err = ((TH1F*)hfolded_mc->Clone("mcsysErr_folded"));
                }
            }

            if(data_over_mc)
            {
                hunfolded_ratio = (TH1F*)hunfolded_data->Clone("hunfolded_ratio");
                hfolded_ratio = (TH1F*)hfolded_data->Clone("hfolded_ratio");
                hunfolded_ratio_sys_err= (TH1F*)hunfolded_data->Clone("ratio_sys");
                hunfolded_ratio->Divide(hunfolded_mc);
                hfolded_ratio->Divide(hfolded_mc);
                hunfolded_ratio_sys_err->Divide(hunfolded_mc);
                hunfolded_ratio_sys_err_mc= (TH1F*)hunfolded_ratio_sys_err->Clone("ratio_sys_mc");
            }
            else
            {
                hunfolded_ratio = (TH1F*)hunfolded_mc->Clone("hunfolded_ratio");
                hfolded_ratio = (TH1F*)hfolded_mc->Clone("hfolded_ratio");
                hunfolded_ratio_sys_err= (TH1F*)hunfolded_mc->Clone("ratio_sys");
                hfolded_ratio_sys_err= (TH1F*)hfolded_mc->Clone("ratio_sys_folded");

                hunfolded_ratio->Divide(hunfolded_data);
                hfolded_ratio->Divide(hfolded_data);
                hunfolded_ratio_sys_err->Divide(hunfolded_data);
                hfolded_ratio_sys_err->Divide(hfolded_data);

                // just to create histograms with the same bin definition
                hunfolded_ratio_sys_err_mc= (TH1F*)hunfolded_ratio_sys_err->Clone("ratio_sys_mc");
                hfolded_ratio_sys_err_mc= (TH1F*)hfolded_ratio_sys_err->Clone("ratio_sys_mc_folded");
            }
        }

        if(firstIter || dilepPtShapeComparison)
        {

            if(firstIter)
            {
                c1=new TCanvas("c1", "c1", 50, 50, 800, 800);
                gStyle->SetOptStat(0);
                gStyle->SetOptFit(0);
                c1->cd();

                pad1 = new TPad("pad1","pad1",0,0.4,1,1);
                pad1->SetBottomMargin(0.);
                pad1->SetTopMargin(0.1);
                pad1->SetTicks(1);
                pad1->SetLogy();
                pad1->Draw();
                pad1->cd();
            }
            else
            {
                c1->cd();
                pad1->cd();
            }

            TString draw_option = "p9e";
            int color = kRed;
            if(dilepPtShapeComparison)
            {
                draw_option = "p9samee";
                if(nthMassBin==1)
                    color = 632;
                if(nthMassBin==2)
                    color = 600;
                if(nthMassBin==3)
                    color = 616;
                if(nthMassBin==4)
                    color = 800;
            }

            hunfolded_data->SetTitle("");
            hunfolded_data->Draw(draw_option);
            hunfolded_mc->Draw("histsamee");

            hunfolded_mc->SetMarkerStyle(20);
            hunfolded_mc->SetMarkerSize(1.);
            hunfolded_mc->SetLineColor(color);
            hunfolded_mc->SetMarkerColor(color);

            if(!dilepPtShapeComparison)
            {
                hfolded_data->Draw("p9samee");
                hfolded_mc->Draw("histsamee");
                hfolded_mc->SetMarkerStyle(21);
                hfolded_mc->SetMarkerSize(1.);
                hfolded_mc->SetLineColor(kBlue);
                hfolded_mc->SetMarkerColor(kBlue);

                hfolded_data->SetMarkerColor(kBlack);
                hfolded_data->SetLineColor(kBlack);
                hfolded_data->SetMarkerSize(1.5);
                hfolded_data->SetMarkerStyle(25);
            }
            hunfolded_data->SetMarkerStyle(24);
            hunfolded_data->SetMarkerSize(1.5);
            hunfolded_data->SetLineColor(kBlack);
            hunfolded_data->GetYaxis()->SetTitle("Events/bin");

            if(var == "Pt" && doNormalisation)
            {
                TF1* landau = new TF1("landau","[0]*TMath::Landau(x,[1],[2]) + [3] * expo(x * [4])", 0., 100);
                landau->SetParameters(2., 6., 2., 0.05, -0.01);
                landau->SetParLimits(0, 0.0, 5.);
                landau->SetParLimits(1, 2.0, 9.);
                landau->SetParLimits(2, 1., 4.0);
                landau->SetParLimits(3, -0.05, 0.12);
                landau->SetParLimits(4, -0.2, -0.001);

                TF1* landau_mc = new TF1("landau_mc","[0]*TMath::Landau(x,[1],[2]) + [3] * expo(x * [4])", 0., 100);
                landau_mc->SetParameters(2., 6., 2., 0.05, -0.01);
                landau_mc->SetParLimits(0, 0.0, 5.);
                landau_mc->SetParLimits(1, 2.0, 9.);
                landau_mc->SetParLimits(2, 1., 4.0);
                landau_mc->SetParLimits(3, -0.05, 0.12);
                landau_mc->SetParLimits(4, -0.2, -0.001);

                landau->SetLineColor(kBlack);
                landau->SetLineStyle(1);
                hunfolded_data->Fit(landau);

                landau_mc->SetLineColor(kBlack);
                landau_mc->SetLineStyle(2);
                hunfolded_mc->Fit(landau_mc);

                TH1* h = landau->GetHistogram();
                TH1* h_mc = landau_mc->GetHistogram();

                TString mean_fit;
                TString mean_hist;
                mean_fit.Form("%.2f", h->GetMean());
                mean_hist.Form("%.2f", meanPt_data_pre_fsr.at(nthMassBin_));

                TString mean_fit_mc;
                TString mean_hist_mc;
                TString mean_true_mc;

                TString err_fit_mc;
                TString err_hist_mc;

                mean_fit_mc.Form("%.2f", h_mc->GetMean());
                mean_hist_mc.Form("%.2f", meanPt_mc_pre_fsr.at(nthMassBin_));
                mean_true_mc.Form("%.2f", true_mc_mean);

                err_fit_mc.Form("%.2f", fabs(h_mc->GetMean()-true_mc_mean)/ true_mc_mean * 100.);
                err_hist_mc.Form("%.2f",fabs(meanPt_mc_pre_fsr.at(nthMassBin_)-true_mc_mean)/ true_mc_mean * 100.);

                if(!dilepPtShapeComparison)
                {
                TLatex mean_true;
                mean_true.SetTextFont(63);
                mean_true.SetTextSize(20);
                mean_true.SetTextColor(kRed);
                mean_true.DrawLatexNDC(0.2, 0.45, "True Mean MC: " + mean_true_mc);

                TLatex mean_values;
                mean_values.SetTextFont(63);
                mean_values.SetTextSize(20);
                mean_values.DrawLatexNDC(0.2, 0.15, "#splitline{Mean (hist): " + mean_hist + "}{Mean (fit): " + mean_fit + "}");

                TLatex mean_values_mc;
                mean_values_mc.SetTextFont(63);
                mean_values_mc.SetTextSize(20);
                mean_values_mc.DrawLatexNDC(0.2, 0.35, "#splitline{Mean (hist): " + mean_hist_mc+ "(" + err_hist_mc +"%)}{Mean (fit): " + mean_fit_mc + "(" + err_fit_mc +"%)}");
                }
            }


            if(doNormalisation)
            {
                hunfolded_data->SetMinimum(5e-4);
                hunfolded_data->SetMaximum(0.9);
            }
            else
            {
                hunfolded_data->SetMinimum(5.);
                hunfolded_data->SetMaximum(1e10);
            }


            TString mean_nom;
            mean_nom.Form("%.5f", hunfolded_data->GetMean());

            TLegend* leg_nom = new TLegend(0.55, 0.65, 0.9, 0.85,"","brNDC");
            leg_nom->SetTextFont(63);
            leg_nom->SetTextSize(23);
            //leg_nom->SetTextSize(0.04);
            leg_nom->SetFillStyle(0);
            leg_nom->SetBorderSize(0);
            //leg_nom->SetNColumns(2);
            //leg_nom->AddEntry(hunfolded_data, "#splitline{Unfolded data}{(mean: " + mean_nom + ")}", "pl");
            {
                leg_nom->AddEntry(hunfolded_data, "FSR unfolded data", "pl");
                leg_nom->AddEntry(hunfolded_mc, "MC", "pl");

                if(!dilepPtShapeComparison)
                {
                    leg_nom->AddEntry(hfolded_data, "Detector unfolded data ", "pl");
                    leg_nom->AddEntry(hfolded_mc, "MC", "pl");
                }
            }
            leg_nom->Draw();

            if(channel_name=="electron") lepton_type = "ee";
            else lepton_type = "#mu#mu";

            if(var=="Pt" && !dilepPtShapeComparison)
            {
                TString low_bound_, upper_bound_;
                low_bound_.Form("%d", (int)massBins[nthMassBin]);
                upper_bound_.Form("%d", (int)massBins[nthMassBin+1]);
                cut_info = low_bound_ + " < M(" + lepton_type + ") < " + upper_bound_ + " (GeV)";
            }
            else
            {
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
            double chi2_without_correlation = Chi2Test(hfolded_data, hfolded_mc);
            bool printChi2 = false;

            if(printChi2)
            {
                {
                    after_unfolding_chi2.Form("%.2f", DoFit(var, nthMassBin, true));
                    before_unfolding_chi2.Form("%.2f", DoFit(var, nthMassBin));
                    chi2_.DrawLatexNDC(0.2, 0.15, "#splitline{Norm. #chi_{FSR unfolded}^{2}: " + after_unfolding_chi2 + "}{Norm. #chi_{detector unfolded}^{2}: " + before_unfolding_chi2 + "}");
                }
            }
        }

        /////////////////////////// systematics ////////////////////////////
        if(systematic)
        {
            // initialize errors the first iteration
            if(firstIter)
            {
                for(int ibin = 1; ibin<hunfolded_sys_err->GetNbinsX()+1;ibin++)
                {

                    hunfolded_sys_err->SetBinError(ibin, 0.);
                    hunfolded_mc_sys_err->SetBinError(ibin, 0.);
                    hunfolded_ratio_sys_err->SetBinContent(ibin, 1.);
                    hunfolded_ratio_sys_err->SetBinError(ibin, 0.);
                    hunfolded_ratio_sys_err_mc->SetBinError(ibin, 0.);

                    hfolded_sys_err->SetBinError(ibin, 1e-7);
                    hfolded_mc_sys_err->SetBinError(ibin, 1e-7);
                    hfolded_ratio_sys_err->SetBinContent(ibin, 1.);
                    hfolded_ratio_sys_err->SetBinError(ibin, 1e-7);
                    hfolded_ratio_sys_err_mc->SetBinError(ibin, 1e-7);
                }
            }
            //
            bool isFoldedSys = false;
            makeSystBand(var, nthMassBin, sysName, fullSys, data_over_mc,
                      hunfolded_data, hunfolded_mc, hunfolded_ratio,
                      hunfolded_sys_err, hunfolded_mc_sys_err, hunfolded_ratio_sys_err, hunfolded_ratio_sys_err_mc, isFoldedSys, combineMassBins, doNormalisation);

            isFoldedSys = true;
            if(sysName != "QED_FSR")
            {
                makeSystBand(var, nthMassBin, sysName, fullSys, data_over_mc,
                          hfolded_data, hfolded_mc, hfolded_ratio,
                          hfolded_sys_err, hfolded_mc_sys_err, hfolded_ratio_sys_err, hfolded_ratio_sys_err_mc, isFoldedSys, combineMassBins, doNormalisation);
            }

            // draw systematic envelope for a systematic source or all systematic sources
            hunfolded_sys_err->SetLineColor(12);
            hunfolded_sys_err->SetFillColor(12);
            hunfolded_sys_err->SetLineWidth(5);
            hunfolded_sys_err->SetFillStyle(3003);
            hunfolded_sys_err->SetMarkerSize(0);
            hunfolded_sys_err->Draw("E2 same");

            if(!dilepPtShapeComparison)
            {
                hfolded_sys_err->SetLineColor(12);
                hfolded_sys_err->SetFillColor(12);
                hfolded_sys_err->SetLineWidth(5);
                hfolded_sys_err->SetFillStyle(3003);
                hfolded_sys_err->SetMarkerSize(0);
                hfolded_sys_err->Draw("E2 same");
            }
        }
        //////////////////////////////// systematic /////////////////////////////////////

        c1->cd();

        // measurement distribution
        if(firstIter && !dilepPtShapeComparison)
        {
            pad2 = new TPad("pad2","pad2",0,0.25,1,0.4);
            pad2->SetTopMargin(0.05);
            pad2->SetBottomMargin(0.05);
            pad2->SetTicks(1);
            //pad2->SetGridy(1);
            pad2->Draw();
            pad2->cd();

            if(data_over_mc)
            {
                hfolded_ratio->Draw("pe");
                hfolded_ratio->SetLineColor(kBlack);
                hfolded_ratio->GetYaxis()->SetTitle("Data/ MC");
            }
            else
            {
                hfolded_ratio->Draw("pe");
                hfolded_ratio->SetMarkerColor(kBlue);
                hfolded_ratio->SetLineColor(kBlue);
                hfolded_ratio->SetMarkerSize(1.2);
                hfolded_ratio->SetMarkerStyle(21);
                hfolded_ratio->GetYaxis()->SetTitle("#splitline{  MC/Data}{(post FSR)}");
                hfolded_ratio->GetYaxis()->CenterTitle();
            }

            hfolded_ratio->SetMinimum(0.6);
            hfolded_ratio->SetMaximum(1.4);
            hfolded_ratio->SetLabelSize(0.);
            hfolded_ratio->GetYaxis()->SetLabelFont(43);
            hfolded_ratio->GetYaxis()->SetLabelSize(25);
            hfolded_ratio->GetYaxis()->SetNdivisions(504);
            hfolded_ratio->GetYaxis()->SetRangeUser(0.61, 1.39);
            hfolded_ratio->GetYaxis()->SetTitleFont(43);
            hfolded_ratio->GetYaxis()->SetTitleSize(25);
            hfolded_ratio->GetYaxis()->SetTitleOffset(2.0);

            if(systematic)
            {
                hfolded_ratio_sys_err->SetLineColor(12);
                hfolded_ratio_sys_err->SetFillColor(12);
                hfolded_ratio_sys_err->SetLineWidth(10);
                hfolded_ratio_sys_err->SetFillStyle(3004);
                hfolded_ratio_sys_err->SetMarkerSize(0);
                hfolded_ratio_sys_err->Draw("E2 same");
                hfolded_ratio_sys_err_mc->SetFillColorAlpha(kBlue,0.3);
                hfolded_ratio_sys_err_mc->SetMarkerSize(0);
                hfolded_ratio_sys_err_mc->Draw("E2 same");
            }

            l_ = new TLine(hunfolded_ratio->GetXaxis()->GetXmin(),1,hunfolded_ratio->GetXaxis()->GetXmax(),1);
            if(var=="Mass")
            {
                if(!combineMassBins)
                    l_ = new TLine(massBins[nthMassBin],1,massBins[nthMassBin+1],1);
                else
                    l_ = new TLine(hunfolded_ratio->GetXaxis()->GetXmin(),1,hunfolded_ratio->GetXaxis()->GetXmax(),1);
            }
            l_->Draw("same");
            l_->SetLineStyle(3);
            if(data_over_mc)
                l_->SetLineColor(kRed);
            else
                l_->SetLineColor(kBlack);
        }


        // unfolded distribution
        c1->cd();

        if(firstIter || dilepPtShapeComparison)
        {
            if(firstIter)
            {
                double y = 0.25;
                if(dilepPtShapeComparison)
                    y = 0.35;

                pad3 = new TPad("pad3","pad3",0,0.0,1, y);
                pad3->SetTopMargin(0.);
                pad3->SetBottomMargin(0.46); // 1. - 0.15 * 0.9 / 0.25
                pad3->SetTicks(1);
                //pad3->SetGridy(1);
                pad3->Draw();
                pad3->cd();
            }
            else
            {
                pad3->cd();
            }

            TString draw_option = "pe";
            if(dilepPtShapeComparison)
                draw_option = "pesame";

            int color = 1;
            if(dilepPtShapeComparison)
            {
                draw_option = "p9samee";
                if(nthMassBin==1)
                    color = 632;
                if(nthMassBin==2)
                    color = 600;
                if(nthMassBin==3)
                    color = 616;
                if(nthMassBin==4)
                    color = 800;
            }

            if(data_over_mc)
            {
                hunfolded_ratio->Draw(draw_option);
            }
            else
            {
                hunfolded_ratio->Draw(draw_option);
                hunfolded_ratio->SetMarkerColor(kRed);
                hunfolded_ratio->SetLineColor(kRed);
                hunfolded_ratio->SetMarkerSize(1.2);
                hunfolded_ratio->SetMarkerStyle(20);
                hunfolded_ratio->GetYaxis()->SetTitle("#splitline{ MC/Data}{(pre FSR)}");
                hunfolded_ratio->GetYaxis()->CenterTitle();
            }

            if(systematic)
            {
                hunfolded_ratio_sys_err->SetLineColor(color);
                hunfolded_ratio_sys_err->SetFillColor(color);
                hunfolded_ratio_sys_err->SetLineWidth(5);
                hunfolded_ratio_sys_err->SetFillStyle(3004);
                hunfolded_ratio_sys_err->SetMarkerSize(0);
                hunfolded_ratio_sys_err->Draw("E2 same");
                //hunfolded_ratio_sys_err->SetFillColorAlpha(kBlack,0.3);
                hunfolded_ratio_sys_err_mc->SetFillColorAlpha(kRed,0.3);
                hunfolded_ratio_sys_err_mc->SetMarkerSize(0);
                hunfolded_ratio_sys_err_mc->Draw("E2 same");
            }

            if(var=="Pt") hunfolded_ratio->GetXaxis()->SetTitle("p_{T}(" + lepton_type + ") (GeV)");
            if(var=="Mass") hunfolded_ratio->GetXaxis()->SetTitle("mass(" + lepton_type + ") (GeV)");

            hunfolded_ratio->SetMinimum(0.6);
            hunfolded_ratio->SetMaximum(1.4);
            hunfolded_ratio->GetYaxis()->SetLabelFont(43);
            hunfolded_ratio->GetYaxis()->SetLabelSize(25);
            hunfolded_ratio->GetYaxis()->SetNdivisions(504);
            hunfolded_ratio->GetYaxis()->SetRangeUser(0.61, 1.39);
            hunfolded_ratio->GetYaxis()->SetTitleFont(43);
            hunfolded_ratio->GetYaxis()->SetTitleSize(25);
            hunfolded_ratio->GetYaxis()->SetTitleOffset(2.0);
            hunfolded_ratio->GetXaxis()->SetLabelFont(43);
            hunfolded_ratio->GetXaxis()->SetLabelSize(25);
            hunfolded_ratio->GetXaxis()->SetTitleOffset(5.);
            hunfolded_ratio->GetXaxis()->SetTitleFont(43);
            hunfolded_ratio->GetXaxis()->SetTitleSize(25);

            l_->Draw("same");
            l_->SetLineStyle(2);
            if(data_over_mc) l_->SetLineColor(kRed);
            else l_->SetLineColor(kBlack);
            if(firstIter)
                CMS_lumi( pad1, 4, 0 );
        }

        c1->cd();

        if((firstIter && !combineMassBins && !dilepPtShapeComparison) || (combineMassBins && nthMassBin + 1 == lastBin) ||
            (dilepPtShapeComparison && nthMassBin + 1 == lastBin))
        {
            // create pdf if every things are done
            if(doNormalisation && !dilepPtShapeComparison)
                c1->SaveAs(outpdf+"_"+ibinMass+"_"+var+"_FSRUnfold_Norm.pdf");
            else if(doNormalisation && dilepPtShapeComparison)
                c1->SaveAs(outpdf+"_"+ibinMass+"_"+var+"_FSRUnfold_PtComparison.pdf");
            else
                c1->SaveAs(outpdf+"_"+ibinMass+"_"+var+"_FSRUnfold.pdf");

            //filein->Close();
            delete hfolded_mc;
            delete hunfolded_mc_sys_err;

            delete hfolded_data;
            delete hunfolded_ratio;
            delete hfolded_ratio;
            delete hunfolded_ratio_sys_err;
            delete hunfolded_ratio_sys_err_mc;

            delete hunfolded_data;
            delete hunfolded_sys_err;
            delete hunfolded_mc;
            delete pad1;
            delete pad2;
            delete c1;
            cout << "delete done" << endl;
        }

        firstIter = false;
    }
}

void ISRUnfold::makeSystBand(const TString var, const int nthMassBin, const TString sysName, const bool fullSys, const bool data_over_mc,
                  const TH1* hunfolded_data, const TH1* hunfolded_mc, const TH1* hunfolded_ratio,
                  TH1* hunfolded_sys_err, TH1* hunfolded_mc_sys_err, TH1* hunfolded_ratio_sys_err, TH1* hunfolded_ratio_sys_err_mc, bool isFoldedSys, bool isMassCombined, bool doNormalisation)
{
    const TUnfoldBinningV17* temp_binning_gen_pt = nomPtUnfold->GetOutputBinning("Gen_Pt");
    const TUnfoldBinningV17* temp_binning_rec_pt = nomPtUnfold->GetInputBinning("Rec_Pt");
    const TUnfoldBinningV17* temp_binning_gen_mass = nomMassUnfold->GetOutputBinning("Gen_Mass");
    const TUnfoldBinningV17* temp_binning_rec_mass = nomMassUnfold->GetInputBinning("Rec_Mass");

    // get mass bin definition from (pt, mass) bin definition
    const TVectorD* temp_tvecd = temp_binning_gen_pt->GetDistributionBinning(1);
    const Double_t* massBins = temp_tvecd->GetMatrixArray();

    // get number of pT bins
    const TVectorD* temp_tvecd_ = temp_binning_gen_pt->GetDistributionBinning(0);
    int nPtBin = temp_tvecd_->GetNrows() - 1;

    TString ibinMass;
    ibinMass.Form("%d", nthMassBin);

    std::map<TString, std::vector<TUnfoldDensityV17*>>::iterator it;
    std::map<TString, std::vector<TUnfoldDensityV17*>>::iterator it_end;

    if(!isFoldedSys)
    {
        it = sysPtFSRUnfold.begin();
        it_end = sysPtFSRUnfold.end();
    }
    else
    {
        it = sysPtUnfold.begin();
        it_end = sysPtUnfold.end();
    }

    while(it != it_end)
    {
        if(!fullSys && it->first != sysName)
        {
            it++;
            continue;
        }

        TH1 * hdatasys0_temp = NULL;
        TH1 * hdatasys1_temp = NULL;

        TH1 * hmcsys0_temp = NULL;
        TH1 * hmcsys1_temp = NULL;

        int systematic_variation_index = 0;
        int systematic_variation_index_mc = 0;

        std::vector<std::map<TString, Int_t>>* p_meanPtErrIdx_data;
        std::vector<std::map<TString, Int_t>>* p_meanPtErrIdx_mc;
        std::map<TString, std::vector<TUnfoldDensityV17*>> temp_sysTUnfDensity;

        // get distribution which gives the maximum variation on mean value
        if(var == "Pt")
        {
            if(!isFoldedSys)
            {
                p_meanPtErrIdx_data = &meanPtErrIdx_sysdata_pre_fsr;
                p_meanPtErrIdx_mc = &meanPtErrIdx_sysmc_pre_fsr;
                temp_sysTUnfDensity = sysPtFSRUnfold;
            }
            else
            {
                p_meanPtErrIdx_data = &meanPtErrIdx_sysdata_det_unf;
                p_meanPtErrIdx_mc = &meanPtErrIdx_sysmc_det_unf;
                temp_sysTUnfDensity = sysPtUnfold;
            }
            // note except QED_FSR, use the systematic histogram giving maximum variation on mean value
            // to calculate systematic error on each bin.
            // this is not correct for some systematic sources, for example PDF error.
            if(it->first != "QED_FSR")
            {
                systematic_variation_index = p_meanPtErrIdx_data->at(nthMassBin)[it->first];
                systematic_variation_index_mc = p_meanPtErrIdx_data->at(nthMassBin)[it->first];

                if(!isMassCombined)
                {
                    hdatasys0_temp = temp_sysTUnfDensity[it->first].at(systematic_variation_index)->GetOutput("hunfolded_pt_systemp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
                    hmcsys0_temp = temp_sysTUnfDensity[it->first].at(systematic_variation_index_mc)->GetBias("hmc_pt_systemp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
                    if(doNormalisation)
                    {
                        doNorm(hdatasys0_temp);
                        doNorm(hmcsys0_temp);
                    }
                }
                else
                {
                    hdatasys0_temp = temp_sysTUnfDensity[it->first].at(systematic_variation_index)->GetOutput("hunfolded_pt_systemp",0,0,"pt[UO];mass[UO]",kFALSE);
                    hmcsys0_temp = temp_sysTUnfDensity[it->first].at(systematic_variation_index_mc)->GetBias("hmc_pt_systemp",0,0,"pt[UO];mass[UO]",kFALSE);
                }

            }
            else
            {
                if(!isMassCombined)
                {
                    hdatasys0_temp = temp_sysTUnfDensity[it->first].at(0)->GetOutput("hunfolded_pt_systemp_fsr0",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
                    hdatasys1_temp = temp_sysTUnfDensity[it->first].at(1)->GetOutput("hunfolded_pt_systemp_fsr1",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);

                    hmcsys0_temp = temp_sysTUnfDensity[it->first].at(0)->GetBias("hmc_pt_systemp_fsr0",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
                    hmcsys1_temp = temp_sysTUnfDensity[it->first].at(1)->GetBias("hmc_pt_systemp_fsr1",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);


                    if(doNormalisation)
                    {
                        doNorm(hdatasys0_temp);
                        doNorm(hmcsys0_temp);
                        doNorm(hdatasys1_temp);
                        doNorm(hmcsys1_temp);
                    }
                }
                else
                {
                    hdatasys0_temp = temp_sysTUnfDensity[it->first].at(0)->GetOutput("hunfolded_pt_systemp_fsr0",0,0,"pt[UO];mass[UO]",kFALSE);
                    hdatasys1_temp = temp_sysTUnfDensity[it->first].at(1)->GetOutput("hunfolded_pt_systemp_fsr1",0,0,"pt[UO];mass[UO]",kFALSE);

                    hmcsys0_temp = temp_sysTUnfDensity[it->first].at(0)->GetBias("hmc_pt_systemp_fsr0",0,0,"pt[UO];mass[UO]",kFALSE);
                    hmcsys1_temp = temp_sysTUnfDensity[it->first].at(1)->GetBias("hmc_pt_systemp_fsr1",0,0,"pt[UO];mass[UO]",kFALSE);

                }
            }
        }

        if(var == "Mass")
        {
            if(!isFoldedSys)
            {
                p_meanPtErrIdx_data = &meanMassErrIdx_sysdata_pre_fsr;
                p_meanPtErrIdx_mc = &meanMassErrIdx_sysmc_pre_fsr;
                temp_sysTUnfDensity = sysMassFSRUnfold;
            }
            else
            {
                p_meanPtErrIdx_data = &meanMassErrIdx_sysdata_det_unf;
                p_meanPtErrIdx_mc = &meanMassErrIdx_sysmc_det_unf;
                temp_sysTUnfDensity = sysMassUnfold;
            }
            if(it->first != "QED_FSR")
            {
                systematic_variation_index = p_meanPtErrIdx_data->at(nthMassBin)[it->first];
                hdatasys0_temp = temp_sysTUnfDensity[it->first].at(systematic_variation_index)->GetOutput("hunfolded_mass_systemp",0,0,"mass[UO];pt[UOC0]",kTRUE);
                hdatasys0_temp->GetXaxis()->SetRange(hdatasys0_temp->GetXaxis()->FindBin(massBins[nthMassBin]+0.01), hdatasys0_temp->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));

                systematic_variation_index_mc = p_meanPtErrIdx_mc->at(nthMassBin)[it->first];
                hmcsys0_temp = temp_sysTUnfDensity[it->first].at(systematic_variation_index_mc)->GetBias("hmc_mass_systemp",0,0,"mass[UO];pt[UOC0]",kTRUE);
                hmcsys0_temp->GetXaxis()->SetRange(hdatasys0_temp->GetXaxis()->FindBin(massBins[nthMassBin]+0.01), hdatasys0_temp->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));

            }
            else
            {
                hdatasys0_temp = temp_sysTUnfDensity[it->first].at(0)->GetOutput("hunfolded_mass_systemp_fsr0",0,0,"mass[UO];pt[UOC0]",kTRUE);
                hdatasys1_temp = temp_sysTUnfDensity[it->first].at(1)->GetOutput("hunfolded_mass_systemp_fsr1",0,0,"mass[UO];pt[UOC0]",kTRUE);
                hdatasys0_temp->GetXaxis()->SetRange(hdatasys0_temp->GetXaxis()->FindBin(massBins[nthMassBin]+0.01), hdatasys0_temp->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));
                hdatasys1_temp->GetXaxis()->SetRange(hdatasys1_temp->GetXaxis()->FindBin(massBins[nthMassBin]+0.01), hdatasys1_temp->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));

                hmcsys0_temp = temp_sysTUnfDensity[it->first].at(0)->GetBias("hmc_mass_systemp_fsr0",0,0,"mass[UO];pt[UOC0]",kTRUE);
                hmcsys1_temp = temp_sysTUnfDensity[it->first].at(1)->GetBias("hmc_mass_systemp_fsr1",0,0,"mass[UO];pt[UOC0]",kTRUE);
                hmcsys0_temp->GetXaxis()->SetRange(hdatasys0_temp->GetXaxis()->FindBin(massBins[nthMassBin]+0.01), hdatasys0_temp->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));
                hmcsys1_temp->GetXaxis()->SetRange(hdatasys1_temp->GetXaxis()->FindBin(massBins[nthMassBin]+0.01), hdatasys1_temp->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));

            }
        }

        TH1F* ratio0_temp = NULL;
        TH1F* ratio1_temp = NULL;

        TH1F* ratio0_temp_mc = NULL;
        TH1F* ratio1_temp_mc = NULL;

        if(it->first != "QED_FSR")
        {
            if(data_over_mc)
            {
                ratio0_temp = (TH1F*)hdatasys0_temp->Clone("ratio");
                ratio0_temp->Divide(hunfolded_mc);

                ratio0_temp_mc = (TH1F*)hunfolded_data->Clone("ratio_mc");
                ratio0_temp_mc->Divide(hmcsys0_temp);
            }
            else
            {
                ratio0_temp = (TH1F*)hunfolded_mc->Clone("ratio");
                ratio0_temp->Divide(hdatasys0_temp);

                ratio0_temp_mc = (TH1F*)hmcsys0_temp->Clone("ratio_mc");
                ratio0_temp_mc->Divide(hunfolded_data);
            }
        }
        else
        {
            if(data_over_mc)
            {
                ratio0_temp = (TH1F*)hdatasys0_temp->Clone("ratio");
                ratio0_temp->Divide(hunfolded_mc);

                ratio1_temp = (TH1F*)hdatasys1_temp->Clone("ratio");
                ratio1_temp->Divide(hunfolded_mc);

                ratio0_temp_mc = (TH1F*)hunfolded_data->Clone("ratio_mc");
                ratio0_temp_mc->Divide(hmcsys0_temp);

                ratio1_temp_mc = (TH1F*)hunfolded_data->Clone("ratio_mc");
                ratio1_temp_mc->Divide(hmcsys1_temp);
            }
            else
            {
                ratio0_temp = (TH1F*)hunfolded_mc->Clone("ratio");
                ratio0_temp->Divide(hdatasys0_temp);

                ratio1_temp = (TH1F*)hunfolded_mc->Clone("ratio");
                ratio1_temp->Divide(hdatasys1_temp);

                ratio0_temp_mc = (TH1F*)hmcsys0_temp->Clone("ratio");
                ratio0_temp_mc->Divide(hunfolded_data);

                ratio1_temp_mc = (TH1F*)hmcsys1_temp->Clone("ratio");
                ratio1_temp_mc->Divide(hunfolded_data);

            }
        }

        // loop over each bin of unfolded histogram and set error due to the systematic source
        // note errors of hunfolded_sys_err initialized as zero before looping over all systematic sources
        // for mass, update only the specific bins in the the current mass range

        int firstBin = 1;
        int lastBin = hunfolded_sys_err->GetNbinsX()+1;
        if(var == "Mass")
        {
            firstBin = hdatasys0_temp->GetXaxis()->FindBin(massBins[nthMassBin]+0.01);
            lastBin = hdatasys0_temp->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01)+1;
        }
        if(var == "Pt" && isMassCombined)
        {
            firstBin = firstBin + nPtBin * nthMassBin;
            lastBin = firstBin + nPtBin;
        }

        for(int ibin = firstBin; ibin < lastBin; ibin++)
        {

            Double_t temp_err = 0.;
            Double_t temp_err_mc = 0.;
            Double_t temp_ratio_err = 0.;
            Double_t temp_ratio_err_mc = 0.;

            Double_t previous_err = hunfolded_sys_err->GetBinError(ibin);
            Double_t previous_err_mc = hunfolded_mc_sys_err->GetBinError(ibin);
            Double_t previous_ratio_err = hunfolded_ratio_sys_err->GetBinError(ibin);
            Double_t previous_ratio_err_mc = hunfolded_ratio_sys_err_mc->GetBinError(ibin);

            if(it->first != "QED_FSR")
            {
                temp_err =  fabs(hunfolded_data->GetBinContent(ibin) - hdatasys0_temp->GetBinContent(ibin));
                temp_ratio_err = fabs(hunfolded_ratio->GetBinContent(ibin) - ratio0_temp->GetBinContent(ibin));

                temp_err_mc =  fabs(hunfolded_mc->GetBinContent(ibin) - hmcsys0_temp->GetBinContent(ibin));
                temp_ratio_err_mc = fabs(hunfolded_ratio->GetBinContent(ibin) - ratio0_temp_mc->GetBinContent(ibin));
            }
            else
            {
                temp_err =  fabs(hdatasys0_temp->GetBinContent(ibin) - hdatasys1_temp->GetBinContent(ibin));
                temp_ratio_err = fabs(ratio0_temp->GetBinContent(ibin) - ratio1_temp->GetBinContent(ibin));

                temp_err_mc =  fabs(hmcsys0_temp->GetBinContent(ibin) - hmcsys1_temp->GetBinContent(ibin));
                temp_ratio_err_mc = fabs(ratio0_temp_mc->GetBinContent(ibin) - ratio1_temp_mc->GetBinContent(ibin));
            }

            hunfolded_sys_err->SetBinError(ibin, sqrt(pow(previous_err,2)+pow(temp_err,2)));
            hunfolded_mc_sys_err->SetBinError(ibin, sqrt(pow(previous_err_mc,2)+pow(temp_err_mc,2)));
            hunfolded_ratio_sys_err->SetBinContent(ibin, 1.);
            if(temp_ratio_err < 5.e-6) temp_ratio_err = 1.e-6;
            hunfolded_ratio_sys_err->SetBinError(ibin, sqrt(pow(previous_ratio_err,2)+pow(temp_ratio_err,2)));
            if(temp_ratio_err_mc < 5.e-6) temp_ratio_err_mc = 1.e-6;
            hunfolded_ratio_sys_err_mc->SetBinError(ibin, sqrt(pow(previous_ratio_err_mc,2)+pow(temp_ratio_err_mc,2)));
        }// loop for bin contents

        delete ratio0_temp;
        delete ratio1_temp;
        delete ratio0_temp_mc;
        delete ratio1_temp_mc;
        delete hdatasys0_temp;
        delete hdatasys1_temp;
        delete hmcsys0_temp;
        delete hmcsys1_temp;

        it++;
    }// end of while
}
