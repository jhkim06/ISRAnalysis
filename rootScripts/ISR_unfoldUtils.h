#ifndef UNFOLDUTILS_H
#define UNFOLDUTILS_H

#include <iostream>
#include <map>
#include <cmath>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TExec.h>
#include <TDOMParser.h>
#include <TXMLDocument.h>

#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TSpline.h>
#include <TGaxis.h>
#include <TMathText.h>
#include <TLatex.h>
#include <THStack.h>
#include <TMatrixDSymEigen.h>
#include "TUnfoldBinningXML.h"
#include "TUnfoldBinning.h"
#include "TUnfoldDensity.h"
#include "TRandom.h"
#include "TF1Convolution.h"

#include "TVectorD.h"
#include "TUnfoldBinningXML.h"
#include "TUnfoldBinning.h"
#include "TUnfoldDensity.h"

using namespace std; 

class ISRUnfold{

private:

    // Bin definitions
    TUnfoldBinning* pt_binning_Rec = NULL;
    TUnfoldBinning* pt_binning_Gen = NULL;
    TUnfoldBinning* mass_binning_Rec = NULL;
    TUnfoldBinning* mass_binning_Gen = NULL;

    // For nominal results
    TUnfoldDensityV17* nomPtUnfold;
    TUnfoldDensityV17* nomMassUnfold;

    TUnfoldDensityV17* nomPtFSRUnfold;
    TUnfoldDensityV17* nomMassFSRUnfold;
   
    std::map<TString, std::vector<TUnfoldDensityV17*>> sysPtUnfold;
    std::map<TString, std::vector<TUnfoldDensityV17*>> sysMassUnfold;
    
    std::map<TString, std::vector<TUnfoldDensityV17*>> sysPtFSRUnfold;
    std::map<TString, std::vector<TUnfoldDensityV17*>> sysMassFSRUnfold;
    
    // Results: mean values 
    // Detector level 
    vector<Double_t> meanMass_data_detector, meanMassStatErr_data_detector, meanMassSysErr_data_detector, meanMassTotErr_data_detector;
    vector<Double_t> meanPt_data_detector,   meanPtStatErr_data_detector,   meanPtSysErr_data_detector, meanPtTotErr_data_detector;
    
    // Detector unfolded results
    vector<Double_t> meanMass_data_det_unf, meanMassStatErr_data_det_unf, meanMassSysErr_data_det_unf, meanMassTotErr_data_det_unf;
    vector<Double_t> meanPt_data_det_unf,   meanPtStatErr_data_det_unf,   meanPtSysErr_data_det_unf, meanPtTotErr_data_det_unf;
    
    // Pre FSR unfolded results
    vector<Double_t> meanMass_data_pre_fsr, meanMassStatErr_data_pre_fsr, meanMassSysErr_data_pre_fsr, meanMassTotErr_data_pre_fsr;
    vector<Double_t> meanPt_data_pre_fsr,   meanPtStatErr_data_pre_fsr,   meanPtSysErr_data_pre_fsr, meanPtTotErr_data_pre_fsr;
    
    // Nominal mean mass and pt for MC
    vector<Double_t> meanMass_mc_det_unf, meanMassErr_mc_det_unf, meanMassStatErr_mc_det_unf, meanMassSysErr_mc_det_unf;
    vector<Double_t> meanPt_mc_det_unf,   meanPtErr_mc_det_unf, meanPtStatErr_mc_det_unf, meanPtSysErr_mc_det_unf;
    
    vector<Double_t> meanMass_mc_pre_fsr, meanMassStatErr_mc_pre_fsr, meanMassSysErr_mc_pre_fsr;
    vector<Double_t> meanPt_mc_pre_fsr,   meanPtStatErr_mc_pre_fsr, meanPtSysErr_mc_pre_fsr;
    
    TCanvas* c1;
    
    // 
    Int_t nScan;
    TSpline *rhoLogTau;
    TGraph *lCurve;
    Int_t iBest;
    
    Int_t nScan_mass;
    TSpline *rhoLogTau_mass;
    TGraph *lCurve_mass;
    Int_t iBest_mass;
    
    // Conditions for unfolding
    TUnfold::ERegMode regMode_detector;
    TUnfold::ERegMode regMode_FSR;
    double nominal_bias;
    
    TString output_baseDir;
    TString channel_name;
    int year;
        
public:
    
    // Constructor
    ISRUnfold(TString channel, int year_ = 2016, int regMode_detector_ = 0, int regMode_FSR_ = 0)
    {
        channel_name = channel;
        year = year_;

        nominal_bias = 1.; // 
        if(regMode_detector_ == 0)
            regMode_detector = TUnfold::kRegModeNone;
        if(regMode_detector_ == 1)
            regMode_detector = TUnfold::kRegModeSize;
        if(regMode_detector_ == 2)
            regMode_detector = TUnfold::kRegModeDerivative; 
        if(regMode_detector_ == 3)
            regMode_detector = TUnfold::kRegModeCurvature; 

        if(regMode_FSR_ == 0)
            regMode_FSR = TUnfold::kRegModeNone;
    }
    ~ISRUnfold(){}

    void setOutputBaseDir(TString outPath);
    void setBias(double bias);

    // Set nominal TUnfoldDensity 
    void setNomResMatrix(TString var, TString filepath, TString matrixName);
    // Set nominal TUnfoldDensity 
    void setNomFSRResMatrix(TString var, TString filepath, TString migrationName, TString phaseSpace);

    // Set input histogram
    void setFSRUnfInput(bool isSys = false, TString sysName = "", int nth = 0);
    void setUnfInput(TString var, TString filepath, TString dirName, bool isSys = false, TString sysName = "", int nth = 0);

    // Set background histograms
    void subBkgs(TString var, TString filepath, TString bkgName, bool isSys = false, TString sysName = "", int totSysN = -1, int nth = 0, TString phase_name = "full_phase");

    // Set systematic TUnfoldDensity
    void setSysTUnfoldDensity(TString var, TString filepath, TString sysName, int totSysN, int nth, TString phase_name = "full_phase", TString fsr_correction_name = "dressed_dRp1");
    void setSysFSRTUnfoldDensity(TString var, TString filepath, TString sysName, int totSysN, int nth, TString phase_name = "full_phase", TString fsr_correction_name = "dressed_dRp1");

    // Do unfold 
    void doISRUnfold(int detOrFSR_unfold = 0, bool doSys = false);

    // Get histograms
    TH1* getDetUnfoldedHists(TString var, TString outHistName = "", TString steering = "", bool useAxis = true);
    TH1* getFSRUnfoldedHists(TString var, TString outHistName = "", TString steering = "", bool useAxis = true);
    TH1* getMCHists(TString var, TString outHistName, TString steering, bool useAxis = true);
    TH1* getDetHists(TString var, TString outHistName = "", TString steering = "", bool useAxis = true);
    TH1* getRawHist(TString filePath, TString histName, TString outHistName, TString steering);

    // Helper functions
    void doNorm(TH1* hist, bool norm = true); 
    void drawtext(TGraph* g);

    int setMeanPt(bool isDet = true);
    int setMeanMass(bool isDet = true);
   
    // Get mean values 
    double getDetMeanPt(int ibin);
    double getDetMeanMass(int ibin);
    double getUnfMeanPt(int ibin);
    double getUnfMeanMass(int ibin);
    double getUnfMeanPtError(int ibin);
    double getUnfMeanMassError(int ibin);
    double getFSRUnfMeanPt(int ibin);
    double getFSRUnfMeanMass(int ibin);
    double getFSRUnfMeanPtError(int ibin);
    double getFSRUnfMeanMassError(int ibin);

    double getMCGenMeanMass(int ibin);
    double getMCGenMeanPt(int ibin);

    double DoFit(TString var = "Pt", int nthMassBin = 0, bool isFSRUnfold = false); // Chi2 fit for unfolded distribution

};

#endif
