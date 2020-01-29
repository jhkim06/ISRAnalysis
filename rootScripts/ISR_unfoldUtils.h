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

#include "TVectorD.h"
#include "TUnfoldBinningXML.h"
#include "TUnfoldBinning.h"
#include "TUnfoldDensity.h"

using namespace std; 

class ISRUnfold{

private:
	TUnfoldDensityV17* nomPtUnfold;
	TUnfoldDensityV17* nomMassUnfold;

	TUnfoldDensityV17* nomPtUnfold_closure;
	TUnfoldDensityV17* nomMassUnfold_closure;

	TUnfoldDensityV17* nomPtFSRUnfold;
	TUnfoldDensityV17* nomMassFSRUnfold;

	TUnfoldDensityV17* nomPtFSRUnfold_closure;
	TUnfoldDensityV17* nomMassFSRUnfold_closure;

	std::map<TString, std::vector<TUnfoldDensityV17*>> sysPtUnfold;
	std::map<TString, std::vector<TUnfoldDensityV17*>> sysMassUnfold;
    
	std::map<TString, std::vector<TUnfoldDensityV17*>> sysPtFSRUnfold;
	std::map<TString, std::vector<TUnfoldDensityV17*>> sysMassFSRUnfold;

        // results
        // detector level 
	vector<Double_t> meanMass_data_detector, meanMassStatErr_data_detector, meanMassSysErr_data_detector, meanMassTotErr_data_detector;
	vector<Double_t> meanPt_data_detector,   meanPtStatErr_data_detector,   meanPtSysErr_data_detector, meanPtTotErr_data_detector;

        // detector unfolded results
	vector<Double_t> meanMass_data_det_unf, meanMassStatErr_data_det_unf, meanMassSysErr_data_det_unf, meanMassTotErr_data_det_unf;
	vector<Double_t> meanPt_data_det_unf,   meanPtStatErr_data_det_unf,   meanPtSysErr_data_det_unf, meanPtTotErr_data_det_unf;

        // pre FSR unfolded results
	vector<Double_t> meanMass_data_pre_fsr, meanMassStatErr_data_pre_fsr, meanMassSysErr_data_pre_fsr, meanMassTotErr_data_pre_fsr;
	vector<Double_t> meanPt_data_pre_fsr,   meanPtStatErr_data_pre_fsr,   meanPtSysErr_data_pre_fsr, meanPtTotErr_data_pre_fsr;

	// nominal mean mass and pt for MC at pre FSR
	vector<Double_t> meanMass_mc_det_unf, meanMassErr_mc_det_unf, meanMassStatErr_mc_det_unf, meanMassSysErr_mc_det_unf;
	vector<Double_t> meanPt_mc_det_unf,   meanPtErr_mc_det_unf, meanPtStatErr_mc_det_unf, meanPtSysErr_mc_det_unf;

	vector<Double_t> meanMass_mc_pre_fsr, meanMassStatErr_mc_pre_fsr, meanMassSysErr_mc_pre_fsr;
	vector<Double_t> meanPt_mc_pre_fsr,   meanPtStatErr_mc_pre_fsr, meanPtSysErr_mc_pre_fsr;

	vector<Double_t> meanMass_mcAlt, meanMassErr_mcAlt;
	vector<Double_t> meanPt_mcAlt,   meanPtErr_mcAlt;

	// map to save systematic uncertainties on mass and pt for each source 
	std::vector<std::map<TString, Double_t>>  meanMassErr_sysdata_det_unf;
	std::vector<std::map<TString, Double_t>>  meanPtErr_sysdata;

	std::vector<std::map<TString, Double_t>>  meanMassErr_sysmc_det_unf;
	std::vector<std::map<TString, Double_t>>  meanPtErr_sysmc_det_unf;

	std::vector<std::map<TString, Double_t>>  meanMassErr_sysdata_pre_fsr;
	std::vector<std::map<TString, Double_t>>  meanPtErr_sysdata_pre_fsr;

	std::vector<std::map<TString, Double_t>>  meanMassErr_sysmc_pre_fsr;
	std::vector<std::map<TString, Double_t>>  meanPtErr_sysmc_pre_fsr;

        // save systematic variation index result in maximum deviation from the nominal measurement
	std::vector<std::map<TString, Int_t>>  meanMassErrIdx_sysdata_det_unf;
	std::vector<std::map<TString, Int_t>>  meanPtErrIdx_sysdata_det_unf;

	std::vector<std::map<TString, Int_t>>  meanMassErrIdx_sysmc_det_unf;
	std::vector<std::map<TString, Int_t>>  meanPtErrIdx_sysmc_det_unf;

	std::vector<std::map<TString, Int_t>>  meanMassErrIdx_sysdata_pre_fsr;
	std::vector<std::map<TString, Int_t>>  meanPtErrIdx_sysdata_pre_fsr;

	std::vector<std::map<TString, Int_t>>  meanMassErrIdx_sysmc_pre_fsr;
	std::vector<std::map<TString, Int_t>>  meanPtErrIdx_sysmc_pre_fsr;

	// map to save systematic variation of mean mass and pt for each source
	std::vector<std::map<TString, std::vector<Double_t>>>  meanMass_sysdata_det_unf;
	std::vector<std::map<TString, std::vector<Double_t>>>  meanPt_sysdata_det_unf;

	std::vector<std::map<TString, std::vector<Double_t>>>  meanMass_sysmc_det_unf;
	std::vector<std::map<TString, std::vector<Double_t>>>  meanPt_sysmc_det_unf;

	std::vector<std::map<TString, std::vector<Double_t>>>  meanMass_sysdata_pre_fsr;
	std::vector<std::map<TString, std::vector<Double_t>>>  meanPt_sysdata_pre_fsr;

	std::vector<std::map<TString, std::vector<Double_t>>>  meanMass_sysmc_pre_fsr;
	std::vector<std::map<TString, std::vector<Double_t>>>  meanPt_sysmc_pre_fsr;

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

        // conditions for unfolding
        TUnfold::ERegMode regMode_detector;
        TUnfold::ERegMode regMode_FSR;
        double nominal_bias;
        
        TString hist_file_path;
        TString hist_file_path_DYDetResM;
        TString hist_file_path_DYFSRResM;
        TString output_baseDir;
        TString channel_name;
        bool do_normalization;
        int year;
        
public:

        // constructor
	ISRUnfold(TString channel, TString filepath, TString filepath_DYDetResM = "", TString filepath_DYFSRResM = "", bool norm = true, int year_ = 2016, 
                int regMode_detector_ = 0, int regMode_FSR_ = 0)
        {
            channel_name = channel;
            hist_file_path = filepath;
            hist_file_path_DYDetResM = filepath_DYDetResM;
            hist_file_path_DYFSRResM = filepath_DYFSRResM;
            do_normalization = norm;
            year = year_;

            nominal_bias = 1.;
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

	// set nominal TUnfoldDensity 
	void SetNomTUnfoldDensity(TString var, TString filepath, TString phase_name = "full_phase", TString fsr_correction_name = "dressed_dRp1");

	// set nominal TUnfoldDensity 
	void setNomFSRTUnfoldDensity(TString var, TString filepath, TString phase_name = "full_phase", TString fsr_correction_name = "dressed_dRp1");

        // do closure test: use the nominal probability matrix 
        void doClosureTest(int detOrFSR_unfold, TString filepath, TString phase_name = "full_phase");

	// set systematic TUnfoldDensity
	void setSysTUnfoldDensity(TString var, TString filepath, TString sysName, int totSysN, int nth, TString phase_name = "full_phase", TString fsr_correction_name = "dressed_dRp1");
        void setSysFSRTUnfoldDensity(TString var, TString filepath, TString sysName, int totSysN, int nth, TString phase_name = "full_phase", TString fsr_correction_name = "dressed_dRp1");

	// set input histogram
	void setFSRUnfoldInput(TString filepath, bool isSys = false, TString sysName = "", int nth = 0, TString phase_name = "full_phase");
	void setInput(TString var, TString filepath, bool isSys = false, TString sysName = "", int nth = 0, double bias = 1., TString phase_name = "full_phase");

	// set background histograms
	void subBkgs(TString var, TString filepath, TString bkgName, bool isSys = false, TString sysName = "", int totSysN = -1, int nth = 0, TString phase_name = "full_phase");

	// set input data and background histograms for systematic TUnfoldDensity
	void initSysTUnfoldDensity();

        // draw migration probability and efficiency plot
        void drawISRMatrixInfo(TString var, TString outpdf, bool detector_unfold = true, bool fsr_systematic = false);

	// do unfold 
	void doISRUnfold(int detOrFSR_unfold = 0, bool doSys = false);

        void drawClosurePlots(TString outpdf, TString var, int nthMassBin);

	// 
	void drawUnfoldedHists(TString outpdf, TString var = "Pt", int nthMassBin = 0, TString sysName = "", bool systematic = false, bool isFSRUnfold = false, bool fullSys = false);
        void makeSystBand(const TString var, const int nthMassBin, const TString sysName, const bool fullSys, const bool data_over_mc, 
                          const TH1* hunfolded_data, const TH1* hunfolded_mc, const TH1* hunfolded_ratio, 
                          TH1* hunfolded_sys_err, TH1* hunfolded_mc_sys_err, TH1* hunfolded_ratio_sys_err, TH1* hunfolded_ratio_sys_err_mc, bool isFoldedSys = false);

	// draw input histograms using GetInput	
	void drawInputPlots(TString outpdf, TString var = "Pt", int nthMassBin = 0, TString sysName = "");

	// draw reco histograms
	void drawNominalRecoPlots(TString outpdf, TString filepath, TString var = "Pt", int nthMassBin = 0, TString sysName = "");
        void drawSysComparionPlots(TString outpdf, TString var, int nthMassBin, TString sysName, bool isFSRUnfold);

	// draw nominal result with systematics 
	void drawSysPlots(TString outpdf, int nthMassBin = 0, TString sysName = "", bool detector_unfold = true);

	void drawtext(TGraph* g);

	void setMeanPt(bool doSys = false, bool altMC = false, bool detector_unfold = false);
	void setMeanMass(bool doSys = false, bool altMC = false, bool detector_unfold = false);

	// need unfolded hist, rho matrix (GetRhoIJtotal), MC truth
	double DoFit(TString var = "Pt", int nthMassBin = 0, bool isFSRUnfold = false); // chi2 fit for unfolded distribution
	double Chi2Test(TH1 *data, TH1 *mc);

	// draw ISR result
	TCanvas* drawISRresult(TString outpdf, bool altMC = false, bool doFit = false);
        void drawISRRun2results(TString outpdf, TCanvas* c_2017, TCanvas* c_2018, TCanvas* c_muon_2016, TCanvas* c_muon_2017, TCanvas* c_muon_2018);
	void drawLCurve(TString outpdf, TString var);
	void drawRhoLog(TString outpdf, TString var);

        void SavePtMassHists();
};

#endif
