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

	TUnfoldDensityV17* nomPtFSRUnfold;
	TUnfoldDensityV17* nomMassFSRUnfold;

	std::map<TString, std::vector<TUnfoldDensityV17*>> sysPtUnfold;
	std::map<TString, std::vector<TUnfoldDensityV17*>> sysMassUnfold;
    
	std::map<TString, std::vector<TUnfoldDensityV17*>> sysPtFSRUnfold;
	std::map<TString, std::vector<TUnfoldDensityV17*>> sysMassFSRUnfold;

	// nominal mean mass and pt for data and mc, and statistical and systematic errors 
	vector<Double_t> meanMass_data, meanMassStatErr_data, meanMassSysErr_data, meanMassTotErr_data;
	vector<Double_t> meanPt_data,   meanPtStatErr_data,   meanPtSysErr_data, meanPtTotErr_data;

	vector<Double_t> meanMass_data_pre_fsr, meanMassStatErr_data_pre_fsr, meanMassSysErr_data_pre_fsr, meanMassTotErr_data_pre_fsr;
	vector<Double_t> meanPt_data_pre_fsr,   meanPtStatErr_data_pre_fsr,   meanPtSysErr_data_pre_fsr, meanPtTotErr_data_pre_fsr;

	vector<Double_t> meanMass_data_pre_fsr_dRp1, meanMassStatErr_data_pre_fsr_dRp1 ;
	vector<Double_t> meanPt_data_pre_fsr_dRp1,   meanPtStatErr_data_pre_fsr_dRp1;

	// nominal mean mass and pt for MC at pre FSR
	vector<Double_t> meanMass_mc, meanMassErr_mc;
	vector<Double_t> meanPt_mc,   meanPtErr_mc;

	vector<Double_t> meanMass_mc_pre_fsr, meanMassErr_mc_pre_fsr;
	vector<Double_t> meanPt_mc_pre_fsr,   meanPtErr_mc_pre_fsr;

	vector<Double_t> meanMass_mc_pre_fsr_, meanMassErr_mc_pre_fsr_;
	vector<Double_t> meanPt_mc_pre_fsr_,   meanPtErr_mc_pre_fsr_;

	vector<Double_t> meanMass_mcAlt, meanMassErr_mcAlt;
	vector<Double_t> meanPt_mcAlt,   meanPtErr_mcAlt;

	// map to save systematic uncertainties on mass and pt for each source 
	std::vector<std::map<TString, Double_t>>  meanMassErr_sysdata;
	std::vector<std::map<TString, Double_t>>  meanPtErr_sysdata;

	std::vector<std::map<TString, Double_t>>  meanMassErr_sysdata_pre_fsr;
	std::vector<std::map<TString, Double_t>>  meanPtErr_sysdata_pre_fsr;

	// map to save systematic variation of mean mass and pt for each source
	std::vector<std::map<TString, std::vector<Double_t>>>  meanMass_sysdata;
	std::vector<std::map<TString, std::vector<Double_t>>>  meanPt_sysdata;

	std::vector<std::map<TString, std::vector<Double_t>>>  meanMass_sysdata_pre_fsr;
	std::vector<std::map<TString, std::vector<Double_t>>>  meanPt_sysdata_pre_fsr;

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

        TString hist_file_path;
        TString channel_name;
public:
	ISRUnfold(TString channel, TString filepath){
            channel_name = channel;
            hist_file_path = filepath;
        }
	~ISRUnfold(){}

	// set nominal TUnfoldDensity 
	void setNomTUnfoldDensity(TString var, TString filepath, TString phase_name = "full_phase", TString fsr_correction_name = "dressed_dRp1");

	// set nominal TUnfoldDensity 
	void setNomFSRTUnfoldDensity(TString var, TString filepath, TString phase_name = "full_phase", TString fsr_correction_name = "dressed_dRp1");

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
	void doISRUnfold(bool doSys = false);
        void doISRQEDFSRUnfold(bool doSys = false);

        void drawClosurePlots(TString outpdf, TString var, int nthMassBin);

	// draw nominal detector level, unfolded plot, and final result plot without systematic
	void drawNominalPlots(TString outpdf, TString var = "Pt", int nthMassBin = 0, TString sysName = "", bool systematic = false, bool isFSRUnfold = false);

	// draw input histograms using GetInput	
	void drawInputPlots(TString outpdf, TString var = "Pt", int nthMassBin = 0, TString sysName = "");

	// draw reco histograms
	void drawNominalRecoPlots(TString outpdf, TString filepath, TString var = "Pt", int nthMassBin = 0, TString sysName = "");

	// draw nominal result with systematics 
	void drawSysPlots(TString outpdf, int nthMassBin = 0, TString sysName = "", bool detector_unfold = true);

	void studyFSRDRPlots(TString outpdf, TString var, int nthMassBin);

	void drawtext(TGraph* g);

        void setMCPreFSRMeanValues(TString filepath);
	void setMeanPt(bool doSys = false, bool altMC = false, bool detector_unfold = false);
	void setMeanMass(bool doSys = false, bool altMC = false, bool detector_unfold = false);

	// need unfolded hist, rho matrix (GetRhoIJtotal), MC truth
	double DoFit(TString var = "Pt", int nthMassBin = 0, bool isFSRUnfold = false); // chi2 fit for unfolded distribution
	double Chi2Test(TH1 *data, TH1 *mc);

	// draw ISR result
	void drawISRresult(TString outpdf, bool altMC = false, bool doFit = false);
	void drawLCurve(TString outpdf, TString var);
	void drawRhoLog(TString outpdf, TString var);
};

#endif
