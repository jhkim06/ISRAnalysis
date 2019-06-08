#ifndef UNFOLDUTILS_H
#define UNFOLDUTILS_H

#include <iostream>
#include <map>
#include <cmath>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TDOMParser.h>
#include <TXMLDocument.h>

#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
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

	std::map<TString, std::vector<TUnfoldDensityV17*>> sysPtUnfold;
	std::map<TString, std::vector<TUnfoldDensityV17*>> sysMassUnfold;

	// nominal mean mass and pt for data and mc, and statistical and systematic errors 
	vector<Double_t> meanMass_data, meanMassStatErr_data, meanMassSysErr_data, meanMassTotErr_data;
	vector<Double_t> meanPt_data,   meanPtStatErr_data,   meanPtSysErr_data, meanPtTotErr_data;

	// nominal mean mass and pt for MC at pre FSR
	vector<Double_t> meanMass_mc, meanMassErr_mc;
	vector<Double_t> meanPt_mc,   meanPtErr_mc;

	vector<Double_t> meanMass_mcAlt, meanMassErr_mcAlt;
	vector<Double_t> meanPt_mcAlt,   meanPtErr_mcAlt;

	// map to save systematic uncertainties on mass and pt for each source 
	std::vector<std::map<TString, Double_t>>  meanMassErr_sysdata;
	std::vector<std::map<TString, Double_t>>  meanPtErr_sysdata;

	// map to save systematic variation of mean mass and pt for each source
	std::vector<std::map<TString, std::vector<Double_t>>>  meanMass_sysdata;
	std::vector<std::map<TString, std::vector<Double_t>>>  meanPt_sysdata;

	TCanvas* c1;

public:
	ISRUnfold() {}
	~ISRUnfold(){}

	// set nominal TUnfoldDensity 
	void setNomTUnfoldDensity(TString filepath, TString var, TString matrixName, bool isfsr = false);

	// set systematic TUnfoldDensity
	void setSysTUnfoldDensity(TString filepath, TString var, TString sysName, int nth, bool isfsr = false);

	// set input histogram
	void setInput(TString var, TString postfix, TString filepath, int nth = 0, bool isSys = false, double bias = 1.);

	// set background histograms
	void subBkgs(TString var, TString postfix, TString filepath, TString bkgName, int nth = 0, bool isSys = false);

	// set input data and background histograms for systematic TUnfoldDensity
	void initSysTUnfoldDensity();

	// do unfold 
	void doISRUnfold();

	// draw nominal detector level, unfolded plot, and final result plot without systematic
	void drawNominalPlots(TString outpdf, TString var = "Pt", int nthMassBin = 0, TString sysName = "");

	// draw input histograms using GetInput	
	void drawInputPlots(TString outpdf, TString var = "Pt", int nthMassBin = 0, TString sysName = "");

	void setMeanPt();
	void setMeanMass();

	// need unfolded hist, rho matrix (GetRhoIJtotal), MC truth
	void DoFit(TString var = "Pt", int nthMassBin = 0); // chi2 fit for unfolded distribution

	// draw ISR result
	void drawISRresult(TString outpdf);
};

#endif
