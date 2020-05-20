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
const int statSize = 1000;

class ISRUnfold{

private:

    std::vector<double> massBinEdges;

    std::vector<TString> bkgNames; // Backgroud Name
    std::vector<TString> bkgTypes; // Backgroud Type

    std::map<TString, int> bkgColors;

    // Bin definitions
    TUnfoldBinning* pt_binning_Rec = NULL;
    TUnfoldBinning* pt_binning_Gen = NULL;
    TUnfoldBinning* mass_binning_Rec = NULL;
    TUnfoldBinning* mass_binning_Gen = NULL;

    std::map<TString, vector<TString>> sysMap;
    std::map<TString, TH1*> sysAbsHist;
    std::map<TString, TH1*> sysRelHist;

    // For nominal results
    TUnfoldDensityV17* nomPtUnfold;
    TUnfoldDensityV17* nomMassUnfold;

    // For statistical uncertainty
    std::vector<TUnfoldDensityV17*> statPtUnfold;
    std::vector<TUnfoldDensityV17*> statMassUnfold;

    // Histogram to save statistical variations for each mass bin
    std::vector<TH1*> meanPtStatVariation;
    std::vector<TH1*> meanMassStatVariation;

    // For systematic uncertainty
    std::map<TString, std::vector<TUnfoldDensityV17*>> sysPtUnfold;
    std::map<TString, std::vector<TUnfoldDensityV17*>> sysMassUnfold;
    
    // Results: mean values 
    // Detector level 
    vector<Double_t> meanMass_data_detector, meanMassStatErr_data_detector, meanMassSysErr_data_detector, meanMassTotErr_data_detector;
    vector<Double_t> meanPt_data_detector,   meanPtStatErr_data_detector,   meanPtSysErr_data_detector, meanPtTotErr_data_detector;
    
    // Detector unfolded results
    vector<Double_t> meanMass_data_unfoled, meanMassStatErr_data_unfoled, meanMassSysErr_data_unfoled, meanMassTotErr_data_unfoled;
    vector<Double_t> meanPt_data_unfoled,   meanPtStatErr_data_unfoled,   meanPtSysErr_data_unfoled, meanPtTotErr_data_unfoled;
    
    // Nominal mean mass and pt for MC
    vector<Double_t> meanMass_mc_unfoled, meanMassStatErr_mc_unfoled, meanMassSysErr_mc_unfoled, meanMassTotErr_mc_unfoled;
    vector<Double_t> meanPt_mc_unfoled, meanPtStatErr_mc_unfoled, meanPtSysErr_mc_unfoled, meanPtTotErr_mc_unfoled;
    
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
    TUnfold::ERegMode regMode;
    double nominal_bias;
    
    TString output_baseDir;
    TString channel_name;
    int year;

    void setMassBindEdges();
    bool makeStatUnfold; 
        
public:
    
    // Constructor
    ISRUnfold(TString channel, int year_ = 2016, int regMode_ = 0, bool makeStatUnfold_ = true)
    {
        cout << "ISRUnfold set!" << endl;

        channel_name = channel;
        year = year_;

        nominal_bias = 1.; // 

        if(regMode_ == 0)
            regMode = TUnfold::kRegModeNone;
        if(regMode_ == 1)
            regMode = TUnfold::kRegModeSize;
        if(regMode_ == 2)
            regMode = TUnfold::kRegModeDerivative; 
        if(regMode_ == 3)
            regMode = TUnfold::kRegModeCurvature; 

        makeStatUnfold = makeStatUnfold_;

        // Fill colors for backgrounds
        bkgColors["WJets"] = kViolet+1;
        bkgColors["EWK"] = kYellow+2;
        bkgColors["Top"] = kBlue;

    }
    ~ISRUnfold(){}

    void setOutputBaseDir(TString outPath);
    void setBias(double bias);

    // Set nominal TUnfoldDensity 
    void setNomResMatrix(TString var, TString filepath, TString dirName, TString histName, bool isSquareMatrix = false);

    // Set input histogram
    void setUnfInput(TString var, TString filepath, TString dirName, TString histName, bool isSys = false, TString sysName = "", int nth = 0);
    void setUnfInput(ISRUnfold* unfold, TString var, bool isSys = false, TString sysName = "", int nth = 0);

    // Set background histograms
    void subBkgs(TString filepath, std::pair<TString, TString>& bkgInfo, 
                 bool isSys = false, TString sysName = "", int totSysN = -1, int nth = 0, TString dirName = "full_phase", double bkgScale = 1.);

    // Set systematic TUnfoldDensity
    void setSysTUnfoldDensity(TString var, TString filepath, TString dirName, TString histName, TString sysName, TString sysPostfix, int nth);

    void setSystematics(TString sysName, TString sysHistName);
    void inline printSystematics()
    {
        std::map<TString, std::vector<TString>>::iterator it = sysMap.begin();
        while(it != sysMap.end())
        {
            cout << "Systematic name: " << it->first << endl;
            it++;
        }
    }
    // Draw folded distribution(before unfolding) using histograms saved in TUnfoldDensity
    TCanvas* drawFoldedHists(TString var, TString filePath, TString steering, bool useAxis, TString sysName = "");
    void setTHStack(TString var, TString filePath, THStack& hs, TH1& hMCtotal, TString steering, bool useAxis, TString sysName = "");

    // Do unfold 
    void doISRUnfold( bool doSys = false);
    void doStatUnfold(); 

    void setStatError();

    // Get histograms
    TH1* getDetUnfoldedHists(TString var, TString outHistName = "", TString steering = "", bool useAxis = true);
    TH1* getMCHists(TString var, TString outHistName, TString steering, bool useAxis = true);
    TH1* getDetHists(TString var, TString outHistName = "", TString steering = "", bool useAxis = true);
    TH1* getRawHist(TString var, TString filePath, TString dirName, TString histName, TString outHistName, TString steering, bool useAxis);

    // Helper functions
    void doNorm(TH1* hist, bool norm = true); 
    void drawtext(TGraph* g);

    int setMeanPt();
    int setMeanMass();

    void fillPtStatVariationHist(int istat);
    void fillMassStatVariationHist(int istat);
   
    // Get mean values 
    double getDetMeanPt(int ibin);
    double getDetMeanMass(int ibin);
    double getDetMeanPtError(int ibin);
    double getDetMeanMassError(int ibin);
    double getUnfMeanPt(int ibin);
    double getUnfMeanMass(int ibin);
    double getUnfMeanPtError(int ibin);
    double getUnfMeanMassError(int ibin);

    double getMCGenMeanMass(int ibin);
    double getMCGenMeanPt(int ibin);

    void drawStatVariation(bool isPt = true, int massBin = 0);

};

#endif
