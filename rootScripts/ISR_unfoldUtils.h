#ifndef UNFOLDUTILS_H
#define UNFOLDUTILS_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <cmath>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TProfile.h>
#include <TExec.h>
#include <TDOMParser.h>
#include <TXMLDocument.h>
#include <TSystem.h>

#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TFrame.h>
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
#include "TDecompSVD.h"

#include "TVectorD.h"
#include "TUnfoldBinningXML.h"
#include "TUnfoldBinning.h"
#include "TUnfoldDensity.h"
#include "TUnfoldIterativeEM.h"

using namespace std;
const int statSize = 100;

class ISRUnfold{

private:
    bool verbose;
    
    TString unfoldName;
    TString channel; // electron or muon
    TString var; // dilepton mass or pt
    TString year;

    TString folded_bin_name;
    TString unfolded_bin_name;

    // Conditions for unfolding
    TUnfold::ERegMode regMode;
    TUnfoldDensity::EDensityMode densityMode;
    double nominal_bias;
    double tau;

    // Bin definitions
    TUnfoldBinning* binningFine = NULL;
    TUnfoldBinning* binningCoarse = NULL;

    // Unfolding
    // For nominal results
    TUnfoldDensity* nominalTUnfold;
    TH2* hResponseM;

    // For systematic uncertainty
    vector<TString> sysVector;
    std::map<TString, TUnfoldDensity*> systematicTUnfold;

    int iBest;
    const int NITER_Iterative = 1000;
    TUnfoldIterativeEM* iterEMTUnfold;
    TGraph *graph_SURE_IterativeSURE,*graph_DFdeviance_IterativeSURE;

    Int_t iBest_nominal;
    TSpline *logTauX,*logTauY;
    TGraph *lCurve;
    TGraph*bestLcurve;
    TGraph*bestLogTauLogChi2;

    // Acceptance correction
    TH1* hFullPhaseData;
    TH1* hFullPhaseMC;

    std::map<TString, TH1*> hSysFullPhaseData;
    std::map<TString, TH1*> hSysFullPhaseMC;

    TH1* hAcceptance;
    std::map<TString, std::map<TString, TH1*>> hSysAcceptance;

    TString outputBaseDir;
    TFile* fUnfoldOut;
    TString fUnfoldOutPath;

public:

    // Constructor
    ISRUnfold(TString unfoldName_, 
              TString channel_, 
              TString year_ = "2016", 
              TString outputBaseDir_ = "", 
              bool verbose_ = false,
              TString var_ = "Mass", 
              TString folded_bin_name_ = "[folded_nominal:folded_nominal]",
              TString unfolded_bin_name_ = "[unfolded_nominal:unfolded_nominal]",
              TUnfold::ERegMode regMode_ = TUnfold::kRegModeNone, 
              TUnfoldDensity::EDensityMode densityMode_ = TUnfoldDensity::kDensityModeNone)
    {
        cout << "ISRUnfold set! " << var_ << endl;

        outputBaseDir = outputBaseDir_;
        unfoldName = unfoldName_;
        channel = channel_;
        year = year_;
        var = var_;

        folded_bin_name = folded_bin_name_;
        unfolded_bin_name = unfolded_bin_name_;

        nominal_bias = 1.;
        tau = 0.;

        regMode = regMode_;
        densityMode = densityMode_;

        verbose = verbose_;

        // Make output root files
        fUnfoldOutPath = outputBaseDir+unfoldName+"_"+channel+"_"+year+"_"+var+".root";
        fUnfoldOut = new TFile(fUnfoldOutPath, "RECREATE");

        fUnfoldOut->mkdir("matrix");
        fUnfoldOut->mkdir("folded");
        fUnfoldOut->mkdir("unfolded");
        fUnfoldOut->mkdir("acceptance");

        fUnfoldOut->mkdir("matrix/" + var);
        fUnfoldOut->mkdir("folded/" + var);
        fUnfoldOut->mkdir("unfolded/" + var);
        fUnfoldOut->mkdir("acceptance/" + var);
    }
    // Destructor
    ~ISRUnfold()
    {
        fUnfoldOut->Close(); 
    }

    inline void closeOutFile()
    {
        fUnfoldOut->Close();
    }

    void checkIterEMUnfold(void);

    TMatrixD makeMatrixFromHist(TH2F*hist);

    void setBias(double bias);

    // Set nominal TUnfoldDensity
    void setNominalRM(TString filepath, TString dirName);
    void save_hists_from_responseM();
    // Set systematic TUnfoldDensity
    void setSystematicRM(TString filepath, TString dirName,
                         TString sysName,
                         TString histPostfix);

    // Set input histogram
    void setUnfInput(TString filepath = "",
                     TString dirName = "",
                     TString sysType = "",
                     TString sysName = "",
                     TString sys_hist_postfix = "");

    // Set background histograms
    void subBkgs(TString filepath,
                 TString dirName = "",
                 TString bkgName = "",
                 TString sysType = "",
                 TString sysName = "",
                 TString histPostfix = "");

    void set_systematics(TString sysHistName);
    
    void inline printSystematics()
    {
        std::vector<TString>::iterator it = sysVector.begin();
        while(it != sysVector.end())
        {
            cout << "Systematic name: " << *it << endl;
            it++;
        }
    }
    inline std::vector<TString>& getSystematicVector()
    {
        return sysVector;
    }

    // Do unfold
    void doISRUnfold(bool partialReg);
    void setPartialRegularize2D(TUnfold::ERegMode partialRegMode,
                                double startMass, double startPt,
                                double endMass, double endPt);

    // Acceptance correction
    void doAcceptCorr(TString filePath, TString top_dir);

    // Get histograms
    TH1* getUnfoldedHists(TString outHistName = "", TString steering = "", bool useAxis = true);
};

#endif
