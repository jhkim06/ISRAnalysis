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

    // Bin definitions
    TUnfoldBinning* binning_Rec = NULL;
    TUnfoldBinning* binning_Gen = NULL;


    // Unfolding
    // For nominal results
    TUnfoldDensity* nominalTUnfold;
    TH2* hResponseM;

    // For statistical uncertainty
    std::vector<TUnfoldDensity*> UnfoldingInputStatTUnfold;
    std::vector<TUnfoldDensity*> UnfoldingMatrixStatTUnfold;
    TUnfoldDensity* modelUncertaintyTUnfold;
    TUnfoldDensity* ignoreBinZeroTUnfold;

    // For systematic uncertainty
    vector<TString> sysVector;
    std::map<TString, TUnfoldDensity*> systematicTUnfold;

    TUnfoldIterativeEM* iterEMTUnfold;

    int iBest;
    const int NITER_Iterative = 1000;

    TGraph *graph_SURE_IterativeSURE,*graph_DFdeviance_IterativeSURE;

    // Acceptance correction
    TH1* hFullPhaseData;
    TH1* hFullPhaseMC;

    std::map<TString, TH1*> hSysFullPhaseData;
    std::map<TString, TH1*> hSysFullPhaseMC;

    TH1* hAcceptance;
    std::map<TString, std::map<TString, TH1*>> hSysAcceptance;

    TFile* fUnfoldOut;
    TFile* fUnfoldReweight;
    TH2* hReweightSF;

    TString fUnfoldOutPath;

    // Conditions for unfolding
    TUnfold::ERegMode regMode;
    double nominal_bias;

    double tau;

    Int_t iBest_nominal;
    TSpline *logTauX,*logTauY;
    TGraph *lCurve;
    TGraph*bestLcurve;
    TGraph*bestLogTauLogChi2;

    TString output_baseDir;
    TString unfold_name;
    TString channel_name;
    TString var;
    int year;

    bool doInputStatUnc;
    bool doRMStatUnc;
    bool ignoreBinZero;
    bool doModelUnc;

public:

    // Constructor
    ISRUnfold(TString unfold_name_, TString channel, int year_ = 2016, int regMode_ = 0, bool doInputStatUnc_ = true, bool doRMStatUnc_ = false, bool ignoreBinZero_ = false, 
              bool verbose_ = false, TString var_ = "Mass", TString output_baseDir_ = "", TString matrix_reweight_file_ = "", bool doModelUnc_ = false)
    {
        cout << "ISRUnfold set! " << var_ << endl;

        output_baseDir = output_baseDir_;
        unfold_name = unfold_name_;
        channel_name = channel;
        year = year_;
        var = var_;

        nominal_bias = 1.;
        tau = 0.;

        if(regMode_ == 0)
            regMode = TUnfold::kRegModeNone;
        if(regMode_ == 1)
            regMode = TUnfold::kRegModeSize;
        if(regMode_ == 2)
            regMode = TUnfold::kRegModeDerivative;
        if(regMode_ == 3)
            regMode = TUnfold::kRegModeCurvature;

        doInputStatUnc = doInputStatUnc_;
        doRMStatUnc = doRMStatUnc_; 
        ignoreBinZero = ignoreBinZero_;
        doModelUnc = doModelUnc_;

        verbose = verbose_;

        // Make output root files
        TString yearStr;
        yearStr.Form("%d", (int)year);

        //output_baseDir = "output/" + yearStr + "/" + channel_name + "/";

        fUnfoldOutPath = output_baseDir+unfold_name+"_"+channel_name+"_"+yearStr+"_"+var+".root";
        fUnfoldOut = new TFile(fUnfoldOutPath, "RECREATE");

        fUnfoldOut->mkdir("matrix");
        fUnfoldOut->mkdir("folded");
        fUnfoldOut->mkdir("unfolded");
        fUnfoldOut->mkdir("acceptance");

        fUnfoldOut->mkdir("matrix/" + var);
        fUnfoldOut->mkdir("folded/" + var);
        fUnfoldOut->mkdir("unfolded/" + var);
        fUnfoldOut->mkdir("acceptance/" + var);

        // open reweight file
        if(doModelUnc && matrix_reweight_file_ != "")
        {
            fUnfoldReweight = new TFile(matrix_reweight_file_);
            if(var == "Pt")
            {
                hReweightSF = (TH2*) fUnfoldReweight->Get("reweighted_pt_matrix");
            }
            if(var == "Mass")
            {
                hReweightSF = (TH2*) fUnfoldReweight->Get("reweighted_mass_matrix");
            }
        }
    }
    // Destructor
    ~ISRUnfold()
    {
        fUnfoldOut->Close(); 
        if(doModelUnc)
        {
            fUnfoldReweight->Close();
        }
    }

    inline void closeOutFile()
    {
        fUnfoldOut->Close();
    }

    void checkIterEMUnfold(void);

    TMatrixD makeMatrixFromHist(TH2F*hist);
    void checkMatrixCond();
    double getSmearedChi2(TString filePath, TString dirName, TString steering, bool useAxis);
    double getUnfoldedChi2(TString steering, bool useAxis);

    void setBias(double bias);

    // Set nominal TUnfoldDensity
    void setNominalRM(TString filepath, TString dirName, TString binDef = "");
    void setFromPrevUnfResult(ISRUnfold* unfold, bool useAccept = false);

    // Set input histogram
    void setUnfInput(TString filepath = "", TString dirName = "", TString binDef = "", TString histName = "", TString sysType = "", TString sysName = "", TString histPostfix = "", bool isFSR = false);
    void setUnfInput(ISRUnfold* unfold,  TString thisSysType = "", TString sysName = "", bool useAccept = false);

    // Set background histograms
    void subBkgs(TString filepath, TString dirName = "",  TString binDef = "", TString bkgName = "", TString sysType = "", TString sysName = "", TString histPostfix = "");

    // Set systematic TUnfoldDensity
    void setSystematicRM(TString filepath, TString dirName,  TString binDef,  TString sysName, TString histPostfix);

    void setSystematics(TString sysHistName);

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

    void varyHistWithStatError(TH1* hist, int sys);

    // Do unfold
    void doISRUnfold();

    // Acceptance correction
    void doAcceptCorr(TString filePath, TString binDef, TString filePath_for_accept = "");

    // Get histograms
    TH1* getUnfoldedHists(TString outHistName = "", TString steering = "", bool useAxis = true);
    TH1* getRawHist(TString filePath, TString dirName, TString histName, TString outHistName, TString steering, bool useAxis);

};

#endif
