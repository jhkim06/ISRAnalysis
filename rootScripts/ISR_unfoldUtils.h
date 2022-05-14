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
    
    // 
    bool doInputStatUnc;
    bool doRMStatUnc;
    bool ignoreBinZero;
    bool doModelUnc;

    TString unfoldName;
    TString channelName; // electron or muon
    TString var; // dilepton mass or pt
    int year;

    // Conditions for unfolding
    TUnfold::ERegMode regMode;
    TUnfoldDensity::EDensityMode densityMode;
    double nominalBias;
    double tau;

    // Bin definitions
    TUnfoldBinning* binningFine = NULL;
    TUnfoldBinning* binningCoarse = NULL;

    // Unfolding
    // For nominal results
    TUnfoldDensity* nominalTUnfold;
    TH2* hResponseM;

    // For unfolding uncertainty
    std::vector<TUnfoldDensity*> unfInputStatTUnfold;
    std::vector<TUnfoldDensity*> unfMatrixStatTUnfold;
    TUnfoldDensity* modelUncertaintyTUnfold;
    TUnfoldDensity* ignoreBinZeroTUnfold;

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

    TFile* fUnfoldReweight;
    TH2* hReweightSF;

    TString outputBaseDir;
    TFile* fUnfoldOut;
    TString fUnfoldOutPath;

public:

    // Constructor
    ISRUnfold(TString unfoldName_, TString channel, int year_ = 2016, TUnfold::ERegMode regMode_ = TUnfold::kRegModeNone, bool doInputStatUnc_ = true, bool doRMStatUnc_ = false, bool ignoreBinZero_ = false, 
              bool verbose_ = false, TString var_ = "Mass", TString outputBaseDir_ = "", TString matrix_reweight_file_ = "", bool doModelUnc_ = false, 
              TUnfoldDensity::EDensityMode densityMode_ = TUnfoldDensity::kDensityModeNone)
    {
        // modified in Xcode
        cout << "ISRUnfold set! " << var_ << endl;

        outputBaseDir = outputBaseDir_;
        unfoldName = unfoldName_;
        channelName = channel;
        year = year_;
        var = var_;

        nominalBias = 1.;
        tau = 0.;

        regMode = regMode_;
        densityMode = densityMode_;

        doInputStatUnc = doInputStatUnc_;
        doRMStatUnc = doRMStatUnc_; 
        ignoreBinZero = ignoreBinZero_;
        doModelUnc = doModelUnc_;

        verbose = verbose_;

        // Make output root files
        TString yearStr;
        yearStr.Form("%d", (int)year);

        //outputBaseDir = "output/" + yearStr + "/" + channelName + "/";
        fUnfoldOutPath = outputBaseDir+unfoldName+"_"+channelName+"_"+yearStr+"_"+var+".root";
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
    void setUnfInputUnfSys();

    // Set background histograms
    void subBkgs(TString filepath, TString dirName = "",  TString binDef = "", TString bkgName = "", TString sysType = "", TString sysName = "", TString histPostfix = "");

    // Set systematic TUnfoldDensity
    void setSystematicRM(TString filepath, TString dirName,  TString binDef,  TString sysName, TString histPostfix);
    // TODO TUnfoldSYS?

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
    void doISRUnfold(bool partialReg);
    void setPartialRegularize2D(TUnfold::ERegMode partialRegMode, double startMass, double startPt, double endMass, double endPt);

    // Acceptance correction
    void doAcceptCorr(TString filePath, TString binDef, TString filePath_for_accept = "");

    // Get histograms
    TH1* getUnfoldedHists(TString outHistName = "", TString steering = "", bool useAxis = true);
    TH1* getRawHist(TString filePath, TString dirName, TString histName, TString outHistName, TString steering, bool useAxis);

    // Bottom Line Test
    void drawChi2Reco();
    void drawChi2Truth();
};

#endif
