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

const int statSize = 1000;

class ISRUnfold{

private:

    // Bin definitions
    TUnfoldBinning* binning_Rec = NULL;
    TUnfoldBinning* binning_Gen = NULL;

    vector<TString> sysVector;

    // Unfolding
    // For nominal results
    TUnfoldDensity* nominalTUnfold;
    TUnfoldDensity* nomMassUnfold;

    TH2* hResponseM;

    // For statistical uncertainty
    std::vector<TUnfoldDensity*> UnfoldingInputStatTUnfold;
    std::vector<TUnfoldDensity*> UnfoldingMatrixStatTUnfold;

    // For systematic uncertainty
    //std::map<TString, std::map<TString, TUnfoldDensity*>> systematicTUnfold;
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
    TH1* hAcceptanceFraction;
    std::map<TString, std::map<TString, TH1*>> hSysAcceptance;
    std::map<TString, std::map<TString, TH1*>> hSysAcceptanceFraction;

    TFile* fUnfoldOut;

    TString fUnfoldOutPath;

    bool verbose;

    // Conditions for unfolding
    TUnfold::ERegMode regMode;
    double nominal_bias;

    TString output_baseDir;
    TString unfold_name;
    TString channel_name;
    TString var;
    int year;

    bool doInputStatUnc;
    bool doRMStatUnc;

public:

    // Constructor
    ISRUnfold(TString unfold_name_, TString channel, int year_ = 2016, int regMode_ = 0, bool doInputStatUnc_ = true, bool doRMStatUnc_ = false, bool verbose_ = false, TString var_ = "Mass", TString output_baseDir_ = "")
    {
        cout << "ISRUnfold set! " << var_ << endl;

        output_baseDir = output_baseDir_;
        unfold_name = unfold_name_;
        channel_name = channel;
        year = year_;
        var = var_;

        nominal_bias = 1.;

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

        verbose = verbose_;

        // Make output root files
        TString yearStr;
        yearStr.Form("%d", (int)year);

        //output_baseDir = "output/" + yearStr + "/" + channel_name + "/";

        fUnfoldOutPath = output_baseDir+unfold_name+"_"+channel_name+"_"+yearStr+"_"+var+".root";
        fUnfoldOut = new TFile(fUnfoldOutPath, "RECREATE");

        fUnfoldOut->mkdir("unfolded");
        fUnfoldOut->mkdir("acceptance");

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

    const TVectorD& checkMatrixCond();
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
    void doAcceptCorr(TString filePath, TString binDef, bool isAccept = false);

    // Get histograms
    TH1* getUnfoldedHists(TString outHistName = "", TString steering = "", bool useAxis = true);
    TH1* getRawHist(TString filePath, TString dirName, TString histName, TString outHistName, TString steering, bool useAxis);

    TH1* getUnfInput(TString steering, bool useAxis, int massBin);

    // Helper functions
    void doNorm(TH1* hist, bool norm = true);
};

#endif
