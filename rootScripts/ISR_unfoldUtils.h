#ifndef UNFOLDUTILS_H
#define UNFOLDUTILS_H

#include <iostream>
#include <iomanip> 
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
    std::map<TString, TH1*> sysAbsPtHist_detector;

    std::map<int, std::map<TString, TH1*>> sysRelPtHist_detectorData;
    std::map<int, std::map<TString, TH1*>> sysRelPtHist_detectorMC;
    std::map<TString, TH1*> sysRelMassHist_detectorData;
    std::map<TString, TH1*> sysRelMassHist_detectorMC;

    std::map<int, std::map<TString, TH1*>> sysRelPtHist_unfoldedData;
    std::map<int, std::map<TString, TH1*>> sysRelPtHist_unfoldedMC;
    std::map<TString, TH1*> sysRelMassHist_unfoldedData;
    std::map<TString, TH1*> sysRelMassHist_unfoldedMC;

    std::map<int, std::map<TString, TH1*>> sysRelPtHist_unfoldedAcceptData;
    std::map<int, std::map<TString, TH1*>> sysRelPtHist_unfoldedAcceptMC;
    std::map<TString, TH1*> sysRelMassHist_unfoldedAcceptData;
    std::map<TString, TH1*> sysRelMassHist_unfoldedAcceptMC;
    std::map<int, std::map<TString, TH1*>> sysAbsPtHist_unfoldedAcceptMC;
    std::map<TString, TH1*> sysAbsMassHist_unfoldedAcceptMC;


    // Unfolding
    // For nominal results
    TUnfoldDensityV17* nomPtUnfold;
    TUnfoldDensityV17* nomMassUnfold;

    TH2* hPtResponseM;
    TH2* hMassResponseM;

    // For statistical uncertainty
    std::vector<TUnfoldDensityV17*> statPtUnfold;
    std::vector<TUnfoldDensityV17*> statMassUnfold;

    // Histogram to save statistical variations for each mass bin
    std::vector<TH1*> meanPtStatVariation;
    std::vector<TH1*> meanMassStatVariation;

    std::vector<TH1*> meanPtPDFVariation;
    std::vector<TH1*> meanMassPDFVariation;

    std::vector<TH1*> meanPtPDFVariation_Accept;
    std::vector<TH1*> meanMassPDFVariation_Accept;

    // For systematic uncertainty
    std::map<TString, std::map<TString, TUnfoldDensityV17*>> sysPtUnfold;
    std::map<TString, std::map<TString, TUnfoldDensityV17*>> sysMassUnfold;
    
    /*------------------------------------------------------------------------------------- Results: mean values ---------------------------------------------------------------------------------------*/ 
    // Folded level: from TUnfoldDensity (mean from the binned histogram) 
    vector<Double_t> meanMass_data_folded, meanMassStatErr_data_folded, meanMassSysErr_data_folded, meanMassTotErr_data_folded;
    vector<Double_t> meanPt_data_folded,   meanPtStatErr_data_folded,   meanPtSysErr_data_folded, meanPtTotErr_data_folded;
    // Mean from raw histograms 
    vector<Double_t> unbinnedMeanMass_data_folded, unbinnedMeanMassStatErr_data_folded, unbinnedMeanMassSysErr_data_folded, unbinnedMeanMassTotErr_data_folded;
    vector<Double_t> unbinnedMeanPt_data_folded,   unbinnedMeanPtStatErr_data_folded,   unbinnedMeanPtSysErr_data_folded, unbinnedMeanPtTotErr_data_folded;
    // Unfolded results
    vector<Double_t> meanMass_data_unfolded, meanMassStatErr_data_unfolded, meanMassSysErr_data_unfolded, meanMassTotErr_data_unfolded;
    vector<Double_t> meanPt_data_unfolded,   meanPtStatErr_data_unfolded,   meanPtSysErr_data_unfolded, meanPtTotErr_data_unfolded;

    // Nominal mean mass and pt for MC
    vector<Double_t> meanMass_mc_folded, meanMassStatErr_mc_folded, meanMassSysErr_mc_folded, meanMassTotErr_mc_folded;
    vector<Double_t> meanPt_mc_folded, meanPtStatErr_mc_folded, meanPtSysErr_mc_folded, meanPtTotErr_mc_folded;

    vector<Double_t> meanMass_mc_unfolded, meanMassStatErr_mc_unfolded, meanMassSysErr_mc_unfolded, meanMassTotErr_mc_unfolded;
    vector<Double_t> meanPt_mc_unfolded, meanPtStatErr_mc_unfolded, meanPtSysErr_mc_unfolded, meanPtTotErr_mc_unfolded;

    std::map<TString, std::map<TString, vector<double>>> meanMass_data_folded_sysVariation;
    std::map<TString, std::map<TString, vector<double>>> meanPt_data_folded_sysVariation;

    std::map<TString, std::map<TString, vector<double>>> meanMass_data_unfolded_sysVariation;
    std::map<TString, std::map<TString, vector<double>>> meanPt_data_unfolded_sysVariation;

    std::map<TString, vector<double>> meanMass_data_folded_systematic;
    std::map<TString, vector<double>> meanPt_data_folded_systematic;

    std::map<TString, vector<double>> meanMass_data_folded_rel_systematic;
    std::map<TString, vector<double>> meanPt_data_folded_rel_systematic;
    
    vector<Double_t> meanMassStatRelErr_data_unfolded, meanPtStatRelErr_data_unfolded;
    /*----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
   
    // Acceptance correction
    TH1* hFullPhaseMassData; // data after acceptance correction 
    TH1* hFullPhasePtData;  

    std::map<TString, std::map<TString, TH1*>> hSysFullPhaseMassData;
    std::map<TString, std::map<TString, TH1*>> hSysFullPhasePtData;
    std::map<TString, std::map<TString, TH1*>> hSysFullPhaseMassMC;
    std::map<TString, std::map<TString, TH1*>> hSysFullPhasePtMC;
    
    TH1* hAcceptanceMass;  
    TH1* hAcceptancePt;  

    vector<Double_t> meanMass_data_acc_corrected, meanMassStatErr_data_acc_corrected, meanMassSysErr_data_acc_corrected, meanMassTotErr_data_acc_corrected;
    vector<Double_t> meanPt_data_acc_corrected,   meanPtStatErr_data_acc_corrected,   meanPtSysErr_data_acc_corrected, meanPtTotErr_data_acc_corrected;

    std::map<TString, std::map<TString, vector<double>>> meanMass_data_accept_sysVariation;
    std::map<TString, std::map<TString, vector<double>>> meanPt_data_accept_sysVariation;

    std::map<TString, vector<double>> meanMass_data_accept_systematic;
    std::map<TString, vector<double>> meanPt_data_accept_systematic;

    std::map<TString, vector<double>> meanMass_data_accept_rel_systematic;
    std::map<TString, vector<double>> meanPt_data_accept_rel_systematic;

    bool verbose;

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
    ISRUnfold(TString channel, int year_ = 2016, int regMode_ = 0, bool makeStatUnfold_ = true, bool verbose_ = false)
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
        bkgColors["Fake"] = kViolet+6;
        bkgColors["WJets"] = kViolet+1;
        bkgColors["EWK"] = kOrange+1;
        bkgColors["Top"] = kAzure+1;

        verbose = verbose_;

    }
    // Destructor
    ~ISRUnfold(){}

    void drawResponseM(TString var, TString sysName = "", TString sysPostfix = "", bool isDetector = true);

    const TVectorD& checkMatrixCond(TString var = "Mass");
    double getSmearedChi2(TString var, TString filePath, TString dirName, TString steering, bool useAxis, bool divBinWidth = false);
    double getUnfoldedChi2(TString var, TString steering, bool useAxis, bool divBinWidth = false);

    void setOutputBaseDir(TString outPath);
    void setBias(double bias);

    // Set nominal TUnfoldDensity 
    void setNomResMatrix(TString var, TString filepath, TString dirName, TString histName, TString binDef = "");

    void setFromPrevUnfResult(ISRUnfold* unfold, bool useAccept = false);

    // Set input histogram
    void setUnfInput(TString var, TString varPostfix = "", TString filepath = "", TString dirName ="", TString histName = "", bool isSys = false, TString sysName = "", TString sysPostfix = "", bool isFSR = false);
    void setUnfInput(ISRUnfold* unfold, TString var, bool isSys = false, TString sysName = "", TString sysPostfix = "", bool useAccept = false);

    // Set background histograms
    void subBkgs(TString filepath, std::pair<TString, TString>& bkgInfo, 
                 bool isSys = false, TString binDef = "", TString dirName = "", TString sysName = "", TString sysPostfix = "", bool isFSR = false);

    // Set systematic TUnfoldDensity
    void setSysTUnfoldDensity(TString var, TString filepath, TString dirName, TString histName, TString sysName, TString sysPostfix, TString binDef);

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
    inline std::map<TString, std::vector<TString>>& getSystematicMap()
    {
        return sysMap;
    }
    // Draw folded distribution(before unfolding) using histograms saved in TUnfoldDensity
    TCanvas* drawFoldedHists(TString var, TString filePath, TString dirName, TString steering, bool useAxis, TString sysName = "", TString outName = "", int nthMassBin = 0, bool divBinWidth = false, TString sysFilePath = "");
    TCanvas* drawUnfoldedHists(TString var, TString steering, bool useAxis, TString sysName = "", TString outName = "", int nthMassBin = 0, bool divBinWidth = false, bool isType3Closure = false);
    TH1* getDetectorSystematicBand(TString var, TString filePath, TString dirName, TString steering, bool useAxis, TString sysName, TH1* hData, TH1* hDY, TH1* hMCtotal, TH1* hRatio, bool divBinWidth, bool isRatio, bool forMC, TString sysFilePath = "", int nthMassBin = 0);
    TH1* getUnfoldedSystematicBand(TString var, TString steering, bool useAxis, TString sysName, TH1* hData, TH1* hDY, TH1* hRatio, bool divBinWidth, bool isRatio, bool forMC, int nthMassBin = 0);
    TH1* getUnfAcceptSystematicBand(TString var, TString steering, bool useAxis, TString sysName, TH1* hData, TH1* hDY, TH1* hRatio, bool divBinWidth, bool isRatio, bool forMC, int nthMassBin = 0);
    void setTHStack(TString var, TString filePath, TString dirName, THStack& hs, TH1& hMCtotal, TLegend& leg, TString steering, bool useAxis, TString sysName = "", bool divBinWidth = false);
    void setXaxisTitle(TH1* hist, TString var, bool useAxis, TString title = "");
    void writeCutInfo(TPad* pad, TString var, int nthMassBin = 0);
    double getBinnedMean(TH1* hist);
    void divideByBinWidth(TH1* hist, bool norm = false);

    // Do unfold 
    void doISRUnfold( bool doSys = false, bool doReg = false);
    void doStatUnfold(); 

    void setStatError();
    void setSysError();
    void setSysError_Accept();
    void setTotSysError();
    void setTotSysError_Accept();

    void doAcceptCorr(TString filePath, TString binDef, bool doSys = false, TString outName = "");
    void drawAcceptance(TString var, TH1* hMC, TString outName);
    TCanvas* drawAcceptCorrHists(TString var, TString filePath, TString binDef, TString steering, bool useAxis, TString sysName, TString outName, int nthMassBin, bool divBinWidth); 
    void drawCorrelation(TString var, TString steering, bool useAis, TString outName = "");
    void drawComparisonPlot(TString var, TString plotName, TString topYaxisName, TString bottomYaxisName, TString bottomXaxisName, TH1* h1, TH1* h2, TH1* hratio, TString outName, int nthMassBin = 0);

    // Get histograms
    TH1* getUnfoldedHists(TString var, TString outHistName = "", TString steering = "", bool useAxis = true, bool divBinWidth = false);
    TH1* getGenMCHist(TString var, TString steering, bool useAxis = true, int massBin = 0, bool binWidth = false);
    TH1* getDetHists(TString var, TString outHistName = "", TString steering = "", bool useAxis = true);
    TH1* getRawHist(TString var, TString filePath, TString dirName, TString histName, TString outHistName, TString steering, bool useAxis, bool divBinWidth = false);

    TH1* getUnfInput(TString var, TString steering, bool useAxis, int massBin, bool binWidth);

    // Helper functions
    void doNorm(TH1* hist, bool norm = true); 
    void drawtext(TGraph* g);

    int setMeanPt(TString filePath = "", TString dirName = "");
    int setMeanMass(TString filePath = "", TString dirName = "");
    void setMeanPt_Accept();
    void setMeanMass_Accept();
    int setSysMeanPt();
    int setSysMeanMass();
    void setSysMeanPt_Accept();
    void setSysMeanMass_Accept();

    void fillPtStatVariationHist(int istat);
    void fillMassStatVariationHist(int istat);

    void fillPtPDFVariationHist(int istat);
    void fillPtPDFVariationHist_Accept(int istat);
    void fillMassPDFVariationHist(int istat);
    void fillMassPDFVariationHist_Accept(int istat);
   
    // Get mean values 
    double getDetMeanPt(int ibin);
    double getDetMeanMass(int ibin);
    double getDetMeanPtError(int ibin);
    double getDetMeanMassError(int ibin);
    double getMCDetMeanPt(int ibin);
    double getMCDetMeanMass(int ibin);
    double getMCDetMeanPtError(int ibin);
    double getMCDetMeanMassError(int ibin);
    double getUnfMeanPt(int ibin);
    double getUnfMeanMass(int ibin);
    double getUnfMeanPtError(int ibin);
    double getUnfMeanMassError(int ibin);
    double getUnfMeanPtSysError(int ibin);
    double getUnfMeanMassSysError(int ibin);
    double getAccMeanPt(int ibin);
    double getAccMeanPtError(int ibin);
    double getAccMeanMass(int ibin);
    double getAccMeanMassError(int ibin);
    double getAccMeanPtSysError(int ibin);
    double getAccMeanMassSysError(int ibin);
    double getAccMeanPtTotError(int ibin);
    double getAccMeanMassTotError(int ibin);

    double getMCGenMeanMass(int ibin);
    double getMCGenMeanPt(int ibin);

    void printMeanValues(bool printSys);
    void printMeanValues_Accept(bool printSys = false);
    void drawStatVariation(bool isPt = true, int massBin = 0);
    void drawPDFVariation(bool isPt = true, int massBin = 0);
    void drawSysVariation(TString sysName, TString var, int massBin);
    void drawSystematics(TString var);
    void drawSystematics_Acceptance(TString var);
};

#endif
