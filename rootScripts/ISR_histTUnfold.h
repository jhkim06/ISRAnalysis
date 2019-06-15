#ifndef HISTUNFOLD_H
#define HISTUNFOLD_H

#include <iostream>
#include <map>
#include <cmath>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TDOMParser.h>
#include <TXMLDocument.h>
#include "TUnfoldBinningXML.h"
#include "TUnfoldBinning.h"

#include "TLorentzVector.h"
using namespace std;

enum ptOrMass{
  PT = 1, MASS
};

const int nmassBins_fine_muon = 58; 
const int nmassBins_wide_muon = 29; 
const Double_t massBins_fine_muon[nmassBins_fine_muon+1] =     {40,42.5,45,47.5,50,52.5,55,57.5,60,62.5,65,67.5,70,72.5,75,77.5,80,82.5,85,87.5,90,92.5,95,97.5,100,102.5,105,107.5,110,112.5,115,117.5,120,123,126,129.5,133,137,141,145.5,150,155,160,165.5,171,178,185,192.5,200,209,218,229,240,254,268,284,300,325,350};
const Double_t massBins_wide_muon[nmassBins_wide_muon+1] =     {40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,126,133,141,150,160,171,185,200,218,240,268,300,350};

const int nmassBins_fine_electron = 54; 
const int nmassBins_wide_electron = 27; 
const Double_t massBins_fine_electron[nmassBins_fine_electron+1] = {50,52.5,55,57.5,60,62.5,65,67.5,70,72.5,75,77.5,80,82.5,85,87.5,90,92.5,95,97.5,100,102.5,105,107.5,110,112.5,115,117.5,120,123,126,129.5,133,137,141,145.5,150,155,160,165.5,171,178,185,192.5,200,209,218,229,240,254,268,284,300,325,350};
const Double_t massBins_wide_electron[nmassBins_wide_electron+1] = {50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,126,133,141,150,160,171,185,200,218,240,268,300,350}; 

const int nmassBins_forPt = 5;
const Double_t massBins_forPt_muon[nmassBins_forPt+1] =     {40,60,80,100,200,350};
const Double_t massBins_forPt_electron[nmassBins_forPt+1] = {50,65,80,100,200,350};

class histTUnfold {

private:

	// binning definition for pt
	TUnfoldBinning* ptBinningRec;
	TUnfoldBinning* ptBinningGen;

	// binning definition for mass
	TUnfoldBinning* massBinningRec;
	TUnfoldBinning* massBinningGen;

  	// detector level histograms
  	std::map<TString, TH1*> histMaps; 
	// response matrix
  	std::map<TString, TH2*> hist2DMaps;
	// map for systematic sources: systematic name, the number of systematic variations
  	std::map<TString, int > sysMaps;

	TTree *tree;

        // variables in the input tree
        Double_t weightGen,weightRec;

        int ispassRec, DYtautau, isBveto;
        int isfiducialPostFSR, isfiducialPreFSR;
        int nentries;

        vector<Double_t> *ptPreFSR;
        vector<Double_t> *mPreFSR;
        vector<Double_t> *ptPostFSR;
        vector<Double_t> *mPostFSR;
        vector<Double_t> *ptRec;
        vector<Double_t> *mRec;
        vector<Double_t> *ptRec_momentumUp;
        vector<Double_t> *mRec_momentumUp;
        vector<Double_t> *ptRec_momentumDown;
        vector<Double_t> *mRec_momentumDown;
        vector<Double_t> *ptRec_momentumResUp;
        vector<Double_t> *mRec_momentumResUp;
        vector<Double_t> *ptRec_momentumResDown;
        vector<Double_t> *mRec_momentumResDown;
        vector<Double_t> *AlphaS;
        vector<Double_t> *Scale;
        vector<Double_t> *PDFerror;
        Int_t nVtx;

        vector<TLorentzVector> *particleFSR;
        vector<TLorentzVector> *anparticleFSR;
        TLorentzVector *particlePostFSR;
        TLorentzVector *anparticlePostFSR;

        Double_t        PUweight;
        Double_t        PUweight_Up;
        Double_t        PUweight_Dn;
        Double_t        trgSF;
        Double_t        trgSF_Up;
        Double_t        trgSF_Dn;
        Double_t        recoSF;
        Double_t        recoSF_Up;
        Double_t        recoSF_Dn;
        Double_t        IdSF;
        Double_t        IdSF_Up;
        Double_t        IdSF_Dn;
        Double_t        IsoSF;
        Double_t        IsoSF_Up;
        Double_t        IsoSF_Dn;
        Double_t        L1Prefire;
        Double_t        L1Prefire_Up;
        Double_t        L1Prefire_Dn;

        TBranch        *b_particleFSR;   //!
        TBranch        *b_anparticleFSR;   //!
        TBranch        *b_particlePostFSR;   //!
        TBranch        *b_anparticlePostFSR;   //!

        Int_t IsElEl, IsMuMu;
        Double_t ZPtCor;
        Double_t bTagReweight;
        Int_t isdielectron, isdimuon;


public:

  	std::map<TString, Double_t> sysNamesWeights;
  	Bool_t isInc, isSig, isAlt;

  	// binning
	void SetPtBinningRec(TString channel);
	void SetPtBinningGen(TString channel);
	void SetMassBinningRec(TString channel);
	void SetMassBinningGen(TString channel);

	inline TUnfoldBinning* GetPtBinningRec(){ return ptBinningRec;}
	inline TUnfoldBinning* GetMassBinningRec(){ return massBinningRec;}
	inline TUnfoldBinning* GetPtBinningGen(){ return ptBinningGen;}
	inline TUnfoldBinning* GetMassBinningGen(){ return massBinningGen;}

	// set histograms: get histogram names from python script and create histograms
	void CreateHistMap(int which_unfold, TString hname, TString postfix = "", bool isRec = true); // postfix for DY to tautau events
	void CreateHist2DMap(int which_unfold, TString hname);
	void SetsysMap(TString sysName, int nVariations);

	Double_t GetSysWeights(TString sysName, bool isReco, int nthSys);
  	
	// save histograms into root files  
	void saveRecoHists(TFile *filein, TFile *fileout1, TString channel); 
        void saveSigHists(TFile *filein, TFile *fileout1, TFile *fileout2, TString channel, Double_t temp_kfactor);

  	void FillMigration2DM(ptOrMass which_unfold, bool selected, TString hname, Double_t recoPt, Double_t RecoMass, Double_t truthPt, Double_t truthMass, Double_t wreco, Double_t wgen, Double_t corr = 1.);
  	void FillHistogram(ptOrMass which_unfold, TString hname, Double_t recoPt, Double_t recoMass, Double_t wreco, bool isRec = true);

  	histTUnfold(std::map<TString, TH1*> histMaps_, std::map<TString, TH2*> hist2DMaps_, Int_t isInc_, Int_t isAlt_):
    		histMaps(std::move(histMaps_)), hist2DMaps(std::move(hist2DMaps_)), isInc(std::move(isInc_)), isAlt(std::move(isAlt_)) {}

  	histTUnfold(std::map<TString, TH1*> histMaps_):
    		histMaps(std::move(histMaps_)) {}

	histTUnfold() {}

	~histTUnfold(){

        }
};

#endif
