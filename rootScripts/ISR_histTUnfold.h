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

const int pt_unfold = 1;
const int mass_unfold = 2;

class histTUnfold {

private:

	// binning definition for pt
	TUnfoldBinning* ptBinningRec;
	TUnfoldBinning* ptBinningGen;

	// binning definition for mass
	TUnfoldBinning* massBinningRec;
	TUnfoldBinning* massBinningGen;

  	std::map<TString, TH1*> histMaps;
  	std::map<TString, TH2*> hist2DMaps;

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
        vector<Double_t> *etaRec;
        vector<Double_t> *mRec;
        vector<Double_t> *l1PreFire;
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
  	Bool_t isInc, isSig;

  	// binning
	void SetPtBinningRec();
	void SetPtBinningGen();
	void SetMassBinningRec();
	void SetMassBinningGen();

	inline TUnfoldBinning* GetPtBinningRec(){ return ptBinningRec;}
	inline TUnfoldBinning* GetMassBinningRec(){ return massBinningRec;}
	inline TUnfoldBinning* GetPtBinningGen(){ return ptBinningGen;}
	inline TUnfoldBinning* GetMassBinningGen(){ return massBinningGen;}

	// set histograms: get histogram names from python script and create histograms
	void CreateHistMap(const int which_unfold, TString hname, TString postfix = ""); // postfix for DY to tautau events
	void CreateHist2DMap(const int which_unfold, TString hname);
  	
	// save histograms into root files  
	void saveRecoHists(TFile *filein, TFile *fileout1, TString channel); 
        void saveSigHists(TFile *filein, TFile *fileout1, TFile *fileout2, TString channel, Double_t temp_kfactor);

  	void FillMigration2DM(const int which_unfold, bool selected, TString hname, Double_t recoPt, Double_t RecoMass, Double_t truthPt, Double_t truthMass, Double_t wreco, Double_t wgen, Double_t corr = 1.);
  	void FillHistogram(const int which_unfold, TString hname, Double_t recoPt, Double_t recoMass, Double_t wreco);

  	histTUnfold(std::map<TString, TH1*> histMaps_, std::map<TString, TH2*> hist2DMaps_, Int_t isInc_):
    		histMaps(std::move(histMaps_)), hist2DMaps(std::move(hist2DMaps_)), isInc(std::move(isInc_)) {}

  	histTUnfold(std::map<TString, TH1*> histMaps_):
    		histMaps(std::move(histMaps_)) {}

	histTUnfold() {}

	~histTUnfold(){

        }
};

#endif
