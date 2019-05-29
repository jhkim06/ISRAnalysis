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

// FIXME change name since this class will be used for all input root files
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

     	// weightMaps
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
	void CreateHistMap(const int which_unfold, TString hname, TString postfix = "");
	void CreateHist2DMap(const int which_unfold, TString hname);
 
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
