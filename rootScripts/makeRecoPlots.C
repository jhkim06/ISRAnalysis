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
#include "tree.h"

// FIXME may be it's better to set this in the unfoldingISR to select do only norminal unfolding or include systematic unfolding
// TODO use map for systematics ex) std::vector<std::map<TString,Double_t>> string for systematic name double for weight
struct recoHistsinfo {
  std::vector<TH1*> hists;
  std::vector<std::map<TString,Double_t>> sysNamesWeights; // will be initialze in the function makeRecoHists
  std::vector<TString> sysNames;
  
  recoHistsinfo(std::vector<TH1*> ihists, std::vector<TString> sysNames_):
    hists(std::move(ihists)), sysNames(std::move(sysNames_)) {}
};

struct sigHistsinfo {
  std::vector<TH1*> hists;
  std::vector<TH2*> matrixs;
  std::vector<std::map<TString,Double_t>> sysNamesWeights;
  std::vector<TString> sysNames;
  Bool_t isInc;

  sigHistsinfo(std::vector<TH1*> ihists, std::vector<TH2*> imatrixs, std::vector<TString> sysNames_, Int_t isInc_):
    hists(std::move(ihists)), matrixs(std::move(imatrixs)), sysNames(std::move(sysNames_)), isInc(std::move(isInc_)) {}
};

TUnfoldBinningV17* binning_rec(){

 // FIXME save binning information in other place and read from there
 const int nmassbin_fine=5;
 double massbin_fine[nmassbin_fine+1]={40,60,80,100,200,350};

 // pt bins for reco
 const int nptbin_fine=50;
 double ptbin_fine[nptbin_fine+1];

 for(int i = 0; i < nptbin_fine + 1; i++){
    ptbin_fine[i] = i*2;
 }

 TUnfoldBinningV17 *binning_Rec=new TUnfoldBinningV17("Rec");
 binning_Rec->AddAxis("pt",nptbin_fine,ptbin_fine,false,true);
 binning_Rec->AddAxis("mass",nmassbin_fine,massbin_fine,true,true);

 return binning_Rec;
}

TUnfoldBinningV17* binning_gen(){

 // FIXME save binning information in other place and read from there
 const int nmassbin_wide=5;
 double massbin_wide[nmassbin_wide+1]={40,60,80,100,200,350};

 // pt bins for gen
 const int nptbin_wide=20;
 double ptbin_wide[nptbin_wide+1];

 for(int i = 0; i < nptbin_wide + 1; i++){
     ptbin_wide[i] = i*5;
 }

 TUnfoldBinningV17 *binning_Gen=new TUnfoldBinningV17("Rec");
 binning_Gen->AddAxis("pt",nptbin_wide,ptbin_wide,false,true);
 binning_Gen->AddAxis("mass",nmassbin_wide,massbin_wide,true,true);

 return binning_Gen;
}

// create hisgotram using TUnfoldBinningV17
TH1* histogram(TString name){
 return binning_rec()->CreateHistogram("hRec"+name);
}

TH2* matrix(TString name){
 return TUnfoldBinningV17::CreateHistogramOfMigrations(binning_gen(),binning_rec(),"hmcGenRec" + name);
}

void recoHists(TFile *filein, TFile *fileout1, const recoHistsinfo &recoHist){ // TODO add list of systematics

 gROOT->SetBatch();
 TH1::SetDefaultSumw2();

 TTree *treco=(TTree *)filein->Get("tree");
 treco->SetBranchAddress("qLep",&qLep);
 treco->SetBranchAddress("ptRec",&ptRec);
 treco->SetBranchAddress("etaRec",&etaRec);
 treco->SetBranchAddress("phiRec",&phiRec);
 treco->SetBranchAddress("mRec",&mRec);
 treco->SetBranchAddress("ispassRec",&ispassRec);
 treco->SetBranchAddress("isBveto",&isBveto);
 treco->SetBranchAddress("weightRec",&weightRec);
 treco->SetBranchAddress("nVtx",&nVtx);
 treco->SetBranchAddress("DYtautau",&DYtautau);
 nentries=treco->GetEntries();

 TUnfoldBinningV17 *bin = binning_rec();

 // TODO based on the info in recoHistsinfo make map for systematics
 for(int i=0;i<nentries;i++){
 //for(int i=0;i<10000;i++){
   if(i%10000000==0) cout<<i<<endl;
   treco->GetEntry(i);
    if(ispassRec && isBveto && ptRec->at(2) < 100){
          fileout1->cd();
          recoHist.hists.at(0)->Fill(bin->GetGlobalBinNumber(ptRec->at(2),mRec->at(2)), weightRec);
    }
 }// event loop

 delete bin;
 delete treco;
 //fileout->cd();
 
 // seems histograms automatically written 
 // recoHist.hists.at(0)->Write();
 // delete outputFile;
}

void sigHists(TFile *filein, TFile *fileout1, TFile *fileout2, const sigHistsinfo &sigHist, const recoHistsinfo &recoHist){ // TODO add list of systematics

 gROOT->SetBatch();
 TH1::SetDefaultSumw2();

 TTree *treco=(TTree *)filein->Get("tree");
 treco->SetBranchAddress("qLep",&qLep);
 treco->SetBranchAddress("ptRec",&ptRec);
 treco->SetBranchAddress("etaRec",&etaRec);
 treco->SetBranchAddress("phiRec",&phiRec);
 treco->SetBranchAddress("mRec",&mRec);
 treco->SetBranchAddress("ispassRec",&ispassRec);
 treco->SetBranchAddress("isBveto",&isBveto);
 treco->SetBranchAddress("weightRec",&weightRec);
 treco->SetBranchAddress("nVtx",&nVtx);
 treco->SetBranchAddress("DYtautau",&DYtautau);
 nentries=treco->GetEntries();

 TUnfoldBinningV17 *bin_rec = binning_rec();
 TUnfoldBinningV17 *bin_gen = binning_gen();

 // TODO based on the info in recoHistsinfo make map for systematics
 for(int i=0;i<nentries;i++){
 //for(int i=0;i<10000;i++){
   if(i%10000000==0) cout<<i<<endl;
   treco->GetEntry(i);

    if(!sigHist.isInc){
       if(ispassRec && isBveto && ptRec->at(2) < 100){
             fileout1->cd();
             sigHist.hists.at(0)->Fill(bin_rec->GetGlobalBinNumber(ptRec->at(2),mRec->at(2)), weightRec); //
       }
    }
    else{ // for inclusive DY MC

         if(ispassRec && isBveto && ptRec->at(2) < 100){
            if(DYtautau){
               fileout2->cd();
               recoHist.hists.at(0)->Fill(bin_rec->GetGlobalBinNumber(ptRec->at(2),mRec->at(2)), weightRec);
            }
            else{
               fileout1->cd();
               sigHist.hists.at(0)->Fill(bin_rec->GetGlobalBinNumber(ptRec->at(2),mRec->at(2)), weightRec);
            }
         }
    }

 }// event loop

 delete bin_rec;
 delete bin_gen;
 delete treco;
 //fileout->cd();

 // seems histograms automatically written 
 // recoHist.hists.at(0)->Write();
 // delete outputFile;
}

