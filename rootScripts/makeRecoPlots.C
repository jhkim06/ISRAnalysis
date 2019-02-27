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
struct recoTH1info {
  std::vector<TH1*> hists;
  std::vector<std::map<TString,Double_t>> sysNamesWeights;
  std::vector<TString> sysNames;
  Bool_t isSig;
  
  recoTH1info(std::vector<TH1*> ihists, std::vector<TString> sysNames_, Bool_t isSig_):
    hists(std::move(ihists)), sysNames(std::move(sysNames_)), isSig(std::move(isSig_)) {}
};

TUnfoldBinningV17* binning(){
 const int nmassbin_fine=5;
 double massbin_fine[nmassbin_fine+1]={40,60,80,100,200,350};
 const int nmassbin_wide=5;
 double massbin_wide[nmassbin_wide+1]={40,60,80,100,200,350};

 // pt bins for reco
 const int nptbin_fine=50;
 double ptbin_fine[nptbin_fine+1];

 for(int i = 0; i < nptbin_fine + 1; i++){
    ptbin_fine[i] = i*2;
 }

 // pt bins for gen
 const int nptbin_wide=20;
 double ptbin_wide[nptbin_wide+1];

 for(int i = 0; i < nptbin_wide + 1; i++){
     ptbin_wide[i] = i*5;
 }

 TUnfoldBinningV17 *binning_Rec=new TUnfoldBinningV17("Rec");
 binning_Rec->AddAxis("pt",nptbin_fine,ptbin_fine,false,true);
 binning_Rec->AddAxis("mass",nmassbin_fine,massbin_fine,true,true);

 return binning_Rec;
}

// create hisgotram using TUnfoldBinningV17
TH1* histogram(TString name){
 return binning()->CreateHistogram("h"+ name+ "Rec");
}

void recoHists(TFile *filein, TFile *fileout1, TFile *fileout2, const recoTH1info &recoHist, const recoTH1info &recoHist_){ // TODO add list of systematics

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

 TUnfoldBinningV17 *bin = binning();

 // TODO based on the info in recoTH1info make map for systematics
 for(int i=0;i<nentries;i++){
 //for(int i=0;i<10000;i++){
   if(i%10000000==0) cout<<i<<endl;
   treco->GetEntry(i);
    if(ispassRec && isBveto && ptRec->at(2) < 100){
       if(!recoHist.isSig){ // for data and MC except DY
          fileout1->cd();
          recoHist.hists.at(0)->Fill(bin->GetGlobalBinNumber(ptRec->at(2),mRec->at(2)), weightRec);
       }
       else{
            if(DYtautau){
               fileout2->cd();
               recoHist_.hists.at(0)->Fill(bin->GetGlobalBinNumber(ptRec->at(2),mRec->at(2)), weightRec);
            }
            else{
               fileout1->cd();
               recoHist.hists.at(0)->Fill(bin->GetGlobalBinNumber(ptRec->at(2),mRec->at(2)), weightRec);
            }
       }
    }
 }// event loop

 delete bin;
 delete treco;
 //fileout->cd();
 
 // seems histograms automatically written 
 // recoHist.hists.at(0)->Write();
 // delete outputFile;
}

