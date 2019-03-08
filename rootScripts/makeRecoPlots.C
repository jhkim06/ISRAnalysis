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

#include "makeRecoPlots.h"

// FIXME may be it's better to set this in the unfoldingISR to select do only norminal unfolding or include systematic unfolding
// TODO use map for systematics ex) std::vector<std::map<TString,Double_t>> string for systematic name double for weight
struct recoHistsinfo {
  std::vector<TH1*> ptHists;
  std::vector<TH1*> massHists;
  std::vector<std::map<TString,Double_t>> sysNamesWeights; // will be initialze in the function makeRecoHists
  std::vector<TString> sysNames;
  
  recoHistsinfo(std::vector<TH1*> iptHists, std::vector<TH1*> imassHists, std::vector<TString> sysNames_):
    ptHists(std::move(iptHists)), massHists(std::move(imassHists)), sysNames(std::move(sysNames_)) {}
};

struct sigHistsinfo {
  std::vector<TH1*> ptHists;
  std::vector<TH1*> massHists;
  std::vector<TH2*> ptMatrixs;
  std::vector<TH2*> massMatrixs;
  std::vector<std::map<TString,Double_t>> sysNamesWeights;
  std::vector<TString> sysNames;
  Bool_t isInc;

  sigHistsinfo(std::vector<TH1*> iptHists, std::vector<TH1*> imassHists, std::vector<TH2*> iptMatrixs, std::vector<TH2*> imassMatrixs, std::vector<TString> sysNames_, Int_t isInc_):
    ptHists(std::move(iptHists)), massHists(std::move(imassHists)), ptMatrixs(std::move(iptMatrixs)), massMatrixs(std::move(imassMatrixs)), sysNames(std::move(sysNames_)), isInc(std::move(isInc_)) {}
};

TUnfoldBinningV17* ptBinning_rec(){

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

TUnfoldBinningV17* ptBinning_gen(){

 // FIXME save binning information in other place and read from there
 const int nmassbin_wide=5;
 double massbin_wide[nmassbin_wide+1]={40,60,80,100,200,350};

 // pt bins for gen
 const int nptbin_wide=20;
 double ptbin_wide[nptbin_wide+1];

 for(int i = 0; i < nptbin_wide + 1; i++){
     ptbin_wide[i] = i*5;
 }

 TUnfoldBinningV17 *binning_Gen=new TUnfoldBinningV17("Gen");
 binning_Gen->AddAxis("pt",nptbin_wide,ptbin_wide,false,true);
 binning_Gen->AddAxis("mass",nmassbin_wide,massbin_wide,true,true);

 return binning_Gen;
}

TUnfoldBinningV17* massBinning_rec(){

 // FIXME save binning information in other place and read from there
 const int nbin_fine=58;
 double bin_fine[nbin_fine+1]={40,42.5,45,47.5,50,52.5,55,57.5,60,62.5,65,67.5,70,72.5,75,77.5,80,82.5,85,87.5,90,92.5,95,97.5,100,102.5,105,107.5,110,112.5,115,117.5,120,123,126,129.5,133,137,141,145.5,150,155,160,165.5,171,178,185,192.5,200,209,218,229,240,254,268,284,300,325,350};

 TUnfoldBinningV17 *binning_Rec=new TUnfoldBinningV17("Rec");
 binning_Rec->AddAxis("reco mass",nbin_fine,bin_fine,false,false);

 return binning_Rec;
}

TUnfoldBinningV17* massBinning_gen(){

 // FIXME save binning information in other place and read from there
 const int nbin_wide=29;
 double bin_wide[nbin_wide+1]={40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,126,133,141,150,160,171,185,200,218,240,268,300,350};

 TUnfoldBinningV17 *binning_Gen=new TUnfoldBinningV17("Gen");
 binning_Gen->AddAxis("gen mass",nbin_wide,bin_wide,true,true);

 return binning_Gen;
}

// create hisgotram using TUnfoldBinningV17
TH1* ptHistogram(TString name){
 return ptBinning_rec()->CreateHistogram("hPtRec"+name);
}

TH1* massHistogram(TString name){
 return massBinning_rec()->CreateHistogram("hMassRec"+name);
}

TH2* ptMatrix(TString name){
 return TUnfoldBinningV17::CreateHistogramOfMigrations(ptBinning_gen(),ptBinning_rec(),"hmcPtGenRec" + name);
}

TH2* massMatrix(TString name){
 return TUnfoldBinningV17::CreateHistogramOfMigrations(massBinning_gen(),massBinning_rec(),"hmcMassGenRec" + name);
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
 treco->SetBranchAddress("weightGen",&weightGen);
 treco->SetBranchAddress("weightRec",&weightRec);
 treco->SetBranchAddress("l1PreFire",&l1PreFire);
 treco->SetBranchAddress("nVtx",&nVtx);
 treco->SetBranchAddress("DYtautau",&DYtautau);
 nentries=treco->GetEntries();

 TUnfoldBinningV17 *ptbin = ptBinning_rec();
 TUnfoldBinningV17 *massbin = massBinning_rec();

 // TODO based on the info in recoHistsinfo make map for systematics
 for(int i=0;i<nentries;i++){
 //for(int i=0;i<10000;i++){
   if(i%10000000==0) cout<<i<<endl;
   treco->GetEntry(i);
    if(ispassRec && isBveto && ptRec->at(2) < 100){
          fileout1->cd();
          recoHist.ptHists.at(0)->Fill(ptbin->GetGlobalBinNumber(ptRec->at(2),mRec->at(2)), weightGen*weightRec);
          recoHist.massHists.at(0)->Fill(massbin->GetGlobalBinNumber(mRec->at(2)), weightGen*weightRec);
    }
 }// event loop

 delete ptbin;
 delete massbin;
 delete treco;
 //fileout->cd();
 
 // seems ptHistograms automatically written 
 // recoHist.ptHists.at(0)->Write();
 // delete outputFile;
}

void sigHists(TFile *filein, TFile *fileout1, TFile *fileout2, const sigHistsinfo &sigHist, const recoHistsinfo &recoHist){ // TODO add list of systematics

 gROOT->SetBatch();
 TH1::SetDefaultSumw2();

 TTree *tsignal=(TTree *)filein->Get("tree");
 tsignal->SetBranchAddress("qLep",&qLep);
 tsignal->SetBranchAddress("ptPreFSR",&ptPreFSR);
 tsignal->SetBranchAddress("mPreFSR",&mPreFSR);
 tsignal->SetBranchAddress("ptDRp2FSR",&ptDRp2FSR);
 tsignal->SetBranchAddress("mDRp2FSR",&mDRp2FSR);
 tsignal->SetBranchAddress("ptRec",&ptRec);
 tsignal->SetBranchAddress("etaRec",&etaRec);
 tsignal->SetBranchAddress("phiRec",&phiRec);
 tsignal->SetBranchAddress("mRec",&mRec);
 tsignal->SetBranchAddress("ispassRec",&ispassRec);
 tsignal->SetBranchAddress("isBveto",&isBveto);
 tsignal->SetBranchAddress("weightGen",&weightGen);
 tsignal->SetBranchAddress("weightRec",&weightRec);
 tsignal->SetBranchAddress("l1PreFire",&l1PreFire);
 tsignal->SetBranchAddress("nVtx",&nVtx);
 tsignal->SetBranchAddress("DYtautau",&DYtautau);
 tsignal->SetBranchAddress("isfiducialPreFSR",&isfiducialPreFSR);

 // FIXME 
 TFile* fZptWeight = new TFile("/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/etc/ZptWeight/aMCNLO/electron/ZptWeight_electron.root", "r");
 //TFile* fZptWeight = new TFile("/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/etc/ZptWeight/MG/electron/ZptWeight_electron.root", "r");
 TH1 *hZptWeight;
 fZptWeight->GetObject("ZptWeight", hZptWeight);

 nentries=tsignal->GetEntries();

 TUnfoldBinningV17 *ptBin_rec = ptBinning_rec();
 TUnfoldBinningV17 *ptBin_gen = ptBinning_gen();

 TUnfoldBinningV17 *massBin_rec = massBinning_rec();
 TUnfoldBinningV17 *massBin_gen = massBinning_gen();

 // TODO based on the info in recoHistsinfo make map for systematics
 for(int i=0;i<nentries;i++){
 //for(int i=0;i<10000;i++){
   if(i%10000000==0) cout<<i<<endl;
   tsignal->GetEntry(i);

    if(!sigHist.isInc){
       if(ispassRec && isBveto && ptRec->at(2) < 100){
             fileout1->cd();
             sigHist.ptHists.at(0)->Fill(ptBin_rec->GetGlobalBinNumber(ptRec->at(2),mRec->at(2)), weightGen*weightRec); //
             sigHist.massHists.at(0)->Fill(massBin_rec->GetGlobalBinNumber(mRec->at(2)), weightGen*weightRec); //
       }
    }
    else{ // for inclusive DY MC

         if(ispassRec && isBveto && ptRec->at(2) < 100){
            if(DYtautau){
               fileout2->cd();
               recoHist.ptHists.at(0)->Fill(ptBin_rec->GetGlobalBinNumber(ptRec->at(2),mRec->at(2)), weightGen*weightRec*0.97*l1PreFire->at(0));
               recoHist.massHists.at(0)->Fill(massBin_rec->GetGlobalBinNumber(mRec->at(2)), weightGen*weightRec*0.97*l1PreFire->at(0));
            }
            else{
               fileout1->cd();
               sigHist.ptHists.at(0)->Fill(ptBin_rec->GetGlobalBinNumber(ptRec->at(2),mRec->at(2)), weightGen*weightRec*0.97*l1PreFire->at(0));
               sigHist.massHists.at(0)->Fill(massBin_rec->GetGlobalBinNumber(mRec->at(2)), weightGen*weightRec*0.97*l1PreFire->at(0));
            }
         }

        /////////////////////////////////////////// fill migration matrix ////////////////////////////////////////
        //
        if(!DYtautau){
          if(ptPreFSR->size()==3){

            Double_t zptWeight = 1.;
            Double_t ZptGen_ = ptPreFSR->at(2);
            if(ZptGen_ > 100.) ZptGen_ = 99.5;
            zptWeight = hZptWeight->GetBinContent(hZptWeight->FindBin(ZptGen_));

            if( ispassRec && isBveto && ptRec->at(2) < 100) {
               int ptBinZero=0;
               sigHist.ptMatrixs.at(0)->Fill(ptBin_gen->GetGlobalBinNumber(ptPreFSR->at(2), mPreFSR->at(2)), ptBin_rec->GetGlobalBinNumber(ptRec->at(2),mRec->at(2)), weightGen*weightRec*0.97*l1PreFire->at(0));
               sigHist.ptMatrixs.at(0)->Fill(ptBin_gen->GetGlobalBinNumber(ptPreFSR->at(2), mPreFSR->at(2)), ptBinZero,                                               weightGen*(1.-weightRec*0.97*l1PreFire->at(0)));

               // for the migration matrix inside the fiducial region at pre FSR
               sigHist.ptMatrixs.at(1)->Fill(ptBin_gen->GetGlobalBinNumber(ptPreFSR->at(2), mPreFSR->at(2)), ptBin_rec->GetGlobalBinNumber(ptRec->at(2),mRec->at(2)), weightGen*weightRec*0.97*l1PreFire->at(0));
               sigHist.ptMatrixs.at(1)->Fill(ptBin_gen->GetGlobalBinNumber(ptPreFSR->at(2), mPreFSR->at(2)), ptBinZero,                                               weightGen*(1.-weightRec*0.97*l1PreFire->at(0)));

               sigHist.ptMatrixs.at(2)->Fill(ptBin_gen->GetGlobalBinNumber(ptPreFSR->at(2), mPreFSR->at(2)), ptBin_rec->GetGlobalBinNumber(ptRec->at(2),mRec->at(2)), weightGen*zptWeight*weightRec*0.97*l1PreFire->at(0));
               sigHist.ptMatrixs.at(2)->Fill(ptBin_gen->GetGlobalBinNumber(ptPreFSR->at(2), mPreFSR->at(2)), ptBinZero,                                               weightGen*zptWeight*(1.-weightRec*0.97*l1PreFire->at(0)));

               sigHist.ptMatrixs.at(3)->Fill(ptBin_gen->GetGlobalBinNumber(ptDRp2FSR->at(2), mDRp2FSR->at(2)), ptBin_rec->GetGlobalBinNumber(ptRec->at(2),mRec->at(2)), weightGen*weightRec*0.97*l1PreFire->at(0));
               sigHist.ptMatrixs.at(3)->Fill(ptBin_gen->GetGlobalBinNumber(ptDRp2FSR->at(2), mDRp2FSR->at(2)), ptBinZero,                                               weightGen*(1.-weightRec*0.97*l1PreFire->at(0)));

               int massBinZero=0;
               sigHist.massMatrixs.at(0)->Fill(massBin_gen->GetGlobalBinNumber(mPreFSR->at(2)), massBin_rec->GetGlobalBinNumber(mRec->at(2)), weightGen*weightRec*0.97*l1PreFire->at(0));
               sigHist.massMatrixs.at(0)->Fill(massBin_gen->GetGlobalBinNumber(mPreFSR->at(2)), massBinZero,                                  weightGen*(1.-weightRec*0.97*l1PreFire->at(0)));

               // for the migration matrix inside the fiducial region at pre FSR
               sigHist.massMatrixs.at(1)->Fill(massBin_gen->GetGlobalBinNumber(mPreFSR->at(2)), massBin_rec->GetGlobalBinNumber(mRec->at(2)), weightGen*weightRec*0.97*l1PreFire->at(0));
               sigHist.massMatrixs.at(1)->Fill(massBin_gen->GetGlobalBinNumber(mPreFSR->at(2)), massBinZero,                                  weightGen*(1.-weightRec*0.97*l1PreFire->at(0)));

               sigHist.massMatrixs.at(2)->Fill(massBin_gen->GetGlobalBinNumber(mPreFSR->at(2)), massBin_rec->GetGlobalBinNumber(mRec->at(2)), weightGen*zptWeight*weightRec*0.97*l1PreFire->at(0));
               sigHist.massMatrixs.at(2)->Fill(massBin_gen->GetGlobalBinNumber(mPreFSR->at(2)), massBinZero,                                  weightGen*zptWeight*(1.-weightRec*0.97*l1PreFire->at(0)));

               sigHist.massMatrixs.at(3)->Fill(massBin_gen->GetGlobalBinNumber(mDRp2FSR->at(2)), massBin_rec->GetGlobalBinNumber(mRec->at(2)), weightGen*weightRec*0.97*l1PreFire->at(0));
               sigHist.massMatrixs.at(3)->Fill(massBin_gen->GetGlobalBinNumber(mDRp2FSR->at(2)), massBinZero,                                  weightGen*(1.-weightRec*0.97*l1PreFire->at(0)));

            }// events passing reco selection
            else{
                int ptBinZero=0;
                sigHist.ptMatrixs.at(0)->Fill(ptBin_gen->GetGlobalBinNumber(ptPreFSR->at(2), mPreFSR->at(2)), ptBinZero, weightGen); 
                sigHist.ptMatrixs.at(2)->Fill(ptBin_gen->GetGlobalBinNumber(ptPreFSR->at(2), mPreFSR->at(2)), ptBinZero, weightGen*zptWeight); 
                sigHist.ptMatrixs.at(3)->Fill(ptBin_gen->GetGlobalBinNumber(ptDRp2FSR->at(2), mDRp2FSR->at(2)), ptBinZero, weightGen); 

                int massBinZero=0;
                sigHist.massMatrixs.at(0)->Fill(massBin_gen->GetGlobalBinNumber(mPreFSR->at(2)), ptBinZero, weightGen); 
                sigHist.massMatrixs.at(2)->Fill(massBin_gen->GetGlobalBinNumber(mPreFSR->at(2)), ptBinZero, weightGen*zptWeight); 
                sigHist.massMatrixs.at(3)->Fill(massBin_gen->GetGlobalBinNumber(mDRp2FSR->at(2)), ptBinZero, weightGen); 
                
                if(isfiducialPreFSR){
                  int ptBinZero=0;
                  sigHist.ptMatrixs.at(1)->Fill(ptBin_gen->GetGlobalBinNumber(ptPreFSR->at(2), mPreFSR->at(2)), ptBinZero, weightGen); 

                  int massBinZero=0;
                  sigHist.massMatrixs.at(1)->Fill(massBin_gen->GetGlobalBinNumber(mPreFSR->at(2)), ptBinZero, weightGen); 
                }
            }
          }// check if the pre FSR dilepton exist
        }// DY to ee or mumu events only
        //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    }// for DYtoLL MC case
 }// event loop

 delete hZptWeight;
 delete fZptWeight;
 delete ptBin_rec;
 delete ptBin_gen;
 delete massBin_rec;
 delete massBin_gen;
 delete tsignal;
 //fileout->cd();

 // seems ptHistograms automatically written 
 // recoHist.ptHists.at(0)->Write();
 // delete outputFile;
}

