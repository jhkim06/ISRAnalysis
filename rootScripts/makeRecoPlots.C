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
struct recoTH1info {
  std::vector<TH1*> hists;
  recoTH1info(std::vector<TH1*> ihists):
    hists(std::move(ihists)){}
};

void recoHists(TFile *filein){ // TODO add list of systematics

 TH1::SetDefaultSumw2();

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
 binning_Rec->PrintStream(std::cout);

 TUnfoldBinningV17 *binning_Gen=new TUnfoldBinningV17("Gen");
 binning_Gen->AddAxis("pt",nptbin_wide,ptbin_wide,false,true);
 binning_Gen->AddAxis("mass",nmassbin_wide,massbin_wide,true,true);
 binning_Gen->PrintStream(std::cout);

 std::vector<TH1*> hp;
 recoTH1info histsToSave(hp);
 histsToSave.hists.push_back((binning_Rec->CreateHistogram("hdataRec")));

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
 nentries=treco->GetEntries();

 for(int i=0;i<nentries;i++){
   if(i%10000000==0) cout<<i<<endl;
   treco->GetEntry(i);
    if(ispassRec && isBveto && ptRec->at(2) < 100){
       histsToSave.hists.at(0)->Fill(binning_Rec->GetGlobalBinNumber(ptRec->at(2),mRec->at(2)));
    }
 }

 TFile *outputFile=new TFile("test.root","recreate");
 histsToSave.hists.at(0)->Write();
 outputFile->Write();

 delete outputFile;
}

/*
 TFile *outputFile=new TFile("test.root","recreate");

 outputFile->cd();
 TUnfoldBinning *fineBinningRoot,*coarseBinningRoot;

 TDOMParser parser;
 Int_t error=parser.ParseFile("/home/jhkim/ISR2016/unfolding/TUnfold/testUnfold7binning.xml");
  if(error) {
     cout<<"error="<<error<<" from TDOMParser\n";
     cout<<"==============================================================\n";
     cout<<"Maybe the file testUnfold7binning.xml is missing?\n";
     cout<<"The content of the file is included in the comments section\n";
     cout<<"of this macro \"testUnfold7b.C\"\n";
     cout<<"==============================================================\n";
 }

 TXMLDocument const *XMLdocument=parser.GetXMLDocument();
 fineBinningRoot=
    TUnfoldBinningXML::ImportXML(XMLdocument,"fine");
 coarseBinningRoot=
    TUnfoldBinningXML::ImportXML(XMLdocument,"coarse");
 
 // write binnig schemes to output file
 fineBinningRoot->Write();
 coarseBinningRoot->Write();

 if(fineBinningRoot) {
    fineBinningRoot->PrintStream(cout);
 } else {
    cout<<"could not read 'detector' binning\n";
 }
 if(coarseBinningRoot) {
    coarseBinningRoot->PrintStream(cout);
 } else {
    cout<<"could not read 'generator' binning\n";
 }

*/
