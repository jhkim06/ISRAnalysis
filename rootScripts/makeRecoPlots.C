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

void  test(){

 TH1::SetDefaultSumw2();

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

 outputFile->Write();
 delete outputFile;

}
