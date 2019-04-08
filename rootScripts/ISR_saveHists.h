#ifndef SAVEHISTS_H
#define SAVEHISTS_H

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

#include "ISR_histTUnfold.h"

using namespace std;

void saveRecoHists(TFile *filein, TFile *fileout1, histTUnfold &recoHist, TString channel);
void saveSigHists(TFile *filein, TFile *fileout1, TFile *fileout2, histTUnfold &sigHist, histTUnfold &recoHist, TString channel, Double_t temp_kfactor);

#endif
