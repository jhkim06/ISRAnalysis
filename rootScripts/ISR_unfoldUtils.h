#ifndef UNFOLDUTILS_H
#define UNFOLDUTILS_H

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
#include "TUnfoldDensity.h"

TUnfoldDensityV17* setTUnfoldDensity(TFile *filein, TString var, TString matrixName, bool isfsr);
void setVetorSystematic(TFile *filein, TUnfoldDensityV17* unfold, TString var, TString sysMatrixName, int size);
void setInput(TUnfoldDensityV17* unfold, TString var, TString postfix, TFile *filein);
void setInputHist(TUnfoldDensityV17* unfold, TH1* hist);
void subBkgs(TUnfoldDensityV17* unfold, TString var, TString hname, TFile *filein, TString name);
void doUnfold(TUnfoldDensityV17* unfold);
TH1* getHist(TUnfoldDensityV17* unfold, TString postfix);

#endif
