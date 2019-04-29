#ifndef DRAWUTILS_H
#define DRAWUTILS_H

#include <iostream>
#include <map>
#include <cmath>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TGaxis.h>
#include <TMathText.h>
#include <TLatex.h>
#include <THStack.h>
#include <TDOMParser.h>
#include <TXMLDocument.h>
#include "TUnfoldBinningXML.h"
#include "TUnfoldBinning.h"
#include "TUnfoldDensity.h"

using namespace std;

void setYaxisHist(TH1* hist);
void setYaxisGraph(TGraphErrors* gr);
void setXaxisHist(TH1* hist);
void setXaxisGraph(TGraphErrors* gr);
void setTGraphAxis(TGraphErrors* data, Double_t x, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, TGraphErrors* mc = nullptr,  TGraphErrors* sys = nullptr, bool axis = true);
void setPadMargins(TPad* pad);
TLine* drawVerLine(Double_t x1, Double_t y1, Double_t x2, Double_t y2);
void getAveragesMass(vector<Double_t> &mean, vector<Double_t> &err, TH1* hmass);
void getAveragesPt(vector<Double_t> &mean, vector<Double_t> &err, TUnfoldDensity* unfold_pt, bool isData);
void getAveragesSysMass(vector<Double_t> &err, TString sysName, int sysSize, TUnfoldDensity* unfold_mass);
void getAveragesSysPt(vector<Double_t> &err, TString sysName, int sysSize, TUnfoldDensity* unfold_pt, bool relative = false);
void getAveragesSysPt_v2(vector<Double_t> &err, TUnfoldDensity* unfold_pt, TUnfoldDensity* unfold_pt_up, TUnfoldDensity* unfold_pt_dn);
void getAveragesSysPtMC(vector<Double_t> &err, TString sysName, int sysSize, TUnfoldDensity* unfold_pt, TString channel, bool relative = false);
void getRatioSys(TUnfoldDensity* unfold_pt, TString sysName, int sysSize, TH1 *sysRatio);
void getRatioSysMass(TUnfoldDensity* unfold_mass, TString sysName, int sysSize, TH1 *sysRatio);
void drawUnfoldedPtDistWithSys(TString outpdf, TUnfoldDensity* unfold_pt, TString sysName, int sysSize);
void drawRatio(TString outpdf, TUnfoldDensity* unfold_pt, TUnfoldDensity* unfold_mass, TFile *filein, TString channel);
void drawCombinedISR(TString outpdf, TUnfoldDensity* unfold_pt2016, TUnfoldDensity* unfold_mass2016, TUnfoldDensity* unfold_pt2017, TUnfoldDensity* unfold_mass2017);
void drawEMuCombinedISR(TString outpdf, TUnfoldDensity* unfold_ptElectron, TUnfoldDensity* unfold_massElectron, TUnfoldDensity* unfold_ptMuon, TUnfoldDensity* unfold_massMuon);
void drawISRfit(TString outpdf, TUnfoldDensity* unfold_pt, TUnfoldDensity* unfold_mass, TFile *filein);
void drawMassRatio(TString outpdf, TUnfoldDensity* unfold, TFile *filein);
void drawPtReco(TString outpdf, TString postfix, TFile *fdata, TFile *fDYsig, TFile *fDYbkg, TFile *fTTbar, TFile *fVV, TFile *fWjets, TFile *fqcd, TString channel);
void drawMassReco(TString outpdf, TString postfix, TFile *fdata, TFile *fDYsig, TFile *fDYbkg, TFile *fTTbar, TFile *fVV, TFile *fWjets, TFile *fqcd, TString channel);
void responseM(TString outpdf, TUnfoldDensity* unfold);
void efficiency(TString outpdf, TUnfoldDensity* unfold);

TUnfoldDensityV17* getSysMatrix(TString var, TString sysName, TString channel);
void sysMCErrRatio(TUnfoldDensity* unfold_pt, TH1* sysRatio, TString var, TString sysName, int size, TString channel);
void drawSystematicISR(TString outpdf, TUnfoldDensity* unfold_pt, TUnfoldDensity* unfold_mass, TString channel, TUnfoldDensity* unfold_pt_up = NULL, TUnfoldDensity* unfold_pt_dn = NULL);
void drawMCSystematicISR(TString outpdf, TUnfoldDensity* unfold_pt, TUnfoldDensity* unfold_mass, TString channel);

#endif
