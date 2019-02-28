#include <iostream>
#include <map>
#include <cmath>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TDOMParser.h>
#include <TXMLDocument.h>
#include "TUnfoldBinningXML.h"
#include "TUnfoldBinning.h"
#include "TUnfoldDensity.h"

#include "makeRecoPlots.h"

void drawRatio(TString outpdf, TUnfoldDensityV17* unfold, TFile *filein){

	TH1* hunfolded = unfold->GetOutput("hunfolded",0,0,0,kFALSE);
	TH1* ratio=(TH1*)hunfolded->Clone("ratio");

        TH2* hmcGenRec = (TH2*)filein->Get("hmcGenRecnorminal");
        TH1 *histMCTruth=hmcGenRec->ProjectionX("histMCTruth",0,-1,"e");
        ratio->Divide(histMCTruth);

  	TCanvas* c1=new TCanvas("c1", "c1", 50, 50, 700, 550);
  	c1->cd();
  	gStyle->SetOptFit(kFALSE);
  	gStyle->SetOptStat(0);
  	gStyle->SetTitleFont(32);
  	gStyle->SetLineWidth(1);

  	TPad *pad1 = new TPad("pad1","pad1",0,0.45,1,1);
  	pad1->SetBottomMargin(0.00);
  	pad1->SetTopMargin(0.05);
  	pad1->SetLeftMargin(0.12);
  	pad1->SetRightMargin(0.05);
  	pad1->SetTicks(1);
  	pad1->SetLogy();
  	pad1->Draw();
  	pad1->cd();

  	hunfolded->SetTitle("");
  	hunfolded->Draw("p9histe");
  	hunfolded->SetMarkerStyle(20);
  	hunfolded->SetMarkerSize(.7);
  	hunfolded->SetLineColor(kBlack);
  	hunfolded->GetXaxis()->SetLabelSize(0.);
  	hunfolded->GetYaxis()->SetLabelFont(63);
  	hunfolded->GetYaxis()->SetLabelSize(25);
  	hunfolded->GetYaxis()->SetTitle("Number of events per bin");
  	hunfolded->GetYaxis()->SetTitleFont(43);
  	hunfolded->GetYaxis()->SetTitleSize(20);
  	hunfolded->GetYaxis()->SetTitleOffset(1.3);
  	//histMCTruth->SetFillColor(kRed);
  	histMCTruth->Draw("histsames");
  	histMCTruth->SetLineColor(2);
  	hunfolded->GetYaxis()->SetRangeUser(100.,hunfolded->GetMaximum()>histMCTruth->GetMaximum()?10.*hunfolded->GetMaximum():10.*histMCTruth->GetMaximum());

  	TLegend* leg_ = new TLegend(0.7, 0.60, 0.95, 0.9,"","brNDC");
  	leg_->SetTextSize(0.06);
  	leg_->SetFillStyle(0);
  	leg_->SetBorderSize(0);
  	leg_->AddEntry(hunfolded, "Unfolded reco", "p");
  	leg_->AddEntry(histMCTruth, "Truth", "l");
  	leg_->Draw();

  	c1->cd();

  	TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.45);
  	pad2->SetBottomMargin(0.3);
  	pad2->SetLeftMargin(0.12);
  	pad2->SetRightMargin(0.05);
  	pad2->SetTicks(1);
  	pad2->Draw();
  	pad2->cd();

  	ratio->Draw("pe");
  	ratio->SetMarkerStyle(20);
  	ratio->SetMarkerSize(.7);
  	ratio->SetLineColor(kBlack);
  	ratio->SetMinimum(0.8);
  	ratio->SetMaximum(1.2);
  	ratio->SetTitle("");
  	ratio->GetXaxis()->SetLabelFont(63);
  	ratio->GetXaxis()->SetLabelSize(25);
  	ratio->GetYaxis()->SetLabelFont(63);
  	ratio->GetYaxis()->SetLabelSize(25);
  	ratio->GetYaxis()->SetTitleOffset(1.3);
  	ratio->GetXaxis()->SetTitleOffset(3.5);
  	ratio->GetYaxis()->SetNdivisions(505);
  	ratio->GetXaxis()->SetLabelOffset(0.02);
  	ratio->GetYaxis()->SetTitleFont(43);
  	ratio->GetYaxis()->SetTitleSize(20);
  	ratio->GetXaxis()->SetTitleFont(43);
  	ratio->GetXaxis()->SetTitleSize(20);

        TLine *l_;
        l_ = new TLine(ratio->GetXaxis()->GetXmin(),1,ratio->GetXaxis()->GetXmax(),1);
        l_->Draw("same");
        l_->SetLineStyle(1);
        l_->SetLineColor(kRed);

        c1->cd();
        c1->SaveAs(outpdf);
        delete leg_;
        delete pad1;
        delete pad2;
        delete c1;

}
