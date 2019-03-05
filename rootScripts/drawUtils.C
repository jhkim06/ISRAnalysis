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
#include <TMathText.h>
#include <TLatex.h>
#include <TDOMParser.h>
#include <TXMLDocument.h>
#include "TUnfoldBinningXML.h"
#include "TUnfoldBinning.h"
#include "TUnfoldDensity.h"

#include "makeRecoPlots.h"

void setYaxisHist(TH1* hist){

        hist->SetTitle("");
        hist->GetYaxis()->SetLabelFont(63);
        hist->GetYaxis()->SetLabelSize(25);
        hist->GetYaxis()->SetTitleFont(43);
        hist->GetYaxis()->SetTitleSize(20);
        hist->GetYaxis()->SetTitleOffset(2.5);
}

void setYaxisGraph(TGraphErrors* gr){

        gr->SetTitle("");
        gr->GetYaxis()->SetLabelFont(63);
        gr->GetYaxis()->SetLabelSize(25);
        gr->GetYaxis()->SetTitleFont(43);
        gr->GetYaxis()->SetTitleSize(20);
        gr->GetYaxis()->SetTitleOffset(2.5);
}

void setXaxisHist(TH1* hist){

        hist->GetXaxis()->SetLabelFont(63);
        hist->GetXaxis()->SetLabelSize(25);
        hist->GetXaxis()->SetTitleOffset(3.5);
        hist->GetXaxis()->SetLabelOffset(0.02);
        hist->GetXaxis()->SetTitleFont(43);
        hist->GetXaxis()->SetTitleSize(20);
}

void setXaxisGraph(TGraphErrors* gr){

        gr->GetXaxis()->SetLabelFont(63);
        gr->GetXaxis()->SetLabelSize(25);
        gr->GetXaxis()->SetTitleOffset(3.5);
        gr->GetXaxis()->SetLabelOffset(0.02);
        gr->GetXaxis()->SetTitleFont(43);
        gr->GetXaxis()->SetTitleSize(20);
}

void setPadMargins(TPad* pad){

        pad->SetTopMargin(0.05);
        pad->SetLeftMargin(0.12);
        pad->SetRightMargin(0.05);

}

TLine* drawVerLine(Double_t x1, Double_t y1, Double_t x2, Double_t y2){

        TLine *l_;
        l_ = new TLine(x1, y1, x2, y2);
        l_->Draw("same");
        l_->SetLineStyle(2);
        l_->SetLineColor(kBlack);
        return l_;
}

TH1* getHistwAxisBin(TH1* hist, int nthMass){

	int nbinxTot = hist->GetNbinsX();
	int nbinxSub = nbinxTot / 5; // assume the total number of mass bin is 5, and equidistant binning definition here. 
	TH1* temp = new TH1D("temp","temp",nbinxSub, 0., 100.);

	for(int ibin = 1; ibin <=nbinxSub; ibin++){
                int ibin_ = nbinxSub * (nthMass-1) + ibin;
		temp->SetBinContent(ibin, hist->GetBinContent(ibin_));	
		temp->SetBinError(ibin, hist->GetBinError(ibin_));	
	} 

        return temp;
}

void drawRatio(TString outpdf, TUnfoldDensityV17* unfold_pt, TUnfoldDensityV17* unfold_mass, TFile *filein){

        gROOT->SetBatch();

	TH1* hunfolded_pt = unfold_pt->GetOutput("hunfolded_pt",0,0,"*[UO]",kFALSE);
	TH1* ratio=(TH1*)hunfolded_pt->Clone("ratio");

        TH1 *histMCTruth_pt=unfold_pt->GetBias("histMCTruth_pt",0,0,"*[UO]",kFALSE);
        ratio->Divide(histMCTruth_pt);

        TH1* hunfolded_mass = unfold_mass->GetOutput("hunfolded_mass",0,0,"*[UO]",kTRUE);
        TH1 *histMCTruth_mass=unfold_mass->GetBias("histMCTruth_mass",0,0,"*[UO]",kTRUE);

        TH1* hmeans = new TH1D("means", "means", 5, 0., 5.);
        TH1* hmeansMC = new TH1D("meansMC", "meansMC", 5, 0., 5.);

        Double_t meanpt_data[5], meanpterr_data[5];
        Double_t meanpt_mc[5], meanpterr_mc[5];

        Double_t meanmass_data[5], meanmasserr_data[5];
        Double_t meanmass_mc[5], meanmasserr_mc[5];

        // get average pT for each mass region
        for(int i = 1; i < 6; i++){
        	TH1* temp;
                TH1* tempMC;

		temp = getHistwAxisBin(hunfolded_pt, i);
                hmeans->SetBinContent(i, temp->GetMean());
                hmeans->SetBinError(i, temp->GetMeanError());
                meanpt_data[i-1] = temp->GetMean();
                meanpterr_data[i-1] = temp->GetMeanError();
        	delete temp;

		tempMC = getHistwAxisBin(histMCTruth_pt, i);
                hmeansMC->SetBinContent(i, tempMC->GetMean());
                hmeansMC->SetBinError(i, tempMC->GetMeanError());
                meanpt_mc[i-1] = tempMC->GetMean();
                meanpterr_mc[i-1] = tempMC->GetMeanError();
                delete tempMC;
        }
	
	// get average mass
        TUnfoldBinningV17* massBin_gen = massBinning_gen();
	TH1 *hunfolded_mass_bin[5];
	TH1 *hmc_mass_bin[5];
        hunfolded_mass_bin[0] = massBin_gen->CreateHistogram("hunfolded_mass_m1",true);
        hunfolded_mass_bin[1] = massBin_gen->CreateHistogram("hunfolded_mass_m2",true);
        hunfolded_mass_bin[2] = massBin_gen->CreateHistogram("hunfolded_mass_m3",true);
        hunfolded_mass_bin[3] = massBin_gen->CreateHistogram("hunfolded_mass_m4",true);
        hunfolded_mass_bin[4] = massBin_gen->CreateHistogram("hunfolded_mass_m5",true);

        hmc_mass_bin[0] = massBin_gen->CreateHistogram("hmc_mass_m1",true);
        hmc_mass_bin[1] = massBin_gen->CreateHistogram("hmc_mass_m2",true);
        hmc_mass_bin[2] = massBin_gen->CreateHistogram("hmc_mass_m3",true);
        hmc_mass_bin[3] = massBin_gen->CreateHistogram("hmc_mass_m4",true);
        hmc_mass_bin[4] = massBin_gen->CreateHistogram("hmc_mass_m5",true);

        for(int ibin = 1; ibin<hunfolded_mass->GetNbinsX()+1; ibin++){

           int massRegion = 0;
           if(hunfolded_mass->GetBinLowEdge(ibin) >= 40 && hunfolded_mass->GetBinLowEdge(ibin)+hunfolded_mass->GetBinWidth(ibin) <= 60)
              massRegion = 0;
           if(hunfolded_mass->GetBinLowEdge(ibin) >= 60 && hunfolded_mass->GetBinLowEdge(ibin)+hunfolded_mass->GetBinWidth(ibin) <= 80)
              massRegion = 1;
           if(hunfolded_mass->GetBinLowEdge(ibin) >= 80 && hunfolded_mass->GetBinLowEdge(ibin)+hunfolded_mass->GetBinWidth(ibin) <= 100)
              massRegion = 2;
           if(hunfolded_mass->GetBinLowEdge(ibin) >= 100 && hunfolded_mass->GetBinLowEdge(ibin)+hunfolded_mass->GetBinWidth(ibin) <= 200)
              massRegion = 3;
           if(hunfolded_mass->GetBinLowEdge(ibin) >= 200 && hunfolded_mass->GetBinLowEdge(ibin)+hunfolded_mass->GetBinWidth(ibin) <= 350)
              massRegion = 4;

           hunfolded_mass_bin[massRegion]->SetBinContent(ibin, hunfolded_mass->GetBinContent(ibin));
           hunfolded_mass_bin[massRegion]->SetBinError(ibin, hunfolded_mass->GetBinError(ibin));

           hmc_mass_bin[massRegion]->SetBinContent(ibin, histMCTruth_mass->GetBinContent(ibin));
           hmc_mass_bin[massRegion]->SetBinError(ibin, histMCTruth_mass->GetBinError(ibin));
        }

        delete massBin_gen;
	for(int i = 0; i < 5; i++){

           meanmass_data[i] = hunfolded_mass_bin[i]->GetMean();
           meanmasserr_data[i] = hunfolded_mass_bin[i]->GetMeanError();

           meanmass_mc[i] = hmc_mass_bin[i]->GetMean();
           meanmasserr_mc[i] = hmc_mass_bin[i]->GetMeanError();

           delete hunfolded_mass_bin[i];
           delete hmc_mass_bin[i];
        }

	/* print out values
        for(int i = 0; i < 5; i++){

            std::cout <<" pt: " << meanpt_data[i] << " mass: " << meanmass_data[i] << std::endl;

        }
 	*/

  	TCanvas* c1=new TCanvas("c1", "c1", 50, 50, 700*1.5, 1000*1.5);
  	c1->cd();
  	gStyle->SetOptStat(0);

  	TPad *pad1 = new TPad("pad1","pad1",0,0.66,1,1);
        setPadMargins(pad1);
  	pad1->SetBottomMargin(0.01);
  	pad1->SetTicks(1);
  	pad1->SetLogy();
  	pad1->Draw();
  	pad1->cd();

  	hunfolded_pt->SetTitle("");
  	hunfolded_pt->Draw("p9histe");
  	hunfolded_pt->SetMarkerStyle(20);
  	hunfolded_pt->SetMarkerSize(.7);
  	hunfolded_pt->SetLineColor(kBlack);
        setYaxisHist(hunfolded_pt);
  	hunfolded_pt->GetXaxis()->SetLabelSize(0.);
  	hunfolded_pt->GetYaxis()->SetTitle("Number of events per bin");
  	histMCTruth_pt->Draw("histsames");
  	histMCTruth_pt->SetLineColor(2);
  	hunfolded_pt->GetYaxis()->SetRangeUser(100.,hunfolded_pt->GetMaximum()>histMCTruth_pt->GetMaximum()?10.*hunfolded_pt->GetMaximum():10.*histMCTruth_pt->GetMaximum());

  	TLegend* leg_ = new TLegend(0.7, 0.60, 0.95, 0.9,"","brNDC");
  	leg_->SetTextSize(0.06);
  	leg_->SetFillStyle(0);
  	leg_->SetBorderSize(0);
  	leg_->AddEntry(hunfolded_pt, "Unfolded reco", "p");
  	leg_->AddEntry(histMCTruth_pt, "Truth", "l");
  	leg_->Draw();

  	c1->cd();

  	TPad *pad2 = new TPad("pad2","pad2",0,0.33,1,0.66);
        setPadMargins(pad2);
  	pad2->SetBottomMargin(0.2);
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
        setYaxisHist(ratio);
        setXaxisHist(ratio);
  	ratio->GetYaxis()->SetNdivisions(505);

        TLine *l_;
        l_ = new TLine(ratio->GetXaxis()->GetXmin(),1,ratio->GetXaxis()->GetXmax(),1);
        l_->Draw("same");
        l_->SetLineStyle(1);
        l_->SetLineColor(kRed);
	
	TLine* lm1 = drawVerLine(20, 0.8, 20, 1.2);
	TLine* lm2 = drawVerLine(40, 0.8, 40, 1.2);
	TLine* lm3 = drawVerLine(60, 0.8, 60, 1.2);
	TLine* lm4 = drawVerLine(80, 0.8, 80, 1.2);
	TLine* lm5 = drawVerLine(100, 0.8, 100, 1.2);

        c1->cd();

        TPad *pad3 = new TPad("pad3","pad3",0,0.,1,0.33);
        setPadMargins(pad3);
        pad3->SetBottomMargin(0.3);
        pad3->SetTicks(1);
        pad3->SetGrid(0,1);
        pad3->Draw();
        pad3->cd();

        TGraphErrors *grUnfolded = new TGraphErrors(5, meanmass_data, meanpt_data, meanmasserr_data, meanpterr_data);
        TGraphErrors *grMC = new TGraphErrors(5, meanmass_mc, meanpt_mc, meanmasserr_mc, meanpterr_mc);

        grUnfolded->Draw("ape");
        grUnfolded->SetTitle("Average \\ {p_{T}^{\\ell\\ell}} vs M{\\ell\\ell}");
        grUnfolded->SetMarkerStyle(20);
        grUnfolded->SetMarkerSize(.9);
        grUnfolded->SetLineColor(kBlack);
        grUnfolded->SetMinimum(10.);
        grUnfolded->SetMaximum(30.);
     	setYaxisGraph(grUnfolded);
	setXaxisGraph(grUnfolded);
        grUnfolded->GetYaxis()->SetTitle("Average \\ {p_{T}^{\\ell\\ell}}");

        grMC->Draw("samepe");
        grMC->SetMarkerStyle(20);
        grMC->SetMarkerSize(.9);
        grMC->SetLineColor(kRed);
        grMC->SetMarkerColor(kRed);

/*
        hmeans->Draw("pe");
        hmeans->SetTitle("");
        hmeans->SetMarkerStyle(20);
        hmeans->SetMarkerSize(.7);
        hmeans->SetLineColor(kBlack);
        hmeans->SetMinimum(10.);
        hmeans->SetMaximum(30.);
        setYaxisHist(hmeans);
        setXaxisHist(hmeans);
        hmeans->GetYaxis()->SetNdivisions(505);
        hmeans->SetYTitle("Average \\ {p_{T}^{\\ell\\ell}}"); // \\ell not working with pdf format
	hmeans->GetXaxis()->SetTitle("Mass(m_{#font[11]{ee}}) region (GeV)");
        hmeans->GetXaxis()->SetTitleOffset(5.);
        hmeans->GetXaxis()->SetBinLabel(hmeans->GetXaxis()->FindBin(0.), "40~60");
        hmeans->GetXaxis()->SetBinLabel(hmeans->GetXaxis()->FindBin(1.), "60~80");
        hmeans->GetXaxis()->SetBinLabel(hmeans->GetXaxis()->FindBin(2.), "80~100");
        hmeans->GetXaxis()->SetBinLabel(hmeans->GetXaxis()->FindBin(3.), "100~200");
        hmeans->GetXaxis()->SetBinLabel(hmeans->GetXaxis()->FindBin(4.), "200~350");

        hmeansMC->Draw("pesame");
        hmeansMC->SetMarkerStyle(20);
        hmeansMC->SetMarkerSize(.7);
        hmeansMC->SetLineColor(kRed);
        hmeansMC->SetMarkerColor(kRed);
*/
        c1->SaveAs(outpdf);
	
	delete grUnfolded; 
	delete grMC; 
        delete hmeans;
        delete hmeansMC;
        delete l_; 
        delete lm1;
        delete lm2;
	delete lm3;
        delete lm4;
        delete lm5;
        delete leg_;
        delete pad1;
        delete pad2;
        delete pad3;
        delete c1;
}

void drawMassRatio(TString outpdf, TUnfoldDensityV17* unfold, TFile *filein){

        gROOT->SetBatch();
	TH1* hunfolded = unfold->GetOutput("hunfolded",0,0,"*[UO]",kTRUE);
	TH1* ratio=(TH1*)hunfolded->Clone("ratio");

        TH1 *histMCTruth=unfold->GetBias("histMCTruth",0,0,"*[UO]",kTRUE);
        ratio->Divide(histMCTruth);

  	TCanvas* c1=new TCanvas("c1", "c1", 50, 50, 700*1.5, 1000*1.5);
  	c1->cd();
  	gStyle->SetOptStat(0);

  	TPad *pad1 = new TPad("pad1","pad1",0,0.66,1,1);
        setPadMargins(pad1);
  	pad1->SetBottomMargin(0.01);
  	pad1->SetTicks(1);
  	pad1->SetLogy();
  	pad1->Draw();
  	pad1->cd();

  	hunfolded->SetTitle("");
  	hunfolded->Draw("p9histe");
  	hunfolded->SetMarkerStyle(20);
  	hunfolded->SetMarkerSize(.7);
  	hunfolded->SetLineColor(kBlack);
        setYaxisHist(hunfolded);
  	hunfolded->GetXaxis()->SetLabelSize(0.);
  	hunfolded->GetYaxis()->SetTitle("Number of events per bin");
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

  	TPad *pad2 = new TPad("pad2","pad2",0,0.33,1,0.66);
        setPadMargins(pad2);
  	pad2->SetBottomMargin(0.2);
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
        setYaxisHist(ratio);
        setXaxisHist(ratio);
  	ratio->GetYaxis()->SetNdivisions(505);

        TLine *l_;
        l_ = new TLine(ratio->GetXaxis()->GetXmin(),1,ratio->GetXaxis()->GetXmax(),1);
        l_->Draw("same");
        l_->SetLineStyle(1);
        l_->SetLineColor(kRed);
	
	//TLine* lm1 = drawVerLine(20, 0.8, 20, 1.2);
	//TLine* lm2 = drawVerLine(40, 0.8, 40, 1.2);
	//TLine* lm3 = drawVerLine(60, 0.8, 60, 1.2);
	//TLine* lm4 = drawVerLine(80, 0.8, 80, 1.2);
	//TLine* lm5 = drawVerLine(100, 0.8, 100, 1.2);

        c1->cd();
/*
        TPad *pad3 = new TPad("pad3","pad3",0,0.,1,0.33);
        setPadMargins(pad3);
        pad3->SetBottomMargin(0.3);
        pad3->SetTicks(1);
        pad3->SetGrid(0,1);
        pad3->Draw();
        pad3->cd();

        hmeans->Draw("pe");
        hmeans->SetTitle("");
        hmeans->SetMarkerStyle(20);
        hmeans->SetMarkerSize(.7);
        hmeans->SetLineColor(kBlack);
        hmeans->SetMinimum(10.);
        hmeans->SetMaximum(30.);
        setYaxisHist(hmeans);
        setXaxisHist(hmeans);
        hmeans->GetYaxis()->SetNdivisions(505);
        hmeans->SetYTitle("Average \\ {p_{T}^{\\ell\\ell}}"); // \\ell not working with pdf format
	hmeans->GetXaxis()->SetTitle("Mass(m_{#font[11]{ee}}) region (GeV)");
        hmeans->GetXaxis()->SetTitleOffset(5.);
        hmeans->GetXaxis()->SetBinLabel(hmeans->GetXaxis()->FindBin(0.), "40~60");
        hmeans->GetXaxis()->SetBinLabel(hmeans->GetXaxis()->FindBin(1.), "60~80");
        hmeans->GetXaxis()->SetBinLabel(hmeans->GetXaxis()->FindBin(2.), "80~100");
        hmeans->GetXaxis()->SetBinLabel(hmeans->GetXaxis()->FindBin(3.), "100~200");
        hmeans->GetXaxis()->SetBinLabel(hmeans->GetXaxis()->FindBin(4.), "200~350");

        hmeansMC->Draw("pesame");
        hmeansMC->SetMarkerStyle(20);
        hmeansMC->SetMarkerSize(.7);
        hmeansMC->SetLineColor(kRed);
        hmeansMC->SetMarkerColor(kRed);
*/
        c1->SaveAs(outpdf);

        //delete hmeans;
        //delete hmeansMC;
        delete l_; 
        //delete lm1;
        //delete lm2;
	//delete lm3;
        //delete lm4;
        //delete lm5;
        delete leg_;
        delete pad1;
        delete pad2;
        //delete pad3;
        delete c1;
}

void responseM(TString outpdf, TUnfoldDensityV17* unfold){
	gROOT->SetBatch();
	TH2 *histProbability=unfold->GetProbabilityMatrix("histProbability"); // get response matrix

        TCanvas *c1 = new TCanvas("c1", "c1", 50, 50, 750, 750);
        gStyle->SetOptStat(0);
        gStyle->SetLineWidth(1);
        c1->SetBottomMargin(0.15);
        c1->SetTopMargin(0.05);
        c1->SetLeftMargin(0.15);
        c1->SetRightMargin(0.15);
        c1->SetTicks(1);
        c1->Draw();
        c1->cd();

        //gPad->SetLogz();  
        gStyle->SetPalette(1);

        histProbability->GetZaxis()->SetRangeUser(0., 1.);
        histProbability->Draw("colz");

        c1->SaveAs(outpdf);
     
        delete histProbability;
        delete c1;
}

void efficiency(TString outpdf, TUnfoldDensityV17* unfold){

        gROOT->SetBatch();
        TH2 *histProbability=unfold->GetProbabilityMatrix("histProbability"); // get response matrix
	TH1 *histEfficiency=histProbability->ProjectionX("histEfficiency");

        TCanvas *c1 = new TCanvas("c1", "c1", 50, 50, 750, 750);
        gStyle->SetOptStat(0);
        gStyle->SetLineWidth(1);
        c1->SetBottomMargin(0.15);
        c1->SetTopMargin(0.05);
        c1->SetLeftMargin(0.15);
        c1->SetRightMargin(0.15);
        c1->SetTicks(1);
        c1->Draw();
        c1->cd();

        histEfficiency->Draw("hist");

        c1->SaveAs(outpdf);

        delete histProbability;
        delete histEfficiency;
        delete c1;

}
