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

#include "makeRecoPlots.h"

#include "./tdrstyle.C"
#include "./CMS_lumi.C"

void setYaxisHist(TH1* hist){

        hist->SetTitle("");
        hist->GetYaxis()->SetLabelFont(63);
        hist->GetYaxis()->SetLabelSize(25);
        hist->GetYaxis()->SetTitleFont(63);
        hist->GetYaxis()->SetTitleSize(35);
        hist->GetYaxis()->SetTitleOffset(2.3);
}

void setYaxisGraph(TGraphErrors* gr){

        gr->SetTitle("");
        gr->GetYaxis()->SetLabelFont(63);
        gr->GetYaxis()->SetLabelSize(25);
        gr->GetYaxis()->SetTitleFont(63);
        gr->GetYaxis()->SetTitleSize(35);
        gr->GetYaxis()->SetTitleOffset(2.3);
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

void drawTest(TString outpdf, TUnfoldDensityV17* unfold_pt){

        TH1* hunfolded_pt = unfold_pt->GetOutput("hunfolded_pt",0,0,"pt[UO];mass[UOC2]",kTRUE);

        TCanvas* c1=new TCanvas("c1", "c1", 50, 50, 700*1.5, 1000*1.5);
        c1->cd();
        gStyle->SetOptStat(0);
	hunfolded_pt->Draw("hist");

	c1->SaveAs(outpdf);
        delete c1;
        delete hunfolded_pt;
}

void drawRatio(TString outpdf, TUnfoldDensityV17* unfold_pt, TUnfoldDensityV17* unfold_mass, TFile *filein){

        gROOT->SetBatch();

        setTDRStyle();
        writeExtraText = true;       // if extra text
        extraText  = "Preliminary";

	// pt distribution for 5 mass bins in one histogram
	TH1* hunfolded_pt = unfold_pt->GetOutput("hunfolded_pt",0,0,"*[UO]",kFALSE);
	TH1* ratio=(TH1*)hunfolded_pt->Clone("ratio");

        TH1 *histMCTruth_pt=unfold_pt->GetBias("histMCTruth_pt",0,0,"*[UO]",kFALSE);
        ratio->Divide(histMCTruth_pt);

	// mass distribution
        TH1* hunfolded_mass = unfold_mass->GetOutput("hunfolded_mass",0,0,"*[UO]",kTRUE);
        TH1 *histMCTruth_mass=unfold_mass->GetBias("histMCTruth_mass",0,0,"*[UO]",kTRUE);

        TH1* hmeans = new TH1D("means", "means", 5, 0., 5.);
        TH1* hmeansMC = new TH1D("meansMC", "meansMC", 5, 0., 5.);

        Double_t meanpt_data[5], meanpterr_data[5];
        Double_t meanpt_mc[5], meanpterr_mc[5];

        Double_t meanmass_data[5], meanmasserr_data[5];
        Double_t meanmass_mc[5], meanmasserr_mc[5];

        // get average pT for each mass region
	for(int i = 0; i < 5; i++){
	   TString ibinMass;
	   ibinMass.Form("%d", i);
	   TH1* hunfolded_pt_temp = unfold_pt->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);	
           meanpt_data[i] = hunfolded_pt_temp->GetMean();
           meanpterr_data[i] = hunfolded_pt_temp->GetMeanError();

	   TH1 *histMCTruth_pt_temp=unfold_pt->GetBias("histMCTruth_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
           meanpt_mc[i] = histMCTruth_pt_temp->GetMean();
           meanpterr_mc[i] = histMCTruth_pt_temp->GetMeanError();

           delete hunfolded_pt_temp;
           delete histMCTruth_pt_temp; 
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
  	pad1->SetTopMargin(0.1);
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
        pad1->Update();

	TLine grid_;
	grid_.SetLineColor(kRed);
	grid_.SetLineStyle(kSolid);
	for( int ii=0; ii<histMCTruth_pt->GetNbinsX(); ii++ )
	{
		Int_t i_bin = ii+1;
		Double_t binEdge = hunfolded_pt->GetBinLowEdge(i_bin);
		grid_.DrawLine(binEdge, 100, binEdge, histMCTruth_pt->GetBinContent(ii+1) );
	}

        std::cout << "n bins: " << ratio->GetXaxis()->GetNbins() << std::endl;	
	int totoalNbin = ratio->GetXaxis()->GetNbins();

        TLine grid_vertical;
        grid_vertical.SetLineColor(kBlack);
        grid_vertical.SetLineStyle(kSolid);
        for( int ii=0; ii<5; ii++ )
        {
                grid_vertical.DrawLine(totoalNbin/5 * (ii+1) + 0.5, hunfolded_pt->GetMinimum(), totoalNbin/5 * (ii+1) + 0.5, hunfolded_pt->GetMaximum() );
        }


  	TLegend* leg_ = new TLegend(0.7, 0.60, 0.95, 0.9,"","brNDC");
  	leg_->SetTextSize(0.06);
  	leg_->SetFillStyle(0);
  	leg_->SetBorderSize(0);
  	leg_->AddEntry(hunfolded_pt, "Unfolded data", "p");
  	leg_->AddEntry(histMCTruth_pt, "Truth", "l");
  	leg_->Draw();

  	c1->cd();

  	TPad *pad2 = new TPad("pad2","pad2",0,0.33,1,0.66);
        setPadMargins(pad2);
  	pad2->SetBottomMargin(0.2);
  	pad2->SetTicks(1);
  	pad2->SetGridy(1);
  	pad2->Draw();
  	pad2->cd();

  	ratio->Draw("pe");
  	ratio->SetMarkerStyle(20);
  	ratio->SetMarkerSize(.7);
  	ratio->SetLineColor(kBlack);
  	ratio->SetMinimum(0.6);
  	ratio->SetMaximum(1.4);
  	ratio->SetTitle("");
        setYaxisHist(ratio);
        setXaxisHist(ratio);
  	ratio->GetYaxis()->SetNdivisions(515);

        TLine *l_;
        l_ = new TLine(ratio->GetXaxis()->GetXmin(),1,ratio->GetXaxis()->GetXmax(),1);
        l_->Draw("same");
        l_->SetLineStyle(1);
        l_->SetLineColor(kRed);

        for( int ii=0; ii<5; ii++ )
        {
                grid_vertical.DrawLine(totoalNbin/5 * (ii+1) + 0.5, ratio->GetMinimum(), totoalNbin/5 * (ii+1) + 0.5, ratio->GetMaximum() );
        }


        c1->cd();

        TPad *pad3 = new TPad("pad3","pad3",0,0.,1,0.33);
        setPadMargins(pad3);
        pad3->SetBottomMargin(0.3);
        pad3->SetTicks(1);
        pad3->SetGrid(0,1);
        pad3->SetLogx();
        pad3->Draw();
        pad3->cd();

        TGraphErrors *grUnfolded = new TGraphErrors(5, meanmass_data, meanpt_data, meanmasserr_data, meanpterr_data);
        TGraphErrors *grMC = new TGraphErrors(5, meanmass_mc, meanpt_mc, meanmasserr_mc, meanpterr_mc);

        grUnfolded->Draw("ape");
        grUnfolded->GetXaxis()->SetMoreLogLabels();
        grUnfolded->SetTitle("Average \\ {p_{T}^{\\ell\\ell}} vs M{\\ell\\ell}");
        grUnfolded->SetMarkerStyle(20);
        grUnfolded->SetMarkerSize(.9);
        grUnfolded->SetLineColor(kBlack);
        grUnfolded->SetMinimum(10.);
        grUnfolded->SetMaximum(30.);
	grUnfolded->GetYaxis()->SetNdivisions(505);
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

        CMS_lumi( pad1, 4, 0 );
        c1->cd();
        c1->SaveAs(outpdf);
	
	delete grUnfolded; 
	delete grMC; 
        delete hmeans;
        delete hmeansMC;
        delete l_; 
        delete leg_;
        delete pad1;
        delete pad2;
        delete pad3;
        delete c1;
}

void drawMassRatio(TString outpdf, TUnfoldDensityV17* unfold, TFile *filein){

        gROOT->SetBatch();

        setTDRStyle();
        writeExtraText = true;       // if extra text
        extraText  = "Preliminary";

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
  	pad1->SetTopMargin(0.1);
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

        TLine grid_;
        grid_.SetLineColor(kRed);
        grid_.SetLineStyle(kSolid);
        for( int ii=0; ii<histMCTruth->GetNbinsX(); ii++ )
        {
                Int_t i_bin = ii+1;
                Double_t binEdge = hunfolded->GetBinLowEdge(i_bin);
                grid_.DrawLine(binEdge, 100, binEdge, histMCTruth->GetBinContent(ii+1) );
        }

  	TLegend* leg_ = new TLegend(0.7, 0.60, 0.95, 0.9,"","brNDC");
  	leg_->SetTextSize(0.06);
  	leg_->SetFillStyle(0);
  	leg_->SetBorderSize(0);
  	leg_->AddEntry(hunfolded, "Unfolded data", "p");
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
	
        CMS_lumi( pad1, 4, 0 );
        c1->cd();
        c1->SaveAs(outpdf);

        delete l_; 
        delete leg_;
        delete pad1;
        delete pad2;
        delete c1;
}

void drawPtReco(TString outpdf, TString postfix, TFile *fdata, TFile *fDYsig, TFile *fDYbkg, TFile *fTTbar, TFile *fVV, TFile *fWjets, TString channel){

	setTDRStyle();
	writeExtraText = true;       // if extra text
	extraText  = "Preliminary";

	TString hname = "hPtRec" + postfix;
  	Int_t linecolorZ   = kOrange-3;
  	Int_t fillcolorZ   = kOrange-2;
   	Int_t linecolorEWK = kOrange+10;
   	Int_t fillcolorEWK = kOrange+7;
   	Int_t linecolorTop = kGreen+2;
   	Int_t fillcolorTop = kGreen-5;

	TUnfoldBinningV17 *ptbin = ptBinning_rec();

	TH1* hdata = (TH1*)fdata->Get(hname);
  	TH1* hdataNoUO = ptbin->ExtractHistogram("hdata", hdata, 0, kFALSE, "*[UO]");

	TH1* hdysig = (TH1*)fDYsig->Get(hname);
  	TH1* hdysigNoUO = ptbin->ExtractHistogram("hdysig", hdysig, 0, kFALSE, "*[UO]");

	TH1* hdybkg = (TH1*)fDYbkg->Get(hname);
  	TH1* hdybkgNoUO = ptbin->ExtractHistogram("hdybkg", hdybkg, 0, kFALSE, "*[UO]");

	TH1* httbar = (TH1*)fTTbar->Get(hname);
  	TH1* httbarNoUO = ptbin->ExtractHistogram("httbar", httbar, 0, kFALSE, "*[UO]");

	TH1* hvv = (TH1*)fVV->Get(hname);
  	TH1* hvvNoUO = ptbin->ExtractHistogram("hvv", hvv, 0, kFALSE, "*[UO]");

	TH1* hwjets = (TH1*)fWjets->Get(hname);
  	TH1* hwjetsNoUO = ptbin->ExtractHistogram("hwjets", hwjets, 0, kFALSE, "*[UO]");

   	THStack *hsMCs = new THStack("hsMCs","hsMCs");
        hdysigNoUO->SetLineColor(linecolorZ);
        hdysigNoUO->SetFillColor(fillcolorZ);
        hdybkgNoUO->SetLineColor(linecolorZ+1);
        hdybkgNoUO->SetFillColor(fillcolorZ+1);
        httbarNoUO->SetLineColor(linecolorTop);
        httbarNoUO->SetFillColor(fillcolorTop);
        hvvNoUO->SetLineColor(linecolorEWK);
        hvvNoUO->SetFillColor(fillcolorEWK);
        hwjetsNoUO->SetLineColor(linecolorEWK-1);
        hwjetsNoUO->SetFillColor(fillcolorEWK-1);

	hsMCs->Add(hwjetsNoUO);
	hsMCs->Add(hvvNoUO);
	hsMCs->Add(httbarNoUO);
	hsMCs->Add(hdybkgNoUO);
	hsMCs->Add(hdysigNoUO);

        TCanvas* c1=new TCanvas("c1", "c1", 50, 50, 700*1.5, 1000*1.5);
        c1->cd();
        gStyle->SetOptStat(0);

        TPad *pad1 = new TPad("pad1","pad1",0,0.45,1,1);
        setPadMargins(pad1);
        pad1->SetTopMargin(0.1);
        pad1->SetBottomMargin(0.01);
        pad1->SetTicks(1);
        pad1->SetLogy();
        pad1->Draw();
        pad1->cd();
        
        hdataNoUO->SetTitle("");
        hdataNoUO->Draw("p9histe");
        hdataNoUO->SetMarkerStyle(20);
        hdataNoUO->SetMarkerSize(.7);
        hdataNoUO->SetLineColor(kBlack);
        setYaxisHist(hdataNoUO);
	hsMCs->Draw("hist same");
        hdataNoUO->GetXaxis()->SetLabelSize(0.);
        hdataNoUO->GetYaxis()->SetTitle("Number of events per bin");
        //histMCTruth_pt->Draw("histsames");
        //histMCTruth_pt->SetLineColor(2);
        hdataNoUO->GetYaxis()->SetRangeUser(1.,hdataNoUO->GetMaximum() * 1e2);
        hdataNoUO->Draw("p9histsamee");

        int totoalNbin = hdataNoUO->GetXaxis()->GetNbins();

        TLine grid_vertical;
        grid_vertical.SetLineColor(kBlack);
        grid_vertical.SetLineStyle(kSolid);
        for( int ii=0; ii<5; ii++ )
        {
                grid_vertical.DrawLine(totoalNbin/5 * (ii+1) + 0.5, hdataNoUO->GetMinimum(), totoalNbin/5 * (ii+1) + 0.5, hdataNoUO->GetMaximum() );
        }

	TLegend *fLeg = new TLegend(0.65,0.55,0.95,0.88);
        fLeg->SetTextSize(0.05);
        fLeg->SetFillStyle(0);
        fLeg->SetBorderSize(0);
 	fLeg->AddEntry(hdataNoUO, "Data", "pe");
 	if( channel.CompareTo("electron") == 0 ) fLeg->AddEntry(hdysigNoUO, "Z #rightarrow ee", "F");
 	if( channel.CompareTo("muon") == 0 )     fLeg->AddEntry(hdysigNoUO, "Z #rightarrow #mu#mu", "F");
 	fLeg->AddEntry(hdybkgNoUO, "Z #rightarrow #tau#tau", "F");
 	fLeg->AddEntry(httbarNoUO, "ttbar", "F");
 	fLeg->AddEntry(hvvNoUO, "VV", "F");
 	fLeg->AddEntry(hwjetsNoUO, "Wjets", "F");
	fLeg->Draw();

        c1->cd();

        TPad *pad2 = new TPad("pad2","pad2",0,0.2,1,0.45);
        setPadMargins(pad2);
        pad2->SetBottomMargin(0.2);
        pad2->SetTicks(1);
        pad2->SetGridy();
        pad2->Draw();
        pad2->cd();

	TH1F *ratio = (TH1F*)hdataNoUO->Clone("ratio");
	TH1F *hmcs = (TH1F*)hdysigNoUO->Clone("mcs");
	hmcs->Add(hdybkgNoUO);
	hmcs->Add(httbarNoUO);
	hmcs->Add(hvvNoUO);
	hmcs->Add(hwjetsNoUO);

        ratio->Divide(hmcs);
        ratio->Draw("pe");
        ratio->SetMarkerStyle(20);
        ratio->SetMarkerSize(.7);
        ratio->SetLineColor(kBlack);
        ratio->SetMinimum(0.6);
        ratio->SetMaximum(1.4);
        ratio->GetYaxis()->SetTitle("#frac{Data}{MC}");
        ratio->GetYaxis()->CenterTitle();
        setYaxisHist(ratio);
        setXaxisHist(ratio);
        ratio->GetYaxis()->SetNdivisions(515);

        TLine *l_;
        l_ = new TLine(ratio->GetXaxis()->GetXmin(),1,ratio->GetXaxis()->GetXmax(),1);
        l_->Draw("same");
        l_->SetLineStyle(1);
        l_->SetLineColor(kRed);

        for( int ii=0; ii<5; ii++ )
        {
                grid_vertical.DrawLine(totoalNbin/5 * (ii+1) + 0.5, ratio->GetMinimum(), totoalNbin/5 * (ii+1) + 0.5, ratio->GetMaximum() );
        }


        c1->cd();

        TPad *pad3 = new TPad("pad3","pad3",0,0.,1,0.2);
        setPadMargins(pad3);
        pad3->SetBottomMargin(0.2);
        pad3->SetTicks(1);
        pad3->SetGrid(0,1);
        pad3->Draw();
        pad3->cd();

        TH1F *hratiottbar = (TH1F*)httbarNoUO->Clone("ratiottbar");
        TH1F *hratiodybkg = (TH1F*)hdybkgNoUO->Clone("ratiodybkg");

        hratiottbar->Divide(hmcs);
        hratiodybkg->Divide(hmcs);
        hratiottbar->Draw("hist");
        hratiottbar->SetMinimum(0.0);
        hratiottbar->SetMaximum(0.65);
        setYaxisHist(hratiottbar);
        setXaxisHist(hratiottbar);
	hratiottbar->SetFillColorAlpha(fillcolorTop, 0.35);
        hratiodybkg->Draw("histsame");
	hratiodybkg->SetFillColorAlpha(fillcolorZ+1, 0.35);

        CMS_lumi( pad1, 4, 0 );
        c1->cd();
        c1->SaveAs(outpdf);

        delete hdata;
        delete hdataNoUO;
        delete hdysig;
        delete hdysigNoUO;
        delete hsMCs;
        delete c1;
        delete ptbin;
}

void drawMassReco(TString outpdf, TString postfix, TFile *fdata, TFile *fDYsig, TFile *fDYbkg, TFile *fTTbar, TFile *fVV, TFile *fWjets, TString channel){

        setTDRStyle();
        writeExtraText = true;       // if extra text
        extraText  = "Preliminary";

	TString hname = "hMassRec" + postfix;
  	Int_t linecolorZ   = kOrange-3;
  	Int_t fillcolorZ   = kOrange-2;
   	Int_t linecolorEWK = kOrange+10;
   	Int_t fillcolorEWK = kOrange+7;
   	Int_t linecolorTop = kGreen+2;
   	Int_t fillcolorTop = kGreen-5;

	TUnfoldBinningV17 *massbin = massBinning_rec();

	TH1* hdata = (TH1*)fdata->Get(hname);
  	TH1* hdataNoUO = massbin->ExtractHistogram("hdata", hdata, 0, kTRUE, "*[UO]");

	TH1* hdysig = (TH1*)fDYsig->Get(hname);
  	TH1* hdysigNoUO = massbin->ExtractHistogram("hdysig", hdysig, 0, kTRUE, "*[UO]");

	TH1* hdybkg = (TH1*)fDYbkg->Get(hname);
  	TH1* hdybkgNoUO = massbin->ExtractHistogram("hdybkg", hdybkg, 0, kTRUE, "*[UO]");

	TH1* httbar = (TH1*)fTTbar->Get(hname);
  	TH1* httbarNoUO = massbin->ExtractHistogram("httbar", httbar, 0, kTRUE, "*[UO]");

	TH1* hvv = (TH1*)fVV->Get(hname);
  	TH1* hvvNoUO = massbin->ExtractHistogram("hvv", hvv, 0, kTRUE, "*[UO]");

	TH1* hwjets = (TH1*)fWjets->Get(hname);
  	TH1* hwjetsNoUO = massbin->ExtractHistogram("hwjets", hwjets, 0, kTRUE, "*[UO]");

   	THStack *hsMCs = new THStack("hsMCs","hsMCs");
        hdysigNoUO->SetLineColor(linecolorZ);
        hdysigNoUO->SetFillColor(fillcolorZ);
        hdybkgNoUO->SetLineColor(linecolorZ+1);
        hdybkgNoUO->SetFillColor(fillcolorZ+1);
        httbarNoUO->SetLineColor(linecolorTop);
        httbarNoUO->SetFillColor(fillcolorTop);
        hvvNoUO->SetLineColor(linecolorEWK);
        hvvNoUO->SetFillColor(fillcolorEWK);
        hwjetsNoUO->SetLineColor(linecolorEWK-1);
        hwjetsNoUO->SetFillColor(fillcolorEWK-1);

	hsMCs->Add(hwjetsNoUO);
	hsMCs->Add(hvvNoUO);
	hsMCs->Add(httbarNoUO);
	hsMCs->Add(hdybkgNoUO);
	hsMCs->Add(hdysigNoUO);

        TCanvas* c1=new TCanvas("c1", "c1", 50, 50, 700*1.5, 1000*1.5);
        c1->cd();
        gStyle->SetOptStat(0);

        TPad *pad1 = new TPad("pad1","pad1",0,0.45,1,1);
        setPadMargins(pad1);
        pad1->SetBottomMargin(0.01);
        pad1->SetTopMargin(0.1);
        pad1->SetTicks(1);
        pad1->SetLogy();
        pad1->Draw();
        pad1->cd();
        
        hdataNoUO->SetTitle("");
        hdataNoUO->Draw("p9histe");
        hdataNoUO->SetMarkerStyle(20);
        hdataNoUO->SetMarkerSize(.7);
        hdataNoUO->SetLineColor(kBlack);
        setYaxisHist(hdataNoUO);
	hsMCs->Draw("hist same");
        hdataNoUO->GetXaxis()->SetLabelSize(0.);
        hdataNoUO->GetYaxis()->SetTitle("Number of events per bin");
        //histMCTruth_mass->Draw("histsames");
        //histMCTruth_mass->SetLineColor(2);
        hdataNoUO->GetYaxis()->SetRangeUser(1.,hdataNoUO->GetMaximum() * 1e2);
        hdataNoUO->Draw("p9histsamee");
	
        TLegend *fLeg = new TLegend(0.65,0.55,0.95,0.88);
        fLeg->SetTextSize(0.05);
        fLeg->SetFillStyle(0);
        fLeg->SetBorderSize(0);
 	fLeg->AddEntry(hdataNoUO, "Data", "pe");
 	if( channel.CompareTo("electron") == 0 ) fLeg->AddEntry(hdysigNoUO, "Z #rightarrow ee", "F");
 	if( channel.CompareTo("muon") == 0 )     fLeg->AddEntry(hdysigNoUO, "Z #rightarrow #mu#mu", "F");
 	fLeg->AddEntry(hdybkgNoUO, "Z #rightarrow #tau#tau", "F");
 	fLeg->AddEntry(httbarNoUO, "ttbar", "F");
 	fLeg->AddEntry(hvvNoUO, "VV", "F");
 	fLeg->AddEntry(hwjetsNoUO, "Wjets", "F");
	fLeg->Draw();

        c1->cd();

        TPad *pad2 = new TPad("pad2","pad2",0,0.2,1,0.45);
        setPadMargins(pad2);
        pad2->SetBottomMargin(0.2);
        pad2->SetTicks(1);
        pad2->SetGridy();
        pad2->Draw();
        pad2->cd();

	TH1F *ratio = (TH1F*)hdataNoUO->Clone("ratio");
	TH1F *hmcs = (TH1F*)hdysigNoUO->Clone("mcs");
	hmcs->Add(hdybkgNoUO);
	hmcs->Add(httbarNoUO);
	hmcs->Add(hvvNoUO);
	hmcs->Add(hwjetsNoUO);

        ratio->Divide(hmcs);
        ratio->Draw("pe");
        ratio->SetMarkerStyle(20);
        ratio->SetMarkerSize(.7);
        ratio->SetLineColor(kBlack);
        ratio->SetMinimum(0.8);
        ratio->SetMaximum(1.2);
        ratio->GetYaxis()->SetTitle("#frac{Data}{MC}");
        ratio->GetYaxis()->CenterTitle();
        setYaxisHist(ratio);
        setXaxisHist(ratio);
        ratio->GetYaxis()->SetNdivisions(515);

        TLine *l_;
        l_ = new TLine(ratio->GetXaxis()->GetXmin(),1,ratio->GetXaxis()->GetXmax(),1);
        l_->Draw("same");
        l_->SetLineStyle(1);
        l_->SetLineColor(kRed);

        c1->cd();

        TPad *pad3 = new TPad("pad3","pad3",0,0.,1,0.2);
        setPadMargins(pad3);
        pad3->SetBottomMargin(0.2);
        pad3->SetTicks(1);
        pad3->SetGrid(0,1);
        pad3->Draw();
        pad3->cd();

        TH1F *hratiottbar = (TH1F*)httbarNoUO->Clone("ratiottbar");
        TH1F *hratiodybkg = (TH1F*)hdybkgNoUO->Clone("ratiodybkg");

        hratiottbar->Divide(hmcs);
        hratiodybkg->Divide(hmcs);
        hratiottbar->Draw("hist");
        hratiottbar->SetMinimum(0.0);
        hratiottbar->SetMaximum(0.65);
        setYaxisHist(hratiottbar);
        setXaxisHist(hratiottbar);
	hratiottbar->SetFillColorAlpha(fillcolorTop, 0.35);
        hratiodybkg->Draw("histsame");
	hratiodybkg->SetFillColorAlpha(fillcolorZ+1, 0.35);

        CMS_lumi( pad1, 4, 0 );
        c1->cd();
        c1->SaveAs(outpdf);

        delete hdata;
        delete hdataNoUO;
        delete hdysig;
        delete hdysigNoUO;
        delete hsMCs;
        delete c1;
        delete massbin;
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

        gPad->SetLogz();  
        gStyle->SetPalette(1);

        histProbability->GetZaxis()->SetRangeUser(0., 1.);
        histProbability->Draw("colz");

        TLine grid_;
        //grid_.SetLineColor(kGray+2);
        grid_.SetLineColorAlpha(kGray+2, 0.35);;
        grid_.SetLineStyle(kSolid);
        for( int ii=0; ii<histProbability->GetXaxis()->GetNbins(); ii++ )
        {
                Int_t i_bin = ii+1;
                Double_t binEdge = histProbability->GetXaxis()->GetBinUpEdge(i_bin);
                grid_.DrawLine(binEdge, histProbability->GetYaxis()->GetBinUpEdge(0), binEdge, histProbability->GetYaxis()->GetBinUpEdge(histProbability->GetYaxis()->GetNbins()) );
        }

        for( int ii=0; ii<histProbability->GetYaxis()->GetNbins(); ii++ )
        {
                Int_t i_bin = ii+1;
                Double_t binEdge = histProbability->GetYaxis()->GetBinUpEdge(i_bin);
                grid_.DrawLine(histProbability->GetXaxis()->GetBinUpEdge(0), binEdge, histProbability->GetXaxis()->GetBinUpEdge(histProbability->GetXaxis()->GetNbins()),binEdge );
        }

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
