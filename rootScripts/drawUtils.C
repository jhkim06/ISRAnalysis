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

void setTGraphAxis(TGraphErrors* data, TGraphErrors* sys, TGraphErrors* mc, Double_t x, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax){


        TString mass;
        mass.Form("%0.2f", x);

        data->Draw("ape");
        sys->SetFillColor(1);
        sys->SetFillStyle(3004);
        sys->Draw("same2");
        data->Draw("pe");
        mc->Draw("samepe");
        data->SetTitle("Average \\ {p_{T}^{\\ell\\ell}} vs M{\\ell\\ell}");
        data->SetMarkerStyle(20);
        data->SetMarkerSize(1.);
        data->SetLineColor(kBlack);
        data->SetMinimum(ymin);
        data->SetMaximum(ymax);
        data->GetXaxis()->SetLimits(xmin, xmax);
//	data->GetXaxis()->SetBinLabel(data->GetXaxis()->FindBin(x), mass);
        data->GetXaxis()->SetNdivisions(501);
        data->GetYaxis()->SetNdivisions(505);
        data->GetYaxis()->SetLabelOffset(-0.2);
        setYaxisGraph(data);
        setXaxisGraph(data);
        //data->GetYaxis()->SetTitle("Average \\ {p_{T}^{\\ell\\ell}}");

        mc->SetMarkerStyle(20);
        mc->SetMarkerSize(1.);
        mc->SetLineColor(kRed);
        mc->SetMarkerColor(kRed);

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

void getAveragesMass(vector<Double_t> &mean, vector<Double_t> &err, TH1* hmass){

	Double_t massBins[6] = {40., 60., 80., 100., 200., 350.};
	for(int ibin = 0; ibin < 5; ibin++){
		hmass->GetXaxis()->SetRange(hmass->GetXaxis()->FindBin(massBins[ibin]+0.01),hmass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
		mean.push_back(hmass->GetMean());
		err.push_back(hmass->GetMeanError());
	}

}

// TODO make a function to get systematic error on the average pt for each systematic source
void getAveragesPt(vector<Double_t> &mean, vector<Double_t> &err, TUnfoldDensityV17* unfold_pt, bool isData){

	int nMassBin = 5;

        for(int i = 0; i < nMassBin; i++){
           TString ibinMass;
           ibinMass.Form("%d", i);

           TH1* hpt_temp;  

           if(isData) hpt_temp=unfold_pt->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
           else       hpt_temp=unfold_pt->GetBias("histMCTruth_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);

           mean.push_back(hpt_temp->GetMean());
           err.push_back(hpt_temp->GetMeanError());
        
           delete hpt_temp;
        }           
}

void getAveragesSysMass(vector<Double_t> &err, TString sysName, int sysSize, TUnfoldDensityV17* unfold_mass){

        Double_t massBins[6] = {40., 60., 80., 100., 200., 350.};

        TH1* hmass_temp;
        hmass_temp=unfold_mass->GetOutput("hunfolded_temp",0,0,"*[UO]",kTRUE);

        for(int ibin = 0; ibin < 5; ibin++){
                hmass_temp->GetXaxis()->SetRange(hmass_temp->GetXaxis()->FindBin(massBins[ibin]+0.01),hmass_temp->GetXaxis()->FindBin(massBins[ibin+1]-0.01));

           	Double_t defaultMean = hmass_temp->GetMean();
		std::cout << "mass " << ibin << " mean : " << defaultMean << std::endl;

           	Double_t err_ = -999.;
           	for(int i = 0; i < sysSize; i++){
           	        if(i==5 || i==7) continue;

           	        TString isys;
           	        isys.Form("%d", i);

           	        TH1* hsysmass_temp;

           	        hsysmass_temp=(TH1*)unfold_mass->GetDeltaSysSource(sysName+"_"+isys, sysName+"_"+isys, sysName+"_"+isys, "Gen", "*[UO]", kTRUE);
           	        hsysmass_temp->Add(hmass_temp);
			hsysmass_temp->GetXaxis()->SetRange(hsysmass_temp->GetXaxis()->FindBin(massBins[ibin]+0.01),hsysmass_temp->GetXaxis()->FindBin(massBins[ibin+1]-0.01));

           	        Double_t temp_err =  fabs( defaultMean - hsysmass_temp->GetMean() );
           	        if( err_ < temp_err ){
           	          err_ = temp_err;
           	        }

           	       delete hsysmass_temp;
           	}

                err.push_back(err_);
        }

	delete hmass_temp;

}

void getAveragesSysPt(vector<Double_t> &err, TString sysName, int sysSize, TUnfoldDensityV17* unfold_pt){

        int nMassBin = 5;

        for(int i = 0; i < nMassBin; i++){
           TString ibinMass;
           ibinMass.Form("%d", i);

           TH1* hpt_temp;

           hpt_temp=unfold_pt->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);

	   Double_t defaultMean = hpt_temp->GetMean();

	   Double_t err_ = -999.;
           for(int i = 0; i < sysSize; i++){
                   if(i==5 || i==7) continue;

               	   TString isys;
                   isys.Form("%d", i);

           	   TH1* hsyspt_temp;

           	   hsyspt_temp=(TH1*)unfold_pt->GetDeltaSysSource(sysName+"_"+isys, sysName+"_"+isys, sysName+"_"+isys, "Gen", "pt[UO];mass[UOC"+ibinMass+"]", kTRUE);
		   hsyspt_temp->Add(hpt_temp);

                   Double_t temp_err =  fabs( defaultMean - hsyspt_temp->GetMean() );
		   if( err_ < temp_err ){
		     err_ = temp_err;
 		   }

		  delete hsyspt_temp;
	   }
	  
           std::cout << "err: " << err_ << std::endl;
	   err.push_back(err_);

           delete hpt_temp;
        }
}

void getRatioSys(TUnfoldDensityV17* unfold_pt, TString sysName, int sysSize, TH1 *sysRatio){

	vector<TH1*> hDeltas;
	vector<TH1*> hSys;

        TH1* hunfolded_pt = unfold_pt->GetOutput("hunfolded_pt_temp",0,0,"*[UO]",kFALSE);
        TH1* ratio=(TH1*)hunfolded_pt->Clone("ratio");
        TH1 *histMCTruth_pt=unfold_pt->GetBias("histMC_pt",0,0,"*[UO]",kFALSE);
        ratio->Divide(histMCTruth_pt);

	for(int i = 0; i < sysSize; i++){
           	TString isys;
           	isys.Form("%d", i);
		hDeltas.push_back((TH1*)unfold_pt->GetDeltaSysSource(sysName+"_"+isys, sysName+"_"+isys, sysName+"_"+isys, "Gen", "*[UO]", kFALSE));
		hSys.push_back(unfold_pt->GetOutput("hunfolded_pt"+sysName+"_"+isys,0,0,"*[UO]",kFALSE));

		hSys.at(i)->Add(hDeltas.at(i));
		hSys.at(i)->Divide(histMCTruth_pt);
	}

        for(int ibin = 1; ibin<histMCTruth_pt->GetNbinsX()+1;ibin++){

           sysRatio->SetBinContent(ibin, ratio->GetBinContent(ibin));

           Double_t err = -999.;
	   for(int i = 0; i < sysSize; i++){
		   if(i==5 || i==7) continue;
		   // get "envelope"
	   	   Double_t temp_err =  fabs(hSys.at(i)->GetBinContent(ibin) - ratio->GetBinContent(ibin));
	           if( temp_err > err){
	   		  err = temp_err;
	           }
	   }
	   sysRatio->SetBinError(ibin, err);
        }

	for(int i = 0; i < sysSize; i++){
	   delete hDeltas.at(i);
           delete hSys.at(i);

	}
	delete hunfolded_pt;
	hDeltas.clear();
	hDeltas.shrink_to_fit();

	hSys.clear();
	hSys.shrink_to_fit();
}

void drawRatio(TString outpdf, TUnfoldDensityV17* unfold_pt, TUnfoldDensityV17* unfold_mass, TFile *filein){

        gROOT->SetBatch();

        setTDRStyle();
        writeExtraText = true;       // if extra text
        extraText  = "Preliminary";

	// pt distribution for 5 mass bins in one histogram
	TH1* hunfolded_pt = unfold_pt->GetOutput("hunfolded_pt",0,0,"*[UO]",kFALSE);
	TH1* ratio=(TH1*)hunfolded_pt->Clone("ratio");

	// get systematic for ratio plot
        TH1* sysErrRatio = (TH1*)hunfolded_pt->Clone("hsysErrRatio");
	//getRatioSys(unfold_pt, "AlphaS", (int)2, sysErrRatio);	
	getRatioSys(unfold_pt, "Scale", (int)9, sysErrRatio);	

        TH1 *histMCTruth_pt=unfold_pt->GetBias("histMCTruth_pt",0,0,"*[UO]",kFALSE);
        ratio->Divide(histMCTruth_pt);

	// mass distribution
        TH1* hunfolded_mass = unfold_mass->GetOutput("hunfolded_mass",0,0,"*[UO]",kTRUE);
        TH1 *histMCTruth_mass=unfold_mass->GetBias("histMCTruth_mass",0,0,"*[UO]",kTRUE);

        vector<Double_t> meanmass_data, meanmasserr_data;
        vector<Double_t> meanmass_mc, meanmasserr_mc;

	vector<Double_t> meanpt_data, meanpterr_data;
	vector<Double_t> meanpt_mc, meanpterr_mc;

	// get average pt for each mass bin
	getAveragesPt(meanpt_data, meanpterr_data, unfold_pt, true);
	getAveragesPt(meanpt_mc, meanpterr_mc, unfold_pt, false);

	vector<Double_t> ScaleErr;
	getAveragesSysPt(ScaleErr, "Scale", (int)9, unfold_pt);

	// get average mass
	//
	getAveragesMass(meanmass_data, meanmasserr_data, hunfolded_mass);
	getAveragesMass(meanmass_mc, meanmasserr_mc, histMCTruth_mass);

	vector<Double_t> ScaleMassErr;
	getAveragesSysMass(ScaleMassErr, "Scale", (int)9, unfold_mass);


  	TCanvas* c1=new TCanvas("c1", "c1", 50, 50, 700*1.5, 1000*1.5);
  	c1->cd();
  	gStyle->SetOptStat(0);

  	TPad *pad1 = new TPad("pad1","pad1",0,0.6,1,1);
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

  	TPad *pad2 = new TPad("pad2","pad2",0,0.33,1,0.6);
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

	// draw systematic error
  	sysErrRatio->Draw("E2same");
  	sysErrRatio->SetMarkerSize(0);
  	sysErrRatio->SetFillColorAlpha(kBlack,0.3);

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

	Double_t initial_padx = 0.05;
	for(int imass = 0; imass < 5; imass++){

        	TString imass_;
        	imass_.Form("%d", imass);

        	c1->cd();

        	TPad *pad = new TPad("mass_" + imass_,"mass_" + imass_, initial_padx + (0.15 * imass) + (0.05 * imass), 0., initial_padx + (0.15 * (imass + 1)) + (0.05 * imass), 0.33);
        	setPadMargins(pad);
        	pad->SetLeftMargin(0.1);
        	pad->SetRightMargin(0.1);
        	pad->SetBottomMargin(0.3);
        	pad->SetTicks(1);
        	pad->SetGrid(0,1);
        	//pad->SetLogx();
        	pad->Draw();
        	pad->cd();

		// FIXME how to handle below memory allocation
        	TGraphErrors *grUnfolded = new TGraphErrors(1, &meanmass_data[imass], &meanpt_data[imass], &meanmasserr_data[imass], &meanpterr_data[imass]);
        	TGraphErrors *sysData = new TGraphErrors(1, &meanmass_data[imass], &meanpt_data[imass], &ScaleMassErr[imass], &ScaleErr[imass]);
        	TGraphErrors *grMC = new TGraphErrors(1, &meanmass_mc[imass], &meanpt_mc[imass], &meanmasserr_mc[imass], &meanpterr_mc[imass]);

        	Double_t xmin = meanmass_data[imass] > meanmass_mc[imass] ? meanmass_mc[imass] * 0.99: meanmass_data[imass] * 0.99,
        	         xmax = meanmass_data[imass] < meanmass_mc[imass] ? meanmass_mc[imass] * 1.01: meanmass_data[imass] * 1.01,
        	         ymin = meanpt_data[imass] > meanpt_mc[imass] ? meanpt_mc[imass] * 0.95: meanpt_data[imass] * 0.95,
        	         ymax = meanpt_data[imass] < meanpt_mc[imass] ? meanpt_mc[imass] * 1.05: meanpt_data[imass] * 1.05;
        	setTGraphAxis(grUnfolded, sysData, grMC, meanmass_data[imass], xmin, xmax, ymin, ymin + 2.);



	}

        CMS_lumi( pad1, 4, 0 );
        c1->cd();
        c1->SaveAs(outpdf);
	
	//delete grUnfolded; 
	//delete grMC; 
        delete l_; 
        delete leg_;
        delete pad1;
        delete pad2;
        //delete pad3;
        delete c1;
}

void drawMassRatio(TString outpdf, TUnfoldDensityV17* unfold, TFile *filein){

        gROOT->SetBatch();

        setTDRStyle();
        writeExtraText = true;       // if extra text
        extraText  = "Preliminary";

	TH1* hunfolded = unfold->GetOutput("hunfolded",0,0,"*[UO]",kTRUE);
	TH1* ratio=(TH1*)hunfolded->Clone("ratio");

        TH1* ratioUp=(TH1*)hunfolded->Clone("ratioUp");     // default unfoled data
        TH1* ratioDown=(TH1*)hunfolded->Clone("ratioDown"); //
        TH1* sysErrRatio = (TH1*)hunfolded->Clone("hsysErrRatio");

        ratioUp->  Add((TH1*)unfold->GetDeltaSysSource("AlphaS_0","SysAlphaS_0","SysAlphaS_0","Gen","*[UO]",kTRUE)); // add systematic delta 
        ratioDown->Add((TH1*)unfold->GetDeltaSysSource("AlphaS_1","SysAlphaS_1","SysAlphaS_1","Gen","*[UO]",kTRUE)); // add systematic delta 

        TH1 *histMCTruth=unfold->GetBias("histMCTruth",0,0,"*[UO]",kTRUE);
        ratio->Divide(histMCTruth);

        ratioUp->Divide(histMCTruth);
        ratioDown->Divide(histMCTruth);

        for(int ibin = 1; ibin<ratioUp->GetNbinsX()+1;ibin++){

           Double_t errUp = 0., errDown = 0.;
           errUp = fabs(ratioUp->GetBinContent(ibin) - ratio->GetBinContent(ibin));
           errDown = fabs(ratioDown->GetBinContent(ibin) - ratio->GetBinContent(ibin));

           sysErrRatio->SetBinContent(ibin, ratio->GetBinContent(ibin));
           sysErrRatio->SetBinError(ibin, errUp >= errDown ? errUp : errDown);
           if(sysErrRatio->GetBinError(ibin)==0) sysErrRatio->SetBinError(ibin, 10e-5);
        }


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

        // draw systematic error
        sysErrRatio->Draw("E2same");
        sysErrRatio->SetMarkerSize(0);
        sysErrRatio->SetFillColorAlpha(kBlack,0.3);

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

        for( int ii=0; ii<5; ii++ )
        {
                grid_vertical.DrawLine(totoalNbin/5 * (ii+1) + 0.5, hratiottbar->GetMinimum(), totoalNbin/5 * (ii+1) + 0.5, hratiottbar->GetMaximum() );
        }


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
