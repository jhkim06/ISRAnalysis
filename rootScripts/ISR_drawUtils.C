#include "ISR_drawUtils.h"
#include "tdrstyle.C"
#include "CMS_lumi.C"

void setYaxisHist(TH1* hist){

        hist->SetTitle("");
        hist->GetYaxis()->SetLabelFont(63);
        hist->GetYaxis()->SetLabelSize(30);
        hist->GetYaxis()->SetTitleFont(63);
        hist->GetYaxis()->SetTitleSize(40);
        hist->GetYaxis()->SetTitleOffset(1.7);
}

void setYaxisGraph(TGraphErrors* gr){

        gr->SetTitle("");
        gr->GetYaxis()->SetLabelFont(63);
        gr->GetYaxis()->SetLabelSize(30);
        gr->GetYaxis()->SetTitleFont(63);
        gr->GetYaxis()->SetTitleSize(40);
        gr->GetYaxis()->SetTitleOffset(1.7);
}

void setXaxisHist(TH1* hist){

        hist->GetXaxis()->SetLabelFont(63);
        hist->GetXaxis()->SetLabelSize(25);
        hist->GetXaxis()->SetTitleOffset(2.5);
        hist->GetXaxis()->SetLabelOffset(0.02);
        hist->GetXaxis()->SetTitleFont(63);
        hist->GetXaxis()->SetTitleSize(40);
}

void setXaxisGraph(TGraphErrors* gr){

        gr->GetXaxis()->SetLabelFont(63);
        gr->GetXaxis()->SetLabelSize(30);
        gr->GetXaxis()->SetTitleOffset(1.5);
        gr->GetXaxis()->SetLabelOffset(0.02);
        gr->GetXaxis()->SetTitleFont(63);
        gr->GetXaxis()->SetTitleSize(40);
}

void setTGraphAxis(TGraphErrors* data, Double_t x, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, TGraphErrors* mc,  TGraphErrors* sys, bool axis){


        TString mass;
        mass.Form("%0.2f", x);

        if(axis) data->Draw("ape");
	else data->Draw("samepe");
	//data->GetXaxis()->ChangeLabel(2, -45, 2);
	if(sys != nullptr) sys->SetFillColor(kGray+1);
        if(sys != nullptr) sys->SetFillStyle(3001);
        if(sys != nullptr) sys->Draw("same2");
        data->Draw("pe");
        if(mc != nullptr) mc->Draw("samepe");
        data->SetTitle("Average \\ {p_{T}^{\\ell\\ell}} vs M{\\ell\\ell}");
        if(axis) data->SetMarkerStyle(20);
 	else data->SetMarkerStyle(25);
        data->SetMarkerSize(1.);
        data->SetLineColor(kBlack);
        data->SetLineWidth(1);
        data->SetMinimum(ymin);
        data->SetMaximum(ymax);
        data->GetXaxis()->SetLimits(xmin, xmax);
        data->GetXaxis()->SetNdivisions(501);
        data->GetYaxis()->SetNdivisions(505);
        data->GetYaxis()->SetLabelOffset(-0.2);
        setYaxisGraph(data);
        setXaxisGraph(data);
	data->GetXaxis()->SetLabelSize(0);
        //data->GetYaxis()->SetTitle("Average \\ {p_{T}^{\\ell\\ell}}");

        if(mc != nullptr) mc->SetMarkerStyle(20);
        if(mc != nullptr) mc->SetMarkerSize(1.);
        if(mc != nullptr) mc->SetMarkerColor(kRed);
        if(mc != nullptr) mc->SetLineColor(kRed);

}

void setPadMargins(TPad* pad){

        pad->SetTopMargin(0.1);
        pad->SetLeftMargin(0.2);
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


void getAveragesMass(vector<Double_t> &mean, vector<Double_t> &err, TH1* hmass){

	//Double_t massBins[6] = {40., 60., 80., 100., 200., 350.};
	Double_t massBins[6] = {50., 65., 80., 100., 200., 350.};
	for(int ibin = 0; ibin < 5; ibin++){
		hmass->GetXaxis()->SetRange(hmass->GetXaxis()->FindBin(massBins[ibin]+0.01),hmass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
		mean.push_back(hmass->GetMean());
		err.push_back(hmass->GetMeanError());
	}

}

// TODO make a function to get systematic error on the average pt for each systematic source
void getAveragesPt(vector<Double_t> &mean, vector<Double_t> &err, TUnfoldDensity* unfold_pt, bool isData){

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

void getAveragesSysMass(vector<Double_t> &err, TString sysName, int sysSize, TUnfoldDensity* unfold_mass){

        //Double_t massBins[6] = {40., 60., 80., 100., 200., 350.};
        Double_t massBins[6] = {50., 65., 80., 100., 200., 350.};

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

           	        hsysmass_temp=(TH1*)unfold_mass->GetDeltaSysSource(sysName+"_"+isys, sysName+"_"+isys, sysName+"_"+isys, "Gen_Mass", "*[UO]", kTRUE);
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

void getAveragesSysPt(vector<Double_t> &err, TString sysName, int sysSize, TUnfoldDensity* unfold_pt){

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

           	   hsyspt_temp=(TH1*)unfold_pt->GetDeltaSysSource(sysName+"_"+isys, sysName+"_"+isys, sysName+"_"+isys, "Gen_Pt", "pt[UO];mass[UOC"+ibinMass+"]", kTRUE);
		   hsyspt_temp->Add(hpt_temp);

                   Double_t temp_err =  fabs( defaultMean - hsyspt_temp->GetMean() );
		   if( err_ < temp_err ){
		     err_ = temp_err;
 		   }

		  delete hsyspt_temp;
	   }
	  
	   err.push_back(err_);

           delete hpt_temp;
        }
}


void getRatioSysMass(TUnfoldDensity* unfold_mass, TString sysName, int sysSize, TH1 *sysRatio){

	vector<TH1*> hDeltas;
	vector<TH1*> hSys;

        TH1* hunfolded_mass = unfold_mass->GetOutput("hunfolded_mass_temp",0,0,"*[UO]",kTRUE);
        TH1* ratio=(TH1*)hunfolded_mass->Clone("ratio");
        TH1 *histMCTruth_mass=unfold_mass->GetBias("histMC_mass",0,0,"*[UO]",kTRUE);
        ratio->Divide(histMCTruth_mass);

	for(int i = 0; i < sysSize; i++){
           	TString isys;
           	isys.Form("%d", i);
		hDeltas.push_back((TH1*)unfold_mass->GetDeltaSysSource(sysName+"_"+isys, sysName+"_"+isys, sysName+"_"+isys, "Gen_Mass", "*[UO]", kTRUE));
		hSys.push_back(unfold_mass->GetOutput("hunfolded_mass"+sysName+"_"+isys,0,0,"*[UO]",kTRUE));

		hSys.at(i)->Add(hDeltas.at(i));
		hSys.at(i)->Divide(histMCTruth_mass);
	}

        for(int ibin = 1; ibin<histMCTruth_mass->GetNbinsX()+1;ibin++){

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
	   //sysRatio->SetBinError(ibin, err);
	   sysRatio->SetBinError(ibin, sqrt(pow(err,2) + pow(sysRatio->GetBinError(ibin),2)));
        }

	for(int i = 0; i < sysSize; i++){
	   delete hDeltas.at(i);
           delete hSys.at(i);

	}
	delete hunfolded_mass;
	hDeltas.clear();
	hDeltas.shrink_to_fit();

	hSys.clear();
	hSys.shrink_to_fit();
}

void getRatioSys(TUnfoldDensity* unfold_pt, TString sysName, int sysSize, TH1 *sysRatio){

	vector<TH1*> hDeltas;
	vector<TH1*> hSys;

        TH1* hunfolded_pt = unfold_pt->GetOutput("hunfolded_pt_temp",0,0,"*[UO]",kFALSE); // default unfolding result
        TH1* ratio=(TH1*)hunfolded_pt->Clone("ratio");
        TH1 *histMCTruth_pt=unfold_pt->GetBias("histMC_pt",0,0,"*[UO]",kFALSE);
        ratio->Divide(histMCTruth_pt); // default ratio 

	for(int i = 0; i < sysSize; i++){
           	TString isys;
           	isys.Form("%d", i);
		hDeltas.push_back((TH1*)unfold_pt->GetDeltaSysSource(sysName+"_"+isys, sysName+"_"+isys, sysName+"_"+isys, "Gen_Pt", "*[UO]", kFALSE));
		hSys.push_back(unfold_pt->GetOutput("hunfolded_pt"+sysName+"_"+isys,0,0,"*[UO]",kFALSE));// default unfolding result

		hSys.at(i)->Add(hDeltas.at(i));
		hSys.at(i)->Divide(histMCTruth_pt); // systematic / truth
	}

        for(int ibin = 1; ibin<histMCTruth_pt->GetNbinsX()+1;ibin++){

           //sysRatio->SetBinContent(ibin, ratio->GetBinContent(ibin));
           sysRatio->SetBinContent(ibin, 1.);

           Double_t err = -999.;
	   for(int i = 0; i < sysSize; i++){
		   if(i==5 || i==7) continue;
		   // get "envelope"
	   	   Double_t temp_err =  fabs(hSys.at(i)->GetBinContent(ibin) - ratio->GetBinContent(ibin));
	           if( temp_err > err){
	   		  err = temp_err;
	           }
	   }
	   sysRatio->SetBinError(ibin, sqrt(pow(err,2) + pow(sysRatio->GetBinError(ibin),2)));
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

void drawUnfoldedPtDistWithSys(TString outpdf, TUnfoldDensity* unfold_pt){

        gROOT->SetBatch();

        setTDRStyle();
        writeExtraText = true;       // if extra text
        extraText  = "Preliminary";

        for(int imass = 0; imass < 5; imass++){

           	TString ibinMass;
           	ibinMass.Form("%d", imass);

           	TH1* hpt_temp;
           	TH1* hpterr_temp;
        	TLegend* leg_ = new TLegend(0.2, 0.70, 0.95, 0.9,"","brNDC");
		leg_->SetNColumns(2);
        	leg_->SetTextSize(0.04);
        	leg_->SetFillStyle(0);
        	leg_->SetBorderSize(0);

           	hpt_temp=unfold_pt->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
           	hpterr_temp=unfold_pt->GetOutput("hunfolded_pterr_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);

        	TCanvas* c1=new TCanvas("c1", "c1", 50, 50, 700*1.5, 800*1.5);
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

        	hpt_temp->SetTitle("");
        	hpt_temp->Draw("p9histe");
        	hpt_temp->SetMarkerStyle(20);
        	hpt_temp->SetMarkerSize(.7);
        	hpt_temp->SetLineColor(kBlack);


                int sysSize = 2;
                TString sysName = "IsoSF";

        	for(int ibin = 1; ibin<hpterr_temp->GetNbinsX()+1;ibin++){

        	   Double_t err = -999.;
        	   for(int i = 0; i < sysSize; i++){
        	           if(i==5 || i==7) continue;

                           TString isys;
                           isys.Form("%d", i);

			   TH1F * hsyspt_temp_=((TH1F*)unfold_pt->GetDeltaSysSource(sysName+"_"+isys, sysName+"_"+isys+ibinMass, sysName+"_"+isys+ibinMass, "Gen_Pt", "pt[UO];mass[UOC"+ibinMass+"]", kTRUE));
			   hsyspt_temp_->Add(hpt_temp);
        	           // get "envelope"
        	           Double_t temp_err =  fabs(hpt_temp->GetBinContent(ibin) - hsyspt_temp_->GetBinContent(ibin));
        	           if( temp_err > err){
        	                  err = temp_err;
        	           }
			   delete hsyspt_temp_;
        	   }
        	   hpterr_temp->SetBinError(ibin, err);
        	}

	        hpterr_temp->Draw("E2same");
	        hpterr_temp->SetMarkerSize(0);
	        hpterr_temp->SetFillColorAlpha(kRed,0.3);

		c1->cd();

        	TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.6);
        	setPadMargins(pad2);
        	pad2->SetBottomMargin(0.2);
        	pad2->SetTicks(1);
        	pad2->SetGridy(1);
        	pad2->Draw();

           	TH1F *hsyspt_temp;
		TH1F *ratio;
           	for(int i = 0; i < sysSize; i++){
           	        if(i==5 || i==7) continue;

           	        TString isys;
           	        isys.Form("%d", i);


           	        hsyspt_temp=((TH1F*)unfold_pt->GetDeltaSysSource(sysName+"_"+isys, sysName+"_"+isys+ibinMass, sysName+"_"+isys+ibinMass, "Gen_Pt", "pt[UO];mass[UOC"+ibinMass+"]", kTRUE));
           	        hsyspt_temp->Add(hpt_temp);


			ratio= ((TH1F*)hsyspt_temp->Clone("ratio"+isys+ibinMass));
			ratio->Divide(hpt_temp);

			pad2->cd();
                	if(i==0)ratio->Draw("hist");
			else ratio->Draw("samehist");
                	ratio->SetLineColor(i+1);
                	ratio->SetMarkerColor(i+1);
                	ratio->SetMinimum(0.9);
                	ratio->SetMaximum(1.1);
                	ratio->SetTitle("");
                	setYaxisHist(ratio);
                	setXaxisHist(ratio);
                	ratio->GetYaxis()->SetNdivisions(515);

                	TString mean_;
                	mean_.Form("%.5f", hsyspt_temp->GetMean());


        		leg_->AddEntry(ratio, sysName+"_"+isys + " " + mean_, "l");
        		leg_->Draw();

           	}


		c1->cd();
        	c1->SaveAs(outpdf+"_"+ibinMass+".pdf");
		delete hpt_temp;
		delete hpterr_temp;
		delete pad1;
		delete c1;

	}
}

void initHistContent(TH1* hist, Double_t value){

        for(int ibin = 1; ibin<hist->GetNbinsX()+1;ibin++){
                hist->SetBinContent(ibin, value);
        }
}

void initHistError(TH1* hist, Double_t value){

        for(int ibin = 1; ibin<hist->GetNbinsX()+1;ibin++){
                hist->SetBinError(ibin, value);
        }
}

TUnfoldDensityV17* getSysMatrix(TString var, TString sysName, TString channel){

	TFile* filein;
	if( channel.CompareTo("electron") == 0 ) filein = new TFile("/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/output/2016/electron/DYtoEE.root");
	if( channel.CompareTo("muon") == 0 ) filein = new TFile("/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/output/2016/muon/DYtoEE.root");
        TH2* hmcGenRec = (TH2*)filein->Get("hmc" + var + "GenRec" + sysName);
	hmcGenRec->SetDirectory(0);
        TUnfoldBinning* binning_Rec = NULL;
        TUnfoldBinning* binning_Gen = NULL;

        if( var.CompareTo("Pt") == 0 ){

          binning_Rec = (TUnfoldBinning*)filein->Get("Rec_Pt");
          binning_Gen = (TUnfoldBinning*)filein->Get("Gen_Pt");
        }
        if( var.CompareTo("Mass") == 0 ){

          binning_Rec = (TUnfoldBinning*)filein->Get("Rec_Mass");
          binning_Gen = (TUnfoldBinning*)filein->Get("Gen_Mass");

        }

        TUnfoldDensityV17 *unfold;
        unfold = new TUnfoldDensityV17(hmcGenRec,
                                       TUnfold::kHistMapOutputHoriz,
                                       TUnfold::kRegModeNone, // fixed to use no regularisation temporary
                                       TUnfold::kEConstraintArea,
                                       TUnfoldDensityV17::kDensityModeBinWidth,
                                       binning_Gen,binning_Rec);

        TH1* hmcRec = (TH1*)filein->Get("h" + var + "Recnominal");
	unfold->SetInput(hmcRec, 1); // dummy input

	return unfold;
}

void sysMCErrRatio(TUnfoldDensity* unfold_pt, TH1* sysRatio, TString var, TString sysName, int size, TString channel){

	vector<TH1*> hMCSys;
	vector<TUnfoldDensityV17*> hMCSysM;

	TH1* hunfolded_pt = unfold_pt->GetOutput("hunfolded_pt_temp",0,0,"*[UO]",kFALSE);
	TH1* ratio=(TH1*)hunfolded_pt->Clone("ratio");
        TH1 *histMCTruth_pt=unfold_pt->GetBias("histMC_pt",0,0,"*[UO]",kFALSE);
        ratio->Divide(histMCTruth_pt); // default ratio 

	for(int i=0; i<size; i++){
                TString isys;
                isys.Form("%d", i);
		hMCSysM.push_back(getSysMatrix("Pt", sysName+"_"+isys, channel));
		hMCSys.push_back((TH1*)hunfolded_pt->Clone("ratio"+isys)); // default data
		hMCSys.at(i)->Divide(hMCSysM.at(i)->GetBias("histMC_pt"+isys,0,0,"*[UO]",kFALSE));
	}

        for(int ibin = 1; ibin<histMCTruth_pt->GetNbinsX()+1;ibin++){

           sysRatio->SetBinContent(ibin, ratio->GetBinContent(ibin));
           //sysRatio->SetBinContent(ibin, 1.);

           Double_t err = -999.;
           for(int i = 0; i < size; i++){
                   if(i==5 || i==7) continue;
                   // get "envelope"
                   Double_t temp_err =  fabs(hMCSys.at(i)->GetBinContent(ibin) - ratio->GetBinContent(ibin));
		   cout << i << " th sys, bin: " << ibin << " err: " << temp_err << " sys bin: " << hMCSys.at(i)->GetBinContent(ibin) << endl;
                   if( temp_err > err){
                          err = temp_err;
                   }
           }
           sysRatio->SetBinError(ibin, sqrt(pow(err,2) + pow(sysRatio->GetBinError(ibin),2)));
        }

        for(int i = 0; i < size; i++){
           delete hMCSysM.at(i);
           delete hMCSys.at(i);

        }

        hMCSys.clear();
        hMCSys.shrink_to_fit();

        hMCSysM.clear();
        hMCSysM.shrink_to_fit();
}

void drawRatio(TString outpdf, TUnfoldDensity* unfold_pt, TUnfoldDensity* unfold_mass, TFile *filein, TString channel){  

        gROOT->SetBatch();

        setTDRStyle();
        writeExtraText = true;       // if extra text
        extraText  = "Preliminary";

	// pt distribution for 5 mass bins in one histogram
	TH1* hunfolded_pt = unfold_pt->GetOutput("hunfolded_pt",0,0,"*[UO]",kFALSE);
        TH1 *histMCTruth_pt=unfold_pt->GetBias("histMCTruth_pt",0,0,"*[UO]",kFALSE);
	TH1* ratio=(TH1*)hunfolded_pt->Clone("ratio");
	initHistError(ratio, 1.); // TODO need to 
        ratio->Divide(histMCTruth_pt);

	// get systematic for ratio plot
        TH1* sysErrRatio = (TH1*)hunfolded_pt->Clone("hsysErrRatio");
	sysErrRatio->Reset();
	getRatioSys(unfold_pt, "AlphaS", (int)2, sysErrRatio);	
	getRatioSys(unfold_pt, "Scale", (int)9, sysErrRatio);	
	getRatioSys(unfold_pt, "unfoldsys", (int)1, sysErrRatio);	
	getRatioSys(unfold_pt, "PU", (int)2, sysErrRatio);	
	getRatioSys(unfold_pt, "trgSF", (int)2, sysErrRatio);	
	getRatioSys(unfold_pt, "recoSF", (int)2, sysErrRatio);	
	getRatioSys(unfold_pt, "IdSF", (int)2, sysErrRatio);	
	getRatioSys(unfold_pt, "IsoSF", (int)2, sysErrRatio);	

        // stat error for measurement
        TH1* statErrRatio = (TH1*)hunfolded_pt->Clone("hstatErrRatio");
        statErrRatio->Divide(histMCTruth_pt);
        initHistContent(statErrRatio, 1.);

        TH1* sysErrRatioMC = (TH1*)hunfolded_pt->Clone("hsysErrRatio");
	sysErrRatioMC->Reset();
        sysMCErrRatio(unfold_pt, sysErrRatioMC, "Pt", "AlphaS", (int)2, channel);
        sysMCErrRatio(unfold_pt, sysErrRatioMC, "Pt", "Scale", (int)9, channel);
        sysMCErrRatio(unfold_pt, sysErrRatioMC, "Pt", "unfoldsys", (int)1, channel);


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
	getAveragesSysPt(ScaleErr, "Scale",   (int)9, unfold_pt);

	vector<Double_t> UnfoldErr;
	getAveragesSysPt(UnfoldErr, "unfoldsys",   (int)1, unfold_pt);

	vector<Double_t> AlphaSErr;
	getAveragesSysPt(AlphaSErr, "AlphaS", (int)2, unfold_pt);

	// get average mass
	//
	getAveragesMass(meanmass_data, meanmasserr_data, hunfolded_mass);
	getAveragesMass(meanmass_mc, meanmasserr_mc, histMCTruth_mass);

	vector<Double_t> ScaleMassErr;
	getAveragesSysMass(ScaleMassErr, "Scale",   (int)9, unfold_mass);

	vector<Double_t> UnfoldMassErr;
	getAveragesSysMass(UnfoldMassErr, "unfoldsys",   (int)1, unfold_mass);

	vector<Double_t> AlphaSMassErr;
	getAveragesSysMass(AlphaSMassErr, "AlphaS", (int)2, unfold_mass);


  	TCanvas* c1=new TCanvas("c1", "c1", 50, 50, 800, 800);
  	c1->cd();
  	gStyle->SetOptStat(0);

  	TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
        setPadMargins(pad1);
        pad1->SetBottomMargin(0.0);
  	pad1->SetTicks(1);
  	pad1->SetLogy();
  	pad1->Draw();
  	pad1->cd();

  	hunfolded_pt->SetTitle("");
  	hunfolded_pt->Draw("p9histe");
  	hunfolded_pt->SetMarkerStyle(20);
  	hunfolded_pt->SetMarkerSize(1.2);
  	hunfolded_pt->SetLineColor(kBlack);
        setYaxisHist(hunfolded_pt);
  	hunfolded_pt->GetXaxis()->SetLabelSize(0.);
  	hunfolded_pt->GetYaxis()->SetTitle("Events/ bin");
  	//histMCTruth_pt->Draw("histsames");
  	histMCTruth_pt->Draw("p9sames");
  	histMCTruth_pt->SetMarkerStyle(20);
  	histMCTruth_pt->SetMarkerColor(2);
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
	//	grid_.DrawLine(binEdge, 100, binEdge, histMCTruth_pt->GetBinContent(ii+1) );
	}

	int totoalNbin = ratio->GetXaxis()->GetNbins();

	// draw vertical lines to divide each mass retion
        TLine grid_vertical;
        grid_vertical.SetLineColor(kBlack);
        grid_vertical.SetLineStyle(kSolid);
        for( int ii=0; ii<5; ii++ )
        {
                grid_vertical.DrawLine(totoalNbin/5 * (ii+1) + 0.5, hunfolded_pt->GetMinimum(), totoalNbin/5 * (ii+1) + 0.5, hunfolded_pt->GetMaximum() );
        }


  	TLegend* leg_ = new TLegend(0.7, 0.60, 0.9, 0.88,"","brNDC");
  	leg_->SetTextSize(0.04);
  	//leg_->SetFillStyle(1);
  	leg_->SetBorderSize(0);
  	leg_->AddEntry(hunfolded_pt, "Unfolded data", "p");
  	leg_->AddEntry(histMCTruth_pt, "Truth", "l");
  	leg_->Draw();

  	c1->cd();

  	TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.4);
        setPadMargins(pad2);
        pad2->SetTopMargin(0.0);
        pad2->SetBottomMargin(0.3);
  	pad2->SetTicks(1);
  	pad2->SetGridy(1);
  	pad2->Draw();
  	pad2->cd();

  	ratio->Draw("pe");
  	ratio->SetMarkerStyle(20);
  	ratio->SetMarkerSize(1.0);
  	ratio->SetLineColor(kRed);
  	ratio->SetMarkerColor(kRed);
  	ratio->SetMinimum(0.6);
  	ratio->SetMaximum(1.4);
  	ratio->SetTitle("");
        setYaxisHist(ratio);
        setXaxisHist(ratio);
  	ratio->GetYaxis()->SetNdivisions(515);
        ratio->GetYaxis()->SetTitle("Unfolded/MC");
        ratio->GetXaxis()->SetTitle("p_{T} bin");
        ratio->GetYaxis()->CenterTitle();

	// draw systematic error
  	sysErrRatio->Draw("E2same");
  	sysErrRatio->SetMarkerSize(0);
  	sysErrRatio->SetFillStyle(3004);
  	sysErrRatio->SetFillColor(kBlack);
  	//sysErrRatio->SetFillColorAlpha(kBlack,0.5);

	statErrRatio->Draw("pesame");
	statErrRatio->SetMarkerSize(0);
  	statErrRatio->SetLineColor(kBlack);

	sysErrRatioMC->Draw("E2same");
	sysErrRatioMC->SetMarkerSize(0);
  	//sysErrRatioMC->SetFillStyle(3004);
  	sysErrRatioMC->SetFillColor(kRed);
  	sysErrRatioMC->SetFillColorAlpha(kRed,0.3);

        TLine *l_;
        l_ = new TLine(ratio->GetXaxis()->GetXmin(),1,ratio->GetXaxis()->GetXmax(),1);
        l_->Draw("same");
        l_->SetLineStyle(1);
        l_->SetLineColor(kBlack);

        for( int ii=0; ii<5; ii++ )
        {
                grid_vertical.DrawLine(totoalNbin/5 * (ii+1) + 0.5, ratio->GetMinimum(), totoalNbin/5 * (ii+1) + 0.5, ratio->GetMaximum() );
        }


        c1->cd();
/*
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
        	//TGraphErrors *sysData = new TGraphErrors(1, &meanmass_data[imass], &meanpt_data[imass], &ScaleMassErr[imass], &ScaleErr[imass]);
        	TGraphErrors *sysData = new TGraphErrors(1, &meanmass_data[imass], &meanpt_data[imass], &UnfoldMassErr[imass], &UnfoldErr[imass]);
        	//TGraphErrors *sysData = new TGraphErrors(1, &meanmass_data[imass], &meanpt_data[imass], &AlphaSMassErr[imass], &AlphaSErr[imass]);
        	TGraphErrors *grMC = new TGraphErrors(1, &meanmass_mc[imass], &meanpt_mc[imass], &meanmasserr_mc[imass], &meanpterr_mc[imass]);

        	Double_t xmin = meanmass_data[imass] > meanmass_mc[imass] ? meanmass_mc[imass] * 0.99: meanmass_data[imass] * 0.99,
        	         xmax = meanmass_data[imass] < meanmass_mc[imass] ? meanmass_mc[imass] * 1.01: meanmass_data[imass] * 1.01,
        	         ymin = meanpt_data[imass] > meanpt_mc[imass] ? meanpt_mc[imass] * 0.95: meanpt_data[imass] * 0.95,
        	         ymax = meanpt_data[imass] < meanpt_mc[imass] ? meanpt_mc[imass] * 1.05: meanpt_data[imass] * 1.05;
        	setTGraphAxis(grUnfolded, meanmass_data[imass], xmin, xmax, ymin, ymin + 3., grMC, sysData);



	}
*/
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

void drawCombinedISR(TString outpdf, TUnfoldDensity* unfold_pt2016, TUnfoldDensity* unfold_mass2016, TUnfoldDensity* unfold_pt2017, TUnfoldDensity* unfold_mass2017){

        gROOT->SetBatch();

        setTDRStyle();
        writeExtraText = true;       // if extra text
        extraText  = "Preliminary";

        // mass distribution
        TH1* hunfolded_mass2016 = unfold_mass2016->GetOutput("hunfolded_mass2016",0,0,"*[UO]",kTRUE);
        TH1 *histMCTruth_mass2016=unfold_mass2016->GetBias("histMCTruth_mass2016",0,0,"*[UO]",kTRUE);

        vector<Double_t> meanmass_data2016, meanmasserr_data2016;
        vector<Double_t> meanmass_mc2016, meanmasserr_mc2016;

        vector<Double_t> meanpt_data2016, meanpterr_data2016;
        vector<Double_t> meanpt_mc2016, meanpterr_mc2016;

        // get average pt for each mass bin
        getAveragesPt(meanpt_data2016, meanpterr_data2016, unfold_pt2016, true);
        getAveragesPt(meanpt_mc2016, meanpterr_mc2016, unfold_pt2016, false);

        vector<Double_t> ScaleErr2016;
        getAveragesSysPt(ScaleErr2016, "Scale", (int)9, unfold_pt2016);

        vector<Double_t> AlphaSErr2016;
        getAveragesSysPt(AlphaSErr2016, "AlphaS", (int)2, unfold_pt2016);

        vector<Double_t> totalPtSys2016;
        vector<Double_t> totalPtErr2016;
        for(unsigned int i=0;i<ScaleErr2016.size();i++){
                totalPtErr2016.push_back(sqrt(pow(ScaleErr2016.at(i),2)+pow(AlphaSErr2016.at(i),2)+pow(meanpterr_data2016.at(i),2)));
                totalPtSys2016.push_back(sqrt(pow(ScaleErr2016.at(i),2)+pow(AlphaSErr2016.at(i),2)));
        }

        // get average mass
        //
        getAveragesMass(meanmass_data2016, meanmasserr_data2016, hunfolded_mass2016);
        getAveragesMass(meanmass_mc2016, meanmasserr_mc2016, histMCTruth_mass2016);

        vector<Double_t> ScaleMassErr2016;
        getAveragesSysMass(ScaleMassErr2016, "Scale", (int)9, unfold_mass2016);

        vector<Double_t> AlphaSMassErr2016;
        getAveragesSysMass(AlphaSMassErr2016, "AlphaS", (int)2, unfold_mass2016);

        vector<Double_t> totalMassSys2016;
        vector<Double_t> totalMassErr2016;
        for(unsigned int i=0;i<ScaleMassErr2016.size();i++){
                totalMassErr2016.push_back(sqrt(pow(ScaleMassErr2016.at(i),2)+pow(AlphaSMassErr2016.at(i),2)+pow(meanmasserr_data2016.at(i),2)));
                totalMassSys2016.push_back(sqrt(pow(ScaleMassErr2016.at(i),2)+pow(AlphaSMassErr2016.at(i),2)));
        }


        TH1* hunfolded_mass2017 = unfold_mass2017->GetOutput("hunfolded_mass2017",0,0,"*[UO]",kTRUE);
        TH1 *histMCTruth_mass2017=unfold_mass2017->GetBias("histMCTruth_mass2017",0,0,"*[UO]",kTRUE);

        vector<Double_t> meanmass_data2017, meanmasserr_data2017;
        vector<Double_t> meanmass_mc2017, meanmasserr_mc2017;

        vector<Double_t> meanpt_data2017, meanpterr_data2017;
        vector<Double_t> meanpt_mc2017, meanpterr_mc2017;

        // get average pt for each mass bin
        getAveragesPt(meanpt_data2017, meanpterr_data2017, unfold_pt2017, true);
        getAveragesPt(meanpt_mc2017, meanpterr_mc2017, unfold_pt2017, false);

        vector<Double_t> ScaleErr2017;
        getAveragesSysPt(ScaleErr2017, "Scale", (int)9, unfold_pt2017);

        vector<Double_t> AlphaSErr2017;
        getAveragesSysPt(AlphaSErr2017, "AlphaS", (int)2, unfold_pt2017);

        vector<Double_t> totalPtSys2017;
        vector<Double_t> totalPtErr2017;
        for(unsigned int i=0;i<ScaleErr2017.size();i++){
                totalPtErr2017.push_back(sqrt(pow(ScaleErr2017.at(i),2)+pow(AlphaSErr2017.at(i),2)+pow(meanpterr_data2017.at(i),2)));
                totalPtSys2017.push_back(sqrt(pow(ScaleErr2017.at(i),2)+pow(AlphaSErr2017.at(i),2)));
        }

        // get average mass
        //
        getAveragesMass(meanmass_data2017, meanmasserr_data2017, hunfolded_mass2017);
        getAveragesMass(meanmass_mc2017, meanmasserr_mc2017, histMCTruth_mass2017);

        vector<Double_t> ScaleMassErr2017;
        getAveragesSysMass(ScaleMassErr2017, "Scale", (int)9, unfold_mass2017);

        vector<Double_t> AlphaSMassErr2017;
        getAveragesSysMass(AlphaSMassErr2017, "AlphaS", (int)2, unfold_mass2017);

        vector<Double_t> totalMassSys2017;
        vector<Double_t> totalMassErr2017;
        for(unsigned int i=0;i<ScaleMassErr2017.size();i++){
                totalMassErr2017.push_back(sqrt(pow(ScaleMassErr2017.at(i),2)+pow(AlphaSMassErr2017.at(i),2)+pow(meanmasserr_data2017.at(i),2)));
                totalMassSys2017.push_back(sqrt(pow(ScaleMassErr2017.at(i),2)+pow(AlphaSMassErr2017.at(i),2)));
        }

        TCanvas * c1 = new TCanvas("","", 900*1.5, 1000*1.5);
        c1->cd();

        TPad *pad1 = new TPad("pad1","pad1",0,0.6,1,1);
        setPadMargins(pad1);
  	gStyle->SetOptStat(0);
        pad1->SetBottomMargin(0.01);
        pad1->SetTopMargin(0.1);
        pad1->SetTicks(1);
        pad1->SetLogx();
        pad1->SetGridy();
        pad1->Draw();
        pad1->cd();

        TGraphErrors *grUnfolded = new TGraphErrors(5, &meanmass_data2016[0], &meanpt_data2016[0], &totalMassErr2016[0], &totalPtErr2016[0]);
        grUnfolded->SetLineColor(kBlack);
        grUnfolded->SetMarkerColor(kBlack);
        grUnfolded->SetMarkerStyle(20);
        grUnfolded->SetMarkerSize(1.);
        grUnfolded->SetLineStyle(1);
        grUnfolded->SetLineWidth(1);
        grUnfolded->Draw("ape");
        setYaxisGraph(grUnfolded);

        TGraphErrors *grUnfolded2017 = new TGraphErrors(5, &meanmass_data2017[0], &meanpt_data2017[0], &totalMassErr2017[0], &totalPtErr2017[0]);
        grUnfolded2017->SetLineColor(kBlack);
        grUnfolded2017->SetMarkerColor(kBlack);
        grUnfolded2017->SetMarkerStyle(21);
        grUnfolded2017->SetMarkerSize(1.);
        grUnfolded2017->SetLineStyle(1);
        grUnfolded2017->SetLineWidth(1);
        grUnfolded2017->Draw("samepe");
	
        vector<Double_t> meanmass_datacombined, meanmasserr_datacombined;
        vector<Double_t> meanpt_datacombined, meanpterr_datacombined;
        vector<Double_t> meanptsyserr_datacombined, meanmasssyserr_datacombined;

	for(int i = 0; i < 5; i++){
		meanpterr_datacombined.push_back( 1./sqrt( 1./pow(meanpterr_data2016.at(i),2) + 1./pow(meanpterr_data2017.at(i),2) ) );
		meanptsyserr_datacombined.push_back( (totalPtSys2016.at(i) + totalPtSys2017.at(i))/2. );
		meanpt_datacombined.push_back( (meanpt_data2016.at(i)/pow(totalPtErr2016.at(i),2) + meanpt_data2017.at(i)/pow(totalPtErr2017.at(i),2)) / ( 1./pow(totalPtErr2016.at(i),2) + 1./pow(totalPtErr2017.at(i),2) )  );

		meanmasserr_datacombined.push_back( 1./sqrt( 1./pow(meanmasserr_data2016.at(i),2) + 1./pow(meanmasserr_data2017.at(i),2) ) );
		meanmasssyserr_datacombined.push_back( (totalMassSys2016.at(i) + totalMassSys2017.at(i))/2. );
		meanmass_datacombined.push_back( (meanmass_data2016.at(i)/pow(totalPtErr2016.at(i),2) + meanmass_data2017.at(i)/pow(totalPtErr2017.at(i),2)) / ( 1./pow(totalPtErr2016.at(i),2) + 1./pow(totalPtErr2017.at(i),2) )  );
	}

        TGraphErrors *grUnfoldedCombined = new TGraphErrors(5, &meanmass_datacombined[0], &meanpt_datacombined[0], &meanmasserr_datacombined[0], &meanpterr_datacombined[0]);
        grUnfoldedCombined->SetLineColor(kRed);
        grUnfoldedCombined->SetMarkerColor(kRed);
        grUnfoldedCombined->SetMarkerStyle(31);
        grUnfoldedCombined->SetMarkerSize(1.);
        grUnfoldedCombined->SetLineStyle(1);
        grUnfoldedCombined->SetLineWidth(1);
        grUnfoldedCombined->Draw("samepe");
	

        TF1 *f1 = new TF1("f1", "[0]+[1]*log(x)", 40., 350.);
        f1->GetXaxis()->SetRangeUser(40., 350.);
        f1->SetLineColor(kBlack);
        grUnfolded->Fit(f1, "R"); // R: fitting sub range
        f1->Draw("same");

        TF1 *f2 = new TF1("f2", "[0]+[1]*log(x)", 40., 350.);
        f2->GetXaxis()->SetRangeUser(40., 350.);
        f2->SetLineColor(kBlack);
        grUnfolded2017->Fit(f2, "R"); // R: fitting sub range
        f2->Draw("same");

        c1->cd();

        Double_t initial_padx = 0.015;
        TLatex ptstring;
        for(int imass = 0; imass < 5; imass++){

		// set dummy mass mean position
		meanmass_data2017[imass] = meanmass_data2016[imass] + 0.5;
		meanmass_datacombined[imass] = meanmass_data2016[imass] + 0.25;
		meanmasserr_data2016[imass] = 0.1;
		totalMassSys2016[imass] = 0.1;
		meanmasserr_data2017[imass] = 0.1;
		totalMassSys2017[imass] = 0.1;
		meanmasserr_datacombined[imass] = 0.1;
		meanmasssyserr_datacombined[imass] = 0.1;

                TString imass_;
                imass_.Form("%d", imass);

                c1->cd();

                TPad *pad = new TPad("mass_" + imass_,"mass_" + imass_, initial_padx + (0.19 * imass) + (0.005 * imass), 0., initial_padx + (0.19 * (imass + 1)) + (0.005 * imass), 0.6);
                setPadMargins(pad);
                pad->SetLeftMargin(0.01);
                pad->SetRightMargin(0.01);
                pad->SetBottomMargin(0.1);
                pad->SetTicks(1);
                pad->SetGrid(0,1);
                //pad->SetLogx();
                pad->Draw();
                pad->cd();

                TLegend* leg_ = new TLegend(0.03, 0.6, 0.7, 0.92,"","brNDC");
                leg_->SetTextSize(0.08);
                //leg_->SetFillStyle(0);
                leg_->SetBorderSize(0);

                // FIXME how to handle below memory allocation
                TGraphErrors *grUnfolded = new TGraphErrors(1, &meanmass_data2016[imass], &meanpt_data2016[imass], &meanmasserr_data2016[imass], &meanpterr_data2016[imass]);
                TGraphErrors *sysData = new TGraphErrors(1, &meanmass_data2016[imass], &meanpt_data2016[imass], &totalMassSys2016[imass], &totalPtSys2016[imass]);

                Double_t xmin = meanmass_data2016[imass] > meanmass_mc2016[imass] ? meanmass_mc2016[imass] - 1.: meanmass_data2016[imass] - 1.,
                         xmax = meanmass_data2016[imass] < meanmass_mc2016[imass] ? meanmass_mc2016[imass] + 1.: meanmass_data2016[imass] + 1.,
                         ymin = meanpt_data2016[imass] > meanpt_mc2016[imass] ? meanpt_mc2016[imass] - 0.7: meanpt_data2016[imass] - 0.7,
                         ymax = meanpt_data2016[imass] < meanpt_mc2016[imass] ? meanpt_mc2016[imass] + 0.7: meanpt_data2016[imass] + 0.7;
                setTGraphAxis(grUnfolded, meanmass_data2016[imass], xmin, xmax, ymin, ymin + 3., nullptr, sysData, true);

                TString pt;
                TString stat;
                TString sys;
                TString total;
                pt.Form("%.2f", meanpt_data2016[imass]);
                stat.Form("%.2f", meanpterr_data2016[imass]);
                sys.Form("%.2f", totalPtSys2016[imass]);
                total.Form("%.2f", sqrt(pow(totalPtSys2016[imass],2)+pow( meanpterr_data2016[imass],2) ));
                ptstring.SetTextSize(0.12);
		leg_->AddEntry(grUnfolded, "#splitline{2016: "+pt+"#pm"+total+"}{        (#pm"+stat+"#pm"+sys+")}");
                //ptstring.DrawLatex(meanmass_data2016[imass]*0.99, meanpt_data2016[imass]*1.02, pt+"#pm"+stat+"#pm"+sys);


                TGraphErrors *grUnfolded2017 = new TGraphErrors(1, &meanmass_data2017[imass], &meanpt_data2017[imass], &meanmasserr_data2017[imass], &meanpterr_data2017[imass]);
                TGraphErrors *sysData2017 = new TGraphErrors(1, &meanmass_data2017[imass], &meanpt_data2017[imass], &totalMassSys2017[imass], &totalPtSys2017[imass]);

                setTGraphAxis(grUnfolded2017, meanmass_data2017[imass], xmin, xmax, ymin, ymin + 3., nullptr, sysData2017, false); 

                pt.Form("%.2f", meanpt_data2017[imass]);
                stat.Form("%.2f", meanpterr_data2017[imass]);
                sys.Form("%.2f", totalPtSys2017[imass]);
                total.Form("%.2f", sqrt(pow(totalPtSys2017[imass],2)+pow( meanpterr_data2017[imass],2) ));
                ptstring.SetTextSize(0.12);
                //ptstring.DrawLatex(meanmass_data2017[imass]*0.99, meanpt_data2017[imass]*1.02, pt+"#pm"+stat+"#pm"+sys);

		//leg_->AddEntry(grUnfolded2017, "#splitline{2017: }{"+pt+"#pm"+stat+"#pm"+sys+"}");
		leg_->AddEntry(grUnfolded2017, "#splitline{2017: "+pt+"#pm"+total+"}{        (#pm"+stat+"#pm"+sys+")}");

                TGraphErrors *grUnfoldedcombined = new TGraphErrors(1, &meanmass_datacombined[imass], &meanpt_datacombined[imass], &meanmasserr_datacombined[imass], &meanpterr_datacombined[imass]);
                TGraphErrors *sysDataCombined = new TGraphErrors(1, &meanmass_datacombined[imass], &meanpt_datacombined[imass], &meanmasssyserr_datacombined[imass], &meanptsyserr_datacombined[imass]);

                setTGraphAxis(grUnfoldedcombined, meanmass_datacombined[imass], xmin, xmax, ymin, ymin + 3., nullptr, sysDataCombined, false);
		grUnfoldedcombined->SetMarkerStyle(31);
		grUnfoldedcombined->SetMarkerColor(kRed);
		grUnfoldedcombined->SetLineColor(kRed);

                pt.Form("%.2f", meanpt_datacombined[imass]);
                stat.Form("%.2f", meanpterr_datacombined[imass]);
                sys.Form("%.2f", meanptsyserr_datacombined[imass]);
                total.Form("%.2f", sqrt(pow(meanptsyserr_datacombined[imass],2)+pow( meanpterr_datacombined[imass],2) ));
                ptstring.SetTextSize(0.12);
                //ptstring.DrawLatex(meanmass_datacombined[imass]*0.99, meanpt_datacombined[imass]*1.02, pt+"#pm"+stat+"#pm"+sys);
		leg_->AddEntry(grUnfoldedcombined, "#splitline{Combined: "+pt+"#pm"+total+"}{        (#pm"+stat+"#pm"+sys+")}");

		leg_->Draw();


        }


        CMS_lumi( c1, 4, 0 );
        c1->cd();
        c1->SaveAs(outpdf);



}

void drawISRfit(TString outpdf, TUnfoldDensity* unfold_pt, TUnfoldDensity* unfold_mass, TFile *filein){

        gROOT->SetBatch();

        setTDRStyle();
        writeExtraText = true;       // if extra text
        extraText  = "Preliminary";

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

        vector<Double_t> AlphaSErr;
        getAveragesSysPt(AlphaSErr, "AlphaS", (int)2, unfold_pt);

        vector<Double_t> UnfoldErr;
        getAveragesSysPt(UnfoldErr, "unfoldsys", (int)1, unfold_pt);

 	vector<Double_t> totalPtSys;
 	vector<Double_t> totalPtErr;
	for(unsigned int i=0;i<ScaleErr.size();i++){
		totalPtErr.push_back(sqrt(pow(ScaleErr.at(i),2)+pow(AlphaSErr.at(i),2)+pow(UnfoldErr.at(i),2)+pow(meanpterr_data.at(i),2)));
		totalPtSys.push_back(sqrt(pow(ScaleErr.at(i),2)+pow(AlphaSErr.at(i),2)+pow(UnfoldErr.at(i),2)));
	}

        // get average mass
        //
        getAveragesMass(meanmass_data, meanmasserr_data, hunfolded_mass);
        getAveragesMass(meanmass_mc, meanmasserr_mc, histMCTruth_mass);

        vector<Double_t> ScaleMassErr;
        getAveragesSysMass(ScaleMassErr, "Scale", (int)9, unfold_mass);

        vector<Double_t> AlphaSMassErr;
        getAveragesSysMass(AlphaSMassErr, "AlphaS", (int)2, unfold_mass);

        vector<Double_t> UnfoldMassErr;
        getAveragesSysMass(UnfoldMassErr, "unfoldsys", (int)1, unfold_mass);

        vector<Double_t> totalMassSys;
        vector<Double_t> totalMassErr;
        for(unsigned int i=0;i<ScaleMassErr.size();i++){
                totalMassErr.push_back(sqrt(pow(ScaleMassErr.at(i),2)+pow(AlphaSMassErr.at(i),2)+pow(UnfoldMassErr.at(i),2)+pow(meanmasserr_data.at(i),2)));
                totalMassSys.push_back(sqrt(pow(ScaleMassErr.at(i),2)+pow(AlphaSMassErr.at(i),2)+pow(UnfoldMassErr.at(i),2)));
        }

 	TCanvas * c1 = new TCanvas("c1","c1", 50, 50, 850, 700);
	c1->cd();
	gStyle->SetOptFit(0);

        TPad *pad1 = new TPad("pad1","pad1",0,0.0,1,1);
        setPadMargins(pad1);
        pad1->SetBottomMargin(0.2);
        pad1->SetTopMargin(0.08);
        pad1->SetTicks(1);
        pad1->SetLogx();
        pad1->SetGridy();
        pad1->Draw();
        pad1->cd();

 	TGraphErrors *grUnfolded = new TGraphErrors(5, &meanmass_data[0], &meanpt_data[0], &totalMassErr[0], &totalPtErr[0]);
 	grUnfolded->SetLineColor(kBlack);
 	grUnfolded->SetMarkerColor(kBlack);
 	grUnfolded->SetMarkerStyle(20);
 	grUnfolded->SetMarkerSize(1.);
 	grUnfolded->SetLineStyle(1);
 	grUnfolded->Draw("ape");
 	grUnfolded->GetYaxis()->SetRangeUser(10.,30.);
 	grUnfolded->GetXaxis()->SetLimits(30.,800.);
	grUnfolded->GetYaxis()->SetTitle("Average p_{T} (GeV)");
	grUnfolded->GetXaxis()->SetTitle("Average Mass (GeV)");
	setYaxisGraph(grUnfolded);
	setXaxisGraph(grUnfolded);

 	TF1 *f1 = new TF1("f1", "[0]+[1]*log(x)", 40., 350.);
 	f1->GetXaxis()->SetRangeUser(40., 350.);
 	f1->SetLineColor(kBlack);
 	grUnfolded->Fit(f1, "R0"); // R: fitting sub range
 	f1->Draw("same");

        c1->cd();
/*
        Double_t initial_padx = 0.01;
 	TLatex ptstring;
        for(int imass = 0; imass < 5; imass++){

                TString imass_;
                imass_.Form("%d", imass);

                c1->cd();

                TPad *pad = new TPad("mass_" + imass_,"mass_" + imass_, initial_padx + (0.18 * imass) + (0.01 * imass), 0., initial_padx + (0.18 * (imass + 1)) + (0.01 * imass), 0.55);
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
                TGraphErrors *sysData = new TGraphErrors(1, &meanmass_data[imass], &meanpt_data[imass], &totalMassSys[imass], &totalPtSys[imass]);
                TGraphErrors *grMC = new TGraphErrors(1, &meanmass_mc[imass], &meanpt_mc[imass], &meanmasserr_mc[imass], &meanpterr_mc[imass]);

                Double_t xmin = meanmass_data[imass] > meanmass_mc[imass] ? meanmass_mc[imass] - 1.: meanmass_data[imass] - 1.,
                         xmax = meanmass_data[imass] < meanmass_mc[imass] ? meanmass_mc[imass] + 1.: meanmass_data[imass] + 1.,
                         ymin = meanpt_data[imass] > meanpt_mc[imass] ? meanpt_mc[imass] - 0.7: meanpt_data[imass] - 0.7,
                         ymax = meanpt_data[imass] < meanpt_mc[imass] ? meanpt_mc[imass] + 0.7: meanpt_data[imass] + 0.7;
                setTGraphAxis(grUnfolded, meanmass_data[imass], xmin, xmax, ymin, ymin + 3., grMC, sysData);

		TString pt;
		TString stat;
		TString sys;
 		pt.Form("%.2f", meanpt_data[imass]);
 		stat.Form("%.2f", meanpterr_data[imass]);
 		sys.Form("%.2f", totalPtSys[imass]);
 		ptstring.SetTextSize(0.12);
 		ptstring.DrawLatex(meanmass_data[imass]*0.99, meanpt_data[imass]*1.02, pt+"#pm"+stat+"#pm"+sys);

        }
*/

        CMS_lumi( pad1, 4, 11 );
        c1->cd();
        c1->SaveAs(outpdf);

}

void drawMassRatio(TString outpdf, TUnfoldDensity* unfold, TFile *filein){

        gROOT->SetBatch();

        setTDRStyle();
        writeExtraText = true;       // if extra text
        extraText  = "Preliminary";

	TH1* hunfolded = unfold->GetOutput("hunfolded",0,0,"*[UO]",kTRUE);
        TH1 *histMCTruth=unfold->GetBias("histMCTruth",0,0,"*[UO]",kTRUE);
	TH1* ratio=(TH1*)hunfolded->Clone("ratio");

	// systematic error for measurement
        TH1* sysErrRatio = (TH1*)hunfolded->Clone("hsysErrRatio");
        sysErrRatio->Reset();

	// stat error for measurement
        TH1* statErrRatio = (TH1*)hunfolded->Clone("hstatErrRatio");
        statErrRatio->Divide(histMCTruth);
	initHistContent(statErrRatio, 1.);
        //statErrRatio->Reset();
	
        getRatioSysMass(unfold, "AlphaS", (int)2, sysErrRatio);
        getRatioSysMass(unfold, "Scale", (int)9, sysErrRatio);
	getRatioSysMass(unfold, "unfoldsys", (int)1, sysErrRatio);

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
  	hunfolded->GetYaxis()->SetTitle("Events/ bin");
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
  	ratio->SetMinimum(0.0);
  	ratio->SetMaximum(2.0);
  	ratio->SetTitle("");
        setYaxisHist(ratio);
        setXaxisHist(ratio);
  	ratio->GetYaxis()->SetNdivisions(505);

        // draw systematic error
        sysErrRatio->Draw("E2same");
        sysErrRatio->SetMarkerSize(0);
        sysErrRatio->SetFillColorAlpha(kBlack,0.3);

	statErrRatio->Draw("E2same");
        sysErrRatio->SetMarkerSize(0);

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

void drawPtReco(TString outpdf, TString postfix, TFile *fdata, TFile *fDYsig, TFile *fDYbkg, TFile *fTTbar, TFile *fVV, TFile *fWjets, TFile *fqcd, TString channel){

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

	//TUnfoldBinning *ptbin = ptBinning_rec();
	TUnfoldBinning* ptbin=(TUnfoldBinning*)fdata->Get("Rec_Pt");

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

        TH1* hqcd;
	TH1* hqcdNoUO;
        if(fqcd != NULL){ 
		hqcd = (TH1*)fqcd->Get(hname);
		hqcdNoUO = ptbin->ExtractHistogram("hqcd", hqcd, 0, kFALSE, "*[UO]");
	}

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
        if(fqcd != NULL) hqcdNoUO->SetLineColor(linecolorEWK-2);
        if(fqcd != NULL) hqcdNoUO->SetFillColor(fillcolorEWK-2);

	if(fqcd != NULL) hsMCs->Add(hqcdNoUO);
	hsMCs->Add(hwjetsNoUO);
	hsMCs->Add(hvvNoUO);
	hsMCs->Add(httbarNoUO);
	hsMCs->Add(hdybkgNoUO);
	hsMCs->Add(hdysigNoUO);

        TCanvas* c1=new TCanvas("c1", "c1", 50, 50, 950*1.5, 1000*1.5);
        c1->cd();
        gStyle->SetOptStat(0);

        TPad *pad1 = new TPad("pad1","pad1",0,0.35,1,1);
        setPadMargins(pad1);
        pad1->SetTopMargin(0.1);
        pad1->SetBottomMargin(0.01);
        pad1->SetTicks(0,1);
        pad1->SetLogy();
        pad1->Draw();
        pad1->cd();
        
        hdataNoUO->SetTitle("");
        hdataNoUO->Draw("p9histe");
        hdataNoUO->SetMarkerStyle(20);
        hdataNoUO->SetMarkerSize(1.5);
        hdataNoUO->SetLineColor(kBlack);
        setYaxisHist(hdataNoUO);
	hsMCs->Draw("hist same");
        hdataNoUO->GetXaxis()->SetLabelSize(0.);
        hdataNoUO->GetYaxis()->SetTitle("Events/ Bin");
        hdataNoUO->GetYaxis()->SetRangeUser(5.,hdataNoUO->GetMaximum() * 1e2);
        hdataNoUO->Draw("p9histsamee");
        pad1->RedrawAxis();

        int totoalNbin = hdataNoUO->GetXaxis()->GetNbins();

        TLine grid_vertical;
        grid_vertical.SetLineColor(kBlack);
        grid_vertical.SetLineStyle(kSolid);
        for( int ii=0; ii<5; ii++ )
        {
                grid_vertical.DrawLine(totoalNbin/5 * (ii+1) + 0.5, hdataNoUO->GetMinimum(), totoalNbin/5 * (ii+1) + 0.5, hdataNoUO->GetMaximum() );
        }

	TLegend *fLeg = new TLegend(0.6,0.6,0.95,0.9);
        fLeg->SetTextSize(50);
        fLeg->SetTextFont(63);
        fLeg->SetFillStyle(0);
        fLeg->SetBorderSize(0);
 	fLeg->AddEntry(hdataNoUO, "Data", "lep");
 	if( channel.CompareTo("electron") == 0 ) fLeg->AddEntry(hdysigNoUO, "Z/#gamma* #rightarrow ee", "F");
 	if( channel.CompareTo("muon") == 0 )     fLeg->AddEntry(hdysigNoUO, "Z/#gamma* #rightarrow #mu#mu", "F");
 	fLeg->AddEntry(hdybkgNoUO, "Z/#gamma* #rightarrow #tau#tau", "F");
 	fLeg->AddEntry(httbarNoUO, "ttbar", "F");
 	fLeg->AddEntry(hvvNoUO, "VV", "F");
 	fLeg->AddEntry(hwjetsNoUO, "Wjets", "F");
 	if(fqcd != NULL) fLeg->AddEntry(hqcdNoUO, "QCD", "F");
	fLeg->Draw();

        c1->cd();

        TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.35);
        setPadMargins(pad2);
        pad2->SetBottomMargin(0.22);
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
	if(fqcd != NULL) hmcs->Add(hqcdNoUO);

        ratio->Divide(hmcs);
        ratio->Draw("pe");
        ratio->SetMarkerStyle(20);
        ratio->SetMarkerSize(1.5);
        ratio->SetLineColor(kBlack);
        ratio->SetMinimum(0.6);
        ratio->SetMaximum(1.4);
        ratio->GetYaxis()->SetTitle("Data/MC");
        ratio->GetYaxis()->CenterTitle();
        setYaxisHist(ratio);
        setXaxisHist(ratio);
        ratio->GetYaxis()->SetNdivisions(515);
        if( channel.CompareTo("electron") == 0 ) ratio->GetXaxis()->SetTitle("Mass, p_{T}(ee) bin");
        if( channel.CompareTo("muon") == 0 ) ratio->GetXaxis()->SetTitle("Mass, p_{T}(#mu#mu) bin");

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
/*
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
*/

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

void drawMassReco(TString outpdf, TString postfix, TFile *fdata, TFile *fDYsig, TFile *fDYbkg, TFile *fTTbar, TFile *fVV, TFile *fWjets, TFile *fqcd, TString channel){

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

	//TUnfoldBinning *massbin = massBinning_rec();
	TUnfoldBinning* massbin=(TUnfoldBinning*)fdata->Get("Rec_Mass");

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

        TH1* hqcd;
        TH1* hqcdNoUO;
        if(fqcd != NULL){
                hqcd = (TH1*)fqcd->Get(hname);
                hqcdNoUO = massbin->ExtractHistogram("hqcd", hqcd, 0, kTRUE, "*[UO]");
        }

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
        if(fqcd != NULL) hqcdNoUO->SetLineColor(linecolorEWK-2);
        if(fqcd != NULL) hqcdNoUO->SetFillColor(fillcolorEWK-2);

        if(fqcd != NULL) hsMCs->Add(hqcdNoUO);
	hsMCs->Add(hwjetsNoUO);
	hsMCs->Add(hvvNoUO);
	hsMCs->Add(httbarNoUO);
	hsMCs->Add(hdybkgNoUO);
	hsMCs->Add(hdysigNoUO);

        TCanvas* c1=new TCanvas("c1", "c1", 950*1.5, 1000*1.5);
        c1->cd();
        gStyle->SetOptStat(0);

        TPad *pad1 = new TPad("pad1","pad1",0,0.35,1,1);
        setPadMargins(pad1);
        pad1->SetBottomMargin(0.01);
        pad1->SetTopMargin(0.1);
        pad1->SetTicks(0,1);
        pad1->SetLogy();
        pad1->Draw();
        pad1->cd();
        
        hdataNoUO->SetTitle("");
        hdataNoUO->Draw("p9histe");
        hdataNoUO->SetMarkerStyle(20);
        hdataNoUO->SetMarkerSize(1.5);
        hdataNoUO->SetLineColor(kBlack);
        setYaxisHist(hdataNoUO);
	hsMCs->Draw("hist same");
        hdataNoUO->GetXaxis()->SetLabelSize(0.);
        hdataNoUO->GetYaxis()->SetTitle("Events/ bin");
        hdataNoUO->GetYaxis()->SetRangeUser(5.,hdataNoUO->GetMaximum() * 0.8e1);
        hdataNoUO->Draw("p9histsamee");
	pad1->RedrawAxis();
	
        TLegend *fLeg = new TLegend(0.6,0.5,0.95,0.9);
        fLeg->SetTextSize(50);
        fLeg->SetTextFont(63);
        fLeg->SetFillStyle(0);
        fLeg->SetBorderSize(0);
 	fLeg->AddEntry(hdataNoUO, "data", "lep");
 	if( channel.CompareTo("electron") == 0 ) fLeg->AddEntry(hdysigNoUO, "Z/#gamma* #rightarrow ee", "F");
 	if( channel.CompareTo("muon") == 0 )     fLeg->AddEntry(hdysigNoUO, "Z/#gamma* #rightarrow #mu#mu", "F");
 	fLeg->AddEntry(hdybkgNoUO, "Z/#gamma* #rightarrow #tau#tau", "F");
 	fLeg->AddEntry(httbarNoUO, "ttbar", "F");
 	fLeg->AddEntry(hvvNoUO, "VV", "F");
 	fLeg->AddEntry(hwjetsNoUO, "Wjets", "F");
        if(fqcd != NULL) fLeg->AddEntry(hqcdNoUO, "QCD", "F");
	fLeg->Draw();

        c1->cd();

        TPad *pad2 = new TPad("pad2","pad2",0,0.0,1,0.35);
        setPadMargins(pad2);
        pad2->SetBottomMargin(0.22);
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
        if(fqcd != NULL) hmcs->Add(hqcdNoUO);

        ratio->Divide(hmcs);
        ratio->Draw("pe");
        ratio->SetMarkerStyle(20);
        ratio->SetMarkerSize(1.5);
        ratio->SetLineColor(kBlack);
        ratio->SetMinimum(0.0);
        ratio->SetMaximum(2.0);
        ratio->GetYaxis()->SetTitle("Data/MC");
        ratio->GetYaxis()->CenterTitle();
        setYaxisHist(ratio);
        setXaxisHist(ratio);
        ratio->GetYaxis()->SetNdivisions(510);
        if( channel.CompareTo("electron") == 0 ) ratio->GetXaxis()->SetTitle("Mass(ee) (GeV)");
        if( channel.CompareTo("muon") == 0 ) ratio->GetXaxis()->SetTitle("Mass(#mu#mu) (GeV)");

        TLine *l_;
        l_ = new TLine(ratio->GetXaxis()->GetXmin(),1,ratio->GetXaxis()->GetXmax(),1);
        l_->Draw("same");
        l_->SetLineStyle(1);
        l_->SetLineColor(kRed);

        c1->cd();
/*
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
*/
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

void responseM(TString outpdf, TUnfoldDensity* unfold){
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

void efficiency(TString outpdf, TUnfoldDensity* unfold){

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

void drawEMuCombinedISR(TString outpdf, TUnfoldDensity* unfold_ptElectron, TUnfoldDensity* unfold_massElectron, TUnfoldDensity* unfold_ptMuon, TUnfoldDensity* unfold_massMuon){

        gROOT->SetBatch();

        setTDRStyle();
        writeExtraText = true;       // if extra text
        extraText  = "Preliminary";

        // mass distribution
        TH1* hunfolded_massElectron = unfold_massElectron->GetOutput("hunfolded_massElectron",0,0,"*[UO]",kTRUE);
        TH1 *histMCTruth_massElectron=unfold_massElectron->GetBias("histMCTruth_massElectron",0,0,"*[UO]",kTRUE);

        vector<Double_t> meanmass_dataElectron, meanmasserr_dataElectron;
        vector<Double_t> meanmass_mcElectron, meanmasserr_mcElectron;

        vector<Double_t> meanpt_dataElectron, meanpterr_dataElectron;
        vector<Double_t> meanpt_mcElectron, meanpterr_mcElectron;

        // get average pt for each mass bin
        getAveragesPt(meanpt_dataElectron, meanpterr_dataElectron, unfold_ptElectron, true);
        getAveragesPt(meanpt_mcElectron, meanpterr_mcElectron, unfold_ptElectron, false);

        vector<Double_t> ScaleErrElectron;
        getAveragesSysPt(ScaleErrElectron, "Scale", (int)9, unfold_ptElectron);

        vector<Double_t> AlphaSErrElectron;
        getAveragesSysPt(AlphaSErrElectron, "AlphaS", (int)2, unfold_ptElectron);

        vector<Double_t> unfoldsysErrElectron;
        getAveragesSysPt(unfoldsysErrElectron, "unfoldsys", (int)1, unfold_ptElectron);

        vector<Double_t> totalPtSysElectron;
        vector<Double_t> totalPtErrElectron;
        for(unsigned int i=0;i<ScaleErrElectron.size();i++){
                totalPtErrElectron.push_back(sqrt(pow(ScaleErrElectron.at(i),2)+pow(AlphaSErrElectron.at(i),2)+pow(meanpterr_dataElectron.at(i),2)+pow(unfoldsysErrElectron.at(i),2)));
                totalPtSysElectron.push_back(sqrt(pow(ScaleErrElectron.at(i),2)+pow(AlphaSErrElectron.at(i),2)));
        }

        // get average mass
        //
        getAveragesMass(meanmass_dataElectron, meanmasserr_dataElectron, hunfolded_massElectron);
        getAveragesMass(meanmass_mcElectron, meanmasserr_mcElectron, histMCTruth_massElectron);

        vector<Double_t> ScaleMassErrElectron;
        getAveragesSysMass(ScaleMassErrElectron, "Scale", (int)9, unfold_massElectron);

        vector<Double_t> AlphaSMassErrElectron;
        getAveragesSysMass(AlphaSMassErrElectron, "AlphaS", (int)2, unfold_massElectron);

        vector<Double_t> totalMassSysElectron;
        vector<Double_t> totalMassErrElectron;
        for(unsigned int i=0;i<ScaleMassErrElectron.size();i++){
                totalMassErrElectron.push_back(sqrt(pow(ScaleMassErrElectron.at(i),2)+pow(AlphaSMassErrElectron.at(i),2)+pow(meanmasserr_dataElectron.at(i),2)));
                totalMassSysElectron.push_back(sqrt(pow(ScaleMassErrElectron.at(i),2)+pow(AlphaSMassErrElectron.at(i),2)));
        }


        TH1* hunfolded_massMuon = unfold_massMuon->GetOutput("hunfolded_massMuon",0,0,"*[UO]",kTRUE);
        TH1 *histMCTruth_massMuon=unfold_massMuon->GetBias("histMCTruth_massMuon",0,0,"*[UO]",kTRUE);

        vector<Double_t> meanmass_dataMuon, meanmasserr_dataMuon;
        vector<Double_t> meanmass_mcMuon, meanmasserr_mcMuon;

        vector<Double_t> meanpt_dataMuon, meanpterr_dataMuon;
        vector<Double_t> meanpt_mcMuon, meanpterr_mcMuon;

        // get average pt for each mass bin
        getAveragesPt(meanpt_dataMuon, meanpterr_dataMuon, unfold_ptMuon, true);
        getAveragesPt(meanpt_mcMuon, meanpterr_mcMuon, unfold_ptMuon, false);

        vector<Double_t> ScaleErrMuon;
        getAveragesSysPt(ScaleErrMuon, "Scale", (int)9, unfold_ptMuon);

        vector<Double_t> AlphaSErrMuon;
        getAveragesSysPt(AlphaSErrMuon, "AlphaS", (int)2, unfold_ptMuon);

        vector<Double_t> totalPtSysMuon;
        vector<Double_t> totalPtErrMuon;
        for(unsigned int i=0;i<ScaleErrMuon.size();i++){
                totalPtErrMuon.push_back(sqrt(pow(ScaleErrMuon.at(i),2)+pow(AlphaSErrMuon.at(i),2)+pow(meanpterr_dataMuon.at(i),2)));
                totalPtSysMuon.push_back(sqrt(pow(ScaleErrMuon.at(i),2)+pow(AlphaSErrMuon.at(i),2)));
        }

        // get average mass
        //
        getAveragesMass(meanmass_dataMuon, meanmasserr_dataMuon, hunfolded_massMuon);
        getAveragesMass(meanmass_mcMuon, meanmasserr_mcMuon, histMCTruth_massMuon);

        vector<Double_t> ScaleMassErrMuon;
        getAveragesSysMass(ScaleMassErrMuon, "Scale", (int)9, unfold_massMuon);

        vector<Double_t> AlphaSMassErrMuon;
        getAveragesSysMass(AlphaSMassErrMuon, "AlphaS", (int)2, unfold_massMuon);

        vector<Double_t> totalMassSysMuon;
        vector<Double_t> totalMassErrMuon;
        for(unsigned int i=0;i<ScaleMassErrMuon.size();i++){
                totalMassErrMuon.push_back(sqrt(pow(ScaleMassErrMuon.at(i),2)+pow(AlphaSMassErrMuon.at(i),2)+pow(meanmasserr_dataMuon.at(i),2)));
                totalMassSysMuon.push_back(sqrt(pow(ScaleMassErrMuon.at(i),2)+pow(AlphaSMassErrMuon.at(i),2)));
        }

        TCanvas * c1 = new TCanvas("c1","c1", 50, 50, 850, 700);
        c1->cd();
	gStyle->SetOptFit(0);

        TPad *pad1 = new TPad("pad1","pad1",0,0.0,1,1);
        setPadMargins(pad1);
        pad1->SetBottomMargin(0.2);
        pad1->SetTopMargin(0.08);
        pad1->SetTicks(1);
        pad1->SetLogx();
        pad1->SetGridy();
        pad1->Draw();
        pad1->cd();

        TGraphErrors *grUnfolded = new TGraphErrors(5, &meanmass_dataElectron[0], &meanpt_dataElectron[0], &totalMassErrElectron[0], &totalPtErrElectron[0]);
        grUnfolded->SetLineColor(kBlack);
        grUnfolded->SetMarkerColor(kBlack);
        grUnfolded->SetMarkerStyle(24);
        grUnfolded->SetMarkerSize(1.5);
        grUnfolded->SetLineStyle(1);
        grUnfolded->SetLineWidth(1);
        grUnfolded->Draw("ape");
        grUnfolded->GetYaxis()->SetTitle("Average p_{T} (GeV)");
        grUnfolded->GetXaxis()->SetTitle("Average Mass (GeV)");
        setYaxisGraph(grUnfolded);
        setXaxisGraph(grUnfolded);
        //grUnfolded->GetYaxis()->SetRangeUser(0.,30.);
        //grUnfolded->GetXaxis()->SetLimits(30.,800.);

        TGraphErrors *grUnfoldedMuon = new TGraphErrors(5, &meanmass_dataMuon[0], &meanpt_dataMuon[0], &totalMassErrMuon[0], &totalPtErrMuon[0]);
        grUnfoldedMuon->SetLineColor(kBlack);
        grUnfoldedMuon->SetMarkerColor(kBlack);
        grUnfoldedMuon->SetMarkerStyle(25);
        grUnfoldedMuon->SetMarkerSize(1.5);
        grUnfoldedMuon->SetLineStyle(1);
        grUnfoldedMuon->SetLineWidth(1);
        grUnfoldedMuon->Draw("samepe");

	TGraphErrors *grMC = new TGraphErrors(5, &meanmass_mcElectron[0], &meanpt_mcElectron[0], &meanmasserr_mcElectron[0], &meanpterr_mcElectron[0]);
        grMC->SetLineColor(kRed);
        grMC->SetMarkerColor(kRed);
        grMC->SetMarkerStyle(20);
        grMC->SetMarkerSize(1.5);
        grMC->SetLineStyle(1);
        grMC->SetLineWidth(1);
        grMC->Draw("samepe");
	
        TF1 *f1 = new TF1("f1", "[0]+2.*[1]*log(x)", 40., 350.);
        f1->GetXaxis()->SetRangeUser(40., 350.);
        f1->SetLineColor(kBlack);
        grUnfolded->Fit(f1, "R"); // R: fitting sub range
        f1->Draw("same");

        TF1 *f2 = new TF1("f2", "[0]+2.*[1]*log(x)", 40., 350.);
        f2->GetXaxis()->SetRangeUser(40., 350.);
        f2->SetLineColor(kBlack);
        grUnfoldedMuon->Fit(f2, "R"); // R: fitting sub range
        f2->Draw("same");

        TLegend* leg_ = new TLegend(0.55, 0.25, 0.9, 0.45,"","brNDC");
        leg_->SetTextSize(0.03);
        //leg_->SetFillStyle(1);
        leg_->SetBorderSize(0);
        leg_->AddEntry(grUnfolded, "Unfolded electron data", "pe");
        leg_->AddEntry(grUnfoldedMuon, "Unfolded muon data", "pe");
        leg_->AddEntry(grMC, "DY MC at pre FSR (aMC@NLO)", "pe");
        leg_->Draw();

        c1->cd();

        CMS_lumi( pad1, 4, 11 );
        c1->cd();
        c1->SaveAs(outpdf);

}
