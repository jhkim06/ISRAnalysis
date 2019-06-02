#include "ISR_unfoldUtils.h"

#include "tdrstyle.C"
#include "CMS_lumi.C"

void ISRUnfold::setNomTUnfoldDensity(TString filepath, TString var, TString matrixName, bool isfsr){

	TFile* filein = new TFile(filepath);

        TH2* hmcGenRec = (TH2*)filein->Get("hmc" + var + "GenRec" + matrixName);
        TUnfoldBinning* binning_Rec = NULL;
        TUnfoldBinning* binning_Gen = NULL;

        if( var.CompareTo("Pt") == 0 ){

          if(!isfsr) binning_Rec = (TUnfoldBinning*)filein->Get("Rec_Pt");
          else binning_Rec = (TUnfoldBinning*)filein->Get("Gen_Pt");

          binning_Gen = (TUnfoldBinning*)filein->Get("Gen_Pt");
        }
        if( var.CompareTo("Mass") == 0 ){

          if(!isfsr) binning_Rec = (TUnfoldBinning*)filein->Get("Rec_Mass");
          else binning_Rec = (TUnfoldBinning*)filein->Get("Gen_Mass");

          binning_Gen = (TUnfoldBinning*)filein->Get("Gen_Mass");

        }

	if( var.CompareTo("Pt") == 0 ){ 
        	nomPtUnfold = new TUnfoldDensityV17(hmcGenRec,
        	                               TUnfold::kHistMapOutputHoriz,
        	                               TUnfold::kRegModeNone, // fixed to use no regularisation temporary
        	                               TUnfold::kEConstraintArea,
        	                               TUnfoldDensityV17::kDensityModeBinWidth,
        	                               binning_Gen,binning_Rec);
	}

        if( var.CompareTo("Mass") == 0 ){
                nomMassUnfold = new TUnfoldDensityV17(hmcGenRec,
                                               TUnfold::kHistMapOutputHoriz,
                                               TUnfold::kRegModeNone, // fixed to use no regularisation temporary
                                               TUnfold::kEConstraintArea,
                                               TUnfoldDensityV17::kDensityModeBinWidth,
                                               binning_Gen,binning_Rec);
        }
}

void ISRUnfold::setSysTUnfoldDensity(TString filepath, TString var, TString sysName, int nth, bool isfsr){

        TFile* filein = new TFile(filepath);
	TString nth_;
	nth_.Form("%d", nth);

	TString matrixName = sysName + "_" + nth_;

        TH2* hmcGenRec = (TH2*)filein->Get("hmc" + var + "GenRec" + matrixName);
        //TH2* hmcGenRec = (TH2*)filein->Get("hmc" + var + "GenRecnominal");
        TUnfoldBinning* binning_Rec = NULL;
        TUnfoldBinning* binning_Gen = NULL;

        if( var.CompareTo("Pt") == 0 ){

          if(!isfsr) binning_Rec = (TUnfoldBinning*)filein->Get("Rec_Pt");
          else binning_Rec = (TUnfoldBinning*)filein->Get("Gen_Pt");

          binning_Gen = (TUnfoldBinning*)filein->Get("Gen_Pt");
        }
        if( var.CompareTo("Mass") == 0 ){

          if(!isfsr) binning_Rec = (TUnfoldBinning*)filein->Get("Rec_Mass");
          else binning_Rec = (TUnfoldBinning*)filein->Get("Gen_Mass");

          binning_Gen = (TUnfoldBinning*)filein->Get("Gen_Mass");

        }

        if( var.CompareTo("Pt") == 0 ){
                sysPtUnfold[sysName].push_back( new TUnfoldDensityV17(hmcGenRec,
                                               TUnfold::kHistMapOutputHoriz,
                                               TUnfold::kRegModeNone, // fixed to use no regularisation temporary
                                               TUnfold::kEConstraintArea,
                                               TUnfoldDensityV17::kDensityModeBinWidth,
                                               binning_Gen,binning_Rec));
        }

        if( var.CompareTo("Mass") == 0 ){
                sysMassUnfold[sysName].push_back( new TUnfoldDensityV17(hmcGenRec,
                                               TUnfold::kHistMapOutputHoriz,
                                               TUnfold::kRegModeNone, // fixed to use no regularisation temporary
                                               TUnfold::kEConstraintArea,
                                               TUnfoldDensityV17::kDensityModeBinWidth,
                                               binning_Gen,binning_Rec));
        }
}

void ISRUnfold::setInput(TString var, TString postfix, TString filepath, int nth, bool isSys, double bias){
	// No effects on the unfolded results respect to bias factor 

	TFile* filein = new TFile(filepath);
        TString nth_;
        nth_.Form("%d", nth);
        TH1* hRec;

	if(!isSys){
        	hRec = (TH1*)filein->Get("h"+var+"Rec"+postfix);
        	if( var.CompareTo("Pt") == 0 )   nomPtUnfold->SetInput(hRec,   bias); // 
        	if( var.CompareTo("Mass") == 0 ) nomMassUnfold->SetInput(hRec, bias); // 
	}
	else{
        	hRec = (TH1*)filein->Get("h"+var+"Rec"+postfix+"_"+nth_);
                hRec->SetDirectory(0);
        	if( var.CompareTo("Pt") == 0 )   sysPtUnfold[postfix].at(nth)  ->SetInput(hRec,   bias); // 
        	if( var.CompareTo("Mass") == 0 ) sysMassUnfold[postfix].at(nth)->SetInput(hRec,   bias); // 
	}
	delete hRec;
	filein->Close();
	delete filein;

}

void ISRUnfold::subBkgs(TString var, TString postfix, TString filepath, TString bkgName, int nth, bool isSys){

	TFile* filein = new TFile(filepath);
        TString nth_;
        nth_.Form("%d", nth);
        TH1* hRec = NULL;

	if(!isSys){
        	hRec = (TH1*)filein->Get("h"+var+"Rec"+postfix);
        	if( var.CompareTo("Pt") == 0 )   nomPtUnfold->  SubtractBackground(hRec, bkgName);
        	if( var.CompareTo("Mass") == 0 ) nomMassUnfold->SubtractBackground(hRec, bkgName);
	}
	else{	
        	hRec = (TH1*)filein->Get("h"+var+"Rec"+postfix+"_"+nth_);
		hRec->SetDirectory(0);
        	if( var.CompareTo("Pt") == 0 )   sysPtUnfold[postfix].at(nth)  ->SubtractBackground(hRec, bkgName);
        	if( var.CompareTo("Mass") == 0 ) sysMassUnfold[postfix].at(nth)->SubtractBackground(hRec, bkgName);
	}
	
	delete hRec;
	filein->Close();
	delete filein;

}

void ISRUnfold::doISRUnfold(){

	nomPtUnfold->DoUnfold(0);
	nomMassUnfold->DoUnfold(0);

	std::map<TString, std::vector<TUnfoldDensityV17*>>::iterator it = sysPtUnfold.begin();

	while(it != sysPtUnfold.end()){
		int nSys = it->second.size();
		for(int i = 0; i < nSys; i++){
			it->second.at(i)->DoUnfold(0);
		}
		it++;
	}

	it = sysMassUnfold.begin();
        while(it != sysMassUnfold.end()){
                int nSys = it->second.size();
                for(int i = 0; i < nSys; i++){
                        it->second.at(i)->DoUnfold(0);
                }
                it++;
        }

}

void ISRUnfold::setMeanMass(){

        Double_t massBins[6] = {50., 65., 80., 100., 200., 350.};

        TH1* hunfolded_mass =  nomMassUnfold->GetOutput("hunfolded_mass",0,0,"*[UO]",kTRUE);
        TH1 *histMCTruth_mass= nomMassUnfold->GetBias("histMCTruth_mass",0,0,"*[UO]",kTRUE);

        for(int ibin = 0; ibin < 5; ibin++){
                hunfolded_mass->GetXaxis()->  SetRange(hunfolded_mass->GetXaxis()->  FindBin(massBins[ibin]+0.01),hunfolded_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
                histMCTruth_mass->GetXaxis()->SetRange(histMCTruth_mass->GetXaxis()->FindBin(massBins[ibin]+0.01),histMCTruth_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));

                meanMass_data.   push_back(hunfolded_mass->GetMean());
                meanMassErr_data.push_back(hunfolded_mass->GetMeanError());

                meanMass_mc.   push_back(histMCTruth_mass->GetMean());
                meanMassErr_mc.push_back(histMCTruth_mass->GetMeanError());
        }

	delete hunfolded_mass;
	delete histMCTruth_mass;
}


// set mean pt from mass and DY mc
void ISRUnfold::setMeanPt(){

        int nMassBin = 5;

        for(int i = 0; i < nMassBin; i++){
           TString ibinMass;
           ibinMass.Form("%d", i);

           TH1* hpt_temp_data;
           TH1* hpt_temp_mc;

           hpt_temp_data = nomPtUnfold->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
           hpt_temp_mc   = nomPtUnfold->GetBias("histMCTruth_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);

           meanPt_data.   push_back(hpt_temp_data->GetMean());
           meanPtErr_data.push_back(hpt_temp_data->GetMeanError());

           meanPt_mc.   push_back(hpt_temp_mc->GetMean());
           meanPtErr_mc.push_back(hpt_temp_mc->GetMeanError());

           delete hpt_temp_data;
	   delete hpt_temp_mc;
        }
}


void ISRUnfold::drawInputPlots(TString outpdf, TString var, int nthMassBin, TString sysName){

        gROOT->SetBatch();

	setTDRStyle();
	writeExtraText = true;
	extraText  = "work in progress";

        TString ibinMass;
        ibinMass.Form("%d", nthMassBin);

        TH1* hpt_temp_data;
        TH1F *ratio;

        hpt_temp_data   = nomPtUnfold->GetInput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);

        c1=new TCanvas("c1", "c1", 50, 50, 800, 800);
        c1->cd();
        gStyle->SetOptStat(0);

        TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
        pad1->SetBottomMargin(0.01);
        pad1->SetTopMargin(0.1);
        pad1->SetTicks(1);
        pad1->SetLogy();
        pad1->Draw();
        pad1->cd();

        hpt_temp_data->SetTitle("");
        hpt_temp_data->Draw("p9histe");
        hpt_temp_data->SetMarkerStyle(20);
        hpt_temp_data->SetMarkerSize(.7);
        hpt_temp_data->SetLineColor(kBlack);
        hpt_temp_data->GetYaxis()->SetTitle("Events/bin");


	TH1* hpt_sys_temp;
	int sysSize = sysPtUnfold[sysName].size();
        Double_t err = -999.;
        Double_t ratio_err = -999.;
        for(int i = 0; i < sysSize; i++){
                if((i==5 || i==7) && sysName.CompareTo("Scale") == 0) continue;

                TString isys;
                isys.Form("%d", i);

                TH1 * hsyspt_temp = sysPtUnfold[sysName].at(i)->GetInput("hunfolded_pt_systemp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
	        hpt_sys_temp = ((TH1F*)hsyspt_temp->Clone("pt_temp"));
	        hpt_sys_temp->Draw("histsame");
	        hpt_sys_temp->SetLineColor(kBlack);
	        hpt_sys_temp->SetLineStyle(2);

                delete hsyspt_temp;
        }

        TString mean_nom;
        mean_nom.Form("%.5f", hpt_temp_data->GetMean());

        TLegend* leg_nom = new TLegend(0.45, 0.70, 0.75, 0.9,"","brNDC");
        leg_nom->SetNColumns(2);
        leg_nom->SetTextSize(0.055);
        leg_nom->SetFillStyle(0);
        leg_nom->SetBorderSize(0);

        leg_nom->AddEntry(hpt_temp_data, "Unfolded data (mean: " + mean_nom + ")", "pl");
        leg_nom->Draw();

	c1->cd();

        TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.4);
        pad2->SetTopMargin(0.05);
        pad2->SetBottomMargin(0.2);
        pad2->SetTicks(1);
        pad2->SetGridy(1);
        pad2->Draw();

        pad2->cd();

        for(int i = 0; i < sysSize; i++){
                if((i==5 || i==7) && sysName.CompareTo("Scale") == 0) continue;

                TString isys;
                isys.Form("%d", i); 

                TH1 * hsyspt_temp = sysPtUnfold[sysName].at(i)->GetInput("hunfolded_pt_systemp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
                ratio = ((TH1F*)hsyspt_temp->Clone("pt_temp"));
		ratio->Divide(hpt_temp_data);
                if(i==0 ){
 			ratio->Draw("hist");
        		ratio->GetYaxis()->SetTitle("Data/ MC");
        		ratio->GetXaxis()->SetTitle("p_{T} at pre FSR(GeV)");
        		ratio->SetMinimum(0.8);
        		ratio->SetMaximum(1.2);
        		ratio->SetTitle("");
        		ratio->GetXaxis()->SetTitleOffset(1.5);
        		ratio->GetYaxis()->SetNdivisions(515);
                	ratio->SetLineColor(kBlack);
                	ratio->SetLineStyle(2);
		}
                else ratio->Draw("histsame");

                delete hsyspt_temp;
        }   

	CMS_lumi( c1, 4, 0 );
        c1->cd();
        c1->SaveAs(outpdf+"_input_"+ibinMass+".pdf");

        delete hpt_temp_data;
	delete ratio;
        delete pad1;
        delete pad2;
        delete c1;
}


void ISRUnfold::drawNominalPlots(TString outpdf, TString var, int nthMassBin, TString sysName){

        gROOT->SetBatch();

	setTDRStyle();
	writeExtraText = true;
	extraText  = "work in progress";

        TString ibinMass;
        ibinMass.Form("%d", nthMassBin);

        TH1* hpt_temp_data;
	TH1* hpt_sys_err_temp;
        TH1* hpt_temp_mc;
        TH1F *ratio;
        TH1F *ratio_sys_err;

	// get nominal unfoled result
        hpt_temp_data   = nomPtUnfold->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
	hpt_sys_err_temp= ((TH1F*)hpt_temp_data->Clone("sysErr")); 
	hpt_temp_mc   = nomPtUnfold->GetBias("histMCTruth_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);

        ratio= ((TH1F*)hpt_temp_data->Clone("ratio"));
        ratio_sys_err= ((TH1F*)hpt_temp_data->Clone("ratio_sys"));
        ratio->Divide(hpt_temp_mc);
        ratio_sys_err->Divide(hpt_temp_mc);

        c1=new TCanvas("c1", "c1", 50, 50, 800, 800);
        c1->cd();
        gStyle->SetOptStat(0);

        TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
        pad1->SetBottomMargin(0.01);
        pad1->SetTopMargin(0.1);
        pad1->SetTicks(1);
        pad1->SetLogy();
        pad1->Draw();
        pad1->cd();

        hpt_temp_data->SetTitle("");
        hpt_temp_data->Draw("p9histe");
        hpt_temp_mc->Draw("histsame");
        hpt_temp_data->SetMarkerStyle(20);
        hpt_temp_data->SetMarkerSize(.7);
        hpt_temp_data->SetLineColor(kBlack);
        hpt_temp_mc->SetLineColor(kRed);
        hpt_temp_data->GetYaxis()->SetTitle("Events/bin");

	TH1* hpt_sys_temp;
        for(int ibin = 1; ibin<hpt_sys_err_temp->GetNbinsX()+1;ibin++){
	   int sysSize = sysPtUnfold[sysName].size();
           Double_t err = -999.;
           Double_t ratio_err = -999.;
           for(int i = 0; i < sysSize; i++){
                   if((i==5 || i==7) && sysName.CompareTo("Scale") == 0) continue;

                   TString isys;
                   isys.Form("%d", i);

                   TH1 * hsyspt_temp = sysPtUnfold[sysName].at(i)->GetOutput("hunfolded_pt_systemp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
		   hpt_sys_temp = ((TH1F*)hsyspt_temp->Clone("pt_temp"));
		   hpt_sys_temp->Draw("histsame");
		   hpt_sys_temp->SetLineColor(kBlack);
		   hpt_sys_temp->SetLineStyle(2);
		   TH1F * ratio_temp = ((TH1F*)hsyspt_temp->Clone("ratio"));
		   ratio_temp->Divide(hpt_temp_mc);
                   // get "envelope"
                   Double_t temp_err =  fabs(hpt_temp_data->GetBinContent(ibin) - hsyspt_temp->GetBinContent(ibin));
		   Double_t temp_sys_err = fabs(ratio->GetBinContent(ibin) - ratio_temp->GetBinContent(ibin));
                   if( temp_err > err){
                          err = temp_err;
                   }
	           if( temp_sys_err > ratio_err){
                          ratio_err = temp_sys_err;
                   }
                   delete hsyspt_temp;
                   delete ratio_temp;
           }
           hpt_sys_err_temp->SetBinError(ibin, err);
	   ratio_sys_err->SetBinContent(ibin, 1.);
           ratio_sys_err->SetBinError(ibin, ratio_err);
        }

	hpt_sys_err_temp->Draw("E2same");
        hpt_sys_err_temp->SetMarkerSize(0);
        hpt_sys_err_temp->SetFillColorAlpha(kBlack,0.3);

        TString mean_nom;
        mean_nom.Form("%.5f", hpt_temp_data->GetMean());

        TLegend* leg_nom = new TLegend(0.45, 0.70, 0.75, 0.9,"","brNDC");
        leg_nom->SetNColumns(2);
        leg_nom->SetTextSize(0.055);
        leg_nom->SetFillStyle(0);
        leg_nom->SetBorderSize(0);

        leg_nom->AddEntry(hpt_temp_data, "Unfolded data (mean: " + mean_nom + ")", "pl");
        leg_nom->Draw();

	c1->cd();

        TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.4);
        pad2->SetTopMargin(0.05);
        pad2->SetBottomMargin(0.2);
        pad2->SetTicks(1);
        pad2->SetGridy(1);
        pad2->Draw();

        pad2->cd();
        ratio->Draw("hist");
        ratio->GetYaxis()->SetTitle("Data/ MC");
        ratio->GetXaxis()->SetTitle("p_{T} at pre FSR(GeV)");
        ratio->SetMinimum(0.8);
        ratio->SetMaximum(1.2);
        ratio->SetTitle("");
        ratio->GetXaxis()->SetTitleOffset(1.5);
        ratio->GetYaxis()->SetNdivisions(515);

	ratio_sys_err->Draw("E2same");
        ratio_sys_err->SetMarkerSize(0);
        ratio_sys_err->SetFillColorAlpha(kBlack,0.3);

	CMS_lumi( c1, 4, 0 );
        c1->cd();
        c1->SaveAs(outpdf+"_"+ibinMass+".pdf");

        delete hpt_temp_data;
        delete hpt_temp_mc;
	delete hpt_sys_temp;
        delete pad1;
        delete pad2;
        delete c1;
}
