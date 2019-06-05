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

        TH2* hmcGenRec = NULL;
	if(sysName.CompareTo("Alt") == 0) hmcGenRec = (TH2*)filein->Get("hmc" + var + "GenRecnominal");
        else hmcGenRec = (TH2*)filein->Get("hmc" + var + "GenRec" + matrixName);

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
		if(postfix.CompareTo("Alt") == 0) hRec = (TH1*)filein->Get("h"+var+"Recnominal");
        	//hRec = (TH1*)filein->Get("h"+var+"Recnominal");
        	if( var.CompareTo("Pt") == 0 )   sysPtUnfold[postfix].at(nth)  ->SetInput(hRec,   bias); // 
        	if( var.CompareTo("Mass") == 0 ) sysMassUnfold[postfix].at(nth)->SetInput(hRec,   bias); // 
	}
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
                if(postfix.CompareTo("Alt") == 0 ) hRec = (TH1*)filein->Get("h"+var+"Recnominal");
        	//hRec = (TH1*)filein->Get("h"+var+"Recnominal");
        	if( var.CompareTo("Pt") == 0 )   sysPtUnfold[postfix].at(nth)  ->SubtractBackground(hRec, bkgName);
        	if( var.CompareTo("Mass") == 0 ) sysMassUnfold[postfix].at(nth)->SubtractBackground(hRec, bkgName);
	}
	
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

        //Double_t massBins[6] = {50., 65., 80., 100., 200., 350.};

        const TUnfoldBinningV17* temp_binning = nomPtUnfold->GetOutputBinning("Gen_Pt");
	const TVectorD* temp_tvecd = temp_binning->GetDistributionBinning(1);
	const Double_t* massBins = temp_tvecd->GetMatrixArray();


        TH1* hunfolded_mass =  nomMassUnfold->GetOutput("hunfolded_mass",0,0,"*[UO]",kTRUE);
        TH1 *histMCTruth_mass= nomMassUnfold->GetBias("histMCTruth_mass",0,0,"*[UO]",kTRUE);
        TH1 *histMCTruth_massAlt= sysMassUnfold["Alt"].at(0)->GetBias("histMCTruth_massAlt",0,0,"*[UO]",kTRUE);

        for(int ibin = 0; ibin < 5; ibin++){
                hunfolded_mass->GetXaxis()->  SetRange(hunfolded_mass->GetXaxis()->  FindBin(massBins[ibin]+0.01),hunfolded_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
                histMCTruth_mass->GetXaxis()->SetRange(histMCTruth_mass->GetXaxis()->FindBin(massBins[ibin]+0.01),histMCTruth_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
                histMCTruth_massAlt->GetXaxis()->SetRange(histMCTruth_mass->GetXaxis()->FindBin(massBins[ibin]+0.01),histMCTruth_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));

                meanMass_data.   push_back(hunfolded_mass->GetMean());
                meanMassStatErr_data.push_back(hunfolded_mass->GetMeanError());

                meanMass_mc.   push_back(histMCTruth_mass->GetMean());
                meanMassErr_mc.push_back(histMCTruth_mass->GetMeanError());

                meanMass_mcAlt.   push_back(histMCTruth_massAlt->GetMean());
                meanMassErr_mcAlt.push_back(histMCTruth_massAlt->GetMeanError());

	        std::map<TString, std::vector<TUnfoldDensityV17*>>::iterator it = sysMassUnfold.begin();
                std::map<TString, std::vector<Double_t>> temp_map; // temp map to save systematic results for a mass bin

           	while(it != sysMassUnfold.end()){
           	        int nSys = it->second.size();
           	        TH1* hdatasys_temp;
           	        for(int i = 0; i < nSys; i++){
           	             hdatasys_temp = sysMassUnfold[it->first].at(i)->GetOutput("hunfolded_mass_systemp",0,0,"*[UO]",kTRUE);
		             hdatasys_temp->GetXaxis()->SetRange(hdatasys_temp->GetXaxis()->FindBin(massBins[ibin]+0.01), hdatasys_temp->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
           	             temp_map[it->first].push_back(hdatasys_temp->GetMean());

           	             delete hdatasys_temp;
           	        }
           	        it++;
           	}
           	meanMass_sysdata.push_back(temp_map);
           	temp_map.clear();

        }

	delete hunfolded_mass;
	delete histMCTruth_mass;
	delete histMCTruth_massAlt;

        // check systematic mean values
        int size = meanMass_sysdata.size();
        for(int i = 0; i < size; i++){

           std::map<TString, Double_t> temp_map_;
           std::map<TString, std::vector<Double_t>>::iterator it = meanMass_sysdata.at(i).begin();
           while(it != meanMass_sysdata.at(i).end()){
                int size_ = it->second.size();

                TH1F *hpdfsys = NULL;
                if((it->first).CompareTo("PDFerror") == 0) hpdfsys = new TH1F("pdfsys", "pdfsys", 100, meanMass_data.at(i)-0.2, meanMass_data.at(i)+0.2); // temp histogram to contain PDF variations

                Double_t err = -999.; // 
                for(int j = 0; j < size_; j++){
                        if( (i==5 || i==7) && (it->first).CompareTo("Scale") == 0) continue;

                        if((it->first).CompareTo("PDFerror") == 0){
                                hpdfsys->Fill(it->second.at(j));
                        }

                        Double_t temp_err =  fabs(meanMass_data.at(i) - it->second.at(j));
                        if( temp_err > err){
                               err = temp_err;
                        }

                        //cout << i << " th mass bin, " << it->first << j << " th sys value: " << it->second.at(j) << endl;
                }
                if((it->first).CompareTo("PDFerror") == 0){
                        err = hpdfsys->GetRMS();
                        delete hpdfsys;
                }

                temp_map_[it->first] = err;
                it++;
           }
           meanMassErr_sysdata.push_back(temp_map_);
           temp_map_.clear();
        }// loop for mass bins

        size = meanMassErr_sysdata.size();
        for(int i = 0; i < size; i++){
           std::map<TString, Double_t>::iterator it = meanMassErr_sysdata.at(i).begin();
           Double_t totalSys = 0.;
           while(it != meanMassErr_sysdata.at(i).end()){
                //cout << i << " th mass bin, mass" << it->first << " " << it->second << endl;
                totalSys += pow(it->second, 2);
                it++;
           }

           cout << i << " th mass bin, total mass systematic uncertainty: " << sqrt(totalSys) << " statistical error: " << meanMassStatErr_data.at(i) << endl;
           meanMassSysErr_data.push_back(sqrt(totalSys));
	   meanMassTotErr_data.push_back(sqrt(totalSys + pow(meanMassStatErr_data.at(i),2)));
        }
}


// set mean pt from mass and DY mc
void ISRUnfold::setMeanPt(){

        int nMassBin = 5;

        for(int i = 0; i < nMassBin; i++){
           TString ibinMass;
           ibinMass.Form("%d", i);

           TH1* hpt_temp_data;
           TH1* hpt_temp_mc;
           TH1* hpt_temp_mcAlt;

           hpt_temp_data = nomPtUnfold->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
           hpt_temp_mc   = nomPtUnfold->GetBias("histMCTruth_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
           hpt_temp_mcAlt   = sysPtUnfold["Alt"].at(0)->GetBias("histMCTruth_pt_tempAlt",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);

           meanPt_data.   push_back(hpt_temp_data->GetMean());
           meanPtStatErr_data.push_back(hpt_temp_data->GetMeanError());

           meanPt_mc.   push_back(hpt_temp_mc->GetMean());
           meanPtErr_mc.push_back(hpt_temp_mc->GetMeanError());

           meanPt_mcAlt.   push_back(hpt_temp_mcAlt->GetMean());
           meanPtErr_mcAlt.push_back(hpt_temp_mcAlt->GetMeanError());

           std::map<TString, std::vector<TUnfoldDensityV17*>>::iterator it = sysPtUnfold.begin();
	   std::map<TString, std::vector<Double_t>> temp_map; // temp map to save systematic results for a mass bin

           while(it != sysPtUnfold.end()){
                   int nSys = it->second.size();
		   TH1* hdatasys_temp;
                   for(int i = 0; i < nSys; i++){
			hdatasys_temp = sysPtUnfold[it->first].at(i)->GetOutput("hunfolded_pt_systemp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
			temp_map[it->first].push_back(hdatasys_temp->GetMean());
	
			delete hdatasys_temp;
                   }
                   it++;
           }
	   meanPt_sysdata.push_back(temp_map);
	   temp_map.clear();

           delete hpt_temp_data;
	   delete hpt_temp_mc;
	   delete hpt_temp_mcAlt;
        }

	// check systematic mean values
	int size = meanPt_sysdata.size();
	for(int i = 0; i < size; i++){

           std::map<TString, Double_t> temp_map_;
           std::map<TString, std::vector<Double_t>>::iterator it = meanPt_sysdata.at(i).begin();
           while(it != meanPt_sysdata.at(i).end()){
		int size_ = it->second.size(); // size of systematic variations
		
		TH1F *hpdfsys = NULL;
		if((it->first).CompareTo("PDFerror") == 0) hpdfsys = new TH1F("pdfsys", "pdfsys", 100, meanPt_data.at(i)-0.2, meanPt_data.at(i)+0.2); // temp histogram to contain PDF variations

		Double_t err = -999.; // 
		for(int j = 0; j < size_; j++){
			if( (i==5 || i==7) && (it->first).CompareTo("Scale") == 0) continue;

			if((it->first).CompareTo("PDFerror") == 0){
	 			hpdfsys->Fill(it->second.at(j));
			}

                        Double_t temp_err =  fabs(meanPt_data.at(i) - it->second.at(j));
                        if( temp_err > err){
                               err = temp_err;
                        }
			//cout << i << " th mass bin, " << it->first << j << " th sys value: " << it->second.at(j) << endl; 
		}// loop for systematic variations
		if((it->first).CompareTo("PDFerror") == 0){
			err = hpdfsys->GetRMS();
			delete hpdfsys;
		}

		temp_map_[it->first] = err;
		it++;
	   }// loop for systematic sources
	   meanPtErr_sysdata.push_back(temp_map_);
	   temp_map_.clear();
	}// loop for mass binss

	size = meanPtErr_sysdata.size(); // size of mass bins
	for(int i = 0; i < size; i++){
           std::map<TString, Double_t>::iterator it = meanPtErr_sysdata.at(i).begin();
	   Double_t totalSys = 0.;
           while(it != meanPtErr_sysdata.at(i).end()){
		//cout << i << " th mass bin, pt" << it->first << " " << it->second << endl;
	 	totalSys += pow(it->second, 2);	
		it++;
	   }
	   
           cout << i << " th mass bin, total pt systematic uncertainty: " << sqrt(totalSys) << " statistical error: " << meanPtStatErr_data.at(i) << endl;
	   meanPtSysErr_data.push_back(sqrt(totalSys));
	   meanPtTotErr_data.push_back(sqrt(totalSys + pow(meanPtStatErr_data.at(i),2)));
	}// loop for mass bins
}

void ISRUnfold::drawISRresult(TString outpdf){

        gROOT->SetBatch();

        setTDRStyle();
        writeExtraText = true;       // if extra text
        extraText  = "work in progress";

        c1 = new TCanvas("c1","c1", 50, 50, 850, 700);
        c1->cd();
        gStyle->SetOptFit(0);

        c1->SetBottomMargin(0.2);
        c1->SetTopMargin(0.08);
        c1->SetTicks(1);
        c1->SetLogx();
        c1->SetGridy();

        TGraphErrors *grUnfolded = new TGraphErrors(5, &meanMass_data[0], &meanPt_data[0], &meanMassTotErr_data[0], &meanPtTotErr_data[0]);
        grUnfolded->SetLineColor(kBlack);
        grUnfolded->SetMarkerColor(kBlack);
        grUnfolded->SetMarkerStyle(20);
        grUnfolded->SetMarkerSize(1.);
        grUnfolded->SetLineStyle(1);
        grUnfolded->Draw("ape");
        grUnfolded->GetYaxis()->SetRangeUser(10.,30.);
        grUnfolded->GetXaxis()->SetLimits(30.,500.);
        grUnfolded->GetYaxis()->SetTitle("Average p_{T} (GeV)");
        grUnfolded->GetXaxis()->SetTitle("Average Mass (GeV)");

        TGraphErrors *grMC = new TGraphErrors(5, &meanMass_mc[0], &meanPt_mc[0], &meanMassErr_mc[0], &meanPtErr_mc[0]);
        grMC->SetLineColor(kRed);
        grMC->SetMarkerColor(kRed);
        grMC->SetMarkerStyle(20);
        grMC->SetMarkerSize(1.);
        grMC->SetLineStyle(1);
        grMC->Draw("pe same");

        TGraphErrors *grMCAlt = new TGraphErrors(5, &meanMass_mcAlt[0], &meanPt_mcAlt[0], &meanMassErr_mcAlt[0], &meanPtErr_mcAlt[0]);
        grMCAlt->SetLineColor(kRed);
        grMCAlt->SetMarkerColor(kRed);
        grMCAlt->SetMarkerStyle(24);
        grMCAlt->SetMarkerSize(1.);
        grMCAlt->SetLineStyle(1);
        grMCAlt->Draw("pe same");

        TF1 *f1 = new TF1("f1", "[0]+[1]*log(x)", 40., 350.);
        f1->GetXaxis()->SetRangeUser(40., 350.);
        f1->SetLineColor(kBlack);
        grUnfolded->Fit(f1, "R0"); // R: fitting sub range
        f1->Draw("same");

        CMS_lumi( c1, 4, 11 );
        c1->SaveAs(outpdf);
	delete grUnfolded;
	delete grMC;
	delete f1;
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

        leg_nom->AddEntry(hpt_temp_data, "Bkg subtracted data (mean: " + mean_nom + ")", "pl");
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
        		ratio->GetYaxis()->SetTitle("Systematic/ Nominal input");
        		ratio->GetXaxis()->SetTitle("p_{T} at pre FSR(GeV)");
        		ratio->SetMinimum(0.8);
        		ratio->SetMaximum(1.2);
        		ratio->SetTitle("");
        		ratio->GetXaxis()->SetTitleOffset(1.5);
        		ratio->GetYaxis()->SetNdivisions(515);
                	ratio->SetLineColor(kBlack);
                	ratio->SetLineStyle(2);
		}
                else{
		      ratio->Draw("histsame");
                      ratio->SetLineStyle(2);
		}

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

	//const Double_t massBins[6] = {50., 65., 80., 100., 200., 350.};

        const TUnfoldBinningV17* temp_binning = nomPtUnfold->GetOutputBinning("Gen_Pt");
        const TVectorD* temp_tvecd = temp_binning->GetDistributionBinning(1);
        const Double_t* massBins = temp_tvecd->GetMatrixArray();

        gROOT->SetBatch();

	setTDRStyle();
	writeExtraText = true;
	extraText  = "work in progress";

        TString ibinMass;
        ibinMass.Form("%d", nthMassBin);

        TH1* hunfolded_data;
	TH1* hunfolded_sys_err;
        TH1* hpreFSR_mc;
        TH1F *ratio;
        TH1F *ratio_sys_err;

	// get nominal unfoled result
	if(var.CompareTo("Pt") == 0 ){
        	hunfolded_data  = nomPtUnfold->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
		hunfolded_sys_err= ((TH1F*)hunfolded_data->Clone("sysErr")); 
		hpreFSR_mc   = nomPtUnfold->GetBias("histMCTruth_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
	}
	if(var.CompareTo("Mass") == 0 ){
        	hunfolded_data  = nomMassUnfold->GetOutput("hunfolded_mass_temp",0,0,"*[UO]",kTRUE);
		hunfolded_data->GetXaxis()->SetRange(hunfolded_data->GetXaxis()->FindBin(massBins[nthMassBin]+0.01),hunfolded_data->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));
		hunfolded_sys_err= ((TH1F*)hunfolded_data->Clone("sysErr")); 
		hpreFSR_mc   = nomMassUnfold->GetBias("histMCTruth_mass_temp",0,0,"*[UO]",kTRUE);
	}

        ratio= ((TH1F*)hunfolded_data->Clone("ratio"));
        ratio_sys_err= ((TH1F*)hunfolded_data->Clone("ratio_sys"));
        ratio->Divide(hpreFSR_mc);
        ratio_sys_err->Divide(hpreFSR_mc);

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

        hunfolded_data->SetTitle("");
        hunfolded_data->Draw("p9histe");
        hpreFSR_mc->Draw("histsame");
        hunfolded_data->SetMarkerStyle(20);
        hunfolded_data->SetMarkerSize(.7);
        hunfolded_data->SetLineColor(kBlack);
        hpreFSR_mc->SetLineColor(kRed);
        hunfolded_data->GetYaxis()->SetTitle("Events/bin");
        hunfolded_data->SetMinimum(10.);

	TH1* hpt_sys_temp;
	int sysSize = sysPtUnfold[sysName].size();
        for(int i = 0; i < sysSize; i++){
                if((i==5 || i==7) && sysName.CompareTo("Scale") == 0) continue;

                TString isys;
                isys.Form("%d", i);

                TH1 * hdatasys_temp;
                TH1 * hmcsys_temp = NULL;
       
                if(var.CompareTo("Pt") == 0 )   hdatasys_temp = sysPtUnfold[sysName].at(i)->GetOutput("hunfolded_pt_systemp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
		if(var.CompareTo("Pt") == 0 && sysName.CompareTo("Alt") == 0 ){
	 		hmcsys_temp = sysPtUnfold[sysName].at(i)->GetBias("histMCTruth_pt_tempAlt",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE); 
			hmcsys_temp->SetDirectory(0);
		}

                if(var.CompareTo("Mass") == 0 ){
              		hdatasys_temp = sysMassUnfold[sysName].at(i)->GetOutput("hunfolded_mass_systemp",0,0,"*[UO]",kTRUE);
             		hdatasys_temp->GetXaxis()->SetRange(hdatasys_temp->GetXaxis()->FindBin(massBins[nthMassBin]+0.01), hdatasys_temp->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));
                }

                for(int ibin = 1; ibin<hunfolded_sys_err->GetNbinsX()+1;ibin++){
                   Double_t err = -999.;
                   Double_t ratio_err = -999.;

                   TH1F * ratio_temp = ((TH1F*)hdatasys_temp->Clone("ratio"));

                   if(hdatasys_temp->GetBinContent(ibin) < 0 ){ 
                	cout << " negative bin exist: " << ibin << " value: " << hdatasys_temp->GetBinContent(ibin) << " prv bin value: " << hdatasys_temp->GetBinContent(ibin - 1) << " next bin value: " << hdatasys_temp->GetBinContent(ibin + 1) <<std::endl;
		   }

                   ratio_temp->Divide(hpreFSR_mc);
                   // get "envelope"
                   Double_t temp_err =  fabs(hunfolded_data->GetBinContent(ibin) - hdatasys_temp->GetBinContent(ibin));
                   Double_t temp_sys_err = fabs(ratio->GetBinContent(ibin) - ratio_temp->GetBinContent(ibin));

                   if( temp_err > err){
                             err = temp_err;
                   }
                   if( temp_sys_err > ratio_err){
                             ratio_err = temp_sys_err;
                   }
                   hunfolded_sys_err->SetBinError(ibin, err);
                   ratio_sys_err->SetBinContent(ibin, 1.);
                   ratio_sys_err->SetBinError(ibin, ratio_err);

                   delete ratio_temp;
              }// loop for bin contents

                hpt_sys_temp = ((TH1F*)hdatasys_temp->Clone("pt_temp"));
                hpt_sys_temp->Draw("histsame");
                hpt_sys_temp->SetLineColor(kBlack);
                hpt_sys_temp->SetLineStyle(2);
		if(hmcsys_temp!=NULL){
			cout << "draw " << sysName << " MC histogram " << " 1 bin content: " << hmcsys_temp->GetBinContent(1) << endl;
			hmcsys_temp->Draw("histsame");
			hmcsys_temp->SetLineColor(kBlue);
		}


              delete hdatasys_temp;
	      //delete hmcsys_temp;
          }

	hunfolded_sys_err->Draw("E2same");
        hunfolded_sys_err->SetMarkerSize(0);
        hunfolded_sys_err->SetFillColorAlpha(kBlack,0.3);

        TString mean_nom;
        mean_nom.Form("%.5f", hunfolded_data->GetMean());

        TLegend* leg_nom = new TLegend(0.45, 0.70, 0.75, 0.9,"","brNDC");
        leg_nom->SetNColumns(2);
        leg_nom->SetTextSize(0.055);
        leg_nom->SetFillStyle(0);
        leg_nom->SetBorderSize(0);

        leg_nom->AddEntry(hunfolded_data, "Unfolded data (mean: " + mean_nom + ")", "pl");
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
        c1->SaveAs(outpdf+"_"+ibinMass+"_"+var+".pdf");

        delete hunfolded_data;
	delete hunfolded_sys_err;
        delete hpreFSR_mc;
	delete hpt_sys_temp;
        delete pad1;
        delete pad2;
        delete c1;
}
