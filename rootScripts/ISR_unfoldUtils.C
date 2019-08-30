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
	if(sysName.CompareTo("Alt") == 0 || sysName.CompareTo("unfoldBias") == 0 || sysName.CompareTo("unfoldScan") == 0 || sysName.CompareTo("Closure") == 0){
            if(sysName.CompareTo("Alt") == 0){
                if(var.CompareTo("Pt") == 0) hmcGenRec = (TH2*)filein->Get("eventSel/ptll TUnfold matrix/hmc" + var + "GenRecnominal");
                if(var.CompareTo("Mass") == 0) hmcGenRec = (TH2*)filein->Get("eventSel/mll TUnfold matrix/hmc" + var + "GenRecnominal");
            }
            else
                hmcGenRec = (TH2*)filein->Get("hmc" + var + "GenRecnominal");
        }
        else hmcGenRec = (TH2*)filein->Get("hmc" + var + "GenRec" + matrixName);

        //TH2* hmcGenRec = (TH2*)filein->Get("hmc" + var + "GenRecnominal");
        TUnfoldBinning* binning_Rec = NULL;
        TUnfoldBinning* binning_Gen = NULL;

        if( var.CompareTo("Pt") == 0 ){

          TString Rec_Pt = "Rec_Pt";
          TString Gen_Pt = "Gen_Pt";

          if (sysName.CompareTo("Alt") == 0){

            Rec_Pt = "eventSel/ptll TUnfold matrix/" + Rec_Pt;
            Gen_Pt = "eventSel/ptll TUnfold matrix/" + Gen_Pt;

          }

          if(!isfsr) binning_Rec = (TUnfoldBinning*)filein->Get(Rec_Pt);
          else binning_Rec = (TUnfoldBinning*)filein->Get(Rec_Pt);

          binning_Gen = (TUnfoldBinning*)filein->Get(Gen_Pt);
        }
        if( var.CompareTo("Mass") == 0 ){

          TString Rec_Mass = "Rec_Mass";
          TString Gen_Mass = "Gen_Mass";

          if (sysName.CompareTo("Alt") == 0){

            Rec_Mass = "eventSel/mll TUnfold matrix/" + Rec_Mass;
            Gen_Mass = "eventSel/mll TUnfold matrix/" + Gen_Mass;

          }

          if(!isfsr) binning_Rec = (TUnfoldBinning*)filein->Get(Rec_Mass);
          else binning_Rec = (TUnfoldBinning*)filein->Get(Rec_Mass);

          binning_Gen = (TUnfoldBinning*)filein->Get(Gen_Mass);

        }

	TUnfold::ERegMode mode=TUnfold::kRegModeNone;
	if( sysName.CompareTo("unfoldScan") == 0 || sysName.CompareTo("unfoldBias") == 0) mode = TUnfold::kRegModeCurvature;

        if( var.CompareTo("Pt") == 0 ){
                sysPtUnfold[sysName].push_back( new TUnfoldDensityV17(hmcGenRec,
                                               TUnfold::kHistMapOutputHoriz,
                                               mode, // fixed to use no regularisation temporary
                                               TUnfold::kEConstraintArea,
                                               TUnfoldDensityV17::kDensityModeBinWidth,
                                               binning_Gen,binning_Rec));
        }

        if( var.CompareTo("Mass") == 0 ){
                sysMassUnfold[sysName].push_back( new TUnfoldDensityV17(hmcGenRec,
                                               TUnfold::kHistMapOutputHoriz,
                                               mode, // fixed to use no regularisation temporary
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
		if(postfix.CompareTo("Alt") == 0 || postfix.CompareTo("unfoldBias") == 0 || postfix.CompareTo("unfoldScan") == 0 || postfix.CompareTo("Closure") == 0) hRec = (TH1*)filein->Get("h"+var+"Recnominal");
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
                if(postfix.CompareTo("Alt") == 0 || postfix.CompareTo("unfoldBias") == 0 || postfix.CompareTo("unfoldScan") == 0) hRec = (TH1*)filein->Get("h"+var+"Recnominal");
        	//hRec = (TH1*)filein->Get("h"+var+"Recnominal");
        	if( var.CompareTo("Pt") == 0 )   sysPtUnfold[postfix].at(nth)  ->SubtractBackground(hRec, bkgName);
        	if( var.CompareTo("Mass") == 0 ) sysMassUnfold[postfix].at(nth)->SubtractBackground(hRec, bkgName);
	}
	
	filein->Close();
	delete filein;

}

void ISRUnfold::drawLCurve(TString outpdf, TString var){

        gROOT->SetBatch();

	TGraph *lCurve_temp;
	Int_t iBest_temp = 0;
	if( var.CompareTo("Pt") == 0){
		lCurve_temp = lCurve;
		iBest_temp = iBest;
	}

	else if( var.CompareTo("Mass") == 0){
                lCurve_temp = lCurve_mass;
                iBest_temp = iBest_mass;
        }
	
	

	Double_t x[1],y[1];
	lCurve_temp->GetPoint(iBest_temp,x[0],y[0]);
	TGraph *bestLCurve=new TGraph(1,x,y);

        c1 = new TCanvas("c1","c1", 50, 50, 850, 700);
        c1->cd();
        gStyle->SetOptFit(0);

        c1->SetBottomMargin(0.2);
        c1->SetTopMargin(0.08);
        c1->SetTicks(1);

	lCurve_temp->Draw("AL");
  	bestLCurve->SetMarkerColor(kRed);
  	bestLCurve->Draw("*");
	c1->SaveAs(outpdf);	
	delete c1;
	

}

void ISRUnfold::drawRhoLog(TString outpdf, TString var){

        gROOT->SetBatch();

	TSpline *rhoLogTau_temp;
        Int_t iBest_temp = 0;
	Int_t nScan_temp = 0;
        if( var.CompareTo("Pt") == 0){
                rhoLogTau_temp = rhoLogTau;
                iBest_temp = iBest;
		nScan_temp = nScan;
        }

        else if( var.CompareTo("Mass") == 0){
                rhoLogTau_temp = rhoLogTau_mass;
                iBest_temp = iBest_mass;
		nScan_temp = nScan_mass;
        }


        Double_t t[1],rho[1];
        rhoLogTau_temp->GetKnot(iBest_temp,t[0],rho[0]);
	TGraph *bestRhoLogTau=new TGraph(1,t,rho);

  	Double_t *tAll=new Double_t[nScan_temp],*rhoAll=new Double_t[nScan_temp];
  	for(Int_t i=0;i<nScan_temp;i++) {
  	   rhoLogTau_temp->GetKnot(i,tAll[i],rhoAll[i]);
	   cout << i << " tAll[i]: " << tAll[i] << " rhoAll[i]: " << rhoAll[i] << endl;
  	}
  	TGraph *knots=new TGraph(nScan_temp,tAll,rhoAll);

        c1 = new TCanvas("c1","c1", 50, 50, 850, 700);
        c1->cd();
        gStyle->SetOptFit(0);

        c1->SetBottomMargin(0.2);
        c1->SetTopMargin(0.08);
        c1->SetTicks(1);

  	//rhoLogTau_temp->Draw();
  	knots->Draw("AP");
  	bestRhoLogTau->SetMarkerColor(kRed);
  	bestRhoLogTau->Draw("*");

        c1->SaveAs(outpdf);
        delete c1;


}

void ISRUnfold::doISRUnfold(bool doSys){

	double tauMin=1.e-4;
	double tauMax=1.e-3;

  	nScan=30;
  	rhoLogTau=0;
  	lCurve=0;
	iBest = 0;

        nScan_mass=30;
        rhoLogTau_mass=0;
        lCurve_mass=0;
        iBest_mass = 0;

	nomPtUnfold->DoUnfold(0);
	nomMassUnfold->DoUnfold(0);


        if(doSys){
	    std::map<TString, std::vector<TUnfoldDensityV17*>>::iterator it = sysPtUnfold.begin();

	    while(it != sysPtUnfold.end()){
	    	int nSys = it->second.size();
	    	for(int i = 0; i < nSys; i++){
	    		if((it->first).CompareTo("unfoldScan") == 0 || (it->first).CompareTo("unfoldBias") == 0){
	    			//it->second.at(i)->ScanLcurve(50,tauMin,tauMax,0);
  	    			iBest=it->second.at(i)->ScanTau(nScan,0.,0.,&rhoLogTau,
  	    			                           TUnfoldDensity::kEScanTauRhoAvgSys,
  	    			                           0,0,
  	    			                           &lCurve);
	    		}
	    		else it->second.at(i)->DoUnfold(0);

	    	}
	    	it++;
	    }

	    it = sysMassUnfold.begin();
            while(it != sysMassUnfold.end()){
                    int nSys = it->second.size();
                    for(int i = 0; i < nSys; i++){
                            if((it->first).CompareTo("unfoldScan") == 0 || (it->first).CompareTo("unfoldBias") == 0){
	    			//it->second.at(i)->ScanLcurve(50,tauMin,tauMax,0);

                                    iBest_mass=it->second.at(i)->ScanTau(nScan_mass,0.,0.,&rhoLogTau_mass,
                                                               TUnfoldDensity::kEScanTauRhoAvgSys,
                                                               0,0,
                                                               &lCurve_mass);
	    		}
	    		else it->second.at(i)->DoUnfold(0);

                    }
                    it++;
            }
        }

}

double ISRUnfold::Chi2Test(TH1 *data, TH1 *mc){ // check if this is right way to get chi square value

 gROOT->SetBatch();

 double chi2 = 0.;
 double ndf  = 0.;

 for(int i=1;i<=mc->GetNbinsX();i++) {
    ndf += 1.;
    if(data->GetBinError(i) == 0){
      std::cout << "unfolded " << i << " bin: " << data->GetBinContent(i) << std::endl;
      std::cout << "error is zero in the " << i << " bin..." << std::endl;
      std::cout << "so skip this bin" << std::endl;
      continue;
    }
    double pull=(data->GetBinContent(i)-mc->GetBinContent(i))/data->GetBinError(i);
    chi2+= pull*pull;
 }
 return chi2/ndf;

}


// need unfolded hist, rho matrix (GetRhoIJtotal), MC truth
double ISRUnfold::DoFit(TString var, int nthMassBin){

	TH1 *g_fcnHist=0;
	TMatrixD *g_fcnMatrix=0;

	TH1* hpt_temp_data;
	TH1* hpt_temp_mc;
	TH2 *hrho;

        TString ibinMass;
        ibinMass.Form("%d", nthMassBin);

	hpt_temp_data = nomPtUnfold->GetOutput("hunfolded_pt_temp_chi",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
	hpt_temp_mc   = nomPtUnfold->GetBias("histMCTruth_pt_temp_chi",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
	hrho          = nomPtUnfold->GetRhoIJtotal("histRho_chi", 0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE); 
	
	//cout << "rho matrix dimention: nbinx " << hrho->GetXaxis()->GetNbins() << " nbiny: " << hrho->GetXaxis()->GetNbins() << endl;
	
        int n=hpt_temp_data->GetNbinsX(); 

        TMatrixDSym v(n);
        //g_fcnMatrix=new TMatrixD(n,n);
        for(int i=0;i<n;i++) {
           for(int j=0;j<n;j++) {
              v(i,j)=hrho->GetBinContent(i+1,j+1)*(hpt_temp_data->GetBinError(i+1)*hpt_temp_data->GetBinError(j+1));
           }
        }

        TMatrixDSymEigen ev(v);
        TMatrixD d(n,n);
        TVectorD di(ev.GetEigenValues());
        for(int i=0;i<n;i++) {
           if(di(i)>0.0) {
              d(i,i)=1./di(i);
           } else {
              cout<<"bad eigenvalue i="<<i<<" di="<<di(i)<<"\n";
              exit(0);
           }
        }

        TMatrixD O(ev.GetEigenVectors());
        TMatrixD DOT(d,TMatrixD::kMultTranspose,O);
        g_fcnMatrix=new TMatrixD(O,TMatrixD::kMult,DOT);
        TMatrixD test(*g_fcnMatrix,TMatrixD::kMult,v);
        int error=0;

        for(int i=0;i<n;i++) {
           if(TMath::Abs(test(i,i)-1.0)>1.E-7) {
              error++;
           }
           for(int j=0;j<n;j++) {
              if(i==j) continue;
              if(TMath::Abs(test(i,j)>1.E-7)) error++;
           }
        }

	// calculate chi2
        double chi2 = 0.;
        double ndf = 0.;
        for(int i=0;i<hpt_temp_data->GetNbinsX();i++) {
           double di=hpt_temp_data->GetBinContent(i+1)-hpt_temp_mc->GetBinContent(i+1);
           if(g_fcnMatrix) {
              for(int j=0;j<hpt_temp_data->GetNbinsX();j++) {
                 double dj=hpt_temp_data->GetBinContent(j+1)-hpt_temp_mc->GetBinContent(j+1);
                 chi2+=di*dj*(*g_fcnMatrix)(i,j);
              }
           } else {
              double pull=di/hpt_temp_data->GetBinError(i+1);
              chi2+=pull*pull;
           }
           ndf+=1.0;
        }


	delete g_fcnHist;
        delete g_fcnMatrix;

        delete hpt_temp_data;
        delete hpt_temp_mc;
        delete hrho;

	return chi2/ndf;

}

void ISRUnfold::setMeanMass(bool doSys, bool altMC){

        //Double_t massBins[6] = {50., 65., 80., 100., 200., 350.};
        int nMassBin = 5;

        const TUnfoldBinningV17* temp_binning = nomPtUnfold->GetOutputBinning("Gen_Pt");
	const TVectorD* temp_tvecd = temp_binning->GetDistributionBinning(1);
	const Double_t* massBins = temp_tvecd->GetMatrixArray();


        TH1* hunfolded_mass =  nomMassUnfold->GetOutput("hunfolded_mass",0,0,"*[UO]",kTRUE);
        TH1 *histMCTruth_mass= nomMassUnfold->GetBias("histMCTruth_mass",0,0,"*[UO]",kTRUE);

        for(int ibin = 0; ibin < nMassBin; ibin++){
                hunfolded_mass->GetXaxis()->  SetRange(hunfolded_mass->GetXaxis()->  FindBin(massBins[ibin]+0.01),hunfolded_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
                histMCTruth_mass->GetXaxis()->SetRange(histMCTruth_mass->GetXaxis()->FindBin(massBins[ibin]+0.01),histMCTruth_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));

                meanMass_data.   push_back(hunfolded_mass->GetMean());
                meanMassStatErr_data.push_back(hunfolded_mass->GetMeanError());

                meanMass_mc.   push_back(histMCTruth_mass->GetMean());
                meanMassErr_mc.push_back(histMCTruth_mass->GetMeanError());

                if(altMC){
                    TH1 *histMCTruth_massAlt= sysMassUnfold["Alt"].at(0)->GetBias("histMCTruth_massAlt",0,0,"*[UO]",kTRUE);
                    histMCTruth_massAlt->GetXaxis()->SetRange(histMCTruth_mass->GetXaxis()->FindBin(massBins[ibin]+0.01),histMCTruth_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
                    meanMass_mcAlt.   push_back(histMCTruth_massAlt->GetMean());
                    meanMassErr_mcAlt.push_back(histMCTruth_massAlt->GetMeanError());

	            delete histMCTruth_massAlt;
                }

                if(doSys){
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

        }

	delete hunfolded_mass;
	delete histMCTruth_mass;

        // calculate systematic uncertaintiess
        if(doSys){
            for(int i = 0; i < nMassBin; i++){

               std::map<TString, Double_t> temp_map_;
               std::map<TString, std::vector<Double_t>>::iterator it = meanMass_sysdata.at(i).begin();
               while(it != meanMass_sysdata.at(i).end()){
                    int size_ = it->second.size();
                    if((it->first).CompareTo("FSRDR") == 0 || (it->first).CompareTo("Closure") == 0){ it++; continue;}

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
        }

        for(int i = 0; i < nMassBin; i++){
           Double_t totalSys = 0.;

           if(doSys){
           std::map<TString, Double_t>::iterator it = meanMassErr_sysdata.at(i).begin();
            while(it != meanMassErr_sysdata.at(i).end()){
                 //cout << i << " th mass bin, mass" << it->first << " " << it->second << endl;
                 totalSys += pow(it->second, 2);
                 it++;
            }

            cout << i << " th mass bin, total mass systematic uncertainty: " << sqrt(totalSys) << " statistical error: " << meanMassStatErr_data.at(i) << endl;
           }
           meanMassSysErr_data.push_back(sqrt(totalSys));
	   meanMassTotErr_data.push_back(sqrt(totalSys + pow(meanMassStatErr_data.at(i),2)));
        }
}


// set mean pt from mass and DY mc
void ISRUnfold::setMeanPt(bool doSys, bool altMC){

        int nMassBin = 5;

        // save mean pt for each systematic variation
        for(int i = 0; i < nMassBin; i++){
           TString ibinMass;
           ibinMass.Form("%d", i);

           TH1* hpt_temp_data;
           TH1* hpt_temp_mc;
           TH1* hpt_temp_mcAlt;

           hpt_temp_data = nomPtUnfold->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
           hpt_temp_mc   = nomPtUnfold->GetBias("histMCTruth_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);

           meanPt_data.   push_back(hpt_temp_data->GetMean());
           meanPtStatErr_data.push_back(hpt_temp_data->GetMeanError());

           meanPt_mc.   push_back(hpt_temp_mc->GetMean());
           meanPtErr_mc.push_back(hpt_temp_mc->GetMeanError());

           if(altMC){
            hpt_temp_mcAlt   = sysPtUnfold["Alt"].at(0)->GetBias("histMCTruth_pt_tempAlt",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
            meanPt_mcAlt.   push_back(hpt_temp_mcAlt->GetMean());
            meanPtErr_mcAlt.push_back(hpt_temp_mcAlt->GetMeanError());
           }

           if(doSys){
            std::map<TString, std::vector<TUnfoldDensityV17*>>::iterator it = sysPtUnfold.begin();
	    std::map<TString, std::vector<Double_t>> temp_map; // temp map to save systematic results for a mass bin
                                                               // format: systematic name, mean pt
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
        }

	// calculate systematic uncertainty for each systematic variation
	if(doSys){
	    for(int i = 0; i < nMassBin; i++){

               std::map<TString, Double_t> temp_map_;
               std::map<TString, std::vector<Double_t>>::iterator it = meanPt_sysdata.at(i).begin();
               while(it != meanPt_sysdata.at(i).end()){
	    	int size_ = it->second.size(); // size of systematic variations
	    	if((it->first).CompareTo("FSRDR") == 0 || (it->first).CompareTo("Closure") == 0){ it++; continue;}
	    	
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
        }

	for(int i = 0; i < nMassBin; i++){
	   Double_t totalSys = 0.;
           if(doSys){
           std::map<TString, Double_t>::iterator it = meanPtErr_sysdata.at(i).begin();
            while(it != meanPtErr_sysdata.at(i).end()){
	         //cout << i << " th mass bin, pt" << it->first << " " << it->second << endl;
	         totalSys += pow(it->second, 2);	
	         it++;
	    }
	    
            cout << i << " th mass bin, total pt systematic uncertainty: " << sqrt(totalSys) << " statistical error: " << meanPtStatErr_data.at(i) << endl;
            it = meanPtErr_sysdata.at(i).begin();
            while(it != meanPtErr_sysdata.at(i).end()){
	         cout << i << " th mass bin, pt" << it->first << " " << it->second/ meanPt_data.at(i) * 100.<< endl;
	         totalSys += pow(it->second, 2);	
	         it++;
	    }
           }
	   meanPtSysErr_data.push_back(sqrt(totalSys));
	   meanPtTotErr_data.push_back(sqrt(totalSys + pow(meanPtStatErr_data.at(i),2)));
	}// loop for mass bins
}

void ISRUnfold::drawISRresult(TString outpdf, bool altMC){

        gROOT->SetBatch();

	Double_t zero_err[5] = {0.};

	// temporary map to show variations respect to dR to add FSR photons
	//std::map<int, std::vector<Double_t>> FSR_temp_pt;
	//std::map<int, std::vector<Double_t>> FSR_temp_mass;

	//for(int i = 0; i < 5; i++){

	//	int size = meanPt_sysdata.at(i)["FSRDR"].size();
	//	for(int j = 0; j < size; j++){
	//		FSR_temp_pt[j].  push_back((meanPt_sysdata.at(i)["FSRDR"]).at(j));
	//		FSR_temp_mass[j].push_back((meanMass_sysdata.at(i)["FSRDR"]).at(j));
	//	}
	//}

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

        //TGraphErrors *grFSRtest;
        //int map_size = FSR_temp_mass.size();
        //for(int i = 0; i < map_size; i++){
        //        grFSRtest = new TGraphErrors(5, &(FSR_temp_mass[i])[0], &(FSR_temp_pt[i])[0], zero_err, zero_err);
        //        grFSRtest->SetLineColor(kBlue);
        //        grFSRtest->SetMarkerColor(kBlue);
        //        grFSRtest->SetMarkerStyle(24);
        //        grFSRtest->SetMarkerSize(0.5);
        //        grFSRtest->SetLineStyle(1);
        //        grFSRtest->Draw("ple same");
        //}

        TGraphErrors *grMC = new TGraphErrors(5, &meanMass_mc[0], &meanPt_mc[0], &meanMassErr_mc[0], &meanPtErr_mc[0]);
        grMC->SetLineColor(kRed);
        grMC->SetMarkerColor(kRed);
        grMC->SetMarkerStyle(20);
        grMC->SetMarkerSize(1.);
        grMC->SetLineStyle(1);
        grMC->SetLineColor(kRed);
        grMC->Draw("lpe same");

        TGraphErrors *grMCAlt;
        if(altMC){
            grMCAlt = new TGraphErrors(5, &meanMass_mcAlt[0], &meanPt_mcAlt[0], &meanMassErr_mcAlt[0], &meanPtErr_mcAlt[0]);
            grMCAlt->SetLineColor(kRed);
            grMCAlt->SetMarkerColor(kRed);
            grMCAlt->SetMarkerStyle(24);
            grMCAlt->SetMarkerSize(1.);
            grMCAlt->SetLineStyle(2);
            grMCAlt->SetLineColor(kRed);
            grMCAlt->Draw("lpe same");
        }

        grUnfolded->Draw("pe same");

        TLegend* leg_ = new TLegend(0.22, 0.6, 0.55, 0.75,"","brNDC");
        leg_->SetTextSize(0.03);
        //leg_->SetFillStyle(1);
        leg_->SetBorderSize(0);
        leg_->AddEntry(grUnfolded, "Unfolded electron data", "pe");
        leg_->AddEntry(grMC, "DY MC at pre FSR (aMC@NLO)", "pe");
        if(altMC) leg_->AddEntry(grMCAlt, "DY MC at pre FSR (Madgraph)", "pe");
        leg_->Draw();

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
	delete c1;
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

// https://root.cern.ch/root/html/tutorials/graphs/graphtext.C.html
void ISRUnfold::drawtext(TGraph *g)
{
   Int_t i,n;
   Double_t x,y;
   TLatex *l;

   n = g->GetN();
   for (i=0; i<n; i++) {
       g->GetPoint(i,x,y);
       l = new TLatex(x-0.1,y+0.02,Form("%d",i+1));
       l->SetTextSize(0.015);
       l->SetTextFont(42);
       l->SetTextAlign(21);
       l->Draw();
   }
}

void ISRUnfold::drawSysPlots(TString outpdf, int nthMassBin, TString sysName){

        gROOT->SetBatch();

        TString ibinMass;
        ibinMass.Form("%d", nthMassBin);

        c1 = new TCanvas("c1","c1", 50, 50, 850, 700);
        c1->cd();
        gStyle->SetOptFit(0);

        c1->SetBottomMargin(0.2);
        c1->SetTopMargin(0.08);
        c1->SetTicks(1);
        //c1->SetGridx();
        //c1->SetGridy();

	int sysSize = sysPtUnfold[sysName].size();
        double amaxMass=meanMass_data[nthMassBin]+ 2. * meanMassTotErr_data[nthMassBin];
        double aminMass=meanMass_data[nthMassBin]- 2. * meanMassTotErr_data[nthMassBin];
        double amaxPt=meanPt_data[nthMassBin]+ 2. * meanPtTotErr_data[nthMassBin];
        double aminPt=meanPt_data[nthMassBin]- 2. * meanPtTotErr_data[nthMassBin];
        for(int i=0;i<sysSize;i++) {
           amaxMass=TMath::Max(amaxMass, (meanMass_sysdata.at(nthMassBin)[sysName])[i]);
           aminMass=TMath::Min(aminMass, (meanMass_sysdata.at(nthMassBin)[sysName])[i]);

           amaxPt=TMath::Max(amaxPt, (meanPt_sysdata.at(nthMassBin)[sysName])[i]);
           aminPt=TMath::Min(aminPt, (meanPt_sysdata.at(nthMassBin)[sysName])[i]);
        }


        TGraphErrors *grUnfolded = new TGraphErrors(1, &meanMass_data[nthMassBin], &meanPt_data[nthMassBin], &meanMassTotErr_data[nthMassBin], &meanPtTotErr_data[nthMassBin]);
        grUnfolded->SetLineColor(kBlack);
        grUnfolded->SetMarkerColor(kBlack);
        grUnfolded->SetMarkerStyle(20);
        grUnfolded->SetMarkerSize(1.);
        grUnfolded->SetLineStyle(1);
        grUnfolded->Draw("ape");
        grUnfolded->GetYaxis()->SetRangeUser(aminPt*0.995,amaxPt*1.005);
        grUnfolded->GetXaxis()->SetLimits(aminMass*0.995,amaxMass*1.005);
        grUnfolded->GetYaxis()->SetTitle("Average p_{T} (GeV)");
        grUnfolded->GetXaxis()->SetTitle("Average Mass (GeV)");

	// meanPtStatErr_data
        TGraphErrors *grStatErr = new TGraphErrors(1, &meanMass_data[nthMassBin], &meanPt_data[nthMassBin], &meanMassStatErr_data[nthMassBin], &meanPtStatErr_data[nthMassBin]);
        grStatErr->SetFillColor(kBlue);;
	grStatErr->SetFillColorAlpha(kBlue,0.3);
        grStatErr->Draw("sameE2");

        TGraph *grSys = new TGraph(sysSize, &(meanMass_sysdata.at(nthMassBin)[sysName])[0], &(meanPt_sysdata.at(nthMassBin)[sysName])[0] );
        grSys->SetLineColor(kRed);
        grSys->SetMarkerColor(kRed);
        grSys->SetMarkerStyle(20);
        grSys->SetMarkerSize(0.8);
        grSys->SetLineStyle(1);
        grSys->Draw("pe same ");
	drawtext(grSys);


        grUnfolded->Draw("pe same");

        TLegend* leg = new TLegend(0.6, 0.70, 0.85, 0.9,"","brNDC");
        leg->SetTextSize(0.04);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);

        leg->AddEntry(grUnfolded, "Nominal point", "pl");
        leg->AddEntry(grSys, "Systematic (" + sysName + ")" , "pl");
        leg->Draw();


        CMS_lumi( c1, 4, 11 );
        c1->SaveAs(outpdf+"_"+ibinMass+"_"+sysName+".pdf");
        delete grUnfolded;
	delete c1;

}

void ISRUnfold::drawNominalRecoPlots(TString outpdf, TString filepath, TString var, int nthMassBin, TString sysName){

        TString ibinMass;
        ibinMass.Form("%d", nthMassBin);

        gROOT->SetBatch();

        setTDRStyle();
        writeExtraText = true;
        extraText  = "work in progress";

	// get binning definition
        const TUnfoldBinning* temp_binning = nomPtUnfold->GetInputBinning("Rec_Pt");

	TFile* filein = new TFile(filepath);
        TH1* hdysig = (TH1*)filein->Get("hPtRecnominal"); // get DY 
        TH1* hdysigNoUO = temp_binning->ExtractHistogram("hdysig", hdysig, 0, kTRUE, "pt[UO];mass[UOC"+ibinMass+"]");

        TH1* hpt_temp_data;
        TH1F *ratio;

	// bkg subracted data
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
        hdysigNoUO->Draw("hist same");
        hdysigNoUO->SetLineColor(kRed);


        TString mean_nom;
        mean_nom.Form("%.5f", hpt_temp_data->GetMean());

        TLegend* leg_nom = new TLegend(0.45, 0.70, 0.75, 0.9,"","brNDC");
        leg_nom->SetNColumns(2);
        leg_nom->SetTextSize(0.055);
        leg_nom->SetFillStyle(0);
        leg_nom->SetBorderSize(0);

        leg_nom->AddEntry(hpt_temp_data, "Bkg subtracted data (mean: " + mean_nom + ")", "pl");
        leg_nom->Draw();

        TLatex chi2_norm;
        TString chi2_;

        // TODO add Mass option 
        if(var.CompareTo("Pt") == 0 ){
                chi2_.Form("%f", Chi2Test(hpt_temp_data, hdysigNoUO));
                chi2_norm.DrawLatexNDC(0.2, 0.3, "#chi^{2}/NDOF= " + chi2_);
        }

        c1->cd();

        TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.4);
        pad2->SetTopMargin(0.05);
        pad2->SetBottomMargin(0.2);
        pad2->SetTicks(1);
        pad2->SetGridy(1);
        pad2->Draw();

        pad2->cd();

	ratio = ((TH1F*)hpt_temp_data->Clone("pt_temp"));
	ratio->Divide(hdysigNoUO);
        ratio->Draw("hist");
        ratio->GetYaxis()->SetTitle("Data (bkg sub)/ DY MC");
        ratio->GetXaxis()->SetTitle("p_{T} (GeV)");
        ratio->SetMinimum(0.8);
        ratio->SetMaximum(1.2);
        ratio->SetTitle("");
        ratio->GetXaxis()->SetTitleOffset(1.5);
        ratio->GetYaxis()->SetNdivisions(515);
        ratio->SetLineColor(kBlack);
        //ratio->SetLineStyle(2);

        CMS_lumi( c1, 4, 0 );
        c1->cd();
        c1->SaveAs(outpdf+"_input_"+ibinMass+".pdf");

        delete hpt_temp_data;
	delete hdysigNoUO;
        delete ratio;
        delete pad1;
        delete pad2;
        delete c1;

}

void ISRUnfold::studyFSRDRPlots(TString outpdf, TString var, int nthMassBin){

        gROOT->SetBatch();

        setTDRStyle();
        writeExtraText = true;
        extraText  = "work in progress";

        TString ibinMass;
        ibinMass.Form("%d", nthMassBin);

        TH1* hunfolded_data;
        TH1* hpreFSR_dressed0p1;
        TH1F *ratio;

        // get nominal unfoled result
        if(var.CompareTo("Pt") == 0 ){
                hunfolded_data  = nomPtUnfold->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
                hpreFSR_dressed0p1   =    sysPtUnfold["FSRDR"].at(0)->GetBias("histMCTruth_pt_tempAlt",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
        }
        if(var.CompareTo("Mass") == 0 ){
                hunfolded_data  = nomMassUnfold->GetOutput("hunfolded_mass_temp",0,0,"*[UO]",kTRUE);
                hpreFSR_dressed0p1   =    sysMassUnfold["FSRDR"].at(0)->GetOutput("hunfolded_mass_systemp",0,0,"*[UO]",kTRUE);
        }

        ratio= ((TH1F*)hunfolded_data->Clone("ratio"));
        ratio->Divide(hpreFSR_dressed0p1);

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
        hpreFSR_dressed0p1->Draw("histsame");
        hunfolded_data->SetMarkerStyle(20);
        hunfolded_data->SetMarkerSize(.7);
        hunfolded_data->SetLineColor(kBlack);
        hpreFSR_dressed0p1->SetLineColor(kRed);
        hunfolded_data->GetYaxis()->SetTitle("Events/bin");
        hunfolded_data->SetMinimum(10.);

        c1->cd();

        TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.4);
        pad2->SetTopMargin(0.05);
        pad2->SetBottomMargin(0.2);
        pad2->SetTicks(1);
        pad2->SetGridy(1);
        pad2->Draw();

        pad2->cd();
        ratio->Draw("hist");
        ratio->GetYaxis()->SetTitle("#DeltaR=X/ #DeltaR=0.1");
        ratio->GetXaxis()->SetTitle("p_{T} at pre FSR(GeV)");
        ratio->SetMinimum(0.5);
        ratio->SetMaximum(1.5);
        ratio->SetTitle("");
        ratio->GetXaxis()->SetTitleOffset(1.5);
        ratio->GetYaxis()->SetNdivisions(515);
        ratio->SetLineColor(kRed);

        TH1* hratio_dr_temp;
        int sysSize = sysPtUnfold["FSRDR"].size();
        for(int i = 1; i < sysSize; i++){

		TH1* htemp;
        	if(var.CompareTo("Pt") == 0 ){
        	        htemp   =    sysPtUnfold["FSRDR"].at(i)->GetBias("histMCTruth_pt_tempAlt_",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
        	}
        	if(var.CompareTo("Mass") == 0 ){
        	        htemp   =    sysMassUnfold["FSRDR"].at(i)->GetOutput("hunfolded_mass_systemp_",0,0,"*[UO]",kTRUE);
        	}
		hratio_dr_temp = ((TH1F*)htemp->Clone("htemp"));
		hratio_dr_temp->Divide(hpreFSR_dressed0p1);
                hratio_dr_temp->Draw("histsame");
                hratio_dr_temp->SetLineColor(i);
                hratio_dr_temp->SetLineStyle(2);
                delete htemp;
	}

        CMS_lumi( c1, 4, 0 );
        c1->cd();
        c1->SaveAs(outpdf+"_"+ibinMass+"_"+var+".pdf");

	//delete hratio_dr_temp;
        delete hunfolded_data;
        delete hpreFSR_dressed0p1;
        delete pad1;
        delete pad2;
        delete c1;
}


void ISRUnfold::drawNominalPlots(TString outpdf, TString var, int nthMassBin, TString sysName){

        const TUnfoldBinningV17* temp_binning = nomPtUnfold->GetOutputBinning("Gen_Pt");
        const TUnfoldBinningV17* temp_binning_mass = nomMassUnfold->GetOutputBinning("Gen_Mass");
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
        TH1F *ratio_MG_aMCNLO;
        TH1F *ratio_sys_err;


        TH1F *ratio_closure;
	// get nominal unfoled result
	if(var.CompareTo("Pt") == 0 ){
        	if(sysName.CompareTo("Closure")!=0)hunfolded_data  = nomPtUnfold->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
		else hunfolded_data  = sysPtUnfold["Closure"].at(0)->GetOutput("hunfolded_pt_closure",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
		hunfolded_sys_err= ((TH1F*)hunfolded_data->Clone("sysErr")); 
		hpreFSR_mc   = nomPtUnfold->GetBias("histMCTruth_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
	}
	if(var.CompareTo("Mass") == 0 ){
        	if(sysName.CompareTo("Closure")!=0) hunfolded_data  = nomMassUnfold->GetOutput("hunfolded_mass_temp",0,0,"*[UO]",kTRUE);
		else hunfolded_data  = sysMassUnfold["Closure"].at(0)->GetOutput("hunfolded_mass_closure",0,0,"*[UO]",kTRUE);
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
       
                if(var.CompareTo("Pt") == 0 && (sysName.CompareTo("Closure") != 0 ))   hdatasys_temp = sysPtUnfold[sysName].at(i)->GetOutput("hunfolded_pt_systemp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);

                if(var.CompareTo("Pt") == 0 && (sysName.CompareTo("Closure") == 0 )){
        		TFile* filein = new TFile("/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/output/2016/electron/DYtoEE.root");

        		//TH1* hdatasys_temp_ = (TH1*)filein->Get("h" + var + "Gennominal");
			//hdatasys_temp = temp_binning->ExtractHistogram("hunfolded_pt_systemp",hdatasys_temp_, 0, kTRUE,"pt[UO];mass[UOC"+ibinMass+"]");

        		TString ibinMass_;
        		ibinMass_.Form("%d", nthMassBin+1);
        		hdatasys_temp = (TH1*)filein->Get("h" + var + "_m"+ibinMass_); //hPt_m1
	
                        ratio_closure= ((TH1F*)hdatasys_temp->Clone("ratio_closure"));
			
        		ratio_closure->Divide(hpreFSR_mc);
		}

		if(var.CompareTo("Pt") == 0 && (sysName.CompareTo("Alt") == 0 || sysName.CompareTo("FSRDR") == 0 ) && i == 0){
	 		hmcsys_temp = sysPtUnfold[sysName].at(i)->GetBias("histMCTruth_pt_tempAlt",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE); // get alternative DY MC 
			hmcsys_temp->SetDirectory(0);
                        ratio_MG_aMCNLO= ((TH1F*)hmcsys_temp->Clone("ratio_MG_aMCNLO"));
                        ratio_MG_aMCNLO->Divide(hpreFSR_mc);
		}

                if(var.CompareTo("Mass") == 0 && (sysName.CompareTo("Closure") != 0 )){
              		hdatasys_temp = sysMassUnfold[sysName].at(i)->GetOutput("hunfolded_mass_systemp",0,0,"*[UO]",kTRUE);
             		hdatasys_temp->GetXaxis()->SetRange(hdatasys_temp->GetXaxis()->FindBin(massBins[nthMassBin]+0.01), hdatasys_temp->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));
                }
		if(var.CompareTo("Mass") == 0 && (sysName.CompareTo("Closure") == 0 )){
                        TFile* filein = new TFile("/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/output/2016/electron/DYtoEE.root");

                        //TH1* hdatasys_temp_ = (TH1*)filein->Get("h" + var + "Gennominal");
                        //hdatasys_temp = temp_binning_mass->ExtractHistogram("hunfolded_mass_systemp",hdatasys_temp_, 0, kTRUE,"*[UO]");
             		//hdatasys_temp->GetXaxis()->SetRange(hdatasys_temp->GetXaxis()->FindBin(massBins[nthMassBin]+0.01), hdatasys_temp->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));

                        hdatasys_temp = (TH1*)filein->Get("h" + var);
             		hdatasys_temp->GetXaxis()->SetRange(hdatasys_temp->GetXaxis()->FindBin(massBins[nthMassBin]+0.01), hdatasys_temp->GetXaxis()->FindBin(massBins[nthMassBin+1]-0.01));

                        ratio_closure= ((TH1F*)hdatasys_temp->Clone("ratio_closure"));
        		ratio_closure->Divide(hpreFSR_mc);
			

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
		   if(ratio_err < 5.e-6) ratio_err = 1.e-6;
                   ratio_sys_err->SetBinError(ibin, ratio_err);

                   delete ratio_temp;
              }// loop for bin contents

                hpt_sys_temp = ((TH1F*)hdatasys_temp->Clone("pt_temp"));
                hpt_sys_temp->Draw("histsame");
                hpt_sys_temp->SetLineColor(kBlack);
                hpt_sys_temp->SetLineStyle(2);

		if( (sysName.CompareTo("Closure") == 0 ) ){

        		TString mean_mc_;
        		mean_mc_.Form("%.5f", hdatasys_temp->GetMean());

        		TLegend* leg_mc_ = new TLegend(0.45, 0.30, 0.75, 0.5,"","brNDC");
        		leg_mc_->SetNColumns(2);
        		leg_mc_->SetTextSize(0.055);
        		leg_mc_->SetFillStyle(0);
        		leg_mc_->SetBorderSize(0);

        		leg_mc_->AddEntry(hpt_sys_temp, "Truth MC (mean: " + mean_mc_ + ")", "pl");
        		leg_mc_->Draw();
		}

		if(hmcsys_temp!=NULL){
			cout << "draw " << sysName << " MC histogram " << " 1 bin content: " << hmcsys_temp->GetBinContent(1) << endl;
			hmcsys_temp->Draw("histsame");
			hmcsys_temp->SetLineColor(kBlue);
		}


              delete hdatasys_temp;
	      //delete hmcsys_temp;
          }

	if(sysName.CompareTo("Closure")!=0){
	hunfolded_sys_err->Draw("E2same");
        hunfolded_sys_err->SetMarkerSize(0);
        hunfolded_sys_err->SetFillColorAlpha(kBlack,0.3);
	}

        TString mean_nom;
        mean_nom.Form("%.5f", hunfolded_data->GetMean());

        TLegend* leg_nom = new TLegend(0.45, 0.4, 0.75, 0.6,"","brNDC");
        //leg_nom->SetNColumns(2);
        leg_nom->SetTextSize(0.055);
        leg_nom->SetFillStyle(0);
        leg_nom->SetBorderSize(0);

        if(sysName.CompareTo("Closure")!=0){
            leg_nom->AddEntry(hunfolded_data, "Unfolded data (mean: " + mean_nom + ")", "pl");
            leg_nom->AddEntry(hpreFSR_mc, "aMC@NLO", "l");
	    //if(sysName.CompareTo("Alt") == 0){
            //    leg_nom->AddEntry(hmcsys_temp, "MG (LO)", "l");
            //}
        }
	else leg_nom->AddEntry(hunfolded_data, "Unfolded DY MC (mean: " + mean_nom + ")", "pl");
        leg_nom->Draw();

	TLatex chi2_norm;
        TString chi2_;

	// TODO add Mass option	
	if(var.CompareTo("Pt") == 0 ){
		chi2_.Form("%f",DoFit("Pt", nthMassBin));
 		chi2_norm.DrawLatexNDC(0.2, 0.3, "#chi^{2}/NDOF= " + chi2_);
	}

	c1->cd();

        TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.4);
        pad2->SetTopMargin(0.05);
        pad2->SetBottomMargin(0.2);
        pad2->SetTicks(1);
        pad2->SetGridy(1);
        pad2->Draw();

        pad2->cd();
        ratio->Draw("histe");
        ratio->GetYaxis()->SetTitle("Data/ MC");
        ratio->GetXaxis()->SetTitle("p_{T} at pre FSR(GeV)");
        ratio->SetMinimum(0.8);
        ratio->SetMaximum(1.2);
        ratio->SetTitle("");
        ratio->GetXaxis()->SetTitleOffset(1.5);
        ratio->GetYaxis()->SetNdivisions(515);

        TLegend* leg_ratios = new TLegend(0.45, 0.7, 0.75, 0.9,"","brNDC");
        //leg_ratios->SetNColumns(2);
        leg_ratios->SetTextSize(0.055);
        leg_ratios->SetFillStyle(0);
        leg_ratios->SetBorderSize(0);
        leg_ratios->AddEntry(ratio, "unfolded data/ aMC@NLO" , "pl");

        if(sysName.CompareTo("Alt") == 0){
            ratio_MG_aMCNLO->Draw("histsamee");
            ratio_MG_aMCNLO->SetLineStyle(2);
            ratio_MG_aMCNLO->SetLineColor(kRed);
            leg_ratios->AddEntry(ratio_MG_aMCNLO, "MG/ aMC@NLO", "l");
        }

	if(sysName.CompareTo("Closure")!=0){
	    ratio_sys_err->Draw("E2same");
            ratio_sys_err->SetMarkerSize(0);
            ratio_sys_err->SetFillColorAlpha(kBlack,0.3);
            leg_ratios->AddEntry(ratio_sys_err, "systematic " + sysName , "F");

	}
	else{
		ratio_closure->Draw("histsamee");
		ratio_closure->SetLineColor(kRed);

	}
        leg_ratios->Draw();



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
