#include "histTUnfold.h"

void histTUnfold::CreateHistMap(const int which_unfold, TString hname){

        if(which_unfold == pt_unfold){
                 histMaps.insert(std::pair<TString, TH1*>("pt_"+hname, ptBinningRec->CreateHistogram("hPtRec"+hname)));
        }

        if(which_unfold == mass_unfold){
                 histMaps.insert(std::pair<TString, TH1*>("mass_"+hname, massBinningRec->CreateHistogram("hMassRec"+hname)));
        }
}

void histTUnfold::CreateHist2DMap(const int which_unfold, TString hname){

	if(which_unfold == pt_unfold){
		 hist2DMaps.insert(std::pair<TString, TH2*>("pt_"+hname, TUnfoldBinningV17::CreateHistogramOfMigrations(ptBinningGen, ptBinningRec,"hmcPtGenRec" + hname)));
	}

	if(which_unfold == mass_unfold){
		 hist2DMaps.insert(std::pair<TString, TH2*>("mass_"+hname, TUnfoldBinningV17::CreateHistogramOfMigrations(massBinningGen, massBinningRec,"hmcMassGenRec" + hname)));
	}
}

void histTUnfold::FillHistogram(const int which_unfold, TString hname, Double_t recoPt, Double_t recoMass, Double_t wreco){

	if(which_unfold == pt_unfold){
		(histMaps.find(hname)->second)->Fill(ptBinningRec->GetGlobalBinNumber(recoPt, recoMass),  wreco);
	}

	if(which_unfold == mass_unfold){
		if(recoPt < 100.) (histMaps.find(hname)->second)->Fill(massBinningRec->GetGlobalBinNumber(recoMass),  wreco);
	}
}

void histTUnfold::FillMigration2DM(const int which_unfold, bool selected, TString hname, Double_t recoPt, Double_t RecoMass, Double_t truthPt, Double_t truthMass, Double_t wreco, Double_t wgen){

	int binZero=0;

	if(selected){
		if( which_unfold == pt_unfold){
			(hist2DMaps.find(hname)->second)->Fill(ptBinningGen->GetGlobalBinNumber(truthPt, truthMass), ptBinningRec->GetGlobalBinNumber(recoPt, RecoMass), wgen*wreco);
			(hist2DMaps.find(hname)->second)->Fill(ptBinningGen->GetGlobalBinNumber(truthPt, truthMass), binZero, wgen*(1.-wreco));
		}
                if( which_unfold == mass_unfold){
                        if(recoPt < 100.) (hist2DMaps.find(hname)->second)->Fill(massBinningGen->GetGlobalBinNumber(truthMass), massBinningRec->GetGlobalBinNumber(RecoMass), wgen*wreco);
                        if(recoPt < 100.) (hist2DMaps.find(hname)->second)->Fill(massBinningGen->GetGlobalBinNumber(truthMass), binZero, wgen*(1.-wreco));
                }
	}
	else{
		if( which_unfold == pt_unfold) (hist2DMaps.find(hname)->second)->Fill(ptBinningGen->GetGlobalBinNumber(truthPt, truthMass), binZero, wgen);
		if( which_unfold == mass_unfold) (hist2DMaps.find(hname)->second)->Fill(massBinningGen->GetGlobalBinNumber(truthMass), binZero, wgen);
	}
}

void histTUnfold::SetPtBinningRec(){

 // FIXME save binning information in other place and read from there
 const int nmassbin_fine=5;
 double massbin_fine[nmassbin_fine+1]={50,65,80,100,200,350};

 // pt bins for reco
 //const int nptbin_fine=50;
 //double ptbin_fine[nptbin_fine+1];
 //for(int i = 0; i < nptbin_fine + 1; i++){
 //   ptbin_fine[i] = i*2;
 //}

 const int nptbin_fine=17;
 double ptbin_fine[nptbin_fine+1]={0., 2., 4., 6., 8., 10., 12., 14., 18., 22., 28., 35., 45., 55., 65., 75., 85., 100.};

 ptBinningRec= new TUnfoldBinningV17("Rec_Pt");
 ptBinningRec->AddAxis("pt",nptbin_fine,ptbin_fine,false,true);
 ptBinningRec->AddAxis("mass",nmassbin_fine,massbin_fine,true,true);
}

void histTUnfold::SetPtBinningGen(){

 // FIXME save binning information in other place and read from there
 const int nmassbin_wide=5;
 double massbin_wide[nmassbin_wide+1]={50,65,80,100,200,350};

 // pt bins for gen
 //const int nptbin_wide=20;
 //double ptbin_wide[nptbin_wide+1];
 //for(int i = 0; i < nptbin_wide + 1; i++){
 //    ptbin_wide[i] = i*5;
 //}

 //const int nptbin_wide=9;
 //double ptbin_wide[nptbin_wide+1]={0., 4., 8., 12., 18., 28., 40., 55., 80., 100.};

 const int nptbin_wide=17;
 double ptbin_wide[nptbin_wide+1]={0., 2., 4., 6., 8., 10., 12., 14., 18., 22., 28., 35., 45., 55., 65., 75., 85., 100.};

 ptBinningGen=(new TUnfoldBinningV17("Gen_Pt"));
 ptBinningGen->AddAxis("pt",nptbin_wide,ptbin_wide,false,true);
 ptBinningGen->AddAxis("mass",nmassbin_wide,massbin_wide,true,true);

}

void histTUnfold::SetMassBinningRec(){

 // FIXME save binning information in other place and read from there
 //const int nbin_fine=58;
 const int nbin_fine=54;
 //double bin_fine[nbin_fine+1]={40,42.5,45,47.5,50,52.5,55,57.5,60,62.5,65,67.5,70,72.5,75,77.5,80,82.5,85,87.5,90,92.5,95,97.5,100,102.5,105,107.5,110,112.5,115,117.5,120,123,126,129.5,133,137,141,145.5,150,155,160,165.5,171,178,185,192.5,200,209,218,229,240,254,268,284,300,325,350};
 double bin_fine[nbin_fine+1]={50,52.5,55,57.5,60,62.5,65,67.5,70,72.5,75,77.5,80,82.5,85,87.5,90,92.5,95,97.5,100,102.5,105,107.5,110,112.5,115,117.5,120,123,126,129.5,133,137,141,145.5,150,155,160,165.5,171,178,185,192.5,200,209,218,229,240,254,268,284,300,325,350};

 massBinningRec=(new TUnfoldBinningV17("Rec_Mass"));
 massBinningRec->AddAxis("reco mass",nbin_fine,bin_fine,false,false);
}

void histTUnfold::SetMassBinningGen(){

 // FIXME save binning information in other place and read from there
 //const int nbin_wide=29;
 const int nbin_wide=27;
 //double bin_wide[nbin_wide+1]={40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,126,133,141,150,160,171,185,200,218,240,268,300,350};
 double bin_wide[nbin_wide+1]={50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,126,133,141,150,160,171,185,200,218,240,268,300,350};

 massBinningGen=(new TUnfoldBinningV17("Gen_Mass"));
 massBinningGen->AddAxis("gen mass",nbin_wide,bin_wide,true,true);
}
/*
void saveRecoHists(TFile *filein, TFile *fileout1, histTUnfold &recoHist, TString channel){ // TODO add list of systematics

	gROOT->SetBatch();
	TH1::SetDefaultSumw2();
	
	TTree *treco=(TTree *)filein->Get("tree");
	treco->SetBranchAddress("ptRec",&ptRec);
	treco->SetBranchAddress("mRec",&mRec);
	treco->SetBranchAddress("IsElEl",&IsElEl);
	treco->SetBranchAddress("IsMuMu",&IsMuMu);
	treco->SetBranchAddress("ispassRec",&ispassRec);
	treco->SetBranchAddress("isBveto",&isBveto);
	treco->SetBranchAddress("weightGen",&weightGen);
	treco->SetBranchAddress("weightRec",&weightRec);
	treco->SetBranchAddress("DYtautau",&DYtautau);
	treco->SetBranchAddress("bTagReweight",&bTagReweight); // FIXME may it is better to change name bTagReweight to bTagSF
	nentries=treco->GetEntries();
	
	// TODO based on the info in histTUnfold make map for systematics
	//for(int i=0;i<nentries;i++){
	for(int i=0;i<10000;i++){
	  if(i%10000000==0) cout<<i<<endl;
	  treco->GetEntry(i);

          int isdilep = 0;
          if( channel.CompareTo("electron") == 0 ) isdilep = IsElEl;
          if( channel.CompareTo("muon") == 0 ) isdilep = IsMuMu;

	   if(isdilep && ispassRec && isBveto && ptRec->at(0) > 25 && ptRec->at(1) > 15){
	         fileout1->cd();

		 recoHist.FillHistogram(pt_unfold, "pt_nominal", ptRec->at(2), mRec->at(2), weightGen*weightRec*bTagReweight);
		 recoHist.FillHistogram(mass_unfold, "mass_nominal", ptRec->at(2), mRec->at(2), weightGen*weightRec*bTagReweight);
	   }
	
	}// event loop

	fileout1->cd();
	recoHist.GetPtBinningRec()->Write();
	recoHist.GetMassBinningRec()->Write();
	
	delete treco;
	//fileout->cd();
	
	// seems ptHistograms automatically written 
	// recoHist.ptHists.at(0)->Write();
	// delete outputFile;
}

void saveSigHists(TFile *filein, TFile *fileout1, TFile *fileout2, histTUnfold &sigHist, histTUnfold &recoHist, TString channel, Double_t temp_kfactor){ // TODO add list of systematics

	gROOT->SetBatch();
	TH1::SetDefaultSumw2();
	
	TTree *tsignal=(TTree *)filein->Get("tree");
	tsignal->SetBranchAddress("ptPreFSR",&ptPreFSR);
	tsignal->SetBranchAddress("mPreFSR",&mPreFSR);
	tsignal->SetBranchAddress("ptPostFSR",&ptPostFSR);
	tsignal->SetBranchAddress("mPostFSR",&mPostFSR);
	tsignal->SetBranchAddress("ptRec",&ptRec);
	tsignal->SetBranchAddress("mRec",&mRec);
	tsignal->SetBranchAddress("IsElEl",&IsElEl);
	tsignal->SetBranchAddress("IsMuMu",&IsMuMu);
	tsignal->SetBranchAddress("ispassRec",&ispassRec);
	tsignal->SetBranchAddress("isBveto",&isBveto);
	tsignal->SetBranchAddress("isdimuon",&isdimuon);
	tsignal->SetBranchAddress("isdielectron",&isdielectron);
	tsignal->SetBranchAddress("weightGen",&weightGen);
	tsignal->SetBranchAddress("weightRec",&weightRec);
	tsignal->SetBranchAddress("ZPtCor",&ZPtCor);
	tsignal->SetBranchAddress("bTagReweight",&bTagReweight);
	tsignal->SetBranchAddress("DYtautau",&DYtautau);
	tsignal->SetBranchAddress("isfiducialPreFSR",&isfiducialPreFSR);
	tsignal->SetBranchAddress("isfiducialPostFSR",&isfiducialPostFSR);

	tsignal->SetBranchAddress("AlphaS",&AlphaS);
	tsignal->SetBranchAddress("Scale",&Scale);
	
	// FIXME 
	TFile* fZptWeight = new TFile("/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/etc/ZptWeight/aMCNLO/electron/ZptWeight_electron.root", "r");
	//TFile* fZptWeight = new TFile("/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/etc/ZptWeight/MG/electron/ZptWeight_electron.root", "r");
	TH1 *hZptWeight;
	fZptWeight->GetObject("ZptWeight", hZptWeight);
	
	nentries=tsignal->GetEntries();
	
	// TODO based on the info in histTUnfold make map for systematics
	//for(int i=0;i<nentries;i++){
	for(int i=0;i<10000;i++){
	  if(i%10000000==0) cout<<i<<endl;
	  tsignal->GetEntry(i);
	
	  weightGen *= temp_kfactor;

          int isdilep = 0, issignal = 0;
          if( channel.CompareTo("electron") == 0 ){
                  isdilep = IsElEl;
                  issignal = isdielectron;
          }
          if( channel.CompareTo("muon") == 0 ){
                  isdilep = IsMuMu;
                  issignal = isdimuon;
          }

	  if( std::isinf(bTagReweight) ){
		 bTagReweight = 1.; //FIXME check analyzer later
		 std::cout << "nan bTagReweight... check later" << std::endl;
          }
	
	   if(!sigHist.isInc){
	      if(isdilep && ispassRec){
	            fileout1->cd();

		    sigHist.FillHistogram(pt_unfold, "pt_nominal", ptRec->at(2), mRec->at(2), weightGen*weightRec*bTagReweight);
		    sigHist.FillHistogram(mass_unfold, "mass_nominal", ptRec->at(2), mRec->at(2), weightGen*weightRec*bTagReweight);

	      }
	   }
	   else{ // for inclusive DY MC
	
	        Double_t zptWeight = 1.;
	        zptWeight = ZPtCor;
	
	        if(isdilep && ispassRec && isBveto && ptRec->at(0) > 25 && ptRec->at(1) > 15){
	           if(DYtautau){
	              fileout2->cd();

		      recoHist.FillHistogram(pt_unfold, "pt_nominal", ptRec->at(2), mRec->at(2), weightGen*weightRec*bTagReweight);
		      recoHist.FillHistogram(mass_unfold, "mass_nominal", ptRec->at(2), mRec->at(2), weightGen*weightRec*bTagReweight);

	           }
	           else{
	              fileout1->cd();

		      sigHist.FillHistogram(pt_unfold, "pt_nominal", ptRec->at(2), mRec->at(2), weightGen*weightRec*bTagReweight);
		      sigHist.FillHistogram(mass_unfold, "mass_nominal", ptRec->at(2), mRec->at(2), weightGen*weightRec*bTagReweight);

	           }
	        }

	       /////////////////////////////////////////// fill migration matrix ////////////////////////////////////////
	       //
	       if(issignal){
	           fileout1->cd();

		   Double_t diptrec = -999., dimassrec = -999.;
		   Double_t diptgen = -999., dimassgen = -999.;

	           bool selected = (isdilep && ispassRec && isBveto && ptRec->at(0) > 25 && ptRec->at(1) > 15);	
		   // check if there are selected dilepton at detector to avoid memory error
		   if(ptRec->size() == 3){ diptrec = ptRec->at(2), dimassrec = mRec->at(2);}

		   sigHist.FillMigration2DM( pt_unfold, selected, "pt_nominal",  diptrec, dimassrec, ptPreFSR->at(2), mPreFSR->at(2), weightRec*bTagReweight, weightGen); 
	           sigHist.FillMigration2DM( mass_unfold, selected, "mass_nominal", diptrec, dimassrec, ptPreFSR->at(2), mPreFSR->at(2), weightRec*bTagReweight, weightGen); 

		   // scale systematic
		   for(unsigned int i=0; i<Scale->size(); i++){
                        TString is;
                        is.Form("%d", i);

		        sigHist.FillMigration2DM( pt_unfold, selected, "pt_Scale_"+is,  diptrec, dimassrec, ptPreFSR->at(2), mPreFSR->at(2), weightRec*bTagReweight, weightGen*Scale->at(i)); 
	                sigHist.FillMigration2DM( mass_unfold, selected, "mass_Scale_"+is, diptrec, dimassrec, ptPreFSR->at(2), mPreFSR->at(2), weightRec*bTagReweight, weightGen*Scale->at(i)); 
		   }

		   // alpha s systematic
                   for(unsigned int i=0; i<AlphaS->size(); i++){
                        TString is;
                        is.Form("%d", i);

                        sigHist.FillMigration2DM( pt_unfold, selected, "pt_AlphaS_"+is,  diptrec, dimassrec, ptPreFSR->at(2), mPreFSR->at(2), weightRec*bTagReweight, weightGen*AlphaS->at(i));
                        sigHist.FillMigration2DM( mass_unfold, selected, "mass_AlphaS_"+is, diptrec, dimassrec, ptPreFSR->at(2), mPreFSR->at(2), weightRec*bTagReweight, weightGen*AlphaS->at(i));
                   }
	
	       }// DY to ee or mumu events only
	       //
	       //////////////////////////////////////////////////// fill migration matrix ////////////////////////////////////////
	
	   }// for DYtoLL MC case
	}// event loop

        fileout2->cd();
        recoHist.GetPtBinningRec()->Write();
        recoHist.GetMassBinningRec()->Write();

        fileout1->cd();
        sigHist.GetPtBinningRec()->Write();
        sigHist.GetMassBinningRec()->Write();
        sigHist.GetPtBinningGen()->Write();
        sigHist.GetMassBinningGen()->Write();
	
	delete hZptWeight;
	delete fZptWeight;
	delete tsignal;
	//fileout->cd();
	
	// seems ptHistograms automatically written 
	// recoHist.ptHists.at(0)->Write();
	// delete outputFile;
}

*/
