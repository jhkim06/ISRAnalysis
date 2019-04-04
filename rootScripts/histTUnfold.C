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
