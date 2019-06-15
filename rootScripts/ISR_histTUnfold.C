#include "ISR_histTUnfold.h"

void histTUnfold::SetsysMap(TString sysName, int nVariations){

	sysMaps.insert(std::pair<TString, int>(sysName, nVariations));
}


void histTUnfold::CreateHistMap(int which_unfold, TString hname, TString postfix, bool isRec){


        if(which_unfold == ptOrMass::PT){
                 if(isRec) histMaps.insert(std::pair<TString, TH1*>("pt_"+hname+postfix, ptBinningRec->CreateHistogram("hPtRec"+hname)));
		 else histMaps.insert(std::pair<TString, TH1*>("ptGen_"+hname+postfix, 	 ptBinningGen->CreateHistogram("hPtGen"+hname)));
        }

        if(which_unfold == ptOrMass::MASS){
                 if(isRec) histMaps.insert(std::pair<TString, TH1*>("mass_"+hname+postfix, massBinningRec->CreateHistogram("hMassRec"+hname)));
		 else histMaps.insert(std::pair<TString, TH1*>("massGen_"+hname+postfix, massBinningGen->CreateHistogram("hMassGen"+hname)));
        }
}

void histTUnfold::CreateHist2DMap(int which_unfold, TString hname){

	if(which_unfold == ptOrMass::PT){
		 hist2DMaps.insert(std::pair<TString, TH2*>("pt_"+hname, TUnfoldBinning::CreateHistogramOfMigrations(ptBinningGen, ptBinningRec,"hmcPtGenRec" + hname)));
	}

	if(which_unfold == ptOrMass::MASS){
		 hist2DMaps.insert(std::pair<TString, TH2*>("mass_"+hname, TUnfoldBinning::CreateHistogramOfMigrations(massBinningGen, massBinningRec,"hmcMassGenRec" + hname)));
	}
}

void histTUnfold::FillHistogram(ptOrMass which_unfold, TString hname, Double_t recoPt, Double_t recoMass, Double_t wreco, bool isRec){

	if(which_unfold == ptOrMass::PT){
		if(isRec) (histMaps.find(hname)->second)->Fill(ptBinningRec->GetGlobalBinNumber(recoPt, recoMass),  wreco);
		else (histMaps.find(hname)->second)->Fill(ptBinningGen->GetGlobalBinNumber(recoPt, recoMass),  wreco);
	}

	if(which_unfold == ptOrMass::MASS){
		if(isRec) {
		   if(recoPt < 100.) {
			(histMaps.find(hname)->second)->Fill(massBinningRec->GetGlobalBinNumber(recoMass),  wreco);
		   }
		}
                else {
		   if(recoPt < 100.) {
                        (histMaps.find(hname)->second)->Fill(massBinningGen->GetGlobalBinNumber(recoMass),  wreco);
		   }
                }
			
	}
}

void histTUnfold::FillMigration2DM(ptOrMass which_unfold, bool selected, TString hname, Double_t recoPt, Double_t RecoMass, Double_t truthPt, Double_t truthMass, Double_t wreco, Double_t wgen, Double_t corr){

	int binZero=0;

	if(selected){
		if( which_unfold == ptOrMass::PT){
			(hist2DMaps.find(hname)->second)->Fill(ptBinningGen->GetGlobalBinNumber(truthPt, truthMass), ptBinningRec->GetGlobalBinNumber(recoPt, RecoMass), wgen*wreco*corr);
			(hist2DMaps.find(hname)->second)->Fill(ptBinningGen->GetGlobalBinNumber(truthPt, truthMass), binZero, wgen*(1.-wreco));
		}
                if( which_unfold == ptOrMass::MASS){
                        if(recoPt < 100.) (hist2DMaps.find(hname)->second)->Fill(massBinningGen->GetGlobalBinNumber(truthMass), massBinningRec->GetGlobalBinNumber(RecoMass), wgen*wreco*corr);
                        if(recoPt < 100.) (hist2DMaps.find(hname)->second)->Fill(massBinningGen->GetGlobalBinNumber(truthMass), binZero, wgen*(1.-wreco));
                }
	}
	else{
		if( which_unfold == ptOrMass::PT) (hist2DMaps.find(hname)->second)->Fill(ptBinningGen->GetGlobalBinNumber(truthPt, truthMass), binZero, wgen);
		if( which_unfold == ptOrMass::MASS) (hist2DMaps.find(hname)->second)->Fill(massBinningGen->GetGlobalBinNumber(truthMass), binZero, wgen);
	}
}

void histTUnfold::saveRecoHists(TFile *filein, TFile *fileout1, TString channel){ // 

        gROOT->SetBatch();
        gInterpreter->GenerateDictionary("vector<TLorentzVector>","TLorentzVector.h;vector");
        TH1::SetDefaultSumw2();

 	ptRec = 0;
        mRec = 0; 	
 	ptRec_momentumUp = 0;
        mRec_momentumUp = 0; 	
 	ptRec_momentumDown = 0;
        mRec_momentumDown = 0; 	
 	ptRec_momentumResUp = 0;
        mRec_momentumResUp = 0; 	
 	ptRec_momentumResDown = 0;
        mRec_momentumResDown = 0; 	
        AlphaS = 0;
        Scale = 0;
        PDFerror = 0;

        ZPtCor = -999.;

        tree=(TTree *)filein->Get("tree");
        tree->SetBranchAddress("ptRec",&ptRec);
        tree->SetBranchAddress("mRec",&mRec);
        tree->SetBranchAddress("ptRec_momentumUp",&ptRec_momentumUp);
        tree->SetBranchAddress("mRec_momentumUp",&mRec_momentumUp);
        tree->SetBranchAddress("ptRec_momentumDown",&ptRec_momentumDown);
        tree->SetBranchAddress("mRec_momentumDown",&mRec_momentumDown);
        tree->SetBranchAddress("ptRec_momentumResUp",&ptRec_momentumResUp);
        tree->SetBranchAddress("mRec_momentumResUp",&mRec_momentumResUp);
        tree->SetBranchAddress("ptRec_momentumResDown",&ptRec_momentumResDown);
        tree->SetBranchAddress("mRec_momentumResDown",&mRec_momentumResDown);
        tree->SetBranchAddress("IsElEl",&IsElEl);
        tree->SetBranchAddress("IsMuMu",&IsMuMu);
        tree->SetBranchAddress("ispassRec",&ispassRec);
        tree->SetBranchAddress("isBveto",&isBveto);
        tree->SetBranchAddress("weightGen",&weightGen);
        tree->SetBranchAddress("weightRec",&weightRec);
        tree->SetBranchAddress("DYtautau",&DYtautau);
        tree->SetBranchAddress("bTagReweight",&bTagReweight); // FIXME may it is better to change name bTagReweight to bTagSF

        tree->SetBranchAddress("PUweight", &PUweight);
        tree->SetBranchAddress("PUweight_Up", &PUweight_Up);
        tree->SetBranchAddress("PUweight_Dn", &PUweight_Dn);
        tree->SetBranchAddress("trgSF", &trgSF);
        tree->SetBranchAddress("trgSF_Up", &trgSF_Up);
        tree->SetBranchAddress("trgSF_Dn", &trgSF_Dn);
        tree->SetBranchAddress("recoSF", &recoSF);
        tree->SetBranchAddress("recoSF_Up", &recoSF_Up);
        tree->SetBranchAddress("recoSF_Dn", &recoSF_Dn);
        tree->SetBranchAddress("IdSF", &IdSF);
        tree->SetBranchAddress("IdSF_Up", &IdSF_Up);
        tree->SetBranchAddress("IdSF_Dn", &IdSF_Dn);
        tree->SetBranchAddress("IsoSF", &IsoSF);
        tree->SetBranchAddress("IsoSF_Up", &IsoSF_Up);
        tree->SetBranchAddress("IsoSF_Dn", &IsoSF_Dn);
        tree->SetBranchAddress("L1Prefire",    &L1Prefire);
        tree->SetBranchAddress("L1Prefire_Up", &L1Prefire_Up);
        tree->SetBranchAddress("L1Prefire_Dn", &L1Prefire_Dn);

	TBranch* b_alphas =(TBranch*) tree->GetListOfBranches()->FindObject("AlphaS");
	TBranch* b_scale =(TBranch*) tree->GetListOfBranches()->FindObject("Scale");
	TBranch* b_pdferror =(TBranch*) tree->GetListOfBranches()->FindObject("PDFerror");
	TBranch* b_zptcorr =(TBranch*) tree->GetListOfBranches()->FindObject("ZPtCor");

        if(b_alphas)   tree->SetBranchAddress("AlphaS",&AlphaS);
        if(b_scale)    tree->SetBranchAddress("Scale",&Scale);
        if(b_pdferror) tree->SetBranchAddress("PDFerror",&PDFerror);
        if(b_zptcorr)  tree->SetBranchAddress("ZPtCor",&ZPtCor);	

        nentries=tree->GetEntries();

	double lep1_ptCut, lep2_ptCut;

	for(int i=0;i<nentries;i++){
        //for(int i=0;i<10000;i++){
          if(i%10000000==0) cout<<i<<endl;
          tree->GetEntry(i);

          int isdilep = 0;
          if( channel.CompareTo("electron") == 0 ){ 
		isdilep = IsElEl;
		lep1_ptCut = 25.;
		lep2_ptCut = 15.;
	  }
          if( channel.CompareTo("muon") == 0 ){ 
		isdilep = IsMuMu;
		lep1_ptCut = 25.;
		lep2_ptCut = 15.;
	  }

	   // nominal dilepton event selection 
           if(isdilep && ispassRec && isBveto && ptRec->at(0) > lep1_ptCut && ptRec->at(1) > lep2_ptCut){
                 fileout1->cd();

		 Double_t wreco = weightRec*bTagReweight;
		 Double_t wgen = weightGen;
                 FillHistogram(ptOrMass::PT, "pt_nominal", ptRec->at(2), mRec->at(2), wgen*wreco);
                 FillHistogram(ptOrMass::MASS, "mass_nominal", ptRec->at(2), mRec->at(2), wgen*wreco);

                 // loop to create systematic response matrix
                 std::map<TString, int>::iterator it = sysMaps.begin();
                 while(it != sysMaps.end()){
                      if( it->first != "LepScale" && it->first != "LepRes"){ // TODO also skip lepton scale/resolution systematic 

                         for(int nSys = 0; nSys < it->second; nSys++){
                             TString is;
                             is.Form("%d", nSys);

                             wreco = GetSysWeights(it->first, true,  nSys);
                             wgen = GetSysWeights(it->first, false,  nSys);

                             FillHistogram(ptOrMass::PT,   "pt_"  + it->first+"_"+is,    ptRec->at(2), mRec->at(2), wgen*wreco);
                             FillHistogram(ptOrMass::MASS, "mass_"  + it->first+"_"+is,  ptRec->at(2), mRec->at(2), wgen*wreco);

                         }
                      }
                      it++;
                 }// loop for systematic

           }// 

	   // dilepton event selection with lepton momentum scale up
	   if(isdilep && ispassRec && isBveto && ptRec_momentumUp->at(0) > lep1_ptCut && ptRec_momentumUp->at(1) > lep2_ptCut){
                 fileout1->cd();

                 Double_t wreco = weightRec*bTagReweight;
                 Double_t wgen = weightGen;
                 FillHistogram(ptOrMass::PT,   "pt_LepScale_0",     ptRec_momentumUp->at(2), mRec_momentumUp->at(2), wgen*wreco);
                 FillHistogram(ptOrMass::MASS, "mass_LepScale_0",   ptRec_momentumUp->at(2), mRec_momentumUp->at(2), wgen*wreco);
	   }

           // dilepton event selection with lepton momentum scale down 
           if(isdilep && ispassRec && isBveto && ptRec_momentumDown->at(0) > lep1_ptCut && ptRec_momentumDown->at(1) > lep2_ptCut){
                 fileout1->cd();

                 Double_t wreco = weightRec*bTagReweight;
                 Double_t wgen = weightGen;
                 FillHistogram(ptOrMass::PT,   "pt_LepScale_1",     ptRec_momentumDown->at(2), mRec_momentumDown->at(2), wgen*wreco);
                 FillHistogram(ptOrMass::MASS, "mass_LepScale_1",   ptRec_momentumDown->at(2), mRec_momentumDown->at(2), wgen*wreco);
           }

           // dilepton event selection with lepton momentum resolution up
           if(isdilep && ispassRec && isBveto && ptRec_momentumResUp->at(0) > lep1_ptCut && ptRec_momentumResUp->at(1) > lep2_ptCut){
                 fileout1->cd();

                 Double_t wreco = weightRec*bTagReweight;
                 Double_t wgen = weightGen;
                 FillHistogram(ptOrMass::PT,   "pt_LepRes_0",     ptRec_momentumResUp->at(2), mRec_momentumResUp->at(2), wgen*wreco);
                 FillHistogram(ptOrMass::MASS, "mass_LepRes_0",   ptRec_momentumResUp->at(2), mRec_momentumResUp->at(2), wgen*wreco);
           }   

           // dilepton event selection with lepton momentum resolution down
           if(isdilep && ispassRec && isBveto && ptRec_momentumResDown->at(0) > lep1_ptCut && ptRec_momentumResDown->at(1) > lep2_ptCut){
                 fileout1->cd();

                 Double_t wreco = weightRec*bTagReweight;
                 Double_t wgen = weightGen;
                 FillHistogram(ptOrMass::PT,   "pt_LepRes_1",     ptRec_momentumResDown->at(2), mRec_momentumResDown->at(2), wgen*wreco);
                 FillHistogram(ptOrMass::MASS, "mass_LepRes_1",   ptRec_momentumResDown->at(2), mRec_momentumResDown->at(2), wgen*wreco);
           }   



        }// event loop

        fileout1->cd();
        GetPtBinningRec()->Write();
        GetMassBinningRec()->Write();

        delete tree;
        //fileout->cd();
}

Double_t histTUnfold::GetSysWeights(TString sysName, bool isReco, int nthSys){

	if(sysName == "PU"){
		if(isReco){
			if(nthSys == 0) return weightRec*bTagReweight*PUweight_Up/PUweight;
			if(nthSys == 1) return weightRec*bTagReweight*PUweight_Dn/PUweight;
		}
		else{
			return weightGen;
		}
	}// PU

        else if(sysName == "trgSF"){
                if(isReco){
                        if(nthSys == 0) return weightRec*bTagReweight*trgSF_Up/trgSF;
                        if(nthSys == 1) return weightRec*bTagReweight*trgSF_Dn/trgSF;
                }
                else{
                        return weightGen;
                }
        }// trgSF

        else if(sysName == "recoSF"){
                if(isReco){
                        if(nthSys == 0) return weightRec*bTagReweight*recoSF_Up/recoSF;
                        if(nthSys == 1) return weightRec*bTagReweight*recoSF_Dn/recoSF;
                }   
                else{
                        return weightGen;
                }   
        }// recoSF

        else if(sysName == "IdSF"){
                if(isReco){
                        if(nthSys == 0) return weightRec*bTagReweight*IdSF_Up/IdSF;
                        if(nthSys == 1) return weightRec*bTagReweight*IdSF_Dn/IdSF;
                }   
                else{
                        return weightGen;
                }   
        }// IdSF

        else if(sysName == "IsoSF"){
                if(isReco){
                        if(nthSys == 0) return weightRec*bTagReweight*IsoSF_Up/IsoSF;
                        if(nthSys == 1) return weightRec*bTagReweight*IsoSF_Dn/IsoSF;
                }   
                else{
                        return weightGen;
                }   
        }// IsoSF

        else if(sysName == "L1Prefire"){
                if(isReco){
                        if(nthSys == 0) return weightRec*bTagReweight*L1Prefire_Up/L1Prefire;
                        if(nthSys == 1) return weightRec*bTagReweight*L1Prefire_Dn/L1Prefire;
                }
                else{
                        return weightGen;
                }
        }// L1Prefire

        else if(sysName == "AlphaS"){
                if(isReco){
			return weightRec*bTagReweight;
                }   
                else{
                        if(AlphaS != 0 && AlphaS->size() != 0) return weightGen*AlphaS->at(nthSys);
			else return weightGen;
                }   
        }// AlphaS
        else if(sysName == "Scale"){
                if(isReco){
			return weightRec*bTagReweight;
                }
                else{
                        if(Scale != 0 && Scale->size() != 0 && nthSys < Scale->size()) return weightGen*Scale->at(nthSys);
			else return weightGen;
                }
        }// Scale
        else if(sysName == "PDFerror"){
                if(isReco){
			return weightRec*bTagReweight;
                }
                else{
                        if(PDFerror != 0 && PDFerror->size() != 0) return weightGen*PDFerror->at(nthSys);
			else return weightGen;
                }
        }// PDFerror
        else if(sysName == "unfoldsys"){
                if(isReco){
			return weightRec*bTagReweight;
                }
                else{
                        if(ZPtCor != -999) return weightGen*ZPtCor;
			else return weightGen;
                }
        }// unfoldsys
        else if(sysName == "FSRDR"){
                if(isReco){
                        return weightRec*bTagReweight;
                }
                else{
                        return weightGen;
                }
        }// FSRDR
	else{
		return 1.;
        }

}

void histTUnfold::saveSigHists(TFile *filein, TFile *fileout1, TFile *fileout2, TString channel, Double_t temp_kfactor){ // 

        gROOT->SetBatch();
        gInterpreter->GenerateDictionary("vector<TLorentzVector>","TLorentzVector.h;vector");
        TH1::SetDefaultSumw2();

        ptPreFSR = 0;
        mPreFSR = 0;
        ptPostFSR = 0;
        mPostFSR = 0;
        ptRec = 0;
        mRec = 0;
        ptRec_momentumUp = 0;
        mRec_momentumUp = 0;
        ptRec_momentumDown = 0;
        mRec_momentumDown = 0;
        ptRec_momentumResUp = 0;
        mRec_momentumResUp = 0;
        ptRec_momentumResDown = 0;
        mRec_momentumResDown = 0;

        AlphaS = 0;
        Scale = 0;
        PDFerror = 0;

        particleFSR = 0;
        anparticleFSR = 0;
        particlePostFSR = 0;
        anparticlePostFSR = 0;

        tree=(TTree *)filein->Get("tree");
        tree->SetBranchAddress("ptPreFSR",&ptPreFSR);
        tree->SetBranchAddress("mPreFSR",&mPreFSR);
        tree->SetBranchAddress("ptPostFSR",&ptPostFSR);
        tree->SetBranchAddress("mPostFSR",&mPostFSR);
        tree->SetBranchAddress("ptRec",&ptRec);
        tree->SetBranchAddress("mRec",&mRec);
        tree->SetBranchAddress("ptRec_momentumUp",&ptRec_momentumUp);
        tree->SetBranchAddress("mRec_momentumUp",&mRec_momentumUp);
        tree->SetBranchAddress("ptRec_momentumDown",&ptRec_momentumDown);
        tree->SetBranchAddress("mRec_momentumDown",&mRec_momentumDown);
        tree->SetBranchAddress("ptRec_momentumResUp",&ptRec_momentumResUp);
        tree->SetBranchAddress("mRec_momentumResUp",&mRec_momentumResUp);
        tree->SetBranchAddress("ptRec_momentumResDown",&ptRec_momentumResDown);
        tree->SetBranchAddress("mRec_momentumResDown",&mRec_momentumResDown);
        tree->SetBranchAddress("IsElEl",&IsElEl);
        tree->SetBranchAddress("IsMuMu",&IsMuMu);
        tree->SetBranchAddress("ispassRec",&ispassRec);
        tree->SetBranchAddress("isBveto",&isBveto);
        tree->SetBranchAddress("isdimuon",&isdimuon);
        tree->SetBranchAddress("isdielectron",&isdielectron);
        tree->SetBranchAddress("weightGen",&weightGen);
        tree->SetBranchAddress("weightRec",&weightRec);
        tree->SetBranchAddress("ZPtCor",&ZPtCor);
        tree->SetBranchAddress("bTagReweight",&bTagReweight);
        tree->SetBranchAddress("DYtautau",&DYtautau);
        tree->SetBranchAddress("isfiducialPreFSR",&isfiducialPreFSR);
        tree->SetBranchAddress("isfiducialPostFSR",&isfiducialPostFSR);

	if(!isAlt){
        	tree->SetBranchAddress("particleFSR", &particleFSR, &b_particleFSR);
        	tree->SetBranchAddress("anparticleFSR", &anparticleFSR, &b_anparticleFSR);
        	tree->SetBranchAddress("particlePostFSR", &particlePostFSR, &b_particlePostFSR);
        	tree->SetBranchAddress("anparticlePostFSR", &anparticlePostFSR, &b_anparticlePostFSR);

        	tree->SetBranchAddress("AlphaS",&AlphaS);
        	tree->SetBranchAddress("Scale",&Scale);
        	tree->SetBranchAddress("PDFerror",&PDFerror);

        	tree->SetBranchAddress("PUweight", &PUweight);
        	tree->SetBranchAddress("PUweight_Up", &PUweight_Up);
        	tree->SetBranchAddress("PUweight_Dn", &PUweight_Dn);
        	tree->SetBranchAddress("trgSF", &trgSF);
        	tree->SetBranchAddress("trgSF_Up", &trgSF_Up);
        	tree->SetBranchAddress("trgSF_Dn", &trgSF_Dn);
        	tree->SetBranchAddress("recoSF", &recoSF);
        	tree->SetBranchAddress("recoSF_Up", &recoSF_Up);
        	tree->SetBranchAddress("recoSF_Dn", &recoSF_Dn);
        	tree->SetBranchAddress("IdSF", &IdSF);
        	tree->SetBranchAddress("IdSF_Up", &IdSF_Up);
        	tree->SetBranchAddress("IdSF_Dn", &IdSF_Dn);
        	tree->SetBranchAddress("IsoSF", &IsoSF);
        	tree->SetBranchAddress("IsoSF_Up", &IsoSF_Up);
        	tree->SetBranchAddress("IsoSF_Dn", &IsoSF_Dn);
        	tree->SetBranchAddress("L1Prefire",    &L1Prefire);
        	tree->SetBranchAddress("L1Prefire_Up", &L1Prefire_Up);
        	tree->SetBranchAddress("L1Prefire_Dn", &L1Prefire_Dn);
	}

	TH2* hMCGenGenPt =   TUnfoldBinning::CreateHistogramOfMigrations(ptBinningGen,ptBinningGen,    "histMCGenGenPt");
	TH2* hMCGenGenMass = TUnfoldBinning::CreateHistogramOfMigrations(massBinningGen,massBinningGen,"histMCGenGenMass");

        nentries=tree->GetEntries();

        double lep1_ptCut, lep2_ptCut;

        for(int i=0;i<nentries;i++){
        //for(int i=0;i<10000;i++){
          if(i%10000000==0) cout<<i<<endl;
          tree->GetEntry(i);

          weightGen *= temp_kfactor;

          int isdilep = 0, issignal = 0;
          if( channel.CompareTo("electron") == 0 ){
                  isdilep = IsElEl;
                  issignal = isdielectron;
                  lep1_ptCut = 25.;
                  lep2_ptCut = 15.;
          }
          if( channel.CompareTo("muon") == 0 ){
                  isdilep = IsMuMu;
                  issignal = isdimuon;
		  lep1_ptCut = 25.;
		  lep2_ptCut = 15.;
          }

          if( std::isinf(bTagReweight) ){
                 bTagReweight = 1.; //FIXME check analyzer later
                 std::cout << "nan bTagReweight... check later" << std::endl;
          }

           if(!isInc){
              if(isdilep && ispassRec){
                    fileout1->cd();

                    FillHistogram(ptOrMass::PT, "pt_nominal", ptRec->at(2), mRec->at(2), weightGen*weightRec*bTagReweight);
                    FillHistogram(ptOrMass::MASS, "mass_nominal", ptRec->at(2), mRec->at(2), weightGen*weightRec*bTagReweight);

              }
           }
           else{ // for inclusive DY MC

                Double_t zptWeight = 1.;
                zptWeight = ZPtCor;

                if(isdilep && ispassRec && isBveto && ptRec->at(0) > lep1_ptCut && ptRec->at(1) > lep2_ptCut){

		   TString postfix = "";
		   
                   if(DYtautau){
		      if(isAlt) continue;
                      fileout2->cd();
	  	      postfix = "_tau";
		   }
		   else{
			fileout1->cd();
	           }

		   Double_t wreco = weightRec*bTagReweight;
		   Double_t wgen = weightGen;

		   if( wreco == 0 ) continue; // FIXME check why zero value is saved (check L1Prefire for muon channel)

		   if(!isAlt){
						    //map name of histogram	
                   	FillHistogram(ptOrMass::PT, "pt_nominal"+postfix,     ptRec->at(2), mRec->at(2), wgen*wreco);
                   	FillHistogram(ptOrMass::MASS, "mass_nominal"+postfix, ptRec->at(2), mRec->at(2), wgen*wreco);

                   	// loop to create systematic response matrix
                   	std::map<TString, int>::iterator it = sysMaps.begin();
                   	while(it != sysMaps.end()){
                   	     if( it->first != "LepScale" && it->first != "LepRes"){ // 

                   	        for(int nSys = 0; nSys < it->second; nSys++){
                   	            TString is;
                   	            is.Form("%d", nSys);

                   	            wreco = GetSysWeights(it->first, true,  nSys);
                   	            wgen = GetSysWeights(it->first, false,  nSys);

                   	            FillHistogram(ptOrMass::PT,   "pt_"  + it->first+"_"+is+postfix,    ptRec->at(2), mRec->at(2), wgen*wreco);
                   	            FillHistogram(ptOrMass::MASS, "mass_"  + it->first+"_"+is+postfix,  ptRec->at(2), mRec->at(2), wgen*wreco);

                   	        }
                   	     }
                   	     it++;
                   	}// loop for systematic
	           }

               }
		if(!isAlt){

                   TString postfix = "";

                   if(DYtautau){
                      if(isAlt) continue;
                      fileout2->cd();
                      postfix = "_tau";
                   }
                   else{
                        fileout1->cd();
                   }

                   Double_t wreco = weightRec*bTagReweight;
                   Double_t wgen = weightGen;

                   // dilepton event selection with lepton momentum scale up
                   if(isdilep && ispassRec && isBveto && ptRec_momentumUp->at(0) > lep1_ptCut && ptRec_momentumUp->at(1) > lep2_ptCut){
                         fileout1->cd();

                         Double_t wreco = weightRec*bTagReweight;
                         Double_t wgen = weightGen;
                         FillHistogram(ptOrMass::PT,   "pt_LepScale_0"+postfix,     ptRec_momentumUp->at(2), mRec_momentumUp->at(2), wgen*wreco);
                         FillHistogram(ptOrMass::MASS, "mass_LepScale_0"+postfix,   ptRec_momentumUp->at(2), mRec_momentumUp->at(2), wgen*wreco);
                   }   

                   // dilepton event selection with lepton momentum scale down 
                   if(isdilep && ispassRec && isBveto && ptRec_momentumDown->at(0) > lep1_ptCut && ptRec_momentumDown->at(1) > lep2_ptCut){
                         fileout1->cd();

                         Double_t wreco = weightRec*bTagReweight;
                         Double_t wgen = weightGen;
                         FillHistogram(ptOrMass::PT,   "pt_LepScale_1"+postfix,     ptRec_momentumDown->at(2), mRec_momentumDown->at(2), wgen*wreco);
                         FillHistogram(ptOrMass::MASS, "mass_LepScale_1"+postfix,   ptRec_momentumDown->at(2), mRec_momentumDown->at(2), wgen*wreco);
                   }   

                   // dilepton event selection with lepton momentum resolution up
                   if(isdilep && ispassRec && isBveto && ptRec_momentumResUp->at(0) > lep1_ptCut && ptRec_momentumResUp->at(1) > lep2_ptCut){
                         fileout1->cd();

                         Double_t wreco = weightRec*bTagReweight;
                         Double_t wgen = weightGen;
                         FillHistogram(ptOrMass::PT,   "pt_LepRes_0"+postfix,     ptRec_momentumResUp->at(2), mRec_momentumResUp->at(2), wgen*wreco);
                         FillHistogram(ptOrMass::MASS, "mass_LepRes_0"+postfix,   ptRec_momentumResUp->at(2), mRec_momentumResUp->at(2), wgen*wreco);
                   }   

                   // dilepton event selection with lepton momentum resolution down
                   if(isdilep && ispassRec && isBveto && ptRec_momentumResDown->at(0) > lep1_ptCut && ptRec_momentumResDown->at(1) > lep2_ptCut){
                         fileout1->cd();

                         Double_t wreco = weightRec*bTagReweight;
                         Double_t wgen = weightGen;
                         FillHistogram(ptOrMass::PT,   "pt_LepRes_1"+postfix,     ptRec_momentumResDown->at(2), mRec_momentumResDown->at(2), wgen*wreco);
                         FillHistogram(ptOrMass::MASS, "mass_LepRes_1"+postfix,   ptRec_momentumResDown->at(2), mRec_momentumResDown->at(2), wgen*wreco);
                   }   
		}


                
               /////////////////////////////////////////// fill migration matrix ////////////////////////////////////////
               //
               if(issignal){
                   fileout1->cd();

                   Double_t diptrec = -999., dimassrec = -999.;
                   Double_t diptgen = -999., dimassgen = -999.;
			
		   // systematic dilepton pt and mass
                   Double_t diptrec_ScaleUp = -999., dimassrec_ScaleUp = -999.;
                   Double_t diptrec_ScaleDn = -999., dimassrec_ScaleDn = -999.;
                   Double_t diptrec_ResUp = -999., dimassrec_ResUp = -999.;
                   Double_t diptrec_ResDn = -999., dimassrec_ResDn = -999.;


                   bool selected = (isdilep && ispassRec && isBveto && ptRec->at(0) > lep1_ptCut && ptRec->at(1) > lep2_ptCut);

                   // check if there are selected dilepton at detector to avoid memory error
                   if(ptRec->size() == 3){ 
			diptrec   = ptRec->at(2);
	                dimassrec = mRec->at(2);
		   }

                   diptgen =   ptPreFSR->at(2);
                   dimassgen = mPreFSR->at(2);

                   FillMigration2DM( ptOrMass::PT, selected,   "pt_nominal",   diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight, weightGen);
                   FillMigration2DM( ptOrMass::MASS, selected, "mass_nominal", diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight, weightGen);
		   
                   FillHistogram(ptOrMass::PT,   "ptGen_nominal",    diptgen, dimassgen, weightGen, false);
                   FillHistogram(ptOrMass::MASS, "massGen_nominal",  diptgen, dimassgen, weightGen, false);

                   if(!isAlt){			
                   	bool selected_ScaleUp = (isdilep && ispassRec && isBveto && ptRec_momentumUp->at(0) 	  >   lep1_ptCut && ptRec_momentumUp->at(1) > lep2_ptCut);
                   	bool selected_ScaleDn = (isdilep && ispassRec && isBveto && ptRec_momentumDown->at(0)    >   lep1_ptCut && ptRec_momentumDown->at(1) > lep2_ptCut);
                   	bool selected_ResUp =   (isdilep && ispassRec && isBveto && ptRec_momentumResUp->at(0) >   lep1_ptCut && ptRec_momentumResUp->at(1) > lep2_ptCut);
                   	bool selected_ResDn =   (isdilep && ispassRec && isBveto && ptRec_momentumResDown->at(0) >   lep1_ptCut && ptRec_momentumResDown->at(1) > lep2_ptCut);

                   	if(ptRec->size() == 3){ 
		   	     diptrec_ScaleUp   = ptRec_momentumUp->at(2);
	           	     dimassrec_ScaleUp = mRec_momentumUp->at(2);

		   	     diptrec_ScaleDn   = ptRec_momentumDown->at(2);
	           	     dimassrec_ScaleDn = mRec_momentumDown->at(2);

		   	     diptrec_ResUp   = ptRec_momentumResUp->at(2);
	           	     dimassrec_ResUp = mRec_momentumResUp->at(2);

		   	     diptrec_ResDn   = ptRec_momentumResDown->at(2);
	           	     dimassrec_ResDn = mRec_momentumResDown->at(2);
		   	}

                   	FillMigration2DM( ptOrMass::PT, selected,   "pt_LepScale_0",   diptrec_ScaleUp, dimassrec_ScaleUp, diptgen, dimassgen, weightRec*bTagReweight, weightGen);
                   	FillMigration2DM( ptOrMass::MASS, selected, "mass_LepScale_0", diptrec_ScaleUp, dimassrec_ScaleUp, diptgen, dimassgen, weightRec*bTagReweight, weightGen);

                   	FillMigration2DM( ptOrMass::PT, selected,   "pt_LepScale_1",   diptrec_ScaleDn, dimassrec_ScaleDn, diptgen, dimassgen, weightRec*bTagReweight, weightGen);
                   	FillMigration2DM( ptOrMass::MASS, selected, "mass_LepScale_1", diptrec_ScaleDn, dimassrec_ScaleDn, diptgen, dimassgen, weightRec*bTagReweight, weightGen);

                   	FillMigration2DM( ptOrMass::PT, selected,   "pt_LepRes_0",   diptrec_ResUp, dimassrec_ResUp, diptgen, dimassgen, weightRec*bTagReweight, weightGen);
                   	FillMigration2DM( ptOrMass::MASS, selected, "mass_LepRes_0", diptrec_ResUp, dimassrec_ResUp, diptgen, dimassgen, weightRec*bTagReweight, weightGen);

                   	FillMigration2DM( ptOrMass::PT, selected,   "pt_LepRes_1",   diptrec_ResDn, dimassrec_ResDn, diptgen, dimassgen, weightRec*bTagReweight, weightGen);
                   	FillMigration2DM( ptOrMass::MASS, selected, "mass_LepRes_1", diptrec_ResDn, dimassrec_ResDn, diptgen, dimassgen, weightRec*bTagReweight, weightGen);

                   	TLorentzVector temp_particlePostFSR[30];
                   	TLorentzVector temp_anparticlePostFSR[30];
                   	Double_t drcut[30] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.};

                   	// add fsr photon in dr < 0.1 for lepton
                   	for(int idrcut = 0; idrcut < 30; idrcut++){
                   	   temp_particlePostFSR[idrcut] = *particlePostFSR;
                   	   for(int gsize = 0; gsize < particleFSR->size(); gsize++){

                   	           Double_t temp_dr = sqrt(pow(particlePostFSR->Phi()-particleFSR->at(gsize).Phi(),2)+pow(particlePostFSR->Eta()-particleFSR->at(gsize).Eta(),2));
                   	           if(temp_dr < drcut[idrcut]) temp_particlePostFSR[idrcut] += particleFSR->at(gsize);
                   	   }
                   	}

                   	for(int idrcut = 0; idrcut < 30; idrcut++){
                   	   temp_anparticlePostFSR[idrcut] = *anparticlePostFSR;
                   	   for(int gsize = 0; gsize < anparticleFSR->size(); gsize++){

                   	     Double_t temp_dr = sqrt(pow(anparticlePostFSR->Phi()-anparticleFSR->at(gsize).Phi(),2)+pow(anparticlePostFSR->Eta()-anparticleFSR->at(gsize).Eta(),2));
                   	     if(temp_dr < drcut[idrcut]) temp_anparticlePostFSR[idrcut] += anparticleFSR->at(gsize);
                   	   }
                   	}

                   	for(int idrcut = 0; idrcut < 30; idrcut++){
                   	   Double_t temp_dipt = (temp_particlePostFSR[idrcut] + temp_anparticlePostFSR[idrcut]).Pt();
                   	   Double_t temp_dimass = (temp_particlePostFSR[idrcut] + temp_anparticlePostFSR[idrcut]).M();
                   	   TString dr_;
                           dr_.Form("%d", idrcut);
                           FillMigration2DM( ptOrMass::PT,   selected, "pt_FSRDR_"+dr_,   diptrec, dimassrec, temp_dipt, temp_dimass, weightRec*bTagReweight, weightGen);
                           FillMigration2DM( ptOrMass::MASS, selected, "mass_FSRDR_"+dr_, diptrec, dimassrec, temp_dipt, temp_dimass, weightRec*bTagReweight, weightGen);

			   if(idrcut ==0){
	 			hMCGenGenPt->  Fill(ptBinningGen->GetGlobalBinNumber(diptgen, dimassgen), ptBinningGen->GetGlobalBinNumber(temp_dipt, temp_dimass), weightGen);
	 			hMCGenGenMass->Fill(massBinningGen->GetGlobalBinNumber(dimassgen),        massBinningGen->GetGlobalBinNumber(temp_dimass), weightGen);
			   }
                   	}

		   	// loop to create systematic response matrix
		   	std::map<TString, int>::iterator it = sysMaps.begin();
		   	while(it != sysMaps.end()){
		   	     if(it->first != "FSRDR"){ // 

		   	        for(int nSys = 0; nSys < it->second; nSys++){
                   	            TString is;
                   	            is.Form("%d", nSys);

		   	            Double_t wreco = GetSysWeights(it->first, true,  nSys);
		   	            Double_t wgen  = GetSysWeights(it->first, false, nSys);

                   	            FillMigration2DM( ptOrMass::PT,  selected,   "pt_"  + it->first+"_"+is,  diptrec, dimassrec, diptgen, dimassgen, wreco, wgen);
                   	            FillMigration2DM( ptOrMass::MASS, selected,  "mass_"+ it->first+"_"+is,  diptrec, dimassrec, diptgen, dimassgen, wreco, wgen);

 		   	        }
 		   	     }
                   	     it++;
 		   	}// loop for systematic
		   }

               }// DY to ee or mumu events only
               //
               //////////////////////////////////////////////////// fill migration matrix ////////////////////////////////////////

           }// for DYtoLL MC case
        }// event loop
        if(!isAlt){
		fileout2->cd();
        	GetPtBinningRec()->Write();
        	GetMassBinningRec()->Write();
	}

        fileout1->cd();
        GetPtBinningRec()->Write();
        GetMassBinningRec()->Write();
        GetPtBinningGen()->Write();
        GetMassBinningGen()->Write();
	hMCGenGenPt->Write();
	hMCGenGenMass->Write();

        delete tree;
        //fileout->cd();

        // seems ptHistograms automatically written 
        // recoHist.ptHists.at(0)->Write();
        // delete outputFile;
}


void histTUnfold::SetPtBinningRec(TString channel){

 const int nptbin_fine=17;
 double ptbin_fine[nptbin_fine+1]={0., 2., 4., 6., 8., 10., 12., 14., 18., 22., 28., 35., 45., 55., 65., 75., 85., 100.};

 ptBinningRec= new TUnfoldBinning("Rec_Pt");
 ptBinningRec->AddAxis("pt",nptbin_fine,ptbin_fine,false,true);
 if( channel.CompareTo("electron") == 0 ) ptBinningRec->AddAxis("mass", nmassBins_forPt, massBins_forPt_electron, true, true);
 if( channel.CompareTo("muon") == 0 )     ptBinningRec->AddAxis("mass", nmassBins_forPt, massBins_forPt_muon, true, true);
}

void histTUnfold::SetPtBinningGen(TString channel){

 const int nptbin_wide=9;
 double ptbin_wide[nptbin_wide+1]={0., 4., 8., 12., 18., 28., 40., 55., 80., 100.};

 //const int nptbin_wide=17;
 //double ptbin_wide[nptbin_wide+1]={0., 2., 4., 6., 8., 10., 12., 14., 18., 22., 28., 35., 45., 55., 65., 75., 85., 100.};

 ptBinningGen=(new TUnfoldBinning("Gen_Pt"));
 ptBinningGen->AddAxis("pt",nptbin_wide,ptbin_wide,false,true);
 if( channel.CompareTo("electron") == 0 ) ptBinningGen->AddAxis("mass", nmassBins_forPt, massBins_forPt_electron, true, true);
 if( channel.CompareTo("muon") == 0 )     ptBinningGen->AddAxis("mass", nmassBins_forPt, massBins_forPt_muon, true, true);

}

void histTUnfold::SetMassBinningRec(TString channel){

 massBinningRec=(new TUnfoldBinning("Rec_Mass"));
 if( channel.CompareTo("electron") == 0 ) massBinningRec->AddAxis("reco mass", nmassBins_fine_electron, massBins_fine_electron, false, false);
 if( channel.CompareTo("muon") == 0 )     massBinningRec->AddAxis("reco mass", nmassBins_fine_muon,     massBins_fine_muon, false, false);
}

void histTUnfold::SetMassBinningGen(TString channel){

 massBinningGen=(new TUnfoldBinning("Gen_Mass"));
 if( channel.CompareTo("electron") == 0 ) massBinningGen->AddAxis("gen mass", nmassBins_wide_electron, massBins_wide_electron, true, true);
 if( channel.CompareTo("muon") == 0 )     massBinningGen->AddAxis("gen mass", nmassBins_wide_muon,     massBins_wide_muon,     true, true);
}
