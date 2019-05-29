#include "ISR_histTUnfold.h"

void histTUnfold::SetsysMap(TString sysName, int nVariations){

	sysMaps.insert(std::pair<TString, int>(sysName, nVariations));
}


void histTUnfold::CreateHistMap(int which_unfold, TString hname, TString postfix){


        if(which_unfold == ptOrMass::PT){
                 histMaps.insert(std::pair<TString, TH1*>("pt_"+hname+postfix, ptBinningRec->CreateHistogram("hPtRec"+hname)));
        }

        if(which_unfold == ptOrMass::MASS){
                 histMaps.insert(std::pair<TString, TH1*>("mass_"+hname+postfix, massBinningRec->CreateHistogram("hMassRec"+hname)));
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

void histTUnfold::FillHistogram(ptOrMass which_unfold, TString hname, Double_t recoPt, Double_t recoMass, Double_t wreco){

	if(which_unfold == ptOrMass::PT){
		(histMaps.find(hname)->second)->Fill(ptBinningRec->GetGlobalBinNumber(recoPt, recoMass),  wreco);
	}

	if(which_unfold == ptOrMass::MASS){
		if(recoPt < 100.) (histMaps.find(hname)->second)->Fill(massBinningRec->GetGlobalBinNumber(recoMass),  wreco);
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

void histTUnfold::saveRecoHists(TFile *filein, TFile *fileout1, TString channel){ // TODO add list of systematics

        gROOT->SetBatch();
        gInterpreter->GenerateDictionary("vector<TLorentzVector>","TLorentzVector.h;vector");
        TH1::SetDefaultSumw2();

 	ptRec = 0;
        mRec = 0; 	

        tree=(TTree *)filein->Get("tree");
        tree->SetBranchAddress("ptRec",&ptRec);
        tree->SetBranchAddress("mRec",&mRec);
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
        nentries=tree->GetEntries();

        // TODO based on the info in histTUnfold make map for systematics

	for(int i=0;i<nentries;i++){
        //for(int i=0;i<10000;i++){
          if(i%10000000==0) cout<<i<<endl;
          tree->GetEntry(i);

          int isdilep = 0;
          if( channel.CompareTo("electron") == 0 ) isdilep = IsElEl;
          if( channel.CompareTo("muon") == 0 ) isdilep = IsMuMu;

           if(isdilep && ispassRec && isBveto && ptRec->at(0) > 25 && ptRec->at(1) > 15){
                 fileout1->cd();

                 FillHistogram(ptOrMass::PT, "pt_nominal", ptRec->at(2), mRec->at(2), weightGen*weightRec*bTagReweight); // read name from histMaps
                 // a function to select weight regarding the histogram name
                 FillHistogram(ptOrMass::MASS, "mass_nominal", ptRec->at(2), mRec->at(2), weightGen*weightRec*bTagReweight);
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

        else if(sysName == "AlphaS"){
                if(isReco){
                        return weightRec*bTagReweight*AlphaS->at(nthSys);
                }   
                else{
                        return weightGen*AlphaS->at(nthSys);
                }   
        }// AlphaS
	else{
		return 1.;
        }

}

void histTUnfold::saveSigHists(TFile *filein, TFile *fileout1, TFile *fileout2, TString channel, Double_t temp_kfactor){ // TODO add list of systematics

        gROOT->SetBatch();
        gInterpreter->GenerateDictionary("vector<TLorentzVector>","TLorentzVector.h;vector");
        TH1::SetDefaultSumw2();

        ptPreFSR = 0;
        mPreFSR = 0;
        ptPostFSR = 0;
        mPostFSR = 0;
        ptRec = 0;
        etaRec = 0;
        mRec = 0;
        l1PreFire = 0;
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

        TFile* fpthist = new TFile("/home/jhkim/ISR2016/unfolding/systematic/hPtRecnominal.root");
        TH1* hptcorr=(TH1*)fpthist->Get("data_");
        hptcorr->SetDirectory(0);
        TUnfoldBinning* binning_ptRec=(TUnfoldBinning*)fpthist->Get("Rec_Pt");

        TFile* fmasshist = new TFile("/home/jhkim/ISR2016/unfolding/systematic/hMassRecnominal.root");
        TH1* hmasscorr=(TH1*)fpthist->Get("data_");
        hmasscorr->SetDirectory(0);
        TUnfoldBinning* binning_massRec=(TUnfoldBinning*)fmasshist->Get("Rec_Mass");

        // FIXME 
        TFile* fZptWeight = new TFile("/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/etc/ZptWeight/aMCNLO/electron/ZptWeight_electron.root", "r");
        //TFile* fZptWeight = new TFile("/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/etc/ZptWeight/MG/electron/ZptWeight_electron.root", "r");
        TH1 *hZptWeight;
        fZptWeight->GetObject("ZptWeight", hZptWeight);

        nentries=tree->GetEntries();

        // TODO based on the info in histTUnfold make map for systematics
        for(int i=0;i<nentries;i++){
        //for(int i=0;i<10000;i++){
          if(i%10000000==0) cout<<i<<endl;
          tree->GetEntry(i);

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

                if(isdilep && ispassRec && isBveto && ptRec->at(0) > 25 && ptRec->at(1) > 15){
                   if(DYtautau){
                      fileout2->cd();

                      FillHistogram(ptOrMass::PT, "pt_nominal_tau", ptRec->at(2), mRec->at(2), weightGen*weightRec*bTagReweight);
                      FillHistogram(ptOrMass::MASS, "mass_nominal_tau", ptRec->at(2), mRec->at(2), weightGen*weightRec*bTagReweight);

                   }
                   else{
                      fileout1->cd();

                      FillHistogram(ptOrMass::PT, "pt_nominal", ptRec->at(2), mRec->at(2), weightGen*weightRec*bTagReweight);
                      FillHistogram(ptOrMass::MASS, "mass_nominal", ptRec->at(2), mRec->at(2), weightGen*weightRec*bTagReweight);

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
                   diptgen = ptPreFSR->at(2);
                   dimassgen = mPreFSR->at(2);

                   TLorentzVector temp_particlePostFSR[19];
                   TLorentzVector temp_anparticlePostFSR[19];
                   Double_t drcut[19] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9};

                   // add fsr photon in dr < 0.1 for lepton
                   for(int idrcut = 0; idrcut < 19; idrcut++){
                      temp_particlePostFSR[idrcut] = *particlePostFSR;
                      for(int gsize = 0; gsize < particleFSR->size(); gsize++){

                              Double_t temp_dr = sqrt(pow(particlePostFSR->Phi()-particleFSR->at(gsize).Phi(),2)+pow(particlePostFSR->Eta()-particleFSR->at(gsize).Eta(),2));
                              if(temp_dr < drcut[idrcut]) temp_particlePostFSR[idrcut] += particleFSR->at(gsize);
                      }
                   }

                   for(int idrcut = 0; idrcut < 19; idrcut++){
                      temp_anparticlePostFSR[idrcut] = *anparticlePostFSR;
                      for(int gsize = 0; gsize < anparticleFSR->size(); gsize++){

                        Double_t temp_dr = sqrt(pow(anparticlePostFSR->Phi()-anparticleFSR->at(gsize).Phi(),2)+pow(anparticlePostFSR->Eta()-anparticleFSR->at(gsize).Eta(),2));
                        if(temp_dr < drcut[idrcut]) temp_anparticlePostFSR[idrcut] += anparticleFSR->at(gsize);
                      }
                   }

                   for(int idrcut = 0; idrcut < 19; idrcut++){
                      Double_t temp_dipt = (temp_particlePostFSR[idrcut] + temp_anparticlePostFSR[idrcut]).Pt();
                      Double_t temp_dimass = (temp_particlePostFSR[idrcut] + temp_anparticlePostFSR[idrcut]).M();
                        TString dr_;
                      if(idrcut<9) {
                        dr_.Form("%d", idrcut+1);
                        FillMigration2DM( ptOrMass::PT, selected, "pt_FSRDR0p"+dr_,  diptrec, dimassrec, temp_dipt, temp_dimass, weightRec*bTagReweight, weightGen);
                        FillMigration2DM( ptOrMass::MASS, selected, "mass_FSRDR0p"+dr_, diptrec, dimassrec, temp_dipt, temp_dimass, weightRec*bTagReweight, weightGen);
                        }
                        else{
                        dr_.Form("%d", (idrcut+1)%10);
                        FillMigration2DM( ptOrMass::PT, selected, "pt_FSRDR1p"+dr_,  diptrec, dimassrec, temp_dipt, temp_dimass, weightRec*bTagReweight, weightGen);
                        FillMigration2DM( ptOrMass::MASS, selected, "mass_FSRDR1p"+dr_, diptrec, dimassrec, temp_dipt, temp_dimass, weightRec*bTagReweight, weightGen);
                        }
                   }

                   FillMigration2DM( ptOrMass::PT, selected, "pt_nominal",  diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight, weightGen);
                   FillMigration2DM( ptOrMass::MASS, selected, "mass_nominal", diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight, weightGen);

                   // unfolding systematic
                   Double_t ptCorr = hptcorr->GetBinContent(binning_ptRec->GetGlobalBinNumber(diptrec, dimassrec));
                   Double_t massCorr = hmasscorr->GetBinContent(binning_massRec->GetGlobalBinNumber(dimassrec));
                   //std::cout << "massCorr: " << massCorr << " mass: " << dimassrec << " selected " << selected << std::endl;

                   FillMigration2DM( ptOrMass::PT, selected, "pt_unfoldsys_0",  diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight, weightGen*ZPtCor);
                   FillMigration2DM( ptOrMass::MASS, selected, "mass_unfoldsys_0", diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight, weightGen*ZPtCor);

                   // scale systematic
                   for(unsigned int i=0; i<Scale->size(); i++){
                        TString is;
                        is.Form("%d", i);

                        FillMigration2DM( ptOrMass::PT, selected, "pt_Scale_"+is,  diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight, weightGen*Scale->at(i));
                        FillMigration2DM( ptOrMass::MASS, selected, "mass_Scale_"+is, diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight, weightGen*Scale->at(i));
                   }

                   // PDF error
                   for(unsigned int i=0; i<PDFerror->size(); i++){
                        TString is;
                        is.Form("%d", i);

                        FillMigration2DM( ptOrMass::PT, selected, "pt_PDFerror_"+is,  diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight, weightGen*PDFerror->at(i));
                        FillMigration2DM( ptOrMass::MASS, selected, "mass_PDFerror_"+is, diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight, weightGen*PDFerror->at(i));
                   }

		   // loop to create systematic response matrix
		   std::map<TString, int>::iterator it = sysMaps.begin();
		   while(it != sysMaps.end()){
			if(it->second == 2){ // FIXME

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

               }// DY to ee or mumu events only
               //
               //////////////////////////////////////////////////// fill migration matrix ////////////////////////////////////////

           }// for DYtoLL MC case
        }// event loop
        fileout2->cd();
        GetPtBinningRec()->Write();
        GetMassBinningRec()->Write();

        fileout1->cd();
        GetPtBinningRec()->Write();
        GetMassBinningRec()->Write();
        GetPtBinningGen()->Write();
        GetMassBinningGen()->Write();

        delete binning_ptRec;
        delete binning_massRec;
        delete hptcorr;
        delete hmasscorr;
        delete fpthist;
        delete fmasshist;

        delete hZptWeight;
        delete fZptWeight;
        delete tree;
        //fileout->cd();

        // seems ptHistograms automatically written 
        // recoHist.ptHists.at(0)->Write();
        // delete outputFile;
}




void histTUnfold::SetPtBinningRec(){

 // FIXME save binning information in other place and read from there
 const int nmassbin_fine=5;
 double massbin_fine[nmassbin_fine+1]={50,60,80,100,200,350};

 // pt bins for reco
 //const int nptbin_fine=50;
 //double ptbin_fine[nptbin_fine+1];
 //for(int i = 0; i < nptbin_fine + 1; i++){
 //   ptbin_fine[i] = i*2;
 //}

 const int nptbin_fine=17;
 double ptbin_fine[nptbin_fine+1]={0., 2., 4., 6., 8., 10., 12., 14., 18., 22., 28., 35., 45., 55., 65., 75., 85., 100.};

 ptBinningRec= new TUnfoldBinning("Rec_Pt");
 ptBinningRec->AddAxis("pt",nptbin_fine,ptbin_fine,false,true);
 ptBinningRec->AddAxis("mass",nmassbin_fine,massbin_fine,true,true);
}

void histTUnfold::SetPtBinningGen(){

 // FIXME save binning information in other place and read from there
 const int nmassbin_wide=5;
 double massbin_wide[nmassbin_wide+1]={50,60,80,100,200,350};

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

 ptBinningGen=(new TUnfoldBinning("Gen_Pt"));
 ptBinningGen->AddAxis("pt",nptbin_wide,ptbin_wide,false,true);
 ptBinningGen->AddAxis("mass",nmassbin_wide,massbin_wide,true,true);

}

void histTUnfold::SetMassBinningRec(){

 // FIXME save binning information in other place and read from there
 //const int nbin_fine=58;
 const int nbin_fine=54;
 //double bin_fine[nbin_fine+1]={40,42.5,45,47.5,50,52.5,55,57.5,60,62.5,65,67.5,70,72.5,75,77.5,80,82.5,85,87.5,90,92.5,95,97.5,100,102.5,105,107.5,110,112.5,115,117.5,120,123,126,129.5,133,137,141,145.5,150,155,160,165.5,171,178,185,192.5,200,209,218,229,240,254,268,284,300,325,350};
 double bin_fine[nbin_fine+1]={50,52.5,55,57.5,60,62.5,65,67.5,70,72.5,75,77.5,80,82.5,85,87.5,90,92.5,95,97.5,100,102.5,105,107.5,110,112.5,115,117.5,120,123,126,129.5,133,137,141,145.5,150,155,160,165.5,171,178,185,192.5,200,209,218,229,240,254,268,284,300,325,350};

 //const int nbin_fine=300;
 //double bin_fine[nbin_fine+1];
 //for(int i = 0; i < nbin_fine + 1; i++){
 //    bin_fine[i] = 50. + (double)i;
 //}

 massBinningRec=(new TUnfoldBinning("Rec_Mass"));
 massBinningRec->AddAxis("reco mass",nbin_fine,bin_fine,false,false);
}

void histTUnfold::SetMassBinningGen(){

 // FIXME save binning information in other place and read from there
 //const int nbin_wide=29;
 const int nbin_wide=27;
 //double bin_wide[nbin_wide+1]={40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,126,133,141,150,160,171,185,200,218,240,268,300,350};
 double bin_wide[nbin_wide+1]={50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,126,133,141,150,160,171,185,200,218,240,268,300,350};

 //const int nbin_wide=150;
 //double bin_wide[nbin_wide+1];
 //for(int i = 0; i < nbin_wide + 1; i++){
 //    bin_wide[i] = 50. + (double)(i*2);
 //}

 massBinningGen=(new TUnfoldBinning("Gen_Mass"));
 massBinningGen->AddAxis("gen mass",nbin_wide,bin_wide,true,true);
}
