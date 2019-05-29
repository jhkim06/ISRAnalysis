#include "ISR_histTUnfold.h"

void histTUnfold::CreateHistMap(const int which_unfold, TString hname, TString postfix){

        if(which_unfold == pt_unfold){
                 histMaps.insert(std::pair<TString, TH1*>("pt_"+hname+postfix, ptBinningRec->CreateHistogram("hPtRec"+hname)));
        }

        if(which_unfold == mass_unfold){
                 histMaps.insert(std::pair<TString, TH1*>("mass_"+hname+postfix, massBinningRec->CreateHistogram("hMassRec"+hname)));
        }
}

void histTUnfold::CreateHist2DMap(const int which_unfold, TString hname){

	if(which_unfold == pt_unfold){
		 hist2DMaps.insert(std::pair<TString, TH2*>("pt_"+hname, TUnfoldBinning::CreateHistogramOfMigrations(ptBinningGen, ptBinningRec,"hmcPtGenRec" + hname)));
	}

	if(which_unfold == mass_unfold){
		 hist2DMaps.insert(std::pair<TString, TH2*>("mass_"+hname, TUnfoldBinning::CreateHistogramOfMigrations(massBinningGen, massBinningRec,"hmcMassGenRec" + hname)));
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

void histTUnfold::FillMigration2DM(const int which_unfold, bool selected, TString hname, Double_t recoPt, Double_t RecoMass, Double_t truthPt, Double_t truthMass, Double_t wreco, Double_t wgen, Double_t corr){

	int binZero=0;

	if(selected){
		if( which_unfold == pt_unfold){
			(hist2DMaps.find(hname)->second)->Fill(ptBinningGen->GetGlobalBinNumber(truthPt, truthMass), ptBinningRec->GetGlobalBinNumber(recoPt, RecoMass), wgen*wreco*corr);
			(hist2DMaps.find(hname)->second)->Fill(ptBinningGen->GetGlobalBinNumber(truthPt, truthMass), binZero, wgen*(1.-wreco));
		}
                if( which_unfold == mass_unfold){
                        if(recoPt < 100.) (hist2DMaps.find(hname)->second)->Fill(massBinningGen->GetGlobalBinNumber(truthMass), massBinningRec->GetGlobalBinNumber(RecoMass), wgen*wreco*corr);
                        if(recoPt < 100.) (hist2DMaps.find(hname)->second)->Fill(massBinningGen->GetGlobalBinNumber(truthMass), binZero, wgen*(1.-wreco));
                }
	}
	else{
		if( which_unfold == pt_unfold) (hist2DMaps.find(hname)->second)->Fill(ptBinningGen->GetGlobalBinNumber(truthPt, truthMass), binZero, wgen);
		if( which_unfold == mass_unfold) (hist2DMaps.find(hname)->second)->Fill(massBinningGen->GetGlobalBinNumber(truthMass), binZero, wgen);
	}
}

void histTUnfold::saveRecoHists(TFile *filein, TFile *fileout1, TString channel){ // TODO add list of systematics

        gROOT->SetBatch();
        gInterpreter->GenerateDictionary("vector<TLorentzVector>","TLorentzVector.h;vector");
        TH1::SetDefaultSumw2();

        // variables in the input tree
        Double_t weightGen,weightRec ;

        Double_t        PUweight;
        Double_t        PUweight_Up;
        Double_t        PUweight_Dn;
        Double_t        trgSF;
        Double_t        trgSF_Up;
        Double_t        trgSF_Dn;
        Double_t        recoSF;
        Double_t        recoSF_Up;
        Double_t        recoSF_Dn;
        Double_t        IdSF;
        Double_t        IdSF_Up;
        Double_t        IdSF_Dn;
        Double_t        IsoSF;
        Double_t        IsoSF_Up;
        Double_t        IsoSF_Dn;

        int ispassRec,DYtautau, isBveto;
        int isfiducialPostFSR, isfiducialPreFSR;
        int nentries;

        Int_t IsElEl, IsMuMu;
        Double_t ZPtCor;
        Double_t bTagReweight;
        Int_t isdielectron, isdimuon;

        vector<Double_t> *ptRec = 0;
	vector<Double_t> *mRec = 0;

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

        treco->SetBranchAddress("PUweight", &PUweight);
        treco->SetBranchAddress("PUweight_Up", &PUweight_Up);
        treco->SetBranchAddress("PUweight_Dn", &PUweight_Dn);
        treco->SetBranchAddress("trgSF", &trgSF);
        treco->SetBranchAddress("trgSF_Up", &trgSF_Up);
        treco->SetBranchAddress("trgSF_Dn", &trgSF_Dn);
        treco->SetBranchAddress("recoSF", &recoSF);
        treco->SetBranchAddress("recoSF_Up", &recoSF_Up);
        treco->SetBranchAddress("recoSF_Dn", &recoSF_Dn);
        treco->SetBranchAddress("IdSF", &IdSF);
        treco->SetBranchAddress("IdSF_Up", &IdSF_Up);
        treco->SetBranchAddress("IdSF_Dn", &IdSF_Dn);
        treco->SetBranchAddress("IsoSF", &IsoSF);
        treco->SetBranchAddress("IsoSF_Up", &IsoSF_Up);
        treco->SetBranchAddress("IsoSF_Dn", &IsoSF_Dn);
        nentries=treco->GetEntries();

        // TODO based on the info in histTUnfold make map for systematics

        for(int i=0;i<10000;i++){
          if(i%10000000==0) cout<<i<<endl;
          treco->GetEntry(i);

          int isdilep = 0;
          if( channel.CompareTo("electron") == 0 ) isdilep = IsElEl;
          if( channel.CompareTo("muon") == 0 ) isdilep = IsMuMu;

           if(isdilep && ispassRec && isBveto && ptRec->at(0) > 25 && ptRec->at(1) > 15){
                 fileout1->cd();

                 FillHistogram(pt_unfold, "pt_nominal", ptRec->at(2), mRec->at(2), weightGen*weightRec*bTagReweight); // read name from histMaps
                 // a function to select weight regarding the histogram name
                 FillHistogram(mass_unfold, "mass_nominal", ptRec->at(2), mRec->at(2), weightGen*weightRec*bTagReweight);
           }

        }// event loop

        fileout1->cd();
        GetPtBinningRec()->Write();
        GetMassBinningRec()->Write();

        delete treco;
        //fileout->cd();
}

void histTUnfold::saveSigHists(TFile *filein, TFile *fileout1, TFile *fileout2, TString channel, Double_t temp_kfactor){ // TODO add list of systematics

        gROOT->SetBatch();
        gInterpreter->GenerateDictionary("vector<TLorentzVector>","TLorentzVector.h;vector");
        TH1::SetDefaultSumw2();

        // variables in the input tree
        Double_t weightGen,weightRec;

        int ispassRec,DYtautau, isBveto;
        int isfiducialPostFSR, isfiducialPreFSR;
        int nentries;

        vector<Double_t> *ptPreFSR = 0;
        vector<Double_t> *mPreFSR = 0;
        vector<Double_t> *ptPostFSR = 0;
        vector<Double_t> *mPostFSR = 0;
        vector<Int_t> *qLep = 0;
        vector<Double_t> *ptRec = 0;
        vector<Double_t> *etaRec = 0;
        vector<Double_t> *phiRec = 0;
        vector<Double_t> *mRec = 0;
        vector<Double_t> *TrigSF = 0;
        vector<Double_t> *Id1SF = 0;
        vector<Double_t> *Id2SF = 0;
        vector<Double_t> *Reco1SF = 0;
        vector<Double_t> *Reco2SF = 0;
        vector<Double_t> *l1PreFire = 0;
        vector<Double_t> *AlphaS = 0;
        vector<Double_t> *Scale = 0;
        vector<Double_t> *PDFerror = 0;
        Int_t nVtx;

        vector<TLorentzVector> *particleFSR;
        vector<TLorentzVector> *anparticleFSR;
        TLorentzVector *particlePostFSR;
        TLorentzVector *anparticlePostFSR;


        Double_t        PUweight;
        Double_t        PUweight_Up;
        Double_t        PUweight_Dn;
        Double_t        trgSF;
        Double_t        trgSF_Up;
        Double_t        trgSF_Dn;
        Double_t        recoSF;
        Double_t        recoSF_Up;
        Double_t        recoSF_Dn;
        Double_t        IdSF;
        Double_t        IdSF_Up;
        Double_t        IdSF_Dn;
        Double_t        IsoSF;
        Double_t        IsoSF_Up;
        Double_t        IsoSF_Dn;

        TBranch        *b_particleFSR;   //!
        TBranch        *b_anparticleFSR;   //!
        TBranch        *b_particlePostFSR;   //!
        TBranch        *b_anparticlePostFSR;   //!

        particleFSR = 0;
        anparticleFSR = 0;
        particlePostFSR = 0;
        anparticlePostFSR = 0;

        Int_t IsElEl, IsMuMu;
        Double_t ZPtCor;
        Double_t bTagReweight;
        Int_t isdielectron, isdimuon;

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
        tsignal->SetBranchAddress("particleFSR", &particleFSR, &b_particleFSR);
        tsignal->SetBranchAddress("anparticleFSR", &anparticleFSR, &b_anparticleFSR);
        tsignal->SetBranchAddress("particlePostFSR", &particlePostFSR, &b_particlePostFSR);
        tsignal->SetBranchAddress("anparticlePostFSR", &anparticlePostFSR, &b_anparticlePostFSR);

        tsignal->SetBranchAddress("AlphaS",&AlphaS);
        tsignal->SetBranchAddress("Scale",&Scale);
        tsignal->SetBranchAddress("PDFerror",&PDFerror);

        tsignal->SetBranchAddress("PUweight", &PUweight);
        tsignal->SetBranchAddress("PUweight_Up", &PUweight_Up);
        tsignal->SetBranchAddress("PUweight_Dn", &PUweight_Dn);
        tsignal->SetBranchAddress("trgSF", &trgSF);
        tsignal->SetBranchAddress("trgSF_Up", &trgSF_Up);
        tsignal->SetBranchAddress("trgSF_Dn", &trgSF_Dn);
        tsignal->SetBranchAddress("recoSF", &recoSF);
        tsignal->SetBranchAddress("recoSF_Up", &recoSF_Up);
        tsignal->SetBranchAddress("recoSF_Dn", &recoSF_Dn);
        tsignal->SetBranchAddress("IdSF", &IdSF);
        tsignal->SetBranchAddress("IdSF_Up", &IdSF_Up);
        tsignal->SetBranchAddress("IdSF_Dn", &IdSF_Dn);
        tsignal->SetBranchAddress("IsoSF", &IsoSF);
        tsignal->SetBranchAddress("IsoSF_Up", &IsoSF_Up);
        tsignal->SetBranchAddress("IsoSF_Dn", &IsoSF_Dn);

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

           if(!isInc){
              if(isdilep && ispassRec){
                    fileout1->cd();

                    FillHistogram(pt_unfold, "pt_nominal", ptRec->at(2), mRec->at(2), weightGen*weightRec*bTagReweight);
                    FillHistogram(mass_unfold, "mass_nominal", ptRec->at(2), mRec->at(2), weightGen*weightRec*bTagReweight);

              }
           }
           else{ // for inclusive DY MC

                Double_t zptWeight = 1.;
                zptWeight = ZPtCor;

                if(isdilep && ispassRec && isBveto && ptRec->at(0) > 25 && ptRec->at(1) > 15){
                   if(DYtautau){
                      fileout2->cd();

                      FillHistogram(pt_unfold, "pt_nominal_tau", ptRec->at(2), mRec->at(2), weightGen*weightRec*bTagReweight);
                      FillHistogram(mass_unfold, "mass_nominal_tau", ptRec->at(2), mRec->at(2), weightGen*weightRec*bTagReweight);

                   }
                   else{
                      fileout1->cd();

                      FillHistogram(pt_unfold, "pt_nominal", ptRec->at(2), mRec->at(2), weightGen*weightRec*bTagReweight);
                      FillHistogram(mass_unfold, "mass_nominal", ptRec->at(2), mRec->at(2), weightGen*weightRec*bTagReweight);

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

                   //TLorentzVector temp_particlePostFSR[19];
                   //TLorentzVector temp_anparticlePostFSR[19];
                   //Double_t drcut[19] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9};

                   //// add fsr photon in dr < 0.1 for lepton
                   //for(int idrcut = 0; idrcut < 19; idrcut++){
                   //   temp_particlePostFSR[idrcut] = *particlePostFSR;
                   //   for(int gsize = 0; gsize < particleFSR->size(); gsize++){

                   //           Double_t temp_dr = sqrt(pow(particlePostFSR->Phi()-particleFSR->at(gsize).Phi(),2)+pow(particlePostFSR->Eta()-particleFSR->at(gsize).Eta(),2));
                   //           if(temp_dr < drcut[idrcut]) temp_particlePostFSR[idrcut] += particleFSR->at(gsize);
                   //   }
                   //}

                   //for(int idrcut = 0; idrcut < 19; idrcut++){
                   //   temp_anparticlePostFSR[idrcut] = *anparticlePostFSR;
                   //   for(int gsize = 0; gsize < anparticleFSR->size(); gsize++){

                   //     Double_t temp_dr = sqrt(pow(anparticlePostFSR->Phi()-anparticleFSR->at(gsize).Phi(),2)+pow(anparticlePostFSR->Eta()-anparticleFSR->at(gsize).Eta(),2));
                   //     if(temp_dr < drcut[idrcut]) temp_anparticlePostFSR[idrcut] += anparticleFSR->at(gsize);
                   //   }
                   //}

                   //for(int idrcut = 0; idrcut < 19; idrcut++){
                   //   Double_t temp_dipt = (temp_particlePostFSR[idrcut] + temp_anparticlePostFSR[idrcut]).Pt();
                   //   Double_t temp_dimass = (temp_particlePostFSR[idrcut] + temp_anparticlePostFSR[idrcut]).M();
                   //     TString dr_;
                   //   if(idrcut<9) {
                   //     dr_.Form("%d", idrcut+1);
                   //     FillMigration2DM( pt_unfold, selected, "pt_FSRDR0p"+dr_,  diptrec, dimassrec, temp_dipt, temp_dimass, weightRec*bTagReweight, weightGen);
                   //     FillMigration2DM( mass_unfold, selected, "mass_FSRDR0p"+dr_, diptrec, dimassrec, temp_dipt, temp_dimass, weightRec*bTagReweight, weightGen);
                   //     }
                   //     else{
                   //     dr_.Form("%d", (idrcut+1)%10);
                   //     FillMigration2DM( pt_unfold, selected, "pt_FSRDR1p"+dr_,  diptrec, dimassrec, temp_dipt, temp_dimass, weightRec*bTagReweight, weightGen);
                   //     FillMigration2DM( mass_unfold, selected, "mass_FSRDR1p"+dr_, diptrec, dimassrec, temp_dipt, temp_dimass, weightRec*bTagReweight, weightGen);
                   //     }
                   //}

                   FillMigration2DM( pt_unfold, selected, "pt_nominal",  diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight, weightGen);
                   FillMigration2DM( mass_unfold, selected, "mass_nominal", diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight, weightGen);

                   // unfolding systematic
                   Double_t ptCorr = hptcorr->GetBinContent(binning_ptRec->GetGlobalBinNumber(diptrec, dimassrec));
                   Double_t massCorr = hmasscorr->GetBinContent(binning_massRec->GetGlobalBinNumber(dimassrec));
                   //std::cout << "massCorr: " << massCorr << " mass: " << dimassrec << " selected " << selected << std::endl;

                   FillMigration2DM( pt_unfold, selected, "pt_unfoldsys_0",  diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight, weightGen*ZPtCor);
                   FillMigration2DM( mass_unfold, selected, "mass_unfoldsys_0", diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight, weightGen*ZPtCor);

                   // scale systematic
                   for(unsigned int i=0; i<Scale->size(); i++){
                        TString is;
                        is.Form("%d", i);

                        FillMigration2DM( pt_unfold, selected, "pt_Scale_"+is,  diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight, weightGen*Scale->at(i));
                        FillMigration2DM( mass_unfold, selected, "mass_Scale_"+is, diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight, weightGen*Scale->at(i));
                   }

                   // PDF error
                   //for(unsigned int i=0; i<PDFerror->size(); i++){
                   //     TString is;
                   //     is.Form("%d", i);

                   //     FillMigration2DM( pt_unfold, selected, "pt_PDFerror_"+is,  diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight, weightGen*PDFerror->at(i));
                   //     FillMigration2DM( mass_unfold, selected, "mass_PDFerror_"+is, diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight, weightGen*PDFerror->at(i));
                   //}

                   // alpha s systematic
                   for(unsigned int i=0; i<AlphaS->size(); i++){
                        TString is;
                        is.Form("%d", i);

                        FillMigration2DM( pt_unfold, selected, "pt_AlphaS_"+is,  diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight, weightGen*AlphaS->at(i));
                        FillMigration2DM( mass_unfold, selected, "mass_AlphaS_"+is, diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight, weightGen*AlphaS->at(i));
                   }

                   // PU
                   for(unsigned int i=0; i<2; i++){
                        TString is;
                        is.Form("%d", i);

                        if(i==0){
                           FillMigration2DM( pt_unfold, selected, "pt_PU_"+is,  diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight*PUweight_Up/PUweight, weightGen);
                           FillMigration2DM( mass_unfold, selected, "mass_PU_"+is, diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight*PUweight_Up/PUweight, weightGen);
                        }
                        else{
                           FillMigration2DM( pt_unfold, selected, "pt_PU_"+is,  diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight*PUweight_Dn/PUweight, weightGen);
                           FillMigration2DM( mass_unfold, selected, "mass_PU_"+is, diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight*PUweight_Dn/PUweight, weightGen);
                        }
                   }
                   // trgSF
                   for(unsigned int i=0; i<2; i++){
                        TString is;
                        is.Form("%d", i);

                        if(i==0){
                           FillMigration2DM( pt_unfold, selected, "pt_trgSF_"+is,  diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight*trgSF_Up/trgSF, weightGen);
                           FillMigration2DM( mass_unfold, selected, "mass_trgSF_"+is, diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight*trgSF_Up/trgSF, weightGen);
                        }
                        else{
                           FillMigration2DM( pt_unfold, selected, "pt_trgSF_"+is,  diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight*trgSF_Dn/trgSF, weightGen);
                           FillMigration2DM( mass_unfold, selected, "mass_trgSF_"+is, diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight*trgSF_Dn/trgSF, weightGen);
                        }
                   }

                   // recoSF
                   for(unsigned int i=0; i<2; i++){
                        TString is;
                        is.Form("%d", i);

                        if(i==0){
                           FillMigration2DM( pt_unfold, selected, "pt_recoSF_"+is,  diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight*recoSF_Up/recoSF, weightGen);
                           FillMigration2DM( mass_unfold, selected, "mass_recoSF_"+is, diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight*recoSF_Up/recoSF, weightGen);
                        }
                        else{
                           FillMigration2DM( pt_unfold, selected, "pt_recoSF_"+is,  diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight*recoSF_Dn/recoSF, weightGen);
                           FillMigration2DM( mass_unfold, selected, "mass_recoSF_"+is, diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight*recoSF_Dn/recoSF, weightGen);
                        }
                   }

                   // IdSF
                   for(unsigned int i=0; i<2; i++){
                        TString is;
                        is.Form("%d", i);

                        if(i==0){
                           FillMigration2DM( pt_unfold, selected, "pt_IdSF_"+is,  diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight*IdSF_Up/IdSF, weightGen);
                           FillMigration2DM( mass_unfold, selected, "mass_IdSF_"+is, diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight*IdSF_Up/IdSF, weightGen);
                        }
                        else{
                           FillMigration2DM( pt_unfold, selected, "pt_IdSF_"+is,  diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight*IdSF_Dn/IdSF, weightGen);
                           FillMigration2DM( mass_unfold, selected, "mass_IdSF_"+is, diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight*IdSF_Dn/IdSF, weightGen);
                        }
                   }

                   // IsoSF
                   for(unsigned int i=0; i<2; i++){
                        TString is;
                        is.Form("%d", i);

                        if(i==0){
                           FillMigration2DM( pt_unfold, selected, "pt_IsoSF_"+is,  diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight*IsoSF_Up/IsoSF, weightGen);
                           FillMigration2DM( mass_unfold, selected, "mass_IsoSF_"+is, diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight*IsoSF_Up/IsoSF, weightGen);
                        }
                        else{
                           FillMigration2DM( pt_unfold, selected, "pt_IsoSF_"+is,  diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight*IsoSF_Dn/IsoSF, weightGen);
                           FillMigration2DM( mass_unfold, selected, "mass_IsoSF_"+is, diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight*IsoSF_Dn/IsoSF, weightGen);
                        }
                   }

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
        delete tsignal;
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
