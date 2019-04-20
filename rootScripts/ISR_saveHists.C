#include "ISR_saveHists.h"
//#include "tree.h"

void saveRecoHists(TFile *filein, TFile *fileout1, histTUnfold &recoHist, TString channel){ // TODO add list of systematics

        gROOT->SetBatch();
	gInterpreter->GenerateDictionary("vector<TLorentzVector>","TLorentzVector.h;vector");
        TH1::SetDefaultSumw2();

 	// variables in the input tree
 	Double_t weightGen,weightRec, weightRecIdUp, weightRecIdDown;
 	Double_t weightRecTriUp, weightRecTriDown;
 	Double_t weightRecRecoUp, weightRecRecoDown;
 	Double_t weightGenPileUp, weightGenPileDown;
 	vector<Double_t> *weightGenScale = 0;
 	vector<Double_t> *weightGenPdf = 0;

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
 	Int_t nVtx;

        Int_t IsElEl, IsMuMu;
        Double_t ZPtCor;
        Double_t bTagReweight;
        Int_t isdielectron, isdimuon;

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
        for(int i=0;i<nentries;i++){
        //for(int i=0;i<10000;i++){
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
	gInterpreter->GenerateDictionary("vector<TLorentzVector>","TLorentzVector.h;vector");
        TH1::SetDefaultSumw2();

        // variables in the input tree
        Double_t weightGen,weightRec, weightRecIdUp, weightRecIdDown;
        Double_t weightRecTriUp, weightRecTriDown;
        Double_t weightRecRecoUp, weightRecRecoDown;
        Double_t weightGenPileUp, weightGenPileDown;
        vector<Double_t> *weightGenScale = 0;
        vector<Double_t> *weightGenPdf = 0;

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
        Int_t nVtx;

        vector<TLorentzVector> *particleFSR;
        vector<TLorentzVector> *anparticleFSR;
        TLorentzVector *particlePostFSR;
        TLorentzVector *anparticlePostFSR;

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
        for(int i=0;i<nentries;i++){
        //for(int i=0;i<10000;i++){
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
                   Double_t diptgenDR0p1 = -999., dimassgenDR0p1 = -999.;

                   bool selected = (isdilep && ispassRec && isBveto && ptRec->at(0) > 25 && ptRec->at(1) > 15);
                   // check if there are selected dilepton at detector to avoid memory error
                   if(ptRec->size() == 3){ diptrec = ptRec->at(2), dimassrec = mRec->at(2);}
		   diptgen = ptPreFSR->at(2);
		   dimassgen = mPreFSR->at(2);

		   // add fsr photon in dr < 0.1 for lepton
		   for(int gsize = 0; gsize < particleFSR->size(); gsize++){
		      Double_t temp_dr = sqrt(pow(particlePostFSR->Phi()-particleFSR->at(gsize).Phi(),2)+pow(particlePostFSR->Eta()-particleFSR->at(gsize).Eta(),2));
		      if(temp_dr < 0.1) *particlePostFSR += particleFSR->at(gsize);
		   }

                   for(int gsize = 0; gsize < anparticleFSR->size(); gsize++){
                      Double_t temp_dr = sqrt(pow(anparticlePostFSR->Phi()-anparticleFSR->at(gsize).Phi(),2)+pow(anparticlePostFSR->Eta()-anparticleFSR->at(gsize).Eta(),2));
                      if(temp_dr < 0.1) *anparticlePostFSR += anparticleFSR->at(gsize);
                   }

		  diptgenDR0p1 = ((*particlePostFSR)+(*anparticlePostFSR)).Pt();
		  dimassgenDR0p1 = ((*particlePostFSR)+(*anparticlePostFSR)).M();


                   sigHist.FillMigration2DM( pt_unfold, selected, "pt_nominal",  diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight, weightGen);
                   sigHist.FillMigration2DM( mass_unfold, selected, "mass_nominal", diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight, weightGen);

                   sigHist.FillMigration2DM( pt_unfold, selected, "pt_FSRDR0p1",  diptrec, dimassrec, diptgenDR0p1, dimassgenDR0p1, weightRec*bTagReweight, weightGen);
                   sigHist.FillMigration2DM( mass_unfold, selected, "mass_FSRDR0p1", diptrec, dimassrec, diptgenDR0p1, dimassgenDR0p1, weightRec*bTagReweight, weightGen);

		   // unfolding systematic
		   Double_t ptCorr = hptcorr->GetBinContent(binning_ptRec->GetGlobalBinNumber(diptrec, dimassrec));
		   Double_t massCorr = hmasscorr->GetBinContent(binning_massRec->GetGlobalBinNumber(dimassrec));
		   //std::cout << "massCorr: " << massCorr << " mass: " << dimassrec << " selected " << selected << std::endl;

                   sigHist.FillMigration2DM( pt_unfold, selected, "pt_unfoldsys_0",  diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight, weightGen*ZPtCor);
                   sigHist.FillMigration2DM( mass_unfold, selected, "mass_unfoldsys_0", diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight, weightGen*ZPtCor);

                   // scale systematic
                   for(unsigned int i=0; i<Scale->size(); i++){
                        TString is;
                        is.Form("%d", i);

                        sigHist.FillMigration2DM( pt_unfold, selected, "pt_Scale_"+is,  diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight, weightGen*Scale->at(i));
                        sigHist.FillMigration2DM( mass_unfold, selected, "mass_Scale_"+is, diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight, weightGen*Scale->at(i));
                   }

                   // alpha s systematic
                   for(unsigned int i=0; i<AlphaS->size(); i++){
                        TString is;
                        is.Form("%d", i);

                        sigHist.FillMigration2DM( pt_unfold, selected, "pt_AlphaS_"+is,  diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight, weightGen*AlphaS->at(i));
                        sigHist.FillMigration2DM( mass_unfold, selected, "mass_AlphaS_"+is, diptrec, dimassrec, diptgen, dimassgen, weightRec*bTagReweight, weightGen*AlphaS->at(i));
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

