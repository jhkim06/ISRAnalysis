
 // variables in the input tree
 Double_t weightGen,weightRec, weightRecIdUp, weightRecIdDown;
 Double_t weightRecTriUp, weightRecTriDown;
 Double_t weightRecRecoUp, weightRecRecoDown;
 Double_t weightGenPileUp, weightGenPileDown;
 vector<Double_t> *weightGenScale = 0;
 vector<Double_t> *weightGenPdf = 0;

 int ispassRec,isfiducialPreFSR,DYtautau, isBveto;
 int nentries;

 vector<Double_t> *ptPreFSR = 0;
 vector<Double_t> *mPreFSR = 0;
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
 Int_t nVtx;

