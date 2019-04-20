#include "ISR_unfoldUtils.h"
#include "TUnfoldBinningXML.h"
#include "TUnfoldBinning.h"
#include "TUnfoldDensity.h"

TUnfoldDensityV17* setTUnfoldDensity(TFile *filein, TString var, TString matrixName, bool isfsr){

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


        TUnfoldDensityV17 *unfold;
        unfold = new TUnfoldDensityV17(hmcGenRec,
                                       TUnfold::kHistMapOutputHoriz,
                                       TUnfold::kRegModeNone, // fixed to use no regularisation temporary
                                       TUnfold::kEConstraintArea,
                                       TUnfoldDensityV17::kDensityModeBinWidth,
                                       binning_Gen,binning_Rec);

        // seems not allowed to delete the inputs of the TUnfoldDensity 
        return unfold;
 
}

void setVetorSystematic(TFile *filein, TUnfoldDensityV17* unfold, TString var, TString sysMatrixName, int size){

        for(int i = 0; i < size; i++){
           TString ithSys; ithSys.Form("%d",i);
           unfold->AddSysError((TH2*)filein->Get("hmc" + var + "GenRec" + sysMatrixName+"_"+ithSys), sysMatrixName+"_"+ithSys, TUnfold::kHistMapOutputHoriz, TUnfoldSys::kSysErrModeMatrix);
        }
}

void setInput(TUnfoldDensityV17* unfold, TString var, TString postfix, TFile *filein){
	
        TH1* hRec = (TH1*)filein->Get("h"+var+"Rec"+postfix);
        unfold->SetInput(hRec, 1); // TODO allow to select option for bias
			           // NEED TO KNOW THE EFFECT OF THE BIAS FACTOR
}

void setInputHist(TUnfoldDensityV17* unfold, TH1* hist){

        unfold->SetInput(hist, 1); // TODO allow to select option for bias
                                   // NEED TO KNOW THE EFFECT OF THE BIAS FACTOR
}

void subBkgs(TUnfoldDensityV17* unfold, TString var, TString hname, TFile *filein, TString name){

        TH1* hRec = (TH1*)filein->Get("h"+var+"Rec"+hname);
        unfold->SubtractBackground(hRec, name); 
}

void doUnfold(TUnfoldDensityV17* unfold){

	unfold->DoUnfold(0); // TODO allow option for regulrization
}

TH1* getHist(TUnfoldDensityV17* unfold, TString postfix){
	TH1* test = unfold->GetOutput("hdataUnfold_"+postfix,0,0,0,kFALSE);
        std::cout << "nbin: " << test->GetXaxis()->GetNbins() << std::endl;
	delete test;
	return unfold->GetOutput("hdataUnfold_"+postfix,0,0,0,kFALSE);
}
