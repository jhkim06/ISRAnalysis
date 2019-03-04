#include <iostream>
#include <map>
#include <cmath>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TDOMParser.h>
#include <TXMLDocument.h>
#include "TUnfoldBinningXML.h"
#include "TUnfoldBinning.h"
#include "TUnfoldDensity.h"

#include "makeRecoPlots.h"

TUnfoldDensityV17* setTUnfoldDensity(TFile *filein, TString var, TString matrixName){

	TH2* hmcGenRec = (TH2*)filein->Get("hmc" + var + "GenRec" + matrixName);
	TUnfoldBinningV17* binning_Rec = NULL;
	TUnfoldBinningV17* binning_Gen = NULL;

 	if( var.CompareTo("Pt") == 0 ){
	  binning_Rec = ptBinning_rec();
	  binning_Gen = ptBinning_gen();
        }
        if( var.CompareTo("Mass") == 0 ){
	  binning_Rec = massBinning_rec();
	  binning_Gen = massBinning_gen();
        }

        TUnfoldDensityV17 *unfold;
        unfold = new TUnfoldDensityV17(hmcGenRec,
                                       TUnfoldV17::kHistMapOutputHoriz,
                                       TUnfold::kRegModeNone, // fixed to use no regularisation temporary
                                       TUnfoldV17::kEConstraintArea,
                                       TUnfoldDensityV17::kDensityModeBinWidth,
                                       binning_Gen,binning_Rec);

        // seems not allowed to delete the inputs of the TUnfoldDensity 
        return unfold;
 
}

void setInput(TUnfoldDensityV17* unfold, TString var, TFile *filein){
	
        TH1* hRec = (TH1*)filein->Get("h"+var+"Recnorminal");
        unfold->SetInput(hRec, 1); // TODO allow to selec option for bias
			           // NEED TO KNOW THE EFFECT OF THE BIAS FACTOR
}

void subBkgs(TUnfoldDensityV17* unfold, TString var, TFile *filein, TString name){

        TH1* hRec = (TH1*)filein->Get("h"+var+"Recnorminal");
        unfold->SubtractBackground(hRec, name); 
}

void doUnfold(TUnfoldDensityV17* unfold){

	unfold->DoUnfold(0); // TODO allow option for regulrization
}
