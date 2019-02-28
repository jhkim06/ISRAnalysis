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

TUnfoldDensityV17* setTUnfoldDensity(TFile *filein){

	TH2* hmcGenRec = (TH2*)filein->Get("hmcGenRecnorminal");
	TUnfoldBinningV17* binning_Rec = binning_rec();
	TUnfoldBinningV17* binning_Gen = binning_gen();

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

void setInput(TUnfoldDensityV17* unfold, TFile *filein){
	
        TH1* hRec = (TH1*)filein->Get("hRecnorminal");
        unfold->SetInput(hRec, 0); // TODO allow to selec option for bias
}

void doUnfold(TUnfoldDensityV17* unfold){

	unfold->DoUnfold(0); // TODO allow option for regulrization
}
