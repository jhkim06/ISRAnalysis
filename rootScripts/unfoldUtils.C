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

TUnfoldDensityV17* setTUnfoldDensity(TFile *filein, TString var, TString matrixName, bool isfsr){

	TH2* hmcGenRec = (TH2*)filein->Get("hmc" + var + "GenRec" + matrixName);
	TUnfoldBinningV17* binning_Rec = NULL;
	TUnfoldBinningV17* binning_Gen = NULL;

 	if( var.CompareTo("Pt") == 0 ){
	  if(!isfsr) binning_Rec = ptBinning_rec();
          else binning_Rec = ptBinning_gen();

	  binning_Gen = ptBinning_gen();
        }
        if( var.CompareTo("Mass") == 0 ){
	  if(!isfsr) binning_Rec = massBinning_rec();
	   else binning_Rec = massBinning_gen();

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
