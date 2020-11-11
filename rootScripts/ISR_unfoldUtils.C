#include "ISR_unfoldUtils.h"
#include "tdrstyle.C"
#include "CMS_lumi.C"

void ISRUnfold::setBias(double bias)
{
   nominal_bias = bias;
}

void ISRUnfold::setOutputBaseDir(TString outPath)
{
    output_baseDir = outPath;
}

const TVectorD& ISRUnfold::checkMatrixCond(TString var)
{
    TH2D* hProb = NULL;
    if(var == "Pt") hProb = (TH2D*) nomPtUnfold->GetProbabilityMatrix("hProb_"+var);
    if(var == "Mass") hProb = (TH2D*) nomMassUnfold->GetProbabilityMatrix("hProb_"+var);

    // TODO check hProb->GetNbinsY()+2,hProb->GetNbinsX()+2 without "F"
    TMatrixD* mProb=new TMatrixD(hProb->GetNbinsX()+2,hProb->GetNbinsY()+2,hProb->GetArray(),"F"); // +2 for under/overflow bins
    TMatrixD* mtProb=new TMatrixD(hProb->GetNbinsY()+2,hProb->GetNbinsX()+2);
    mtProb->Transpose(*mProb);
    TDecompSVD* svdProb=new TDecompSVD(*mtProb);
    //mtProb->Print();
    cout << "Decompose(), successed? " << svdProb->Decompose() << endl;
    const Int_t colLwb = svdProb->GetColLwb();
    const Int_t nCols  = svdProb->GetNcols();
    const TVectorD& singularValues = svdProb->GetSig();

    // Find minimum (>0)
    int i = colLwb+nCols-1;
    while(i >= colLwb)
    {
        if(singularValues[i] > 0)
        {
            double min = singularValues[i];
            cout << var << ", Cond(): " << singularValues[colLwb]/min << endl;
            break;
        }
        else
        {
           i--;
        }
    }
    return svdProb->GetSig();
}

double ISRUnfold::getSmearedChi2(TString var, TString filePath, TString dirName, TString steering, bool useAxis, bool divBinWidth)
{
    double chi2 = 0.;
    double ndf  = 0.;

    TH1* hData; // Data - Bkg
    TH1* hDY; // DY MC

    TString DYHistName_ = "histo_DYJetsToMuMu";
    if(channel_name == "electron") DYHistName_ = "histo_DYJetsToEE";
    if(var.Contains("Mass"))
    {
        hData = nomMassUnfold->GetInput("hData_"+var, 0, 0, steering, useAxis);
        hDY = getRawHist(var, filePath, dirName, DYHistName_, "Signal_"+var, steering, useAxis, divBinWidth); ;
    }
    else
    {
        hData = nomPtUnfold->GetInput("hData_"+var, 0, 0, steering, useAxis);
        hDY = getRawHist(var, filePath, dirName, DYHistName_, "Signal_"+var, steering, useAxis, divBinWidth); ;
    }

    for(int i=1;i<=hDY->GetNbinsX();i++)
    {
        ndf += 1.;
        if(hData->GetBinError(i) == 0)
        {
          std::cout << "unfolded " << i << " bin: " << hData->GetBinContent(i) << std::endl;
          std::cout << "error is zero in the " << i << " bin..." << std::endl;
          std::cout << "so skip this bin" << std::endl;
          continue;
        }
        // TODO add option to use input covariance matrix
        double pull=(hData->GetBinContent(i)-hDY->GetBinContent(i))/hData->GetBinError(i);
        //cout << "data: " << hData->GetBinContent(i) << " mc: " << hDY->GetBinContent(i) << " data error: " << hData->GetBinError(i) << endl;
        chi2+= pull*pull;
    }
    //cout << "chi^{2}, " << chi2 << endl;
    return chi2;
}

double ISRUnfold::getUnfoldedChi2(TString var, TString steering, bool useAxis, bool divBinWidth)
{
    divBinWidth = false;
    double chi2 = 0.;
    double ndf  = 0.;

    TH1* hData; // Data - Bkg
    TH1* hDY; // DY MC
    TH2* hRho;

    //TH1 *g_fcnHist=0;
    TMatrixD *g_fcnMatrix=0;

    if(var.Contains("Mass"))
    {
        hData = nomMassUnfold->GetOutput("hData_"+var, 0, 0, steering, useAxis);
        hDY = nomMassUnfold->GetBias("hData_"+var, 0, 0, steering, useAxis);
        hRho = nomMassUnfold->GetRhoIJtotal("histRho_chi_"+var, 0,0, steering, useAxis);
    }
    else
    {
        hData = nomPtUnfold->GetOutput("hData_"+var, 0, 0, steering, useAxis);
        hDY = nomPtUnfold->GetBias("hData_"+var, 0, 0, steering, useAxis);
        hRho = nomPtUnfold->GetRhoIJtotal("histRho_chi_"+var, 0,0, steering, useAxis);
    }

    // FIXME check if "n" or "n+1"
    int n = hData->GetNbinsX();
    TMatrixDSym v(n);
    for(int i=0;i<n;i++)
    {
       for(int j=0;j<n;j++)
        {
            v(i,j)=hRho->GetBinContent(i+1,j+1)*(hData->GetBinError(i+1)*hData->GetBinError(j+1));
        }
    }

    TMatrixDSymEigen ev(v);
    TMatrixD d(n,n);
    TVectorD di(ev.GetEigenValues());
    for(int i=0;i<n;i++) {
       if(di(i)>0.0) {
          d(i,i)=1./di(i);
       } else {
          cout<<"bad eigenvalue i="<<i<<" di="<<di(i)<<"\n";
          exit(0);
       }
    }

    TMatrixD O(ev.GetEigenVectors());
    TMatrixD DOT(d,TMatrixD::kMultTranspose,O);
    g_fcnMatrix=new TMatrixD(O,TMatrixD::kMult,DOT);
    TMatrixD test(*g_fcnMatrix,TMatrixD::kMult,v);
    int error=0;

    for(int i=0;i<n;i++)
    {
        if(TMath::Abs(test(i,i)-1.0)>1.E-7)
        {
            error++;
        }
        for(int j=0;j<n;j++)
        {
            if(i==j) continue;
            if(TMath::Abs(test(i,j)>1.E-7)) error++;
        }
    }


    // Calculate chi2
    //for(int i=0;i<hData->GetNbinsX();i++)
    //{
    //
    //    double di_=hData->GetBinContent(i+1)-hDY->GetBinContent(i+1);
    //    if(g_fcnMatrix)
    //    {
    //        for(int j=0;j<hData->GetNbinsX();j++)
    //        {
    //            double dj=hData->GetBinContent(j+1)-hDY->GetBinContent(j+1);
    //            chi2+=di_*dj*(*g_fcnMatrix)(i,j);
    //        }
    //    }
    //    else
    //    {
    //        double pull=di_/hData->GetBinError(i+1);
    //        chi2+=pull*pull;
    //    }
    //    ndf+=1.0;
    //}

    for(int i=1;i<=hDY->GetNbinsX();i++)
    {
        ndf += 1.;
        if(hData->GetBinError(i) == 0)
        {
          std::cout << "unfolded " << i << " bin: " << hData->GetBinContent(i) << std::endl;
          std::cout << "error is zero in the " << i << " bin..." << std::endl;
          std::cout << "so skip this bin" << std::endl;
          continue;
        }
        // TODO add option to use input covariance matrix
        double pull=(hData->GetBinContent(i)-hDY->GetBinContent(i))/hData->GetBinError(i);
        //cout << "data: " << hData->GetBinContent(i) << " mc: " << hDY->GetBinContent(i) << " data error: " << hData->GetBinError(i) << endl;
        chi2+= pull*pull;
    }
    //cout << "chi^{2}, " << chi2 << endl;
    return chi2;

}

void ISRUnfold::checkIterEMUnfold(void)
{
    setTDRStyle();
    double yMin=1.;
    double yLine=10.;
    double yMax=graph_SURE_IterativeSURE_pt->GetMaximum();

    gStyle->SetPadBottomMargin(0.2);
    TCanvas *canvas1=new TCanvas("compare","",3600,1200);
    canvas1->Divide(2,1);

    canvas1->cd(1);
    gPad->SetLogy();
    graph_SURE_IterativeSURE_pt->GetYaxis()->SetRangeUser(yMin,yMax);
    graph_SURE_IterativeSURE_pt->GetXaxis()->SetRangeUser(-1.5,100.5);
    graph_SURE_IterativeSURE_pt->GetXaxis()->SetTitle("iteration");
    graph_SURE_IterativeSURE_pt->GetXaxis()->SetTitleOffset(1.2);
    graph_SURE_IterativeSURE_pt->GetXaxis()->SetTitleFont(43);
    graph_SURE_IterativeSURE_pt->GetXaxis()->SetTitleSize(100);
    graph_SURE_IterativeSURE_pt->SetMarkerColor(kBlue);
    graph_SURE_IterativeSURE_pt->SetMarkerStyle(20);
    graph_SURE_IterativeSURE_pt->SetMarkerSize(2);
    graph_SURE_IterativeSURE_pt->DrawClone("APW");
    int n_scanSURE_iterative=graph_SURE_IterativeSURE_pt->GetN();
    double const *nIter_scanSURE_iterative=graph_SURE_IterativeSURE_pt->GetX();
    double const *DF_scanSURE_iterative=graph_DFdeviance_IterativeSURE_pt->GetX();
    double const *deviance_scanSURE=graph_DFdeviance_IterativeSURE_pt->GetY();
    TGraph *DF_iterative=new TGraph
       (n_scanSURE_iterative,nIter_scanSURE_iterative,DF_scanSURE_iterative);
    TGraph *deviance_iterative=new TGraph
       (n_scanSURE_iterative,nIter_scanSURE_iterative,deviance_scanSURE);
    DF_iterative->SetMarkerColor(kRed);
    DF_iterative->SetMarkerStyle(24);
    DF_iterative->SetMarkerSize(2);
    DF_iterative->DrawClone("P");
    deviance_iterative->SetMarkerColor(kMagenta);
    deviance_iterative->SetMarkerStyle(22);
    deviance_iterative->SetMarkerSize(2);
    deviance_iterative->DrawClone("P");
    TLine *line2=new TLine(iBest_pt,yLine,iBest_pt,1e4);
    line2->SetLineStyle(1);
    line2->Draw();
    TLegend *legend3=new TLegend(0.25,0.2,0.9,0.45,"Iterative EM, minimize SURE");
    legend3->SetBorderSize(0);
    legend3->SetFillStyle(0);
    legend3->SetTextSize(0.045);
    legend3->AddEntry(graph_SURE_IterativeSURE_pt,"SURE","p");
    legend3->AddEntry(DF_iterative,"D.F.","p");
    legend3->AddEntry(deviance_iterative,"deviance","p");

    legend3->AddEntry(line2,TString::Format
                      ("min(SURE) at iteration=%d",iBest_pt),"l");
    legend3->AddEntry((TObject *)0,TString::Format
                      ("D.F.=%3g",DF_scanSURE_iterative[iBest_pt]),"");
    legend3->Draw();

    // For mass unfolding
    yMax=graph_SURE_IterativeSURE_mass->GetMaximum();
    yMin=1e-2;
    canvas1->cd(2);
    gPad->SetLogy();
    graph_SURE_IterativeSURE_mass->GetYaxis()->SetRangeUser(yMin,yMax);
    graph_SURE_IterativeSURE_mass->GetXaxis()->SetRangeUser(-1.5,100.5);
    graph_SURE_IterativeSURE_mass->SetTitle(";iteration");
    graph_SURE_IterativeSURE_mass->GetXaxis()->SetTitleOffset(1.2);
    graph_SURE_IterativeSURE_mass->GetXaxis()->SetTitleFont(43);
    graph_SURE_IterativeSURE_mass->GetXaxis()->SetTitleSize(100);
    graph_SURE_IterativeSURE_mass->SetMarkerColor(kBlue);
    graph_SURE_IterativeSURE_mass->SetMarkerStyle(20);
    graph_SURE_IterativeSURE_mass->SetMarkerSize(2);
    graph_SURE_IterativeSURE_mass->DrawClone("APW");
    n_scanSURE_iterative=graph_SURE_IterativeSURE_mass->GetN();
    double const *nIter_scanSURE_iterative_mass=graph_SURE_IterativeSURE_mass->GetX();
    double const *DF_scanSURE_iterative_mass=graph_DFdeviance_IterativeSURE_mass->GetX();
    double const *deviance_scanSURE_mass=graph_DFdeviance_IterativeSURE_mass->GetY();
    TGraph *DF_iterative_mass=new TGraph
       (n_scanSURE_iterative,nIter_scanSURE_iterative_mass,DF_scanSURE_iterative_mass);
    TGraph *deviance_iterative_mass=new TGraph
       (n_scanSURE_iterative,nIter_scanSURE_iterative_mass,deviance_scanSURE_mass);
    DF_iterative_mass->SetMarkerColor(kRed);
    DF_iterative_mass->SetMarkerStyle(24);
    DF_iterative_mass->SetMarkerSize(2);
    DF_iterative_mass->DrawClone("P");
    deviance_iterative_mass->SetMarkerColor(kMagenta);
    deviance_iterative_mass->SetMarkerStyle(22);
    deviance_iterative_mass->SetMarkerSize(2);
    deviance_iterative_mass->DrawClone("P");
    TLine *line2_mass=new TLine(iBest_mass,yLine,iBest_mass,1e3);
    line2_mass->SetLineStyle(1);
    line2_mass->Draw();
    TLegend *legend3_mass=new TLegend(0.25,0.2,0.9,0.45,"Iterative EM, minimize SURE");
    legend3_mass->SetBorderSize(0);
    legend3_mass->SetFillStyle(0);
    legend3_mass->SetTextSize(0.045);
    legend3_mass->AddEntry(graph_SURE_IterativeSURE_mass,"SURE","p");
    legend3_mass->AddEntry(DF_iterative,"D.F.","p");
    legend3_mass->AddEntry(deviance_iterative,"deviance","p");

    legend3_mass->AddEntry(line2,TString::Format
                      ("min(SURE) at iteration=%d",iBest_mass),"l");
    legend3_mass->AddEntry((TObject *)0,TString::Format
                     ("D.F.=%3g",DF_scanSURE_iterative[iBest_mass]),"");
    legend3_mass->Draw();

    canvas1->SaveAs("ISR_scan.pdf");
}

void ISRUnfold::drawResponseM(TString var, TString sysName, TString sysPostfix, bool isDetector)
{
    TString yearStr;
    yearStr.Form("%d", (int)year);
    TFile f(output_baseDir+"Matrix_"+var+yearStr+".root","recreate");

    sysPostfix = "";
    sysName = "";
    const TVectorD* xaxis1_tvecd;
    int xaxis1_nbin;

    const TVectorD* yaxis1_tvecd = NULL;
    int yaxis1_nbin;

    TH2 *histProb = NULL;
    TH2 *histProb_woUO = NULL;
    bool draw_wo_UO = true;

    if(var=="Pt")
    {
        xaxis1_tvecd = pt_binning_Gen->GetDistributionBinning(0);
        xaxis1_nbin = xaxis1_tvecd->GetNrows() - 1; // number of bins without UO

        yaxis1_tvecd = pt_binning_Rec->GetDistributionBinning(0);
        yaxis1_nbin = yaxis1_tvecd->GetNrows() - 1;

        const TVectorD* massBinVector = pt_binning_Rec->GetDistributionBinning(1);
        const double* ptBinArrayRec = yaxis1_tvecd->GetMatrixArray();
        const double* ptBinArrayGen = xaxis1_tvecd->GetMatrixArray();
        //const double* massBinArray = massBinVector->GetMatrixArray();
        int nPtBinRec = yaxis1_tvecd->GetNrows() - 1;
        int nPtBinGen = xaxis1_tvecd->GetNrows() - 1;
        int nMassBin = massBinVector->GetNrows() - 1;

        vector<double> newPtBinVectorRec;
        //int nTotalPtBins = nPtBinRec * nMassBin;
        for(int iMassEdge = 0; iMassEdge < nMassBin; iMassEdge++)
        {
            for(int iPtEdge = 0; iPtEdge < nPtBinRec + 1; iPtEdge++)
            {
                if(iMassEdge == 0)
                {
                    newPtBinVectorRec.push_back(ptBinArrayRec[iPtEdge]);    
                }
                else
                {
                    if(iPtEdge == 0) continue;
                    double newPtEdge = iMassEdge * ptBinArrayRec[nPtBinRec] + ptBinArrayRec[iPtEdge];
                    newPtBinVectorRec.push_back(newPtEdge);
                }
            }
        } 

        vector<double> newPtBinVectorGen;
        //int nTotalPtBinsGen = nPtBinGen * nMassBin;
        for(int iMassEdge = 0; iMassEdge < nMassBin; iMassEdge++)
        {
            for(int iPtEdge = 0; iPtEdge < nPtBinGen + 1; iPtEdge++)
            {
                if(iMassEdge == 0)
                {
                    newPtBinVectorGen.push_back(ptBinArrayGen[iPtEdge]);    
                }
                else
                {
                    if(iPtEdge == 0) continue;
                    double newPtEdge = iMassEdge * ptBinArrayGen[nPtBinGen] + ptBinArrayGen[iPtEdge];
                    newPtBinVectorGen.push_back(newPtEdge);
                }
            }
        } 

        histProb = nomPtUnfold->GetProbabilityMatrix("Migration prob. for pt mass bin",";p_T(gen);p_T(Rec)");

        // TODO make a function
        if(draw_wo_UO)
        {
            histProb_woUO = new TH2D("responseM_woUO","responseM_woUO", xaxis1_nbin * nMassBin, &newPtBinVectorGen[0], yaxis1_nbin * nMassBin, &newPtBinVectorRec[0]);

            int iGenMassBin = 0;
            int iRecMassBin = 0;
            int nGenPtBinTotal = xaxis1_nbin;
            int nRecPtBinTotal = yaxis1_nbin;
            if(pt_binning_Gen->HasUnderflow(0)) nGenPtBinTotal++;
            if(pt_binning_Gen->HasOverflow(0)) nGenPtBinTotal++;
            if(pt_binning_Rec->HasUnderflow(0)) nRecPtBinTotal++;
            if(pt_binning_Rec->HasOverflow(0)) nRecPtBinTotal++;

            int ibinx_woUO = 1;
            for(int ibinx = 1; ibinx < histProb->GetNbinsX()+1; ibinx++)
            {
                iGenMassBin = ibinx / nGenPtBinTotal;
                if(ibinx % nGenPtBinTotal == 0) iGenMassBin = iGenMassBin - 1;

                if(pt_binning_Gen->HasUnderflow(1) && iGenMassBin == 0) continue;
                if(pt_binning_Gen->HasOverflow(1) && iGenMassBin == nMassBin) continue;

                if(pt_binning_Gen->HasUnderflow(0) && ibinx % nGenPtBinTotal == 1) continue; // first pt bin
                if(pt_binning_Gen->HasOverflow(0) && ibinx % nGenPtBinTotal == 0) continue; // last pt bin

                int ibiny_woUO = 1;
                for(int ibiny = 1; ibiny < histProb->GetNbinsY(); ibiny++)
                {

                    iRecMassBin = ibiny / nRecPtBinTotal;
                    if(ibiny % nRecPtBinTotal == 0) iRecMassBin = iRecMassBin - 1;

                    if(pt_binning_Rec->HasUnderflow(1) && iRecMassBin == 0) continue;
                    if(pt_binning_Rec->HasOverflow(1) && iRecMassBin == nMassBin) continue;

                    if(pt_binning_Rec->HasUnderflow(0) && ibiny % nRecPtBinTotal == 1) continue; // first pt bin
                    if(pt_binning_Rec->HasOverflow(0) && ibiny % nRecPtBinTotal == 0) continue; // last pt bin

                    //cout << "ibinx, ibiny " << ibinx << " " << ibiny << " nRecPtBinTotal: " << nRecPtBinTotal << endl;
                    histProb_woUO->SetBinContent(ibinx_woUO, ibiny_woUO, histProb->GetBinContent(ibinx, ibiny));
                    //cout << "x, y " << ibinx_woUO << " " << ibiny_woUO << " content: " << histProb->GetBinContent(ibinx, ibiny) << endl;
                    ibiny_woUO++;
                }
                ibinx_woUO++;
            }
        }
    }
    else
    {
        xaxis1_tvecd = mass_binning_Gen->GetDistributionBinning(0);
        xaxis1_nbin = xaxis1_tvecd->GetNrows() - 1;
        const Double_t* massGenBins = xaxis1_tvecd->GetMatrixArray();

        yaxis1_tvecd = mass_binning_Rec->GetDistributionBinning(0);
        yaxis1_nbin = yaxis1_tvecd->GetNrows() - 1;
        const Double_t* massRecBins = yaxis1_tvecd->GetMatrixArray();

        histProb = nomMassUnfold->GetProbabilityMatrix("Migration prob. for mass bin",";mass(gen);mass(Rec)", true);
        if(draw_wo_UO)
        {
            histProb_woUO = new TH2D("responseM_woUO","responseM_woUO", xaxis1_nbin, massGenBins, yaxis1_nbin, massRecBins);

            int iGenMassBin = 0;
            int iRecMassBin = 0;
            int nGenPtBinTotal = xaxis1_nbin;
            int nRecPtBinTotal = yaxis1_nbin;
            if(mass_binning_Gen->HasUnderflow(0)) nGenPtBinTotal++;
            if(mass_binning_Gen->HasOverflow(0)) nGenPtBinTotal++;
            if(mass_binning_Rec->HasUnderflow(0)) nRecPtBinTotal++;
            if(mass_binning_Rec->HasOverflow(0)) nRecPtBinTotal++;

            int ibinx_woUO = 1;
            for(int ibinx = 1; ibinx < histProb->GetNbinsX()+1; ibinx++)
            {
                iGenMassBin = ibinx / nGenPtBinTotal;
                if(ibinx % nGenPtBinTotal == 0) iGenMassBin = iGenMassBin - 1;

                if(mass_binning_Gen->HasUnderflow(1) && iGenMassBin == 0) continue;
                if(mass_binning_Gen->HasOverflow(1) && iGenMassBin == 1) continue;

                if(mass_binning_Gen->HasUnderflow(0) && ibinx % nGenPtBinTotal == 1) continue; // first mass bin
                if(mass_binning_Gen->HasOverflow(0) && ibinx % nGenPtBinTotal == 0) continue; // last mass bin

                int ibiny_woUO = 1;
                for(int ibiny = 1; ibiny < histProb->GetNbinsY(); ibiny++)
                {

                    iRecMassBin = ibiny / nRecPtBinTotal;
                    if(ibiny % nRecPtBinTotal == 0) iRecMassBin = iRecMassBin - 1;

                    if(mass_binning_Rec->HasUnderflow(1) && iRecMassBin == 0) continue;
                    if(mass_binning_Rec->HasOverflow(1) && iRecMassBin == 1) continue;

                    if(mass_binning_Rec->HasUnderflow(0) && ibiny % nRecPtBinTotal == 1) continue; // first mass bin
                    if(mass_binning_Rec->HasOverflow(0) && ibiny % nRecPtBinTotal == 0) continue; // last mass bin

                    //cout << "ibinx, ibiny " << ibinx << " " << ibiny << " nRecPtBinTotal: " << nRecPtBinTotal << endl;
                    histProb_woUO->SetBinContent(ibinx_woUO, ibiny_woUO, histProb->GetBinContent(ibinx, ibiny));
                    //cout << "x, y " << ibinx_woUO << " " << ibiny_woUO << " content: " << histProb->GetBinContent(ibinx, ibiny) << endl;
                    ibiny_woUO++;
                }
                ibinx_woUO++;
            }
        }

    }

    histProb_woUO->Write();   
    f.Close();
}

// Set the nominal response matrix
void ISRUnfold::setNominalRM(TString var, TString filepath, TString dirName, TString histName, TString binDef)
{
    //cout << "ISRUnfold::setNominalRM set response matrix..." << endl;
    TH1::AddDirectory(kFALSE);
    TFile* filein = new TFile(filepath);

    TString fullDirPath = dirName + "/" + var + "_ResMatrix_" + histName + binDef + "/";

    TString Rec_binName = "Rec_"+var;
    TString Gen_binName = "Gen_"+var;
    Rec_binName = fullDirPath + Rec_binName;
    Gen_binName = fullDirPath + Gen_binName;

    // Set bin definition
    if(var == "Pt")
    {
        pt_binning_Rec = (TUnfoldBinning*)filein->Get(Rec_binName);
        pt_binning_Gen = (TUnfoldBinning*)filein->Get(Gen_binName);
    }
    else if(var == "Mass")
    {
        mass_binning_Rec = (TUnfoldBinning*)filein->Get(Rec_binName);
        mass_binning_Gen = (TUnfoldBinning*)filein->Get(Gen_binName);
    }
    else
    {
        cout << "ISRUnfold::setNominalRM, only Pt and Mass available for var" << endl;
        exit (EXIT_FAILURE);
    }

    // Set mass bin edges
    setMassBindEdges();

    // Set response matrix
    // First, get the response matrix
    TH2* hmcGenRec;
    hmcGenRec = (TH2*)filein->Get(fullDirPath + "hmc" + var + "GenRec");

    if( var == "Pt" )
    {
    	nomPtUnfold = new TUnfoldDensity(hmcGenRec, TUnfold::kHistMapOutputHoriz, regMode, TUnfold::kEConstraintArea, TUnfoldDensity::kDensityModeNone, pt_binning_Gen, pt_binning_Rec);
        cout << "Check TUnfold version " << nomPtUnfold->GetTUnfoldVersion() << endl;
        hPtResponseM = (TH2*) hmcGenRec->Clone("hPtResponseM");

        // For statistical uncertainty
        if(makeStatUnfold)
        {
            // cout << "Create response matrix for statistical uncertainty..." << endl;
            for(int i = 0; i < statSize; i++)
            {
                statPtUnfold.push_back(new TUnfoldDensity(hmcGenRec, TUnfold::kHistMapOutputHoriz, regMode, TUnfold::kEConstraintArea, TUnfoldDensity::kDensityModeNone, pt_binning_Gen, pt_binning_Rec));
            }
        }
    }
    else
    {
        nomMassUnfold = new TUnfoldDensity(hmcGenRec, TUnfold::kHistMapOutputHoriz, TUnfold::kRegModeNone, TUnfold::kEConstraintArea, TUnfoldDensity::kDensityModeNone, mass_binning_Gen, mass_binning_Rec);
        hMassResponseM = (TH2*) hmcGenRec->Clone("hMassResponseM");

        // For statistical uncertainty
        if(makeStatUnfold)
        {
            for(int i = 0; i < statSize; i++)
            {
                statMassUnfold.push_back(new TUnfoldDensity(hmcGenRec, TUnfold::kHistMapOutputHoriz, TUnfold::kRegModeNone, TUnfold::kEConstraintArea, TUnfoldDensity::kDensityModeNone, mass_binning_Gen, mass_binning_Rec));
            }
        }
    }

    filein->Close();
    delete filein;
}

void ISRUnfold::setMassBindEdges()
{
    //cout << "ISRUnfold::setMassBindEdges " << endl;
    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    const Double_t* massBins = temp_tvecd->GetMatrixArray();
    int nMassBinEdges = temp_tvecd->GetNrows();

    if(massBinEdges.size() == 0)
    {
        for(int i = 0 ; i < nMassBinEdges; i++)
        {
            massBinEdges.push_back(massBins[i]);
            //cout << i << " th mass bin edge: " << massBins[i] << endl;
        }
    }
    else
    {
        cout << "ISRUnfold::setMassBindEdges massBinEdges already set." << endl;
    }

}

void ISRUnfold::setFromPrevUnfResult(ISRUnfold* unfold, bool useAccept)
{
    //cout << "setFromPrevUnfResult(), useAccept? " << useAccept << endl;
    // Loop over sytematics considered in the previous unfold class
    // So first get sysMap map object
    std::map<TString, std::vector<TString>> sysMap_previous = unfold->getSystematicMap();
    std::map<TString, std::vector<TString>>::iterator it = sysMap_previous.begin();
    while(it != sysMap_previous.end())
    {
        //cout << "Systematic name: " << it->first << endl;
        if(this->sysMap.find(it->first) == this->sysMap.end())
        {
            // Not found in this ISRUnfold class, but exits in the previous one
            int size = sysMap_previous[it->first].size();
            for(int ith = 0; ith < size; ith++)
            {
                // Create TUnfoldDensity using the DEFAULT response matrix
                //cout << "Systematic variation, " << sysMap_previous[it->first][ith] << endl;
                if((it->first).Contains("Unfolding") && !(sysMap_previous[it->first][ith]).Contains("Nominal"))
                {
                    this->iterEMPtUnfold   = new TUnfoldIterativeEM(hPtResponseM,TUnfoldDensity::kHistMapOutputHoriz,pt_binning_Gen,pt_binning_Rec);
                    this->iterEMMassUnfold = new TUnfoldIterativeEM(hMassResponseM,TUnfoldDensity::kHistMapOutputHoriz,mass_binning_Gen,mass_binning_Rec);

                    if(!useAccept)
                    {
                        this->iterEMPtUnfold->SetInput(unfold->iterEMPtUnfold->GetOutput("hUnfoldedPt_"+ it->first + "_" + sysMap_previous[it->first][ith],0,0,"*[*]",false), nominal_bias);
                        this->iterEMMassUnfold->SetInput(unfold->iterEMMassUnfold->GetOutput("hUnfoldedMass_"+ it->first + "_" + sysMap_previous[it->first][ith],0,0,"*[*]",false), nominal_bias);
                    }
                    else
                    {
                        //cout << "use acceptance corrected output!" << endl;
                        this->iterEMPtUnfold->SetInput(unfold->hSysFullPhasePtData[it->first][sysMap_previous[it->first][ith]], nominal_bias);
                        this->iterEMMassUnfold->SetInput(unfold->hSysFullPhaseMassData[it->first][sysMap_previous[it->first][ith]], nominal_bias);
                    }
                }
                else
                {
                    this->sysPtUnfold[it->first][sysMap_previous[it->first][ith]]   = new TUnfoldDensity(hPtResponseM,TUnfold::kHistMapOutputHoriz,regMode, TUnfold::kEConstraintArea, TUnfoldDensity::kDensityModeNone, pt_binning_Gen,pt_binning_Rec);
                    this->sysMassUnfold[it->first][sysMap_previous[it->first][ith]] = new TUnfoldDensity(hMassResponseM,TUnfold::kHistMapOutputHoriz,regMode, TUnfold::kEConstraintArea, TUnfoldDensity::kDensityModeNone,mass_binning_Gen, mass_binning_Rec);

                    if(!useAccept)
                    {
                        this->sysPtUnfold[it->first][sysMap_previous[it->first][ith]]->SetInput(unfold->sysPtUnfold[it->first][sysMap_previous[it->first][ith]]->GetOutput("hUnfoldedPt_"+ it->first + "_" + sysMap_previous[it->first][ith],0,0,"*[*]",false), nominal_bias);
                        this->sysMassUnfold[it->first][sysMap_previous[it->first][ith]]->SetInput(unfold->sysMassUnfold[it->first][sysMap_previous[it->first][ith]]->GetOutput("hUnfoldedMass_"+ it->first + "_" + sysMap_previous[it->first][ith],0,0,"*[*]",false), nominal_bias);
                    }
                    else
                    {
                        //cout << "use acceptance corrected output!" << endl;
                        this->sysPtUnfold[it->first][sysMap_previous[it->first][ith]]->SetInput(unfold->hSysFullPhasePtData[it->first][sysMap_previous[it->first][ith]], nominal_bias);
                        this->sysMassUnfold[it->first][sysMap_previous[it->first][ith]]->SetInput(unfold->hSysFullPhaseMassData[it->first][sysMap_previous[it->first][ith]], nominal_bias);
                    }
                }

                this->sysMap[it->first].push_back(sysMap_previous[it->first][ith]);
            }
        }
        else
        {
            // Found
            // Systematic both considered in this and previous unfolding
            // Loop over systematic varations
            int size = sysMap_previous[it->first].size();
            for(int ith = 0; ith < size; ith++)
            {
                if(!useAccept)
                {
                    this->sysPtUnfold[it->first][sysMap_previous[it->first][ith]]->SetInput(unfold->sysPtUnfold[it->first][sysMap_previous[it->first][ith]]->GetOutput("hUnfoldedPt_"+ it->first + "_" + sysMap_previous[it->first][ith],0,0,"*[*]",false), nominal_bias);
                    this->sysMassUnfold[it->first][sysMap_previous[it->first][ith]]->SetInput(unfold->sysMassUnfold[it->first][sysMap_previous[it->first][ith]]->GetOutput("hUnfoldedMass_"+ it->first + "_" + sysMap_previous[it->first][ith],0,0,"*[*]",false), nominal_bias);
                }
                else
                {

                    if((it->first).Contains("Unfolding") && !(sysMap_previous[it->first][ith]).Contains("Nominal"))
                    {
                        this->iterEMPtUnfold->SetInput(unfold->hSysFullPhasePtData[it->first][sysMap_previous[it->first][ith]], nominal_bias);
                        this->iterEMMassUnfold->SetInput(unfold->hSysFullPhaseMassData[it->first][sysMap_previous[it->first][ith]], nominal_bias);
                    }
                    else
                    {
                        this->sysPtUnfold[it->first][sysMap_previous[it->first][ith]]->SetInput(unfold->hSysFullPhasePtData[it->first][sysMap_previous[it->first][ith]], nominal_bias);
                        this->sysMassUnfold[it->first][sysMap_previous[it->first][ith]]->SetInput(unfold->hSysFullPhaseMassData[it->first][sysMap_previous[it->first][ith]], nominal_bias);
                    }
                }
            }
        }
        it++;
    }
    // Loop over variations of event selection efficiency correction
}

void ISRUnfold::setSystematicRM(TString var, TString filepath, TString dirName, TString histName, TString sysName, TString sysPostfix, TString histPostfix, TString binDef)
{
    TFile* filein = new TFile(filepath);
    TH2* hmcGenRec = NULL;

    TString fullHistName = "hmc" + var + "GenRec";
    if(histPostfix != "")
        fullHistName = fullHistName + "_" + sysPostfix;

    hmcGenRec = (TH2*)filein->Get(dirName + "/" + var + "_ResMatrix_" + histName + binDef + "/" + fullHistName);

    if( var == "Pt" )
    {
        if(sysName.Contains("Unfolding") && !sysPostfix.Contains("Nominal"))
        {
            iterEMPtUnfold = new TUnfoldIterativeEM(hmcGenRec,TUnfoldDensity::kHistMapOutputHoriz,pt_binning_Gen,pt_binning_Rec);
        }
        else
        {
            sysPtUnfold[sysName][sysPostfix] = new TUnfoldDensity(hmcGenRec, TUnfold::kHistMapOutputHoriz, regMode, TUnfold::kEConstraintArea, TUnfoldDensity::kDensityModeNone, pt_binning_Gen, pt_binning_Rec);
        }
    }
    else if( var == "Mass" )
    {
        if(sysName.Contains("Unfolding") && !sysPostfix.Contains("Nominal"))
        {
            iterEMMassUnfold = new TUnfoldIterativeEM(hmcGenRec,TUnfoldDensity::kHistMapOutputHoriz,mass_binning_Gen,mass_binning_Rec);
        }
        else
        {
            sysMassUnfold[sysName][sysPostfix] = new TUnfoldDensity(hmcGenRec, TUnfold::kHistMapOutputHoriz, TUnfold::kRegModeNone, TUnfold::kEConstraintArea, TUnfoldDensity::kDensityModeNone, mass_binning_Gen, mass_binning_Rec);
        }
    }
    else
    {
        cout << "ISRUnfold::setSystematicRM, only Pt and Mass available for var" << endl;
        exit (EXIT_FAILURE);
    }

    filein->Close();
    delete filein;
}

// Set input histogram from unfolding output
void ISRUnfold::setUnfInput(ISRUnfold* unfold, TString var, bool isSys, TString sysName, TString sysPostfix, bool useAccept)
{
    TH1::AddDirectory(kFALSE);

    if(!useAccept)
    {
        if(!isSys)
        {
            if(var=="Pt")
            {
                nomPtUnfold->SetInput(unfold->getUnfoldedHists("Pt", "UnfoldOut_Pt", "*[*]", false), 1.);
            }
            else
            {
                nomMassUnfold->SetInput(unfold->getUnfoldedHists("Mass", "UnfoldOut_Mass", "*[*]", false), 1.);
            }
        }
        else
        {
            if(var=="Pt")
            {
                sysPtUnfold[sysName][sysPostfix]->SetInput(unfold->getUnfoldedHists("Pt", "UnfoldOut_Pt"+sysName+sysPostfix, "*[*]", false), 1.);
            }
            else
            {
                sysMassUnfold[sysName][sysPostfix]->SetInput(unfold->getUnfoldedHists("Mass", "UnfoldOut_Mass"+sysName+sysPostfix, "*[*]", false), 1.);
            }
        }
    }
    else
    {
        //cout << "set from previous unfold class, isSys " << isSys << endl;
        if(!isSys)
        {
            if(var=="Pt")
            {
                nomPtUnfold->SetInput(unfold->hFullPhasePtData, 1.);
            }
            else
            {
                nomMassUnfold->SetInput(unfold->hFullPhaseMassData, 1.);
            }
        }
        else
        {
            if(var=="Pt")
            {
                if(sysName.Contains("Unfolding") && !sysPostfix.Contains("Nominal"))
                {
                    iterEMPtUnfold->SetInput(unfold->hFullPhasePtData, 1.);
                }
                else
                {
                    sysPtUnfold[sysName][sysPostfix]->SetInput(unfold->hFullPhasePtData, 1.);
                }
            }
            else
            {
                if(sysName.Contains("Unfolding") && !sysPostfix.Contains("Nominal"))
                {
                    iterEMMassUnfold->SetInput(unfold->hFullPhaseMassData, 1.);
                }
                else
                {
                    sysMassUnfold[sysName][sysPostfix]->SetInput(unfold->hFullPhaseMassData, 1.);
                }
            }
        }
    }
}

// Set input histogram from root file
void ISRUnfold::setUnfInput(TString var, TString varPostfix, TString filepath, TString dirName, TString histName, bool isSys, TString sysName, TString sysPostfix, bool isFSR)
{
    TH1::AddDirectory(kFALSE);

    TFile* filein = new TFile(filepath);
    TH1* hRec = NULL;
    //cout << dirName+"/"+var+varPostfix+"/"+histName << endl;
    hRec = (TH1*)filein->Get(dirName+"/"+var+varPostfix+"/"+histName);

    // Use DY MC as unfolding input, i.e. simple closure test
    if(!isFSR)
    {
        if(histName.Contains("DYJetsTo"))
        {
            histName.ReplaceAll("DYJetsTo", "DYJets10to50To");
            hRec->Add((TH1*)filein->Get(dirName+"/"+var+varPostfix+"/"+histName));
        }
    }
    else
    {
        if(histName.Contains("DYJets"))
        {
            histName.ReplaceAll("DYJets", "DYJets10to50");
            hRec->Add((TH1*)filein->Get(dirName+"/"+var+varPostfix+"/"+histName));
        }
    }

    // Very preliminary test for input covariance using ID SF
    //TFile* fcov = new TFile("/home/jhkim/ISR_Run2/unfolding/TUnfoldISR2016/rootScripts/covariance.root");
    //TFile* fcov_pt = new TFile("/home/jhkim/ISR_Run2/unfolding/TUnfoldISR2016/rootScripts/covariance_pt.root");
    //TH2* hCov = (TH2*) fcov->Get("cov");
    //TH2* hCov_pt = (TH2*) fcov_pt->Get("cov");

    // Nominal
    if(!isSys)
    {
        if(var == "Pt")
        {
            nomPtUnfold->SetInput(hRec,   nominal_bias, 0);
        }
        else if(var == "Mass")
        {
            nomMassUnfold->SetInput(hRec, nominal_bias, 0);
        }
        else{
            cout << "ISRUnfold::setUnfInput, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
        }
    }
    else
    // Systematic histograms
    {
        if(var == "Pt")
        {
            if(sysName.Contains("Unfolding") && !sysPostfix.Contains("Nominal"))
            {
                iterEMPtUnfold->SetInput(hRec, nominal_bias);
            }
            else
            {
                sysPtUnfold[sysName][sysPostfix]->SetInput(hRec, nominal_bias);
            }
        }
        else if(var == "Mass")
        {
            if(sysName.Contains("Unfolding") && !sysPostfix.Contains("Nominal"))
            {
                iterEMMassUnfold->SetInput(hRec, nominal_bias);
            }
            else
            {
                sysMassUnfold[sysName][sysPostfix]->SetInput(hRec, nominal_bias);
            }
        }
        else
        {
            cout << "ISRUnfold::setUnfInput, only Pt and Mass available for var" << endl;
            exit (EXIT_FAILURE);
        }
    }

    filein->Close();
    delete filein;
}

void ISRUnfold::subBkgs(TString filepath, std::pair<TString, TString>& bkgInfo, bool isSys, TString binDef, TString dirName, TString sysName, TString sysPostfix, TString histPostfix)
{
    TFile* filein = new TFile(filepath);
    TH1* hPtRec = NULL;
    TH1* hMassRec = NULL;

    // Nominal histograms
    if(!isSys)
    {
        bkgNames.push_back(bkgInfo.first);
        bkgTypes.push_back(bkgInfo.second);

        hPtRec = (TH1*)filein->Get(dirName + "/Pt"+binDef+"/histo_" + bkgInfo.first);
        nomPtUnfold->  SubtractBackground(hPtRec, bkgInfo.first);

        hMassRec = (TH1*)filein->Get(dirName + "/Mass"+binDef+"/histo_" + bkgInfo.first);
        nomMassUnfold->SubtractBackground(hMassRec, bkgInfo.first);
    }
    else
    // Systematic
    {
        TString fullHistName = bkgInfo.first + "_" + sysPostfix;
        if(histPostfix == "")
            fullHistName = bkgInfo.first;

        hPtRec = (TH1*)filein->Get(dirName + "/Pt"+binDef+"/histo_" + fullHistName);
        hMassRec = (TH1*)filein->Get(dirName + "/Mass"+binDef+"/histo_" + fullHistName);

        //cout << "file path: " << filepath << endl;
        //cout << dirName + "/Pt"+binDef+"/histo_" + fullHistName << endl;

        if(sysName.Contains("Unfolding") && !sysPostfix.Contains("Nominal"))
        {
            iterEMPtUnfold->SubtractBackground(hPtRec, bkgInfo.first);
            iterEMMassUnfold->SubtractBackground(hMassRec, bkgInfo.first);
        }
        else
        {
            // FIXME temporary method for background systematic
            if(sysName.Contains("Background"))
            {
                if(sysPostfix.Contains("Up"))
                {
                    sysPtUnfold[sysName][sysPostfix]->SubtractBackground(hPtRec, bkgInfo.first, 1.05);
                    sysMassUnfold[sysName][sysPostfix]->SubtractBackground(hMassRec, bkgInfo.first, 1.05);
                }
                if(sysPostfix.Contains("Down"))
                {
                    sysPtUnfold[sysName][sysPostfix]->SubtractBackground(hPtRec, bkgInfo.first, 0.95);
                    sysMassUnfold[sysName][sysPostfix]->SubtractBackground(hMassRec, bkgInfo.first, 0.95);
                }
            }
            else
            {
                sysPtUnfold[sysName][sysPostfix]->SubtractBackground(hPtRec, bkgInfo.first);
                sysMassUnfold[sysName][sysPostfix]->SubtractBackground(hMassRec, bkgInfo.first);
            }
        }
    }

    filein->Close();
    delete filein;
}

void ISRUnfold::setSystematics(TString sysName, TString sysHistName)
{
    sysMap[sysName].push_back(sysHistName);
}

TProfile* ISRUnfold::cloneHistToTProf(TH1* hist, TString histName)
{

  if (hist == NULL) return NULL;
  TProfile* cloneHist = new TProfile(histName, 
                             histName, 
                             hist->GetNbinsX(),
                             hist->GetXaxis()->GetXbins()->GetArray());
  return cloneHist;

}

TH1* ISRUnfold::cloneEmptyHist(TH1* hist, TString histName)
{

  if (hist == NULL) return NULL;
  TH1F* cloneHist = new TH1F(histName, 
                             histName, 
                             hist->GetNbinsX(),
                             hist->GetXaxis()->GetXbins()->GetArray());
  return cloneHist;

}

TGraphErrors* ISRUnfold::histToTGraphError(TH1* hist, bool zeroXerror)
{
    const int nPoints = hist->GetNbinsX();
    double x[nPoints];
    double y[nPoints];
    double xError[nPoints];
    double yError[nPoints];
    for(int i = 1; i < nPoints + 1; i++)
    {
        x[i-1] = hist->GetBinCenter(i);
        y[i-1] = hist->GetBinContent(i);
        if(zeroXerror)
            xError[i-1] = 0.0;
        else 
            xError[i-1] = hist->GetBinWidth(i)/2.; 
        yError[i-1] = hist->GetBinError(i);
    }
 
    TGraphErrors* out = new TGraphErrors(nPoints, x, y, xError, yError);
    return out;
}

void ISRUnfold::divideByBinWidth(TH1* hist, bool norm)
{
    for(int ibin = 1; ibin < hist->GetXaxis()->GetNbins()+1; ibin++)
    {
        double binWidth = hist->GetBinWidth(ibin);
        hist->SetBinContent(ibin, hist->GetBinContent(ibin)/ binWidth);
        hist->SetBinError(ibin, hist->GetBinError(ibin)/ binWidth);
    }
    if(norm)
    {
        hist->Scale(1./ hist->Integral());
    }
}

double ISRUnfold::getBinnedMean(TH1* hist)
{
    TH1* hBinned = (TH1*) hist->Clone("Binned");
    TH1* hDummy = (TH1*) hist->Clone("Dummy");
    hDummy->Reset("ICES");

    hBinned->Add(hDummy, -1); // After Add() Sumwx information removed, so can get binned mean.
    return hBinned->GetMean();
}

void ISRUnfold::doStatUnfold()
{
    //cout << "ISRUnfold::doStatUnfold() " << endl;
    for(int istat = 0; istat < statSize; istat++)
    {
        //cout << istat << " th stat.." << endl;
        TH1* tempMassInput;
        TH1* tempPtInput;

        TString nth_;
        nth_.Form("%d", istat);
        tempPtInput = nomPtUnfold->GetInput("tempPtHist_" + nth_, 0, 0, 0, false);
        tempMassInput = nomMassUnfold->GetInput("tempMassHist_" + nth_, 0, 0, 0, false);

        // randomize histogram bin content
        for(int ibin = 1; ibin<tempPtInput->GetNbinsX()+1;ibin++)
        {
            double err = tempPtInput->GetBinError(ibin);
            if(err > 0.0)
            {
                tempPtInput->SetBinContent(ibin, tempPtInput->GetBinContent(ibin) + gRandom->Gaus(0,err));
            }
        }
        for(int ibin = 1; ibin<tempMassInput->GetNbinsX()+1;ibin++)
        {
            double err = tempMassInput->GetBinError(ibin);
            if(err > 0.0)
            {
                tempMassInput->SetBinContent(ibin, tempMassInput->GetBinContent(ibin) + gRandom->Gaus(0,err));
            }
        }
        statPtUnfold.at(istat)->SetInput(tempPtInput, nominal_bias);
        statMassUnfold.at(istat)->SetInput(tempMassInput, nominal_bias);

        statPtUnfold.at(istat)->DoUnfold(0);
        statMassUnfold.at(istat)->DoUnfold(0);

        fillPtStatVariationHist(istat);
        fillMassStatVariationHist(istat);

        delete tempMassInput;
        delete tempPtInput;
        delete statPtUnfold.at(istat);
        delete statMassUnfold.at(istat);
    }
    statPtUnfold.clear();
    statMassUnfold.clear();
}

void ISRUnfold::doISRUnfold(bool doSys)
{

    TString yearStr;
    yearStr.Form("%d", (int)year);

    TFile* f;
    TDirectory* topDir;
    TDirectory* massDir;
    TDirectory* ptDir;

    TString fullPath=output_baseDir+unfold_name+"_"+channel_name+"_"+yearStr+".root";
    if(!gSystem->AccessPathName(fullPath, kFileExists))
    {
        f=TFile::Open(fullPath, "UPDATE"); 

        topDir=f->GetDirectory("unfolded");
        massDir=f->GetDirectory("unfolded/Mass");
        ptDir=f->GetDirectory("unfolded/Pt");
    }
    else
    {
        f=new TFile(fullPath, "CREATE");

        // Create directory
        topDir=f->mkdir("unfolded");
        massDir=topDir->mkdir("Mass");
        ptDir=topDir->mkdir("Pt"); 

        ptDir->cd();
        pt_binning_Gen->Write();
        massDir->cd();
        mass_binning_Gen->Write();
    }

    // cout << "ISRUnfold::doISRUnfold!!" << endl;
    // Nominal
    if(!doSys)
    {
        //cout << "Unfold without systematic" << endl;
        // No regularisation
        if(regMode == TUnfold::kRegModeNone)
        {
            // Nominal unfolding
            nomPtUnfold->DoUnfold(0);
            nomMassUnfold->DoUnfold(0);
        }
        else
        {
            nomMassUnfold->DoUnfold(0);

            int istart = pt_binning_Gen->GetGlobalBinNumber(0., 200.);
            int iend = pt_binning_Gen->GetGlobalBinNumber(99., 200.);
            nomPtUnfold->RegularizeBins(istart, 1, iend-istart+1, regMode);
            /*
            double tauMin=1.e-4;
            double tauMax=1.e-1;
            nomPtUnfold->ScanLcurve(100, tauMin, tauMax, 0);
            */
            TH2 *histL= nomPtUnfold->GetL("L");
            if(histL)
            {
                for(Int_t j=1;j<=histL->GetNbinsY();j++)
                {
                    cout<<"L["<<nomPtUnfold->GetLBinning()->GetBinName(j)<<"]";
                    for(Int_t i=1;i<=histL->GetNbinsX();i++) {
                        Double_t c=histL->GetBinContent(i,j);
                        if(c!=0.0) cout<<" ["<<i<<"]="<<c;
                    }
                    cout<<"\n";
                }
            }
        }

        ptDir->cd(); 
        nomPtUnfold->GetOutput("histo_Data",0,0, "*[*]", false)->Write(); 
        nomPtUnfold->GetBias("histo_DY", 0, 0, "*[*]", false)->Write(); 
        massDir->cd();
        nomMassUnfold->GetOutput("histo_Data",0,0, "*[*]", false)->Write(); 
        nomMassUnfold->GetBias("histo_DY", 0, 0, "*[*]", false)->Write(); 
    }
    // For systematic
    else
    {
        std::map<TString, std::vector<TString>>::iterator it = sysMap.begin();
        while(it != sysMap.end())
        {
            int size = (it->second).size();
            for(int i = 0; i < size; i++)
            {
                if( (it->first).Contains("Unfolding") && !((it->second).at(i)).Contains("Nominal"))
                {
                    iBest_pt=iterEMPtUnfold->ScanSURE(NITER_Iterative, &graph_SURE_IterativeSURE_pt, &graph_DFdeviance_IterativeSURE_pt);
                    iBest_mass=iterEMMassUnfold->ScanSURE(NITER_Iterative, &graph_SURE_IterativeSURE_mass, &graph_DFdeviance_IterativeSURE_mass);
                    cout << "iBest pt, Mass: " << iBest_pt << " " << iBest_mass << endl;

                    ptDir->cd();
                    iterEMPtUnfold->GetOutput("histo_Data_"+(it->second).at(i),0,0, "*[*]", false)->Write(); 
                    nomPtUnfold->GetBias("histo_DY_"+(it->second).at(i), 0, 0, "*[*]", false)->Write(); 

                    massDir->cd();
                    iterEMMassUnfold->GetOutput("histo_Data_"+(it->second).at(i),0,0, "*[*]", false)->Write(); 
                    nomMassUnfold->GetBias("histo_DY_"+(it->second).at(i), 0, 0, "*[*]", false)->Write(); 
                }
                else
                {

                    if(regMode == TUnfold::kRegModeNone)
                    {
                        sysPtUnfold[it->first][(it->second).at(i)]->DoUnfold(0);
                        sysMassUnfold[it->first][(it->second).at(i)]->DoUnfold(0);
                    }
                    else
                    {
                        double tauMin=1.e-4;
                        double tauMax=1.e-1;
                        sysPtUnfold[it->first][(it->second).at(i)]->ScanLcurve(100, tauMin, tauMax, 0);
                        sysMassUnfold[it->first][(it->second).at(i)]->DoUnfold(0);
                    }

                    if(it->first == "PDF")
                    {
                        fillMassPDFVariationHist(i+1); // i+1 since it starts from 1
                        fillPtPDFVariationHist(i+1);
                    }
            
                    ptDir->cd();
                    sysPtUnfold[it->first][(it->second).at(i)]->GetOutput("histo_Data_"+(it->second).at(i),0,0, "*[*]", false)->Write(); 
                    sysPtUnfold[it->first][(it->second).at(i)]->GetBias("histo_DY_"+(it->second).at(i), 0, 0, "*[*]", false)->Write(); 

                    massDir->cd();
                    sysMassUnfold[it->first][(it->second).at(i)]->GetOutput("histo_Data_"+(it->second).at(i),0,0, "*[*]", false)->Write(); 
                    sysMassUnfold[it->first][(it->second).at(i)]->GetBias("histo_DY_"+(it->second).at(i), 0, 0, "*[*]", false)->Write(); 
                }
            }
            it++;
        }
    }// Unfold for systematic
    topDir->Write();
    f->Close();
}

void ISRUnfold::drawCorrelation(TString var, TString steering, bool useAxis, TString outName)
{

    setTDRStyle();
    writeExtraText = true;
    extraText  = "Work in progress";
    gStyle->SetLineWidth(5);
    gStyle->SetFrameLineWidth(5);
    gStyle->SetOptStat(0);
    gROOT->ForceStyle();

    TCanvas* c = new TCanvas("c","c", 2400, 2400);
    c->SetTopMargin(0.1);
    c->SetRightMargin(0.15);
    c->cd();

    TH2* hCorrelation = NULL;
    TH2* hCorrelationUseAxis = NULL;

    if(var.Contains("Mass"))
    {
        hCorrelation=nomMassUnfold->GetRhoIJtotal("histRho", 0, 0, steering, useAxis);
    }
    else
    {

        hCorrelation=nomPtUnfold->GetRhoIJtotal("histRho", 0, 0, steering, useAxis);
        if(steering == "pt[UO];mass[UO]")
        {

            const TVectorD* xaxis1_tvecd = pt_binning_Gen->GetDistributionBinning(0);
            int xaxis1_nbin = xaxis1_tvecd->GetNrows() - 1; // number of bins without UO

            const TVectorD* massBinVector = pt_binning_Rec->GetDistributionBinning(1);
            const double* ptBinArrayGen = xaxis1_tvecd->GetMatrixArray();
            //const double* massBinArray = massBinVector->GetMatrixArray();
            int nPtBinGen = xaxis1_tvecd->GetNrows() - 1;
            int nMassBin = massBinVector->GetNrows() - 1;

            vector<double> newPtBinVectorGen;
            //int nTotalPtBinsGen = nPtBinGen * nMassBin;
            for(int iMassEdge = 0; iMassEdge < nMassBin; iMassEdge++)
            {
                for(int iPtEdge = 0; iPtEdge < nPtBinGen + 1; iPtEdge++)
                {
                    if(iMassEdge == 0)
                    {
                        newPtBinVectorGen.push_back(ptBinArrayGen[iPtEdge]);    
                    }
                    else
                    {
                        if(iPtEdge == 0) continue;
                        double newPtEdge = iMassEdge * ptBinArrayGen[nPtBinGen] + ptBinArrayGen[iPtEdge];
                        newPtBinVectorGen.push_back(newPtEdge);
                    }
                }
            } 

            hCorrelationUseAxis = new TH2D("correlation","correlation", xaxis1_nbin * nMassBin, &newPtBinVectorGen[0], xaxis1_nbin * nMassBin, &newPtBinVectorGen[0]);

            for(int ibinx = 1; ibinx < hCorrelation->GetNbinsX()+1; ibinx++)
            {
                for(int ibiny = 1; ibiny < hCorrelation->GetNbinsY()+1; ibiny++)
                {

                    hCorrelationUseAxis->SetBinContent(ibinx, ibiny, hCorrelation->GetBinContent(ibinx, ibiny));
                }
            }
            
        }
    }

    TString channel_name_;
    if(channel_name=="electron") channel_name_ = "ee";
    else channel_name_ = "#mu#mu";

    gStyle->SetPalette(kRainBow);

    if(hCorrelationUseAxis == NULL)
    {
        hCorrelation->SetMinimum(-1.);
        hCorrelation->SetMaximum(1.);
        hCorrelation->SetLineColor(kWhite);
        hCorrelation->Draw("COLZ");

        if(var.Contains("Pt"))
        {
            hCorrelation->GetYaxis()->SetTitle("p_{T}^{" + channel_name_ + "} [GeV]");
            hCorrelation->GetYaxis()->SetTitleOffset(1.0);
            hCorrelation->GetXaxis()->SetTitle("p_{T}^{" + channel_name_ + "} [GeV]");
            hCorrelation->GetXaxis()->SetTitleOffset(1.2);
        }
    }
    else
    {

        hCorrelationUseAxis->SetMinimum(-1.);
        hCorrelationUseAxis->SetMaximum(1.);
        hCorrelationUseAxis->SetLineColor(kWhite);
        hCorrelationUseAxis->Draw("COLZ");

        if(var.Contains("Pt"))
        {
            hCorrelationUseAxis->GetYaxis()->SetTitle("p_{T}^{" + channel_name_ + "} [GeV]");
            hCorrelationUseAxis->GetYaxis()->SetTitleOffset(1.0);
            hCorrelationUseAxis->GetXaxis()->SetTitle("p_{T}^{" + channel_name_ + "} [GeV]");
            hCorrelationUseAxis->GetXaxis()->SetTitleOffset(1.2);
        }

        TLine grid_;
        TLine grid_bin_boundary;
        grid_.SetLineColor(kBlack);
        grid_.SetLineStyle(1);

        const TVectorD* xaxis1_tvecd = pt_binning_Gen->GetDistributionBinning(0); 
        int nPtBinGen = xaxis1_tvecd->GetNrows() - 1;

        int boundarybin_x = 1;
        int countMassBin = 0;
        for( int ii=0; ii<hCorrelationUseAxis->GetXaxis()->GetNbins(); ii++ )
        {
            Int_t i_bin = ii+1;
            Double_t binEdge = hCorrelationUseAxis->GetXaxis()->GetBinLowEdge(i_bin);

            if(boundarybin_x == i_bin)
            {
                //grid_.DrawLine(binEdge, hCorrelationUseAxis->GetYaxis()->GetBinUpEdge(0), binEdge, hCorrelationUseAxis->GetYaxis()->GetBinUpEdge(hCorrelationUseAxis->GetYaxis()->GetNbins()) );
                boundarybin_x += nPtBinGen; // next edge to draw
                countMassBin++;
                if(i_bin == 1) continue;
                grid_.DrawLine(binEdge, hCorrelationUseAxis->GetYaxis()->GetBinUpEdge(nPtBinGen * (countMassBin - 2)), binEdge, hCorrelationUseAxis->GetYaxis()->GetBinUpEdge(nPtBinGen * (countMassBin)) );
            }
        }

        int boundarybin_y = 1;
        countMassBin = 0;
        for( int ii=0; ii<hCorrelationUseAxis->GetYaxis()->GetNbins(); ii++ )
        {
            Int_t i_bin = ii+1;
            Double_t binEdge = hCorrelationUseAxis->GetYaxis()->GetBinLowEdge(i_bin);

            if(boundarybin_y == i_bin)
            {
                boundarybin_y += nPtBinGen; // next edge to draw
                countMassBin++;
                if(i_bin == 1) continue;
                //grid_.DrawLine(hCorrelationUseAxis->GetXaxis()->GetBinUpEdge(0), binEdge, hCorrelationUseAxis->GetXaxis()->GetBinUpEdge(hCorrelationUseAxis->GetXaxis()->GetNbins()), binEdge);
                grid_.DrawLine(hCorrelationUseAxis->GetXaxis()->GetBinUpEdge(nPtBinGen * (countMassBin - 2)), binEdge, hCorrelationUseAxis->GetXaxis()->GetBinUpEdge(nPtBinGen * (countMassBin)), binEdge);
            }
        }

        if(var.Contains("Pt"))
        {
            hCorrelationUseAxis->GetXaxis()->SetNdivisions(506);
            hCorrelationUseAxis->GetYaxis()->SetNdivisions(506);

            hCorrelationUseAxis->GetXaxis()->ChangeLabel(2,-1,-1,-1,-1,-1,"100");
            hCorrelationUseAxis->GetXaxis()->ChangeLabel(3,-1,-1,-1,-1,-1,"100");
            hCorrelationUseAxis->GetXaxis()->ChangeLabel(4,-1,-1,-1,-1,-1,"100");
            hCorrelationUseAxis->GetXaxis()->ChangeLabel(5,-1,-1,-1,-1,-1,"100");
            hCorrelationUseAxis->GetXaxis()->ChangeLabel(6,-1,-1,-1,-1,-1,"100");

            hCorrelationUseAxis->GetYaxis()->ChangeLabel(2,-1,-1,-1,-1,-1,"100");
            hCorrelationUseAxis->GetYaxis()->ChangeLabel(3,-1,-1,-1,-1,-1,"100");
            hCorrelationUseAxis->GetYaxis()->ChangeLabel(4,-1,-1,-1,-1,-1,"100");
            hCorrelationUseAxis->GetYaxis()->ChangeLabel(5,-1,-1,-1,-1,-1,"100");
            hCorrelationUseAxis->GetYaxis()->ChangeLabel(6,-1,-1,-1,-1,-1,"100");

            TGaxis *xaxisMassLabel = new TGaxis(0,0,500,0,0,500,511,"+LS");
            xaxisMassLabel->SetTickLength(0);
            xaxisMassLabel->SetLabelSize(0);
            xaxisMassLabel->ChangeLabel(2,-1,0.015,-1,kGray+2,62,"50<M<64 GeV");
            xaxisMassLabel->ChangeLabel(4,-1,0.015,-1,kGray+2,62,"64<M<60 GeV");
            xaxisMassLabel->ChangeLabel(6,-1,0.015,-1,kGray+2,62,"81<M<101 GeV");
            xaxisMassLabel->ChangeLabel(8,-1,0.015,-1,kGray+2,62,"101<M<200 GeV");
            xaxisMassLabel->ChangeLabel(10,-1,0.015,-1,kGray+2,62,"200<M<320 GeV");
            xaxisMassLabel->SetLabelOffset(0.07);
            xaxisMassLabel->Draw();

            TGaxis *yaxisMassLabel = new TGaxis(0,0,0,500,0,500,511,"+LS");
            yaxisMassLabel->SetTickLength(0);
            yaxisMassLabel->SetLabelSize(0);
            yaxisMassLabel->ChangeLabel(1,90,0.015,-1,kGray+2,62,"50<M<64 GeV");
            yaxisMassLabel->ChangeLabel(3,90,0.015,-1,kGray+2,62,"64<M<60 GeV");
            yaxisMassLabel->ChangeLabel(5,90,0.015,-1,kGray+2,62,"81<M<101 GeV");
            yaxisMassLabel->ChangeLabel(7,90,0.015,-1,kGray+2,62,"101<M<200 GeV");
            yaxisMassLabel->ChangeLabel(9,90,0.015,-1,kGray+2,62,"200<M<320 GeV");
            yaxisMassLabel->SetLabelOffset(-0.07);
            yaxisMassLabel->CenterLabels();
            yaxisMassLabel->Draw();
        }

    }

    int iPeriod_ = 4;
    if(year == 2017)
        iPeriod_ = 5;
    if(year == 2018)
        iPeriod_ = 6;
    CMS_lumi(c, iPeriod_, 11);

    c->SaveAs(output_baseDir+"Correlation_"+var+"_"+outName+".pdf");
    delete hCorrelation;
    if(hCorrelationUseAxis != NULL)
        delete hCorrelationUseAxis;
}

void ISRUnfold::doAcceptCorr(TString filePath, TString binDef, bool doSys, TString outName, bool isAccept)
{
    TString yearStr;
    yearStr.Form("%d", (int)year);

    TFile* f; 
    TDirectory* topDir;
    TDirectory* massDir;
    TDirectory* ptDir;

    TString fullPath=output_baseDir+unfold_name+"_Acceptance_"+channel_name+"_"+yearStr+".root";
    if(!gSystem->AccessPathName(fullPath, kFileExists))
    {

        f=TFile::Open(fullPath, "UPDATE"); 

        topDir=f->GetDirectory("acceptance");
        massDir=f->GetDirectory("acceptance/Mass");
        ptDir=f->GetDirectory("acceptance/Pt");
    }
    else
    {

        f=new TFile(fullPath, "CREATE");

        // Create directory
        topDir=f->mkdir("acceptance");
        massDir=topDir->mkdir("Mass");
        ptDir=topDir->mkdir("Pt"); 

        ptDir->cd();
        pt_binning_Gen->Write();
        massDir->cd();
        mass_binning_Gen->Write();

    }

    TFile* filein = new TFile(filePath);

    TString accepCorrOrEffCorr;
    if(isAccept)
        accepCorrOrEffCorr = "Acceptance";
    else
        accepCorrOrEffCorr = "Efficiency";

    TH1* hFiducialPhaseMassMC = NULL;
    TH1* hFiducialPhasePtMC = NULL;

    // Nominal acceptance
    // Mass
    hFullPhaseMassMC = (TH1*) filein->Get("Acceptance/MassGen" + binDef + "/histo_DYJets");

    if(year==2016)
        hFullPhaseMassMC->Add((TH1*) filein->Get("Acceptance/MassGen" + binDef + "/histo_DYJets10to50"));
    else
        hFullPhaseMassMC->Add((TH1*) filein->Get("Acceptance/MassGen" + binDef + "/histo_DYJets10to50_MG"));

    hFiducialPhaseMassMC = nomMassUnfold->GetBias("hFiducialMass", 0, 0, "*[*]", false);
    hAcceptanceMass = (TH1*) hFullPhaseMassMC->Clone("hAcceptanceMass");
    hAcceptanceMass->Divide(hFiducialPhaseMassMC);

    hAcceptanceFractionMass = (TH1*) hFiducialPhaseMassMC->Clone("hAcceptanceFractionMass");
    hAcceptanceFractionMass->Divide(hFullPhaseMassMC);

    hFullPhaseMassData = nomMassUnfold->GetOutput("histo_Data",0,0, "*[*]", false);
    hFullPhaseMassData->Multiply(hAcceptanceMass); // Bin by bin acceptance correction 

    massDir->cd();
    hFullPhaseMassData->Write();
    hFullPhaseMassMC->SetName("histo_DY");
    hFullPhaseMassMC->Write();
    hAcceptanceMass->Write();

    // Pt
    hFullPhasePtMC = (TH1*) filein->Get("Acceptance/PtGen" + binDef + "/histo_DYJets");
    if(year==2016)
        hFullPhasePtMC->Add((TH1*) filein->Get("Acceptance/PtGen" + binDef + "/histo_DYJets10to50"));
    else
        hFullPhasePtMC->Add((TH1*) filein->Get("Acceptance/PtGen" + binDef + "/histo_DYJets10to50_MG"));

    hFiducialPhasePtMC = nomPtUnfold->GetBias("hFiducialPt", 0, 0, "*[*]", false);
    hAcceptancePt = (TH1*) hFullPhasePtMC->Clone("hAcceptancePt");

    hAcceptanceFractionPt = (TH1*) hFiducialPhasePtMC->Clone("hAcceptanceFractionPt");
    hAcceptanceFractionPt->Divide(hFullPhasePtMC);

    hAcceptancePt->Divide(hFiducialPhasePtMC);
    hFullPhasePtData = nomPtUnfold->GetOutput("histo_Data",0,0, "*[*]", false);
    hFullPhasePtData->Multiply(hAcceptancePt);

    ptDir->cd();
    hFullPhasePtData->Write();
    hFullPhasePtMC->SetName("histo_DY");
    hFullPhasePtMC->Write();
    hAcceptancePt->Write();

    if(doSys)
    {
        std::map<TString, vector<TString>> sysMapForAcceptance;

        std::map<TString, std::vector<TString>>::iterator it = sysMap.begin();
        while(it != sysMap.end())
        {
            int size = (it->second).size();
            for(int i = 0; i < size; i++)
            {
                TH1* hFullPhaseMassMC_raw_sys = NULL;
                TH1* hFullPhasePtMC_raw_sys = NULL;

                TH1* hFiducialPhaseMassMC_sys = NULL;
                TH1* hFiducialPhasePtMC_sys = NULL;
                
                if((it->first).Contains("Unfolding") && !((it->second).at(i)).Contains("Nominal"))
                {
                    hSysFullPhaseMassData[it->first][(it->second).at(i)] = iterEMMassUnfold->GetOutput("histo_Data_"+(it->second).at(i),0,0, "*[*]", false);
                    hSysFullPhasePtData[it->first][(it->second).at(i)]   = iterEMPtUnfold->GetOutput("histo_Data_"+(it->second).at(i),0,0, "*[*]", false);
                    hFiducialPhaseMassMC_sys=hFiducialPhaseMassMC;
                    hFiducialPhaseMassMC_sys->SetName("histo_DY_"+(it->second).at(i));
                    hFiducialPhasePtMC_sys=hFiducialPhasePtMC;
                    hFiducialPhasePtMC_sys->SetName("histo_DY_"+(it->second).at(i));
                }
                else
                {
                    hSysFullPhaseMassData[it->first][(it->second).at(i)] = sysMassUnfold[it->first][(it->second).at(i)]->GetOutput("histo_Data_"+(it->second).at(i),0,0, "*[*]", false);
                    hSysFullPhasePtData[it->first][(it->second).at(i)]   = sysPtUnfold[it->first][(it->second).at(i)]->GetOutput("histo_Data_"+(it->second).at(i),0,0, "*[*]", false);
                    hFiducialPhaseMassMC_sys = sysMassUnfold[it->first][(it->second).at(i)]->GetBias("histo_DY_"+(it->second).at(i), 0, 0, "*[*]", false);
                    hFiducialPhasePtMC_sys = sysPtUnfold[it->first][(it->second).at(i)]->GetBias("hFiducialPt_sys"+(it->second).at(i), 0, 0, "*[*]", false);
                }

                // For PDF, AlphaS, Scale etc, denominator changed
                if( (((it->first).Contains("Scale") && !(it->first).Contains("Lep")) || (it->first).Contains("PDF") || (it->first).Contains("AlphaS")) && !(it->first).Contains("_") )
                {
                    hFullPhaseMassMC_raw_sys = (TH1*) filein->Get("Acceptance/MassGen" + binDef + "/histo_DYJets_"+(it->second).at(i));
                    if(year==2016)
                        hFullPhaseMassMC_raw_sys->Add((TH1*) filein->Get("Acceptance/MassGen" + binDef + "/histo_DYJets10to50_"+(it->second).at(i)));
                    else
                        hFullPhaseMassMC_raw_sys->Add((TH1*) filein->Get("Acceptance/MassGen" + binDef + "/histo_DYJets10to50_MG_"+(it->second).at(i)));

                    hFullPhasePtMC_raw_sys = (TH1*) filein->Get("Acceptance/PtGen" + binDef + "/histo_DYJets_"+(it->second).at(i));
                    if(year==2016)
                        hFullPhasePtMC_raw_sys->Add((TH1*) filein->Get("Acceptance/PtGen" + binDef + "/histo_DYJets10to50_"+(it->second).at(i)));
                    else
                        hFullPhasePtMC_raw_sys->Add((TH1*) filein->Get("Acceptance/PtGen" + binDef + "/histo_DYJets10to50_MG_"+(it->second).at(i)));
                }
                else
                {
                    hFullPhaseMassMC_raw_sys=hFullPhaseMassMC;
                    hFullPhasePtMC_raw_sys=hFullPhasePtMC;
                }

                // For mass
                TH1* hAcceptanceMass_sys = (TH1*) hFullPhaseMassMC_raw_sys->Clone("hAcceptanceMass_sys");
                hAcceptanceMass_sys->Divide(hFiducialPhaseMassMC_sys);

                TH1* hAcceptanceFractionMass_sys = (TH1*) hFiducialPhaseMassMC_sys->Clone("hAcceptanceFractionMass_sys");
                hAcceptanceFractionMass_sys->Divide(hFullPhaseMassMC_raw_sys);
                
                //hSysFullPhaseMassData[accepCorrOrEffCorr + "_" + it->first][(it->second).at(i)] = nomMassUnfold->GetOutput("hAcceptMassData" +it->first+(it->second).at(i),0,0, "*[*]", false);
                hSysFullPhaseMassData[it->first][(it->second).at(i)]->Multiply(hAcceptanceMass_sys);
                hSysFullPhaseMassMC[it->first][(it->second).at(i)] = hFullPhaseMassMC_raw_sys;
                sysMapForAcceptance[it->first].push_back((it->second).at(i)); // Update sysMapForAcceptance 

                //hSysAcceptanceMass[it->first][(it->second).at(i)] = (TH1*) hAcceptanceMass_sys->Clone("Mass_" + it->first + "_" + (it->second).at(i));
                hSysAcceptanceFractionMass[it->first][(it->second).at(i)] = (TH1*) hAcceptanceFractionMass_sys->Clone("FractionMass_" + it->first + "_" + (it->second).at(i));
                delete hAcceptanceMass_sys;
                delete hAcceptanceFractionMass_sys;

                // For pt
                TH1* hAcceptancePt_sys = (TH1*) hFullPhasePtMC_raw_sys->Clone("hAcceptancePt_sys");
                hAcceptancePt_sys->Divide(hFiducialPhasePtMC_sys);

                TH1* hAcceptanceFractionPt_sys = (TH1*) hFiducialPhasePtMC_sys->Clone("hAcceptanceFractionPt_sys");
                hAcceptanceFractionPt_sys->Divide(hFullPhasePtMC_raw_sys);

                //hSysFullPhasePtData[accepCorrOrEffCorr + "_" + it->first][(it->second).at(i)] = nomPtUnfold->GetOutput("hAcceptPtData" +it->first+(it->second).at(i),0,0, "*[*]", false);
                hSysFullPhasePtData[it->first][(it->second).at(i)]->Multiply(hAcceptancePt_sys);
                hSysFullPhasePtMC[it->first][(it->second).at(i)] = hFullPhasePtMC_raw_sys;

                if((it->first).Contains("PDF"))
                {
                    fillMassPDFVariationHist_Accept(i+1);
                    fillPtPDFVariationHist_Accept(i+1);
                }

                //hSysAcceptancePt[it->first][(it->second).at(i)] = (TH1*) hAcceptancePt_sys->Clone("Pt_" + it->first + "_" + (it->second).at(i));
                hSysAcceptanceFractionPt[it->first][(it->second).at(i)] = (TH1*) hAcceptanceFractionPt_sys->Clone("FractionPt_" + it->first + "_" + (it->second).at(i));
                delete hAcceptancePt_sys;
                delete hAcceptanceFractionPt_sys;

                //hSysFullPhaseMassData[it->first][(it->second).at(i)]->Multiply(hAcceptanceMass);
                //hSysFullPhasePtData[it->first][(it->second).at(i)]->Multiply(hAcceptancePt);
                hSysFullPhaseMassMC[it->first][(it->second).at(i)] = hFullPhaseMassMC_raw_sys;
                hSysFullPhasePtMC[it->first][(it->second).at(i)]   = hFullPhasePtMC_raw_sys;

                ptDir->cd();
                hSysFullPhasePtData[it->first][(it->second).at(i)]->Write();
                hFullPhasePtMC_raw_sys->SetName("histo_DY_"+(it->second).at(i));
                hFullPhasePtMC_raw_sys->Write();

                massDir->cd();
                hSysFullPhaseMassData[it->first][(it->second).at(i)]->Write();
                hFullPhaseMassMC_raw_sys->SetName("histo_DY_"+(it->second).at(i));
                hFullPhaseMassMC_raw_sys->Write();

                hSysAcceptanceMass[it->first][(it->second).at(i)] = (TH1*) hAcceptanceMass->Clone("Mass_" + it->first + "_" + (it->second).at(i));
                hSysAcceptancePt[it->first][(it->second).at(i)] = (TH1*) hFullPhasePtMC->Clone("Pt_" + it->first + "_" + (it->second).at(i));
            }
            it++;
        }

        // Update sys map
        // Sys. unc. of acceptance correction
        // PDF, alphaS, Scale
        //it = sysMapForAcceptance.begin();
        //while(it != sysMapForAcceptance.end())
        //{
        //    int size = (it->second).size();
        //    for(int i = 0; i < size; i++)
        //    {
        //       sysMap[it->first].push_back((it->second).at(i)); 
        //    }
        //    it++;
        //}
  
        // Stat. unc. of acceptance correction

        hSysFullPhaseMassData[accepCorrOrEffCorr + "_Stat"]["Up"] = nomMassUnfold->GetOutput("hFullPhaseMassData"+accepCorrOrEffCorr+"StatUp",0,0, "*[*]", false); 
        hSysFullPhaseMassData[accepCorrOrEffCorr + "_Stat"]["Down"] = nomMassUnfold->GetOutput("hFullPhaseMassData"+accepCorrOrEffCorr+"StatDown",0,0, "*[*]", false); 

        TH1* hAcceptanceMass_statUp = (TH1*) hAcceptanceMass->Clone("h"+accepCorrOrEffCorr+"StatUpMass"); 
        TH1* hAcceptanceMass_statDown = (TH1*) hAcceptanceMass->Clone("h"+accepCorrOrEffCorr+"StatDownMass"); 
        varyHistWithStatError(hAcceptanceMass_statUp, 1);
        varyHistWithStatError(hAcceptanceMass_statDown, -1);

        hSysFullPhaseMassData[accepCorrOrEffCorr + "_Stat"]["Up"]->Multiply(hAcceptanceMass_statUp);
        hSysFullPhaseMassData[accepCorrOrEffCorr + "_Stat"]["Down"]->Multiply(hAcceptanceMass_statDown);
        hSysFullPhaseMassMC[accepCorrOrEffCorr + "_Stat"]["Up"] = hFullPhaseMassMC;
        hSysFullPhaseMassMC[accepCorrOrEffCorr + "_Stat"]["Down"] = hFullPhaseMassMC;

        hSysFullPhasePtData[accepCorrOrEffCorr + "_Stat"]["Up"] = nomPtUnfold->GetOutput("hFullPhasePtData"+accepCorrOrEffCorr+"StatUp",0,0, "*[*]", false); 
        hSysFullPhasePtData[accepCorrOrEffCorr + "_Stat"]["Down"] = nomPtUnfold->GetOutput("hFullPhasePtData"+accepCorrOrEffCorr+"StatDown",0,0, "*[*]", false); 

        TH1* hAcceptancePt_statUp = (TH1*) hAcceptancePt->Clone("h"+accepCorrOrEffCorr+"StatUpPt"); 
        TH1* hAcceptancePt_statDown = (TH1*) hAcceptancePt->Clone("h"+accepCorrOrEffCorr+"StatDownPt"); 
        varyHistWithStatError(hAcceptancePt_statUp, 1);
        varyHistWithStatError(hAcceptancePt_statDown, -1);

        hSysFullPhasePtData[accepCorrOrEffCorr + "_Stat"]["Up"]->Multiply(hAcceptancePt_statUp);
        hSysFullPhasePtData[accepCorrOrEffCorr + "_Stat"]["Down"]->Multiply(hAcceptancePt_statDown);
        hSysFullPhasePtMC[accepCorrOrEffCorr + "_Stat"]["Up"] = hFullPhasePtMC;
        hSysFullPhasePtMC[accepCorrOrEffCorr + "_Stat"]["Down"] = hFullPhasePtMC;

        //setSystematics(accepCorrOrEffCorr, "StatUp");
        //setSystematics(accepCorrOrEffCorr, "StatDown");
    }

    topDir->Write();;
    f->Close();

    delete hFiducialPhaseMassMC;
    delete hFiducialPhasePtMC;
}

void ISRUnfold::varyHistWithStatError(TH1* hist, int sys)
{
    for(int ibin = 1; ibin < hist->GetNbinsX()+1; ibin++)  
    {
        hist->SetBinContent(ibin, hist->GetBinContent(ibin) + double(sys) * hist->GetBinError(ibin));
    }
}

void ISRUnfold::setMeanMass_Accept()
{

    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;
    const Double_t* massBins = temp_tvecd->GetMatrixArray();

    TH1* hAccept_corr_mass = mass_binning_Gen->ExtractHistogram("hData_accept", hFullPhaseMassData, 0, kTRUE, "mass[UO];pt[UOC0]");
    for(int ibin = 0; ibin < nMassBin; ibin++)
    {
        hAccept_corr_mass->GetXaxis()->SetRange(hAccept_corr_mass->GetXaxis()->FindBin(massBins[ibin]+0.01), hAccept_corr_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));

        meanMass_data_acc_corrected. push_back(hAccept_corr_mass->GetMean());
        meanMassStatErr_data_acc_corrected.push_back(hAccept_corr_mass->GetMeanError());
    }

    delete hAccept_corr_mass;
}

void ISRUnfold::setMeanPt_Accept()
{

    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;

    for(int i = 0; i < nMassBin; i++)
    {
        TString ibinMass;
        ibinMass.Form("%d", i);

        TH1* hAccept_corr_pt = pt_binning_Gen->ExtractHistogram("hData_accept", hFullPhasePtData, 0, kTRUE, "pt[UO];mass[UOC"+ibinMass+"]");

        //cout << "Detector, " << i << " th mass bin, mean: " << hdetector_data->GetMean() << " +/- " << hdetector_data->GetMeanError() << endl;
        meanPt_data_acc_corrected.push_back(hAccept_corr_pt->GetMean());
        meanPtStatErr_data_acc_corrected.push_back(hAccept_corr_pt->GetMeanError());

        delete hAccept_corr_pt;
    }
}

void ISRUnfold::setTheoryMeanValues(TString filePath, TString binDef)
{

    TFile* filein = new TFile(filePath);
    //cout << "file: " << filePath << endl;
    TH1* hFullPhaseMassMCRaw = NULL;
    TH1* hFullPhasePtMCRaw = NULL;
   
    // Nominal mean value
    // Mass  
    hFullPhaseMassMCRaw = (TH1*) filein->Get("Acceptance/MassGen" + binDef + "/histo_DYJets");
    //cout << "Acceptance/MassGen" + binDef + "/histo_DYJets" << endl; 
    if(year==2016)
        hFullPhaseMassMCRaw->Add((TH1*) filein->Get("Acceptance/MassGen" + binDef + "/histo_DYJets10to50"));
    else
        hFullPhaseMassMCRaw->Add((TH1*) filein->Get("Acceptance/MassGen" + binDef + "/histo_DYJets10to50_MG"));

    // Pt
    hFullPhasePtMCRaw = (TH1*) filein->Get("Acceptance/PtGen" + binDef + "/histo_DYJets");
    if(year==2016)
        hFullPhasePtMCRaw->Add((TH1*) filein->Get("Acceptance/PtGen" + binDef + "/histo_DYJets10to50"));
    else
        hFullPhasePtMCRaw->Add((TH1*) filein->Get("Acceptance/PtGen" + binDef + "/histo_DYJets10to50_MG"));

    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;
    const Double_t* massBins = temp_tvecd->GetMatrixArray();

    // Nominal values
    // Extract histogram with axis definition
    TH1* hFullPhaseMassMC = mass_binning_Gen->ExtractHistogram("hMCMass", hFullPhaseMassMCRaw, 0, true, "mass[UO];pt[UOC0]");
    for(int ibin = 0; ibin < nMassBin; ibin++)
    {
        hFullPhaseMassMC->GetXaxis()->SetRange(hFullPhaseMassMC->GetXaxis()->FindBin(massBins[ibin]+0.01), hFullPhaseMassMC->GetXaxis()->FindBin(massBins[ibin+1]-0.01)); 
        meanMass_theory_acc_corrected. push_back(hFullPhaseMassMC->GetMean());  
        meanMassStatErr_theory_acc_corrected.push_back(hFullPhaseMassMC->GetMeanError()); 

        //cout << ibin << " th mass bin, mean mass " << hFullPhaseMassMC->GetMean() << " +/- " << hFullPhaseMassMC->GetMeanError() << endl;
    }
    delete hFullPhaseMassMC;

    for(int i = 0; i < nMassBin; i++)
    {
        TString ibinMass;
        ibinMass.Form("%d", i);

        TH1* hFullPhasePtMCTemp = pt_binning_Gen->ExtractHistogram("hMCPt", hFullPhasePtMCRaw, 0, kTRUE, "pt[UO];mass[UOC"+ibinMass+"]");

        //cout << "Detector, " << i << " th mass bin, mean: " << hdetector_data->GetMean() << " +/- " << hdetector_data->GetMeanError() << endl;
        meanPt_theory_acc_corrected.push_back(hFullPhasePtMCTemp->GetMean());
        meanPtStatErr_theory_acc_corrected.push_back(hFullPhasePtMCTemp->GetMeanError());

        //cout << i << " th mass bin, mean pt " << hFullPhasePtMCTemp->GetMean() << " +/- " << hFullPhasePtMCTemp->GetMeanError() << endl;

        delete hFullPhasePtMCTemp;
    }

    std::map<TString, std::vector<TString>>::iterator it = sysMap.begin();
    bool initTotalUnc = true;
    while(it != sysMap.end())
    {
        if(it->first == "Scale" || it->first == "AlphaS" || it->first == "PDF")
        {
            int size = (it->second).size();
            for(int i = 0; i < size; i++)
            {
                TH1* hFullPhaseMassMCRaw_ = NULL;
                TH1* hFullPhasePtMCRaw_ = NULL;

                hFullPhaseMassMCRaw_ = (TH1*) filein->Get("Acceptance/MassGen" + binDef + "/histo_DYJets_"+(it->second).at(i));
                //cout << "Acceptance/MassGen" + binDef + "/histo_DYJets" << endl; 
                if(year==2016)
                    hFullPhaseMassMCRaw_->Add((TH1*) filein->Get("Acceptance/MassGen" + binDef + "/histo_DYJets10to50_"+(it->second).at(i)));
                else
                    hFullPhaseMassMCRaw_->Add((TH1*) filein->Get("Acceptance/MassGen" + binDef + "/histo_DYJets10to50_MG_"+(it->second).at(i)));

                // Pt
                hFullPhasePtMCRaw_ = (TH1*) filein->Get("Acceptance/PtGen" + binDef + "/histo_DYJets_"+(it->second).at(i));
                if(year==2016)
                    hFullPhasePtMCRaw_->Add((TH1*) filein->Get("Acceptance/PtGen" + binDef + "/histo_DYJets10to50_"+(it->second).at(i)));
                else
                    hFullPhasePtMCRaw_->Add((TH1*) filein->Get("Acceptance/PtGen" + binDef + "/histo_DYJets10to50_MG_"+(it->second).at(i)));

                hFullPhaseMassMC = mass_binning_Gen->ExtractHistogram("hMCMass", hFullPhaseMassMCRaw_, 0, true, "mass[UO];pt[UOC0]");
                for(int ibin = 0; ibin < nMassBin; ibin++)
                {
                    hFullPhaseMassMC->GetXaxis()->SetRange(hFullPhaseMassMC->GetXaxis()->FindBin(massBins[ibin]+0.01), hFullPhaseMassMC->GetXaxis()->FindBin(massBins[ibin+1]-0.01)); 

                    TString ibinMass;
                    ibinMass.Form("%d", ibin);

                    if(it->first == "PDF")
                    {
                        if(i == 0)
                        {
                            meanMassPDFVariation_theory_Accept.push_back(new TH1F("MeanMassPDFTheory_bin"+ibinMass, "MeanMassPDFThory_bin"+ibinMass, 1000, 5., 500.));
                        }
                        meanMassPDFVariation_theory_Accept.at(ibin)->Fill(hFullPhaseMassMC->GetMean());
                    }
                    else
                    {
                        meanMass_theory_accept_sysVariation[it->first][(it->second).at(i)]. push_back(hFullPhaseMassMC->GetMean());  
                    }
                }
                delete hFullPhaseMassMC;

                for(int ibin = 0; ibin < nMassBin; ibin++)
                {
                    TString ibinMass;
                    ibinMass.Form("%d", ibin);

                    TH1* hFullPhasePtMC = pt_binning_Gen->ExtractHistogram("hMCPt", hFullPhasePtMCRaw_, 0, kTRUE, "pt[UO];mass[UOC"+ibinMass+"]");

                    //cout << "Detector, " << i << " th mass bin, mean: " << hdetector_data->GetMean() << " +/- " << hdetector_data->GetMeanError() << endl;
                    if(it->first == "PDF")
                    {
                        if(i == 0)
                        {
                            meanPtPDFVariation_theory_Accept.push_back(new TH1F("MeanPtPDFTheory_bin"+ibinMass, "MeanPtPDFThory_bin"+ibinMass, 1000, 5., 50.));
                        }
                        meanPtPDFVariation_theory_Accept.at(ibin)->Fill(hFullPhasePtMC->GetMean());
                    }
                    else
                    {
                        meanPt_theory_accept_sysVariation[it->first][(it->second).at(i)].push_back(hFullPhasePtMC->GetMean());
                    }
                    delete hFullPhasePtMC;
                }

                delete hFullPhaseMassMCRaw_;
                delete hFullPhasePtMCRaw_;
            }

            for(int ibin = 0; ibin < nMassBin; ibin++)
            {
                double error_mass = 0.;
                double error_pt = 0.;
                int size = (it->second).size();

                // Take the maximum variation as systematic uncertainty
                if(it->first != "PDF")
                {
                    for(int i = 0; i < size; i++)
                    {
                        double temp_error_mass = fabs(meanMass_theory_accept_sysVariation[it->first][(it->second).at(i)].at(ibin) - meanMass_theory_acc_corrected.at(ibin));
                        double temp_error_pt = fabs(meanPt_theory_accept_sysVariation[it->first][(it->second).at(i)].at(ibin) - meanPt_theory_acc_corrected.at(ibin));
                        error_mass = error_mass < temp_error_mass?temp_error_mass:error_mass;
                        error_pt = error_pt < temp_error_pt?temp_error_pt:error_pt;
                    }
                }
                else
                {
                    error_mass = meanMassPDFVariation_theory_Accept.at(ibin)->GetRMS();
                    error_pt = meanPtPDFVariation_theory_Accept.at(ibin)->GetRMS();
                }
                meanMass_theory_accept_systematic[it->first].push_back(error_mass);
                meanPt_theory_accept_systematic[it->first].push_back(error_pt);

                if(initTotalUnc)
                {
                    meanMass_theory_accept_systematic["total"].push_back(error_mass);
                    meanPt_theory_accept_systematic["total"].push_back(error_pt);
                }
                else
                {
                    meanMass_theory_accept_systematic["total"].at(ibin) = (sqrt(pow(error_mass,2) + pow(meanMass_theory_accept_systematic["total"].at(ibin),2)));
                    meanPt_theory_accept_systematic["total"].at(ibin) =  (sqrt(pow(error_pt,2) + pow(meanPt_theory_accept_systematic["total"].at(ibin),2)));
                }
                //cout << ibin << " th mass bin, " << it->first << " systematic unc. mass " << error_mass  << " pt " << error_pt << endl; 
                //cout << ibin << " th mass bin, " << " total systematic unc. mass " << meanMass_theory_accept_systematic["total"].at(ibin)  << " pt " << meanPt_theory_accept_systematic["total"].at(ibin) << endl; 
            }
            initTotalUnc = false;
            it++;
        }
        else
        {
            it++;
        }
    }
}

int ISRUnfold::setMeanMass(TString filePath)
{
    //cout << "ISRUnfold::setMeanMass()   Save mean of dilepton..." << endl;
    TString DYHistName_ = "histo_DYJetsToMuMu";
    if(channel_name == "electron")
    {
        DYHistName_  = "histo_DYJetsToEE";
    }

    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;
    const Double_t* massBins = temp_tvecd->GetMatrixArray();

    TUnfoldDensity* p_unfold = NULL;
    p_unfold = nomMassUnfold;

    TH1 * hdetector_mass = p_unfold->GetInput("hdetector_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);
    TH1* hunfolded_mass =  p_unfold->GetOutput("hunfolded_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);
    TH1* hMC_mass =  p_unfold->GetBias("hunfolded_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);
    TH1* hDY = NULL;
    if(filePath != "")
    {
        hDY = getRawHist("Mass_FineCoarse", filePath, "Detector", DYHistName_, "Signal", "mass[UO];pt[UOC0]", true, false);
    }

    // Loop over mass bins
    for(int ibin = 0; ibin < nMassBin; ibin++)
    {
        // Set x-axis range
        hdetector_mass->GetXaxis()->SetRange(hdetector_mass->GetXaxis()->FindBin(massBins[ibin]+0.01), hdetector_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
        hunfolded_mass->GetXaxis()->SetRange(hunfolded_mass->GetXaxis()->FindBin(massBins[ibin]+0.01), hunfolded_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
        hMC_mass->GetXaxis()->SetRange(hunfolded_mass->GetXaxis()->FindBin(massBins[ibin]+0.01), hunfolded_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
        if(hDY != NULL)
        {
            hDY->GetXaxis()->SetRange(hDY->GetXaxis()->FindBin(massBins[ibin]+0.01), hDY->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
            meanMass_mc_folded.   push_back(hDY->GetMean());
            meanMassStatErr_mc_folded.push_back(hDY->GetMeanError());
        }

        // Get mean values
        //cout << "Detector, " << ibin << " th mass bin, mean: " << hdetector_mass->GetMean() << " +/- " << hdetector_mass->GetMeanError() << endl;
        meanMass_data_folded. push_back(hdetector_mass->GetMean());
        meanMassStatErr_data_folded.push_back(hdetector_mass->GetMeanError());

        //cout << "Unfolded, " << ibin << " th mass bin, mean: " << hunfolded_mass->GetMean() << " +/- " << hunfolded_mass->GetMeanError() << endl;
        meanMass_data_unfolded.   push_back(hunfolded_mass->GetMean());
        meanMassStatErr_data_unfolded.push_back(hunfolded_mass->GetMeanError());

        //cout << "MC, " << ibin << " th mass bin, mean: " << hMC_mass->GetMean() << " +/- " << hMC_mass->GetMeanError() << endl;
        meanMass_mc_unfolded.   push_back(hMC_mass->GetMean());
        meanMassStatErr_mc_unfolded.push_back(hMC_mass->GetMeanError());

    }// End of mass bin loop

    delete hdetector_mass;
    delete hunfolded_mass;
    delete hMC_mass;
    if(hDY != NULL) delete hDY;

    return nMassBin;
}

int ISRUnfold::setSysMeanMass()
{
    //cout << "ISRUnfold::setSysMeanMass()   Save mean of dilepton..." << endl;

    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;
    const Double_t* massBins = temp_tvecd->GetMatrixArray();

    const TVectorD* temp_tvecd_rec = pt_binning_Rec->GetDistributionBinning(1);
    int nMassBin_rec = temp_tvecd_rec->GetNrows() - 1;
    const Double_t* massBins_rec = temp_tvecd_rec->GetMatrixArray();

    // For systematic
    std::map<TString, std::vector<TString>>::iterator it = sysMap.begin();
    while(it != sysMap.end())
    {
        TUnfoldDensity* p_unfold = NULL;

        //cout << "Unfold for " << it->first << " systematic." << endl;
        int size = (it->second).size();
        //cout << size << " systematic variation exist." << endl;

        for(int i = 0; i < size; i++)
        {
            TH1* hunfolded_mass = NULL;
            TH1* hunfolded_mass_input = NULL;
            if( (it->first).Contains("Unfolding") && !((it->second).at(i)).Contains("Nominal"))
            {
                hunfolded_mass = iterEMMassUnfold->GetOutput("hunfolded_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);
                hunfolded_mass_input = nomMassUnfold->GetInput("hunfolded_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);
            }
            else
            {
                if((it->first).Contains("Unfolding"))
                {
                    p_unfold = nomMassUnfold;;
                    hunfolded_mass = p_unfold->GetOutput("hunfolded_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);
                    hunfolded_mass_input = p_unfold->GetInput("hunfolded_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);
                }
                else
                {
                    p_unfold = sysMassUnfold[it->first][(it->second).at(i)];
                    hunfolded_mass = p_unfold->GetOutput("hunfolded_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);
                    hunfolded_mass_input = p_unfold->GetInput("hunfolded_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);
                }
            }


            // Loop over mass bins
            for(int ibin = 0; ibin < nMassBin; ibin++)
            {
                // Set x-axis range
                hunfolded_mass->GetXaxis()->SetRange(hunfolded_mass->GetXaxis()->FindBin(massBins[ibin]+0.01), hunfolded_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
                hunfolded_mass_input->GetXaxis()->SetRange(hunfolded_mass_input->GetXaxis()->FindBin(massBins_rec[ibin]+0.01), hunfolded_mass_input->GetXaxis()->FindBin(massBins_rec[ibin+1]-0.01));
                meanMass_data_unfolded_sysVariation[it->first][(it->second).at(i)].push_back(hunfolded_mass->GetMean());
                meanMass_data_folded_sysVariation[it->first][(it->second).at(i)].push_back(hunfolded_mass_input->GetMean());
                //cout << it->first << " " << (it->second).at(i) << " " << hunfolded_mass->GetMean() << endl;
            }// End of mass bin loop
            delete hunfolded_mass;
            delete hunfolded_mass_input;
        }
        it++;
    }
    return nMassBin;
}

void ISRUnfold::setSysMeanMass_Accept()
{
    //cout << "ISRUnfold::setSysMeanMass_Accept()   Save mean of dilepton..." << endl;
    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;
    const Double_t* massBins = temp_tvecd->GetMatrixArray();
    //cout << "nMassBin: " << nMassBin << endl;
    // For systematic
    std::map<TString, std::vector<TString>>::iterator it = sysMap.begin();
    while(it != sysMap.end())
    {
        int size = (it->second).size();

        for(int i = 0; i < size; i++)
        {

            TH1* hunfolded_mass_raw = hSysFullPhaseMassData[it->first][(it->second).at(i)];
            TH1* hunfolded_mass = mass_binning_Gen->ExtractHistogram("temp_mass", hunfolded_mass_raw, 0, kTRUE, "mass[UO];pt[UOC0]");
            // Loop over mass bins
            for(int ibin = 0; ibin < nMassBin; ibin++)
            {
                // Set x-axis range
                hunfolded_mass->GetXaxis()->SetRange(hunfolded_mass->GetXaxis()->FindBin(massBins[ibin]+0.01), hunfolded_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
                meanMass_data_accept_sysVariation[it->first][(it->second).at(i)].push_back(hunfolded_mass->GetMean());
                //cout << it->first << " " << (it->second).at(i) << " " << hunfolded_mass->GetMean() << " vector size: " << meanMass_data_accept_sysVariation[it->first][(it->second).at(i)].size() << endl;
            }// End of mass bin loop
            delete hunfolded_mass;
        }
        it++;
    }
}

void ISRUnfold::setSysError()
{
    cout << "------------------------------------- Systematic Uncertainty for ISR analysis ---------------------------------------" << endl;
    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;

    std::streamsize ss = cout.precision();
    cout.precision(2);
    cout.setf( std::ios::fixed, std:: ios::floatfield );

    // For systematic
    std::map<TString, std::vector<TString>>::iterator it = sysMap.begin();
    while(it != sysMap.end())
    {
        cout << "Systematic: " << it->first << endl;
        cout << setw(10) << "Mass bin" ;
        cout << setw(10) << "Mass" ;
        cout << setw(10) << "Pt" << endl;
        for(int ibin = 0; ibin < nMassBin; ibin++)
        {
            double error_mass = 0.;
            double error_pt = 0.;
            int size = (it->second).size();
            cout << setw(10) << ibin ;


            if(!(it->first).Contains("PDF"))
            {
                // Take the maximum variation as systematic uncertainty
                for(int i = 0; i < size; i++)
                {
                    double temp_error_mass = fabs(meanMass_data_unfolded_sysVariation[it->first][(it->second).at(i)].at(ibin) - meanMass_data_unfolded.at(ibin));
                    double temp_error_pt = fabs(meanPt_data_unfolded_sysVariation[it->first][(it->second).at(i)].at(ibin) - meanPt_data_unfolded.at(ibin));
                    error_mass = error_mass < temp_error_mass?temp_error_mass:error_mass;
                    error_pt = error_pt < temp_error_pt?temp_error_pt:error_pt;

                    if((it->first).Contains("FSR"))
                    {
                        error_mass = fabs(meanMass_data_unfolded_sysVariation[it->first][(it->second).at(0)].at(ibin)-meanMass_data_unfolded_sysVariation[it->first][(it->second).at(1)].at(ibin));
                        error_pt = fabs(meanPt_data_unfolded_sysVariation[it->first][(it->second).at(0)].at(ibin)-meanPt_data_unfolded_sysVariation[it->first][(it->second).at(1)].at(ibin));
                        break;
                    }
                }
            }
            else
            {
                //error_mass = meanMassPDFVariation.at(ibin)->GetRMS();
                //error_pt = meanPtPDFVariation.at(ibin)->GetRMS();

                TH1* htemp_mass = new TH1D("pdf_temp_mass", "pdf_temp_mass", 5000, 5, 500);
                TH1* htemp_pt = new TH1D("pdf_temp_pt", "pdf_temp_pt", 5000, 5, 500);
                for(int i = 0; i < size; i++)
                {
                    htemp_mass->Fill(meanMass_data_unfolded_sysVariation[it->first][(it->second).at(i)].at(ibin));
                    htemp_pt->Fill(meanPt_data_unfolded_sysVariation[it->first][(it->second).at(i)].at(ibin));
                }
                error_mass = htemp_mass->GetRMS();
                error_pt = htemp_pt->GetRMS();

                delete htemp_mass;
                delete htemp_pt;
            }

            meanMass_data_folded_systematic[it->first].push_back(error_mass);
            meanPt_data_folded_systematic[it->first].push_back(error_pt);

            meanMass_data_folded_rel_systematic[it->first].push_back(error_mass/ meanMass_data_unfolded.at(ibin) * 100.);
            meanPt_data_folded_rel_systematic[it->first].push_back(error_pt/ meanPt_data_unfolded.at(ibin) * 100.);

            cout << setw(10) << error_mass ;
            cout << setw(10) << error_pt << endl;
        }
        cout << endl;
        it++;
    }
    cout.precision(ss);
}

void ISRUnfold::setSysError_Accept()
{
    cout << "------------------------------------- Systematic Uncertainty for ISR analysis (Acceptance corrected)---------------------------------------" << endl;
    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;

    std::streamsize ss = cout.precision();
    cout.precision(2);
    cout.setf( std::ios::fixed, std:: ios::floatfield );

    // For systematic
    std::map<TString, std::vector<TString>>::iterator it = sysMap.begin();
    while(it != sysMap.end())
    {
        cout << "Systematic: " << it->first << endl;
        cout << setw(10) << "Mass bin" ;
        cout << setw(10) << "Mass" ;
        cout << setw(10) << "Pt" << endl;
        //cout << "Systematic: " << it->first << endl;
        for(int ibin = 0; ibin < nMassBin; ibin++)
        {
            double error_mass = 0.;
            double error_pt = 0.;
            int size = (it->second).size();
            cout << setw(10) << ibin ;

            // Take the maximum variation as systematic uncertainty
            if(!(it->first).Contains("PDF"))
            {
                for(int i = 0; i < size; i++)
                {
                    double temp_error_mass = fabs(meanMass_data_accept_sysVariation[it->first][(it->second).at(i)].at(ibin) - meanMass_data_acc_corrected.at(ibin));
                    double temp_error_pt = fabs(meanPt_data_accept_sysVariation[it->first][(it->second).at(i)].at(ibin) - meanPt_data_acc_corrected.at(ibin));
                    error_mass = error_mass < temp_error_mass?temp_error_mass:error_mass;
                    error_pt = error_pt < temp_error_pt?temp_error_pt:error_pt;

                    if((it->first).Contains("FSR"))
                    {
                        error_mass = fabs(meanMass_data_accept_sysVariation[it->first][(it->second).at(0)].at(ibin)-meanMass_data_accept_sysVariation[it->first][(it->second).at(1)].at(ibin));
                        error_pt = fabs(meanPt_data_accept_sysVariation[it->first][(it->second).at(0)].at(ibin)-meanPt_data_accept_sysVariation[it->first][(it->second).at(1)].at(ibin));
                        break;
                    }
                }
            }
            else
            {
                // PDF uncertainty
                //error_mass = meanMassPDFVariation_Accept.at(ibin)->GetRMS();
                //error_pt = meanPtPDFVariation_Accept.at(ibin)->GetRMS();

                TH1* htemp_mass = new TH1D("pdf_temp_mass", "pdf_temp_mass", 5000, 5, 500);
                TH1* htemp_pt = new TH1D("pdf_temp_pt", "pdf_temp_pt", 5000, 5, 500);
                for(int i = 0; i < size; i++)
                {
                    htemp_mass->Fill(meanMass_data_accept_sysVariation[it->first][(it->second).at(i)].at(ibin));
                    htemp_pt->Fill(meanPt_data_accept_sysVariation[it->first][(it->second).at(i)].at(ibin));
                }
                error_mass = htemp_mass->GetRMS();
                error_pt = htemp_pt->GetRMS();

                delete htemp_mass;
                delete htemp_pt;
            }
            meanMass_data_accept_systematic[it->first].push_back(error_mass);
            meanPt_data_accept_systematic[it->first].push_back(error_pt);

            meanMass_data_accept_rel_systematic[it->first].push_back(error_mass/ meanMass_data_acc_corrected.at(ibin) * 100.);
            meanPt_data_accept_rel_systematic[it->first].push_back(error_pt/ meanPt_data_acc_corrected.at(ibin) * 100.);

            //cout << ibin << " mass bin " << endl;
            //cout << "mass: " << error_mass / meanMass_data_acc_corrected.at(ibin) * 100.  << " pt: " << error_pt / meanPt_data_acc_corrected.at(ibin) * 100.<< endl;
            cout << setw(10) << error_mass ;
            cout << setw(10) << error_pt << endl;
        }
        it++;
        cout << endl;
    }

    setAcceptError();

    cout.precision(ss);
}

void ISRUnfold::setTotSysError()
{

    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;

    // For systematic
    for(int ibin = 0; ibin < nMassBin; ibin++)
    {
        double totSys_mass = 0.;
        double totSys_pt = 0.;
        std::map<TString, std::vector<TString>>::iterator it = sysMap.begin();
        while(it != sysMap.end())
        {
            double sysErr_mass = meanMass_data_folded_systematic[it->first].at(ibin);
            double sysErr_pt = meanPt_data_folded_systematic[it->first].at(ibin);

            totSys_mass += pow(sysErr_mass, 2);
            totSys_pt += pow(sysErr_pt, 2);
            it++;
        }
        meanMassSysErr_data_unfolded.push_back(sqrt(totSys_mass));
        meanPtSysErr_data_unfolded.push_back(sqrt(totSys_pt));
    }
}

void ISRUnfold::printMeanValues(bool printSys)
{
    std::streamsize ss = std::cout.precision();
    std::cout.precision(2);
    std::cout.setf( std::ios::fixed, std:: ios::floatfield );

    int size = meanMass_data_unfolded.size();

    cout << "Mean values" << endl;
    cout << setw(10) << "Mass bin" ;
    cout << setw(30) << "Mass (GeV)" ;
    cout << setw(30) << "Pt (GeV)" << endl;
    for(int i = 0; i < 70; i++) cout << "-";
    cout << endl;
    for(int i = 0; i < size; i++)
    {
        cout << setw(10) << i;
        if(printSys)
        {
            cout << setw(16) << meanMass_data_unfolded.at(i) << "+/-" << meanMassSysErr_data_unfolded.at(i) << "+/-" << meanMassStatErr_data_unfolded.at(i) ;
            cout << setw(16) << meanPt_data_unfolded.at(i) << "+/-" << meanPtSysErr_data_unfolded.at(i) << "+/-" << meanPtStatErr_data_unfolded.at(i) << endl;
        }
    }
    std::cout.precision(ss);
}

void ISRUnfold::printMeanValues_Accept(bool printSys)
{
    cout << "Acceptance (or efficiency) corrected mean values" << endl;
    std::streamsize ss = std::cout.precision();
    std::cout.precision(2);
    std::cout.setf( std::ios::fixed, std:: ios::floatfield );

    int size = meanMass_data_unfolded.size();

    cout << "Mean values" << endl;
    cout << setw(10) << "Mass bin" ;
    cout << setw(30) << "Mass (GeV)" ;
    cout << setw(30) << "Pt (GeV)" << endl;
    for(int i = 0; i < 70; i++) cout << "-";
    cout << endl;
    for(int i = 0; i < size; i++)
    {
        cout << setw(10) << i;
        if(printSys)
        {
            cout << setw(16) << meanMass_data_acc_corrected.at(i) << "+/-" << meanMassSysErr_data_acc_corrected.at(i) << "+/-" << meanMassStatErr_data_acc_corrected.at(i) ;
            cout << setw(16) << meanPt_data_acc_corrected.at(i) << "+/-" << meanPtSysErr_data_acc_corrected.at(i) << "+/-" << meanPtStatErr_data_acc_corrected.at(i) << endl;
        }
    }
    std::cout.precision(ss);
}

void ISRUnfold::setTotSysError_Accept()
{

    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;

    // For systematic
    for(int ibin = 0; ibin < nMassBin; ibin++)
    {
        double totSys_mass = 0.;
        double totSys_pt = 0.;
        std::map<TString, std::vector<TString>>::iterator it = sysMap.begin();
        while(it != sysMap.end())
        {
            double sysErr_mass = meanMass_data_accept_systematic[it->first].at(ibin);
            double sysErr_pt = meanPt_data_accept_systematic[it->first].at(ibin);

            totSys_mass += pow(sysErr_mass, 2);
            totSys_pt += pow(sysErr_pt, 2);
            it++;
        }
        meanMassSysErr_data_acc_corrected.push_back(sqrt(totSys_mass));
        meanPtSysErr_data_acc_corrected.push_back(sqrt(totSys_pt));

        meanMassRelSysErr_data_acc_corrected.push_back(sqrt(totSys_mass)/ meanMass_data_acc_corrected.at(ibin) * 100.);
        meanPtRelSysErr_data_acc_corrected.push_back(sqrt(totSys_pt)/ meanPt_data_acc_corrected.at(ibin) * 100.);

        meanMassTotErr_data_acc_corrected.push_back(sqrt(totSys_mass + pow(meanMassStatErr_data_acc_corrected.at(ibin), 2)));
        meanPtTotErr_data_acc_corrected.push_back(sqrt(totSys_pt + pow(meanPtStatErr_data_acc_corrected.at(ibin), 2)));
    }
}
void ISRUnfold::setStatError()
{
    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;

    // Loop over mass bins
    for(int ibin = 0; ibin < nMassBin; ibin++)
    {
        meanMassStatErr_data_unfolded.push_back(meanMassStatVariation.at(ibin)->GetRMS());
        meanPtStatErr_data_unfolded.push_back(meanPtStatVariation.at(ibin)->GetRMS());

        meanMassStatRelErr_data_unfolded.push_back(meanMassStatVariation.at(ibin)->GetRMS()/ meanMass_data_unfolded.at(ibin) * 100.);
        meanPtStatRelErr_data_unfolded.push_back(meanPtStatVariation.at(ibin)->GetRMS()/ meanPt_data_unfolded.at(ibin) * 100.);
    }
}

void ISRUnfold::setAcceptError()
{
    int nMassBin = massBinEdges.size() - 1;
    for(int ibin = 0; ibin < nMassBin; ibin++)
    {
        meanMassEffRelErr_data.push_back(0);
        meanPtEffRelErr_data.push_back(0);

        meanMassAcceptRelErr_data.push_back(0);
        meanPtAcceptRelErr_data.push_back(0);

        std::map<TString, std::vector<TString>>::iterator it = sysMap.begin();
        while(it != sysMap.end())
        {
            if(!(it->first).Contains("Efficiency_") && !(it->first).Contains("Acceptance_"))
            {
                it++;
                continue;
            }
            if((it->first).Contains("Efficiency_"))
            {
               meanMassEffRelErr_data.at(ibin) = sqrt(pow(meanMassEffRelErr_data.at(ibin), 2) + pow(meanMass_data_accept_systematic[it->first].at(ibin), 2)); 
               meanPtEffRelErr_data.at(ibin) = sqrt(pow(meanPtEffRelErr_data.at(ibin), 2) + pow(meanPt_data_accept_systematic[it->first].at(ibin), 2)); 
            }
            if((it->first).Contains("Acceptance_")) 
            {
               meanMassAcceptRelErr_data.at(ibin) = sqrt(pow(meanMassAcceptRelErr_data.at(ibin), 2) + pow(meanMass_data_accept_systematic[it->first].at(ibin), 2)); 
               meanPtAcceptRelErr_data.at(ibin) = sqrt(pow(meanPtAcceptRelErr_data.at(ibin), 2) + pow(meanPt_data_accept_systematic[it->first].at(ibin), 2)); 
            }
            it++;
        }
        meanMassEffRelErr_data.at(ibin) = meanMassEffRelErr_data.at(ibin) / meanMass_data_acc_corrected.at(ibin) * 100;
        meanPtEffRelErr_data.at(ibin) = meanPtEffRelErr_data.at(ibin) / meanPt_data_acc_corrected.at(ibin) * 100;

        meanMassAcceptRelErr_data.at(ibin) = meanMassAcceptRelErr_data.at(ibin) / meanMass_data_acc_corrected.at(ibin) * 100;
        meanPtAcceptRelErr_data.at(ibin) = meanPtAcceptRelErr_data.at(ibin) / meanPt_data_acc_corrected.at(ibin) * 100;
    }
}

void ISRUnfold::fillPtPDFVariationHist(int ith)
{
    //cout << "ISRUnfold::fillPtPDFVariationHist()  ith: " << ith << endl;

    // Find number of mass bins
    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;

    TString ith_;
    ith_.Form("%03d", ith);

    // Save mean pt
    for(int i = 0; i < nMassBin; i++)
    {
        TString ibinMass;
        ibinMass.Form("%d", i);

        TH1* hpt_temp_nominal = nomPtUnfold->GetOutput("hunfolded_pt_temp_nom",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
        double nominal_mean = hpt_temp_nominal->GetMean();

        if(ith == 1)
        {
            meanPtPDFVariation.push_back(new TH1F("MeanPtPDF_bin"+ibinMass, "MeanPtPDF_bin"+ibinMass, 100, nominal_mean-0.05, nominal_mean+0.05));
        }

        TH1* hpt_temp_data;

        // Get histograms to set mean values
        hpt_temp_data = sysPtUnfold["PDF"]["PDFerror"+ith_]->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
        meanPtPDFVariation.at(i)->Fill(hpt_temp_data->GetMean());

        delete hpt_temp_data;
        delete hpt_temp_nominal;
    }
}

void ISRUnfold::fillPtPDFVariationHist_Accept(int ith)
{
    //cout << "ISRUnfold::fillPtPDFVariationHist()  ith: " << ith << endl;

    // Find number of mass bins
    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;

    TString ith_;
    ith_.Form("%03d", ith);

    // Save mean pt
    for(int i = 0; i < nMassBin; i++)
    {
        TString ibinMass;
        ibinMass.Form("%d", i);

        double nominal_mean = hFullPhasePtData->GetMean();

        if(ith == 1)
        {
            meanPtPDFVariation_Accept.push_back(new TH1F("MeanPtPDF_bin"+ibinMass, "MeanPtPDF_bin"+ibinMass, 100, nominal_mean-0.05, nominal_mean+0.05));
        }

        // Get histograms to set mean values
        meanPtPDFVariation_Accept.at(i)->Fill(hSysFullPhasePtData["PDF"]["PDFerror"+ith_]->GetMean());
    }
}


void ISRUnfold::fillMassPDFVariationHist(int ith)
{
    //cout << "ISRUnfold::setMeanMass()   Save mean of dilepton..." << endl;

    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;
    const Double_t* massBins = temp_tvecd->GetMatrixArray();

    TString ith_;
    ith_.Form("%03d", ith);

    //cout << "PDFerror"+ith_ << endl;
    TH1* hunfolded_mass = sysMassUnfold["PDF"]["PDFerror"+ith_]->GetOutput("hunfolded_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);
    TH1* hmass_temp_nominal = nomMassUnfold->GetOutput("hunfolded_mass_temp_nom",0,0,"mass[UO];pt[UOC0]",kTRUE);

    // Loop over mass bins
    for(int ibin = 0; ibin < nMassBin; ibin++)
    {
        TString ibinMass;
        ibinMass.Form("%d", ibin);

        // set x-axis range
        hunfolded_mass->GetXaxis()->SetRange(hunfolded_mass->GetXaxis()->FindBin(massBins[ibin]+0.01), hunfolded_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
        hmass_temp_nominal->GetXaxis()->SetRange(hunfolded_mass->GetXaxis()->FindBin(massBins[ibin]+0.01), hunfolded_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));

        if(ith == 1)
        {
            double nominal_mean = hmass_temp_nominal->GetMean();
            meanMassPDFVariation.push_back(new TH1F("MeanMassPDF_bin"+ibinMass, "MeanMassPDF_bin"+ibinMass, 100, nominal_mean-0.05, nominal_mean+0.05));
        }
        meanMassPDFVariation.at(ibin)->Fill(hunfolded_mass->GetMean());
    }// end of mass bin loop

    delete hunfolded_mass;
    delete hmass_temp_nominal;
}

void ISRUnfold::fillMassPDFVariationHist_Accept(int ith)
{
    //cout << "ISRUnfold::setMeanMass()   Save mean of dilepton..." << endl;

    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;
    const Double_t* massBins = temp_tvecd->GetMatrixArray();

    TString ith_;
    ith_.Form("%03d", ith);

    // Loop over mass bins
    for(int ibin = 0; ibin < nMassBin; ibin++)
    {
        TString ibinMass;
        ibinMass.Form("%d", ibin);

        // set x-axis range
        hSysFullPhaseMassData["PDF"]["PDFerror"+ith_]->GetXaxis()->SetRange(hFullPhaseMassData->GetXaxis()->FindBin(massBins[ibin]+0.01), hFullPhaseMassData->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
        hFullPhaseMassData->GetXaxis()->SetRange(hFullPhaseMassData->GetXaxis()->FindBin(massBins[ibin]+0.01), hFullPhaseMassData->GetXaxis()->FindBin(massBins[ibin+1]-0.01));

        if(ith == 1)
        {
            double nominal_mean = hFullPhaseMassData->GetMean();
            meanMassPDFVariation_Accept.push_back(new TH1F("MeanMassPDF_bin"+ibinMass, "MeanMassPDF_bin"+ibinMass, 100, nominal_mean-0.05, nominal_mean+0.05));
        }
        meanMassPDFVariation_Accept.at(ibin)->Fill(hSysFullPhaseMassData["PDF"]["PDFerror"+ith_]->GetMean());
    }// end of mass bin loop
}

void ISRUnfold::fillPtStatVariationHist(int istat)
{
    //cout << "ISRUnfold::fillPtStatVariationHist()  istat: " << istat << endl;

    // Find number of mass bins
    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;

    // Save mean pt
    for(int i = 0; i < nMassBin; i++)
    {
        TString ibinMass;
        ibinMass.Form("%d", i);

        if(istat == 0)
        {
            meanPtStatVariation.push_back(new TH1F("MeanPtStat_bin"+ibinMass, "MeanPtStat_bin"+ibinMass, 40, meanPt_data_unfolded.at(i)-1., meanPt_data_unfolded.at(i)+1.));
        }

        TH1* hpt_temp_data;

        // Get histograms to set mean values
        hpt_temp_data = statPtUnfold.at(istat)->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
        meanPtStatVariation.at(i)->Fill(hpt_temp_data->GetMean());

        delete hpt_temp_data;
    }
}

void ISRUnfold::fillMassStatVariationHist(int istat)
{

    //cout << "ISRUnfold::setMeanMass()   Save mean of dilepton..." << endl;

    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;
    const Double_t* massBins = temp_tvecd->GetMatrixArray();

    TH1* hunfolded_mass = statMassUnfold.at(istat)->GetOutput("hunfolded_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);

    // Loop over mass bins
    for(int ibin = 0; ibin < nMassBin; ibin++)
    {
        TString ibinMass;
        ibinMass.Form("%d", ibin);

        // set x-axis range
        hunfolded_mass->GetXaxis()->SetRange(hunfolded_mass->GetXaxis()->FindBin(massBins[ibin]+0.01), hunfolded_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));

        if(istat == 0)
        {
            meanMassStatVariation.push_back(new TH1F("MeanMassStat_bin"+ibinMass, "MeanMassStat_bin"+ibinMass, 80, meanMass_data_unfolded.at(ibin)-2., meanMass_data_unfolded.at(ibin)+2.));
        }
        meanMassStatVariation.at(ibin)->Fill(hunfolded_mass->GetMean());
    }// end of mass bin loop

    delete hunfolded_mass;
}

int ISRUnfold::setMeanPt(TString filePath)
{
    TString DYHistName_ = "histo_DYJetsToMuMu";
    if(channel_name == "electron")
    {
        DYHistName_  = "histo_DYJetsToEE";
    }

    // Find number of mass bins
    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;

    TUnfoldDensity* p_unfold = NULL;
    p_unfold = nomPtUnfold;

    // Save mean pt
    for(int i = 0; i < nMassBin; i++)
    {
        TString ibinMass;
        ibinMass.Form("%d", i);

        // Get detector level MC
        if(filePath != "")
        {
            TH1* hDY = getRawHist("Pt_FineCoarse", filePath, "Detector", DYHistName_, "Signal", "pt[UO];mass[UOC"+ibinMass+"]", true, false);
            //cout << "detector mc, pt: " << hDY->GetMean() << " err: " << hDY->GetMeanError() << endl;
            //cout << "n bins: " << hDY->GetNbinsX() << endl;
            meanPt_mc_folded.push_back(hDY->GetMean());
            meanPtStatErr_mc_folded.push_back(hDY->GetMeanError());
            delete hDY;
        }

        TH1* hpt_temp_data;
        TH1* hpt_temp_mc;

        // get histograms to set mean values
        TH1* hdetector_data = p_unfold->GetInput("h_folded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
        hpt_temp_data = p_unfold->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
        hpt_temp_mc   = p_unfold->GetBias("histMC_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);

        //cout << "Detector, " << i << " th mass bin, mean: " << hdetector_data->GetMean() << " +/- " << hdetector_data->GetMeanError() << endl;
        meanPt_data_folded.push_back(hdetector_data->GetMean());
        meanPtStatErr_data_folded.push_back(hdetector_data->GetMeanError());

        //cout << "Unfolded, " << i << " th mass bin, mean: " << hpt_temp_data->GetMean() << " +/- " << hpt_temp_data->GetMeanError() << endl;
        meanPt_data_unfolded.push_back(hpt_temp_data->GetMean());
        meanPtStatErr_data_unfolded.push_back(hpt_temp_data->GetMeanError());

        meanPt_mc_unfolded.push_back(hpt_temp_mc->GetMean());
        meanPtStatErr_mc_unfolded.push_back(hpt_temp_mc->GetMeanError());

        delete hdetector_data;
        delete hpt_temp_data;
        delete hpt_temp_mc;
    }

    return nMassBin;
}

int ISRUnfold::setSysMeanPt()
{
    //cout << "ISRUnfold::setSysMeanPt()   Save mean of dilepton momentum..." << endl;

    // Find number of mass bins
    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;

    std::map<TString, std::vector<TString>>::iterator it = sysMap.begin();
    while(it != sysMap.end())
    {
        TUnfoldDensity* p_unfold = NULL;

        //cout << "Unfold for " << it->first << " systematic." << endl;
        int size = (it->second).size();
        //cout << size << " systematic variation exist." << endl;

        for(int i = 0; i < size; i++)
        {
            if((it->first).Contains("Unfolding") && !((it->second).at(i)).Contains("Nominal"))
            {
                p_unfold = NULL;
            }
            else
            {
                if((it->first).Contains("Unfolding"))
                {
                    p_unfold = nomPtUnfold;
                }
                else
                {
                    p_unfold = sysPtUnfold[it->first][(it->second).at(i)];
                }
            }
            // Save mean pt
            for(int j = 0; j < nMassBin; j++)
            {
                TString ibinMass;
                ibinMass.Form("%d", j);

                TH1* hunfolded_pt = NULL;
                TH1* hunfolded_pt_input = NULL;

                // Get histograms to set mean values
                if((it->first).Contains("Unfolding") && !((it->second).at(i)).Contains("Nominal"))
                {
                    hunfolded_pt = iterEMPtUnfold->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
                    hunfolded_pt_input = nomPtUnfold->GetInput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
                }
                else
                {
                    hunfolded_pt = p_unfold->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
                    hunfolded_pt_input = p_unfold->GetInput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
                }
                meanPt_data_unfolded_sysVariation[it->first][(it->second).at(i)].push_back(hunfolded_pt->GetMean());
                meanPt_data_folded_sysVariation[it->first][(it->second).at(i)].push_back(hunfolded_pt_input->GetMean());
                //cout << it->first << " " << (it->second).at(i) << " " << hunfolded_pt->GetMean() << endl;
                delete hunfolded_pt;
                delete hunfolded_pt_input;
            }
        }
        it++;
    }

    return nMassBin;
}

void ISRUnfold::setSysMeanPt_Accept()
{

    //cout << "ISRUnfold::setSysMeanPt_Accept()   Save mean of dilepton momentum..." << endl;

    // Find number of mass bins
    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;

    std::map<TString, std::vector<TString>>::iterator it = sysMap.begin();
    while(it != sysMap.end())
    {
        int size = (it->second).size();

        for(int i = 0; i < size; i++)
        {
            // Save mean pt
            for(int j = 0; j < nMassBin; j++)
            {
                TString ibinMass;
                ibinMass.Form("%d", j);

                TH1* hunfolded_pt_raw = hSysFullPhasePtData[it->first][(it->second).at(i)];
                TH1* hunfolded_pt = pt_binning_Gen->ExtractHistogram("temp_pt", hunfolded_pt_raw, 0, kTRUE, "pt[UO];mass[UOC"+ibinMass+"]");

                // Get histograms to set mean values
                meanPt_data_accept_sysVariation[it->first][(it->second).at(i)].push_back(hunfolded_pt->GetMean());
                //cout << it->first << " " << (it->second).at(i) << " " << hunfolded_pt->GetMean() << endl;

                delete hunfolded_pt;
            }
        }
        it++;
    }
}

vector<double> ISRUnfold::getFoldedMeanPtVectors(TString sysName, TString variationName)
{
    if(sysName=="")
    {
        return meanPt_data_folded;
    }
    else
        return meanPt_data_folded_sysVariation[sysName][variationName];
}

vector<double> ISRUnfold::getFoldedMeanMassVectors(TString sysName, TString variationName)
{
    if(sysName=="")
    {
        return meanMass_data_folded;
    }
    else
        return meanMass_data_folded_sysVariation[sysName][variationName];
}


vector<double> ISRUnfold::getUnfoldedMeanPtVectors(TString sysName, TString variationName)
{
    if(sysName=="")
    {
        return meanPt_data_unfolded;
    }
    else
        return meanPt_data_unfolded_sysVariation[sysName][variationName];
}

vector<double> ISRUnfold::getUnfoldedMeanMassVectors(TString sysName, TString variationName)
{
    if(sysName=="")
    {
        return meanMass_data_unfolded;
    }
    else
        return meanMass_data_unfolded_sysVariation[sysName][variationName];
}


vector<double> ISRUnfold::getAccCorrectedMeanPtVectors(TString sysName, TString variationName)
{
    if(sysName=="")
    {
        return meanPt_data_acc_corrected;
    }
    else
        return meanPt_data_accept_sysVariation[sysName][variationName];
}

vector<double> ISRUnfold::getAccCorrectedMeanMassVectors(TString sysName, TString variationName)
{
    if(sysName=="")
    {
        return meanMass_data_acc_corrected;
    }
    else
        return meanMass_data_accept_sysVariation[sysName][variationName];
}

void ISRUnfold::doNorm(TH1* hist, bool norm)
{
    for(int ibin = 1; ibin < hist->GetXaxis()->GetNbins()+1; ibin++)
    {
        double binWidth = hist->GetBinWidth(ibin);
        hist->SetBinContent(ibin, hist->GetBinContent(ibin)/ binWidth);
        hist->SetBinError(ibin, hist->GetBinError(ibin)/ binWidth);
    }
    if(norm)
        hist->Scale(1./ hist->Integral());
}

TH1* ISRUnfold::getUnfoldedHists(TString var, TString outHistName, TString steering, bool useAxis, bool binWidth)
{
    TH1* outHist = NULL;
    if(var == "Mass")
    {
        outHist = nomMassUnfold->GetOutput(outHistName,0,0,steering,useAxis);
        if(binWidth)
        {
            divideByBinWidth(outHist, false);
        }
        return outHist;
    }
    else
    {
        outHist = nomPtUnfold->GetOutput(outHistName,0,0,steering,useAxis);
        if(binWidth)
        {
            divideByBinWidth(outHist, false);
        }
        else
            return outHist;
    }
    return outHist;
}

TH1* ISRUnfold::getRawHist(TString var, TString filePath, TString dirName, TString histName, TString outHistName, TString steering, bool useAxis, bool divBinWidth)
{
    TH1::AddDirectory(kFALSE);
    TFile* filein = new TFile(filePath);
    TH1* hist = NULL;

    if(steering != "")
    {
        TH1* raw_hist = (TH1*)filein->Get(dirName+"/"+var+"/"+histName);
        if(histName.Contains("DYJetsTo") && !histName.Contains("Tau"))
        {
            histName.ReplaceAll("DYJetsTo", "DYJets10to50To");
            raw_hist->Add((TH1*)filein->Get(dirName+"/"+var+"/"+histName));
        }

        if(var.Contains("Pt"))
        {
            hist = pt_binning_Rec->ExtractHistogram(outHistName, raw_hist, 0, useAxis, steering);
        }
        else
        {
            hist = mass_binning_Rec->ExtractHistogram(outHistName, raw_hist, 0, useAxis, steering);
        }

        delete raw_hist;
    }
    else
    {
        //cout << "Steering not provided, get raw histogram." << endl;
        //cout << dirName+"/"+var+"/"+histName << endl;
        hist = (TH1*)filein->Get(dirName+"/"+var+"/"+histName);
        if(histName.Contains("DYJetsTo") && !histName.Contains("Tau"))
        {
            histName.ReplaceAll("DYJetsTo", "DYJets10to50To");
            hist->Add((TH1*)filein->Get(dirName+"/"+var+"/"+histName));
        }
    }

    delete filein;

    if(divBinWidth)
        divideByBinWidth(hist, false);
    return hist;
}

