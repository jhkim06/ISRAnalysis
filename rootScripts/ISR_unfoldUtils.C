#include "ISR_unfoldUtils.h"

void ISRUnfold::setBias(double bias)
{
   nominal_bias = bias;
}

void ISRUnfold::setOutputBaseDir(TString outPath)
{
    output_baseDir = outPath;
}

const TVectorD& ISRUnfold::checkMatrixCond()
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

double ISRUnfold::getSmearedChi2(TString filePath, TString dirName, TString steering, bool useAxis, bool divBinWidth)
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
        hDY = getRawHist(filePath, dirName, DYHistName_, "Signal_"+var, steering, useAxis, divBinWidth); ;
    }
    else
    {
        hData = nomPtUnfold->GetInput("hData_"+var, 0, 0, steering, useAxis);
        hDY = getRawHist(filePath, dirName, DYHistName_, "Signal_"+var, steering, useAxis, divBinWidth); ;
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

double ISRUnfold::getUnfoldedChi2(TString steering, bool useAxis, bool divBinWidth)
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

void ISRUnfold::drawResponseM(TString sysName, TString sysPostfix, bool isDetector)
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
void ISRUnfold::setNominalRM(TString filepath, TString dirName, TString histName, TString binDef)
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
                    if(var=="Pt") this->iterEMPtUnfold   = new TUnfoldIterativeEM(hPtResponseM,TUnfoldDensity::kHistMapOutputHoriz,pt_binning_Gen,pt_binning_Rec);
                    if(var=="Mass") this->iterEMMassUnfold = new TUnfoldIterativeEM(hMassResponseM,TUnfoldDensity::kHistMapOutputHoriz,mass_binning_Gen,mass_binning_Rec);

                    if(!useAccept)
                    {
                        if(var=="Pt") this->iterEMPtUnfold->SetInput(unfold->iterEMPtUnfold->GetOutput("hUnfoldedPt_"+ it->first + "_" + sysMap_previous[it->first][ith],0,0,"*[*]",false), nominal_bias);
                        if(var=="Mass") this->iterEMMassUnfold->SetInput(unfold->iterEMMassUnfold->GetOutput("hUnfoldedMass_"+ it->first + "_" + sysMap_previous[it->first][ith],0,0,"*[*]",false), nominal_bias);
                    }
                    else
                    {
                        //cout << "use acceptance corrected output!" << endl;
                        if(var=="Pt") this->iterEMPtUnfold->SetInput(unfold->hSysFullPhasePtData[it->first][sysMap_previous[it->first][ith]], nominal_bias);
                        if(var=="Mass") this->iterEMMassUnfold->SetInput(unfold->hSysFullPhaseMassData[it->first][sysMap_previous[it->first][ith]], nominal_bias);
                    }
                }
                else
                {
                    if(var=="Pt") this->sysPtUnfold[it->first][sysMap_previous[it->first][ith]]   = new TUnfoldDensity(hPtResponseM,TUnfold::kHistMapOutputHoriz,regMode, TUnfold::kEConstraintArea, TUnfoldDensity::kDensityModeNone, pt_binning_Gen,pt_binning_Rec);
                    if(var=="Mass") this->sysMassUnfold[it->first][sysMap_previous[it->first][ith]] = new TUnfoldDensity(hMassResponseM,TUnfold::kHistMapOutputHoriz,regMode, TUnfold::kEConstraintArea, TUnfoldDensity::kDensityModeNone,mass_binning_Gen, mass_binning_Rec);

                    if(!useAccept)
                    {
                        if(var=="Pt") this->sysPtUnfold[it->first][sysMap_previous[it->first][ith]]->SetInput(unfold->sysPtUnfold[it->first][sysMap_previous[it->first][ith]]->GetOutput("hUnfoldedPt_"+ it->first + "_" + sysMap_previous[it->first][ith],0,0,"*[*]",false), nominal_bias);
                        if(var=="Mass") this->sysMassUnfold[it->first][sysMap_previous[it->first][ith]]->SetInput(unfold->sysMassUnfold[it->first][sysMap_previous[it->first][ith]]->GetOutput("hUnfoldedMass_"+ it->first + "_" + sysMap_previous[it->first][ith],0,0,"*[*]",false), nominal_bias);
                    }
                    else
                    {
                        //cout << "use acceptance corrected output!" << endl;
                        if(var=="Pt") this->sysPtUnfold[it->first][sysMap_previous[it->first][ith]]->SetInput(unfold->hSysFullPhasePtData[it->first][sysMap_previous[it->first][ith]], nominal_bias);
                        if(var=="Mass") this->sysMassUnfold[it->first][sysMap_previous[it->first][ith]]->SetInput(unfold->hSysFullPhaseMassData[it->first][sysMap_previous[it->first][ith]], nominal_bias);
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
                    if(var=="Pt") this->sysPtUnfold[it->first][sysMap_previous[it->first][ith]]->SetInput(unfold->sysPtUnfold[it->first][sysMap_previous[it->first][ith]]->GetOutput("hUnfoldedPt_"+ it->first + "_" + sysMap_previous[it->first][ith],0,0,"*[*]",false), nominal_bias);
                    if(var=="Mass") this->sysMassUnfold[it->first][sysMap_previous[it->first][ith]]->SetInput(unfold->sysMassUnfold[it->first][sysMap_previous[it->first][ith]]->GetOutput("hUnfoldedMass_"+ it->first + "_" + sysMap_previous[it->first][ith],0,0,"*[*]",false), nominal_bias);
                }
                else
                {

                    if((it->first).Contains("Unfolding") && !(sysMap_previous[it->first][ith]).Contains("Nominal"))
                    {
                        if(var=="Pt") this->iterEMPtUnfold->SetInput(unfold->hSysFullPhasePtData[it->first][sysMap_previous[it->first][ith]], nominal_bias);
                        if(var=="Mass") this->iterEMMassUnfold->SetInput(unfold->hSysFullPhaseMassData[it->first][sysMap_previous[it->first][ith]], nominal_bias);
                    }
                    else
                    {
                        if(var=="Pt") this->sysPtUnfold[it->first][sysMap_previous[it->first][ith]]->SetInput(unfold->hSysFullPhasePtData[it->first][sysMap_previous[it->first][ith]], nominal_bias);
                        if(var=="Mass") this->sysMassUnfold[it->first][sysMap_previous[it->first][ith]]->SetInput(unfold->hSysFullPhaseMassData[it->first][sysMap_previous[it->first][ith]], nominal_bias);
                    }
                }
            }
        }
        it++;
    }
    // Loop over variations of event selection efficiency correction
}

void ISRUnfold::setSystematicRM(TString filepath, TString dirName, TString histName, TString sysName, TString sysPostfix, TString histPostfix, TString binDef)
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
void ISRUnfold::setUnfInput(ISRUnfold* unfold, bool isSys, TString sysName, TString sysPostfix, bool useAccept)
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
void ISRUnfold::setUnfInput(TString varPostfix, TString filepath, TString dirName, TString histName, bool isSys, TString sysName, TString sysPostfix, bool isFSR)
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
        if(var == "Pt")
        {
            hPtRec = (TH1*)filein->Get(dirName + "/Pt"+binDef+"/histo_" + bkgInfo.first);
            nomPtUnfold->  SubtractBackground(hPtRec, bkgInfo.first);
        }

        if(var == "Mass")
        {
            hMassRec = (TH1*)filein->Get(dirName + "/Mass"+binDef+"/histo_" + bkgInfo.first);
            nomMassUnfold->SubtractBackground(hMassRec, bkgInfo.first);
        }
    }
    else
    // Systematic
    {
        TString fullHistName = bkgInfo.first + "_" + sysPostfix;
        if(histPostfix == "")
            fullHistName = bkgInfo.first;

        if(var == "Pt") hPtRec = (TH1*)filein->Get(dirName + "/Pt"+binDef+"/histo_" + fullHistName);
        if(var == "Mass") hMassRec = (TH1*)filein->Get(dirName + "/Mass"+binDef+"/histo_" + fullHistName);

        //cout << "file path: " << filepath << endl;
        //cout << dirName + "/Pt"+binDef+"/histo_" + fullHistName << endl;

        if(sysName.Contains("Unfolding") && !sysPostfix.Contains("Nominal"))
        {
            if(var == "Pt") iterEMPtUnfold->SubtractBackground(hPtRec, bkgInfo.first);
            if(var == "Mass") iterEMMassUnfold->SubtractBackground(hMassRec, bkgInfo.first);
        }
        else
        {
            // FIXME temporary method for background systematic
            if(sysName.Contains("Background"))
            {
                if(sysPostfix.Contains("Up"))
                {
                    if(var == "Pt") sysPtUnfold[sysName][sysPostfix]->SubtractBackground(hPtRec, bkgInfo.first, 1.05);
                    if(var == "Mass") sysMassUnfold[sysName][sysPostfix]->SubtractBackground(hMassRec, bkgInfo.first, 1.05);
                }
                if(sysPostfix.Contains("Down"))
                {
                    if(var == "Pt") sysPtUnfold[sysName][sysPostfix]->SubtractBackground(hPtRec, bkgInfo.first, 0.95);
                    if(var == "Mass") sysMassUnfold[sysName][sysPostfix]->SubtractBackground(hMassRec, bkgInfo.first, 0.95);
                }
            }
            else
            {
                if(var == "Pt") sysPtUnfold[sysName][sysPostfix]->SubtractBackground(hPtRec, bkgInfo.first);
                if(var == "Mass") sysMassUnfold[sysName][sysPostfix]->SubtractBackground(hMassRec, bkgInfo.first);
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

        //fillPtStatVariationHist(istat);
        //fillMassStatVariationHist(istat);

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
    TDirectory* varDir;

    TString fullPath=output_baseDir+unfold_name+"_"+channel_name+"_"+yearStr+"_"+var+".root";
    cout << fullPath << endl;
    if(!gSystem->AccessPathName(fullPath, kFileExists))
    {
        f=TFile::Open(fullPath, "UPDATE"); 

        topDir=f->GetDirectory("unfolded");
        varDir=f->GetDirectory("unfolded/"+var);
    }
    else
    {
        f=new TFile(fullPath, "CREATE");

        // Create directory
        topDir=f->mkdir("unfolded");
        varDir=topDir->mkdir(var); 

        varDir->cd();
        if(var=="Pt") pt_binning_Gen->Write();
        if(var=="Mass") mass_binning_Gen->Write();
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
            if(var=="Pt") nomPtUnfold->DoUnfold(0);
            if(var=="Mass") nomMassUnfold->DoUnfold(0);
        }
        else
        {
            if(var=="Mass") nomMassUnfold->DoUnfold(0);

            if(var=="Pt") 
            {
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
        }

        varDir->cd(); 
        if(var=="Pt")
        {
        nomPtUnfold->GetOutput("histo_Data",0,0, "*[*]", false)->Write(); 
        nomPtUnfold->GetBias("histo_DY", 0, 0, "*[*]", false)->Write(); 
        }
        if(var=="Mass")
        {
        nomMassUnfold->GetOutput("histo_Data",0,0, "*[*]", false)->Write(); 
        nomMassUnfold->GetBias("histo_DY", 0, 0, "*[*]", false)->Write(); 
        }
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
                    if(var=="Pt")  iBest_pt=iterEMPtUnfold->ScanSURE(NITER_Iterative, &graph_SURE_IterativeSURE_pt, &graph_DFdeviance_IterativeSURE_pt);
                    if(var=="Mass")  iBest_mass=iterEMMassUnfold->ScanSURE(NITER_Iterative, &graph_SURE_IterativeSURE_mass, &graph_DFdeviance_IterativeSURE_mass);
                    //cout << "iBest pt, Mass: " << iBest_pt << " " << iBest_mass << endl;

                    varDir->cd();
                    if(var=="Pt") 
                    {
                    iterEMPtUnfold->GetOutput("histo_Data_"+(it->second).at(i),0,0, "*[*]", false)->Write(); 
                    nomPtUnfold->GetBias("histo_DY_"+(it->second).at(i), 0, 0, "*[*]", false)->Write(); 
                    }
                    if(var=="Mass") 
                    {
                    iterEMMassUnfold->GetOutput("histo_Data_"+(it->second).at(i),0,0, "*[*]", false)->Write(); 
                    nomMassUnfold->GetBias("histo_DY_"+(it->second).at(i), 0, 0, "*[*]", false)->Write(); 
                    }
                }
                else
                {

                    if(regMode == TUnfold::kRegModeNone)
                    {
                        if(var=="Pt") sysPtUnfold[it->first][(it->second).at(i)]->DoUnfold(0);
                        if(var=="Mass")  sysMassUnfold[it->first][(it->second).at(i)]->DoUnfold(0);
                    }
                    else
                    {
                        double tauMin=1.e-4;
                        double tauMax=1.e-1;
                        if(var=="Pt")  sysPtUnfold[it->first][(it->second).at(i)]->ScanLcurve(100, tauMin, tauMax, 0);
                        if(var=="Mass")  sysMassUnfold[it->first][(it->second).at(i)]->DoUnfold(0);
                    }

                    varDir->cd();
                    if(var=="Pt") 
                    {
                    sysPtUnfold[it->first][(it->second).at(i)]->GetOutput("histo_Data_"+(it->second).at(i),0,0, "*[*]", false)->Write(); 
                    sysPtUnfold[it->first][(it->second).at(i)]->GetBias("histo_DY_"+(it->second).at(i), 0, 0, "*[*]", false)->Write(); 
                    }
                    if(var=="Mass")
                    { 
                    sysMassUnfold[it->first][(it->second).at(i)]->GetOutput("histo_Data_"+(it->second).at(i),0,0, "*[*]", false)->Write(); 
                    sysMassUnfold[it->first][(it->second).at(i)]->GetBias("histo_DY_"+(it->second).at(i), 0, 0, "*[*]", false)->Write(); 
                    }
                }
            }
            it++;
        }
    }// Unfold for systematic
    topDir->Write();
    f->Close();
}

void ISRUnfold::drawCorrelation(TString steering, bool useAxis, TString outName)
{

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
    TDirectory* varDir;

    TString fullPath=output_baseDir+unfold_name+"_Acceptance_"+channel_name+"_"+yearStr+"_"+var+".root";
    if(!gSystem->AccessPathName(fullPath, kFileExists))
    {

        f=TFile::Open(fullPath, "UPDATE"); 

        topDir=f->GetDirectory("acceptance");
        varDir=f->GetDirectory("acceptance/"+var);
    }
    else
    {

        f=new TFile(fullPath, "CREATE");

        // Create directory
        topDir=f->mkdir("acceptance");
        varDir=topDir->mkdir(var);

        varDir->cd();
        if(var=="Pt") pt_binning_Gen->Write();
        if(var=="Mass") mass_binning_Gen->Write();
    }

    TFile* filein = new TFile(filePath);

    TString accepCorrOrEffCorr;
    if(isAccept)
        accepCorrOrEffCorr = "Acceptance";
    else
        accepCorrOrEffCorr = "Efficiency";

    TH1* hFiducialPhaseMassMC = NULL;
    TH1* hFiducialPhasePtMC = NULL;

    if(var=="Mass")
    {
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

    varDir->cd();
    hFullPhaseMassData->Write();
    hFullPhaseMassMC->SetName("histo_DY");
    hFullPhaseMassMC->Write();
    hAcceptanceMass->Write();
    }

    if(var=="Pt")
    {
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

    varDir->cd();
    hFullPhasePtData->Write();
    hFullPhasePtMC->SetName("histo_DY");
    hFullPhasePtMC->Write();
    hAcceptancePt->Write();
    }

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
                    if(var=="Mass")
                    {
                    hSysFullPhaseMassData[it->first][(it->second).at(i)] = iterEMMassUnfold->GetOutput("histo_Data_"+(it->second).at(i),0,0, "*[*]", false);
                    hFiducialPhaseMassMC_sys=hFiducialPhaseMassMC;
                    hFiducialPhaseMassMC_sys->SetName("histo_DY_"+(it->second).at(i));
                    }
    
                    if(var=="Pt")
                    {
                    hSysFullPhasePtData[it->first][(it->second).at(i)]   = iterEMPtUnfold->GetOutput("histo_Data_"+(it->second).at(i),0,0, "*[*]", false);
                    hFiducialPhasePtMC_sys=hFiducialPhasePtMC;
                    hFiducialPhasePtMC_sys->SetName("histo_DY_"+(it->second).at(i));
                    }
                }
                else
                {
                    if(var=="Mass")
                    {
                    hSysFullPhaseMassData[it->first][(it->second).at(i)] = sysMassUnfold[it->first][(it->second).at(i)]->GetOutput("histo_Data_"+(it->second).at(i),0,0, "*[*]", false);
                    hFiducialPhaseMassMC_sys = sysMassUnfold[it->first][(it->second).at(i)]->GetBias("histo_DY_"+(it->second).at(i), 0, 0, "*[*]", false);
                    }
            
                    if(var=="Pt")
                    {
                    hSysFullPhasePtData[it->first][(it->second).at(i)]   = sysPtUnfold[it->first][(it->second).at(i)]->GetOutput("histo_Data_"+(it->second).at(i),0,0, "*[*]", false);
                    hFiducialPhasePtMC_sys = sysPtUnfold[it->first][(it->second).at(i)]->GetBias("hFiducialPt_sys"+(it->second).at(i), 0, 0, "*[*]", false);
                    }
                }

                // For PDF, AlphaS, Scale etc, denominator changed
                if( (((it->first).Contains("Scale") && !(it->first).Contains("Lep")) || (it->first).Contains("PDF") || (it->first).Contains("AlphaS")) && !(it->first).Contains("_") )
                {
                    if(var=="Mass")
                    {
                    hFullPhaseMassMC_raw_sys = (TH1*) filein->Get("Acceptance/MassGen" + binDef + "/histo_DYJets_"+(it->second).at(i));
                    if(year==2016)
                        hFullPhaseMassMC_raw_sys->Add((TH1*) filein->Get("Acceptance/MassGen" + binDef + "/histo_DYJets10to50_"+(it->second).at(i)));
                    else
                        hFullPhaseMassMC_raw_sys->Add((TH1*) filein->Get("Acceptance/MassGen" + binDef + "/histo_DYJets10to50_MG_"+(it->second).at(i)));
                    }
    
                    if(var=="Pt")
                    {
                    hFullPhasePtMC_raw_sys = (TH1*) filein->Get("Acceptance/PtGen" + binDef + "/histo_DYJets_"+(it->second).at(i));
                    if(year==2016)
                        hFullPhasePtMC_raw_sys->Add((TH1*) filein->Get("Acceptance/PtGen" + binDef + "/histo_DYJets10to50_"+(it->second).at(i)));
                    else
                        hFullPhasePtMC_raw_sys->Add((TH1*) filein->Get("Acceptance/PtGen" + binDef + "/histo_DYJets10to50_MG_"+(it->second).at(i)));
                    }
                }
                else
                {
                    if(var=="Mass") hFullPhaseMassMC_raw_sys=hFullPhaseMassMC;
                    if(var=="Pt") hFullPhasePtMC_raw_sys=hFullPhasePtMC;
                }

                if(var=="Mass")
                {
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
                }

                if(var=="Pt")
                {
                // For pt
                TH1* hAcceptancePt_sys = (TH1*) hFullPhasePtMC_raw_sys->Clone("hAcceptancePt_sys");
                hAcceptancePt_sys->Divide(hFiducialPhasePtMC_sys);

                TH1* hAcceptanceFractionPt_sys = (TH1*) hFiducialPhasePtMC_sys->Clone("hAcceptanceFractionPt_sys");
                hAcceptanceFractionPt_sys->Divide(hFullPhasePtMC_raw_sys);

                //hSysFullPhasePtData[accepCorrOrEffCorr + "_" + it->first][(it->second).at(i)] = nomPtUnfold->GetOutput("hAcceptPtData" +it->first+(it->second).at(i),0,0, "*[*]", false);
                hSysFullPhasePtData[it->first][(it->second).at(i)]->Multiply(hAcceptancePt_sys);
                hSysFullPhasePtMC[it->first][(it->second).at(i)] = hFullPhasePtMC_raw_sys;
                sysMapForAcceptance[it->first].push_back((it->second).at(i)); // Update sysMapForAcceptance 

                //hSysAcceptancePt[it->first][(it->second).at(i)] = (TH1*) hAcceptancePt_sys->Clone("Pt_" + it->first + "_" + (it->second).at(i));
                hSysAcceptanceFractionPt[it->first][(it->second).at(i)] = (TH1*) hAcceptanceFractionPt_sys->Clone("FractionPt_" + it->first + "_" + (it->second).at(i));
                delete hAcceptancePt_sys;
                delete hAcceptanceFractionPt_sys;
                }
                //hSysFullPhaseMassData[it->first][(it->second).at(i)]->Multiply(hAcceptanceMass);
                //hSysFullPhasePtData[it->first][(it->second).at(i)]->Multiply(hAcceptancePt);
                if(var=="Mass") hSysFullPhaseMassMC[it->first][(it->second).at(i)] = hFullPhaseMassMC_raw_sys;
                if(var=="Pt") hSysFullPhasePtMC[it->first][(it->second).at(i)]   = hFullPhasePtMC_raw_sys;

                varDir->cd();
                if(var=="Pt")
                {
                hSysFullPhasePtData[it->first][(it->second).at(i)]->Write();
                hFullPhasePtMC_raw_sys->SetName("histo_DY_"+(it->second).at(i));
                hFullPhasePtMC_raw_sys->Write();
                hSysAcceptancePt[it->first][(it->second).at(i)] = (TH1*) hFullPhasePtMC->Clone("Pt_" + it->first + "_" + (it->second).at(i));
                }
        
                if(var=="Mass")
                {
                hSysFullPhaseMassData[it->first][(it->second).at(i)]->Write();
                hFullPhaseMassMC_raw_sys->SetName("histo_DY_"+(it->second).at(i));
                hFullPhaseMassMC_raw_sys->Write();
                hSysAcceptanceMass[it->first][(it->second).at(i)] = (TH1*) hAcceptanceMass->Clone("Mass_" + it->first + "_" + (it->second).at(i));
                }
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

        if(var=="Mass")
        {
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
        }
   
        if(var=="Pt")
        { 
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
        }
        //setSystematics(accepCorrOrEffCorr, "StatUp");
        //setSystematics(accepCorrOrEffCorr, "StatDown");
    }

    topDir->Write();;
    f->Close();

    if(var=="Mass") delete hFiducialPhaseMassMC;
    if(var=="Pt") delete hFiducialPhasePtMC;
}

void ISRUnfold::varyHistWithStatError(TH1* hist, int sys)
{
    for(int ibin = 1; ibin < hist->GetNbinsX()+1; ibin++)  
    {
        hist->SetBinContent(ibin, hist->GetBinContent(ibin) + double(sys) * hist->GetBinError(ibin));
    }
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

TH1* ISRUnfold::getUnfoldedHists(TString outHistName, TString steering, bool useAxis, bool binWidth)
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

TH1* ISRUnfold::getRawHist(TString filePath, TString dirName, TString histName, TString outHistName, TString steering, bool useAxis, bool divBinWidth)
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

