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
    TFile f("./Matrix_"+var+yearStr+".root","recreate");

    sysPostfix = "";
    sysName = "";
    const TVectorD* xaxis1_tvecd;
    int xaxis1_nbin;

    const TVectorD* yaxis1_tvecd = NULL;
    int yaxis1_nbin;

    setTDRStyle();
    writeExtraText = true;
    extraText  = "simulation";
    gStyle->SetLineWidth(5.);
    gStyle->SetFrameLineWidth(5.);
    gROOT->ForceStyle();

    TCanvas* c1 = new TCanvas("c1","c1", 50, 50, 2400, 2400);
    gStyle->SetOptFit(0);
    gStyle->SetPalette(55);
    c1->cd();

    c1->SetTopMargin(0.1);
    c1->SetRightMargin(0.15);
    c1->SetTicks(1);
    c1->SetLogz();
    //if(var=="Mass")
    //{
    //    c1->SetLogx();
    //    c1->SetLogy();
    //}

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

    TString draw_option = "COLZ";

    gStyle->SetPaintTextFormat("1.2f");
    histProb_woUO->SetMarkerSize(histProb_woUO->GetMarkerSize()*0.3);
    histProb_woUO->Draw("COLZ text");
    histProb_woUO->GetZaxis()->SetRangeUser(5e-3, 1.0);
    histProb_woUO->GetYaxis()->SetTitleFont(43);
    histProb_woUO->GetYaxis()->SetTitleSize(100);
    histProb_woUO->GetYaxis()->SetTitleOffset(1.5);
    histProb_woUO->GetXaxis()->SetTitleFont(43);
    histProb_woUO->GetXaxis()->SetTitleSize(100);
    histProb_woUO->GetXaxis()->SetTitleOffset(1.5);

    vector<TGaxis*> v_xaxis, v_yaxis;

    if(var.Contains("Pt"))
    {
        histProb_woUO->GetXaxis()->SetNdivisions(506);
        histProb_woUO->GetYaxis()->SetNdivisions(506);

        histProb_woUO->GetXaxis()->ChangeLabel(2,-1,-1,-1,-1,-1,"100");
        histProb_woUO->GetXaxis()->ChangeLabel(3,-1,-1,-1,-1,-1,"100");
        histProb_woUO->GetXaxis()->ChangeLabel(4,-1,-1,-1,-1,-1,"100");
        histProb_woUO->GetXaxis()->ChangeLabel(5,-1,-1,-1,-1,-1,"100");
        histProb_woUO->GetXaxis()->ChangeLabel(6,-1,-1,-1,-1,-1,"100");

        histProb_woUO->GetYaxis()->ChangeLabel(2,-1,-1,-1,-1,-1,"100");
        histProb_woUO->GetYaxis()->ChangeLabel(3,-1,-1,-1,-1,-1,"100");
        histProb_woUO->GetYaxis()->ChangeLabel(4,-1,-1,-1,-1,-1,"100");
        histProb_woUO->GetYaxis()->ChangeLabel(5,-1,-1,-1,-1,-1,"100");
        histProb_woUO->GetYaxis()->ChangeLabel(6,-1,-1,-1,-1,-1,"100");

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

    if(isDetector)
    {
        if(var=="Pt")
        {
        histProb_woUO->GetYaxis()->SetTitle("Detector p_{T} [GeV]");
        histProb_woUO->GetXaxis()->SetTitle("Generator p_{T} (dressed level) [GeV]");
        }
        else
        {
        histProb_woUO->GetYaxis()->SetTitle("Detector mass [GeV]");
        histProb_woUO->GetXaxis()->SetTitle("Generator mass (dressed level) [GeV]");
        }
    }
    else
    {
        if(var=="Pt")
        {
        histProb_woUO->GetYaxis()->SetTitle("Generator p_{T} (dressed level) [GeV]");
        histProb_woUO->GetXaxis()->SetTitle("Generator p_{T} (parton level) [GeV]");
        }
        else
        {
        histProb_woUO->GetYaxis()->SetTitle("Generator mass (dressed level) [GeV]");
        histProb_woUO->GetXaxis()->SetTitle("Generator mass (parton level) [GeV]");
        }
    }

    TLine grid_;
    TLine grid_bin_boundary;
    grid_.SetLineColor(kBlack);
    grid_.SetLineStyle(1);

    if(var=="Pt")
    {
        int nPtBinRec = yaxis1_tvecd->GetNrows() - 1;
        int nPtBinGen = xaxis1_tvecd->GetNrows() - 1;

        int boundarybin_x = 1;
        int countMassBin = 0;
        for( int ii=0; ii<histProb_woUO->GetXaxis()->GetNbins(); ii++ )
        {
            Int_t i_bin = ii+1;
            Double_t binEdge = histProb_woUO->GetXaxis()->GetBinLowEdge(i_bin);

            if(boundarybin_x == i_bin)
            {
                //grid_.DrawLine(binEdge, histProb_woUO->GetYaxis()->GetBinUpEdge(0), binEdge, histProb_woUO->GetYaxis()->GetBinUpEdge(histProb_woUO->GetYaxis()->GetNbins()) );
                boundarybin_x += xaxis1_nbin; // next edge to draw
                countMassBin++;
                if(i_bin == 1) continue;
                grid_.DrawLine(binEdge, histProb_woUO->GetYaxis()->GetBinUpEdge(nPtBinRec * (countMassBin - 2)), binEdge, histProb_woUO->GetYaxis()->GetBinUpEdge(nPtBinRec * (countMassBin)) );
            }
        }

        int boundarybin_y = 1;
        countMassBin = 0;
        for( int ii=0; ii<histProb_woUO->GetYaxis()->GetNbins(); ii++ )
        {
            Int_t i_bin = ii+1;
            Double_t binEdge = histProb_woUO->GetYaxis()->GetBinLowEdge(i_bin);

            if(boundarybin_y == i_bin)
            {
                boundarybin_y += yaxis1_nbin; // next edge to draw
                countMassBin++;
                if(i_bin == 1) continue;
                //grid_.DrawLine(histProb_woUO->GetXaxis()->GetBinUpEdge(0), binEdge, histProb_woUO->GetXaxis()->GetBinUpEdge(histProb_woUO->GetXaxis()->GetNbins()), binEdge);
                grid_.DrawLine(histProb_woUO->GetXaxis()->GetBinUpEdge(nPtBinGen * (countMassBin - 2)), binEdge, histProb_woUO->GetXaxis()->GetBinUpEdge(nPtBinGen * (countMassBin)), binEdge);
            }
        }
    }

    c1->RedrawAxis();
    CMS_lumi(c1, 7, 11);
    c1->cd();

    histProb_woUO->Write();   
    f.Close();

    //c1->SaveAs(isDetector?var + "_responseM_" + year + ".pdf":var + "_FSRresponseM_" + year + ".pdf");
    delete c1;
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
                if((it->first).Contains("iterEM") && !(sysMap_previous[it->first][ith]).Contains("Nominal"))
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

                    if((it->first).Contains("iterEM") && !(sysMap_previous[it->first][ith]).Contains("Nominal"))
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
        if(sysName.Contains("iterEM") && !sysPostfix.Contains("Nominal"))
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
        if(sysName.Contains("iterEM") && !sysPostfix.Contains("Nominal"))
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
                if(sysName.Contains("iterEM") && !sysPostfix.Contains("Nominal"))
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
                if(sysName.Contains("iterEM") && !sysPostfix.Contains("Nominal"))
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
            if(sysName.Contains("iterEM") && !sysPostfix.Contains("Nominal"))
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
            if(sysName.Contains("iterEM") && !sysPostfix.Contains("Nominal"))
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

        if(sysName.Contains("iterEM") && !sysPostfix.Contains("Nominal"))
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

double ISRUnfold::setHistCosmetics(TH1* hist, bool isLogy)
{

    hist->SetTitle("");
    hist->SetStats(false);
    hist->GetXaxis()->SetMoreLogLabels(true);
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(3.2);
    hist->SetLineColor(kBlack);
    hist->GetYaxis()->SetTitle("Events/Bin");
    if(isLogy)
    {
        hist->SetMaximum(hist->GetMaximum() * 1e5);
        hist->SetMinimum(5e-2);
    }
    else
    {
        hist->SetMaximum(hist->GetMaximum() * 2.);
        hist->SetMinimum(0.);
    }

    return hist->GetMaximum();
}

TLegend* ISRUnfold::createLegend(double xStartPos_, double yStartPos_ )
{
    double xStartPos = xStartPos_;
    double yStartPos = yStartPos_;
    double xEndPos = 0.9;
    double yEndPos = 0.85;

    TLegend* leg = new TLegend(xStartPos, yStartPos, xEndPos, yEndPos,"","brNDC");
    leg->SetTextFont(43);
    leg->SetTextSize(70);
    leg->SetFillStyle(0); // transparent
    leg->SetBorderSize(0);

    return leg;
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

TCanvas* ISRUnfold::drawPtBkgRatio(TString filePath)
{

    const TVectorD* ptBinVector   = pt_binning_Rec->GetDistributionBinning(0);
    const TVectorD* massBinVector = pt_binning_Rec->GetDistributionBinning(1);
    const double* ptBinArray = ptBinVector->GetMatrixArray();
    //const double* massBinArray = massBinVector->GetMatrixArray();
    int nPtBin = ptBinVector->GetNrows() - 1;
    int nMassBin = massBinVector->GetNrows() - 1;

    vector<double> newPtBinVector;
    int nTotalPtBins = nPtBin * nMassBin;
    //int countMassBins = 0;
    for(int iMassEdge = 0; iMassEdge < nMassBin; iMassEdge++)
    {
        for(int iPtEdge = 0; iPtEdge < nPtBin + 1; iPtEdge++)
        {
            if(iMassEdge == 0)
            {
                newPtBinVector.push_back(ptBinArray[iPtEdge]);    
            }
            else
            {
                if(iPtEdge == 0) continue;
                double newPtEdge = iMassEdge * ptBinArray[nPtBin] + ptBinArray[iPtEdge];
                newPtBinVector.push_back(newPtEdge);
            }
        }
    } 

    TH1::AddDirectory(kFALSE);
    //TFile* filein = new TFile(filePath);

    setTDRStyle();
    writeExtraText = true;
    extraText  = "Work in progress";
    gStyle->SetLineWidth(5);
    gStyle->SetFrameLineWidth(5);
    gROOT->ForceStyle();
    
    TH1* hDY = new TH1D("hDYPt", "hDYPt", nTotalPtBins, &newPtBinVector[0]);;
    TH1* hDYtautau = new TH1D("hDYtauPt", "hDYtauPt", nTotalPtBins, &newPtBinVector[0]);;
    TH1* hWW = new TH1D("hWW", "hWW", nTotalPtBins, &newPtBinVector[0]);;
    TH1* hWZ = new TH1D("hWZ", "hWZ", nTotalPtBins, &newPtBinVector[0]);;
    TH1* hZZ = new TH1D("hZZ", "hZZ", nTotalPtBins, &newPtBinVector[0]);;
    TH1* hTop = new TH1D("hTop", "hTop", nTotalPtBins, &newPtBinVector[0]);;
    TH1* hQCD = new TH1D("hQCD", "hQCD", nTotalPtBins, &newPtBinVector[0]);;
    TH1* hWjet = new TH1D("hWjet", "hWjet", nTotalPtBins, &newPtBinVector[0]);;
    TH1* hMCtotal = new TH1D("hMCtotal", "hMCtotal", nTotalPtBins, &newPtBinVector[0]);

    //TCanvas* c_out = new TCanvas("detector_level_", "detector_level_", 10000, 2000);
    TCanvas* c_out = new TCanvas("detector_level_", "detector_level_", 3200, 2800);
    c_out->Draw();
    c_out->cd();

    TPad *pad1 = new TPad("pad1","pad1",0,0.0,1,1);
    pad1->SetBottomMargin(topPadBottomMargin * 2.);
    pad1->SetTopMargin(topPadTopMargin);
    pad1->SetLogy();
    pad1->SetTicks(1);
    pad1->Draw();
    pad1->cd();

    TString DYHistName_ = "histo_DYJetsToMuMu";
    if(channel_name == "electron")
    {
        DYHistName_  = "histo_DYJetsToEE";
    }

    for(int iMassBin = 0; iMassBin < nMassBin; iMassBin++)
    {
        TString massString; 
        massString.Form("%d",iMassBin);
        TH1* hDYTemp = getRawHist("Pt_FineCoarse", filePath, "Detector", DYHistName_, "Signal", "pt[UO];mass[UOC"+massString+"]", true, true);
        TH1* hDYtauTemp = getRawHist("Pt_FineCoarse", filePath, "Detector", "histo_DYJetsToTauTau", "DYtautau", "pt[UO];mass[UOC"+massString+"]", true, true);
        TH1* hWWTemp = getRawHist("Pt_FineCoarse", filePath, "Detector", "histo_WW_pythia", "WW", "pt[UO];mass[UOC"+massString+"]", true, true);
        TH1* hWZTemp = getRawHist("Pt_FineCoarse", filePath, "Detector", "histo_WZ_pythia", "WZ", "pt[UO];mass[UOC"+massString+"]", true, true);
        TH1* hZZTemp = getRawHist("Pt_FineCoarse", filePath, "Detector", "histo_ZZ_pythia", "ZZ", "pt[UO];mass[UOC"+massString+"]", true, true);
        TH1* hTopTemp = getRawHist("Pt_FineCoarse", filePath, "Detector", "histo_TTLL_powheg", "TT", "pt[UO];mass[UOC"+massString+"]", true, true);
        TH1* hQCDTemp = getRawHist("Pt_FineCoarse", filePath, "Detector", "histo_QCD", "TT", "pt[UO];mass[UOC"+massString+"]", true, true);
        TH1* hWjetTemp = getRawHist("Pt_FineCoarse", filePath, "Detector", "histo_WJet", "TT", "pt[UO];mass[UOC"+massString+"]", true, true);

        for(int ibin = 1; ibin < hDYTemp->GetNbinsX() + 1; ibin++)
        {
            hDY->SetBinContent(ibin + nPtBin * iMassBin, hDYTemp->GetBinContent(ibin));        
            hDY->SetBinError(ibin + nPtBin * iMassBin, hDYTemp->GetBinError(ibin));        

            hDYtautau->SetBinContent(ibin + nPtBin * iMassBin, hDYtauTemp->GetBinContent(ibin));        
            hDYtautau->SetBinError(ibin + nPtBin * iMassBin, hDYtauTemp->GetBinError(ibin));        

            hWW->SetBinContent(ibin + nPtBin * iMassBin, hWWTemp->GetBinContent(ibin));        
            hWW->SetBinError(ibin + nPtBin * iMassBin, hWWTemp->GetBinError(ibin));        

            hWZ->SetBinContent(ibin + nPtBin * iMassBin, hWZTemp->GetBinContent(ibin));        
            hWZ->SetBinError(ibin + nPtBin * iMassBin, hWZTemp->GetBinError(ibin));        

            hZZ->SetBinContent(ibin + nPtBin * iMassBin, hZZTemp->GetBinContent(ibin));        
            hZZ->SetBinError(ibin + nPtBin * iMassBin, hZZTemp->GetBinError(ibin));        

            hTop->SetBinContent(ibin + nPtBin * iMassBin, hTopTemp->GetBinContent(ibin));        
            hTop->SetBinError(ibin + nPtBin * iMassBin, hTopTemp->GetBinError(ibin));        

            hQCD->SetBinContent(ibin + nPtBin * iMassBin, hQCDTemp->GetBinContent(ibin));        
            hQCD->SetBinError(ibin + nPtBin * iMassBin, hQCDTemp->GetBinError(ibin));        

            hWjet->SetBinContent(ibin + nPtBin * iMassBin, hWjetTemp->GetBinContent(ibin));        
            hWjet->SetBinError(ibin + nPtBin * iMassBin, hWjetTemp->GetBinError(ibin));        
        }
        delete hDYTemp;
        delete hDYtauTemp;
        delete hWWTemp; 
        delete hWZTemp;
        delete hZZTemp;
        delete hTopTemp;
        delete hQCDTemp;
        delete hWjetTemp;
    }

    THStack* hsMC = new THStack("hsMC", "hsMC");
    hMCtotal->Add(hDY);
    TLegend* leg = new TLegend(0.17, 0.5, 0.2, 0.89,"","brNDC");
    leg->SetTextFont(43);
    leg->SetTextSize(100);
    leg->SetFillStyle(0); // transparent
    leg->SetBorderSize(0);
    leg->AddEntry(hDYtautau, "DY#rightarrow#tau#tau", "l");
    leg->AddEntry(hWW, "VV", "l");
    leg->AddEntry(hTop, "t#bar{t}", "l");
    leg->AddEntry(hQCD, "Fake", "l");

    for(int iMassBin = 0; iMassBin < nMassBin; iMassBin++)
    {
        TString massString; 
        massString.Form("%d",iMassBin);
        setTHStack("Pt_FineCoarse", filePath, "Detector", *hsMC, *hMCtotal, *leg, "pt[UO];mass[UOC"+massString+"]", true, "", true, false, true); 
    }
    
    hDYtautau->Divide(hMCtotal);
    hWW->Divide(hMCtotal);
    hWZ->Divide(hMCtotal);
    hZZ->Divide(hMCtotal);
    hTop->Divide(hMCtotal);
    hQCD->Divide(hMCtotal);
    hWjet->Divide(hMCtotal);

    hQCD->Add(hWjet);
    hWW->Add(hWZ);
    hWW->Add(hZZ);

    //double max = setHistCosmetics(hDYtautau, false);
    hDYtautau->GetYaxis()->SetTitle("Background ratio");
    hDYtautau->GetYaxis()->SetTitleOffset(1.3);
    hDYtautau->GetYaxis()->SetTickLength(hDYtautau->GetYaxis()->GetTickLength()*0.3); 
    hDYtautau->SetMaximum(.5);
    hDYtautau->SetMinimum(5e-3);
    hDYtautau->SetFillStyle(0);
    hDYtautau->Draw("hist");
    hDYtautau->SetLineColor(kOrange+2);
    hDYtautau->SetLineWidth(5);

    hWW->SetFillStyle(0);
    hWW->Draw("hist same");
    hWW->SetLineColor(kOrange+1);
    hWW->SetLineWidth(5);

    hTop->SetFillStyle(0);
    hTop->Draw("hist same");
    hTop->SetLineColor(kAzure+1);
    hTop->SetLineWidth(5);

    hQCD->SetFillStyle(0);
    hQCD->Draw("hist same");
    hQCD->SetLineColor(kViolet+6);
    hQCD->SetLineWidth(5);

    hDYtautau->GetXaxis()->ChangeLabel(2,-1,-1,-1,-1,-1,"100");
    hDYtautau->GetXaxis()->ChangeLabel(3,-1,-1,-1,-1,-1,"100");
    hDYtautau->GetXaxis()->ChangeLabel(4,-1,-1,-1,-1,-1,"100");
    hDYtautau->GetXaxis()->ChangeLabel(5,-1,-1,-1,-1,-1,"100");
    hDYtautau->GetXaxis()->ChangeLabel(6,-1,-1,-1,-1,-1,"100");

    leg->Draw();
    int iPeriod_ = 4;
    if(year == 2017)
        iPeriod_ = 5;
    if(year == 2018)
        iPeriod_ = 6;
    CMS_lumi(pad1, iPeriod_, 0);

    TLine massEdgeLine; 
    massEdgeLine.SetLineStyle(1); 
    for(int iMassBin = 0; iMassBin < nMassBin; iMassBin++) 
    {
        for(int ibin = 1; ibin < nPtBin + 1; ibin++)
        {
            if(iMassBin != 0 && ibin == 1) 
            {
                massEdgeLine.DrawLine(hDYtautau->GetXaxis()->GetBinLowEdge(ibin + nPtBin * iMassBin), 5e-3, hDYtautau->GetXaxis()->GetBinLowEdge(ibin + nPtBin * iMassBin), .5);
            }
        }
    }
    c_out->SaveAs(output_baseDir+"BkgRatio_CombinedPt"+".pdf");
    return c_out;
}

TCanvas* ISRUnfold::drawPtDistributions(TString filePath)
{
    const TVectorD* ptBinVector    = pt_binning_Rec->GetDistributionBinning(0);
    const TVectorD* massBinVector = pt_binning_Rec->GetDistributionBinning(1);
    const double* ptBinArray = ptBinVector->GetMatrixArray();
    //const double* massBinArray = massBinVector->GetMatrixArray();
    int nPtBin = ptBinVector->GetNrows() - 1;
    int nMassBin = massBinVector->GetNrows() - 1;

    vector<double> newPtBinVector;
    int nTotalPtBins = nPtBin * nMassBin;
    //int countMassBins = 0;
    for(int iMassEdge = 0; iMassEdge < nMassBin; iMassEdge++)
    {
        for(int iPtEdge = 0; iPtEdge < nPtBin + 1; iPtEdge++)
        {
            if(iMassEdge == 0)
            {
                newPtBinVector.push_back(ptBinArray[iPtEdge]);    
            }
            else
            {
                if(iPtEdge == 0) continue;
                double newPtEdge = iMassEdge * ptBinArray[nPtBin] + ptBinArray[iPtEdge];
                newPtBinVector.push_back(newPtEdge);
            }
        }
    } 

    TH1::AddDirectory(kFALSE);
    //TFile* filein = new TFile(filePath);

    setTDRStyle();
    writeExtraText = true;
    extraText  = "Work in progress";
    gStyle->SetLineWidth(1);
    gStyle->SetFrameLineWidth(5);
    gROOT->ForceStyle();
    
    TH1* hData = new TH1D("hDataPt", "hDataPt", nTotalPtBins, &newPtBinVector[0]);
    TH1* hDY = new TH1D("hDYPt", "hDYPt", nTotalPtBins, &newPtBinVector[0]);;
    TH1* hMCtotal = new TH1D("hMCtotal", "hMCtotal", nTotalPtBins, &newPtBinVector[0]);

    //TCanvas* c_out = new TCanvas("detector_level_", "detector_level_", 10000, 2800);
    TCanvas* c_out = new TCanvas("detector_level_", "detector_level_", 3200, 2800);
    c_out->Draw();
    c_out->cd();

    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetBottomMargin(topPadBottomMargin);
    pad1->SetTopMargin(topPadTopMargin);
    pad1->SetTicks(1);
    pad1->SetLogy();
    pad1->Draw();
    pad1->cd();

    TString dataHistName_ = "histo_DoubleMuon";
    TString DYHistName_ = "histo_DYJetsToMuMu";
    if(channel_name == "electron")
    {
        dataHistName_ = "histo_DoubleEG";
        if(year==2018)
            dataHistName_ = "histo_EGamma";
        DYHistName_  = "histo_DYJetsToEE";
    }

    for(int iMassBin = 0; iMassBin < nMassBin; iMassBin++)
    {
        TString massString; 
        massString.Form("%d",iMassBin);
        TH1* hDataTemp = getRawHist("Pt_FineCoarse", filePath, "Detector", dataHistName_, "Data", "pt[UO];mass[UOC"+massString+"]", true, true);
        TH1* hDYTemp = getRawHist("Pt_FineCoarse", filePath, "Detector", DYHistName_, "Signal", "pt[UO];mass[UOC"+massString+"]", true, true);

        for(int ibin = 1; ibin < hDataTemp->GetNbinsX() + 1; ibin++)
        {
            hData->SetBinContent(ibin + nPtBin * iMassBin, hDataTemp->GetBinContent(ibin));        
            hDY->SetBinContent(ibin + nPtBin * iMassBin, hDYTemp->GetBinContent(ibin));        

            hData->SetBinError(ibin + nPtBin * iMassBin, hDataTemp->GetBinError(ibin));        
            hDY->SetBinError(ibin + nPtBin * iMassBin, hDYTemp->GetBinError(ibin));        
        }

        delete hDataTemp;
        delete hDYTemp;
    }

    THStack* hsMC = new THStack("hsMC", "hsMC");
    hMCtotal->Add(hDY);
    TLegend* leg = new TLegend(0.17, 0.5, 0.2, 0.89,"","brNDC");
    leg->SetTextFont(43);
    leg->SetTextSize(70);
    leg->SetFillStyle(0); // transparent
    leg->SetBorderSize(0);
    leg->AddEntry(hData, "Data", "pl");
    leg->AddEntry(hDY, "Drell-Yan", "F");
    bool updateLegend = true;
    for(int iMassBin = 0; iMassBin < nMassBin; iMassBin++)
    {
        TString massString; 
        massString.Form("%d",iMassBin);
        setTHStack("Pt_FineCoarse", filePath, "Detector", *hsMC, *hMCtotal, *leg, "pt[UO];mass[UOC"+massString+"]", true, "", true, updateLegend, true); 
        updateLegend = false;
    }
    hsMC->Add(hDY);

    double max = setHistCosmetics(hData, true);
    hData->GetYaxis()->SetTickLength(hData->GetYaxis()->GetTickLength()*0.3); 
    hData->Draw("p9e");
    hDY->SetFillColor(kOrange);
    hsMC->Draw("hist same");
    leg->Draw();

    TLine massEdgeLine; 
    massEdgeLine.SetLineStyle(1); 
    for(int iMassBin = 0; iMassBin < nMassBin; iMassBin++) 
    {
        for(int ibin = 1; ibin < nPtBin + 1; ibin++)
        {
            if(iMassBin != 0 && ibin == 1) 
            {
                massEdgeLine.DrawLine(hData->GetXaxis()->GetBinLowEdge(ibin + nPtBin * iMassBin), 5e-2, hData->GetXaxis()->GetBinLowEdge(ibin + nPtBin * iMassBin), max);
            }
        }
    }

    TH1 *last = (TH1*)hsMC->GetStack()->Last();
    for(int ibin = 1; ibin < hData->GetNbinsX()+1; ibin++)
    {
        massEdgeLine.SetLineStyle(1);
        massEdgeLine.SetLineColor(kWhite);
        if(ibin == 1 || ibin % nPtBin == 1)
        {
        massEdgeLine.DrawLine(hData->GetXaxis()->GetBinLowEdge(ibin+1), 5e-2, hData->GetXaxis()->GetBinLowEdge(ibin+1), last->GetBinContent(ibin));
        }
        else if(ibin % nPtBin == 0)
        {
        massEdgeLine.DrawLine(hData->GetXaxis()->GetBinLowEdge(ibin), 5e-2, hData->GetXaxis()->GetBinLowEdge(ibin), last->GetBinContent(ibin));
        }
        else
        {
        massEdgeLine.DrawLine(hData->GetXaxis()->GetBinLowEdge(ibin), 5e-2, hData->GetXaxis()->GetBinLowEdge(ibin), last->GetBinContent(ibin));
        massEdgeLine.DrawLine(hData->GetXaxis()->GetBinLowEdge(ibin+1), 5e-2, hData->GetXaxis()->GetBinLowEdge(ibin+1), last->GetBinContent(ibin));
        }
    }
    hData->GetXaxis()->Draw();
    hData->GetXaxis()->SetLabelSize(0);
    hData->Draw("p9e same");

    int iPeriod_ = 4;
    if(year == 2017)
        iPeriod_ = 5;
    if(year == 2018)
        iPeriod_ = 6;
    CMS_lumi(pad1, iPeriod_, 0);

    writeCutInfo(pad1, "Pt", 0, 0.22, 0.8, true);
    writeCutInfo(pad1, "Pt", 1, 0.375, 0.8, false);
    writeCutInfo(pad1, "Pt", 2, 0.530, 0.8, false);
    writeCutInfo(pad1, "Pt", 3, 0.685, 0.8, false);
    writeCutInfo(pad1, "Pt", 4, 0.840, 0.8, false);

    c_out->cd(); 
    
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.35);
    pad2->SetTicks(1);
    pad2->Draw();
    pad2->cd();

    TH1* hRatio = (TH1*) hData->Clone("hRatio");
    hRatio->GetXaxis()->SetLabelSize(70);
    hRatio->SetStats(false);
    hRatio->Divide(hMCtotal);

    hRatio->Draw("p9e");
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(3.2);
    hRatio->SetLineColor(kBlack);
    hRatio->GetYaxis()->SetTitle("Data/MC");

    hRatio->SetMinimum(0.65);
    hRatio->SetMaximum(1.35);
    hRatio->GetYaxis()->SetDecimals();

    massEdgeLine.SetLineColor(kBlack);
    for(int iMassBin = 0; iMassBin < nMassBin; iMassBin++) 
    {
        for(int ibin = 1; ibin < nPtBin + 1; ibin++)
        {
            if(iMassBin != 0 && ibin == 1) 
            {
                massEdgeLine.DrawLine(hData->GetXaxis()->GetBinLowEdge(ibin + nPtBin * iMassBin), 0.65, hData->GetXaxis()->GetBinLowEdge(ibin + nPtBin * iMassBin), 1.35);
            }
        }
    }

    hRatio->GetXaxis()->ChangeLabel(2,-1,-1,-1,-1,-1,"100");
    hRatio->GetXaxis()->ChangeLabel(3,-1,-1,-1,-1,-1,"100");
    hRatio->GetXaxis()->ChangeLabel(4,-1,-1,-1,-1,-1,"100");
    hRatio->GetXaxis()->ChangeLabel(5,-1,-1,-1,-1,-1,"100");
    hRatio->GetXaxis()->ChangeLabel(6,-1,-1,-1,-1,-1,"100");

    setXaxisTitle(hRatio, "Pt", true);

    TLine* l_ = new TLine(hRatio->GetXaxis()->GetXmin(), 1, hRatio->GetXaxis()->GetXmax(), 1);
    l_->SetLineColor(kRed);
    l_->Draw("same");
    l_->SetLineStyle(1);

    c_out->cd();
    c_out->SaveAs(output_baseDir+"Det_CombinedPt"+".pdf");
    return c_out;
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


// FIXME Too complecate, simpify this function!
// Draw detector distributions using input root file
TCanvas& ISRUnfold::drawFoldedHists(TString var, TString filePath, TString dirName, TString steering, bool useAxis, TString sysName, TString outName, int nthMassBin, bool divBinWidth, TString sysFilePath, bool isBkgSubData)
{
    // If steering == "", then usual TH1 histogram
    // If seering != "", TH1 from TUnfold
    TH1::AddDirectory(kFALSE);
    TFile* filein = new TFile(filePath);

    double meanDipt = 0.;
    double meanDipt_bkgsub = 0.;
    double histMin = 5e-2;

    setTDRStyle();
    writeExtraText = true;
    extraText  = "Work in progress";
    gStyle->SetLineWidth(1);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetImageScaling(3.);
    gROOT->ForceStyle();

    // For nominal histogram
    TH1* hData = NULL;
    TH1* hDY = NULL;
    TH1* hMCtotal = NULL;
    TH1* hRatio = NULL;

    TString dataHistName_ = "histo_DoubleMuon";
    TString DYHistName_ = "histo_DYJetsToMuMu";
    if(channel_name == "electron")
    {
        dataHistName_ = "histo_DoubleEG";
        if(year==2018)
            dataHistName_ = "histo_EGamma";
        DYHistName_  = "histo_DYJetsToEE";
    }

    hData = getRawHist(var, filePath, dirName, dataHistName_, "Data", steering, useAxis, divBinWidth);
    TGraphErrors* gData = histToTGraphError(hData);
    if(sysName == "ZptCorr")
        hDY = getRawHist(var, sysFilePath, dirName, DYHistName_, "Signal", steering, useAxis, divBinWidth);
    else
        hDY = getRawHist(var, filePath, dirName, DYHistName_, "Signal", steering, useAxis, divBinWidth);

    hMCtotal = (TH1*) hDY->Clone("hMCtotal");

    // Create canvas
    TCanvas* c_out = new TCanvas("detector_level_"+var, "detector_level_"+var, 3200, 2800);
    c_out->Draw();
    c_out->cd();

    // Create top pad
    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetBottomMargin(topPadBottomMargin);
    pad1->SetTopMargin(topPadTopMargin);
    pad1->SetTicks(1);
    pad1->SetLogy();
    if(var.Contains("Mass"))
        pad1->SetLogx();
    pad1->Draw();
    pad1->cd();

    THStack* hsMC = new THStack("hsMC", "hsMC");
    TH1* hDataBkgSubtracted = (TH1*) hData->Clone("DataBKGsubtracted");
    double xStartPos = 0.8;
    if(isBkgSubData) xStartPos = 0.75;
    TLegend* leg = createLegend(xStartPos, 0.5);
    if(isBkgSubData)
        leg->AddEntry(hDataBkgSubtracted, "Data (bkg. sub.)", "pl");
    else
        leg->AddEntry(hData, "Data", "pl");
    leg->AddEntry(hDY, "Drell-Yan", "F");

    if(sysName == "ZptCorr")
        setTHStack(var, sysFilePath, dirName, *hsMC, *hMCtotal, *leg, steering, useAxis, "", divBinWidth, !isBkgSubData);
    else
        setTHStack(var, filePath, dirName, *hsMC, *hMCtotal, *leg, steering, useAxis, "", divBinWidth, !isBkgSubData);

    hDY->SetLineColor(kWhite);
    hsMC->Add(hDY);
    hDY->SetFillColor(kOrange);

    // Get average transeverse momentum values
    TH1* hMCtotalDYSubtracted = (TH1*) hMCtotal->Clone("MCtotalDYsubtracted");
    hMCtotalDYSubtracted->Add(hDY, -1);
    hDataBkgSubtracted->Add(hMCtotalDYSubtracted, -1);
    meanDipt = hData->GetMean();
    meanDipt_bkgsub = hDataBkgSubtracted->GetMean();
    //double binnedMean = getBinnedMean(hData); // To check how average value changes with binned histogram

    if(isBkgSubData)
    {
        hRatio = (TH1*) hDataBkgSubtracted->Clone("hRatio");
        setHistCosmetics(hDataBkgSubtracted, false);
        hDataBkgSubtracted->GetXaxis()->SetLabelSize(0);
        hDataBkgSubtracted->Draw("p9e");

        hDY->Draw("hist same");
        hDataBkgSubtracted->Draw("p9e same");
        pad1->RedrawAxis();

        TLine meanLine;
        meanLine.SetLineStyle(1);
        meanLine.SetLineColor(kRed);
        meanLine.DrawLine(meanDipt_bkgsub, histMin, meanDipt_bkgsub, hDataBkgSubtracted->GetMaximum());

        TLine meanLine_;
        meanLine_.SetLineStyle(2);
        meanLine_.SetLineColor(kBlack);
        meanLine_.DrawLine(meanDipt, histMin, meanDipt, hDataBkgSubtracted->GetMaximum());
    }
    else
    {
        // histogram to tgraph 
        hRatio = (TH1*) hData->Clone("hRatio");
        setHistCosmetics(hData, true);
        hData->GetXaxis()->SetLabelSize(0);
        hData->Draw("hist p");
        hData->SetMarkerSize(4);
        hData->SetLineWidth(1);

        hsMC->Draw("hist same");
        //hData->Draw("hist same");
        gData->Draw("pe same");
        gData->SetMarkerSize(4);
        gData->SetLineWidth(1);
        pad1->RedrawAxis();
    }

    //TLatex meanDipt_;
    //TLatex meanDipt_bkgsub_;
    //meanDipt_.SetTextFont(43);
    //meanDipt_.SetTextSize(40);
    //meanDipt_bkgsub_.SetTextFont(43);
    //meanDipt_bkgsub_.SetTextSize(40);

    //if(var.Contains("Pt"))
    //{
    //    TString itos, itos_;
    //    itos.Form ("%.2f", meanDipt);
    //    itos_.Form ("%.2f", binnedMean);
    //    meanDipt_.DrawLatexNDC(0.2, 0.6, "avg. p_{T}^{Data}: "+itos+"("+itos_+") GeV");
    //    itos.Form ("%.2f", meanDipt_bkgsub);
    //    meanDipt_bkgsub_.DrawLatexNDC(0.2, 0.6-0.07, "avg. p_{T}^{Data-Bkg}: "+itos);
    //}

    //getSmearedChi2(TString var, TString filePath, TString dirName, TString steering, bool useAxis, bool divBinWidth)
    //TLatex smearedChi2;
    //smearedChi2.SetTextFont(43);
    //smearedChi2.SetTextSize(100);
    //TString chi2_;
    //chi2_.Form("%.2f", getSmearedChi2(var, filePath, dirName, steering,  useAxis, false));
    //smearedChi2.DrawLatexNDC(0.2, 0.6, "#chi^{2}: " + chi2_);

    leg->Draw();

    int iPeriod_ = 4;
    if(year == 2017)
        iPeriod_ = 5;
    if(year == 2018)
        iPeriod_ = 6;
    CMS_lumi(pad1, iPeriod_, 11);
    writeCutInfo(pad1, var, nthMassBin);

    TLine massEdgeLine;
    //if(var.Contains("Mass"))
    {
        //for(unsigned int i = 0; i < massBinEdges.size(); i++)
        //{
        //    if(i==0) continue;
        //    if(i==massBinEdges.size()-1) continue;
        //    massEdgeLine.SetLineStyle(2);
        //    massEdgeLine.DrawLine(massBinEdges[i], hData->GetMinimum(), massBinEdges[i], hData->GetMaximum());
        //}

        TH1 *last = (TH1*)hsMC->GetStack()->Last();
        for(int ibin = 1; ibin < hData->GetNbinsX()+1; ibin++)
        {
            massEdgeLine.SetLineStyle(1);
            massEdgeLine.SetLineColor(kWhite);
            massEdgeLine.SetLineWidth(1);
            if(ibin == 1)
            {
                massEdgeLine.DrawLine(hData->GetXaxis()->GetBinLowEdge(ibin+1), histMin, hData->GetXaxis()->GetBinLowEdge(ibin+1), last->GetBinContent(ibin));
            }
            else
            {
                massEdgeLine.DrawLine(hData->GetXaxis()->GetBinLowEdge(ibin), histMin, hData->GetXaxis()->GetBinLowEdge(ibin), last->GetBinContent(ibin));
                massEdgeLine.DrawLine(hData->GetXaxis()->GetBinLowEdge(ibin+1), histMin, hData->GetXaxis()->GetBinLowEdge(ibin+1), last->GetBinContent(ibin));
            }
        }
    }
    c_out->cd();

    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    pad2->SetTopMargin(bottomPadTopMargin);
    pad2->SetBottomMargin(bottomPadBottomMargin);
    if(var.Contains("Mass"))
    {
        pad2->SetLogx();
        hRatio->GetXaxis()->SetMoreLogLabels(true);
    }
    pad2->SetTicks(1);
    //pad2->SetGridy(1);
    pad2->Draw();
    pad2->cd();

    hRatio->SetStats(false);
    if(isBkgSubData)
    {
        hRatio->Divide(hDY);
    }
    else
    {
        hRatio->Divide(hMCtotal);
    }

    TGraphErrors* gRatio = histToTGraphError(hRatio);
    hRatio->Draw("hist p");
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(4);
    hRatio->SetLineWidth(1);
    hRatio->SetLineColor(kBlack);
    hRatio->GetYaxis()->SetTitle("Data/MC");

    hRatio->SetMinimum(0.65);
    hRatio->SetMaximum(1.35);
    hRatio->GetYaxis()->SetDecimals();

    setXaxisTitle(hRatio, var, useAxis);

    TH1* sysBand_ratio_forData = NULL;
    TH1* sysBand_ratio_forMC = NULL;
    TLegend* leg_sys = NULL;
    if(sysName != "")
    {
        leg_sys = new TLegend(0.7, 0.8, 0.95, 0.9,"","brNDC");

        sysBand_ratio_forData = getDetectorSystematicBand(var, filePath, dirName, steering, useAxis,  sysName, hData, hDY, hMCtotal, hRatio, divBinWidth, true, false, sysFilePath, nthMassBin);
        sysBand_ratio_forData->SetFillColorAlpha(kBlack,0.8);

        sysBand_ratio_forMC = getDetectorSystematicBand(var, filePath, dirName, steering, useAxis,  sysName, hData, hDY, hMCtotal, hRatio, divBinWidth, true, true, sysFilePath, nthMassBin);

        TH1* upBound = cloneEmptyHist(sysBand_ratio_forMC, "up");
        TH1* downBound = cloneEmptyHist(sysBand_ratio_forMC, "down");

        if(var.Contains("Pt"))
        {
            sysRelPtHist_detectorMC[nthMassBin][sysName] = sysBand_ratio_forMC;
            sysRelPtHist_detectorData[nthMassBin][sysName]= sysBand_ratio_forData;
        }
        else
        {
            sysRelMassHist_detectorMC[sysName] = sysBand_ratio_forMC;
            sysRelMassHist_detectorData[sysName] = sysBand_ratio_forData;
        }


        for(int i = 1; i < sysBand_ratio_forMC->GetNbinsX()+1; i++)
        {
            upBound->SetBinContent(i, 1. + sysBand_ratio_forMC->GetBinError(i)); 
            upBound->SetBinError(i, 0);
            downBound->SetBinContent(i, 1. - sysBand_ratio_forMC->GetBinError(i)); 
            downBound->SetBinError(i, 0);
        }
        upBound->Draw("same");
        upBound->SetMarkerSize(0);
        upBound->SetLineColorAlpha(kRed, 0.2);

        downBound->Draw("same");
        downBound->SetMarkerSize(0);
        downBound->SetLineColorAlpha(kRed, 0.2);

        gStyle->SetHatchesSpacing(1);
        sysBand_ratio_forData->SetFillStyle(3004);
        sysBand_ratio_forData->SetMarkerSize(0.);
        sysBand_ratio_forData->Draw("E2 same");
        sysBand_ratio_forMC->SetFillStyle(3004);
        //sysBand_ratio_forMC->SetFillColorAlpha(kRed, 0.1);
        sysBand_ratio_forMC->SetFillColor(kRed-7);
        sysBand_ratio_forMC->SetMarkerSize(0.);
        sysBand_ratio_forMC->SetLineColor(kRed);
        sysBand_ratio_forMC->SetLineWidth(1);
        sysBand_ratio_forMC->Draw("E2 same");

        for(int ibin = 2; ibin < sysBand_ratio_forMC->GetNbinsX()+1; ibin++)
        {
            massEdgeLine.SetLineStyle(2);
            massEdgeLine.SetLineColor(kWhite);
            massEdgeLine.DrawLine(hData->GetXaxis()->GetBinLowEdge(ibin), downBound->GetBinContent(ibin), hData->GetXaxis()->GetBinLowEdge(ibin), upBound->GetBinContent(ibin));
        }

        leg_sys->SetTextFont(43);
        leg_sys->SetTextSize(100);
        leg_sys->SetFillStyle(0); // transparent
        leg_sys->SetBorderSize(0);
        TString sysName_ = getSysNameToShow(sysName);
        leg_sys->AddEntry(sysBand_ratio_forMC, sysName_, "F");
        leg_sys->Draw();

    }

    TLine* l_ = new TLine(hRatio->GetXaxis()->GetXmin(), 1, hRatio->GetXaxis()->GetXmax(), 1);
    l_->SetLineColor(kRed-7);
    l_->Draw("same");
    l_->SetLineStyle(1);

    gRatio->Draw("p same");
    gRatio->SetLineWidth(1);
    gRatio->SetMarkerSize(4);

    // Save canvas
    c_out->cd();
    c_out->SaveAs(outName!=""?output_baseDir+outName+var+sysName+".png":output_baseDir+"detector_"+var+sysName+".png");

    //delete filein;
    //delete gRatio;
    //delete gData;
    return *c_out;
}

TString ISRUnfold::getSysNameToShow(TString sysName)
{
    if(sysName == "ID")
        return "ID SF";
    else if(sysName == "Reco")
        return "Reco SF";
    else if(sysName == "TRG")
        return "Trigger SF";
    else if(sysName == "AlphaS")
        return "#alpha_{S}";
    else if(sysName == "Scale")
        return "#mu_{F}, #mu_{R}";
    else if(sysName == "iterEM")
        return "Unfolding";
    else if(sysName == "LepScale")
        return "Lepton momentum scale";
    else if(sysName == "LepRes")
        return "Lepton momentum resolution";
    else if(sysName == "totalSys")
        return "sys. unc.";
    else return sysName;
}

TCanvas* ISRUnfold::drawAcceptVarHists(TString var, TString steering, bool useAxis, TString sysName, TString outName, int nthMassBin, bool isAccept)
{
    double histMin = 0.;

    setTDRStyle();
    writeExtraText = true;
    extraText  = "Simulation";
    gStyle->SetLineWidth(globalLinedWidth);
    gStyle->SetFrameLineWidth(globalFrameWidth);
    gROOT->ForceStyle();

    TH1::AddDirectory(kFALSE);

    // For nominal histogram
    TH1* hData = NULL;

    if(var.Contains("Pt"))
    {
        hData = pt_binning_Gen->ExtractHistogram("hData", hAcceptanceFractionPt, 0, useAxis, steering);
    }
    else
    {
        hData = mass_binning_Gen->ExtractHistogram("hData", hAcceptanceFractionMass, 0, useAxis, steering);
    }

    // Create canvas
    TCanvas* c_out ;
    c_out = new TCanvas("SysUnfolded_level_"+var, "unfolded_level_"+var, 3200, 2800);

    c_out->Draw();
    c_out->cd();

    TPad *pad1 = new TPad("pad1","pad1",0,0.,1,1);
    pad1->SetBottomMargin(0.12);
    pad1->SetTopMargin(0.1);
    pad1->SetTicks(1);
    //pad1->SetLogy();
    if(var.Contains("Mass"))
    {
        pad1->SetLogx();
        hData->GetXaxis()->SetMoreLogLabels(true);
    }
    pad1->Draw();
    pad1->cd();

    TGraphErrors* gData = histToTGraphError(hData, false);
    hData->SetTitle("");
    hData->SetStats(false);
    hData->Draw("hist p");
    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(4);
    hData->SetLineColor(kBlack);
    hData->GetYaxis()->SetTitle("Fraction of Events");
    hData->SetMaximum(1.);
    hData->SetMinimum(histMin);
    hData->GetXaxis()->SetTitleOffset(1.2);
    setXaxisTitle(hData, var, useAxis);

    gStyle->SetHatchesSpacing(1);
    gData->SetFillStyle(3354); 
    gData->SetFillColor(kBlack); 
    gData->SetMarkerSize(0);
    gData->SetLineWidth(1);
    gData->SetLineColor(kBlack);
    gData->Draw("E2 same");

    TLegend* leg = new TLegend(0.55, 0.75, 0.9, 0.95,"","brNDC");
    leg->SetTextFont(43);
    leg->SetTextSize(100);
    leg->SetFillStyle(0); // transparent
    leg->SetBorderSize(0);
    if(isAccept)
        leg->AddEntry(hData, "Acceptance", "pl");
    else
        leg->AddEntry(hData, "Selection efficiency", "pl");

    pad1->RedrawAxis();

    TH1* hAcceptSys = cloneEmptyHist(hData, "acceptanceSystematic");
    for(int ibin = 1; ibin < hAcceptSys->GetNbinsX()+1; ibin++)
    {
        hAcceptSys->SetBinContent(ibin, hData->GetBinContent(ibin));
        hAcceptSys->SetBinError(ibin, 0.);
    }
    // Draw systematic results
    // PDF, Scale, AlphaS
    std::map<TString, std::vector<TString>>::iterator it = sysMap.begin();
    int sysSize = sysMap[sysName].size();
    for(int ith = 0; ith < sysSize; ith++)
    {
        TH1* htemp = NULL;
        if(var.Contains("Pt"))
        {
            htemp = pt_binning_Gen->ExtractHistogram("hData", hSysAcceptanceFractionPt[sysName][sysMap[sysName].at(ith)], 0, useAxis, steering);
        }
        else
        {
            htemp = mass_binning_Gen->ExtractHistogram("hData", hSysAcceptanceFractionMass[sysName][sysMap[sysName].at(ith)], 0, useAxis, steering);
        }

        for(int ibin = 1; ibin < htemp->GetNbinsX()+1; ibin++)
        {
            double tempBinError = fabs(hData->GetBinContent(ibin)-htemp->GetBinContent(ibin));
            // FIXME For PDF use RMS as uncertainty
            if(tempBinError > hAcceptSys->GetBinError(ibin))
            {
                hAcceptSys->SetBinError(ibin, tempBinError);
            }
        }
    }

    hAcceptSys->SetFillColorAlpha(kRed, 0.5);    
    hAcceptSys->SetLineColor(kRed);
    hAcceptSys->SetMarkerSize(0.); 
    hAcceptSys->Draw("E2 same");

    leg->Draw();

    int iPeriod_ = 4;
    if(year == 2017)
        iPeriod_ = 5;
    if(year == 2018)
        iPeriod_ = 6;
    CMS_lumi(pad1, iPeriod_, 0);
    writeCutInfo(pad1, var, nthMassBin, 0.55, 0.2);

    c_out->cd();
    c_out->SaveAs(outName!=""?output_baseDir+outName+var+sysName+"_var.pdf":output_baseDir+"unfolded_"+var+sysName+".pdf");

    return c_out;
}

TCanvas* ISRUnfold::drawUnfoldedVarHists(TString var, TString steering, bool useAxis, TString sysName, TString outName, int nthMassBin, bool divBinWidth)
{

    double histMin = 5e-1;

    setTDRStyle();
    writeExtraText = true;
    extraText  = "Work in progress";
    gStyle->SetLineWidth(3.);
    gStyle->SetFrameLineWidth(3.);
    gROOT->ForceStyle();

    TH1::AddDirectory(kFALSE);

    // For nominal histogram
    TH1* hData = NULL;
    TH1* hRatio = NULL;

    if(var.Contains("Pt"))
    {
        hData = nomPtUnfold->GetOutput("hUnfoldedPt",0,0,steering,useAxis);
    }
    else
    {
        hData = nomMassUnfold->GetOutput("hUnfoldedMass",0,0,steering,useAxis);
    }
    if(divBinWidth)
    {
        divideByBinWidth(hData, false);
    }
    hRatio = (TH1*) hData->Clone("hRatio");

    // Create canvas
    TCanvas* c_out ;
    c_out = new TCanvas("SysUnfolded_level_"+var, "unfolded_level_"+var, 3200, 2800);

    c_out->Draw();
    c_out->cd();

    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetBottomMargin(0.01);
    pad1->SetTopMargin(0.1);
    pad1->SetTicks(1);
    //pad1->SetLogy();
    if(var.Contains("Mass"))
        pad1->SetLogx();
    pad1->Draw();
    pad1->cd();

    setHistCosmetics(hData, false);
    hData->Draw("p9e");

    TLegend* leg = new TLegend(0.55, 0.7, 0.9, 0.92,"","brNDC");
    //leg->SetNColumns(2);
    leg->SetTextFont(43);
    leg->SetTextSize(70);
    leg->SetFillStyle(0); // transparent
    leg->SetBorderSize(0);
    if(outName.Contains("Closure"))
    {
        leg->AddEntry(hData, "Unfolded MC", "pl");
    }
    else
    {
        leg->AddEntry(hData, "Unfolded data", "pl");
    }

    hData->Draw("p9e same");
    pad1->RedrawAxis();

    // Draw systematic results
    std::map<TString, std::vector<TString>>::iterator it = sysMap.begin();
    int sysSize = sysMap[sysName].size();
    vector<TH1*> v_hData_temp ;
    vector<TH1*> v_hRatio_temp ;
    for(int ith = 0; ith < sysSize; ith++)
    {
        if(sysName.Contains("iterEM") && sysMap[sysName].at(ith) != "Nominal")
        {
            if(var.Contains("Pt"))
            {
                v_hData_temp.push_back(iterEMPtUnfold->GetOutput("hUnfoldedPt",0,0,steering,useAxis));
                v_hRatio_temp.push_back(iterEMPtUnfold->GetOutput("hRatioPt",0,0,steering,useAxis));
            }
            else
            {
                v_hData_temp.push_back(iterEMMassUnfold->GetOutput("hUnfoldedMass_",0,0,steering,useAxis));
                v_hRatio_temp.push_back(iterEMMassUnfold->GetOutput("hRatioMass_",0,0,steering,useAxis));
            }
        }
        else
        {
            if(var.Contains("Pt"))
            {
                v_hData_temp.push_back(sysPtUnfold[sysName][sysMap[sysName].at(ith)]->GetOutput("hUnfoldedPt",0,0,steering,useAxis));
                v_hRatio_temp.push_back(sysPtUnfold[sysName][sysMap[sysName].at(ith)]->GetOutput("hRatioPt",0,0,steering,useAxis));
            }
            else
            {
                v_hData_temp.push_back(sysMassUnfold[sysName][sysMap[sysName].at(ith)]->GetOutput("hUnfoldedMass_",0,0,steering,useAxis));
                v_hRatio_temp.push_back(sysMassUnfold[sysName][sysMap[sysName].at(ith)]->GetOutput("hRatioMass_",0,0,steering,useAxis));
            }
        }

        divideByBinWidth(v_hData_temp.at(ith), false);
        divideByBinWidth(v_hRatio_temp.at(ith), false);
        v_hData_temp.at(ith)->SetLineStyle(1);
        v_hData_temp.at(ith)->SetLineWidth(5);
        v_hData_temp.at(ith)->SetLineColor(2 + ith);
        v_hData_temp.at(ith)->Draw("hist same");
        //delete hData_temp;
    }

    leg->Draw();

    int iPeriod_ = 4;
    if(year == 2017)
        iPeriod_ = 5;
    if(year == 2018)
        iPeriod_ = 6;
    CMS_lumi(pad1, iPeriod_, 11);
    // writeCutInfo(pad, var, nthMassBin);
    if(!outName.Contains("Closure"))
        writeCutInfo(pad1, var, nthMassBin);

    // getUnfoldedChi2(TString var, TString steering, bool useAxis, bool divBinWidth)
    TLatex unfoldedChi2;
    unfoldedChi2.SetTextFont(43);
    unfoldedChi2.SetTextSize(50);
    TString chi2_;
    chi2_.Form("%.2f", getUnfoldedChi2(var, steering, useAxis, divBinWidth));
    unfoldedChi2.DrawLatexNDC(0.2, 0.6, "#chi^{2}: " + chi2_);

    c_out->cd();

    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.35);
    if(var=="Mass")
    {
        pad2->SetLogx();
        hRatio->GetXaxis()->SetMoreLogLabels(true);
    }
    pad2->SetTicks(1);
    pad2->SetGridy(1);
    pad2->Draw();
    pad2->cd();

    hRatio->Divide(hData);
    hRatio->Draw("hist");
    hRatio->SetMinimum(0.9);
    hRatio->SetMaximum(1.1);

    for(int ith = 0; ith < sysSize; ith++)
    {
        v_hRatio_temp.at(ith)->Divide(hData);
        v_hRatio_temp.at(ith)->SetLineStyle(1);
        v_hRatio_temp.at(ith)->SetLineWidth(5);
        v_hRatio_temp.at(ith)->SetLineColor(2 + ith);
        v_hRatio_temp.at(ith)->Draw("hist same");
    }

    TLine* l_ = new TLine(hRatio->GetXaxis()->GetXmin(),1,hRatio->GetXaxis()->GetXmax(),1);
    l_->SetLineColor(kRed);
    l_->Draw("same");
    l_->SetLineStyle(3);

    // Save canvas
    c_out->cd();
    c_out->SaveAs(outName!=""?output_baseDir+outName+var+sysName+"_var.pdf":output_baseDir+"unfolded_"+var+sysName+".pdf");

    for(int ith = 0; ith < sysSize; ith++)
    {
        delete v_hData_temp.at(ith);
        delete v_hRatio_temp.at(ith);
    }
    return c_out;
}

TCanvas* ISRUnfold::drawUnfoldedHists(TString var, TString steering, bool useAxis, TString sysName, TString outName, int nthMassBin, bool divBinWidth, bool isType3Closure)
{
    // If steering == "", then usual TH1 histogram
    // If seering != "", TH1 from TUnfold
    double histMin = 2.;

    setTDRStyle();
    writeExtraText = true;
    extraText  = "Work in progress";
    gStyle->SetLineWidth(globalLinedWidth);
    gStyle->SetFrameLineWidth(globalFrameWidth);
    gROOT->ForceStyle();

    TH1::AddDirectory(kFALSE);

    // For nominal histogram
    TH1* hData = NULL;
    TH1* hData_input = NULL;
    TH1* hDY = NULL;
    TH1* hRatio = NULL;

    if(var.Contains("Pt"))
    {
        hData = nomPtUnfold->GetOutput("hUnfoldedPt",0,0,steering,useAxis);
        if(!isType3Closure)
            hDY = nomPtUnfold->GetBias("hDYMCPt",0,0,steering,useAxis);
        else
            hDY = sysPtUnfold["Closure"]["Nominal"]->GetBias("hDYMCPt",0,0,steering,useAxis);
    }
    else
    {
        hData = nomMassUnfold->GetOutput("hUnfoldedMass",0,0,steering,useAxis);
        if(!isType3Closure)
            hDY = nomMassUnfold->GetBias("hDYMCMass",0,0,steering,useAxis);
        else
            hDY = sysMassUnfold["Closure"]["Nominal"]->GetBias("hDYMCMass",0,0,steering,useAxis);
    }
    if(divBinWidth)
    {
        divideByBinWidth(hData, false);
        divideByBinWidth(hDY, false);
    }
    TGraphErrors* gData = histToTGraphError(hData);
    hRatio = (TH1*) hData->Clone("hRatio");

    // Create canvas
    TCanvas* c_out ;
    c_out = new TCanvas("unfolded_level_"+var, "unfolded_level_"+var, 3200, 2800);

    c_out->Draw();
    c_out->cd();

    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetBottomMargin(topPadBottomMargin);
    pad1->SetTopMargin(topPadTopMargin);
    pad1->SetTicks(1);
    pad1->SetLogy();
    if(var.Contains("Mass"))
        pad1->SetLogx();
    pad1->Draw();
    pad1->cd();

    setHistCosmetics(hData, true);
    hData->GetXaxis()->SetLabelSize(0);
    hData->Draw("hist p");
    hData->SetMinimum(histMin);
    hData->SetMarkerSize(4);
    hData->SetLineWidth(4);

    hDY->SetFillColor(kOrange);
    hDY->SetLineColor(kWhite);
    hDY->Draw("hist same");

    TLine massEdgeLine;
    for(int ibin = 2; ibin < hData->GetNbinsX()+1; ibin++)
    {
        massEdgeLine.SetLineStyle(1);
        massEdgeLine.SetLineWidth(2);
        massEdgeLine.SetLineColor(kWhite);
        massEdgeLine.DrawLine(hData->GetXaxis()->GetBinLowEdge(ibin), histMin, hData->GetXaxis()->GetBinLowEdge(ibin), hDY->GetBinContent(ibin));
    }

    TLegend* leg = createLegend(0.6, 0.7);

    if(outName.Contains("Closure"))
    {
        leg->AddEntry(hData, "Unfolded MC", "pl");
    }
    else
    {
        leg->AddEntry(hData, "Unfolded data", "pl");
    }
    leg->AddEntry(hDY, "Drell-Yan (MG5_aMC@NLO)", "F");

    if(hData_input != NULL)
    {
        divideByBinWidth(hData_input, false);
        hData_input->Draw("p9e same");
        hData_input->SetMarkerStyle(20);
        hData_input->SetMarkerSize(3.2);
        hData_input->SetLineColor(kRed);
        hData_input->SetMarkerColor(kRed);

    }
    gData->Draw("pe same");
    gData->SetMarkerSize(4);
    gData->SetLineWidth(5);
    pad1->RedrawAxis();

    leg->Draw();

    int iPeriod_ = 4;
    if(year == 2017)
        iPeriod_ = 5;
    if(year == 2018)
        iPeriod_ = 6;
    CMS_lumi(pad1, iPeriod_, 11);
    // writeCutInfo(pad, var, nthMassBin);
    if(!outName.Contains("Closure"))
        writeCutInfo(pad1, var, nthMassBin);

    // getUnfoldedChi2(TString var, TString steering, bool useAxis, bool divBinWidth)
    TLatex unfoldedChi2;
    unfoldedChi2.SetTextFont(43);
    unfoldedChi2.SetTextSize(50);
    TString chi2_;
    chi2_.Form("%.2f", getUnfoldedChi2(var, steering, useAxis, divBinWidth));
    unfoldedChi2.DrawLatexNDC(0.2, 0.6, "#chi^{2}: " + chi2_);

    //TLine massEdgeLine;
    //if(var=="Mass")
    //{
    //    for(unsigned int i = 0; i < massBinEdges.size(); i++)
    //    {
    //        if(i==0) continue;
    //        if(i==massBinEdges.size()-1) continue;
    //        massEdgeLine.SetLineStyle(2);
    //        massEdgeLine.DrawLine(massBinEdges[i], hData->GetMinimum(), massBinEdges[i], hData->GetMaximum());
    //    }
    //}
    c_out->cd();

    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    pad2->SetTopMargin(bottomPadTopMargin);
    pad2->SetBottomMargin(bottomPadBottomMargin);
    if(var=="Mass")
    {
        pad2->SetLogx();
        hRatio->GetXaxis()->SetMoreLogLabels(true);
    }
    pad2->SetTicks(1);
    //pad2->SetGridy(1);
    pad2->Draw();
    pad2->cd();

    hRatio->SetStats(false);
    hRatio->Divide(hDY);
    hRatio->Draw("hist p");
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(4);
    hRatio->SetLineColor(kBlack);
    hRatio->GetYaxis()->SetTitle("Unfolded/MC");
    hRatio->SetMinimum(0.65);
    hRatio->SetMaximum(1.35);
    hRatio->GetYaxis()->SetDecimals();

    setXaxisTitle(hRatio, var, useAxis);

    TH1* sysBand_ratio_forData = NULL;
    TH1* sysBand_ratio_forMC = NULL;
    TLegend* leg_sys = NULL;
    if(sysName != "")
    {
        leg_sys = new TLegend(0.7, 0.8, 0.95, 0.9,"","brNDC");

        sysBand_ratio_forData = getUnfoldedSystematicBand(var, steering, true, sysName, hData, hDY, hRatio, divBinWidth, true, false, nthMassBin);
        sysBand_ratio_forMC = getUnfoldedSystematicBand(var, steering, true, sysName, hData, hDY, hRatio, divBinWidth, true, true, nthMassBin);

        if(var.Contains("Pt"))
        {
            sysRelPtHist_unfoldedMC[nthMassBin][sysName] = sysBand_ratio_forMC;
            sysRelPtHist_unfoldedData[nthMassBin][sysName] = sysBand_ratio_forData;
        }
        else
        {
            sysRelMassHist_unfoldedMC[sysName] = sysBand_ratio_forMC;
            sysRelMassHist_unfoldedData[sysName] = sysBand_ratio_forData;
        }

        gStyle->SetHatchesSpacing(1);
        //sysBand_ratio_forData->SetFillStyle(3345);
        sysBand_ratio_forData->SetFillColorAlpha(kGray, 0.5);
        sysBand_ratio_forData->SetLineColor(kGray);
        sysBand_ratio_forData->SetLineWidth(1);
        sysBand_ratio_forData->SetMarkerSize(0.);
        sysBand_ratio_forData->Draw("E2 same");

        sysBand_ratio_forMC->SetFillStyle(3354);
        sysBand_ratio_forMC->SetFillColor(kOrange);
        sysBand_ratio_forMC->SetMarkerSize(0.);
        sysBand_ratio_forMC->SetLineColor(kOrange);
        sysBand_ratio_forMC->SetLineWidth(1);
        sysBand_ratio_forMC->Draw("E2 same");

        leg_sys->SetTextFont(43);
        leg_sys->SetTextSize(100);
        leg_sys->SetFillStyle(0); // transparent
        leg_sys->SetBorderSize(0);
        TString sysName_ = getSysNameToShow(sysName);
        leg_sys->AddEntry(sysBand_ratio_forMC, sysName_, "F");
        leg_sys->Draw();
    }

    TLine* l_ = new TLine(hRatio->GetXaxis()->GetXmin(),1,hRatio->GetXaxis()->GetXmax(),1);
    l_->SetLineColor(kOrange);
    l_->Draw("same");
    l_->SetLineStyle(1);

    TGraphErrors* gRatio = histToTGraphError(hRatio); 
    gRatio->Draw("p same");
    gRatio->SetLineWidth(5);
    gRatio->SetMarkerSize(4);

    // Save canvas
    c_out->cd();
    c_out->SaveAs(outName!=""?output_baseDir+outName+var+sysName+".pdf":output_baseDir+"unfolded_"+var+sysName+".pdf");

    return c_out;
}

TCanvas* ISRUnfold::drawAcceptCorrHists(TString var, TString filePath, TString binDef, TString steering, bool useAxis, TString sysName, TString outName, int nthMassBin, bool divBinWidth)
{
    double histMin = 5e-1;

    setTDRStyle();
    writeExtraText = true;
    extraText  = "Work in progress";
    gStyle->SetLineWidth(globalLinedWidth);
    gStyle->SetFrameLineWidth(globalFrameWidth);
    gROOT->ForceStyle();

    TH1::AddDirectory(kFALSE);

    TFile* filein = new TFile(filePath);

    TH1* hData = NULL;
    TH1* hDY_raw = NULL;
    TH1* hDY = NULL;
    TH1* hRatio = NULL;

    if(var.Contains("Pt"))
    {
        hData = pt_binning_Gen->ExtractHistogram("hData", hFullPhasePtData, 0, useAxis, steering);

        hDY_raw = (TH1*) filein->Get("Acceptance/PtGen" + binDef + "/histo_DYJets");
        if(year==2016)
            hDY_raw->Add((TH1*) filein->Get("Acceptance/PtGen" + binDef + "/histo_DYJets10to50"));
        else
            hDY_raw->Add((TH1*) filein->Get("Acceptance/PtGen" + binDef + "/histo_DYJets10to50_MG"));
        hDY = pt_binning_Gen->ExtractHistogram("hDY", hDY_raw, 0, useAxis, steering);
    }
    else
    {
        hData = mass_binning_Gen->ExtractHistogram("hData", hFullPhaseMassData, 0, useAxis, steering);
        hDY_raw = (TH1*) filein->Get("Acceptance/MassGen" + binDef + "/histo_DYJets");
        if(year==2016)
            hDY_raw->Add((TH1*) filein->Get("Acceptance/MassGen" + binDef + "/histo_DYJets10to50"));
        else
            hDY_raw->Add((TH1*) filein->Get("Acceptance/MassGen" + binDef + "/histo_DYJets10to50_MG"));
        hDY = mass_binning_Gen->ExtractHistogram("hDY", hDY_raw, 0, useAxis, steering);
    }
    if(divBinWidth)
    {
        divideByBinWidth(hData, false);
        divideByBinWidth(hDY, false);
    }
    TGraphErrors* gData = histToTGraphError(hData);
    hRatio = (TH1*) hData->Clone("hRatio");

    TCanvas* c_out = new TCanvas("acceptance_corrected_"+var, "acceptance_corrected_"+var, 3200, 2800);
    c_out->Draw();
    c_out->cd();

    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetBottomMargin(topPadBottomMargin);
    pad1->SetTopMargin(topPadTopMargin);
    pad1->SetTicks(1);
    pad1->SetLogy();
    if(var.Contains("Mass"))
        pad1->SetLogx();
    pad1->Draw();
    pad1->cd();

    hData->SetTitle("");
    hData->SetStats(false);
    hData->GetXaxis()->SetMoreLogLabels(true);
    hData->Draw("hist p");
    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(4);
    hData->SetLineColor(kBlack);
    hData->GetYaxis()->SetTitle("Events/Bin");
    hData->SetMaximum(hDY->GetMaximum() * 1e5);
    hData->SetMinimum(histMin);

    hDY->SetFillColor(kOrange);
    hDY->Draw("hist same");

    TH1* sysBand_forMC = NULL;
    if(sysName != "")
    {
        sysBand_forMC = getUnfAcceptSystematicBand(var, steering, true, sysName, hData, hDY, hDY, divBinWidth, false, true, nthMassBin);
        sysBand_forMC->SetFillColorAlpha(kRed,0.6);

        if(var.Contains("Pt"))
        {
            sysAbsPtHist_unfoldedAcceptMC[nthMassBin][sysName] = sysBand_forMC;
        }
        else
        {
            sysAbsMassHist_unfoldedAcceptMC[sysName] = sysBand_forMC;
        }

        sysBand_forMC->SetFillStyle(3001);
        sysBand_forMC->SetMarkerSize(0.);
        //sysBand_forMC->Draw("E2 same");
    }

    TLine massEdgeLine;
    for(int ibin = 2; ibin < hData->GetNbinsX()+1; ibin++)
    {
        massEdgeLine.SetLineStyle(1);
        massEdgeLine.SetLineWidth(1);
        massEdgeLine.SetLineColor(kWhite);
        massEdgeLine.DrawLine(hData->GetXaxis()->GetBinLowEdge(ibin), histMin, hData->GetXaxis()->GetBinLowEdge(ibin), hDY->GetBinContent(ibin));
    }

    TLegend* leg = createLegend(0.6, 0.7);
    leg->AddEntry(hData, "Unfolded data", "pl");
    leg->AddEntry(hDY, "Drell-Yan (MG5_aMC@NLO)", "F");

    gData->Draw("pe same");
    gData->SetMarkerSize(4);
    gData->SetLineWidth(5);
    pad1->RedrawAxis();

    leg->Draw();

    int iPeriod_ = 4;
    if(year == 2017)
        iPeriod_ = 5;
    if(year == 2018)
        iPeriod_ = 6;
    CMS_lumi(pad1, iPeriod_, 11);

    writeCutInfo(pad1, var, nthMassBin);

    c_out->cd();

    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    pad2->SetTopMargin(bottomPadTopMargin);
    pad2->SetBottomMargin(bottomPadBottomMargin);
    if(var=="Mass")
    {
        pad2->SetLogx();
        hRatio->GetXaxis()->SetMoreLogLabels(true);
    }
    pad2->SetTicks(1);
    //pad2->SetGridy(1);
    pad2->Draw();
    pad2->cd();

    hRatio->SetStats(false);
    hRatio->Divide(hDY);
    hRatio->Draw("hist p");
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(4);
    hRatio->SetLineColor(kBlack);
    hRatio->GetYaxis()->SetTitle("Unfolded/MC");

    hRatio->SetMinimum(0.65);
    hRatio->SetMaximum(1.35);
    hRatio->GetYaxis()->SetDecimals();

    setXaxisTitle(hRatio, var, useAxis);

    TH1* sysBand_ratio_forData = NULL;
    TH1* sysBand_ratio_forMC = NULL;
    TLegend* leg_sys = NULL;

    if(sysName != "")
    {
        leg_sys = new TLegend(0.7, 0.8, 0.95, 0.9,"","brNDC");

        // Caution! Systematic band keep increasing if getUnfAcceptSystematicBand is called repeatedly 
        sysBand_ratio_forData = getUnfAcceptSystematicBand(var, steering, true, sysName, hData, hDY, hRatio, divBinWidth, true, false, nthMassBin);
        sysBand_ratio_forMC = getUnfAcceptSystematicBand(var, steering, true, sysName, hData, hDY, hRatio, divBinWidth, true, true, nthMassBin);

        if(var.Contains("Pt"))
        {
            sysRelPtHist_unfoldedAcceptMC[nthMassBin][sysName] = sysBand_ratio_forMC;
            sysRelPtHist_unfoldedAcceptData[nthMassBin][sysName] = sysBand_ratio_forData;
        }
        else
        {
            sysRelMassHist_unfoldedAcceptMC[sysName] = sysBand_ratio_forMC;
            sysRelMassHist_unfoldedAcceptData[sysName] = sysBand_ratio_forData;
        }
        
        gStyle->SetHatchesSpacing(1);
        sysBand_ratio_forData->SetFillColorAlpha(kGray, 0.5);
        sysBand_ratio_forData->SetLineColor(kGray);
        sysBand_ratio_forData->SetLineWidth(1);
        sysBand_ratio_forData->SetMarkerSize(0.);
        sysBand_ratio_forData->Draw("E2 same");

        sysBand_ratio_forMC->SetFillStyle(3354);
        sysBand_ratio_forMC->SetFillColor(kOrange);
        sysBand_ratio_forMC->SetMarkerSize(0.);
        sysBand_ratio_forMC->SetLineColor(kOrange);
        sysBand_ratio_forMC->SetLineWidth(1);
        sysBand_ratio_forMC->Draw("E2 same");

        leg_sys->SetTextFont(43);
        leg_sys->SetTextSize(100);
        leg_sys->SetFillStyle(0); // transparent
        leg_sys->SetBorderSize(0);
        leg_sys->AddEntry(sysBand_ratio_forMC, sysName, "F");
        leg_sys->Draw();

    }

    TLine* l_ = new TLine(hRatio->GetXaxis()->GetXmin(),1,hRatio->GetXaxis()->GetXmax(),1);
    l_->SetLineColor(kOrange);
    l_->Draw("same");
    l_->SetLineStyle(1);

    TGraphErrors* gRatio = histToTGraphError(hRatio); 
    gRatio->Draw("p same");
    gRatio->SetLineWidth(5);
    gRatio->SetMarkerSize(4);

    c_out->cd();
    if( ((sysName).Contains("Scale") && !(sysName).Contains("Lep"))) sysName = "Scale";
    //c_out->SaveAs(outName!=""?output_baseDir+outName+var+sysName+".pdf":output_baseDir+"unfoldedAccept_"+var+sysName+".pdf");
    return 0;
}

TH1* ISRUnfold::getUnfoldedSystematicBand(TString var, TString steering, bool useAxis, TString sysName, TH1* hData, TH1* hDY, TH1* hRatio, bool divBinWidth, bool isRatio, bool forMC, int nthMassBin)
{
    TH1* sysBand_ratio = NULL;
    TProfile* hprofForPDF = NULL;

    if(!sysName.Contains("totalSys"))
    {
        int variationSize = sysMap[sysName].size();
        vector<TH1*> hRatio_temp;
        for(int ith = 0; ith < variationSize; ith++)
        {
            TH1* hSYS_temp = NULL;
            //TH1* hRatio_temp = NULL;

            if(var.Contains("Pt"))
            {
                if(sysName.Contains("iterEM") && !sysMap[sysName][ith].Contains("Nominal"))
                {
                    if(forMC)
                        hSYS_temp = nomPtUnfold->GetBias("hDYMCPt_temp",0,0,steering,useAxis);
                    else
                    {
                        hSYS_temp = iterEMPtUnfold->GetOutput("hUnfoldedPt_temp",0,0,steering,useAxis);
                    }
                }
                else
                {
                    if(forMC)
                        hSYS_temp = sysPtUnfold[sysName][sysMap[sysName][ith]]->GetBias("hDYMCPt_temp",0,0,steering,useAxis);
                    else
                        hSYS_temp = sysPtUnfold[sysName][sysMap[sysName][ith]]->GetOutput("hUnfoldedPt_temp",0,0,steering,useAxis);
                }
            }
            else
            {
                if(sysName.Contains("iterEM") && !sysMap[sysName][ith].Contains("Nominal"))
                {
                    if(forMC)
                        hSYS_temp = nomMassUnfold->GetBias("hDYMCMass_temp",0,0,steering,useAxis);
                    else
                    {
                        hSYS_temp = iterEMMassUnfold->GetOutput("hUnfoldedMass_temp",0,0,steering,useAxis);
                    }
                }
                else
                {
                    if(forMC)
                        hSYS_temp = sysMassUnfold[sysName][sysMap[sysName][ith]]->GetBias("hDYMCMass_temp",0,0,steering,useAxis);
                    else
                        hSYS_temp = sysMassUnfold[sysName][sysMap[sysName][ith]]->GetOutput("hUnfoldedMass_temp",0,0,steering,useAxis);
                }
            }
            if(divBinWidth)
            {
                divideByBinWidth(hSYS_temp, false);
            }

            if(forMC)
            {
                hRatio_temp.push_back((TH1*) hData->Clone("hRatio_temp"));
                hRatio_temp.at(ith)->Divide(hSYS_temp);
                if(ith==0)
                {
                    //sysBand_ratio = (TH1*)hRatio_temp.at(ith)->Clone("sysUnfBand_ratio"+sysName);
                    sysBand_ratio = cloneEmptyHist(hRatio_temp.at(ith), "sysUnfBand_ratio"+sysName);
                    if(sysName.Contains("PDF")) hprofForPDF = (TProfile*) cloneHistToTProf(sysBand_ratio, "PDFUncertainty");
                }
            }
            else
            {
                hRatio_temp.push_back((TH1*) hSYS_temp->Clone("hRatio_temp"));
                hRatio_temp.at(ith)->Divide(hDY);
                if(ith==0)
                {
                    //sysBand_ratio = (TH1*)hRatio_temp.at(ith)->Clone("sysUnfBand_ratio"+sysName);
                    sysBand_ratio = cloneEmptyHist(hRatio_temp.at(ith), "sysUnfBand_ratio"+sysName);
                    if(sysName.Contains("PDF")) hprofForPDF = (TProfile*) cloneHistToTProf(sysBand_ratio, "PDFUncertainty");
                }
            }
            for(int ibin = 1; ibin < sysBand_ratio->GetNbinsX()+1; ibin++)
            {
                double delta = fabs(hRatio_temp.at(ith)->GetBinContent(ibin) - hRatio->GetBinContent(ibin));
                if(sysName.Contains("FSR"))
                {
                    if(ith == 0) break;
                    delta = fabs(hRatio_temp.at(0)->GetBinContent(ibin) - hRatio_temp.at(1)->GetBinContent(ibin));
                }
                else if(sysName.Contains("PDF"))
                {
                    // PDF
                    hprofForPDF->Fill(sysBand_ratio->GetXaxis()->GetBinCenter(ibin), delta);
                }
                else
                {
                    if(ith != 0)
                    {
                        delta = delta > sysBand_ratio->GetBinError(ibin) ? delta : sysBand_ratio->GetBinError(ibin);
                    }
                }

                if(sysName.Contains("PDF"))
                {
                    sysBand_ratio->SetBinError(ibin, hprofForPDF->GetBinError(ibin));
                    if(forMC)
                        sysBand_ratio->SetBinContent(ibin, 1.);
                    else
                        sysBand_ratio->SetBinContent(ibin, hRatio->GetBinContent(ibin)); 
                }
                else
                {
                    if(delta < 1e-5)
                        sysBand_ratio->SetBinError(ibin, 1e-5);
                    else
                        sysBand_ratio->SetBinError(ibin, delta);

                    if(forMC)
                    {
                        sysBand_ratio->SetBinContent(ibin, 1.);
                    }
                    else
                    {
                        sysBand_ratio->SetBinContent(ibin, hRatio->GetBinContent(ibin));
                    }
                }
            }

            delete hSYS_temp; // This could be data or MC systematic variation
        }// end of variation loop
        delete hprofForPDF;
        hRatio_temp.clear();
    }
    else
    {
            map<TString, TH1*>::iterator it;
            map<TString, TH1*>::iterator end;

            if(var.Contains("Pt"))
            {
                if(forMC)
                {
                    it = sysRelPtHist_unfoldedMC[nthMassBin].begin();
                    end = sysRelPtHist_unfoldedMC[nthMassBin].end();
                }
                else
                {
                    it = sysRelPtHist_unfoldedData[nthMassBin].begin();
                    end = sysRelPtHist_unfoldedData[nthMassBin].end();
                }
            }
            else
            {
                if(forMC)
                {
                    it = sysRelMassHist_unfoldedMC.begin();
                    end = sysRelMassHist_unfoldedMC.end();
                }
                else
                {
                    it = sysRelMassHist_unfoldedData.begin();
                    end = sysRelMassHist_unfoldedData.end();
                }
            }
            bool firstSys = true;
            while(it != end)
            {
                if(firstSys)
                {
                    sysBand_ratio = (TH1*)(it->second)->Clone("sysBand_ratioTotal");
                    for(int ibin = 1; ibin < sysBand_ratio->GetNbinsX()+1; ibin++)
                    {
                        sysBand_ratio->SetBinError(ibin, pow(sysBand_ratio->GetBinError(ibin), 2));
                    }
                    it++;
                    firstSys = false;
                    continue;
                }
                else
                {
                    for(int ibin = 1; ibin < sysBand_ratio->GetNbinsX()+1; ibin++)
                    {
                        sysBand_ratio->SetBinError(ibin, sysBand_ratio->GetBinError(ibin) + pow((it->second)->GetBinError(ibin), 2));
                    }
                    it++;
                }
            }
            for(int ibin = 1; ibin < sysBand_ratio->GetNbinsX()+1; ibin++)
            {
                sysBand_ratio->SetBinError(ibin, sqrt(sysBand_ratio->GetBinError(ibin)));
            }
    }

    return sysBand_ratio;
}

TH1* ISRUnfold::getUnfAcceptSystematicBand(TString var, TString steering, bool useAxis, TString sysName, TH1* hData, TH1* hDY, TH1* hRatio, bool divBinWidth, bool isRatio, bool forMC, int nthMassBin)
{
    TH1* sysBand_ratio = NULL;
    TProfile* hprofForPDF = NULL;

    if(!sysName.Contains("totalSys"))
    {
        int variationSize = sysMap[sysName].size();
        //cout << "sys name: " << sysName << " size:" << variationSize << endl;
        vector<TH1*> hRatio_temp;
        for(int ith = 0; ith < variationSize; ith++)
        {
            TH1* hSYS_raw_temp = NULL;
            TH1* hSYS_temp = NULL;
            //TH1* hRatio_temp = NULL; // TODO use vector to save each ratio histogram

            if(var.Contains("Pt"))
            {
                if(forMC)
                    hSYS_raw_temp = hSysFullPhasePtMC[sysName][sysMap[sysName][ith]];
                else
                    hSYS_raw_temp = hSysFullPhasePtData[sysName][sysMap[sysName][ith]];

                hSYS_temp = pt_binning_Gen->ExtractHistogram("UnfAcceptPt_"+sysMap[sysName][ith], hSYS_raw_temp, 0, useAxis, steering);
            }
            else
            {
                if(forMC)
                    hSYS_raw_temp = hSysFullPhaseMassMC[sysName][sysMap[sysName][ith]];
                else
                    hSYS_raw_temp = hSysFullPhaseMassData[sysName][sysMap[sysName][ith]];

                hSYS_temp = mass_binning_Gen->ExtractHistogram("UnfAcceptMass_"+sysMap[sysName][ith], hSYS_raw_temp, 0, useAxis, steering);
            }

            if(divBinWidth)
            {
                divideByBinWidth(hSYS_temp, false);
            }

            if(forMC)
            {

                if(isRatio)
                {
                    hRatio_temp.push_back((TH1*) hData->Clone("hRatio_temp"));
                    hRatio_temp.at(ith)->Divide(hSYS_temp);
                }
                else
                {
                    hRatio_temp.push_back((TH1*) hSYS_temp->Clone("hRatio_temp"));
                }
                if(ith==0)
                {
                    sysBand_ratio = (TH1*)hRatio_temp.at(ith)->Clone("sysUnfBand_ratio"+sysName);
                    if(sysName.Contains("PDF")) hprofForPDF = (TProfile*) cloneHistToTProf(sysBand_ratio, "PDFUncertainty");  
                }
            }
            else
            {
                hRatio_temp.push_back((TH1*) hSYS_temp->Clone("hRatio_temp"));
                hRatio_temp.at(ith)->Divide(hDY);
                if(ith==0)
                {
                    sysBand_ratio = (TH1*)hRatio_temp.at(ith)->Clone("sysUnfBand_ratio"+sysName);
                    if(sysName.Contains("PDF")) hprofForPDF = (TProfile*) cloneHistToTProf(sysBand_ratio, "PDFUncertainty");
                }
            }

            for(int ibin = 1; ibin < sysBand_ratio->GetNbinsX()+1; ibin++)
            {
                double delta = fabs(hRatio_temp.at(ith)->GetBinContent(ibin) - hRatio->GetBinContent(ibin));
                if(sysName.Contains("FSR"))
                {
                    if(ith == 0) break;
                    delta = fabs(hRatio_temp.at(0)->GetBinContent(ibin) - hRatio_temp.at(1)->GetBinContent(ibin));
                }
                else if(sysName.Contains("PDF"))
                {
                    // PDF
                    hprofForPDF->Fill(sysBand_ratio->GetXaxis()->GetBinCenter(ibin), delta);
                }
                else
                {
                    if(ith != 0)
                    {
                        delta = delta > sysBand_ratio->GetBinError(ibin) ? delta : sysBand_ratio->GetBinError(ibin);
                    }
                }

                if(sysName.Contains("PDF"))
                {
                    sysBand_ratio->SetBinError(ibin, hprofForPDF->GetBinError(ibin));
                    if(forMC)
                    {
                        if(isRatio) sysBand_ratio->SetBinContent(ibin, 1.);
                        else sysBand_ratio->SetBinContent(ibin, hDY->GetBinContent(ibin));
                    }
                    else
                    {
                        sysBand_ratio->SetBinContent(ibin, hRatio->GetBinContent(ibin)); 
                    }
                }
                else
                {
                    if(delta < 1e-5)
                        sysBand_ratio->SetBinError(ibin, 1e-5);
                    else
                        sysBand_ratio->SetBinError(ibin, delta);

                    if(forMC)
                    {
                        if(isRatio) sysBand_ratio->SetBinContent(ibin, 1.);
                        else sysBand_ratio->SetBinContent(ibin, hDY->GetBinContent(ibin));
                    }
                    else
                    {
                        sysBand_ratio->SetBinContent(ibin, hRatio->GetBinContent(ibin));
                    }
                }
            }

            delete hSYS_temp; // This could be data or MC systematic variation
        }
        delete hprofForPDF;
        hRatio_temp.clear();
    }
    else
    {
            map<TString, TH1*>::iterator it;
            map<TString, TH1*>::iterator end;

            if(var.Contains("Pt"))
            {
                if(forMC)
                {
                    it = sysRelPtHist_unfoldedAcceptMC[nthMassBin].begin();
                    end = sysRelPtHist_unfoldedAcceptMC[nthMassBin].end();
                }
                else
                {
                    it = sysRelPtHist_unfoldedAcceptData[nthMassBin].begin();
                    end = sysRelPtHist_unfoldedAcceptData[nthMassBin].end();
                }
            }
            else
            {
                if(forMC)
                {
                    it = sysRelMassHist_unfoldedAcceptMC.begin();
                    end = sysRelMassHist_unfoldedAcceptMC.end();
                }
                else
                {
                    it = sysRelMassHist_unfoldedAcceptData.begin();
                    end = sysRelMassHist_unfoldedAcceptData.end();
                }
            }
            bool firstSys = true;
            while(it != end)
            {
                if(firstSys)
                {
                    sysBand_ratio = (TH1*)(it->second)->Clone("sysBand_ratioTotal");
                    for(int ibin = 1; ibin < sysBand_ratio->GetNbinsX()+1; ibin++)
                    {
                        sysBand_ratio->SetBinError(ibin, pow(sysBand_ratio->GetBinError(ibin), 2));
                    }
                    it++;
                    firstSys = false;
                    continue;
                }
                else
                {
                    for(int ibin = 1; ibin < sysBand_ratio->GetNbinsX()+1; ibin++)
                    {
                        sysBand_ratio->SetBinError(ibin, sysBand_ratio->GetBinError(ibin) + pow((it->second)->GetBinError(ibin), 2));
                    }
                    it++;
                }
            }
            for(int ibin = 1; ibin < sysBand_ratio->GetNbinsX()+1; ibin++)
            {
                sysBand_ratio->SetBinError(ibin, sqrt(sysBand_ratio->GetBinError(ibin)));
            }
    }

    return sysBand_ratio;
}

TH1* ISRUnfold::getDetectorSystematicBand(TString var, TString filePath, TString dirName, TString steering, bool useAxis, TString sysName, TH1* hData, TH1* hDY, TH1* hMCtotal, TH1* hRatio,
                                            bool divBinWidth, bool isRatio, bool forMC, TString sysFilePath, int nthMassBin, bool isBkgSubData)
{

    TH1* sysBand_ratio = NULL;
    TProfile* hprofForPDF = NULL;

    if(sysFilePath != "")
    {
        filePath = sysFilePath;
    }

    TString dataHistName_ = "histo_DoubleMuon";
    TString DYHistName_ = "histo_DYJetsToMuMu";
    if(channel_name == "electron")
    {
        dataHistName_ = "histo_DoubleEG";
        if(year==2018)
            dataHistName_ = "histo_EGamma";
        DYHistName_  = "histo_DYJetsToEE";
    }

    if(!sysName.Contains("totalSys"))
    {
        int variationSize = sysMap[sysName].size();
        for(int ith = 0; ith < variationSize; ith++)
        {
            TH1* hDYSYS_temp = NULL; // This could be data or MC systematic variation
            TH1* hDataSYS_temp = NULL; // This could be data or MC systematic variation
            TH1* hMCtotal_temp = NULL;
            TH1* hRatio_temp = NULL;

            // Dummy THStack, TLegend
            THStack* hsMC_temp;
            TLegend* leg = new TLegend(0.5, 0.7, 0.9, 0.9,"","brNDC");;
            hsMC_temp = new THStack("hsMC_temp", "hsMC_temp");

            TString tempDYHistName = DYHistName_;
            TString tempDYDirName = dirName;
            TString tempDYVar = var;
            TString tempDYNewHistName = "Signal"+sysMap[sysName][ith];

            TString tempDataHistName = dataHistName_;
            TString tempDataDirName = dirName;
            TString tempDataVar = var;
            TString tempDataNewHistName = "Data"+sysMap[sysName][ith];

            if(sysName.Contains("LepScale") || sysName.Contains("LepRes"))
            {
                tempDYHistName = DYHistName_; 
                tempDYDirName = dirName + "_" + sysMap[sysName][ith];
                tempDataDirName = dirName + "_" + sysMap[sysName][ith];
                tempDYVar = var;  
            }
            else
            {
                tempDYHistName = DYHistName_ + "_" + sysMap[sysName][ith];
                if(sysName.Contains("Unfold") || sysName.Contains("ZptCorr"))
                    tempDYHistName = DYHistName_;
            }

            if(sysMap[sysName][ith] != "Nominal")
            {
                hDataSYS_temp = getRawHist(tempDataVar, filePath, tempDataDirName, tempDataHistName, tempDataNewHistName, steering, useAxis, divBinWidth);
                hDYSYS_temp = getRawHist(tempDYVar, filePath, tempDYDirName, tempDYHistName, tempDYNewHistName, steering, useAxis, divBinWidth);

                hMCtotal_temp = (TH1*) hDYSYS_temp->Clone("hMCtotal_temp");

                if(forMC)
                    hRatio_temp = (TH1*) hData->Clone("hRatio_temp");
                else
                    hRatio_temp = (TH1*) hDataSYS_temp->Clone("hRatio_temp");

                setTHStack(tempDYVar, filePath, tempDYDirName, *hsMC_temp, *hMCtotal_temp, *leg, steering, useAxis, sysMap[sysName][ith], divBinWidth);

                delete hsMC_temp;
                delete leg;
            }
            else
            {
                if(forMC)
                {
                    hRatio_temp = (TH1*) hData->Clone("hRatio_down");
                }
                else
                {
                    hRatio_temp = (TH1*) hData->Clone("hRatio_down");
                }
                if(forMC)
                {
                    hMCtotal_temp = (TH1*) hMCtotal->Clone("hMCtotal_down");
                }
            }

            // Set systematic band
            // Create histogram
            if(forMC)
            {
                if(isBkgSubData)
                {
                    TH1* hMCtotalBkgSub = (TH1*) hMCtotal->Clone("hMCtotalBkgSub");
                    hRatio_temp->Add(hMCtotalBkgSub, -1);
                    hRatio_temp->Divide(hDYSYS_temp);
                    if(ith==0){
                        sysBand_ratio = (TH1*)hRatio_temp->Clone("sysBand_ratio"+sysName);
                        if(sysName.Contains("PDF")) hprofForPDF = (TProfile*) cloneHistToTProf(sysBand_ratio, "PDFUncertainty");
                    }

                    delete hMCtotalBkgSub;
                }
                else
                {
                    hRatio_temp->Divide(hMCtotal_temp);
                    if(ith==0){
                        sysBand_ratio = (TH1*)hRatio_temp->Clone("sysBand_ratio"+sysName);
                        if(sysName.Contains("PDF")) hprofForPDF = (TProfile*) cloneHistToTProf(sysBand_ratio, "PDFUncertainty");
                    }
                }
            }
            else
            {
                if(isBkgSubData)
                {
                    TH1* hMCtotalBkgSub = (TH1*) hMCtotal_temp->Clone("hMCtotalBkgSub");
                    hMCtotalBkgSub->Add(hDYSYS_temp, -1);
                    hRatio_temp->Add(hMCtotalBkgSub, -1);
                    hRatio_temp->Divide(hDY);
                    if(ith==0){
                        sysBand_ratio = (TH1*)hRatio_temp->Clone("sysBand_ratio"+sysName);
                        if(sysName.Contains("PDF")) hprofForPDF = (TProfile*) cloneHistToTProf(sysBand_ratio, "PDFUncertainty");
                    }

                    delete hMCtotalBkgSub;
                }
                else
                {
                    hRatio_temp->Divide(hMCtotal);
                    if(ith==0){
                        sysBand_ratio = (TH1*)hRatio_temp->Clone("sysBand_ratio"+sysName);
                        if(sysName.Contains("PDF")) hprofForPDF = (TProfile*) cloneHistToTProf(sysBand_ratio, "PDFUncertainty");
                    }
                }
            }

            // Update error
            for(int ibin = 1; ibin < sysBand_ratio->GetNbinsX()+1; ibin++)
            {
                double delta = fabs(hRatio_temp->GetBinContent(ibin) - hRatio->GetBinContent(ibin));
                if(sysName.Contains("PDF"))
                {
                    hprofForPDF->Fill(sysBand_ratio->GetXaxis()->GetBinCenter(ibin), delta);

                    if(hprofForPDF->GetBinError(ibin) < 1e-5)
                        sysBand_ratio->SetBinError(ibin, 1e-5);
                    else
                        sysBand_ratio->SetBinError(ibin, hprofForPDF->GetBinError(ibin));
        
                    if(forMC)
                    {
                        if(isRatio) sysBand_ratio->SetBinContent(ibin, 1.);
                        else sysBand_ratio->SetBinContent(ibin, hDY->GetBinContent(ibin));
                    }
                    else
                    {
                        sysBand_ratio->SetBinContent(ibin, hRatio->GetBinContent(ibin)); 
                    }
                }
                else
                {
                    if(ith != 0)
                    {
                        delta = delta > sysBand_ratio->GetBinError(ibin) ? delta : sysBand_ratio->GetBinError(ibin);
                    }

                    if(delta < 1e-5)
                        sysBand_ratio->SetBinError(ibin, 1e-5);
                    else
                        sysBand_ratio->SetBinError(ibin, delta);

                    if(forMC)
                    {
                        sysBand_ratio->SetBinContent(ibin, 1.);
                    }
                    else
                    {
                        sysBand_ratio->SetBinContent(ibin, hRatio->GetBinContent(ibin));
                    }
                }

            }

            delete hDataSYS_temp; // This could be data or MC systematic variation
            delete hDYSYS_temp; // This could be data or MC systematic variation
            delete hMCtotal_temp;
            delete hRatio_temp;
        }
    }
    else
    {
            map<TString, TH1*>::iterator it;
            map<TString, TH1*>::iterator end;

            if(var.Contains("Pt"))
            {
                if(forMC)
                {
                    it = sysRelPtHist_detectorMC[nthMassBin].begin();
                    end = sysRelPtHist_detectorMC[nthMassBin].end();
                }
                else
                {
                    it = sysRelPtHist_detectorData[nthMassBin].begin();
                    end = sysRelPtHist_detectorData[nthMassBin].end();
                }
            }
            else
            {
                if(forMC)
                {
                    it = sysRelMassHist_detectorMC.begin();
                    end = sysRelMassHist_detectorMC.end();
                }
                else
                {
                    it = sysRelMassHist_detectorData.begin();
                    end = sysRelMassHist_detectorData.end();
                }
            }
            bool firstSys = true;
            while(it != end)
            {
                if(firstSys)
                {
                    sysBand_ratio = cloneEmptyHist(it->second, "sysBand_ratioTotal");
                    for(int ibin = 1; ibin < sysBand_ratio->GetNbinsX()+1; ibin++)
                    {
                        sysBand_ratio->SetBinError(ibin, pow((it->second)->GetBinError(ibin), 2));
                        if(forMC)
                            sysBand_ratio->SetBinContent(ibin, 1.);
                        else
                            sysBand_ratio->SetBinContent(ibin, hRatio->GetBinContent(ibin));
                    }
                    it++;
                    firstSys = false;
                    continue;
                }
                else
                {
                    for(int ibin = 1; ibin < sysBand_ratio->GetNbinsX()+1; ibin++)
                    {
                        sysBand_ratio->SetBinError(ibin, sysBand_ratio->GetBinError(ibin) + pow((it->second)->GetBinError(ibin), 2));
                    }
                    it++;
                }
            }
            for(int ibin = 1; ibin < sysBand_ratio->GetNbinsX()+1; ibin++)
            {
                sysBand_ratio->SetBinError(ibin, sqrt(sysBand_ratio->GetBinError(ibin)));
            }
    }

    return sysBand_ratio;
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

void ISRUnfold::setXaxisTitle(TH1* hist, TString var, bool useAxis, TString title)
{
    TString channel_name_;
    if(channel_name=="electron") channel_name_ = "ee";
    else channel_name_ = "#mu#mu";

    if(useAxis)
    {
        TString title_;
        if(var.Contains("Mass"))
        {
            title_ = "Mass^{" + channel_name_  + "} [GeV]";
            if(title != "") title_ = title;
            hist->GetXaxis()->SetTitle(title_);
        }
        if(var.Contains("Pt"))
        {
            title_ = "p_{T}^{" + channel_name_  + "} [GeV]";
            if(title != "") title_ = title;
            hist->GetXaxis()->SetTitle(title_);
        }
    }
    else
    {
        if(var.Contains("Mass"))
        {
            hist->GetXaxis()->SetTitle("Mass bin index");
        }
        if(var.Contains("Pt"))
        {
            hist->GetXaxis()->SetTitle("p_{T} bin index");
        }
    }
}

void ISRUnfold::writeCutInfo(TPad* pad, TString var, int nthMassBin, double x_, double y_, bool showLepCuts)
{
    pad->cd();

    TString mass_cut_info;
    TString lepton_cut_info;
    TString lepton_type;
    if(channel_name=="electron"){
        lepton_type = "ee";
        lepton_cut_info = "p_{T}>25(20) GeV, |#eta|<2.5";
    }
    else
    {
        lepton_type = "#mu#mu";
        lepton_cut_info = "p_{T}>20(10) GeV, |#eta|<2.4";
    }

    TString low_bound_, upper_bound_;
    low_bound_.Form("%d", (int)massBinEdges[nthMassBin]);
    upper_bound_.Form("%d", (int)massBinEdges[nthMassBin+1]);
    mass_cut_info = low_bound_ + " < M(" + lepton_type + ") < " + upper_bound_ + " (GeV)";

    TLatex dimass_cut;
    TLatex lepton_cut;
    dimass_cut.SetTextFont(43);
    dimass_cut.SetTextSize(40 * 2.5);
    lepton_cut.SetTextFont(43);
    lepton_cut.SetTextSize(40 * 2.5);

    if(var.Contains("Pt"))
    {
        dimass_cut.DrawLatexNDC(x_, y_, mass_cut_info);
        if(showLepCuts) lepton_cut.DrawLatexNDC(x_, y_-0.05, lepton_cut_info);
    }
    if(var.Contains("Mass"))
    {
        //dimass_cut.DrawLatexNDC(x_, y_, mass_cut_info);
        lepton_cut.DrawLatexNDC(x_, y_, lepton_cut_info);
    }
}

void ISRUnfold::setTHStack(TString var, TString filePath, TString dirName, THStack& hs, TH1& hMCtotal, TLegend& leg, TString steering, bool useAxis, TString sysName, bool divBinWidth, bool updateLegend, bool doPtCombine)
{
    TH1::AddDirectory(kFALSE);
    int bkgSize = bkgNames.size();
    map<TString, TH1*> tempMap_legend;
    vector<TString> temp_bkgName;

    // Count total number of each background type N
    map<TString, int> bkgTypeN;
    //cout << "N bkg: " << bkgSize << endl;
    for(int i = 0; i < bkgSize; i++)
    {
        if(bkgTypes[i] == "DY") continue;
        map<TString, int>::iterator it = bkgTypeN.find(bkgTypes[i]);
        if(it != bkgTypeN.end())
        {
            bkgTypeN[bkgTypes[i]]++;
        }
        else
        {
            //cout << bkgTypes[i] << " first found" << endl;
            bkgTypeN[bkgTypes[i]] = 1;
        }
    }

    TH1* htemp = NULL;
    bool isFirstBkg = true;
    int nthBkg = 0;

    for(int i = 0; i < bkgSize; i++)
    {
        if(bkgTypes[i] == "DY") continue;
        TString histName_ = "histo_" + bkgNames[i];;
        TString dirName_ = dirName;

        if(isFirstBkg)
        {
            if(sysName == "")
            {
                htemp = getRawHist(var, filePath, dirName_, histName_, "h"+bkgNames[i], steering, useAxis, divBinWidth);
            }
            else
            {
                TString histPostfix = "histo_"+bkgNames[i] + "_" + sysName;
                if(sysName.Contains("LepScale") || sysName.Contains("LepRes")) // FIXME sysName -> sysPostfix
                {
                    histPostfix = "histo_"+bkgNames[i];
                    dirName_ = dirName;
                }
                if(bkgTypes[i] == "Fake")
                {
                    dirName_ = "Detector";
                    histPostfix = "histo_"+bkgNames[i];
                }

                htemp = getRawHist(var, filePath, dirName_, histPostfix, "h"+bkgNames[i], steering, useAxis, divBinWidth);
            }
            isFirstBkg = false;
            nthBkg++;
        }
        else
        {
            if(sysName == "")
            {
                htemp->Add(getRawHist(var, filePath, dirName_, histName_, "h"+bkgNames[i], steering, useAxis, divBinWidth));
            }
            else
            {
                TString histPostfix = "histo_"+bkgNames[i] + "_" + sysName;
                if(sysName.Contains("LepScale") || sysName.Contains("LepRes"))
                {
                    histPostfix = "histo_"+bkgNames[i];
                    dirName_ = dirName;
                }
                if(bkgTypes[i] == "Fake")
                {
                    dirName_ = "Detector";
                    histPostfix = "histo_"+bkgNames[i];
                }
                htemp->Add(getRawHist(var, filePath, dirName_, histPostfix, "h"+bkgNames[i], steering, useAxis, divBinWidth));
            }
            nthBkg++;
        }

        // This type of backgrounds all added, so add them to THStack
        if(nthBkg == bkgTypeN[bkgTypes[i]])
        {
            //cout << bkgTypes[i] << " " << bkgTypeN[bkgTypes[i]] << endl;
            htemp->SetFillColor(bkgColors[bkgTypes[i]]);
            htemp->SetLineColor(kWhite);
            if(doPtCombine)
            {
                TString massString = steering[15]; 
                int iMassBin = massString.Atoi();
                TH1* htempTotal = cloneEmptyHist(&hMCtotal, "hMCtotal_"+bkgTypes[i]+"_" + massString); 
                for(int ibin = 1; ibin < htemp->GetNbinsX() + 1; ibin++)
                {
                    htempTotal->SetBinContent(ibin + 17 * iMassBin, htemp->GetBinContent(ibin));        
                }
                htempTotal->SetFillColor(bkgColors[bkgTypes[i]]);
                htempTotal->SetLineColor(kWhite);
                hs.Add(htempTotal);
                hMCtotal.Add(htempTotal);
            }
            else
            {
                hs.Add(htemp);
                hMCtotal.Add(htemp);
            }

            if(sysName == "")
            {
                tempMap_legend[bkgTypes[i]] = htemp;
                temp_bkgName.push_back(bkgTypes[i]);
            }

            isFirstBkg = true;
            nthBkg = 0;
        }
    }

    for(unsigned int i = temp_bkgName.size(); i > 0; i--)
    {
        if(updateLegend)
            leg.AddEntry(tempMap_legend[temp_bkgName.at(i-1)], temp_bkgName.at(i-1), "F");
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

    TString fullPath="./UnfoldedHist_"+channel_name+yearStr+".root";
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
                if(it->first == "iterEM" && (it->second).at(i) != "Nominal")
                {
                    iBest_pt=iterEMPtUnfold->ScanSURE(NITER_Iterative, &graph_SURE_IterativeSURE_pt, &graph_DFdeviance_IterativeSURE_pt);
                    iBest_mass=iterEMMassUnfold->ScanSURE(NITER_Iterative, &graph_SURE_IterativeSURE_mass, &graph_DFdeviance_IterativeSURE_mass);
                    cout << "iBest pt, Mass: " << iBest_pt << " " << iBest_mass << endl;
                    continue;
                }

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
    TFile f("./AcceptHist_"+channel_name+yearStr+".root","UPDATE"); 

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

    hFullPhaseMassData = nomMassUnfold->GetOutput("hMassFullPhaseNominal",0,0, "*[*]", false);
    hFullPhaseMassData->Multiply(hAcceptanceMass); // Bin by bin acceptance correction 

    f.cd();
    hFullPhaseMassData->Write();
    hFullPhaseMassMC->SetName("hMassDYMCNominal");
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
    hFullPhasePtData = nomPtUnfold->GetOutput("hPtFullPhaseNominal",0,0, "*[*]", false);
    hFullPhasePtData->Multiply(hAcceptancePt);

    f.cd();
    hFullPhasePtData->Write();
    hFullPhasePtMC->SetName("hPtDYMCNominal");
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
                if(it->first == "iterEM" && (it->second).at(i) != "Nominal")
                {
                    hSysFullPhaseMassData[it->first][(it->second).at(i)] = iterEMMassUnfold->GetOutput("hMassFullPhase_"+it->first+"_"+(it->second).at(i),0,0, "*[*]", false);
                    hSysFullPhasePtData[it->first][(it->second).at(i)]   = iterEMPtUnfold->GetOutput("hPtFullPhase_"+it->first+"_"+(it->second).at(i),0,0, "*[*]", false);
                }
                else
                {
                    hSysFullPhaseMassData[it->first][(it->second).at(i)] = sysMassUnfold[it->first][(it->second).at(i)]->GetOutput("hMassFullPhase_"+it->first+"_"+(it->second).at(i),0,0, "*[*]", false);
                    hSysFullPhasePtData[it->first][(it->second).at(i)]   = sysPtUnfold[it->first][(it->second).at(i)]->GetOutput("hPtFullPhase_"+it->first+"_"+(it->second).at(i),0,0, "*[*]", false);
                }

                // Use different acceptance for PDF, AlphaS, Scale etc
                if( (((it->first).Contains("Scale") && !(it->first).Contains("Lep")) || (it->first).Contains("PDF") || (it->first).Contains("AlphaS")) && !(it->first).Contains("_") )
                {
                    // For mass
                    TH1* hFullPhaseMassMC_raw_sys = (TH1*) filein->Get("Acceptance/MassGen" + binDef + "/histo_DYJets_"+(it->second).at(i));
                    if(year==2016)
                        hFullPhaseMassMC_raw_sys->Add((TH1*) filein->Get("Acceptance/MassGen" + binDef + "/histo_DYJets10to50_"+(it->second).at(i)));
                    else
                        hFullPhaseMassMC_raw_sys->Add((TH1*) filein->Get("Acceptance/MassGen" + binDef + "/histo_DYJets10to50_MG_"+(it->second).at(i)));

                    TH1* hFiducialPhaseMassMC_sys = sysMassUnfold[it->first][(it->second).at(i)]->GetBias("hFiducialMass_sys", 0, 0, "*[*]", false);
                    TH1* hAcceptanceMass_sys = (TH1*) hFullPhaseMassMC_raw_sys->Clone("hAcceptanceMass_sys");
                    hAcceptanceMass_sys->Divide(hFiducialPhaseMassMC_sys);

                    TH1* hAcceptanceFractionMass_sys = (TH1*) hFiducialPhaseMassMC_sys->Clone("hAcceptanceFractionMass_sys");
                    hAcceptanceFractionMass_sys->Divide(hFullPhaseMassMC_raw_sys);
                    
                    hSysFullPhaseMassData[accepCorrOrEffCorr + "_" + it->first][(it->second).at(i)] = nomMassUnfold->GetOutput("hAcceptMassData" +it->first+(it->second).at(i),0,0, "*[*]", false);
                    hSysFullPhaseMassData[accepCorrOrEffCorr + "_" + it->first][(it->second).at(i)]->Multiply(hAcceptanceMass_sys);
                    hSysFullPhaseMassMC[accepCorrOrEffCorr + "_" + it->first][(it->second).at(i)] = hFullPhaseMassMC_raw_sys;
                    sysMapForAcceptance[accepCorrOrEffCorr + "_" + it->first].push_back((it->second).at(i)); // Update sysMapForAcceptance 

                    //hSysAcceptanceMass[it->first][(it->second).at(i)] = (TH1*) hAcceptanceMass_sys->Clone("Mass_" + it->first + "_" + (it->second).at(i));
                    hSysAcceptanceFractionMass[it->first][(it->second).at(i)] = (TH1*) hAcceptanceFractionMass_sys->Clone("FractionMass_" + it->first + "_" + (it->second).at(i));
                    delete hAcceptanceMass_sys;
                    delete hAcceptanceFractionMass_sys;

                    // For pt
                    TH1* hFullPhasePtMC_raw_sys = (TH1*) filein->Get("Acceptance/PtGen" + binDef + "/histo_DYJets_"+(it->second).at(i));
                    if(year==2016)
                        hFullPhasePtMC_raw_sys->Add((TH1*) filein->Get("Acceptance/PtGen" + binDef + "/histo_DYJets10to50_"+(it->second).at(i)));
                    else
                        hFullPhasePtMC_raw_sys->Add((TH1*) filein->Get("Acceptance/PtGen" + binDef + "/histo_DYJets10to50_MG_"+(it->second).at(i)));

                    TH1* hFiducialPhasePtMC_sys = sysPtUnfold[it->first][(it->second).at(i)]->GetBias("hFiducialPt_sys", 0, 0, "*[*]", false);

                    TH1* hAcceptancePt_sys = (TH1*) hFullPhasePtMC_raw_sys->Clone("hAcceptancePt_sys");
                    hAcceptancePt_sys->Divide(hFiducialPhasePtMC_sys);

                    TH1* hAcceptanceFractionPt_sys = (TH1*) hFiducialPhasePtMC_sys->Clone("hAcceptanceFractionPt_sys");
                    hAcceptanceFractionPt_sys->Divide(hFullPhasePtMC_raw_sys);

                    hSysFullPhasePtData[accepCorrOrEffCorr + "_" + it->first][(it->second).at(i)] = nomPtUnfold->GetOutput("hAcceptPtData" +it->first+(it->second).at(i),0,0, "*[*]", false);
                    hSysFullPhasePtData[accepCorrOrEffCorr + "_" + it->first][(it->second).at(i)]->Multiply(hAcceptancePt_sys);
                    hSysFullPhasePtMC[accepCorrOrEffCorr + "_" + it->first][(it->second).at(i)] = hFullPhasePtMC_raw_sys;

                    if((it->first).Contains("PDF"))
                    {
                        fillMassPDFVariationHist_Accept(i+1);
                        fillPtPDFVariationHist_Accept(i+1);
                    }

                    //hSysAcceptancePt[it->first][(it->second).at(i)] = (TH1*) hAcceptancePt_sys->Clone("Pt_" + it->first + "_" + (it->second).at(i));
                    hSysAcceptanceFractionPt[it->first][(it->second).at(i)] = (TH1*) hAcceptanceFractionPt_sys->Clone("FractionPt_" + it->first + "_" + (it->second).at(i));
                    delete hAcceptancePt_sys;
                    delete hAcceptanceFractionPt_sys;

                    hSysFullPhaseMassData[it->first][(it->second).at(i)]->Multiply(hAcceptanceMass);
                    hSysFullPhasePtData[it->first][(it->second).at(i)]->Multiply(hAcceptancePt);
                    hSysFullPhaseMassMC[it->first][(it->second).at(i)] = hFullPhaseMassMC_raw_sys;
                    hSysFullPhasePtMC[it->first][(it->second).at(i)]   = hFullPhasePtMC_raw_sys;

                    f.cd();
                    hSysFullPhaseMassData[it->first][(it->second).at(i)]->Write();
                    hSysFullPhasePtData[it->first][(it->second).at(i)]->Write();
                    hFullPhaseMassMC_raw_sys->SetName("hMassDYMC_" +it->first+"_"+(it->second).at(i));
                    hFullPhasePtMC_raw_sys->SetName("hPtDYMC_" +it->first+"_"+(it->second).at(i));
                    hFullPhaseMassMC_raw_sys->Write();
                    hFullPhasePtMC_raw_sys->Write();

                    hSysAcceptanceMass[it->first][(it->second).at(i)] = (TH1*) hAcceptanceMass->Clone("Mass_" + it->first + "_" + (it->second).at(i));
                    hSysAcceptancePt[it->first][(it->second).at(i)] = (TH1*) hFullPhasePtMC->Clone("Pt_" + it->first + "_" + (it->second).at(i));
                }
                else
                {
                    // Use nominal acceptance
                    hSysFullPhaseMassData[it->first][(it->second).at(i)]->Multiply(hAcceptanceMass);
                    hSysFullPhasePtData[it->first][(it->second).at(i)]->Multiply(hAcceptancePt);
                    hSysFullPhaseMassMC[it->first][(it->second).at(i)] = hFullPhaseMassMC;
                    hSysFullPhasePtMC[it->first][(it->second).at(i)]   = hFullPhasePtMC;

                    hSysAcceptanceMass[it->first][(it->second).at(i)] = (TH1*) hAcceptanceMass->Clone("Mass_" + it->first + "_" + (it->second).at(i));
                    hSysAcceptancePt[it->first][(it->second).at(i)] = (TH1*) hFullPhasePtMC->Clone("Pt_" + it->first + "_" + (it->second).at(i));
                }
            }
            it++;
        }

        // Update sys map
        // Sys. unc. of acceptance correction
        // PDF, alphaS, Scale
        it = sysMapForAcceptance.begin();
        while(it != sysMapForAcceptance.end())
        {
            int size = (it->second).size();
            for(int i = 0; i < size; i++)
            {
               sysMap[it->first].push_back((it->second).at(i)); 
            }
            it++;
        }
  
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

    f.cd();
    f.Close();

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

void ISRUnfold::drawAcceptance(TString var, TH1* hMC, TString outName)
{
    setTDRStyle();
    writeExtraText = true;
    extraText  = "Simulation";
    gStyle->SetLineWidth(globalLinedWidth);
    gStyle->SetFrameLineWidth(globalLinedWidth);
    gROOT->ForceStyle();

    TH1* hFullPhase = NULL;
    TH1* hFiducialPhase = NULL;
    TH1* hAcceptance = NULL;

    if(var.Contains("Mass"))
    {
        hFullPhase = mass_binning_Gen->ExtractHistogram("hData", hFullPhaseMassMC, 0, true, "mass[UO];pt[UOC0]");
        hFiducialPhase = mass_binning_Gen->ExtractHistogram("hData", hMC, 0, true, "mass[UO];pt[UOC0]");

        divideByBinWidth(hFullPhase, false);
        divideByBinWidth(hFiducialPhase, false);

        hAcceptance = (TH1*) hFiducialPhase->Clone("hAcceptance");
        hAcceptance->Divide(hFullPhase);
        drawComparisonPlot(var, "Acceptance", "Events/bin", "Acceptance", "", hFullPhase, hFiducialPhase, hAcceptance, outName);
    }
    if(var.Contains("Pt"))
    {
        const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
        int nMassBin = temp_tvecd->GetNrows() - 1;

        for(int ibin = 0; ibin < nMassBin; ibin++)
        {
            TString ibinMass;
            ibinMass.Form("%d", ibin);

            hFullPhase = pt_binning_Gen->ExtractHistogram("hData", hFullPhasePtMC, 0, true, "pt[UO];mass[UOC"+ibinMass+"]");
            hFiducialPhase = pt_binning_Gen->ExtractHistogram("hData", hMC, 0, true, "pt[UO];mass[UOC"+ibinMass+"]");

            divideByBinWidth(hFullPhase, false);
            divideByBinWidth(hFiducialPhase, false);

            hAcceptance = (TH1*) hFiducialPhase->Clone("hAcceptance");
            hAcceptance->Divide(hFullPhase);
            drawComparisonPlot(var, "Acceptance", "Events/bin", "Acceptance", "", hFullPhase, hFiducialPhase, hAcceptance, "AcceptancePt_M_"+ibinMass+"_"+outName, ibin);

            delete hFullPhase;
            delete hFiducialPhase;
            delete hAcceptance;
        }
    }
}

void ISRUnfold::drawComparisonPlot(TString var, TString plotName, TString topYaxisName, TString bottomYaxisName, TString bottomXaxisName, TH1* h1, TH1* h2, TH1* hratio, TString outName, int nthMassBin)
{

    // TODO make below lines as a function
    TCanvas* c_out = new TCanvas(plotName + "_" + var, plotName + "_" + var, 3200, 2800);
    c_out->Draw();
    c_out->cd();

    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetBottomMargin(0.01);
    pad1->SetTopMargin(0.1);
    pad1->SetTicks(1);
    pad1->SetLogy();
    if(var.Contains("Mass"))
        pad1->SetLogx();
    pad1->Draw();
    pad1->cd();

    h1->SetTitle("");
    h1->SetStats(false);
    h1->GetXaxis()->SetMoreLogLabels(true);
    h1->Draw("p9e");
    h1->SetMarkerStyle(20);
    h1->SetMarkerSize(3.2);
    h1->SetLineColor(kBlack);
    h1->GetYaxis()->SetTitle(topYaxisName);
    h1->SetMaximum(5e9);
    h1->SetMinimum(2e-1);

    h2->SetMarkerStyle(20);
    h2->SetMarkerSize(3.2);
    h2->SetLineColor(kRed);
    h2->SetMarkerColor(kRed);
    h2->Draw("p9e same");

    int iPeriod_ = 4;
    if(year == 2017)
        iPeriod_ = 5;
    if(year == 2018)
        iPeriod_ = 6;
    CMS_lumi(pad1, iPeriod_, 11);


    if(var.Contains("Pt"))
        writeCutInfo(pad1, var, nthMassBin);

    c_out->cd();

    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.35);
    if(var=="Mass")
    {
        pad2->SetLogx();
        hratio->GetXaxis()->SetMoreLogLabels(true);
    }
    pad2->SetTicks(1);
    pad2->SetGridy(1);
    //if(var.Contains("Mass"))
    //    pad2->SetLogy();
    pad2->Draw();
    pad2->cd();

    hratio->SetStats(false);
    hratio->Draw("p9e");
    hratio->SetMarkerStyle(20);
    hratio->SetMarkerSize(3.2);
    hratio->SetLineColor(kBlack);
    hratio->GetYaxis()->SetTitle(bottomYaxisName);
    bottomXaxisName = "";

    hratio->SetMinimum(1e-3);
    hratio->SetMaximum(1.0);

    setXaxisTitle(hratio, var, true);
    hratio->Draw("p9e same");

    c_out->cd();
    c_out->SaveAs(outName!=""?output_baseDir+plotName+"_"+var+"_"+outName+".pdf":output_baseDir+plotName+"_"+var+".pdf");
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
            if( (it->first).Contains("iterEM") && !((it->second).at(i)).Contains("Nominal"))
            {
                hunfolded_mass = iterEMMassUnfold->GetOutput("hunfolded_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);
            }
            else
            {
                p_unfold = sysMassUnfold[it->first][(it->second).at(i)];
                hunfolded_mass = p_unfold->GetOutput("hunfolded_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);
            }


            // Loop over mass bins
            for(int ibin = 0; ibin < nMassBin; ibin++)
            {
                // Set x-axis range
                hunfolded_mass->GetXaxis()->SetRange(hunfolded_mass->GetXaxis()->FindBin(massBins[ibin]+0.01), hunfolded_mass->GetXaxis()->FindBin(massBins[ibin+1]-0.01));
                meanMass_data_unfolded_sysVariation[it->first][(it->second).at(i)].push_back(hunfolded_mass->GetMean());
                //cout << it->first << " " << (it->second).at(i) << " " << hunfolded_mass->GetMean() << endl;
            }// End of mass bin loop
            delete hunfolded_mass;
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
                //cout << it->first << " " << (it->second).at(i) << " " << hunfolded_mass->GetMean() << endl;
            }// End of mass bin loop
            delete hunfolded_mass;
        }
        it++;
    }
}

void ISRUnfold::setSysError()
{
    ofstream output;
    output.open (output_baseDir + "meanValuesAfterUnfolding.txt");
    output << "------------------------------------- Systematic Uncertainty for ISR analysis ---------------------------------------" << endl;
    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;

    std::streamsize ss = output.precision();
    output.precision(2);
    output.setf( std::ios::fixed, std:: ios::floatfield );

    // For systematic
    std::map<TString, std::vector<TString>>::iterator it = sysMap.begin();
    while(it != sysMap.end())
    {
        output << "Systematic: " << it->first << endl;
        output << setw(10) << "Mass bin" ;
        output << setw(10) << "Mass" ;
        output << setw(10) << "Pt" << endl;
        for(int ibin = 0; ibin < nMassBin; ibin++)
        {
            double error_mass = 0.;
            double error_pt = 0.;
            int size = (it->second).size();
            output << setw(10) << ibin ;


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

            output << setw(10) << error_mass ;
            output << setw(10) << error_pt << endl;
        }
        output << endl;
        it++;
    }
    output.precision(ss);
    output.close();
}

void ISRUnfold::setSysError_Accept()
{
    ofstream output;
    output.open (output_baseDir + "meanValuesAfterAcceptanceCorrection.txt");
    output << "------------------------------------- Systematic Uncertainty for ISR analysis (Acceptance corrected)---------------------------------------" << endl;
    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;

    std::streamsize ss = output.precision();
    output.precision(2);
    output.setf( std::ios::fixed, std:: ios::floatfield );

    // For systematic
    std::map<TString, std::vector<TString>>::iterator it = sysMap.begin();
    while(it != sysMap.end())
    {
        output << "Systematic: " << it->first << endl;
        output << setw(10) << "Mass bin" ;
        output << setw(10) << "Mass" ;
        output << setw(10) << "Pt" << endl;
        //output << "Systematic: " << it->first << endl;
        for(int ibin = 0; ibin < nMassBin; ibin++)
        {
            double error_mass = 0.;
            double error_pt = 0.;
            int size = (it->second).size();
            output << setw(10) << ibin ;

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

            //output << ibin << " mass bin " << endl;
            //output << "mass: " << error_mass / meanMass_data_acc_corrected.at(ibin) * 100.  << " pt: " << error_pt / meanPt_data_acc_corrected.at(ibin) * 100.<< endl;
            output << setw(10) << error_mass ;
            output << setw(10) << error_pt << endl;
        }
        it++;
        output << endl;
    }

    setAcceptError();

    output.precision(ss);
    output.close();
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

double ISRUnfold::getDetMeanMass(int ibin)
{

    int size = meanMass_data_folded.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanMass_data_folded.at(ibin);
    }
}

double ISRUnfold::getDetMeanMassError(int ibin)
{

    int size = meanMass_data_folded.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanMassStatErr_data_folded.at(ibin);
    }
}

double ISRUnfold::getMCDetMeanMass(int ibin)
{

    int size = meanMass_data_folded.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanMass_mc_folded.at(ibin);
    }
}

double ISRUnfold::getMCDetMeanMassError(int ibin)
{

    int size = meanMass_data_folded.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanMassStatErr_mc_folded.at(ibin);
    }
}

double ISRUnfold::getUnfMeanMass(int ibin)
{

    int size = meanMass_data_unfolded.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanMass_data_unfolded.at(ibin);
    }
}

double ISRUnfold::getUnfMeanMassError(int ibin)
{

    int size = meanMassStatErr_data_unfolded.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanMassStatErr_data_unfolded.at(ibin);
    }
}

double ISRUnfold::getUnfMeanMassSysError(int ibin)
{

    int size = meanMassStatErr_data_unfolded.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanMassSysErr_data_unfolded.at(ibin);
    }
}

double ISRUnfold::getAccMeanMassSysError(int ibin)
{

    int size = meanMassStatErr_data_unfolded.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanMassSysErr_data_acc_corrected.at(ibin);
    }
}

double ISRUnfold::getAccMeanMassTotError(int ibin)
{

    int size = meanMassStatErr_data_unfolded.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanMassTotErr_data_acc_corrected.at(ibin);
    }
}

double ISRUnfold::getAccMeanMass(int ibin)
{

    int size = meanMass_data_unfolded.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanMass_data_acc_corrected.at(ibin);
    }
}

double ISRUnfold::getAccMeanMassError(int ibin)
{

    int size = meanMassStatErr_data_unfolded.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanMassStatErr_data_acc_corrected.at(ibin);
    }
}

double ISRUnfold::getMCGenMeanMass(int ibin)
{

    int size = meanMass_mc_unfolded.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanMass_mc_unfolded.at(ibin);
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
            if((it->first).Contains("iterEM") && !((it->second).at(i)).Contains("Nominal"))
            {
                p_unfold = NULL;
            }
            else
            {
                p_unfold = sysPtUnfold[it->first][(it->second).at(i)];
            }
            // Save mean pt
            for(int j = 0; j < nMassBin; j++)
            {
                TString ibinMass;
                ibinMass.Form("%d", j);

                TH1* hunfolded_pt = NULL;

                // Get histograms to set mean values
                if((it->first).Contains("iterEM") && !((it->second).at(i)).Contains("Nominal"))
                {
                    hunfolded_pt = iterEMPtUnfold->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
                }
                else
                {
                    hunfolded_pt = p_unfold->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
                }
                meanPt_data_unfolded_sysVariation[it->first][(it->second).at(i)].push_back(hunfolded_pt->GetMean());
                //cout << it->first << " " << (it->second).at(i) << " " << hunfolded_pt->GetMean() << endl;
                delete hunfolded_pt;
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

double ISRUnfold::getDetMeanPt(int ibin)
{

    int size = meanPt_data_folded.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanPt_data_folded.at(ibin);
    }
}

double ISRUnfold::getDetMeanPtError(int ibin)
{

    int size = meanPt_data_folded.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanPtStatErr_data_folded.at(ibin);
    }
}

double ISRUnfold::getMCFullPhaseMeanPt(int ibin, TString sysName)
{

    int size = meanPt_data_folded.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        if(sysName == "")
        {
            return meanPt_theory_acc_corrected.at(ibin);
        }
        else if(sysName == "stat")
        {
            return meanPtStatErr_theory_acc_corrected.at(ibin);
        }
        else
        {
            return meanPt_theory_accept_systematic[sysName].at(ibin);
        }
    }
}

double ISRUnfold::getMCFullPhaseMeanMass(int ibin, TString sysName)
{

    int size = meanMass_data_folded.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        if(sysName == "")
        {
            return meanMass_theory_acc_corrected.at(ibin);
        }
        else if(sysName == "stat")
        {
            return meanMassStatErr_theory_acc_corrected.at(ibin);
        }
        else
        {
            return meanMass_theory_accept_systematic[sysName].at(ibin);
        }
    }
}

double ISRUnfold::getMCDetMeanPt(int ibin)
{

    int size = meanPt_data_folded.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanPt_mc_folded.at(ibin);
    }
}

double ISRUnfold::getMCDetMeanPtError(int ibin)
{

    int size = meanPt_data_folded.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanPtStatErr_mc_folded.at(ibin);
    }
}


double ISRUnfold::getUnfMeanPt(int ibin)
{

    int size = meanPt_data_unfolded.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanPt_data_unfolded.at(ibin);
    }
}

double ISRUnfold::getUnfMeanPtError(int ibin)
{

    int size = meanPtStatErr_data_unfolded.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanPtStatErr_data_unfolded.at(ibin);
    }
}

double ISRUnfold::getUnfMeanPtSysError(int ibin)
{

    int size = meanPtStatErr_data_unfolded.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanPtSysErr_data_unfolded.at(ibin);
    }
}

double ISRUnfold::getAccMeanPtSysError(int ibin)
{

    int size = meanPtStatErr_data_unfolded.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanPtSysErr_data_acc_corrected.at(ibin);
    }
}

double ISRUnfold::getAccMeanPtTotError(int ibin)
{

    int size = meanPtStatErr_data_unfolded.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanPtTotErr_data_acc_corrected.at(ibin);
    }
}

double ISRUnfold::getAccMeanPt(int ibin)
{

    int size = meanPt_data_unfolded.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanPt_data_acc_corrected.at(ibin);
    }
}

double ISRUnfold::getAccMeanPtError(int ibin)
{

    int size = meanPtStatErr_data_unfolded.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanPtStatErr_data_acc_corrected.at(ibin);
    }
}

double ISRUnfold::getMCGenMeanPt(int ibin)
{

    int size = meanPt_mc_unfolded.size();
    if(ibin >= size)
    {
        exit (EXIT_FAILURE);
    }
    else
    {
        return meanPt_mc_unfolded.at(ibin);
    }
}

// https://root.cern.ch/root/html/tutorials/graphs/graphtext.C.html
void ISRUnfold::drawtext(TGraph *g)
{
    Int_t i,n;
    Double_t x,y;
    TLatex *l;

    n = g->GetN();
    for (i=0; i<n; i++)
    {
       g->GetPoint(i,x,y);
       l = new TLatex(x-0.1,y+0.02,Form("%d",i+1));
       l->SetTextSize(0.015);
       l->SetTextFont(42);
       l->SetTextAlign(21);
       l->Draw();
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

TH1* ISRUnfold::getGenMCHist(TString var, TString steering, bool useAxis, int massBin, bool binWidth)
{
    TH1* outHist = NULL;
    if(var == "Mass")
    {
        outHist = nomMassUnfold->GetBias("GenMass",0,0,steering,useAxis);

        if(binWidth)
        {
            divideByBinWidth(outHist, false);
        }
    }
    else
    {
        TString nth;
        nth.Form("%d", massBin);
        outHist =  nomPtUnfold->GetBias("GenPt_"+nth,0,0,steering,useAxis);

        if(binWidth)
        {
            divideByBinWidth(outHist, false);
        }
    }

    return outHist;
}

TH1* ISRUnfold::getDetHists(TString var, TString outHistName, TString steering, bool useAxis)
{
    // Return background subtracted data

    if(var == "Mass")
        return nomMassUnfold->GetInput(outHistName,0,0,steering,useAxis);
    else
        return nomPtUnfold->GetInput(outHistName,0,0,steering,useAxis);
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

void ISRUnfold::drawStatVariation(bool isPt, int massBin)
{
    gStyle->SetOptStat(1101);
    TCanvas* c = new TCanvas("c","c", 2400, 2400);
    //c->SetLogy();
    c->cd();

    TString nth;
    nth.Form("%d", massBin);
    TString year_;
    year_.Form("%d", year);

    if(isPt)
    {
        meanPtStatVariation.at(massBin)->Draw("pe");
        c->Update();

        double y_max = c->GetFrame()->GetY2();
        double y_min = c->GetFrame()->GetY1();
        double x_nominal = meanPt_data_unfolded.at(massBin);

        TLine grid;
        grid.SetLineColor(kRed);
        grid.SetLineStyle(2);
        grid.DrawLine(x_nominal, y_min, x_nominal, y_max);

        c->SaveAs(output_baseDir+"MeanPtStat_" + nth + year_ + ".pdf");
    }
    else
    {
        meanMassStatVariation.at(massBin)->Draw("pe");
        c->Update();

        double y_max = c->GetFrame()->GetY2();
        double y_min = c->GetFrame()->GetY1();
        double x_nominal = meanMass_data_unfolded.at(massBin);

        TLine grid;
        grid.SetLineColor(kRed);
        grid.SetLineStyle(2);
        grid.DrawLine(x_nominal, y_min, x_nominal, y_max);

        c->SaveAs(output_baseDir+"MeanMassStat_" + nth + year_ + ".pdf");
    }
    delete c;
}


void ISRUnfold::drawPDFVariation(bool isPt, int massBin)
{
    gStyle->SetOptStat(1101);
    TCanvas* c = new TCanvas("c","c", 2400, 2400);
    //c->SetLogy();
    c->cd();

    TString nth;
    nth.Form("%d", massBin);
    TString year_;
    year_.Form("%d", year);

    if(isPt)
    {
        meanPtPDFVariation.at(massBin)->Draw("pe");
        c->Update();

        double y_max = c->GetFrame()->GetY2();
        double y_min = c->GetFrame()->GetY1();
        double x_nominal = meanPt_data_unfolded.at(massBin);

        TLine grid;
        grid.SetLineColor(kRed);
        grid.SetLineStyle(2);
        grid.DrawLine(x_nominal, y_min, x_nominal, y_max);

        c->SaveAs(output_baseDir+"MeanPtPDF_" + nth + year_ + ".pdf");
    }
    else
    {
        meanMassPDFVariation.at(massBin)->Draw("pe");
        c->Update();

        double y_max = c->GetFrame()->GetY2();
        double y_min = c->GetFrame()->GetY1();
        double x_nominal = meanMass_data_unfolded.at(massBin);

        TLine grid;
        grid.SetLineColor(kRed);
        grid.SetLineStyle(2);
        grid.DrawLine(x_nominal, y_min, x_nominal, y_max);

        c->SaveAs(output_baseDir+"MeanMassPDF_" + nth + year_ + ".pdf");
    }
    delete c;
}

void ISRUnfold::drawSysVariation(TString sysName, TString var, int massBin)
{
    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("c","c", 2400, 2400);
    c->cd();

    TH1* hout = NULL;

    TString nth;
    nth.Form("%d", massBin);
    TString year_;
    year_.Form("%d", year);

    int size = sysMap[sysName].size();

    if(var.Contains("Pt"))
    {
        double x_nominal = meanPt_data_unfolded.at(massBin);
        double x_err = meanPtSysErr_data_unfolded.at(massBin);
        hout = new TH1D("hSysVar", "hSysVar", 10, x_nominal-x_err, x_nominal+x_err); // # of bins doesn't matter here
        hout->GetYaxis()->SetTickLength(0);
        hout->GetYaxis()->SetLabelSize(0);
        hout->GetXaxis()->SetTitle("<p_{T}^{\ell\ell}> [GeV]");
        hout->Draw();
        c->Update();

        double y_max = c->GetFrame()->GetY2();
        double y_min = c->GetFrame()->GetY1();

        TLatex *l = new TLatex(x_nominal,y_max/2., "Nominal");
        l->SetTextAngle(90);
        l->SetTextFont(43);
        l->SetTextSize(50);
        l->Draw();

        TLine grid;
        grid.SetLineColor(kRed);
        grid.SetLineStyle(1);
        grid.DrawLine(x_nominal, y_min, x_nominal, y_max);

        TLine grid_var;
        grid_var.SetLineColor(kBlack);
        grid_var.SetLineStyle(2);

        for(int i = 0; i < size; i++)
        {
            double x_variation = meanPt_data_unfolded_sysVariation[sysName][sysMap[sysName][i]].at(massBin);
            grid_var.DrawLine(x_variation, y_min, x_variation, y_max);
            l = new TLatex(x_variation,y_max/2., sysMap[sysName][i]);
            l->SetTextAngle(90);
            l->SetTextFont(43);
            l->SetTextSize(50);
            l->Draw();
        }

        c->SaveAs(output_baseDir+"MeanPt_"+sysName+"_"+nth+"_"+year_+".pdf");
    }
    else
    {
        //
        double x_nominal = meanMass_data_unfolded.at(massBin);
        double x_err = meanMassSysErr_data_unfolded.at(massBin);
        hout = new TH1D("hSysVar", "hSysVar", 10, x_nominal-x_err, x_nominal+x_err); // # of bins doesn't matter here
        hout->GetYaxis()->SetTickLength(0);
        hout->GetYaxis()->SetLabelSize(0);
        hout->GetXaxis()->SetTitle("<Mass^{\ell\ell}> [GeV]");
        hout->Draw();
        c->Update();

        double y_max = c->GetFrame()->GetY2();
        double y_min = c->GetFrame()->GetY1();

        TLatex *l = new TLatex(x_nominal,y_max/2., "Nominal");
        l->SetTextAngle(90);
        l->SetTextFont(43);
        l->SetTextSize(50);
        l->Draw();

        TLine grid;
        grid.SetLineColor(kRed);
        grid.SetLineStyle(1);
        grid.DrawLine(x_nominal, y_min, x_nominal, y_max);

        TLine grid_var;
        grid_var.SetLineColor(kBlack);
        grid_var.SetLineStyle(2);

        for(int i = 0; i < size; i++)
        {
            double x_variation = meanMass_data_unfolded_sysVariation[sysName][sysMap[sysName][i]].at(massBin);
            grid_var.DrawLine(x_variation, y_min, x_variation, y_max);
            l = new TLatex(x_variation,y_max/2., sysMap[sysName][i]);
            l->SetTextAngle(90);
            l->SetTextFont(43);
            l->SetTextSize(50);
            l->Draw();
        }

        c->SaveAs(output_baseDir+"MeanMass_"+sysName+"_"+nth+"_"+year_+".pdf");
    }
}

void ISRUnfold::drawSysVariation_Accept(TString sysName, TString var, int massBin)
{
    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("c","c", 2400, 2400);
    c->cd();

    TH1* hout = NULL;

    TString nth;
    nth.Form("%d", massBin);
    TString year_;
    year_.Form("%d", year);

    int size = sysMap[sysName].size();

    if(var.Contains("Pt"))
    {
        double x_nominal = meanPt_data_acc_corrected.at(massBin);
        double x_err = meanPtSysErr_data_acc_corrected.at(massBin);
        hout = new TH1D("hSysVar", "hSysVar", 10, x_nominal-x_err, x_nominal+x_err); // # of bins doesn't matter here
        hout->GetYaxis()->SetTickLength(0);
        hout->GetYaxis()->SetLabelSize(0);
        hout->GetXaxis()->SetTitle("<p_{T}^{\ell\ell}> [GeV]");
        hout->Draw();
        c->Update();

        double y_max = c->GetFrame()->GetY2();
        double y_min = c->GetFrame()->GetY1();

        TLatex *l = new TLatex(x_nominal,y_max/2., "Nominal");
        l->SetTextAngle(90);
        l->SetTextFont(43);
        l->SetTextSize(50);
        l->Draw();

        TLine grid;
        grid.SetLineColor(kRed);
        grid.SetLineStyle(1);
        grid.DrawLine(x_nominal, y_min, x_nominal, y_max);

        TLine grid_var;
        grid_var.SetLineColor(kBlack);
        grid_var.SetLineStyle(2);

        for(int i = 0; i < size; i++)
        {
            double x_variation = meanPt_data_accept_sysVariation[sysName][sysMap[sysName][i]].at(massBin);
            grid_var.DrawLine(x_variation, y_min, x_variation, y_max);
            l = new TLatex(x_variation,y_max/2., sysMap[sysName][i]);
            l->SetTextAngle(90);
            l->SetTextFont(43);
            l->SetTextSize(50);
            l->Draw();
        }

        c->SaveAs(output_baseDir+"AcceptMeanPt_"+sysName+"_"+nth+"_"+year_+".pdf");
    }
    else
    {
        //
        double x_nominal = meanMass_data_acc_corrected.at(massBin);
        double x_err = meanMassSysErr_data_acc_corrected.at(massBin);
        hout = new TH1D("hSysVar", "hSysVar", 10, x_nominal-x_err, x_nominal+x_err); // # of bins doesn't matter here
        hout->GetYaxis()->SetTickLength(0);
        hout->GetYaxis()->SetLabelSize(0);
        hout->GetXaxis()->SetTitle("<Mass^{\ell\ell}> [GeV]");
        hout->Draw();
        c->Update();

        double y_max = c->GetFrame()->GetY2();
        double y_min = c->GetFrame()->GetY1();

        TLatex *l = new TLatex(x_nominal,y_max/2., "Nominal");
        l->SetTextAngle(90);
        l->SetTextFont(43);
        l->SetTextSize(50);
        l->Draw();

        TLine grid;
        grid.SetLineColor(kRed);
        grid.SetLineStyle(1);
        grid.DrawLine(x_nominal, y_min, x_nominal, y_max);

        TLine grid_var;
        grid_var.SetLineColor(kBlack);
        grid_var.SetLineStyle(2);

        for(int i = 0; i < size; i++)
        {
            double x_variation = meanMass_data_accept_sysVariation[sysName][sysMap[sysName][i]].at(massBin);
            grid_var.DrawLine(x_variation, y_min, x_variation, y_max);
            l = new TLatex(x_variation,y_max/2., sysMap[sysName][i]);
            l->SetTextAngle(90);
            l->SetTextFont(43);
            l->SetTextSize(50);
            l->Draw();
        }

        c->SaveAs(output_baseDir+"AcceptMeanMass_"+sysName+"_"+nth+"_"+year_+".pdf");
    }
}


void ISRUnfold::drawSystematics(TString var, bool isHistStye)
{

    // Loop over
    // meanMass_data_folded_systematic, meanPt_data_folded_systematic
    setTDRStyle();
    writeExtraText = true;
    extraText  = "Work in progress";
    gStyle->SetLineWidth(3);
    gStyle->SetFrameLineWidth(3);
    gROOT->ForceStyle();

    // Prepare TGraph for each systematic
    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;
    vector<double> v_bins;
    vector<double> v_bins_error;
    vector<double> v_dummy_error;
    for(int i = 0; i < nMassBin; i++)
    {
        v_bins.push_back(i+1);
        v_bins_error.push_back(0.5);
        v_dummy_error.push_back(0.0);
    }

    map<TString, TGraph*> map_sys_graph;
    map<TString, TGraphErrors*> map_sys_grapherrors;
    std::map<TString, std::vector<TString>>::iterator it = sysMap.begin();
    while(it != sysMap.end())
    {
        if((it->first).Contains("Acceptance") || (it->first).Contains("Efficiency"))
        {
            it++;
            continue;
        }
        if(var.Contains("Pt"))
        {
            if(isHistStye == false)
                map_sys_graph[it->first] = new TGraph(nMassBin, &v_bins[0], &meanPt_data_folded_rel_systematic[it->first][0]);
            else
                map_sys_grapherrors[it->first] = new TGraphErrors(nMassBin, &v_bins[0], &meanPt_data_folded_rel_systematic[it->first][0], &v_bins_error[0], &v_dummy_error[0]);
        }
        else
        {
            if(isHistStye == false)
                map_sys_graph[it->first] = new TGraph(nMassBin, &v_bins[0], &meanMass_data_folded_rel_systematic[it->first][0]);
            else
                map_sys_grapherrors[it->first] = new TGraphErrors(nMassBin, &v_bins[0], &meanMass_data_folded_rel_systematic[it->first][0], &v_bins_error[0], &v_dummy_error[0]);
        }
        it++;
    }

    // Create canvas
    TCanvas* c_out = new TCanvas("relative_uncertainty_"+var, "relative_uncertainty_"+var, 3000, 1800);
    c_out->SetGridy(1);
    //c_out->SetGridx(1);
    c_out->SetTopMargin(0.08);
    c_out->Draw();
    c_out->cd();

    TLegend* leg = new TLegend(0.45, 0.5, 0.95, 0.9,"","brNDC");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(43);
    leg->SetTextSize(70);

    int marker = 20;
    int markerSize = 3;
    int markerColor = 1;
    bool first_draw = true;

    it = sysMap.begin();
    while(it != sysMap.end())
    {

        if((it->first).Contains("Acceptance") || (it->first).Contains("Efficiency"))
        {
            it++;
            continue;
        }

        if(markerColor == 10) markerColor = 1;
        if(marker == 24) marker = 29;
        if(marker == 30) marker = 33;
        if(marker == 35) marker = 39;
        if(marker == 40) marker = 41;
        if(marker == 42) marker = 47;
        if(marker > 47) marker = 20;

        if(isHistStye == false)
        {
            map_sys_graph[it->first]->SetMarkerStyle(marker);
            map_sys_graph[it->first]->SetMarkerSize(markerSize);
            map_sys_graph[it->first]->SetMarkerColor(markerColor==5?46:markerColor);
            map_sys_graph[it->first]->SetLineColor(markerColor==5?46:markerColor);
        }
        else
        {
            map_sys_grapherrors[it->first]->SetLineWidth(5);
            map_sys_grapherrors[it->first]->SetMarkerStyle(marker);
            map_sys_grapherrors[it->first]->SetMarkerSize(markerSize);
            map_sys_grapherrors[it->first]->SetMarkerColor(markerColor==5?46:markerColor);
            map_sys_grapherrors[it->first]->SetLineColor(markerColor==5?46:markerColor);
        }
        markerColor++;
        marker++;

        if(first_draw)
        {
            //map_sys_graph[it->first]->SetTitleOffset(0.02);
            if(isHistStye == false)
            {
                if(var.Contains("Pt"))
                    map_sys_graph[it->first]->SetTitle("Relative uncertainty on <p_{T}^{\ell\ell}> [%]");
                else
                    map_sys_graph[it->first]->SetTitle("Relative uncertainty on <Mass^{\ell\ell}> [%]");

                map_sys_graph[it->first]->GetYaxis()->SetRangeUser(0., 3.);
                map_sys_graph[it->first]->GetXaxis()->SetLimits(0.5,5.5);
                map_sys_graph[it->first]->GetYaxis()->SetTitleOffset(1.0);
                map_sys_graph[it->first]->GetYaxis()->SetTitle("Relative uncertainty [%]");
                map_sys_graph[it->first]->GetXaxis()->SetTitleOffset(0.7);
                map_sys_graph[it->first]->GetXaxis()->SetTitle("mass bin");
                map_sys_graph[it->first]->Draw("APC");
            }
            else
            {
                if(var.Contains("Pt"))
                    map_sys_grapherrors[it->first]->SetTitle("Relative uncertainty on <p_{T}^{\ell\ell}> [%]");
                else
                    map_sys_grapherrors[it->first]->SetTitle("Relative uncertainty on <Mass^{\ell\ell}> [%]");

                map_sys_grapherrors[it->first]->GetYaxis()->SetRangeUser(0., 3.);
                map_sys_grapherrors[it->first]->GetXaxis()->SetLimits(0.5,5.5);
                map_sys_grapherrors[it->first]->GetYaxis()->SetTitleOffset(1.0);
                map_sys_grapherrors[it->first]->GetYaxis()->SetTitle("Relative uncertainty [%]");
                map_sys_grapherrors[it->first]->GetXaxis()->SetTitleOffset(0.7);
                map_sys_grapherrors[it->first]->GetXaxis()->SetTitle("mass bin");
                map_sys_grapherrors[it->first]->Draw("AP");
            }
            first_draw = false;
        }
        else
        {
            if(isHistStye == false)
                map_sys_graph[it->first]->Draw("PC SAME");
            else
                map_sys_grapherrors[it->first]->Draw("P SAME");
        }
        if(isHistStye == false)
            leg->AddEntry(map_sys_graph[it->first], it->first, "pl");
        else
        {
            if(it->first == "Scale")
                leg->AddEntry(map_sys_grapherrors[it->first], "#mu_{F}/#mu_{R} variations", "pl");
            else if (it->first == "iterEM")
                leg->AddEntry(map_sys_grapherrors[it->first], "Unfolding method", "pl");
            else if (it->first == "AlphaS")
                leg->AddEntry(map_sys_grapherrors[it->first], "#alpha_{S}", "pl");
            else
                leg->AddEntry(map_sys_grapherrors[it->first], it->first, "pl");
        }
        it++;
    }
    if(isHistStye == false)
    {
        // Stat
        if(var.Contains("Pt"))
        {
            map_sys_graph["Stat."] = new TGraph(nMassBin, &v_bins[0], &meanPtStatRelErr_data_unfolded[0]);
        }
        else
        {
            map_sys_graph["Stat."] = new TGraph(nMassBin, &v_bins[0], &meanMassStatRelErr_data_unfolded[0]);
        }
        map_sys_graph["Stat."]->SetMarkerSize(0);
        map_sys_graph["Stat."]->SetLineColor(1);
        map_sys_graph["Stat."]->SetLineStyle(1);
        map_sys_graph["Stat."]->SetLineWidth(3);
        map_sys_graph["Stat."]->Draw("PC SAME");
        leg->AddEntry(map_sys_graph["Stat."], "Statistical", "l");
    }
    else
    {
        if(var.Contains("Pt"))
        {
            map_sys_grapherrors["Stat."] = new TGraphErrors(nMassBin, &v_bins[0], &meanPtStatRelErr_data_unfolded[0], &v_bins_error[0], &v_dummy_error[0]);
        }
        else
        {
            map_sys_grapherrors["Stat."] = new TGraphErrors(nMassBin, &v_bins[0], &meanMassStatRelErr_data_unfolded[0], &v_bins_error[0], &v_dummy_error[0]);
        }

        map_sys_grapherrors["Stat."]->SetMarkerSize(0);
        map_sys_grapherrors["Stat."]->SetLineColor(1);
        map_sys_grapherrors["Stat."]->SetLineStyle(1);
        map_sys_grapherrors["Stat."]->SetLineWidth(10);
        map_sys_grapherrors["Stat."]->Draw("P SAME");
        leg->AddEntry(map_sys_grapherrors["Stat."], "Statistical", "l");
    }

    if(isHistStye == false)
    {
        // Stat
        if(var.Contains("Pt"))
        {
            map_sys_graph["Accept."] = new TGraph(nMassBin, &v_bins[0], &meanPtRelSysErr_data_acc_corrected[0]);
        }
        else
        {
            map_sys_graph["Accept."] = new TGraph(nMassBin, &v_bins[0], &meanMassRelSysErr_data_acc_corrected[0]);
        }
        map_sys_graph["Accept."]->SetMarkerSize(0);
        map_sys_graph["Accept."]->SetLineColor(2);
        map_sys_graph["Accept."]->SetLineStyle(5);
        map_sys_graph["Accept."]->SetLineWidth(3);
        map_sys_graph["Accept."]->Draw("PC SAME");
        leg->AddEntry(map_sys_graph["Accept."], "Acceptance correction", "l");
    }
    else
    {
        if(var.Contains("Pt"))
        {
            map_sys_grapherrors["Accept."] = new TGraphErrors(nMassBin, &v_bins[0], &meanPtRelSysErr_data_acc_corrected[0], &v_bins_error[0], &v_dummy_error[0]);
        }
        else
        {
            map_sys_grapherrors["Accept."] = new TGraphErrors(nMassBin, &v_bins[0], &meanMassRelSysErr_data_acc_corrected[0], &v_bins_error[0], &v_dummy_error[0]);
        }

        map_sys_grapherrors["Accept."]->SetMarkerSize(0);
        map_sys_grapherrors["Accept."]->SetLineColor(2);
        map_sys_grapherrors["Accept."]->SetLineStyle(5);
        map_sys_grapherrors["Accept."]->SetLineWidth(10);
        map_sys_grapherrors["Accept."]->Draw("P SAME");
        leg->AddEntry(map_sys_grapherrors["Accept."], "Acceptance correction", "l");
    }

    int iPeriod_ = 4;
    if(year == 2017)
        iPeriod_ = 5;
    if(year == 2018)
        iPeriod_ = 6;
    CMS_lumi(c_out, iPeriod_, 11);

    leg->Draw();
    c_out->cd();
    c_out->SaveAs(output_baseDir+"Systematic_"+var+".pdf");

    delete c_out;
}

void ISRUnfold::drawSystematics_Acceptance(TString var, bool isHistStye)
{

    // Loop over
    // meanMass_data_folded_systematic, meanPt_data_folded_systematic
    setTDRStyle();
    writeExtraText = true;
    extraText  = "Work in progress";
    gStyle->SetLineWidth(globalLinedWidth);
    gStyle->SetFrameLineWidth(globalFrameWidth);
    gROOT->ForceStyle();

    // Prepare TGraph for each systematic
    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;
    vector<double> v_dummy_error;
    vector<double> v_bins_error;
    int nSystematics = sysMap.size() - 4;
    for(int i = 0; i < nMassBin; i++)
    {
        v_bins_error.push_back(1./(nSystematics+2)*0.5);
        v_dummy_error.push_back(0.0);
    }

    vector<double> v_bins_guideLine;
    vector<double> v_bins_error_guideLine;
    for(int i = 0; i < nMassBin; i++)
    {
        v_bins_guideLine.push_back(i+1);
        v_bins_error_guideLine.push_back(0.5);
    }

    map<TString, TGraph*> map_sys_graph;
    map<TString, TGraphErrors*> map_sys_grapherrors;
    map<TString, TGraphErrors*> map_sys_grapherrors_guideLine;
    std::map<TString, std::vector<TString>>::iterator it = sysMap.begin();


    int ithSys = 0;
    while(it != sysMap.end())
    {
        if((it->first).Contains("Efficiency_") || (it->first).Contains("Acceptance_"))
        {
            it++;
            continue;
        }

        vector<double> v_bins;
        for(int i = 0; i < nMassBin; i++)
        {
            v_bins.push_back(0.5 + i + 1./(nSystematics+2)*0.5 + 1./(nSystematics+2)*(ithSys));
        }

        if(var.Contains("Pt"))
        {
            if(isHistStye == false)
                map_sys_graph[it->first] = new TGraph(nMassBin, &v_bins[0], &meanPt_data_accept_rel_systematic[it->first][0]);
            else
            {
                map_sys_grapherrors[it->first] = new TGraphErrors(nMassBin, &v_bins[0], &meanPt_data_accept_rel_systematic[it->first][0], &v_bins_error[0], &v_dummy_error[0]);
                map_sys_grapherrors_guideLine[it->first] = new TGraphErrors(nMassBin, &v_bins_guideLine[0], &meanPt_data_accept_rel_systematic[it->first][0], &v_bins_error_guideLine[0], &v_dummy_error[0]);
            }
        }
        else
        {
            if(isHistStye == false)
                map_sys_graph[it->first] = new TGraph(nMassBin, &v_bins[0], &meanMass_data_accept_rel_systematic[it->first][0]);
            else
            {
                map_sys_grapherrors[it->first] = new TGraphErrors(nMassBin, &v_bins[0], &meanMass_data_accept_rel_systematic[it->first][0], &v_bins_error[0], &v_dummy_error[0]);
                map_sys_grapherrors_guideLine[it->first] = new TGraphErrors(nMassBin, &v_bins_guideLine[0], &meanMass_data_accept_rel_systematic[it->first][0], &v_bins_error_guideLine[0], &v_dummy_error[0]);
            }
        }
        it++;
        ithSys++;
        v_bins.clear();
    }

    // Create canvas
    TCanvas* c_out = new TCanvas("relative_uncertainty_"+var, "relative_uncertainty_"+var, 3000, 2500);
    c_out->SetGridy(1);
    //c_out->SetGridx(1);
    c_out->SetLogy();
    c_out->SetTopMargin(0.08);
    c_out->Draw();
    c_out->cd();

    TLegend* leg = new TLegend(0.3, 0.62, 0.95, 0.92,"","brNDC");
    leg->SetNColumns(2);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(43);
    leg->SetTextSize(50);

    int marker = 20;
    int markerSize = 5;
    int markerColor = 1;
    bool first_draw = true;

    it = sysMap.begin();
    while(it != sysMap.end())
    {
        if((it->first).Contains("Efficiency_") || (it->first).Contains("Acceptance_"))
        {
            it++;
            continue;
        }
        if(markerColor == 10) markerColor = 1;
        if(marker == 24) marker = 29;
        if(marker == 30) marker = 33;
        if(marker == 35) marker = 39;
        if(marker == 40) marker = 41;
        if(marker == 42) marker = 47;
        if(marker > 47) marker = 20;

        if(isHistStye == false)
        {
            map_sys_graph[it->first]->SetMarkerStyle(marker);
            map_sys_graph[it->first]->SetMarkerSize(markerSize);
            map_sys_graph[it->first]->SetMarkerColor(markerColor==5?46:markerColor);
            map_sys_graph[it->first]->SetLineColor(markerColor==5?46:markerColor);
        }
        else
        {
            map_sys_grapherrors[it->first]->SetLineWidth(2);
            map_sys_grapherrors[it->first]->SetMarkerStyle(marker);
            map_sys_grapherrors[it->first]->SetMarkerSize(markerSize);
            map_sys_grapherrors[it->first]->SetMarkerColor(markerColor==5?46:markerColor);
            map_sys_grapherrors[it->first]->SetLineColor(markerColor==5?46:markerColor);

            map_sys_grapherrors_guideLine[it->first]->SetLineWidth(2);
            map_sys_grapherrors_guideLine[it->first]->SetMarkerStyle(marker);
            map_sys_grapherrors_guideLine[it->first]->SetMarkerSize(0);
            map_sys_grapherrors_guideLine[it->first]->SetMarkerColor(markerColor==5?46:markerColor);
            map_sys_grapherrors_guideLine[it->first]->SetLineColor(markerColor==5?46:markerColor);
        }
        markerColor++;
        marker++;

        if(first_draw)
        {
            //map_sys_graph[it->first]->SetTitleOffset(0.02);
            if(isHistStye == false)
            {
                if(var.Contains("Pt"))
                    map_sys_graph[it->first]->SetTitle("Relative uncertainty on <p_{T}^{\ell\ell}> [%]");
                else
                    map_sys_graph[it->first]->SetTitle("Relative uncertainty on <Mass^{\ell\ell}> [%]");

                map_sys_graph[it->first]->GetYaxis()->SetRangeUser(0., 5.);
                map_sys_graph[it->first]->GetXaxis()->SetLimits(0.5,5.5);
                map_sys_graph[it->first]->GetYaxis()->SetTitleOffset(1.0);
                map_sys_graph[it->first]->GetYaxis()->SetTitle("Relative uncertainty [%]");
                map_sys_graph[it->first]->GetXaxis()->SetTitleOffset(0.7);
                map_sys_graph[it->first]->GetXaxis()->SetTitle("mass bin");
                map_sys_graph[it->first]->Draw("APC");
            }
            else
            {
                if(var.Contains("Pt"))
                    map_sys_grapherrors[it->first]->SetTitle("Relative uncertainty on <p_{T}^{\ell\ell}> [%]");
                else
                    map_sys_grapherrors[it->first]->SetTitle("Relative uncertainty on <Mass^{\ell\ell}> [%]");

                map_sys_grapherrors[it->first]->GetYaxis()->SetRangeUser(1e-3, 90.);
                map_sys_grapherrors[it->first]->GetXaxis()->SetLimits(0.5,5.5);
                map_sys_grapherrors[it->first]->GetYaxis()->SetTitleOffset(1.0);
                map_sys_grapherrors[it->first]->GetYaxis()->SetTitle("Relative uncertainty [%]");
                map_sys_grapherrors[it->first]->GetXaxis()->SetTitleOffset(0.7);
                map_sys_grapherrors[it->first]->GetXaxis()->SetTitle("mass bin");
                map_sys_grapherrors[it->first]->GetXaxis()->SetNdivisions(5);
                map_sys_grapherrors[it->first]->Draw("AP");
                map_sys_grapherrors_guideLine[it->first]->Draw("P");

                TLine massEdgeLine;
                massEdgeLine.SetLineStyle(1);
                massEdgeLine.DrawLine(1.5, 1e-3, 1.5, 90.);
                massEdgeLine.DrawLine(2.5, 1e-3, 2.5, 90.);
                massEdgeLine.DrawLine(3.5, 1e-3, 3.5, 90.);
                massEdgeLine.DrawLine(4.5, 1e-3, 4.5, 90.);
            }
            first_draw = false;
        }
        else
        {
            if(isHistStye == false)
                map_sys_graph[it->first]->Draw("PC SAME");
            else
            {
                map_sys_grapherrors[it->first]->Draw("P SAME");
                map_sys_grapherrors_guideLine[it->first]->Draw("P SAME");
            }
        }

        TString systematicName = getSysNameToShow(it->first);
        if(isHistStye == false)
        {
            leg->AddEntry(map_sys_graph[it->first], systematicName, "p");
        }
        else
        {
            leg->AddEntry(map_sys_grapherrors[it->first], systematicName, "p");
        }
        it++;
    }

    vector<double> v_stat;
    vector<double> v_stat_err;
    for(int i = 0; i < nMassBin; i++)
    {
        v_stat.push_back(1 + i);
        v_stat_err.push_back(0.5);
    }

    vector<double> v_eff;
    for(int i = 0; i < nMassBin; i++)
    {
        v_eff.push_back(0.5 + i + 1./(nSystematics+2)*0.5 + 1./(nSystematics+2)*(ithSys));
    }
    ithSys++;

    vector<double> v_accept;
    for(int i = 0; i < nMassBin; i++)
    {
        v_accept.push_back(0.5 + i + 1./(nSystematics+2)*0.5 + 1./(nSystematics+2)*(ithSys));
    }
    ithSys++;
    // Stat
    if(isHistStye == false)
    {
        if(var.Contains("Pt"))
        {
            map_sys_graph["Stat."] = new TGraph(nMassBin, &v_stat[0], &meanPtStatRelErr_data_unfolded[0]);
        }
        else
        {
            map_sys_graph["Stat."] = new TGraph(nMassBin, &v_stat[0], &meanMassStatRelErr_data_unfolded[0]);
        }

        map_sys_graph["Stat."]->SetMarkerSize(0);
        map_sys_graph["Stat."]->SetLineColor(1);
        map_sys_graph["Stat."]->SetLineStyle(1);
        map_sys_graph["Stat."]->SetLineWidth(10);
        map_sys_graph["Stat."]->Draw("PC SAME");
        leg->AddEntry(map_sys_graph["Stat."], "Statistical", "l");
    }
    else
    {

        if(var.Contains("Pt"))
        {
            map_sys_grapherrors["Stat."] = new TGraphErrors(nMassBin, &v_stat[0], &meanPtStatRelErr_data_unfolded[0], &v_stat_err[0], &v_dummy_error[0]);
            map_sys_grapherrors["Efficiency"] = new TGraphErrors(nMassBin, &v_eff[0], &meanPtEffRelErr_data[0], &v_bins_error[0], &v_dummy_error[0]);
            map_sys_grapherrors["Acceptance"] = new TGraphErrors(nMassBin, &v_accept[0], &meanPtAcceptRelErr_data[0], &v_bins_error[0], &v_dummy_error[0]);
        }
        else
        {
            map_sys_grapherrors["Stat."] = new TGraphErrors(nMassBin, &v_stat[0], &meanMassStatRelErr_data_unfolded[0], &v_stat_err[0], &v_dummy_error[0]);
            map_sys_grapherrors["Efficiency"] = new TGraphErrors(nMassBin, &v_eff[0], &meanMassEffRelErr_data[0], &v_bins_error[0], &v_dummy_error[0]);
            map_sys_grapherrors["Acceptance"] = new TGraphErrors(nMassBin, &v_accept[0], &meanMassAcceptRelErr_data[0], &v_bins_error[0], &v_dummy_error[0]);
        }

        map_sys_grapherrors["Stat."]->SetMarkerSize(0);
        map_sys_grapherrors["Stat."]->SetLineColor(1);
        map_sys_grapherrors["Stat."]->SetLineStyle(1);
        map_sys_grapherrors["Stat."]->SetLineWidth(10);
        map_sys_grapherrors["Stat."]->Draw("P SAME");
        leg->AddEntry(map_sys_grapherrors["Stat."], "Statistical", "l");

        map_sys_grapherrors["Efficiency"]->SetMarkerSize(4);
        map_sys_grapherrors["Efficiency"]->SetMarkerStyle(49);
        map_sys_grapherrors["Efficiency"]->SetMarkerColor(8);
        map_sys_grapherrors["Efficiency"]->SetLineColor(8);
        map_sys_grapherrors["Efficiency"]->SetLineStyle(1);
        map_sys_grapherrors["Efficiency"]->Draw("P SAME");
        leg->AddEntry(map_sys_grapherrors["Efficiency"], "Efficiency", "p");

        map_sys_grapherrors["Acceptance"]->SetMarkerSize(4);
        map_sys_grapherrors["Acceptance"]->SetMarkerStyle(49);
        map_sys_grapherrors["Acceptance"]->SetMarkerColor(4);
        map_sys_grapherrors["Acceptance"]->SetLineColor(4);
        map_sys_grapherrors["Acceptance"]->SetLineStyle(1);
        map_sys_grapherrors["Acceptance"]->Draw("P SAME");
        leg->AddEntry(map_sys_grapherrors["Acceptance"], "Acceptance", "p");
    }

    v_stat.clear();
    v_accept.clear();

    double x_ = 0.18;
    double y_ = 0.85;
    TString lepton_type;
    if(channel_name=="electron"){
        lepton_type = "ee channel";
    }
    else
    {
        lepton_type = "#mu#mu channel";
    }

    TLatex printChannel;
    printChannel.SetTextFont(43);
    printChannel.SetTextSize(40 * 2);

    printChannel.DrawLatexNDC(x_, y_, lepton_type);

    int iPeriod_ = 4;
    if(year == 2017)
        iPeriod_ = 5;
    if(year == 2018)
        iPeriod_ = 6;
    CMS_lumi(c_out, iPeriod_, 0);

    leg->Draw();
    c_out->cd();
    c_out->SaveAs(output_baseDir+"SystematicAccept_"+var+".pdf");

    delete c_out;
}

TH1* ISRUnfold::getUnfInput(TString var, TString steering, bool useAxis, int massBin, bool binWidth)
{
    TH1* outHist = NULL;
    if(var.Contains("Mass"))
    {
        outHist = nomMassUnfold->GetInput("unfold_input_mass", 0, 0, steering, useAxis);
        if(binWidth)
        {
            divideByBinWidth(outHist, false);
        }
    }
    else
    {
        TString nth;
        nth.Form("%d", massBin);
        outHist = nomPtUnfold->GetInput("unfold_input_pt_"+nth, 0, 0, steering, useAxis);
        if(binWidth)
        {
            divideByBinWidth(outHist, false);
        }
    }
    return outHist;
}
