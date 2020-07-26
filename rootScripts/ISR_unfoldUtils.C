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
    TH2D* hProb;
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

    double chi2 = 0.;
    double ndf  = 0.;

    TH1* hData; // Data - Bkg
    TH1* hDY; // DY MC
    TH2* hRho;

    TH1 *g_fcnHist=0;
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

void ISRUnfold::drawResponseM(TString var, TString sysName, TString sysPostfix, bool isDetector)
{
    const TVectorD* xaxis1_tvecd;
    int xaxis1_nbin;

    const TVectorD* yaxis1_tvecd = NULL;
    int yaxis1_nbin;

    int nMassBin = massBinEdges.size() - 1;

    setTDRStyle();
    writeExtraText = true;
    extraText  = "simulation";

    TCanvas* c1 = new TCanvas("c1","c1", 50, 50, 2400, 3000);
    gStyle->SetOptFit(0);
    gStyle->SetPalette(55);
    c1->cd();

    c1->SetBottomMargin(0.2);
    c1->SetRightMargin(0.15);
    c1->SetLeftMargin(0.2);
    c1->SetTopMargin(0.08);
    c1->SetTicks(1);
    c1->SetLogz();

    TH2 *histProb = NULL;
    TH2 *histProb_woUO = NULL;
    bool draw_wo_UO = true;

    if(var=="Pt")
    {
        xaxis1_tvecd = pt_binning_Gen->GetDistributionBinning(0);
        xaxis1_nbin = xaxis1_tvecd->GetNrows() - 1; // number of bins without UO

        yaxis1_tvecd = pt_binning_Rec->GetDistributionBinning(0);
        yaxis1_nbin = yaxis1_tvecd->GetNrows() - 1;

        histProb = nomPtUnfold->GetProbabilityMatrix("Migration prob. for pt mass bin",";p_T(gen);p_T(Rec)");
        if(draw_wo_UO)
        {
            histProb_woUO = new TH2D("responseM_woUO","responseM_woUO", xaxis1_nbin * nMassBin, 0, xaxis1_nbin * nMassBin, yaxis1_nbin * nMassBin, 0, yaxis1_nbin * nMassBin);

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

    histProb_woUO->Draw("COLZ");
    histProb_woUO->GetZaxis()->SetRangeUser(1e-3, 1.0);
    histProb_woUO->GetYaxis()->SetTitleFont(43);
    histProb_woUO->GetYaxis()->SetTitleSize(100);
    histProb_woUO->GetYaxis()->SetTitleOffset(1.5);
    histProb_woUO->GetXaxis()->SetTitleFont(43);
    histProb_woUO->GetXaxis()->SetTitleSize(100);
    histProb_woUO->GetXaxis()->SetTitleOffset(1.5);

    vector<TGaxis*> v_xaxis; 

    if(var.Contains("Pt"))
    {
        histProb_woUO->GetXaxis()->SetLabelSize(0);
        histProb_woUO->GetXaxis()->SetTickSize(0);
        for(int ibin = 0; ibin < nMassBin; ibin++)
        {
            v_xaxis.push_back( new TGaxis(9 * ibin, 0, 9 * (ibin+1), 0, 0, 100, 10,""));
            v_xaxis.at(ibin)->SetLabelFont(43);
            v_xaxis.at(ibin)->SetLabelSize(0);
            v_xaxis.at(ibin)->Draw();
            v_xaxis.at(ibin)->ChangeLabel(10,-1,100,-1,-1,-1,"100");  
            v_xaxis.at(ibin)->SetLabelOffset(0.05);  
            v_xaxis.at(ibin)->Draw();
        }
    }
    //TGaxis *axis1 = new TGaxis(0, 0, 9, 0, 0, 100, 10,"");
    //axis1->SetLabelFont(43);
    //axis1->SetLabelSize(0);
    //axis1->Draw();
    //axis1->ChangeLabel(10,-1,100,-1,-1,-1,"100");  
    //axis1->SetLabelOffset(0.05);  
    //axis1->Draw();

    //TGaxis *axis2 = new TGaxis(9, 0, 18, 0, 0, 100, 10,"");
    //axis2->SetLabelFont(43);
    //axis2->SetLabelSize(0);
    //axis2->Draw();
    //axis2->ChangeLabel(10,-1,100,-1,-1,-1,"100");  
    //axis2->SetLabelOffset(0.05);  
    //axis2->Draw();

    //if(var == "Pt")
    //{

    //    ticks_ = new TH2D("tick", "tick", histProb_woUO->GetNbinsX(), histProb_woUO->GetXaxis()->GetXmin(), histProb_woUO->GetXaxis()->GetXmax(),
    //                                      histProb_woUO->GetNbinsY(), histProb_woUO->GetYaxis()->GetXmin(), histProb_woUO->GetYaxis()->GetXmax());


    //    int center = (xaxis1_nbin) / 2;
    //    int totalBins = xaxis1_nbin;
    //    int center_y = (yaxis1_nbin) / 2;
    //    int totalBins_y = yaxis1_nbin;

    //    //if(pt_binning_Gen->HasUnderflow(1))
    //    //{
    //    //    ticks_->GetXaxis()->SetBinLabel(histProb_woUO->GetXaxis()->FindBin(center), "Underflow");  // 5 = # of pt bins / 2
    //    //    ticks_->GetYaxis()->SetBinLabel(histProb_woUO->GetYaxis()->FindBin(center_y), "Underflow");  // 5 = # of pt bins / 2

    //    //    center += totalBins;
    //    //    center_y += totalBins_y;
    //    //}

    //    for(int ibin = 0; ibin < nMassBin; ibin++)
    //    {
    //        TString lowMassEdge;
    //        TString highMassEdge;

    //        lowMassEdge.Form("%d", (int)massBinEdges[ibin]);
    //        highMassEdge.Form("%d", (int)massBinEdges[ibin+1]);
    //        TString var_name = "M";
    //        if(var == "Mass") var_name = "p_{T}";
    //            ticks_->GetXaxis()->SetBinLabel(histProb_woUO->GetXaxis()->FindBin(center), lowMassEdge+"<"+var_name+"<"+highMassEdge+" GeV");  // TODO option to set mass region

    //        ticks_->GetYaxis()->SetBinLabel(histProb_woUO->GetYaxis()->FindBin(center_y), lowMassEdge+"<"+var_name+"<"+highMassEdge+" GeV");

    //        center += totalBins;
    //        center_y += totalBins_y;
    //    }
    //    if(pt_binning_Gen->HasOverflow(1))
    //    {
    //        ticks_->GetXaxis()->SetBinLabel(histProb_woUO->GetXaxis()->FindBin(center), "Overflow");
    //        ticks_->GetYaxis()->SetBinLabel(histProb_woUO->GetYaxis()->FindBin(center_y), "Overflow");
    //    }

    //    //ticks_->GetYaxis()->SetTitleFont(43);
    //    //ticks_->GetYaxis()->SetTitleSize(40);
    //    //ticks_->GetYaxis()->SetTitleOffset(1.2);
    //    //ticks_->GetXaxis()->SetTitleFont(43);
    //    //ticks_->GetXaxis()->SetTitleSize(40);
    //    //ticks_->GetXaxis()->SetTitleOffset(1.2);

    //    ticks_->GetXaxis()->SetLabelColor(kGray+2);
    //    ticks_->GetYaxis()->SetLabelColor(kGray+2);
    //    ticks_->GetXaxis()->SetLabelSize(30);
    //    ticks_->GetYaxis()->SetLabelSize(30);
    //    //ticks_->GetXaxis()->LabelsOption("v");
    //    ticks_->GetXaxis()->SetTickSize(0);
    //    ticks_->GetYaxis()->SetTickSize(0);
    //    ticks_->Draw("SAME COLZ"); //
    //    //histProb_woUO->Draw("same axis COLZ"); // "colz same" not working
    //    //draw_option = "same COLZ";
    //}

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
    grid_.SetLineStyle(2);

    if(var=="Pt")
    {
    int boundarybin_x = 1;
    for( int ii=0; ii<histProb_woUO->GetXaxis()->GetNbins(); ii++ )
    {
        Int_t i_bin = ii+1;
        Double_t binEdge = histProb_woUO->GetXaxis()->GetBinLowEdge(i_bin);

        if(boundarybin_x == i_bin)
        {
            grid_.DrawLine(binEdge, histProb_woUO->GetYaxis()->GetBinUpEdge(0), binEdge, histProb_woUO->GetYaxis()->GetBinUpEdge(histProb_woUO->GetYaxis()->GetNbins()) );
            boundarybin_x += xaxis1_nbin; // next edge to draw
        }
    }

    int boundarybin_y = 1;
    for( int ii=0; ii<histProb_woUO->GetYaxis()->GetNbins(); ii++ )
    {
        Int_t i_bin = ii+1;
        Double_t binEdge = histProb_woUO->GetYaxis()->GetBinLowEdge(i_bin);

        if(boundarybin_y == i_bin)
        {
            grid_.DrawLine(histProb_woUO->GetXaxis()->GetBinUpEdge(0), binEdge, histProb_woUO->GetXaxis()->GetBinUpEdge(histProb_woUO->GetXaxis()->GetNbins()), binEdge);
            boundarybin_y += yaxis1_nbin; // next edge to draw
        }
    }
    }

    c1->RedrawAxis();
    CMS_lumi(c1, 7, 11);
    c1->cd();

    c1->SaveAs(isDetector?var + "_responseM_" + year + ".png":var + "_FSRresponseM_" + year + ".png");
    delete c1;
}

void ISRUnfold::setNomResMatrix(TString var, TString filepath, TString dirName, TString histName, TString binDef)
{
    //cout << "ISRUnfold::setNomResMatrix set response matrix..." << endl;
    TFile* filein = new TFile(filepath);

    TString Rec_binName = "Rec_"+var;
    TString Gen_binName = "Gen_"+var;
    Rec_binName = dirName + "/" + var + "_ResMatrix_" + histName + binDef + "/" + Rec_binName;
    Gen_binName = dirName + "/" + var + "_ResMatrix_" + histName + binDef + "/" + Gen_binName;

    //cout << "Rec_binName: " << Rec_binName << endl;
    //cout << "Gen_binName: " << Gen_binName << endl;
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
        cout << "ISRUnfold::setNomResMatrix, only Pt and Mass available for var" << endl;
        exit (EXIT_FAILURE);
    }

    // Set mass bin edges
    // setBassBinEdges
    setMassBindEdges();

    // Set response matrix
    TH2* hmcGenRec;
    hmcGenRec = (TH2*)filein->Get(dirName + "/" + var + "_ResMatrix_" + histName + binDef + "/hmc" + var + "GenRec");
    //cout << dirName + "/" + var + "_ResMatrix_" + histName + binDef + "/hmc" + var + "GenRec" << endl;

    if( var == "Pt" )
    {
    	nomPtUnfold = new TUnfoldDensityV17(hmcGenRec,
    	                                    TUnfold::kHistMapOutputHoriz,
    	                                    regMode,
    	                                    TUnfold::kEConstraintArea,
    	                                    TUnfoldDensityV17::kDensityModeBinWidth,
    	                                    pt_binning_Gen,pt_binning_Rec);

        hPtResponseM = (TH2*) hmcGenRec->Clone("hPtResponseM");

        // For statistical uncertainty
        if(makeStatUnfold)
        {
            cout << "Create response matrix for statistical uncertainty..." << endl;
            for(int i = 0; i < statSize; i++)
            {
                statPtUnfold.push_back(new TUnfoldDensityV17(hmcGenRec,
    	                                                    TUnfold::kHistMapOutputHoriz,
    	                                                    regMode,
    	                                                    TUnfold::kEConstraintArea,
    	                                                    TUnfoldDensityV17::kDensityModeBinWidth,
    	                                                    pt_binning_Gen,pt_binning_Rec));
            }
        }
    }
    else
    {
        nomMassUnfold = new TUnfoldDensityV17(hmcGenRec,
                                              TUnfold::kHistMapOutputHoriz,
                                              regMode,
                                              TUnfold::kEConstraintArea,
                                              TUnfoldDensityV17::kDensityModeBinWidth,
                                              mass_binning_Gen,mass_binning_Rec);

        hMassResponseM = (TH2*) hmcGenRec->Clone("hMassResponseM");

        // For statistical uncertainty
        if(makeStatUnfold)
        {
            for(int i = 0; i < statSize; i++)
            {
                statMassUnfold.push_back(new TUnfoldDensityV17(hmcGenRec,
                                                               TUnfold::kHistMapOutputHoriz,
                                                               regMode,
                                                               TUnfold::kEConstraintArea,
                                                               TUnfoldDensityV17::kDensityModeBinWidth,
                                                               mass_binning_Gen,mass_binning_Rec));
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
            // Systematic not considered in this ISRUnfold class, create TUnfoldDensity using the DEFAULT response matrix
            int size = sysMap_previous[it->first].size();
            for(int ith = 0; ith < size; ith++)
            {
                //cout << "Systematic variation, " << sysMap_previous[it->first][ith] << endl;
                this->sysPtUnfold[it->first][sysMap_previous[it->first][ith]] = new TUnfoldDensityV17(hPtResponseM,TUnfold::kHistMapOutputHoriz,
                                                                                                regMode, TUnfold::kEConstraintArea, TUnfoldDensityV17::kDensityModeBinWidth,
                                                                                                pt_binning_Gen,pt_binning_Rec);

                this->sysMassUnfold[it->first][sysMap_previous[it->first][ith]] = new TUnfoldDensityV17(hMassResponseM,TUnfold::kHistMapOutputHoriz,
                                                                                                regMode, TUnfold::kEConstraintArea, TUnfoldDensityV17::kDensityModeBinWidth,
                                                                                                mass_binning_Gen, mass_binning_Rec);

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


                this->sysMap[it->first].push_back(sysMap_previous[it->first][ith]);
            }
        }
        else
        {
            // Found
            // Systematic only for this ISRUnfold class
            // Loop over systematic varations

            // If previous ISRUnfold class have this systematic result, then use the output as input for this ISRUnfold class
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
                    this->sysPtUnfold[it->first][sysMap_previous[it->first][ith]]->SetInput(unfold->hSysFullPhasePtData[it->first][sysMap_previous[it->first][ith]], nominal_bias);
                    this->sysMassUnfold[it->first][sysMap_previous[it->first][ith]]->SetInput(unfold->hSysFullPhaseMassData[it->first][sysMap_previous[it->first][ith]], nominal_bias);
                }
            }
            // Else, use the default input
        }
        it++;
    }
}

void ISRUnfold::setSysTUnfoldDensity(TString var, TString filepath, TString dirName, TString histName, TString sysName, TString sysPostfix, TString binDef)
{
    TFile* filein = new TFile(filepath);
    TH2* hmcGenRec = NULL;

    if(sysName == "Fake")
    {
        hmcGenRec = (TH2*)filein->Get(dirName + "/" + var + "_ResMatrix_" + histName + binDef +"/hmc" + var + "GenRec");
    }
    else if(sysName.Contains("LepMom") || sysName.Contains("Unfold") || sysName.Contains("FSR") || sysName.Contains("ZptCorr"))
    {
        //cout << dirName + "/" + var + "_ResMatrix_" + histName + binDef +"/hmc" + var + "GenRec" << endl;
        hmcGenRec = (TH2*)filein->Get(dirName + "/" + var + "_ResMatrix_" + histName + binDef +"/hmc" + var + "GenRec");
    }
    else
    {
        hmcGenRec = (TH2*)filein->Get(dirName + "/" + var + "_ResMatrix_" + histName + binDef +"/hmc" + var + "GenRec_" + sysPostfix);
    }

    if( var == "Pt" )
    {
        //cout << "sys: " << sysName << " postfix: " << sysPostfix << endl;
        sysPtUnfold[sysName][sysPostfix] = new TUnfoldDensityV17(hmcGenRec,
                                                                 TUnfold::kHistMapOutputHoriz,
                                                                 regMode,
                                                                 TUnfold::kEConstraintArea,
                                                                 TUnfoldDensityV17::kDensityModeBinWidth,
                                                                 pt_binning_Gen,pt_binning_Rec);
    }
    else if( var == "Mass" )
    {
        sysMassUnfold[sysName][sysPostfix] = new TUnfoldDensityV17(hmcGenRec,
                                                                   TUnfold::kHistMapOutputHoriz,
                                                                   regMode,
                                                                   TUnfold::kEConstraintArea,
                                                                   TUnfoldDensityV17::kDensityModeBinWidth,
                                                                   mass_binning_Gen,mass_binning_Rec);
    }
    else
    {
        cout << "ISRUnfold::setSysTUnfoldDensity, only Pt and Mass available for var" << endl;
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
                sysPtUnfold[sysName][sysPostfix]->SetInput(unfold->hFullPhasePtData, 1.);
            }
            else
            {
                sysMassUnfold[sysName][sysPostfix]->SetInput(unfold->hFullPhaseMassData, 1.);
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

    // Nominal histograms
    // First input covariance test (ID)
    TFile* fcov = new TFile("/home/jhkim/ISR_Run2/unfolding/TUnfoldISR2016/rootScripts/covariance.root");
    TFile* fcov_pt = new TFile("/home/jhkim/ISR_Run2/unfolding/TUnfoldISR2016/rootScripts/covariance_pt.root");
    TH2* hCov = (TH2*) fcov->Get("cov");
    TH2* hCov_pt = (TH2*) fcov_pt->Get("cov");

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
            sysPtUnfold[sysName][sysPostfix]->SetInput(hRec, nominal_bias);
        }
        else if(var == "Mass")
        {
            //cout << "sysName: " << sysName << " postfix: " << sysPostfix << endl;
            sysMassUnfold[sysName][sysPostfix]->SetInput(hRec, nominal_bias);
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

void ISRUnfold::subBkgs(TString filepath, std::pair<TString, TString>& bkgInfo, bool isSys, TString binDef, TString dirName, TString sysName, TString sysPostfix, bool isFSR)
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
        if(!isFSR)
        {
            if(!(sysName.Contains("LepMom") || sysName.Contains("Unfold") || sysName.Contains("ZptCorr")))
            {
                TString histPostfix = bkgInfo.first + "_" + sysPostfix;
                if(bkgInfo.second != "Fake" && sysName == "Fake") histPostfix = bkgInfo.first;
                if(bkgInfo.second == "Fake" && sysName != "Fake") histPostfix = bkgInfo.first;

                //cout << dirName + "/Pt"+binDef+"/histo_" + histPostfix << endl;
                hPtRec = (TH1*)filein->Get(dirName + "/Pt"+binDef+"/histo_" + histPostfix);
                sysPtUnfold[sysName][sysPostfix]->SubtractBackground(hPtRec, bkgInfo.first);

                hMassRec = (TH1*)filein->Get(dirName + "/Mass"+binDef+"/histo_" + histPostfix);
                sysMassUnfold[sysName][sysPostfix]->SubtractBackground(hMassRec, bkgInfo.first);
            }
            else
            {
                TString histPostfix = bkgInfo.first;

                //cout << dirName + "/Pt"+binDef+"/histo_" + histPostfix << endl;
                hPtRec = (TH1*)filein->Get(dirName + "/Pt"+binDef+"/histo_" + histPostfix);
                sysPtUnfold[sysName][sysPostfix]->SubtractBackground(hPtRec, bkgInfo.first);

                hMassRec = (TH1*)filein->Get(dirName + "/Mass"+binDef+"/histo_" + histPostfix);
                sysMassUnfold[sysName][sysPostfix]->SubtractBackground(hMassRec, bkgInfo.first);
            }
        }
        else
        {
            if((sysName.Contains("Scale") || sysName.Contains("AlphaS")))
            {
                TString histPostfix = bkgInfo.first + "_" + sysPostfix;

                //cout << dirName + "/Pt"+binDef+"/histo_" + histPostfix << endl;
                hPtRec = (TH1*)filein->Get(dirName + "/Pt"+binDef+"/histo_" + histPostfix);
                sysPtUnfold[sysName][sysPostfix]->SubtractBackground(hPtRec, bkgInfo.first);

                hMassRec = (TH1*)filein->Get(dirName + "/Mass"+binDef+"/histo_" + histPostfix);
                sysMassUnfold[sysName][sysPostfix]->SubtractBackground(hMassRec, bkgInfo.first);
            }
            else
            {
                TString histPostfix = bkgInfo.first;

                //cout << dirName + "/Pt"+binDef+"/histo_" + histPostfix << endl;
                hPtRec = (TH1*)filein->Get(dirName + "/Pt"+binDef+"/histo_" + histPostfix);
                sysPtUnfold[sysName][sysPostfix]->SubtractBackground(hPtRec, bkgInfo.first);

                hMassRec = (TH1*)filein->Get(dirName + "/Mass"+binDef+"/histo_" + histPostfix);
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

// Draw detector distributions using input root file
TCanvas* ISRUnfold::drawFoldedHists(TString var, TString filePath, TString dirName, TString steering, bool useAxis, TString sysName, TString outName, int nthMassBin, bool divBinWidth, TString sysFilePath)
{
    // If steering == "", then usual TH1 histogram
    // If seering != "", TH1 from TUnfold

    double meanDipt = 0.;
    double meanDipt_bkgsub = 0.;
    double histMin = 5e-1;

    setTDRStyle();
    writeExtraText = true;
    extraText  = "Work in progress";
    gStyle->SetLineWidth(3);
    gStyle->SetFrameLineWidth(3);
    gROOT->ForceStyle();

    TH1::AddDirectory(kFALSE);
    //cout << "ISRUnfold::drawFoldedHists, Draw plot!" << endl;

    TFile* filein = new TFile(filePath);

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
    hDY = getRawHist(var, filePath, dirName, DYHistName_, "Signal", steering, useAxis, divBinWidth);

    hMCtotal = (TH1*) hDY->Clone("hMCtotal");
    hRatio = (TH1*) hData->Clone("hRatio");

    // Create canvas
    TCanvas* c_out = new TCanvas("detector_level_"+var, "detector_level_"+var, 3200, 2800);
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

    if(var.Contains("Mass"))
    {
        hData->GetXaxis()->SetRangeUser(massBinEdges.at(0), massBinEdges.at(5));
        hDY->GetXaxis()->SetRangeUser(massBinEdges.at(0), massBinEdges.at(5));
    }

    hData->SetTitle("");
    hData->SetStats(false);
    hData->GetXaxis()->SetMoreLogLabels(true);
    hData->Draw("p9e");
    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(3.2);
    hData->SetLineColor(kBlack);
    hData->GetYaxis()->SetTitle("Events/Bin");
    hData->SetMinimum(histMin);

    hDY->SetFillColor(kOrange);

    TLegend* leg = new TLegend(0.5, 0.6, 0.95, 0.85,"","brNDC");
    //leg->SetNColumns(2);
    leg->SetTextFont(43);
    leg->SetTextSize(100);
    leg->SetFillStyle(0); // transparent
    leg->SetBorderSize(0);
    leg->AddEntry(hData, "Data", "pl");
    leg->AddEntry(hDY, "Drell-Yan", "F");

    THStack* hsMC = new THStack("hsMC", "hsMC");
    setTHStack(var, filePath, dirName, *hsMC, *hMCtotal, *leg, steering, useAxis, "", divBinWidth);
    hsMC->Add(hDY);
    hData->SetMaximum(hsMC->GetMaximum() * 1e5);

    hsMC->Draw("hist same");
    hData->Draw("p9e same");
    pad1->RedrawAxis();

    // Get average transeverse momentum values
    meanDipt = hData->GetMean();
    TH1* hData_ = (TH1*) hData->Clone("DataBKGsubtracted");
    TH1* hMCtotal_ = (TH1*) hMCtotal->Clone("MCtotalDYsubtracted");
    hMCtotal_->Add(hDY, -1);
    hData_->Add(hMCtotal_, -1);
    meanDipt_bkgsub = hData_->GetMean();
    double binnedMean = getBinnedMean(hData); // To check how average value changes with binned histogram

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
    TLatex smearedChi2;
    smearedChi2.SetTextFont(43);
    smearedChi2.SetTextSize(100);
    TString chi2_;
    chi2_.Form("%.2f", getSmearedChi2(var, filePath, dirName, steering,  useAxis, false));
    smearedChi2.DrawLatexNDC(0.2, 0.6, "#chi^{2}: " + chi2_);

    leg->Draw();

    int iPeriod_ = 4;
    if(year == 2017)
        iPeriod_ = 5;
    if(year == 2018)
        iPeriod_ = 6;
    CMS_lumi(pad1, iPeriod_, 11);
    // writeCutInfo(pad, var, nthMassBin);
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
        for(int ibin = 2; ibin < hData->GetNbinsX()+1; ibin++)
        {
            massEdgeLine.SetLineStyle(1);
            massEdgeLine.SetLineColor(kWhite);
            massEdgeLine.DrawLine(hData->GetXaxis()->GetBinLowEdge(ibin), histMin, hData->GetXaxis()->GetBinLowEdge(ibin), last->GetBinContent(ibin));
        }
    }
    c_out->cd();

    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.35);
    if(var.Contains("Mass"))
        pad2->SetLogx();
    pad2->SetTicks(1);
    pad2->SetGridy(1);
    pad2->Draw();
    pad2->cd();

    if(var.Contains("Mass"))
    {
        hRatio->GetXaxis()->SetRangeUser(massBinEdges.at(0), massBinEdges.at(5));
        hMCtotal->GetXaxis()->SetRangeUser(massBinEdges.at(0), massBinEdges.at(5));
    }

    hRatio->SetStats(false);
    hRatio->Divide(hMCtotal);
    hRatio->Draw("p9e");
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(3.2);
    hRatio->SetLineColor(kBlack);
    hRatio->GetYaxis()->SetTitle("Data/MC");

    hRatio->SetMinimum(0.5);
    hRatio->SetMaximum(1.5);

    setXaxisTitle(hRatio, var, useAxis);

    // TODO Save systematic histograms
    TH1* sysBand_ratio_forData = NULL;
    TH1* sysBand_ratio_forMC = NULL;
    TLegend* leg_sys = NULL;
    if(sysName != "")
    {
        leg_sys = new TLegend(0.5, 0.75, 0.95, 0.85,"","brNDC");

        sysBand_ratio_forData = getDetectorSystematicBand(var, filePath, dirName, steering, useAxis,  sysName, hData, hDY, hMCtotal, hRatio, divBinWidth, true, false, sysFilePath, nthMassBin);
        sysBand_ratio_forData->SetFillColorAlpha(kBlack,0.8);

        sysBand_ratio_forMC = getDetectorSystematicBand(var, filePath, dirName, steering, useAxis,  sysName, hData, hDY, hMCtotal, hRatio, divBinWidth, true, true, sysFilePath, nthMassBin);
        sysBand_ratio_forMC->SetFillColorAlpha(kRed,0.2);


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

        sysBand_ratio_forData->SetFillStyle(3001);
        sysBand_ratio_forData->SetMarkerSize(0.);
        sysBand_ratio_forData->Draw("E2 same");
        sysBand_ratio_forMC->SetFillStyle(3001);
        sysBand_ratio_forMC->SetMarkerSize(0.);
        sysBand_ratio_forMC->Draw("E2 same");

        leg_sys->SetTextFont(43);
        leg_sys->SetTextSize(100);
        leg_sys->SetFillStyle(0); // transparent
        leg_sys->SetBorderSize(0);
        leg_sys->AddEntry(sysBand_ratio_forMC, sysName, "F");
        leg_sys->Draw();
            
    }
    hRatio->Draw("p9e same");

    TLine* l_ = new TLine(hRatio->GetXaxis()->GetXmin(), 1, hRatio->GetXaxis()->GetXmax(), 1);
    l_->SetLineColor(kRed);
    l_->Draw("same");
    l_->SetLineStyle(1);

    // Save canvas
    c_out->cd();
    c_out->SaveAs(outName!=""?output_baseDir+outName+var+sysName+".png":output_baseDir+"detector_"+var+sysName+".png");

    delete filein;
    return c_out;
}

TCanvas* ISRUnfold::drawUnfoldedHists(TString var, TString steering, bool useAxis, TString sysName, TString outName, int nthMassBin, bool divBinWidth)
{
    // If steering == "", then usual TH1 histogram
    // If seering != "", TH1 from TUnfold

    double histMin = 5e-1;

    setTDRStyle();
    writeExtraText = true;
    extraText  = "Work in progress";
    gStyle->SetLineWidth(3.);
    gStyle->SetFrameLineWidth(3.);
    gROOT->ForceStyle();

    TH1::AddDirectory(kFALSE);
    //cout << "ISRUnfold::drawFoldedHists, Draw plot!" << endl;

    // For nominal histogram
    TH1* hData = NULL;
    TH1* hDY = NULL;
    TH1* hRatio = NULL;

    if(var.Contains("Pt"))
    {
        hData = nomPtUnfold->GetOutput("hUnfoldedPt",0,0,steering,useAxis);
        hDY = nomPtUnfold->GetBias("hDYMCPt",0,0,steering,useAxis);;
    }
    else
    {
        hData = nomMassUnfold->GetOutput("hUnfoldedMass",0,0,steering,useAxis);
        hDY = nomMassUnfold->GetBias("hDYMCMass",0,0,steering,useAxis);;
    }
    if(divBinWidth)
    {
        divideByBinWidth(hData, false);
        divideByBinWidth(hDY, false);
    }
    hRatio = (TH1*) hData->Clone("hRatio");

    // Create canvas
    TCanvas* c_out ;
    c_out = new TCanvas("unfolded_level_"+var, "unfolded_level_"+var, 3200, 2800);

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

    hData->SetTitle("");
    hData->SetStats(false);
    hData->GetXaxis()->SetMoreLogLabels(true);
    hData->Draw("p9e");
    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(3.2);
    hData->SetLineColor(kBlack);
    hData->GetYaxis()->SetTitle("Events/Bin");
    hData->SetMaximum(hDY->GetMaximum() * 1e5);
    hData->SetMinimum(histMin);

    hDY->SetFillColor(kOrange);
    hDY->Draw("hist same");

    TLine massEdgeLine;
    for(int ibin = 2; ibin < hData->GetNbinsX()+1; ibin++)
    {
        massEdgeLine.SetLineStyle(1);
        massEdgeLine.SetLineColor(kWhite);
        massEdgeLine.DrawLine(hData->GetXaxis()->GetBinLowEdge(ibin), histMin, hData->GetXaxis()->GetBinLowEdge(ibin), hDY->GetBinContent(ibin));
    }

    TLegend* leg = new TLegend(0.5, 0.7, 0.9, 0.9,"","brNDC");
    //leg->SetNColumns(2);
    leg->SetTextFont(43);
    leg->SetTextSize(80);
    leg->SetFillStyle(0); // transparent
    leg->SetBorderSize(0);
    if(outName.Contains("Closure"))
        leg->AddEntry(hData, "Unfolded MC", "pl");
    else
        leg->AddEntry(hData, "Unfolded data", "pl");
    leg->AddEntry(hDY, "Drell-Yan (MG5_aMC@NLO)", "F");

    hData->Draw("p9e same");
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

    hRatio->SetStats(false);
    hRatio->Divide(hDY);
    hRatio->Draw("p9e");
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(3.2);
    hRatio->SetLineColor(kBlack);
    hRatio->GetYaxis()->SetTitle("Unfolded/MC");

    hRatio->SetMinimum(0.5);
    hRatio->SetMaximum(1.5);

    setXaxisTitle(hRatio, var, useAxis);

    hRatio->Draw("p9e same");
    TH1* sysBand_ratio_forData = NULL;
    TH1* sysBand_ratio_forMC = NULL;
    TLegend* leg_sys = NULL;
    if(sysName != "")
    {
        leg_sys = new TLegend(0.5, 0.75, 0.95, 0.85,"","brNDC");

        sysBand_ratio_forData = getUnfoldedSystematicBand(var, steering, true, sysName, hData, hDY, hRatio, divBinWidth, true, false, nthMassBin);
        sysBand_ratio_forData->SetFillColorAlpha(kBlack,0.8);
        sysBand_ratio_forMC = getUnfoldedSystematicBand(var, steering, true, sysName, hData, hDY, hRatio, divBinWidth, true, true, nthMassBin);
        sysBand_ratio_forMC->SetFillColorAlpha(kRed,0.2);

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

        sysBand_ratio_forData->SetFillStyle(3001);
        sysBand_ratio_forData->SetMarkerSize(0.);
        sysBand_ratio_forData->Draw("E2 same");
        sysBand_ratio_forMC->SetFillStyle(3001);
        sysBand_ratio_forMC->SetMarkerSize(0.);
        sysBand_ratio_forMC->Draw("E2 same");

        leg_sys->SetTextFont(43);
        leg_sys->SetTextSize(100);
        leg_sys->SetFillStyle(0); // transparent
        leg_sys->SetBorderSize(0);
        leg_sys->AddEntry(sysBand_ratio_forMC, sysName, "F");
        leg_sys->Draw();

    }

    TLine* l_ = new TLine(hRatio->GetXaxis()->GetXmin(),1,hRatio->GetXaxis()->GetXmax(),1);
    l_->SetLineColor(kRed);
    l_->Draw("same");
    l_->SetLineStyle(3);

    // Save canvas
    c_out->cd();
    c_out->SaveAs(outName!=""?output_baseDir+outName+var+".png":output_baseDir+"unfolded_"+var+sysName+".png");

    return c_out;
}

TCanvas* ISRUnfold::drawAcceptCorrHists(TString var, TString filePath, TString binDef, TString steering, bool useAxis, TString sysName, TString outName, int nthMassBin, bool divBinWidth)
{
    double histMin = 5e-1; 

    setTDRStyle();
    writeExtraText = true;
    extraText  = "Work in progress";
    gStyle->SetLineWidth(3.);
    gStyle->SetFrameLineWidth(3.);
    gROOT->ForceStyle();

    TH1::AddDirectory(kFALSE);
    // cout << "ISRUnfold::drawFoldedHists, Draw plot!" << endl;
    binDef = "_FineCoarse"; // FIXME

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
    hRatio = (TH1*) hData->Clone("hRatio");

    TCanvas* c_out = new TCanvas("acceptance_corrected_"+var, "acceptance_corrected_"+var, 3200, 2800);
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

    hData->SetTitle("");
    hData->SetStats(false);
    hData->GetXaxis()->SetMoreLogLabels(true);
    hData->Draw("p9e");
    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(3.2);
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
        sysBand_forMC->Draw("E2 same");
    }

    TLine massEdgeLine;
    for(int ibin = 2; ibin < hData->GetNbinsX()+1; ibin++)
    {
        massEdgeLine.SetLineStyle(1);
        massEdgeLine.SetLineColor(kWhite);
        massEdgeLine.DrawLine(hData->GetXaxis()->GetBinLowEdge(ibin), histMin, hData->GetXaxis()->GetBinLowEdge(ibin), hDY->GetBinContent(ibin));
    }

    TLegend* leg = new TLegend(0.5, 0.7, 0.9, 0.9,"","brNDC");
    //leg->SetNColumns(2);
    leg->SetTextFont(43);
    leg->SetTextSize(80);
    leg->SetFillStyle(0); // transparent
    leg->SetBorderSize(0);
    leg->AddEntry(hData, "Unfolded data", "pl");
    leg->AddEntry(hDY, "Drell-Yan (MG5_aMC@NLO)", "F");

    hData->Draw("p9e same");
    pad1->RedrawAxis();

    leg->Draw();

    int iPeriod_ = 4;
    if(year == 2017)
        iPeriod_ = 5;
    if(year == 2018)
        iPeriod_ = 6;
    CMS_lumi(pad1, iPeriod_, 11);

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

    hRatio->SetStats(false);
    hRatio->Divide(hDY);
    hRatio->Draw("p9e");
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(3.2);
    hRatio->SetLineColor(kBlack);
    hRatio->GetYaxis()->SetTitle("Unfolded/MC");

    hRatio->SetMinimum(0.5);
    hRatio->SetMaximum(1.5);

    setXaxisTitle(hRatio, var, useAxis);

    hRatio->Draw("p9e same");

    TH1* sysBand_ratio_forData = NULL;
    TH1* sysBand_ratio_forMC = NULL;
    TLegend* leg_sys = NULL;

    if(sysName != "")
    {
        leg_sys = new TLegend(0.5, 0.75, 0.95, 0.85,"","brNDC");

        sysBand_ratio_forData = getUnfAcceptSystematicBand(var, steering, true, sysName, hData, hDY, hRatio, divBinWidth, true, false, nthMassBin);
        sysBand_ratio_forData->SetFillColorAlpha(kBlack,0.8);
        sysBand_ratio_forMC = getUnfAcceptSystematicBand(var, steering, true, sysName, hData, hDY, hRatio, divBinWidth, true, true, nthMassBin);
        sysBand_ratio_forMC->SetFillColorAlpha(kRed,0.2);

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

        sysBand_ratio_forData->SetFillStyle(3001);
        sysBand_ratio_forData->SetMarkerSize(0.);
        sysBand_ratio_forData->Draw("E2 same");
        sysBand_ratio_forMC->SetFillStyle(3001);
        sysBand_ratio_forMC->SetMarkerSize(0.);
        sysBand_ratio_forMC->Draw("E2 same");

        leg_sys->SetTextFont(43);
        leg_sys->SetTextSize(100);
        leg_sys->SetFillStyle(0); // transparent
        leg_sys->SetBorderSize(0);
        leg_sys->AddEntry(sysBand_ratio_forMC, sysName, "F");
        leg_sys->Draw();

    }

    c_out->cd();
    if(sysName.Contains("Scale")) sysName = "Scale";
    c_out->SaveAs(outName!=""?output_baseDir+outName+var+sysName+".png":output_baseDir+"unfoldedAccept_"+var+sysName+".png");
    return c_out;
}

TH1* ISRUnfold::getUnfoldedSystematicBand(TString var, TString steering, bool useAxis, TString sysName, TH1* hData, TH1* hDY, TH1* hRatio, bool divBinWidth, bool isRatio, bool forMC, int nthMassBin)
{
    TH1* sysBand_ratio = NULL;
    
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
                if(forMC)
                    hSYS_temp = sysPtUnfold[sysName][sysMap[sysName][ith]]->GetBias("hDYMCPt_temp",0,0,steering,useAxis);
                else
                    hSYS_temp = sysPtUnfold[sysName][sysMap[sysName][ith]]->GetOutput("hUnfoldedPt_temp",0,0,steering,useAxis);
            }
            else
            {
                if(forMC)
                    hSYS_temp = sysMassUnfold[sysName][sysMap[sysName][ith]]->GetBias("hDYMCMass_temp",0,0,steering,useAxis);
                else
                    hSYS_temp = sysMassUnfold[sysName][sysMap[sysName][ith]]->GetOutput("hUnfoldedMass_temp",0,0,steering,useAxis);
            }
            if(divBinWidth)
            {
                divideByBinWidth(hSYS_temp, false);
            }

            if(forMC)
            {
                hRatio_temp.push_back((TH1*) hData->Clone("hRatio_temp"));
                hRatio_temp.at(ith)->Divide(hSYS_temp);
                if(ith==0) sysBand_ratio = (TH1*)hRatio_temp.at(ith)->Clone("sysUnfBand_ratio"+sysName);
            }
            else
            {
                hRatio_temp.push_back((TH1*) hSYS_temp->Clone("hRatio_temp"));
                hRatio_temp.at(ith)->Divide(hDY);
                if(ith==0) sysBand_ratio = (TH1*)hRatio_temp.at(ith)->Clone("sysUnfBand_ratio"+sysName);
            }

            for(int ibin = 1; ibin < sysBand_ratio->GetNbinsX()+1; ibin++)
            {
                double delta = fabs(hRatio_temp.at(ith)->GetBinContent(ibin) - hRatio->GetBinContent(ibin));
                if(sysName.Contains("FSR"))
                {
                    if(ith == 0) break;
                    delta = fabs(hRatio_temp.at(0)->GetBinContent(ibin) - hRatio_temp.at(1)->GetBinContent(ibin)); 
                }
                else
                {
                    if(ith != 0)
                    {
                        delta = delta > sysBand_ratio->GetBinError(ibin) ? delta : sysBand_ratio->GetBinError(ibin);
                    }
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

            delete hSYS_temp; // This could be data or MC systematic variation
        }
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
                if(ith==0) sysBand_ratio = (TH1*)hRatio_temp.at(ith)->Clone("sysUnfBand_ratio"+sysName);
            }
            else
            {
                hRatio_temp.push_back((TH1*) hSYS_temp->Clone("hRatio_temp"));
                hRatio_temp.at(ith)->Divide(hDY);
                if(ith==0) sysBand_ratio = (TH1*)hRatio_temp.at(ith)->Clone("sysUnfBand_ratio"+sysName);
            }

            for(int ibin = 1; ibin < sysBand_ratio->GetNbinsX()+1; ibin++)
            {
                double delta = fabs(hRatio_temp.at(ith)->GetBinContent(ibin) - hRatio->GetBinContent(ibin));
                if(sysName.Contains("FSR"))
                {
                    if(ith == 0) break;
                    delta = fabs(hRatio_temp.at(0)->GetBinContent(ibin) - hRatio_temp.at(1)->GetBinContent(ibin)); 
                }
                else
                {
                    if(ith != 0)
                    {
                        delta = delta > sysBand_ratio->GetBinError(ibin) ? delta : sysBand_ratio->GetBinError(ibin);
                    }
                }

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

            delete hSYS_temp; // This could be data or MC systematic variation
        }
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

TH1* ISRUnfold::getDetectorSystematicBand(TString var, TString filePath, TString dirName, TString steering, bool useAxis, TString sysName, TH1* hData, TH1* hDY, TH1* hMCtotal, TH1* hRatio, bool divBinWidth, bool isRatio, bool forMC, TString sysFilePath, int nthMassBin)
{

    TH1* sysBand_ratio = NULL;

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
            TH1* hSYS_temp = NULL; // This could be data or MC systematic variation
            TH1* hMCtotal_temp = NULL;
            TH1* hRatio_temp = NULL;

            // Dummy THStack, TLegend
            THStack* hsMC_temp;
            TLegend* leg = new TLegend(0.5, 0.7, 0.9, 0.9,"","brNDC");;
            hsMC_temp = new THStack("hsMC_temp", "hsMC_temp");

            TString tempHistName = DYHistName_;
            TString tempDirName = dirName;
            TString tempNewHistName = "Signal"+sysMap[sysName][ith];

            if(forMC)
            {
                if(!sysName.Contains("LepMom"))
                {
                    tempHistName = DYHistName_ + "_" + sysMap[sysName][ith];
                    // FIXME
                    if(sysName.Contains("Unfold") || sysName.Contains("ZptCorr"))
                        tempHistName = DYHistName_;
                }
                else
                {
                    tempDirName = dirName + "_" + sysMap[sysName][ith];
                }
            }
            else
            {
                tempHistName = dataHistName_;
                if(sysName.Contains("LepMom"))
                {
                    tempDirName = dirName + "_" + sysMap[sysName][ith];
                    tempNewHistName = "Data"+sysMap[sysName][ith];
                }
            }

            if(sysMap[sysName][ith] != "Nominal")
            {
                hSYS_temp = getRawHist(var, filePath, tempDirName, tempHistName, tempNewHistName, steering, useAxis, divBinWidth);
                hMCtotal_temp = (TH1*) hSYS_temp->Clone("hMCtotal_temp");
                hRatio_temp = (TH1*) hData->Clone("hRatio_temp");

                if(forMC)
                    hRatio_temp = (TH1*) hData->Clone("hRatio_temp");
                else
                    hRatio_temp = (TH1*) hSYS_temp->Clone("hRatio_temp");

                if(forMC)
                {
                    //cout << "call setTHStack" << endl;
                    setTHStack(var, filePath, tempDirName, *hsMC_temp, *hMCtotal_temp, *leg, steering, useAxis, sysMap[sysName][ith], divBinWidth);
                    delete hsMC_temp;
                    delete leg;
                }
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
                hRatio_temp->Divide(hMCtotal_temp);
                if(ith==0) sysBand_ratio = (TH1*)hRatio_temp->Clone("sysBand_ratio"+sysName);
            }
            else
            {
                hRatio_temp->Divide(hMCtotal);
                if(ith==0) sysBand_ratio = (TH1*)hRatio_temp->Clone("sysBand_ratio"+sysName);
            }

            // Update error
            for(int ibin = 1; ibin < sysBand_ratio->GetNbinsX()+1; ibin++)
            {
                double delta = fabs(hRatio_temp->GetBinContent(ibin) - hRatio->GetBinContent(ibin));
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

            delete hSYS_temp; // This could be data or MC systematic variation
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

void ISRUnfold::writeCutInfo(TPad* pad, TString var, int nthMassBin)
{
    pad->cd();

    double x_ = 0.2;
    double y_ = 0.7;
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
        lepton_cut.DrawLatexNDC(x_, y_-0.05, lepton_cut_info);
    }
    if(var.Contains("Mass"))
    {
        //dimass_cut.DrawLatexNDC(x_, y_, mass_cut_info);
        lepton_cut.DrawLatexNDC(x_, y_, lepton_cut_info);
    }
}

void ISRUnfold::setTHStack(TString var, TString filePath, TString dirName, THStack& hs, TH1& hMCtotal, TLegend& leg, TString steering, bool useAxis, TString sysName, bool divBinWidth)
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

        if(isFirstBkg)
        {
            if(sysName == "")
            {
                htemp = getRawHist(var, filePath, dirName, histName_, "h"+bkgNames[i], steering, useAxis, divBinWidth);
            }
            else
            {
                TString histPostfix = "histo_"+bkgNames[i] + "_" + sysName;
                if(sysName.Contains("LepMom") || sysName.Contains("ZpTCorrected")) // FIXME sysName -> sysPostfix
                {
                    //cout <<"sysName: " << sysName << endl;
                    histPostfix = "histo_"+bkgNames[i];
                    //cout <<"set histPostfix as " << histPostfix << endl;
                }
                if(bkgTypes[i] != "Fake" && sysName == "Fake") histPostfix = "histo_"+bkgNames[i];
                if(bkgTypes[i] == "Fake" && sysName != "Fake") histPostfix = "histo_"+bkgNames[i];
                htemp = getRawHist(var, filePath, dirName, histPostfix, "h"+bkgNames[i], steering, useAxis, divBinWidth);
            }
            isFirstBkg = false;
            nthBkg++;
        }
        else
        {
            if(sysName == "")
            {
                htemp->Add(getRawHist(var, filePath, dirName, histName_, "h"+bkgNames[i], steering, useAxis, divBinWidth));
            }
            else
            {
                TString histPostfix = "histo_"+bkgNames[i] + "_" + sysName;
                if(sysName.Contains("LepMom") || sysName.Contains("ZpTCorrected"))
                {
                    histPostfix = "histo_"+bkgNames[i];
                }
                if(bkgTypes[i] != "Fake" && sysName == "Fake") histPostfix = "histo_"+bkgNames[i];
                if(bkgTypes[i] == "Fake" && sysName != "Fake") histPostfix = "histo_"+bkgNames[i];
                htemp->Add(getRawHist(var, filePath, dirName, histPostfix, "h"+bkgNames[i], steering, useAxis, divBinWidth));
            }
            nthBkg++;
        }

        // This type of backgrounds all added, so add them to THStack
        if(nthBkg == bkgTypeN[bkgTypes[i]])
        {
            //cout << bkgTypes[i] << " " << bkgTypeN[bkgTypes[i]] << endl;

            if(var.Contains("Mass"))
            {
                htemp->GetXaxis()->SetRangeUser(massBinEdges.at(0), massBinEdges.at(5));
            }
            htemp->SetFillColor(bkgColors[bkgTypes[i]]);
            hs.Add(htemp);
            hMCtotal.Add(htemp);

            if(sysName == "")
            {
                //leg.AddEntry(htemp, bkgTypes[i], "F");
                tempMap_legend[bkgTypes[i]] = htemp;
                temp_bkgName.push_back(bkgTypes[i]);
            }

            isFirstBkg = true;
            nthBkg = 0;
        }
    }

    for(unsigned int i = temp_bkgName.size(); i > 0; i--)
    {
        leg.AddEntry(tempMap_legend[temp_bkgName.at(i-1)], temp_bkgName.at(i-1), "F");
    }

}
void ISRUnfold::doStatUnfold()
{
    //cout << "ISRUnfold::doStatUnfold() " << endl;
    if(regMode == TUnfold::kRegModeNone)
    {
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
    }
}

void ISRUnfold::doISRUnfold(bool doSys, bool doReg)
{
    //cout << "ISRUnfold::doISRUnfold!!" << endl;
    if(!doSys)
    {
        //cout << "Unfold without systematic" << endl;
        if(!doReg)
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
            nomPtUnfold->RegularizeBins(istart,1,iend-istart+1,TUnfoldV17::kRegModeCurvature);

            //TGraph *lcurve;
            //TSpline *logtaux,*logtauy,*logtaucurvature;
            //int ibest=nomPtUnfold->ScanLcurve(100,0,0,&lcurve,&logtaux,&logtauy,&logtaucurvature);
            double tauMin=0.;
            double tauMax=0.;
            nomPtUnfold->ScanLcurve(100,tauMin,tauMax,0);
        }
    }
    else
    {
        //cout << "Do systematic unfold!" << endl;
        //For systematic
        std::map<TString, std::vector<TString>>::iterator it = sysMap.begin();
        while(it != sysMap.end())
        {
            //cout << "Unfold for " << it->first << " systematic." << endl;
            int size = (it->second).size();
            //cout << size << " systematic variation exist." << endl;

            for(int i = 0; i < size; i++)
            {
                if(!doReg)
                {
                    //cout << "posfix: " << (it->second).at(i) << endl;
                    sysPtUnfold[it->first][(it->second).at(i)]->DoUnfold(0);
                    sysMassUnfold[it->first][(it->second).at(i)]->DoUnfold(0);
                }

                if(it->first == "PDF")
                {
                    fillMassPDFVariationHist(i+1); // i+1 since it starts from 1
                    fillPtPDFVariationHist(i+1);
                }
            }
            it++;
        }
    }// Unfold for systematic
}

void ISRUnfold::drawCorrelation(TString var, TString steering, bool useAxis, TString outName)
{
    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("c","c", 2400, 2400);
    c->cd();

    TH2* hCorrelation = NULL;

    if(var.Contains("Mass"))
    {
        hCorrelation=nomMassUnfold->GetRhoIJtotal("histRho", 0, 0, steering, useAxis);
    }
    else
    {
        hCorrelation=nomPtUnfold->GetRhoIJtotal("histRho", 0, 0, steering, useAxis);
    }

    hCorrelation->SetMinimum(-1.);
    hCorrelation->SetMaximum(1.);
    hCorrelation->Draw("COLZ");

    c->SaveAs(output_baseDir+"Correlation_"+var+"_"+outName+".png");
    delete hCorrelation;
}
void ISRUnfold::doAcceptCorr(TString filePath, TString binDef, bool doSys, TString outName)
{
    TFile* filein = new TFile(filePath);

    TH1* hFullPhaseMassMC_raw = NULL;
    TH1* hFiducialPhaseMassMC = NULL;
    TH1* hFullPhasePtMC_raw = NULL;
    TH1* hFiducialPhasePtMC = NULL;

    // Nominal acceptance
    // Mass
    hFullPhaseMassMC_raw = (TH1*) filein->Get("Acceptance/MassGen" + binDef + "/histo_DYJets");
    if(year==2016)
        hFullPhaseMassMC_raw->Add((TH1*) filein->Get("Acceptance/MassGen" + binDef + "/histo_DYJets10to50"));
    else
        hFullPhaseMassMC_raw->Add((TH1*) filein->Get("Acceptance/MassGen" + binDef + "/histo_DYJets10to50_MG"));

    hFiducialPhaseMassMC = nomMassUnfold->GetBias("hFiducialMass", 0, 0, "*[*]", false);
    hAcceptanceMass = (TH1*) hFullPhaseMassMC_raw->Clone("hAcceptanceMass");
    hAcceptanceMass->Divide(hFiducialPhaseMassMC);

    hFullPhaseMassData = nomMassUnfold->GetOutput("hFullPhaseMassData",0,0, "*[*]", false);
    hFullPhaseMassData->Multiply(hAcceptanceMass);

    // Pt
    hFullPhasePtMC_raw = (TH1*) filein->Get("Acceptance/PtGen" + binDef + "/histo_DYJets");
    if(year==2016)
        hFullPhasePtMC_raw->Add((TH1*) filein->Get("Acceptance/PtGen" + binDef + "/histo_DYJets10to50"));
    else
        hFullPhasePtMC_raw->Add((TH1*) filein->Get("Acceptance/PtGen" + binDef + "/histo_DYJets10to50_MG"));

    hFiducialPhasePtMC = nomPtUnfold->GetBias("hFiducialPt", 0, 0, "*[*]", false);
    hAcceptancePt = (TH1*) hFullPhasePtMC_raw->Clone("hAcceptancePt");
    hAcceptancePt->Divide(hFiducialPhasePtMC);

    hFullPhasePtData = nomPtUnfold->GetOutput("hFullPhasePtData",0,0, "*[*]", false);
    hFullPhasePtData->Multiply(hAcceptancePt);

    // Draw acceptance
    drawAcceptance("Mass", hFiducialPhaseMassMC, outName);
    drawAcceptance("Pt", hFiducialPhasePtMC, outName);

    if(doSys)
    {
        std::map<TString, std::vector<TString>>::iterator it = sysMap.begin();
        while(it != sysMap.end())
        {
            int size = (it->second).size();
            for(int i = 0; i < size; i++)
            {
                hSysFullPhaseMassData[it->first][(it->second).at(i)] = sysMassUnfold[it->first][(it->second).at(i)]->GetOutput("hFullPhaseMassData"+it->first+(it->second).at(i),0,0, "*[*]", false);
                hSysFullPhasePtData[it->first][(it->second).at(i)]   = sysPtUnfold[it->first][(it->second).at(i)]->GetOutput("hFullPhasePtData"+it->first+(it->second).at(i),0,0, "*[*]", false);

                // Use different acceptance for PDF AlphaS Scale etc
                if((it->first).Contains("Scale") || (it->first).Contains("PDF") || (it->first).Contains("AlphaS"))
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

                    hSysFullPhaseMassData[it->first][(it->second).at(i)]->Multiply(hAcceptanceMass_sys);
                    hSysFullPhaseMassMC[it->first][(it->second).at(i)] = hFullPhaseMassMC_raw_sys;

                    delete hAcceptanceMass_sys;

                    // For pt
                    TH1* hFullPhasePtMC_raw_sys = (TH1*) filein->Get("Acceptance/PtGen" + binDef + "/histo_DYJets_"+(it->second).at(i));
                    if(year==2016)
                        hFullPhasePtMC_raw_sys->Add((TH1*) filein->Get("Acceptance/PtGen" + binDef + "/histo_DYJets10to50_"+(it->second).at(i)));
                    else
                        hFullPhasePtMC_raw_sys->Add((TH1*) filein->Get("Acceptance/PtGen" + binDef + "/histo_DYJets10to50_MG_"+(it->second).at(i)));

                    TH1* hFiducialPhasePtMC_sys = sysPtUnfold[it->first][(it->second).at(i)]->GetBias("hFiducialPt_sys", 0, 0, "*[*]", false);
                    TH1* hAcceptancePt_sys = (TH1*) hFullPhasePtMC_raw_sys->Clone("hAcceptancePt_sys");
                    hAcceptancePt_sys->Divide(hFiducialPhasePtMC_sys);

                    hSysFullPhasePtData[it->first][(it->second).at(i)]->Multiply(hAcceptancePt_sys);
                    hSysFullPhasePtMC[it->first][(it->second).at(i)] = hFullPhasePtMC_raw_sys;

                    if((it->first).Contains("PDF"))
                    {
                        fillMassPDFVariationHist_Accept(i+1);
                        fillPtPDFVariationHist_Accept(i+1);
                    }

                    delete hAcceptancePt_sys;
                }
                else
                {
                    // Use nominal acceptance 
                    hSysFullPhaseMassData[it->first][(it->second).at(i)]->Multiply(hAcceptanceMass);
                    hSysFullPhasePtData[it->first][(it->second).at(i)]->Multiply(hAcceptancePt);
                    hSysFullPhaseMassMC[it->first][(it->second).at(i)] = hFullPhaseMassMC_raw;
                    hSysFullPhasePtMC[it->first][(it->second).at(i)]   = hFullPhasePtMC_raw;
                }


            }
            it++;
        }
    }

    delete hFiducialPhaseMassMC;
    delete hFiducialPhasePtMC;
}

void ISRUnfold::drawAcceptance(TString var, TH1* hMC, TString outName)
{
    setTDRStyle();
    writeExtraText = true;
    extraText  = "Simulation";
    gStyle->SetLineWidth(3.);
    gStyle->SetFrameLineWidth(3.);
    gROOT->ForceStyle();

    TH1* hFullPhase = NULL;
    TH1* hFiducialPhase = NULL;
    TH1* hAcceptance = NULL;

    if(var.Contains("Mass"))
    {
        hFullPhase = mass_binning_Gen->ExtractHistogram("hData", hFullPhaseMassData, 0, true, "mass[UO];pt[UOC0]");
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

        TH1* hAccept_corr_mass = mass_binning_Gen->ExtractHistogram("hData_accept", hFullPhaseMassData, 0, kTRUE, "mass[UO];pt[UOC0]");
        for(int ibin = 0; ibin < nMassBin; ibin++)
        {
            TString ibinMass;
            ibinMass.Form("%d", ibin);

            hFullPhase = pt_binning_Gen->ExtractHistogram("hData", hFullPhasePtData, 0, true, "pt[UO];mass[UOC"+ibinMass+"]");
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

    hratio->SetMinimum(1e-3);
    hratio->SetMaximum(1.0);

    setXaxisTitle(hratio, var, true);
    hratio->Draw("p9e same");

    c_out->cd();
    c_out->SaveAs(outName!=""?output_baseDir+plotName+"_"+var+"_"+outName+".png":output_baseDir+plotName+"_"+var+".png");
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

int ISRUnfold::setMeanMass(TString filePath, TString dirName)
{
    //cout << "ISRUnfold::setMeanMass()   Save mean of dilepton..." << endl;
    TFile* filein = new TFile(filePath);
    TString DYHistName_ = "histo_DYJetsToMuMu";
    if(channel_name == "electron")
    {
        DYHistName_  = "histo_DYJetsToEE";
    }

    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;
    const Double_t* massBins = temp_tvecd->GetMatrixArray();

    TUnfoldDensityV17* p_unfold = NULL;
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
        TUnfoldDensityV17* p_unfold = NULL;

        //cout << "Unfold for " << it->first << " systematic." << endl;
        int size = (it->second).size();
        //cout << size << " systematic variation exist." << endl;

        for(int i = 0; i < size; i++)
        {
            p_unfold = sysMassUnfold[it->first][(it->second).at(i)];
            TH1* hunfolded_mass =  p_unfold->GetOutput("hunfolded_mass",0,0,"mass[UO];pt[UOC0]",kTRUE);

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
    cout << "------------------------------------- Systematic Uncertainty for ISR analysis ---------------------------------------" << endl;
    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;

    std::streamsize ss = std::cout.precision();
    std::cout.precision(2);
    std::cout.setf( std::ios::fixed, std:: ios::floatfield );

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


            if(it->first != "PDF")
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
                error_mass = meanMassPDFVariation.at(ibin)->GetRMS();
                error_pt = meanPtPDFVariation.at(ibin)->GetRMS();
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
    std::cout.precision(ss);
}

void ISRUnfold::setSysError_Accept()
{

    cout << "------------------------------------- Systematic Uncertainty for ISR analysis (Acceptance corrected)---------------------------------------" << endl;
    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;

    std::streamsize ss = std::cout.precision();
    std::cout.precision(2);
    std::cout.setf( std::ios::fixed, std:: ios::floatfield );

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
            if(it->first != "PDF")
            {
            for(int i = 0; i < size; i++)
            {
                double temp_error_mass = fabs(meanMass_data_accept_sysVariation[it->first][(it->second).at(i)].at(ibin) - meanMass_data_acc_corrected.at(ibin));
                double temp_error_pt = fabs(meanPt_data_accept_sysVariation[it->first][(it->second).at(i)].at(ibin) - meanPt_data_acc_corrected.at(ibin));
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
                error_mass = meanMassPDFVariation_Accept.at(ibin)->GetRMS();
                error_pt = meanPtPDFVariation_Accept.at(ibin)->GetRMS();
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
    std::cout.precision(ss);
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

int ISRUnfold::setMeanPt(TString filePath, TString dirName)
{
    //cout << "ISRUnfold::setMeanPt()   Save mean of dilepton momentum..." << endl;    TFile* filein = new TFile(filePath);
    TFile* filein = new TFile(filePath);
    TString DYHistName_ = "histo_DYJetsToMuMu";
    if(channel_name == "electron")
    {
        DYHistName_  = "histo_DYJetsToEE";
    }

    // Find number of mass bins
    const TVectorD* temp_tvecd = pt_binning_Gen->GetDistributionBinning(1);
    int nMassBin = temp_tvecd->GetNrows() - 1;

    TUnfoldDensityV17* p_unfold = NULL;
    p_unfold = nomPtUnfold;

    // Save mean pt
    for(int i = 0; i < nMassBin; i++)
    {
        TString ibinMass;
        ibinMass.Form("%d", i);

        // Get detector level MC
        if(filePath != "")
        {
            cout << filePath << endl;
            TH1* hDY = getRawHist("Pt_FineCoarse", filePath, "Detector", DYHistName_, "Signal", "pt[UO];mass[UOC"+ibinMass+"]", true, false);
            cout << "detector mc, pt: " << hDY->GetMean() << " err: " << hDY->GetMeanError() << endl;
            cout << "n bins: " << hDY->GetNbinsX() << endl;
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
        TUnfoldDensityV17* p_unfold = NULL;

        //cout << "Unfold for " << it->first << " systematic." << endl;
        int size = (it->second).size();
        //cout << size << " systematic variation exist." << endl;

        for(int i = 0; i < size; i++)
        {
            p_unfold = sysPtUnfold[it->first][(it->second).at(i)];
            // Save mean pt
            for(int j = 0; j < nMassBin; j++)
            {
                TString ibinMass;
                ibinMass.Form("%d", j);

                TH1* hunfolded_pt;

                // Get histograms to set mean values
                hunfolded_pt = p_unfold->GetOutput("hunfolded_pt_temp",0,0,"pt[UO];mass[UOC"+ibinMass+"]",kTRUE);
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
        //cout << dirName+"/"+var+"/"+histName << endl;
        TH1* raw_hist = (TH1*)filein->Get(dirName+"/"+var+"/"+histName);
        //cout << "# of bins: " << raw_hist->GetNbinsX() << endl;
        if(histName.Contains("DYJetsTo") && !histName.Contains("Tau"))
        {
            histName.ReplaceAll("DYJetsTo", "DYJets10to50To");
            raw_hist->Add((TH1*)filein->Get(dirName+"/"+var+"/"+histName));
        }

        if(var.Contains("Pt"))
        {
            hist = pt_binning_Rec->ExtractHistogram(outHistName, raw_hist, 0, useAxis, steering);
            //cout << "steering: " << steering << endl;
            //cout << "# of bins: " << hist->GetNbinsX() << endl;
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

        c->SaveAs(output_baseDir+"MeanPtStat_" + nth + year_ + ".png");
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

        c->SaveAs(output_baseDir+"MeanMassStat_" + nth + year_ + ".png");
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

        c->SaveAs(output_baseDir+"MeanPtPDF_" + nth + year_ + ".png");
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

        c->SaveAs(output_baseDir+"MeanMassPDF_" + nth + year_ + ".png");
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

        c->SaveAs(output_baseDir+"MeanPt_"+sysName+"_"+nth+"_"+year_+".png");
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

        c->SaveAs(output_baseDir+"MeanMass_"+sysName+"_"+nth+"_"+year_+".png");
    }
}

void ISRUnfold::drawSystematics(TString var)
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
    for(int i = 0; i < nMassBin; i++)
    {
        v_bins.push_back(i+1);
    }

    map<TString, TGraph*> map_sys_graph;
    std::map<TString, std::vector<TString>>::iterator it = sysMap.begin();
    while(it != sysMap.end())
    {
        if(var.Contains("Pt"))
        {
            map_sys_graph[it->first] = new TGraph(nMassBin, &v_bins[0], &meanPt_data_folded_rel_systematic[it->first][0]);
        }
        else
        {
            map_sys_graph[it->first] = new TGraph(nMassBin, &v_bins[0], &meanMass_data_folded_rel_systematic[it->first][0]);
        }
        it++;
    }

    // Create canvas
    TCanvas* c_out = new TCanvas("relative_uncertainty_"+var, "relative_uncertainty_"+var, 3000, 1800);
    c_out->SetGridy(1);
    c_out->SetGridx(1);
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
        if(markerColor == 10) markerColor = 1;
        if(marker == 24) marker = 29;
        if(marker == 30) marker = 33;
        if(marker == 35) marker = 39;
        if(marker == 40) marker = 41;
        if(marker == 42) marker = 47;
        if(marker > 47) marker = 20;
        map_sys_graph[it->first]->SetMarkerStyle(marker);
        map_sys_graph[it->first]->SetMarkerSize(markerSize);
        map_sys_graph[it->first]->SetMarkerColor(markerColor==5?46:markerColor);
        map_sys_graph[it->first]->SetLineColor(markerColor==5?46:markerColor);
        markerColor++;
        marker++;

        if(first_draw)
        {
            //map_sys_graph[it->first]->SetTitleOffset(0.02);
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
            first_draw = false;
        }
        else
        {
            map_sys_graph[it->first]->Draw("PC SAME");
        }
        leg->AddEntry(map_sys_graph[it->first], it->first, "pl");
        it++;
    }
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

    leg->Draw();
    c_out->cd();
    c_out->SaveAs(output_baseDir+"Systematic_"+var+".png");

    delete c_out;
}

void ISRUnfold::drawSystematics_Acceptance(TString var)
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
    for(int i = 0; i < nMassBin; i++)
    {
        v_bins.push_back(i+1);
    }

    map<TString, TGraph*> map_sys_graph;
    std::map<TString, std::vector<TString>>::iterator it = sysMap.begin();
    while(it != sysMap.end())
    {
        if(var.Contains("Pt"))
        {
            map_sys_graph[it->first] = new TGraph(nMassBin, &v_bins[0], &meanPt_data_accept_rel_systematic[it->first][0]);
        }
        else
        {
            map_sys_graph[it->first] = new TGraph(nMassBin, &v_bins[0], &meanMass_data_accept_rel_systematic[it->first][0]);
        }
        it++;
    }

    // Create canvas
    TCanvas* c_out = new TCanvas("relative_uncertainty_"+var, "relative_uncertainty_"+var, 3000, 1800);
    c_out->SetGridy(1);
    c_out->SetGridx(1);
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
        if(markerColor == 10) markerColor = 1;
        if(marker == 24) marker = 29;
        if(marker == 30) marker = 33;
        if(marker == 35) marker = 39;
        if(marker == 40) marker = 41;
        if(marker == 42) marker = 47;
        if(marker > 47) marker = 20;
        map_sys_graph[it->first]->SetMarkerStyle(marker);
        map_sys_graph[it->first]->SetMarkerSize(markerSize);
        map_sys_graph[it->first]->SetMarkerColor(markerColor==5?46:markerColor);
        map_sys_graph[it->first]->SetLineColor(markerColor==5?46:markerColor);
        markerColor++;
        marker++;

        if(first_draw)
        {
            //map_sys_graph[it->first]->SetTitleOffset(0.02);
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
            first_draw = false;
        }
        else
        {
            map_sys_graph[it->first]->Draw("PC SAME");
        }
        leg->AddEntry(map_sys_graph[it->first], it->first, "pl");
        it++;
    }
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
    map_sys_graph["Stat."]->SetLineWidth(10);
    map_sys_graph["Stat."]->Draw("PC SAME");
    leg->AddEntry(map_sys_graph["Stat."], "Statistical", "l");

    leg->Draw();
    c_out->cd();
    c_out->SaveAs(output_baseDir+"SystematicAccept_"+var+".png");

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
