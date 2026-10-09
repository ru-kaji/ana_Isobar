// My header files
#include "CorTreeAnalyzer.h"
#include "EventClass.h"

// ROOT header files 
#include "TSystem.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TBox.h"
#include "TLegend.h"
#include "TLatex.h"

// C/C++ heder
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

void calcPurity_isobar(TString filename, TString particle = "L") {

    TString pidStr;
    TString histName;
    TString plotTitle;
    TString xAxisLabel;
    TString outPdf;
    TString particleNameLatex;

    if (particle == "L") {
        pidStr            = "L";
        plotTitle         = "#Lambda Mass";
        xAxisLabel        = "M_{p#pi^{-}} (GeV/c^{2})";
        outPdf            = "Purity_isobar_L.pdf";
        particleNameLatex = "#Lambda Mass";
    } else if (particle == "A") {
        pidStr             = "A"; 
        plotTitle         = "#bar{#Lambda} Mass"; // #bar{...} で上線が付きます
        xAxisLabel        = "M_{#bar{p}#pi^{+}} (GeV/c^{2})";   // 崩壊娘粒子は反陽子とπ+になります
        outPdf            = "Purity_isobar_A.pdf";
        particleNameLatex = "#bar{#Lambda} Mass";
    } else {
        std::cerr << "Error: Unknown particle type!" << std::endl;
        return;
    }

    // 1. Open the file and get the histogram
    TFile *file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    // 2. Define the fit ranges and signal region
    double pdgMass = 1.115683; // Lambda PDG mass (GeV/c^2)
    double fitMin = 1.09;
    double fitMax = 1.13;
    
    // 3. Set up the Fit Function: Double Gaussian + Pol2 Background
    // p0: Constant of Gaussian 1 (Core)
    // p1: Mean of Gaussians
    // p2: Sigma of Gaussian 1 (Core)
    // p3: Constant of Gaussian 2 (Tail)
    // p4: Sigma of Gaussian 2 (Tail)
    // p5, p6, p7: Pol2 parameters for Background
    TF1 *fTotal = new TF1("fTotal", "gaus(0) + [3]*exp(-0.5*((x-[1])/[4])^2) + [5]+[6]*(x-1.115)", fitMin, fitMax);
    //TF1 *fTotal = new TF1("fTotal", "gaus(0) + [3]*exp(-0.5*((x-[1])/[4])^2) + [5]+[6]*(x-1.115)+[7]*pow(x-1.115, 2)+ [8]*pow(x-1.115, 3)", fitMin, fitMax);
    fTotal->SetNpx(2000); // Set high drawing resolution for narrow peaks

    TCanvas *c_purity = new TCanvas("c_purity", "#Lambda Purity Calculation", 1400, 400);
    c_purity->Divide(3, 1);

    const char* centName[8]  = {"0-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%"};
    TH1F* hMass[8] = {nullptr};

    for(int icent=0; icent<8; icent++){
        hMass[icent] = (TH1F*)file->Get(Form("hinvM_%s_%s", pidStr.Data(), centName[icent]));
        if (!hMass[icent]) {
            std::cerr << "Error: Cannot find histogram hinvM in " << filename << std::endl;
            continue;
        }
    }

    TH1F* hMass_0to80  = (TH1F*)hMass[0]->Clone("hMass_0to80");  hMass_0to80->Reset();
    TH1F* hMass_10to80 = (TH1F*)hMass[0]->Clone("hMass_10to80"); hMass_10to80->Reset();
    TH1F* hMass_20to80 = (TH1F*)hMass[0]->Clone("hMass_20to80"); hMass_20to80->Reset();

    for(int icent=0; icent<8; icent++){
        if(!hMass[icent]) continue;
        hMass_0to80->Add(hMass[icent]);
        if(icent>=1) hMass_10to80->Add(hMass[icent]);
        if(icent>=2) hMass_20to80->Add(hMass[icent]);
    }

    TH1F* hTarget[3] = {hMass_0to80, hMass_10to80, hMass_20to80};
    TString targetName[3] = {"0-80%", "10-80%", "20-80%"};

    for(int i = 0; i < 3; i++){
        TH1F *hFit = (TH1F*)hTarget[i]->Clone(Form("hFit_Purity_%s", targetName[i].Data()));
        hFit->SetDirectory(0); // detach from file

        if (hFit->GetEntries() < 10 || hFit->GetMaximum() <= 0) {
            std::cerr << "Warning: Low statistics for " << targetName[i].Data() << ". Skipping fit." << std::endl;
            delete hFit;
            continue;
        };
    
        // Initial parameter guesses
        double bgHeight = hFit->GetBinContent(hFit->FindBin(fitMin));
        double maxCounts = hFit->GetMaximum();
        double sigHeight = maxCounts - bgHeight;
        
        fTotal->SetParameter(0, sigHeight*0.9);       // Core height
        fTotal->SetParLimits(0, sigHeight*0.85, sigHeight);
        fTotal->SetParameter(1, 1.115683);         // Mean
        fTotal->SetParLimits(1, pdgMass - 0.01, pdgMass + 0.01);
        
        fTotal->SetParameter(2, 0.0007);          // Core sigma (narrow)
        fTotal->SetParLimits(2, 0.0005, 0.003);
        
        // Disable Tail Gaussian by fixing its amplitude (p3) to 0
        // To enable later, change FixParameter to SetParameter/SetParLimits
        fTotal->ReleaseParameter(3);
        fTotal->ReleaseParameter(4);
        fTotal->SetParameter(3, sigHeight*0.3); 
        fTotal->SetParLimits(3, sigHeight*0.1, sigHeight*0.3);             // Tail height
        fTotal->SetParameter(4, 0.002);    
        fTotal->SetParLimits(4, 0.0025, 0.02);
        
        // Initial guess for background
        fTotal->SetParameter(5,  2.5e7); 
        //fTotal->SetParLimits(5, bgHeight*0.8, bgHeight*1.1);
        //fTotal->ReleaseParameter(5);
        fTotal->SetParameter(6, -1.7e3);
        fTotal->SetParLimits(6, -1e+7, 1.0);
        //fTotal->ReleaseParameter(6);

        // 各パラメータの取得
        double p0_core = fTotal->GetParameter(0); // Core Constant
        double p2_sig1 = fTotal->GetParameter(2); // Core Sigma
        double p3_tail = fTotal->GetParameter(3); // Tail Constant
        double p4_sig2 = fTotal->GetParameter(4); // Tail Sigma

        // 面積（Yield）に比例する値を計算
        double area_core = p0_core * p2_sig1; 
        double area_tail = p3_tail * p4_sig2;
        double area_total = area_core + area_tail;

        // 面積で重み付けした分散から、Effective Sigma を計算
        double sgm_eff = std::sqrt((area_core * std::pow(p2_sig1, 2) + area_tail * std::pow(p4_sig2, 2)) / area_total);

        // 4. Perform the Fit
        c_purity->cd(i+1);
        hFit->SetStats(0);
        hFit->SetMarkerStyle(20);
        hFit->SetMarkerSize(0.7);
        hFit->SetLineColor(kBlack);
        hFit->GetXaxis()->SetRangeUser(1.08, 1.16);
        hFit->SetTitle(Form("%s %s", plotTitle.Data(), targetName[i].Data()));
        hFit->GetXaxis()->SetTitle(xAxisLabel);
        hFit->GetXaxis()->SetTitleSize(0.04);
        hFit->GetYaxis()->SetTitle("Counts");
        hFit->GetYaxis()->SetTitleOffset(0.9);
        hFit->GetYaxis()->SetTitleSize(0.05);
        hFit->Draw("PE");

        hFit->Fit(fTotal, "R");
        TFitResultPtr r = hFit->Fit(fTotal, "SR");

        // 5. Extract Signal and Background functions
        TF1 *fSig1 = new TF1("fSig1", "gaus", fitMin, fitMax);
        fSig1->SetParameters(fTotal->GetParameter(0), fTotal->GetParameter(1), fTotal->GetParameter(2));
        fSig1->SetLineColor(kGreen+3);
        fSig1->SetLineStyle(2);
        fSig1->Draw("same");

        TF1 *fSig2 = new TF1("fSig2", "gaus", fitMin, fitMax);
        fSig2->SetParameters(fTotal->GetParameter(3), fTotal->GetParameter(1), fTotal->GetParameter(4));
        fSig2->SetLineColor(kCyan);
        fSig2->SetLineStyle(2);
        fSig2->Draw("same");

        
        TF1 *fBg = new TF1("fBg", "[0]+[1]*(x-1.115)", fitMin, fitMax);
        fBg->SetParameters(fTotal->GetParameter(5), fTotal->GetParameter(6));
        fBg->SetLineColor(kBlue);
        fBg->SetLineStyle(2);
        fBg->Draw("same");


        // Usually signal region is +/- n*sigma around the lambda mass. 
        double sgm = fTotal->GetParameter(2);
        double sigMin = pdgMass - 3.0 * sgm;
        double sigMax = pdgMass + 3.0 * sgm;

        // 6. Calculate Integrals within the Signal Window
        int binMin = hFit->FindBin(sigMin);
        int binMax = hFit->FindBin(sigMax - 1e-6); // ※ sigMax がちょうどBinの境界上にある場合、1つ上のBinにならないように微小値を引く
        double binWidth = hFit->GetBinWidth(1);

        double exactSigMin = hFit->GetBinLowEdge(binMin);
        double exactSigMax = hFit->GetBinLowEdge(binMax) + binWidth;
        
        // Integral of Fit Function
        double S1 = fSig1->Integral(exactSigMin, exactSigMax) / binWidth;
        double S2 = fSig2->Integral(exactSigMin, exactSigMax) / binWidth;
        double S  = S1 + S2;
        double B  = fBg->Integral(exactSigMin, exactSigMax) / binWidth;
        double S_plus_B_fit = S + B;

        // Integral of Data
        double S_plus_B_data = hFit->Integral(binMin, binMax);

        // Calculate Purity
        double purity_fit  = S / S_plus_B_fit;
        double purity_data = S / S_plus_B_data; // Using fit for S, data for total
        
        double chi2 = fTotal->GetChisquare();
        double ndf  = fTotal->GetNDF();
        double chi2ndf = (ndf != 0) ? chi2 / ndf : -1.0;

        // 7. Visualization
        double maxy = hFit->GetMaximum();
        TLine *lMin = new TLine(sigMin, 0, sigMin, maxy*0.8);
        TLine *lMax = new TLine(sigMax, 0, sigMax, maxy*0.8);
        lMin->SetLineStyle(3); lMax->SetLineStyle(3);
        lMin->Draw("same"); lMax->Draw("same");
        
        // Add Legend
        TLegend *leg = new TLegend(0.6, 0.15, 0.9, 0.45);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.03);
        leg->AddEntry(hFit, "Data", "pe");
        leg->AddEntry(fTotal, "Total Fit", "l");
        leg->AddEntry(fSig1, "Core Gaussian", "l");
        leg->AddEntry(fSig2, "Tail Gaussian", "l");
        leg->AddEntry(fBg, "Background (pol1)", "l");
        leg->Draw();
        
        TLatex tex;
        tex.SetNDC();
        tex.SetTextSize(0.03); // テキストサイズを少し調整
        //tex.DrawLatex(0.58, 0.87, particleNameLatex);
        tex.DrawLatex(0.55, 0.83, Form("Signal Window: [%.4f, %.4f]", sigMin, sigMax));
        //tex.DrawLatex(0.58, 0.79, Form("S (Fit) = %.0f", S));
        //tex.DrawLatex(0.58, 0.75, Form("B (Fit) = %.0f", B));
        //tex.DrawLatex(0.58, 0.71, Form("S+B (Data) = %.0f", S_plus_B_data));
        tex.DrawLatex(0.55, 0.78, Form("Purity (S_{fit}/S+B_{fit}) = %.4f", purity_fit));
        tex.DrawLatex(0.55, 0.73, Form("Purity (S_{fit}/S+B_{data}) = %.4f", purity_data));
        tex.DrawLatex(0.55, 0.68, Form("#chi^{2}/ndf = %.1f / %d = %.2f", chi2, (int)ndf, chi2ndf)); // chi2表示を追加
        
    } //end of i roop
    
    c_purity->SaveAs(outPdf);
    file->Close();
}
