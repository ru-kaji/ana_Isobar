#ifndef CFCALCULATOR_H
#define CFCALCULATOR_H

#include "TStyle.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"
#include <iostream>
#include <cmath>

class CFAnalyzer {

private:
    TH1 *hSame, *hMix;
    TH1 *hRaw, *hNorm;
    TH1 *hNorm_SB, *hNorm_BS, *hNorm_BB;
    TH1 *hNorm_Uncorr;
    TH1 *hNorm_Corr;
    TString fName;
    double mP_L, mP_A; 

    void SetStyle() {  
        if(hNorm) {  
            hNorm->SetXTitle("Q (GeV/c)");
            hNorm->SetYTitle("C(Q)");
            hNorm->SetMarkerStyle(20);
            hNorm->SetMarkerColor(kBlack);
            hNorm->SetLineColor(kRed);
            hNorm->SetMarkerSize(0.5);
            hNorm->SetMinimum(0.0);
            hNorm->SetMaximum(2.0);
        }
    }

public:
    // コンストラクタとデストラクタ
    CFAnalyzer(TH1* same, TH1* mix, TString name, double pL = 0.8740, double pA = 0.8894) {//初期化
        fName = name;
        mP_L = pL;
        mP_A = pA;
        hRaw = nullptr;
        hNorm = nullptr;
        hNorm_Uncorr = nullptr;
        hNorm_Corr = nullptr;
        hSame = nullptr;
        hMix = nullptr;

        if(same && mix){
            hSame = (TH1*)same->Clone(fName + "_same");
            hMix  = (TH1*)mix->Clone(fName + "_mix");
            hSame->Sumw2();
            hMix->Sumw2();
        } else {
            //std::cout << "Error: Input histograms are null!" << std::endl;
        }
    }
    virtual ~CFAnalyzer() {
        if(hSame) delete hSame;
        if(hMix)  delete hMix;
        if(hRaw)  delete hRaw;
        if(hNorm) delete hNorm;
        if(hNorm_Uncorr) delete hNorm_Uncorr;
        if(hNorm_Corr) delete hNorm_Corr;
    }

    void Calculate(double normMin, double normMax) {
        if(!hSame || !hMix) return;

        // 複数回呼ばれた時のための初期化
        if(hRaw) delete hRaw;
        if(hNorm) delete hNorm;
        if(hNorm_Uncorr) delete hNorm_Uncorr;
        if(hNorm_Corr) delete hNorm_Corr;

        hRaw = (TH1*)hSame->Clone(fName + "_Raw");
        hRaw->SetDirectory(nullptr);
        hRaw->Divide(hMix);

        int bin_min = hSame->FindBin(normMin);
        int bin_max = hSame->FindBin(normMax);
        double int_Same = hSame->Integral(bin_min, bin_max);
        double int_Mix  = hMix->Integral(bin_min, bin_max);

        hNorm = (TH1*)hRaw->Clone(fName + "_Norm");
        hNorm->SetDirectory(nullptr);
        
        if (int_Same > 0) {
            double scale = int_Mix / int_Same;
            hNorm->Scale(scale);
        }

        hNorm_Uncorr = (TH1*)hNorm->Clone(fName + "_Uncorr");
        hNorm_Corr   = (TH1*)hNorm->Clone(fName + "_Corr");
        hNorm_Uncorr->SetDirectory(nullptr);
        hNorm_Corr->SetDirectory(nullptr);

        SetStyle();
    }

    TH1* CorrectCF(TH1* hNorm_SB, TH1* hNorm_BS, TH1* hNorm_BB){
        if (!hNorm || !hNorm_SB || !hNorm_BB || !hNorm_Corr) return nullptr;

        double P_avg = (mP_L + mP_A) / 2.0;
        double P1 = mP_L; double P2 = mP_L; 

        if (fName.Contains("LL")) { P1 = mP_L; P2 = mP_L; } 
        else if (fName.Contains("AA")) { P1 = mP_A; P2 = mP_A; } 
        else if (fName.Contains("LA")) { P1 = mP_L; P2 = mP_A; } 
        else if (fName.Contains("LL+AA")) { P1 = P_avg; P2 = P_avg; }

        double f_SS = P1 * P2;
        double f_BB = (1.0 - P1) * (1.0 - P2);
        
        for (int i = 1; i <= hNorm->GetNbinsX(); i++) {
            double c_meas = hNorm->GetBinContent(i);  
            double err_meas = hNorm->GetBinError(i);   

            double c_SB   = hNorm_SB->GetBinContent(i);
            double err_SB = hNorm_SB->GetBinError(i);
            double c_BS   = hNorm_BS->GetBinContent(i);
            double err_BS = hNorm_BS->GetBinError(i);
            double c_BB   = hNorm_BB->GetBinContent(i);
            double err_BB = hNorm_BB->GetBinError(i);
            
            double c_SS = (c_meas - (1.0-P2)*c_SB - (1.0-P1)*c_BS + f_BB*c_BB) / f_SS;
            
            double err_SS = std::sqrt(
                (err_meas * err_meas) + 
                ((1.0-P1) * (1.0-P1) * err_SB * err_SB) + 
                ((1.0-P2) * (1.0-P2) * err_BS * err_BS) +
                (f_BB * f_BB * err_BB * err_BB)
            ) / f_SS;

            hNorm_Corr->SetBinContent(i, c_SS);
            hNorm_Corr->SetBinError(i, err_SS);
        }
        return hNorm_Corr;
    }

    TH1* CorrectCF_Flat(){
        if (!hNorm || !hNorm_Corr) return nullptr;

        double P_avg = (mP_L + mP_A) / 2.0;
        double P1 = mP_L; double P2 = mP_L; 

        if (fName.Contains("LL")) { P1 = mP_L; P2 = mP_L; } 
        else if (fName.Contains("AA")) { P1 = mP_A; P2 = mP_A; } 
        else if (fName.Contains("LA")) { P1 = mP_L; P2 = mP_A; } 
        else if (fName.Contains("LL+AA")) { P1 = P_avg; P2 = P_avg; }

        double f_SS = P1 * P2;

        for (int i = 1; i <= hNorm->GetNbinsX(); i++) {
            double c_meas = hNorm->GetBinContent(i);
            double err_meas = hNorm->GetBinError(i);

            if (c_meas == 0 && err_meas == 0) continue;

            double c_SS = 1.0 + (c_meas - 1.0) / f_SS;
            double err_SS = err_meas / f_SS;

            hNorm_Corr->SetBinContent(i, c_SS);
            hNorm_Corr->SetBinError(i, err_SS);
        }
        return hNorm_Corr;
    }

    
    // 計算用メソッド
    void Normalized(double normMin, double normMax);
    TH1* CorrectCF(TH1* hNorm_SB, TH1* hNorm_BB);

    // ゲッター（結果の取得）
    TH1* GetRaw()        const { return hRaw; }
    TH1* GetNorm()       const { return hNorm; }
    TH1* GetNormUncorr() const { return hNorm_Uncorr; }
    TH1* GetNormCorr()   const { return hNorm_Corr; }

    void DrawRaw(TH1* h, TString opt = "PE") {
        if(!h) return;
        h->SetXTitle("Q (GeV/c)");
        h->SetYTitle("C(Q)");
        h->SetMinimum(0.2);
        h->SetMaximum(2.4);
        h->SetMarkerStyle(20);
        h->SetMarkerColor(kBlack);
        h->SetLineColor(kBlack);
        h->SetMarkerSize(0.5);
        h->Draw(opt);
    }

    void DrawNorm(TH1* h, TString opt = "PE") {
        if(!h) return;
        h->SetXTitle("Q (GeV/c)");
        h->SetYTitle("C(Q)");
        h->SetMinimum(0.4);
        h->SetMaximum(1.4);
        h->SetMarkerStyle(20);
        h->SetMarkerColor(kRed);
        h->SetLineColor(kRed);
        h->SetMarkerSize(0.5);
        SetAxisCF(h);
        h->Draw(opt);

        int firstBin = h->GetXaxis()->GetFirst();
        int lastBin  = h->GetXaxis()->GetLast();
        double xmin = h->GetXaxis()->GetBinLowEdge(firstBin);
        double xmax = h->GetXaxis()->GetBinUpEdge(lastBin);
        TLine *line = new TLine(0.0, 1.0, 1.5, 1.0);
        line->SetLineStyle(2); 
        line->SetLineColor(kBlack);
        line->SetLineWidth(2);
        line->Draw("same");
    }
};

inline TH1* MergeHist(TH1* hLL, TH1* hAA) {
    if (!hLL || !hAA ) return nullptr;
    TH1 *hMerged = (TH1*)hLL->Clone(hLL->GetName() + TString("_Merged"));
    hMerged->SetDirectory(nullptr);
    if (!hMerged->GetSumw2N()) hMerged->Sumw2();
    hMerged->Add(hAA, 1.0);
    return hMerged;
}

#endif