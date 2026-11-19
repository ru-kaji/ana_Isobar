#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "AxisStyle.h"
#include <iostream>

// 1. 2Dヒストグラム (qinv vs mass) から重み付き1D mass分布への射影
TH1D* ProjectWeightedMass(TH2D* h2, TH1D* hWeight, const char* name) {
    if (!h2) return nullptr;
    TH1D* hProj = h2->ProjectionY(name, 1, 1);
    hProj->Reset();
    hProj->Sumw2();

    for (int ix = 1; ix <= h2->GetNbinsX(); ++ix) {
        double qinv = h2->GetXaxis()->GetBinCenter(ix);
        double w = 1.0;
        if (hWeight && qinv > 0.0 && qinv < 2.0) {
            w = hWeight->GetBinContent(hWeight->FindBin(qinv));
        }
        TH1D* hSlice = h2->ProjectionY(Form("%s_slice_%d", name, ix), ix, ix);
        hProj->Add(hSlice, w);
        delete hSlice;
    }
    return hProj;
}

// 2. 指定領域でのスケーリング
TH1D* ScaleHist(TH1D* hSame, TH1D* hMix, double xmin, double xmax, const char* name) {
    int b1 = hSame->FindBin(xmin);
    int b2 = hSame->FindBin(xmax);
    TH1D* hScaled = (TH1D*)hMix->Clone(name);
    
    double intMix = hMix->Integral(b1, b2);
    if (intMix > 0) {
        hScaled->Scale(hSame->Integral(b1, b2) / intMix);
    }
    return hScaled;
}

// 3. メイン関数（引数は OutAnaTree_isobar_◯◯.root の1つだけ！）
void DrawWeightedMass(TString fileName = "OutAnaTree_isobar_phi5eta4_cuttracks_R1p4.root") {
    gStyle->SetOptStat(0);
    gStyle->SetPadBottomMargin(0.14);
    gStyle->SetPadLeftMargin(0.15);

    // ① ファイル名から自動で Weight_◯◯.root の名前を生成
    TString weightFile = fileName;
    weightFile.ReplaceAll("OutAnaTree_isobar_", "Weight_");

    std::cout << "📂 Input File  : " << fileName << std::endl;
    std::cout << "⚖️ Weight File : " << weightFile << std::endl;

    // ② データファイルと重みファイルの読み込み
    TFile *fin = TFile::Open(fileName, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "❌ Error: Cannot open " << fileName << std::endl;
        return;
    }

    TFile *fw = TFile::Open(weightFile, "READ");
    TH1D *hWeight = (fw && !fw->IsZombie()) ? (TH1D*)fw->Get("hWeight") : nullptr;
    if (!hWeight) {
        std::cerr << "⚠️ Warning: " << weightFile << " not found! Using weight = 1.0" << std::endl;
    }

    // ③ LL と AA のヒストグラム取得 (icut=0, SS, 90-180[Anti])
    TH1D *hSame_LL = (TH1D*)fin->Get("hmass_same_0_LL_SS_90-180");
    TH1D *hSame_AA = (TH1D*)fin->Get("hmass_same_0_AA_SS_90-180");
    TH1D *hMix_LL  = (TH1D*)fin->Get("hmass_mix_0_LL_SS_90-180");
    TH1D *hMix_AA  = (TH1D*)fin->Get("hmass_mix_0_AA_SS_90-180");
    TH2D *h2_LL    = (TH2D*)fin->Get("h2_q_m_mix_LL");
    TH2D *h2_AA    = (TH2D*)fin->Get("h2_q_m_mix_AA");

    if (!hSame_LL || !hSame_AA || !hMix_LL || !hMix_AA || !h2_LL || !h2_AA) {
        std::cerr << "❌ Error: Missing required histograms in " << fileName << std::endl;
        return;
    }

    // ④ 2Dヒストグラムを合算して重み付き1D mass分布へ射影
    TH2D *h2_LLAA = (TH2D*)h2_LL->Clone("h2_q_m_mix_LLAA");
    h2_LLAA->Add(h2_AA);
    TH1D *hMixW_LLAA = ProjectWeightedMass(h2_LLAA, hWeight, "hMixW_LLAA");

    // ⑤ 1Dヒストグラム (Same と 重みなしMix) の LL+AA 合算
    TH1D *hSame_LLAA = (TH1D*)hSame_LL->Clone("hSame_LLAA");
    if (hSame_LLAA->GetSumw2N() == 0) hSame_LLAA->Sumw2();
    hSame_LLAA->Add(hSame_AA);

    TH1D *hMixNoW_LLAA = (TH1D*)hMix_LL->Clone("hMixNoW_LLAA");
    if (hMixNoW_LLAA->GetSumw2N() == 0) hMixNoW_LLAA->Sumw2();
    hMixNoW_LLAA->Add(hMix_AA);

    // ⑥ スケーリングと差分 (Same - Mix) の計算
    TH1D *hMixW_scaled   = ScaleHist(hSame_LLAA, hMixW_LLAA,   2.24, 2.2799, "hMixW_scaled");
    TH1D *hMixNoW_scaled = ScaleHist(hSame_LLAA, hMixNoW_LLAA, 2.24, 2.2799, "hMixNoW_scaled");

    TH1D *hDiff_W = (TH1D*)hSame_LLAA->Clone("hDiff_W");
    hDiff_W->Add(hMixW_scaled, -1.0);

    TH1D *hDiff_NoW = (TH1D*)hSame_LLAA->Clone("hDiff_NoW");
    hDiff_NoW->Add(hMixNoW_scaled, -1.0);

    // ⑦ 描画
    TCanvas *c1 = new TCanvas("c_weigh_diff", "Same - Mix (LL+AA)", 700, 600);
    
    hDiff_W->SetTitle("#Lambda#Lambda + #bar{#Lambda}#bar{#Lambda} (SS, Anti-Parallel);M_{inv} (GeV/c^{2});Same - Mix Counts");
    hDiff_W->SetMarkerStyle(20);
    hDiff_W->SetMarkerSize(0.9);
    hDiff_W->SetMarkerColor(kRed);
    hDiff_W->SetLineColor(kRed);
    hDiff_W->GetYaxis()->SetRangeUser(-800, 800);
    SetAxis(hDiff_W);
    hDiff_W->Draw("PE");

    hDiff_NoW->SetMarkerStyle(24);
    hDiff_NoW->SetMarkerSize(0.9);
    hDiff_NoW->SetMarkerColor(kBlue);
    hDiff_NoW->SetLineColor(kBlue);
    hDiff_NoW->Draw("PE same");

    TLine *line = new TLine(hDiff_W->GetXaxis()->GetXmin(), 0.0, hDiff_W->GetXaxis()->GetXmax(), 0.0);
    line->SetLineStyle(2);
    line->SetLineColor(kBlack);
    line->SetLineWidth(2);
    line->Draw("same");

    TLegend *leg = new TLegend(0.45, 0.72, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(43);
    leg->SetTextSize(18);
    leg->AddEntry(hDiff_W,   "Same - Mix_{weight (no BW)}", "pe");
    leg->AddEntry(hDiff_NoW, "Same - Mix_{no weight}", "pe");
    leg->Draw();
}