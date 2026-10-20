#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"
#include "TLatex.h"
#include "AxisStyle.h"
#include <iostream>
#include <cmath>

//void SetAxis(TH1* h);

class CFAnalyzer {

private:
    TH1D *hSame, *hMix;
    TH1D *hRaw, *hNorm;
    TH1D *hNorm_SB, *hNorm_BS, *hNorm_BB;
    TH1D *hNorm_Uncorr;    // Purity補正前
    TH1D *hNorm_Corr;
    TString fName;         // 管理用の名前
    double mP_L, mP_A;
    
public:
    CFAnalyzer(TH1D* same, TH1D* mix, TString name, double pL = 0.9833, double pA = 0.9843) {//初期化
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
            hSame = (TH1D*)same->Clone(fName + "_same");
            hMix  = (TH1D*)mix->Clone(fName + "_mix");

            /*int rebin_factor = 2; 
            hSame->Rebin(rebin_factor);
            hMix->Rebin(rebin_factor);*/

            if (!hSame->GetSumw2N()) hSame->Sumw2();
            if (!hMix->GetSumw2N()) hMix->Sumw2();
        } else {
            //std::cout << "Error: Input histograms are null!" << std::endl;
        }
    }

    // --- 2. デストラクタ (後片付け) ---
    // クラスが不要になったらメモリを掃除する
    virtual ~CFAnalyzer() {
        // 必要に応じて delete するが、ROOTの場合、TFileを閉じると
        // 自動で消えたりするので、二重deleteに注意が必要。
        // ここでは安全のため明示的な delete は省略します。
    }


    // --- Caliculate CF ---
    void Calculate(double normMin, double normMax) {
        if(!hSame || !hMix) return;

        // Raw
        hRaw = (TH1D*)hSame->Clone(fName + "_Raw");
        hRaw->Divide(hMix);

        // Norm
        int bin_min = hSame->FindBin(normMin);
        int bin_max = hSame->FindBin(normMax);

        // Normalize
        double int_Same = hSame->Integral(bin_min, bin_max);
        double int_Mix  = hMix->Integral(bin_min, bin_max);

        hNorm = (TH1D*)hRaw->Clone(fName + "_Norm");
        
        if (int_Same > 0) {
            double scale = int_Mix / int_Same;
            hNorm->Scale(scale);
        }

        hNorm_Uncorr = (TH1D*)hNorm->Clone(fName + "_Uncorr");
        hNorm_Corr = (TH1D*)hNorm->Clone(fName + "_Corr");

        SetStyle();
    }


    // -------------------- Purity Correction --------------------
    // C_corr(Q) = (C_meas(Q) - 1) / P_pair + 1
    // Error_corr = Error_meas / P_pair
    TH1D* CorrectCF(TH1D* hNorm_SB, TH1D* hNorm_BS, TH1D* hNorm_BB){

        if (!hNorm || !hNorm_SB || !hNorm_BB || !hNorm_Corr) {
        //std::cout << "Error: Null pointers in CorrectCF!" << std::endl;
        return nullptr;
        }

        double P_avg = (mP_L + mP_A) / 2.0;
        double P1 = mP_L;
        double P2 = mP_L; // デフォルト値

        if (fName == "LL") {
            P1 = mP_L; P2 = mP_L;
        } else if (fName == "AA") {
            P1 = mP_A; P2 = mP_A;
        } else if (fName == "LA") {
            P1 = mP_L; P2 = mP_A;
        } else if (fName == "LL+AA") {
            P1 = P_avg; P2 = P_avg;
        }

        double f_SS = P1 * P2;
        //double f_SB = P1 * (1.0 - P2) + (1.0 - P1) * P2;
        double f_BB = (1.0 - P1) * (1.0 - P2);
        
        for (int i = 1; i <= hNorm->GetNbinsX(); i++) {
            double c_meas = hNorm->GetBinContent(i);  //縦軸の値を取得
            double err_meas = hNorm->GetBinError(i);   //縦軸の誤差を取得

            double c_SB   = hNorm_SB->GetBinContent(i);
            double err_SB = hNorm_SB->GetBinError(i);

            double c_BS   = hNorm_BS->GetBinContent(i);
            double err_BS = hNorm_BS->GetBinError(i);

            double c_BB   = hNorm_BB->GetBinContent(i);
            double err_BB = hNorm_BB->GetBinError(i);
            
            double c_SS = (c_meas - (1.0-P2)*c_SB-(1.0-P1)*c_BS + f_BB*c_BB) /f_SS;
            
            double err_SS = std::sqrt(
                (err_meas * err_meas) + 
                ((1.0-P1) * (1.0-P1) * err_SB * err_SB) + 
                ((1.0-P2) * (1.0-P2) * err_BS * err_BS) +
                (f_BB * f_BB * err_BB * err_BB)
            ) / f_SS;

            hNorm_Corr->SetBinContent(i, c_SS);
            hNorm_Corr->SetBinError(i, err_SS);
        }

        SetStyle(); // View
        return hNorm_Corr;
    }


    // --- Draw CF (hRaw, hNorm)---
    void DrawRaw(TH1D* hRaw, TString opt = "PE") {
        if(!hRaw) return;

        hRaw->SetXTitle("Q_{inv} (GeV/c)");
        hRaw->SetYTitle("C(Q_{inv})");
        hRaw->SetMinimum(-0.1);
        hRaw->SetMaximum(1.6);
        hRaw->SetMarkerStyle(20);
        hRaw->SetMarkerColor(kBlack);
        hRaw->SetLineColor(kBlack);
        hRaw->SetMarkerSize(0.5);
        SetAxisCF(hRaw);
        hRaw->Draw(opt);
    }
       
    void DrawNorm(TH1D* hNorm, TString opt) {
        if(!hNorm) return;

        hNorm->SetXTitle("Q_{inv} (GeV/c)");
        hNorm->SetYTitle("C(Q_{inv})");
        hNorm->SetMarkerStyle(20);
        hNorm->SetMarkerColor(kRed);
        hNorm->SetLineColor(kRed);
        hNorm->SetMarkerSize(0.5);
        SetAxisCF(hNorm);
        hNorm->Draw(opt);

        // 基準線
        double xmin = hNorm->GetXaxis()->GetXmin();
        double xmax = hNorm->GetXaxis()->GetXmax();
        TLine *line = new TLine(xmin, 1.0, xmax, 1.0);
        line->SetLineStyle(2); // dot line
        line->SetLineColor(kBlack);
        line->SetLineWidth(2);
        line->Draw("same");
    }

    // ゲッター (結果のヒストグラムだけ欲しい時に使う)
    TH1D* GetRaw()  { return hRaw; }
    TH1D* GetNorm() { return hNorm; }
    TH1D* GetNormUncorr() { return hNorm_Uncorr; }
    TH1D* GetNormCorr() { return hNorm_Corr; }
    

private:
    void SetStyle() {  //※要修正
        if(hNorm) {  
            hNorm->SetXTitle("Q_{inv} (GeV/c)");
            hNorm->SetYTitle("C(Q_{inv})");
            hNorm->SetMarkerStyle(20);
            hNorm->SetMarkerColor(kBlack);
            hNorm->SetLineColor(kRed);
            hNorm->SetMarkerSize(0.5);
            hNorm->SetMinimum(0.0);
            hNorm->SetMaximum(2.0);
        }
    }
}; // end of class CFAnalyzer


void DrawSpin(TH1D* hPara, TH1D* hAnti, TString title){

    if( !hPara || !hAnti){
        //std::cout << "Error: Input histograms for DrawSpin are null!" << std::endl;
        return;
    }


    TH1D *hNorm_Para = (TH1D*)hPara->Clone(title + "_Para");
    TH1D *hNorm_Anti = (TH1D*)hAnti->Clone(title + "_Anti");

        // Parallel
        hNorm_Para->SetTitle(title);
        hNorm_Para->GetXaxis()->SetTitle("Q_{inv} (GeV/c)");
        hNorm_Para->GetYaxis()->SetTitle("C(Q_{inv})");
        hNorm_Para->SetMarkerStyle(20); // 丸
        hNorm_Para->SetMarkerColor(kRed);
        hNorm_Para->SetLineColor(kRed);
        hNorm_Para->SetMarkerSize(0.6);
        hNorm_Para->SetMinimum(0.4); 
        hNorm_Para->SetMaximum(1.4);
        //hNorm_Para->GetXaxis()->SetRangeUser(0.0, 2.0); 
        hNorm_Para->GetXaxis()->SetRangeUser(0.0, 2.0); 

        // Anti-Parallel
        hNorm_Anti->SetMarkerStyle(21); // 四角
        hNorm_Anti->SetMarkerColor(kBlue);
        hNorm_Anti->SetLineColor(kBlue);
        hNorm_Anti->SetMarkerSize(0.6);
        //hNorm_Anti->GetXaxis()->SetRangeUser(0.0, 2.0); 
        hNorm_Anti->GetXaxis()->SetRangeUser(0.0, 2.0); 

        // Draw same
        hNorm_Para->Draw("PE"); 
        hNorm_Anti->Draw("PE same"); 
        SetAxisCF(hNorm_Para);
        SetAxisCF(hNorm_Anti);
        hNorm_Para->GetYaxis()->SetTitleOffset(1.4);

        SetAxisCF(hNorm_Para);

        // 凡例
        TLegend *leg = new TLegend(0.6, 0.22, 0.8, 0.37);
        leg->SetBorderSize(0);
        leg->SetBorderSize(0); // 枠線をなくす
        leg->SetFillStyle(0);  // 背景を透明にする
        leg->SetTextSize(0.055);
        leg->AddEntry(hNorm_Para, "cos#theta_{12} > 0", "lp"); 
        leg->AddEntry(hNorm_Anti, "cos#theta_{12} < 0", "lp");
        leg->Draw();

        // 基準線
        //TLine *line = new TLine(0.0, 1.0, 2.0, 1.0);
        TLine *line = new TLine(0.0, 1.0, 0.5, 1.0);
        line->SetLineStyle(2);
        line->SetLineColor(kBlack);
        line->Draw();

}

TH1D* MergeHist(TH1D* hLL, TH1D* hAA) {
    if (!hLL || !hAA ) return nullptr;

    TH1D *hMerged = (TH1D*)hLL->Clone(hLL->GetName() + TString("_Merged"));
    hMerged->SetDirectory(nullptr);
    if (!hMerged->GetSumw2N()) hMerged->Sumw2();

    hMerged->Add(hAA, 1.0);

    return hMerged;
}


// =========================================================
//  メイン関数
// =========================================================
void plotCF_isobar(TString fileName) { 

    //std::cout << "DEBUG: plotCF.C is loading file: " << fileName << std::endl;

    TFile *fin = new TFile(fileName, "READ");
    if (!fin || fin->IsZombie()) {
        std::cout << "!!! Error: Cannot open file .root !!!" << std::endl;
        return;
    }

    double pL[3] = {0.9833, 0.9840, 0.9806}; 
    double pA[3] = {0.9843, 0.9822, 0.9788};
    
    TH1D *hSame[5][4][4][7][9] = {nullptr};
    TH1D *hMix[5][4][4][7][9] = {nullptr};
    TH1D *hSame_nocut[4][4][7][9] = {nullptr};
    TH1D *hMix_nocut[4][4][7][9] = {nullptr};

    const char* pidName[4]   = {"LL", "AA", "LA", "LL+AA"};
    const char* pairName[4]  = {"SS", "SB", "BS", "BB"};
    const char* spinName[3]  = {"All", "Para", "Anti"};
    const char* spinName7[7] = {"All", "0-90", "90-180", "0-45", "45-90", "90-135", "135-180"};
    const char* centName[9]  = {"0-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "All"};
    const char* cent3Name[3] = {"0-20%", "20-50%", "50-80%"};
    
    for(int ipid = 0; ipid < 3; ipid++) {
        for(int ipair = 0; ipair < 4; ipair++) {
            for(int ispin = 0; ispin < 7; ispin++) {
                for(int icent = 0; icent < 9; icent++) {

                    for(int icut=0; icut<5; icut++){
                        TString nameQinvS = Form("hqinv_same_%d_%s_%s_%s_%s", icut, pidName[ipid], pairName[ipair], spinName7[ispin], centName[icent]);
                        TString nameQinvM = Form("hqinv_mix_%d_%s_%s_%s_%s", icut, pidName[ipid], pairName[ipair], spinName7[ispin], centName[icent]);

                        hSame[icut][ipid][ipair][ispin][icent] = (TH1D*)fin->Get(nameQinvS);
                        hMix[icut][ipid][ipair][ispin][icent]  = (TH1D*)fin->Get(nameQinvM);

                        if (!hMix[icut][ipid][ipair][ispin][icent]) {
                            std::cout << "⚠️  WARNING: " << nameQinvM << " がファイルに見つかりません！" << std::endl;
                        }
                    }
                    TString nameQinvS_nocut = Form("hqinv_same_%s_%s_%s_%s", pidName[ipid], pairName[ipair], spinName7[ispin], centName[icent]);
                    TString nameQinvM_nocut = Form("hqinv_mix_%s_%s_%s_%s", pidName[ipid], pairName[ipair], spinName7[ispin], centName[icent]);
                    hSame_nocut[ipid][ipair][ispin][icent] = (TH1D*)fin->Get(nameQinvS_nocut+"nocut");
                    hMix_nocut[ipid][ipair][ispin][icent]  = (TH1D*)fin->Get(nameQinvM_nocut+"nocut");

                    if (!hSame_nocut[ipid][ipair][ispin][icent]) {
                        std::cout << "⚠️  WARNING: " << nameQinvS_nocut+"nocut" << " がファイルに見つかりません！" << std::endl;
                    }
                }
            }
        }
    }

    
    for(int ipair = 0; ipair < 4; ipair++) {
        for(int ispin = 0; ispin < 7; ispin++) {
            for(int icent = 0; icent < 9; icent++){
                for(int icut = 0; icut < 5; icut++) {
                    if (!hSame[icut][0][ipair][ispin][icent] || !hSame[icut][1][ipair][ispin][icent]) {
                        std::cout << "⚠️  WARNING: MergeHist ipid=3 失敗 - 元のヒストグラムが存在しない ipair=" << ipair << " ispin=" << ispin << " icent=" << icent << std::endl;
                        hSame[icut][3][ipair][ispin][icent] = nullptr;
                        hMix[icut][3][ipair][ispin][icent] = nullptr;
                    } else {
                        hSame[icut][3][ipair][ispin][icent] = (TH1D*)MergeHist(hSame[icut][0][ipair][ispin][icent], hSame[icut][1][ipair][ispin][icent]);
                        hMix[icut][3][ipair][ispin][icent]  = (TH1D*)MergeHist(hMix[icut][0][ipair][ispin][icent], hMix[icut][1][ipair][ispin][icent]);
                    }
                }
                    hSame_nocut[3][ipair][ispin][icent] = (TH1D*)MergeHist(hSame_nocut[0][ipair][ispin][icent], hSame_nocut[1][ipair][ispin][icent]);
                    hMix_nocut[3][ipair][ispin][icent]  = (TH1D*)MergeHist(hMix_nocut[0][ipair][ispin][icent], hMix_nocut[1][ipair][ispin][icent]);
            }
        }
    }
    

    TH1D *hSame_0to20per[4][4][7]  = {nullptr};   TH1D *hMix_0to20per[4][4][7]  = {nullptr}; // [ipid][ipair][ispin]
    TH1D *hSame_20to40per[4][4][7] = {nullptr};   TH1D *hMix_20to40per[4][4][7] = {nullptr};
    TH1D *hSame_20to50per[4][4][7] = {nullptr};   TH1D *hMix_20to50per[4][4][7] = {nullptr};
    TH1D *hSame_50to70per[4][4][7] = {nullptr};   TH1D *hMix_50to70per[4][4][7] = {nullptr};
    TH1D *hSame_50to80per[4][4][7] = {nullptr};   TH1D *hMix_50to80per[4][4][7] = {nullptr};
    TH1D *hSame_20to80per[4][4][7] = {nullptr};   TH1D *hMix_20to80per[4][4][7] = {nullptr};
    TH1D *hSame_10to80per[4][4][7] = {nullptr};   TH1D *hMix_10to80per[4][4][7] = {nullptr};

    for(int ipid = 0; ipid < 4; ipid++) {
        for(int ispin = 0; ispin < 7; ispin++){
            for(int ipair = 0; ipair < 4; ipair++){
                hSame_0to20per[ipid][ipair][ispin]  =  (TH1D*)MergeHist(hSame[3][ipid][ipair][ispin][0], hSame[3][ipid][ipair][ispin][1]);
                hMix_0to20per[ipid][ipair][ispin]   =  (TH1D*)MergeHist(hMix[3][ipid][ipair][ispin][0], hMix[3][ipid][ipair][ispin][1]);

                hSame_20to40per[ipid][ipair][ispin]    = (TH1D*)MergeHist(hSame[3][ipid][ipair][ispin][2], hSame[3][ipid][ipair][ispin][3]);
                hMix_20to40per[ipid][ipair][ispin]     = (TH1D*)MergeHist(hMix[3][ipid][ipair][ispin][2], hMix[3][ipid][ipair][ispin][3]);
                hSame_20to50per[ipid][ipair][ispin] = (TH1D*)MergeHist(hSame[3][ipid][ipair][ispin][4], hSame_20to40per[ipid][ipair][ispin]);
                hMix_20to50per[ipid][ipair][ispin]  = (TH1D*)MergeHist(hMix[3][ipid][ipair][ispin][4], hMix_20to40per[ipid][ipair][ispin]);

                hSame_50to70per[ipid][ipair][ispin] = (TH1D*)MergeHist(hSame[3][ipid][ipair][ispin][5], hSame[3][ipid][ipair][ispin][6]);
                hMix_50to70per[ipid][ipair][ispin]  = (TH1D*)MergeHist(hMix[3][ipid][ipair][ispin][5], hMix[3][ipid][ipair][ispin][6]);
                hSame_50to80per[ipid][ipair][ispin] = (TH1D*)MergeHist(hSame[3][ipid][ipair][ispin][7], hSame_50to70per[ipid][ipair][ispin]);
                hMix_50to80per[ipid][ipair][ispin]  = (TH1D*)MergeHist(hMix[3][ipid][ipair][ispin][7], hMix_50to70per[ipid][ipair][ispin]);
                
                hSame_20to80per[ipid][ipair][ispin] = (TH1D*)MergeHist(hSame_50to80per[ipid][ipair][ispin], hSame_20to50per[ipid][ipair][ispin]);
                hMix_20to80per[ipid][ipair][ispin]  = (TH1D*)MergeHist(hMix_50to80per[ipid][ipair][ispin],  hMix_20to50per[ipid][ipair][ispin]);
                hSame_10to80per[ipid][ipair][ispin] = (TH1D*)MergeHist(hSame[3][ipid][ipair][ispin][1], hSame_20to80per[ipid][ipair][ispin]);
                hMix_10to80per[ipid][ipair][ispin]  = (TH1D*)MergeHist(hMix[3][ipid][ipair][ispin][1],  hMix_20to80per[ipid][ipair][ispin]);
            }
        }
    }

    TH1D *hPara_LL = nullptr;
    TH1D *hPara_AA = nullptr;
    TH1D *hAnti_LL = nullptr;
    TH1D *hAnti_AA = nullptr;

    TH1D *hSave_CF_All[4]  = {nullptr};
    TH1D *hSave_CF_Para[4] = {nullptr};
    TH1D *hSave_CF_Anti[4] = {nullptr};

    TH1D *hSave_Same_All[4] = {nullptr}; TH1D *hSave_Mix_All[4] = {nullptr};
    TH1D *hSave_Same_Para[4]= {nullptr}; TH1D *hSave_Mix_Para[4]= {nullptr};
    TH1D *hSave_Same_Anti[4]= {nullptr}; TH1D *hSave_Mix_Anti[4]= {nullptr};
    
    //--------------------------------------------------
    //  [pid] = [0: LL, 1: AA, 2: LA, 3:LL+AA]
    //
	//  c_cf1[pid]  : (Raw, Norm) × (SS, SB, BB) 
	//  c_cf2[pid]  : (Raw, Norm) × (All, Para, Anti)
	//  c_cf3[pid]  : Normed (with TLegend) 
	//  c_cf4[pid]  : Overlay Spin (SS) (All range)
	//  c_cf5[pid]  : Overlay Spin (SS, SB, BB) (Low range)
	//  c_cf6[pid]  : Overlay Correction (Para, Anti)
    //  c_cf7[pid]  : Overlay Spin (centrality by 10% step)
    //  c_cf8[pid]  : Overlay Spin (centrality 0-20, 20-50. 50-80%)
	//--------------------------------------------------

    TCanvas *c_cf1, *c_cf2, *c_cf3, *c_cf4, *c_cf5, *c_cf6, *c_cf6_2, *c_cf7[4], *c_cf8, *c_cf9, *c_cf10;
    

        if (!hSame[3][0][0][0][0] || !hSame[3][1][1][0][1] || !hSame[3][2][0][1][8]) {
            // << "Warning: hSame is empty or null!" << std::endl;
        }

        // ---------------------------------------------------------
        // c_cf1 : (Raw, Norm) × (SS, SB, BS, BB)
        //----------------------------------------------------------
        c_cf1 = new TCanvas("c_cf1", "CF1", 1000, 800);
        c_cf1->Divide(4, 4);

        for(int ipair = 0; ipair < 4; ipair++) {
            for(int ipid = 0; ipid < 4; ipid++) {
                
                int ispin = 0;

                CFAnalyzer ana1(hSame[3][ipid][ipair][ispin][8], hMix[3][ipid][ipair][ispin][8], pidName[ipid]);
                ana1.Calculate(0.5, 0.8);

                c_cf1->cd(ipid+ipair*4+1);
                TH1D *hRaw  = ana1.GetRaw();
                TH1D *hNorm = ana1.GetNorm();

                hRaw->SetTitle(Form("%s - %s", pidName[ipid], pairName[ipair]));
                ana1.DrawRaw(hRaw, "PE");
                ana1.DrawNorm(hNorm, "PE same");
            }
        }
      
        // --------------------------------------------------------
        // c_cf2 : (Raw, Norm) × (All, Para, Anti)
        //---------------------------------------------------------
        c_cf2 = new TCanvas("c_cf2", "CF2", 1100, 800);
        c_cf2->Divide(4, 3);

        for(int ispin = 0; ispin < 4; ispin++) {
            for(int ipid = 0; ipid < 4; ipid++) {

                int ipair = 0;

                CFAnalyzer ana2(hSame[3][ipid][ipair][ispin][8], hMix[3][ipid][ipair][ispin][8], pidName[ipid]);
                ana2.Calculate(0.5, 0.8);

                c_cf2->cd(ipid+ispin*4+1);
                TH1D *hRaw  = ana2.GetRaw();
                TH1D *hNorm = ana2.GetNorm();

                hRaw->SetTitle(Form("%s - %s", pidName[ipid], spinName[ispin]));
                ana2.DrawRaw(hRaw, "PE");
                ana2.DrawNorm(hNorm, "PE same");
            }
        }

        // ---------------------------------------------------------
        // c_cf3 : Normed + Correcter (All Spin)
        //----------------------------------------------------------
        c_cf3 = new TCanvas("c_cf3", "CF3",1400, 400);
        c_cf3 ->Divide(4, 1);

        for(int ipid = 0; ipid < 4; ipid++) {
            c_cf3 -> cd(ipid+1);
            CFAnalyzer ana3_SS(hSame[3][ipid][0][0][8], hMix[3][ipid][0][0][8], pidName[ipid], pL[0], pA[0]);
            CFAnalyzer ana3_SB(hSame[3][ipid][1][0][8], hMix[3][ipid][1][0][8], pidName[ipid], pL[0], pA[0]);
            CFAnalyzer ana3_BS(hSame[3][ipid][2][0][8], hMix[3][ipid][2][0][8], pidName[ipid], pL[0], pA[0]);
            CFAnalyzer ana3_BB(hSame[3][ipid][3][0][8], hMix[3][ipid][3][0][8], pidName[ipid], pL[0], pA[0]);

            // 規格化
            ana3_SS.Calculate(0.5, 0.8);
            ana3_SB.Calculate(0.5, 0.8);
            ana3_BS.Calculate(0.5, 0.8);
            ana3_BB.Calculate(0.5, 0.8);

            ana3_SS.CorrectCF(ana3_SB.GetNorm(), ana3_BS.GetNorm(), ana3_BB.GetNorm());
            TH1D *hNormCorr_All = ana3_SS.GetNormCorr();
            if(!hNormCorr_All) continue;

            hSave_CF_All[ipid]   = (TH1D*)hNormCorr_All->Clone(Form("hCF_All_%s", pidName[ipid]));
            hSave_Same_All[ipid] = (TH1D*)hSame[3][ipid][0][0][8]->Clone(Form("hSame_All_%s", pidName[ipid]));
            hSave_Mix_All[ipid]  = (TH1D*)hMix[3][ipid][0][0][8]->Clone(Form("hMix_All_%s", pidName[ipid]));
            
            // 描画設定
            hNormCorr_All->SetTitle(Form("%s;Q_{inv} (GeV/c);C(Q_{inv})", pidName[ipid]));
            hNormCorr_All->GetXaxis()->SetRangeUser(0.0, 0.6); // 低Q領域にフォーカス
            hNormCorr_All->SetMinimum(0.6);
            hNormCorr_All->SetMaximum(1.4);
            
            hNormCorr_All->SetMarkerStyle(20);
            hNormCorr_All->SetMarkerColor(kBlack);
            hNormCorr_All->SetLineColor(kBlack);
            
            hNormCorr_All->Draw("PE");
            SetAxisCF(hNormCorr_All);
            
            // 基準線
            TLine *line_cf3 = new TLine(0.0, 1.0, 0.6, 1.0);
            line_cf3->SetLineStyle(2);
            line_cf3->Draw();
        }


        // ---------------------------------------------------------
        // c_cf4 : Overlay Spin (SS) (All range) (Corrected)
        //----------------------------------------------------------
        c_cf4 = new TCanvas("c_cf4", "CF4", 1400, 350);
        c_cf4->Divide(4, 1);

        for(int ipid = 0; ipid < 4; ipid++) {

            if (!hSame[3][ipid][0][1][8] || !hMix[3][ipid][0][1][8]) {
                std::cout << "❌ ERROR: c_cf4 - hSame[" << ipid << "][0][1][8] または hMix が nullptr です！" << std::endl;
                std::cout << "   hSame[" << ipid << "][0][1][8] = " << hSame[3][ipid][0][1][8] << std::endl;
                std::cout << "   hMix[" << ipid << "][0][1][8] = " << hMix[3][ipid][0][1][8] << std::endl;
                continue;
            }

            CFAnalyzer ana4_Para_SS(hSame[3][ipid][0][1][8], hMix[3][ipid][0][1][8], pidName[ipid]);
            CFAnalyzer ana4_Para_SB(hSame[3][ipid][1][1][8], hMix[3][ipid][1][1][8], pidName[ipid]);
            CFAnalyzer ana4_Para_BS(hSame[3][ipid][2][1][8], hMix[3][ipid][2][1][8], pidName[ipid]);
            CFAnalyzer ana4_Para_BB(hSame[3][ipid][3][1][8], hMix[3][ipid][3][1][8], pidName[ipid]);

            CFAnalyzer ana4_Anti_SS(hSame[3][ipid][0][2][8], hMix[3][ipid][0][2][8], pidName[ipid]);
            CFAnalyzer ana4_Anti_SB(hSame[3][ipid][1][2][8], hMix[3][ipid][1][2][8], pidName[ipid]);
            CFAnalyzer ana4_Anti_BS(hSame[3][ipid][2][2][8], hMix[3][ipid][2][2][8], pidName[ipid]);
            CFAnalyzer ana4_Anti_BB(hSame[3][ipid][3][2][8], hMix[3][ipid][3][2][8], pidName[ipid]);

            std::cout << "✅ c_cf4: pid= " << pidName[ipid] << " - ヒストグラム読み込み成功" << std::endl;

            //CFを計算し、規格化する
            ana4_Para_SS.Calculate(0.5, 0.8);
            ana4_Para_SB.Calculate(0.5, 0.8);
            ana4_Para_BS.Calculate(0.5, 0.8);
            ana4_Para_BB.Calculate(0.5, 0.8);
            ana4_Anti_SS.Calculate(0.5, 0.8);
            ana4_Anti_SB.Calculate(0.5, 0.8);
            ana4_Anti_BS.Calculate(0.5, 0.8);
            ana4_Anti_BB.Calculate(0.5, 0.8);

            ana4_Para_SS.CorrectCF(ana4_Para_SB.GetNorm(), ana4_Para_BS.GetNorm(), ana4_Para_BB.GetNorm());
            ana4_Anti_SS.CorrectCF(ana4_Anti_SB.GetNorm(), ana4_Anti_BS.GetNorm(), ana4_Anti_BB.GetNorm());

            c_cf4->cd(ipid+1);
            DrawSpin(ana4_Para_SS.GetNormCorr(), ana4_Anti_SS.GetNormCorr(), "Spin Dependence");

            TH1D *hPara_Show = ana4_Para_SS.GetNormCorr(); // 表示用
            if(!hPara_Show) {
                cout << "❌ Error: hPara_Show is null for pid " << pidName[ipid] << "!" << std::endl;
                continue;
            }
            hPara_Show->SetTitle(Form("%s ;Q_{inv} (GeV/c);C(Q_{inv})", pidName[ipid]));
            hPara_Show->GetXaxis()->SetRangeUser(0.0, 0.6); 
            hPara_Show->SetMinimum(0.6); 
            hPara_Show->SetMaximum(1.2); 
            hPara_Show->SetMarkerSize(0.9);
            hPara_Show->SetMarkerStyle(20);

            hPara_Show->Draw("PE");
            SetAxisCF(hPara_Show);

            hPara_Show->GetYaxis()->SetTitleOffset(1.4);

            TLine *lFull = new TLine(0.0, 1.0, 0.6, 1.0);
            lFull->SetLineStyle(2);
            lFull->Draw();

            /*
            TLatex *tex = new TLatex();
            tex->SetNDC();          
            tex->SetTextFont(43);    // フォントサイズをピクセル単位で指定するモード
            tex->SetTextSize(22);    
            tex->SetTextAlign(12);   // 左揃え (1), 上下中央揃え (2)
            tex->SetTextColor(kBlack); 

            // テキストの表示位置 (NDC座標: 左下が(0,0), 右上が(1,1))
            double x_pos = 0.38;     
            double y_pos_start = 0.3; // 1行目のY座標
            double y_step = 0.055;     // 行間のスペース

            tex->DrawLatex(x_pos, y_pos_start, "Run14");
            tex->DrawLatex(x_pos, y_pos_start - y_step, "Au+Au 200 GeV 0#font[122]{\55}80%");
            // ※ #font[122]{\55} はマイナス記号(-)を綺麗に出すためのおまじないです。

            // 斜体(#it)
            tex->DrawLatex(x_pos, y_pos_start - 2 * y_step, "#it{p}_{T} > 0.4(GeV/c), |#it{y}| < 1.0");
        */
        }

        // ---------------------------------------------------------
        // c_cf5 : Overlay Spin (SS, SB, BS, BB) (Low range)
        //----------------------------------------------------------
        c_cf5 = new TCanvas("c_cf5", "CF5",1100, 1100);
        c_cf5->Divide(4, 4);

        CFAnalyzer* ana5_Para[4][4];
        CFAnalyzer* ana5_Anti[4][4];

        for(int ipair = 0; ipair < 4; ipair++) {
            for(int ipid = 0; ipid < 4; ipid++) {
            

                ana5_Para[ipid][ipair] = new CFAnalyzer(hSame[3][ipid][ipair][1][8], hMix[3][ipid][ipair][1][8], pidName[ipid], pL[0], pA[0]); 
                ana5_Anti[ipid][ipair] = new CFAnalyzer(hSame[3][ipid][ipair][2][8], hMix[3][ipid][ipair][2][8], pidName[ipid], pL[0], pA[0]);

                ana5_Para[ipid][ipair]->Calculate(0.5, 0.8); 
                ana5_Anti[ipid][ipair]->Calculate(0.5, 0.8);

                TH1D* h5_Para = ana5_Para[ipid][ipair]->GetNormCorr();
                TH1D* h5_Anti = ana5_Anti[ipid][ipair]->GetNormCorr();

                if (!h5_Para || !h5_Anti) {
                  continue; // nullptrの場合はこのループの処理をスキップ
                }

                TString title_5 = Form("%s %s", pairName[ipair], pidName[ipid]);
                h5_Anti->SetTitle(title_5);
                h5_Para->GetXaxis()->SetRangeUser(0.0, 2.0);
                h5_Anti->GetXaxis()->SetRangeUser(0.0, 2.0);
            
            }
        }
            

        for(int ipid = 0; ipid < 4; ipid++) {

            ana5_Para[ipid][0]->CorrectCF(ana5_Para[ipid][1]->GetNorm(), ana5_Para[ipid][2]->GetNorm(), ana5_Para[ipid][3]->GetNorm());
            ana5_Anti[ipid][0]->CorrectCF(ana5_Anti[ipid][1]->GetNorm(), ana5_Anti[ipid][2]->GetNorm(), ana5_Anti[ipid][3]->GetNorm());
            //cout << "DEBUG: ana5, pid=" << pidName[ipid] << ", Para: " << ana5_Para[ipid][0]->GetNormCorr() << ", Anti: " << ana5_Anti[ipid][0]->GetNormCorr() << endl;

            for(int ipair = 0; ipair < 4; ipair++) {
                c_cf5->cd(ipair*4+ipid+1);
                DrawSpin(ana5_Para[ipid][ipair]->GetNormCorr(), ana5_Anti[ipid][ipair]->GetNormCorr(), Form("%s %s", pidName[ipid], pairName[ipair]));
            }
        }

        // ---------------------------------------------------------
        // c_cf6 : Overlay Correction (Para, Anti)
        //----------------------------------------------------------
        c_cf6 = new TCanvas("c_cf6", "CF6", 1200, 600);
        c_cf6->Divide(4, 2);

        CFAnalyzer* ana_nocut_Para[4][4];
        CFAnalyzer* ana_nocut_Anti[4][4];

        for(int ipair = 0; ipair < 4; ipair++) {
            for(int ipid = 0; ipid < 4; ipid++) {
            
                ana_nocut_Para[ipid][ipair] = new CFAnalyzer(hSame_nocut[ipid][ipair][1][8], hMix_nocut[ipid][ipair][1][8], Form("%s_nocut", pidName[ipid]), pL[0], pA[0]); 
                ana_nocut_Anti[ipid][ipair] = new CFAnalyzer(hSame_nocut[ipid][ipair][2][8], hMix_nocut[ipid][ipair][2][8], Form("%s_nocut", pidName[ipid]), pL[0], pA[0]);

                ana_nocut_Para[ipid][ipair]->Calculate(0.5, 0.8); 
                ana_nocut_Anti[ipid][ipair]->Calculate(0.5, 0.8);
            }
        }

        for(int ipid = 0; ipid < 4; ipid++) {
            
            ana_nocut_Para[ipid][0]->CorrectCF(ana_nocut_Para[ipid][1]->GetNorm(), ana_nocut_Para[ipid][2]->GetNorm(), ana_nocut_Para[ipid][3]->GetNorm());
            ana_nocut_Anti[ipid][0]->CorrectCF(ana_nocut_Anti[ipid][1]->GetNorm(), ana_nocut_Anti[ipid][2]->GetNorm(), ana_nocut_Anti[ipid][3]->GetNorm());

            TH1D *hPara_nocut_Uncorr   = ana_nocut_Para[ipid][0]->GetNormCorr();
            TH1D *hAnti_nocut_Uncorr = ana_nocut_Anti[ipid][0]->GetNormUncorr();
            TH1D *hPara_nocut_Corr   = ana_nocut_Para[ipid][0]->GetNormCorr();
            TH1D *hAnti_nocut_Corr   = ana_nocut_Anti[ipid][0]->GetNormCorr();

            TH1D *hPara_cut_Uncorr = ana5_Para[ipid][0]->GetNormUncorr();
            TH1D *hAnti_cut_Uncorr = ana5_Anti[ipid][0]->GetNormUncorr();
            TH1D *hPara_cut_Corr   = ana5_Para[ipid][0]->GetNormCorr();
            TH1D *hAnti_cut_Corr   = ana5_Anti[ipid][0]->GetNormCorr();

            if (!hPara_nocut_Uncorr || !hAnti_nocut_Corr || !hPara_cut_Uncorr || !hAnti_cut_Corr) {
                std::cout << "Warning: Skipping c_cf6 drawing for pid=" << pidName[ipid] << " due to null pointers." << std::endl;
                continue;
            }

        /*
            if(ipid == 0) {
                hPara_LL = (TH1D*)hPara_Corr->Clone("hPara_LL");
                hAnti_LL = (TH1D*)hAnti_Corr->Clone("hAnti_LL");
            }else if(ipid == 1) {
                hPara_AA = (TH1D*)hPara_Corr->Clone("hPara_AA");
                hAnti_AA = (TH1D*)hAnti_Corr->Clone("hAnti_AA");
            }*/

            //Parallel
            c_cf6->cd(ipid+1);

            hPara_nocut_Uncorr->SetTitle(Form("cos#theta_{12} > 0, %s;Q_{inv} (GeV/c);C(Q_{inv})", pidName[ipid]));
            hPara_nocut_Uncorr->GetXaxis()->SetRangeUser(0.0, 0.5); 
            hPara_nocut_Uncorr->SetMinimum(0.4); 
            hPara_nocut_Uncorr->SetMaximum(1.5);

            // nocut + Uncorr(グレー実線、白抜きの丸)
            hPara_nocut_Uncorr->SetMarkerStyle(24);
            hPara_nocut_Uncorr->SetMarkerColor(kBlack);
            hPara_nocut_Uncorr->SetLineColor(kBlack);
            hPara_nocut_Uncorr->SetMarkerSize(0.7);

            // nocut + Corr(黒い実線、塗りの丸)
            hPara_nocut_Corr->SetMarkerStyle(20);
            hPara_nocut_Corr->SetMarkerColor(kBlack);
            hPara_nocut_Corr->SetLineColor(kBlack);
            hPara_nocut_Corr->SetMarkerSize(0.7);
            
            // cut + Corr(グレーの実線、白抜きの丸)
            hPara_cut_Uncorr->SetMarkerStyle(24); 
            hPara_cut_Uncorr->SetMarkerColor(kRed);
            hPara_cut_Uncorr->SetLineColor(kRed);
            hPara_cut_Uncorr->SetMarkerSize(0.7);
            
            // cut + Corr (赤い実線、塗りの丸)
            hPara_cut_Corr->SetMarkerStyle(20);
            hPara_cut_Corr->SetMarkerColor(kRed);
            hPara_cut_Corr->SetLineColor(kRed);
            hPara_cut_Corr->SetMarkerSize(0.7);
        
            hPara_nocut_Uncorr->Draw("PE");
            hPara_cut_Uncorr->Draw("PE same");
            hPara_nocut_Corr->Draw("PE same");
            hPara_cut_Corr->Draw("PE same");
            SetAxisCF(hPara_nocut_Uncorr); 
            
            TLegend *legC1 = new TLegend(0.4, 0.24, 0.75, 0.42);
            legC1->SetBorderSize(0);
            legC1->SetFillStyle(0);
            legC1->SetTextSize(0.045);
            legC1->AddEntry(hPara_nocut_Uncorr, "Nocut + Before Corr.", "lp");
            legC1->AddEntry(hPara_nocut_Corr, "Nocut + Purity Corr.", "lp");
            legC1->AddEntry(hPara_cut_Uncorr, "Cut + Before Corr.", "lp");
            legC1->AddEntry(hPara_cut_Corr, "Cut + Purity Corr.", "lp");
            legC1->Draw();
            
            TLine *line_comp = new TLine(0.0, 1.0, 0.5, 1.0);
            line_comp->SetLineStyle(2);
            line_comp->Draw();

            //Anti-Parallel
            c_cf6->cd(ipid+5);
        
            hAnti_nocut_Uncorr->SetTitle(Form("cos#theta_{12} < 0 (Anti-Parallel), %s;Q_{inv} (GeV/c);C(Q_{inv})", pidName[ipid]));
            hAnti_nocut_Uncorr->GetXaxis()->SetRangeUser(0.0, 0.5); 
            hAnti_nocut_Uncorr->SetMinimum(0.4); 
            hAnti_nocut_Uncorr->SetMaximum(1.5);

            // nocut + Uncorr (グレ０、白抜きの丸)
            hAnti_nocut_Uncorr->SetMarkerStyle(24);
            hAnti_nocut_Uncorr->SetMarkerColor(kBlack);
            hAnti_nocut_Uncorr->SetLineColor(kBlack);
            hAnti_nocut_Uncorr->SetMarkerSize(0.7);
            
            // nocut + Corr (黒、塗りの丸)
            hAnti_nocut_Corr->SetMarkerStyle(20);
            hAnti_nocut_Corr->SetMarkerColor(kBlack);
            hAnti_nocut_Corr->SetLineColor(kBlack);
            hAnti_nocut_Corr->SetMarkerSize(0.7);

            // cut + Uncorr (グレーの実線、白抜きの丸)
            hAnti_cut_Uncorr->SetMarkerStyle(24); 
            hAnti_cut_Uncorr->SetMarkerColor(kBlue);
            hAnti_cut_Uncorr->SetLineColor(kBlue);
            hAnti_cut_Uncorr->SetMarkerSize(0.7);
            
            // cut + Corr (青い実線、塗りの丸)
            hAnti_cut_Corr->SetMarkerStyle(20);
            hAnti_cut_Corr->SetMarkerColor(kBlue);
            hAnti_cut_Corr->SetLineColor(kBlue);
            hAnti_cut_Corr->SetMarkerSize(0.7);
        
            hAnti_nocut_Uncorr->Draw("PE");
            hAnti_cut_Uncorr->Draw("PE same");
            hAnti_nocut_Corr->Draw("PE same");
            hAnti_cut_Corr->Draw("PE same");
            SetAxisCF(hAnti_nocut_Uncorr);
            
            TLegend *legC2 = new TLegend(0.4, 0.24, 0.75, 0.42);
            legC2->SetBorderSize(0);
            legC2->SetFillStyle(0);
            legC2->SetTextSize(0.045);
            legC2->AddEntry(hAnti_nocut_Uncorr, "Nocut + Before Corr.", "lp");
            legC2->AddEntry(hAnti_nocut_Corr, "Nocut + Purity Corr.", "lp");
            legC2->AddEntry(hAnti_cut_Uncorr, "Cut + Before Corr.", "lp");
            legC2->AddEntry(hAnti_cut_Corr, "Cut + Purity Corr.", "lp");
            legC2->Draw();
            line_comp->Draw();

            hPara_cut_Corr->SetName(Form("hCF_Para_%s", pidName[ipid]));
            hAnti_cut_Corr->SetName(Form("hCF_Anti_%s", pidName[ipid]));
            
            hSave_CF_Para[ipid] = (TH1D*)hPara_cut_Corr->Clone(Form("hCF_Para_%s", pidName[ipid]));
            hSave_CF_Anti[ipid] = (TH1D*)hAnti_cut_Corr->Clone(Form("hCF_Anti_%s", pidName[ipid]));

            hSave_Same_Para[ipid] = (TH1D*)hSame[4][ipid][0][1][8]->Clone(Form("hSame_Para_%s", pidName[ipid]));
            hSave_Mix_Para[ipid]  = (TH1D*)hMix[4][ipid][0][1][8]->Clone(Form("hMix_Para_%s", pidName[ipid]));

            hSave_Same_Anti[ipid] = (TH1D*)hSame[4][ipid][0][2][8]->Clone(Form("hSame_Anti_%s", pidName[ipid]));
            hSave_Mix_Anti[ipid]  = (TH1D*)hMix[4][ipid][0][2][8]->Clone(Form("hMix_Anti_%s", pidName[ipid]));
            }

        
        // ---------------------------------------------------------
        // c_cf6_2 : Overlay Cut Correction (Para, Anti)
        //----------------------------------------------------------
        c_cf6_2 = new TCanvas("c_cf6_2", "CF6_2", 1200, 600);
        c_cf6_2->Divide(4, 2);

        CFAnalyzer* ana_cut_Para[4][4][4];   CFAnalyzer* ana_cut_Anti[4][4][4];
        CFAnalyzer* ana_Nocut_Para[4][4];    CFAnalyzer* ana_Nocut_Anti[4][4];

        const char* icutName[4] = {"pp", "#pi#pi", "p#pi", "#pip"};

        for(int ipair = 0; ipair < 4; ipair++) {
            for(int ipid = 0; ipid < 4; ipid++) {
                for(int icut = 0; icut < 4; icut++) {
            
                ana_cut_Para[icut][ipid][ipair] = new CFAnalyzer(hSame[icut][ipid][ipair][1][8], hMix[icut][ipid][ipair][1][8], Form("%s_%s", pidName[ipid], icutName[icut]), pL[0], pA[0]); 
                ana_cut_Anti[icut][ipid][ipair] = new CFAnalyzer(hSame[icut][ipid][ipair][2][8], hMix[icut][ipid][ipair][2][8], Form("%s_%s", pidName[ipid], icutName[icut]), pL[0], pA[0]);

                ana_cut_Para[icut][ipid][ipair]->Calculate(0.5, 0.8); 
                ana_cut_Anti[icut][ipid][ipair]->Calculate(0.5, 0.8);
                }
                ana_Nocut_Para[ipid][ipair] = new CFAnalyzer(hSame_nocut[ipid][ipair][1][8], hMix_nocut[ipid][ipair][1][8], Form("%s", pidName[ipid]), pL[0], pA[0]); 
                ana_Nocut_Anti[ipid][ipair] = new CFAnalyzer(hSame_nocut[ipid][ipair][2][8], hMix_nocut[ipid][ipair][2][8], Form("%s", pidName[ipid]), pL[0], pA[0]);

                ana_Nocut_Para[ipid][ipair]->Calculate(0.5, 0.8); 
                ana_Nocut_Anti[ipid][ipair]->Calculate(0.5, 0.8);

            }
        }
        TH1D *hPara_cut[4], *hAnti_cut[4];
        TH1D *hPara_nocut, *hAnti_nocut;

        for(int ipid = 0; ipid < 4; ipid++) {
            for(int icut = 0; icut < 4; icut++) {
            
                ana_cut_Para[icut][ipid][0]->CorrectCF(ana_cut_Para[icut][ipid][1]->GetNorm(), ana_cut_Para[icut][ipid][2]->GetNorm(), ana_cut_Para[icut][ipid][3]->GetNorm());
                ana_cut_Anti[icut][ipid][0]->CorrectCF(ana_cut_Anti[icut][ipid][1]->GetNorm(), ana_cut_Anti[icut][ipid][2]->GetNorm(), ana_cut_Anti[icut][ipid][3]->GetNorm());

                hPara_cut[icut]   = ana_cut_Para[icut][ipid][0]->GetNormCorr();
                hAnti_cut[icut]   = ana_cut_Anti[icut][ipid][0]->GetNormCorr();

                if (!hPara_cut[icut]) {
                    std::cout << "Warning: Skipping c_cf6_2 drawing for pid=" << pidName[ipid] << " due to null pointers." << std::endl;
                    continue;
                }
            }
            ana_Nocut_Para[ipid][0]->CorrectCF(ana_Nocut_Para[ipid][1]->GetNorm(), ana_Nocut_Para[ipid][2]->GetNorm(), ana_Nocut_Para[ipid][3]->GetNorm());
            ana_Nocut_Anti[ipid][0]->CorrectCF(ana_Nocut_Anti[ipid][1]->GetNorm(), ana_Nocut_Anti[ipid][2]->GetNorm(), ana_Nocut_Anti[ipid][3]->GetNorm()); 
        
            hPara_nocut = ana_Nocut_Para[ipid][0]->GetNormCorr();
            hAnti_nocut = ana_Nocut_Anti[ipid][0]->GetNormCorr();
            if (!hPara_nocut || !hAnti_nocut) {
                std::cout << "Warning: Skipping c_cf6_2 drawing for pid=" << pidName[ipid] << " due to null pointers." << std::endl;
                continue;
            }

            //Parallel
            c_cf6_2->cd(ipid+1);

            hPara_nocut->SetTitle(Form("cos#theta_{12} > 0, %s;Q_{inv} (GeV/c);C(Q_{inv})", pidName[ipid]));
            hPara_nocut->GetXaxis()->SetRangeUser(0.0, 0.5); 
            hPara_nocut->SetMinimum(0.4); 
            hPara_nocut->SetMaximum(1.5);

            // nocut + Corr(黒い実線、塗りの丸)
            hPara_nocut->SetMarkerStyle(24);
            hPara_nocut->SetMarkerColor(kBlack);
            hPara_nocut->SetLineColor(kBlack);
            hPara_nocut->SetMarkerSize(0.6);

            hPara_nocut->Draw("PE");
            SetAxisCF(hPara_nocut);

            int cutColors[4] = {kMagenta, kOrange+2, kBlue, kGreen+2};
            int cutStyles[4] = {25, 26, 27, 28}; // Marker styles for cuts

            for(int icut = 0; icut < 4; icut++) {
                hPara_cut[icut]->SetMarkerStyle(cutStyles[icut]);
                hPara_cut[icut]->SetMarkerSize(0.6);
                hPara_cut[icut]->SetMarkerColor(cutColors[icut]);
                hPara_cut[icut]->SetLineColor(cutColors[icut]);
                hPara_cut[icut]->Draw("PE same");
            }
            
            TLegend *legC1_1 = new TLegend(0.4, 0.24, 0.75, 0.45);
            legC1_1->SetBorderSize(0);
            legC1_1->SetFillStyle(0);
            legC1_1->SetTextSize(0.045);
            legC1_1->AddEntry(hPara_nocut, "Nocut", "lp");
            legC1_1->AddEntry(hPara_cut[0], "pp", "lp");
            legC1_1->AddEntry(hPara_cut[1], "#pi#pi", "lp");
            legC1_1->AddEntry(hPara_cut[2], "p#pi", "lp");
            legC1_1->AddEntry(hPara_cut[3], "#pip", "lp");
            legC1_1->Draw();
            
            TLine *line_1 = new TLine(0.0, 1.0, 0.5, 1.0);
            line_1->SetLineStyle(2);
            line_1->Draw();

            //Anti-Parallel
            c_cf6_2->cd(ipid+5);
        
            hAnti_nocut->SetTitle(Form("cos#theta_{12} < 0 (Anti-Parallel), %s;Q_{inv} (GeV/c);C(Q_{inv})", pidName[ipid]));
            hAnti_nocut->GetXaxis()->SetRangeUser(0.0, 0.5); 
            hAnti_nocut->SetMinimum(0.4); 
            hAnti_nocut->SetMaximum(1.5);
            
            // nocut + Corr (黒、塗りの丸)
            hAnti_nocut->SetMarkerStyle(24);
            hAnti_nocut->SetMarkerColor(kBlack);
            hAnti_nocut->SetLineColor(kBlack);
            hAnti_nocut->SetMarkerSize(0.7);

            hAnti_nocut->Draw("PE");
            SetAxisCF(hAnti_nocut);

            for(int icut = 0; icut < 4; icut++) {
                hAnti_cut[icut]->SetMarkerStyle(cutStyles[icut]);
                hAnti_cut[icut]->SetMarkerSize(0.7);
                hAnti_cut[icut]->SetMarkerColor(cutColors[icut]);
                hAnti_cut[icut]->SetLineColor(cutColors[icut]);
                hAnti_cut[icut]->Draw("PE same");
            }

            TLegend *legC2_1 = new TLegend(0.4, 0.24, 0.75, 0.45);
            legC2_1->SetBorderSize(0);
            legC2_1->SetFillStyle(0);
            legC2_1->SetTextSize(0.045);
            legC2_1->AddEntry(hAnti_nocut, "Nocut ", "lp");
            legC2_1->AddEntry(hAnti_cut[0], "pp", "lp");
            legC2_1->AddEntry(hAnti_cut[1], "#pi#pi", "lp");
            legC2_1->AddEntry(hAnti_cut[2], "p#pi", "lp");
            legC2_1->AddEntry(hAnti_cut[3], "#pip", "lp");
            legC2_1->Draw();
            line_1->Draw();

            
            }
/*
        // ---------------------------------------------------------
        // c_cf7 : Overlay Correction (Centrality)
        //----------------------------------------------------------
        CFAnalyzer* ana7_Para[9][4];   CFAnalyzer* ana7_Anti[9][4];

        for(int ipid = 0; ipid < 4; ipid++) {
            for(int icent = 0; icent < 9; icent++) {
                for(int ipair = 0; ipair < 4; ipair++){

                    ana7_Para[icent][ipair] = new CFAnalyzer(hSame[4][ipid][ipair][1][icent], hMix[4][ipid][ipair][1][icent], Form("Para %s, %s, %d", pidName[ipid], pairName[ipair], icent)); 
                    ana7_Anti[icent][ipair] = new CFAnalyzer(hSame[4][ipid][ipair][2][icent], hMix[4][ipid][ipair][2][icent], Form("Anti %s, %s, %d", pidName[ipid], pairName[ipair], icent));
                    ana7_Para[icent][ipair]->Calculate(0.5, 0.8); 
                    ana7_Anti[icent][ipair]->Calculate(0.5, 0.8);
                }
            
                ana7_Para[icent][0]->CorrectCF(ana7_Para[icent][1]->GetNorm(), ana7_Para[icent][2]->GetNorm(), ana7_Para[icent][3]->GetNorm());
                ana7_Anti[icent][0]->CorrectCF(ana7_Anti[icent][1]->GetNorm(), ana7_Anti[icent][2]->GetNorm(), ana7_Anti[icent][3]->GetNorm());

                TH1* h7_Para = ana7_Para[icent][0]->GetNormCorr(); // 補正後のSSを取得
                TH1* h7_Anti = ana7_Anti[icent][0]->GetNormCorr();

                if (!h7_Para || !h7_Anti) {
                    cout << "h7 does not exist" << endl;
                    continue;
                }

                TString title_7 = Form("%s %d", pidName[ipid], icent);
                h7_Anti->SetTitle(title_7);
                h7_Para->GetXaxis()->SetRangeUser(0.0, 0.6);
                h7_Anti->GetXaxis()->SetRangeUser(0.0, 0.6);
            }

            c_cf7[ipid] = new TCanvas(Form("c_cf7_%s", pidName[ipid]), "CF7", 1200, 600);
            c_cf7[ipid]->Divide(4, 2);

            for(int icent = 0; icent < 8; icent++) {
                c_cf7[ipid]->cd(icent+1);
                DrawSpin(ana7_Para[icent][0]->GetNormCorr(), ana7_Anti[icent][0]->GetNormCorr(), Form("%s, %s", pidName[ipid], centName[icent]));
            }
        }

        // ---------------------------------------------------------
        // c_cf8 : Overlay Correction (Centrality 0-20, 20-50. 50-80%)
        //----------------------------------------------------------
        c_cf8 = new TCanvas("c_cf8", "CF8", 1300, 900);
        c_cf8->Divide(4, 3);

        for(int ipid = 0; ipid < 4; ipid++) {

            CFAnalyzer* ana8_Para[3][4]; // ic(centrality), ipair
            CFAnalyzer* ana8_Anti[3][4];

            for(int ipair = 0; ipair < 4; ipair++){

                // 0-20%
                ana8_Para[0][ipair] = new CFAnalyzer(hSame_0to20per[ipid][ipair][1], hMix_0to20per[ipid][ipair][1], Form("Para %s 0-20%% %s", pidName[ipid], pairName[ipair])); 
                ana8_Anti[0][ipair] = new CFAnalyzer(hSame_0to20per[ipid][ipair][2], hMix_0to20per[ipid][ipair][2], Form("Anti %s 0-20%% %s", pidName[ipid], pairName[ipair]));
                // 20-50%
                ana8_Para[1][ipair] = new CFAnalyzer(hSame_20to50per[ipid][ipair][1], hMix_20to50per[ipid][ipair][1], Form("Para %s 20-50 %s", pidName[ipid], pairName[ipair])); 
                ana8_Anti[1][ipair] = new CFAnalyzer(hSame_20to50per[ipid][ipair][2], hMix_20to50per[ipid][ipair][2], Form("Anti %s 20-50 %s", pidName[ipid], pairName[ipair]));
                // 50-80%
                ana8_Para[2][ipair] = new CFAnalyzer(hSame_50to80per[ipid][ipair][1], hMix_50to80per[ipid][ipair][1], Form("Para %s 50-80 %s", pidName[ipid], pairName[ipair])); 
                ana8_Anti[2][ipair] = new CFAnalyzer(hSame_50to80per[ipid][ipair][2], hMix_50to80per[ipid][ipair][2], Form("Anti %s 50-80 %s", pidName[ipid], pairName[ipair]));

                for(int ic3 = 0; ic3 < 3; ic3++) {
                    ana8_Para[ic3][ipair]->Calculate(0.5, 0.8);
                    ana8_Anti[ic3][ipair]->Calculate(0.5, 0.8);
                }
            }

            for(int ic3 = 0; ic3 < 3; ic3++) {
                ana8_Para[ic3][0]->CorrectCF(ana8_Para[ic3][1]->GetNorm(), ana8_Para[ic3][2]->GetNorm(), ana8_Para[ic3][3]->GetNorm());
                ana8_Anti[ic3][0]->CorrectCF(ana8_Anti[ic3][1]->GetNorm(), ana8_Anti[ic3][2]->GetNorm(), ana8_Anti[ic3][3]->GetNorm());

                TH1D* h8_Para = ana8_Para[ic3][0]->GetNormCorr();
                TH1D* h8_Anti = ana8_Anti[ic3][0]->GetNormCorr();

                if (!h8_Para || !h8_Anti ) {
                std::cout << "⚠️ Skip c_cf8 drawing for " << pidName[ipid] << " due to missing histograms." << std::endl;
                continue;
                }

                h8_Anti->SetTitle(pidName[ipid]);
                h8_Para->GetXaxis()->SetRangeUser(0.0, 0.6);
                h8_Anti->GetXaxis()->SetRangeUser(0.0, 0.6);
                
                c_cf8->cd(ipid + 1 + (ic3 * 4));
                DrawSpin(h8_Para, h8_Anti, Form("%s, %s", pidName[ipid], cent3Name[ic3]));
            }
        }*/

        // ---------------------------------------------------------
        // c_cf9 : Overlay (θ : 0-45, 45-90, 90-135, 135-180)
        //----------------------------------------------------------
        c_cf9 = new TCanvas("c_cf9", "CF9", 1500, 400);
        c_cf9->Divide(4, 1);

        for(int ipid = 0; ipid < 4; ipid++) {

            CFAnalyzer* ana9[7][4];
            c_cf9->cd(ipid+1);

            for(int ispin = 3; ispin < 7; ispin++) {
                for(int ipair = 0; ipair < 4; ipair++){

                ana9[ispin][ipair] = new CFAnalyzer(hSame[3][ipid][ipair][ispin][8], hMix[3][ipid][ipair][ispin][8], Form("%s, %s, %s", pidName[ipid], pairName[ipair], spinName7[ispin]), pL[0], pA[0]); 
                ana9[ispin][ipair]->Calculate(0.5, 0.8); 
                }

                ana9[ispin][0]->CorrectCF(ana9[ispin][1]->GetNorm(), ana9[ispin][2]->GetNorm(), ana9[ispin][3]->GetNorm());
                TH1D* h9 = ana9[ispin][0]->GetNormCorr();

                if (!h9) {
                    std::cout << "⚠️ Skip c_cf9 drawing for " << pidName[ipid] << " due to missing histograms." << std::endl;
                    continue;
                }

                TString title_9 = Form("%s", pidName[ipid]);
                h9->SetTitle(title_9);
                h9->GetXaxis()->SetRangeUser(0.0, 0.6);
                h9->SetMarkerSize(0.7);
                
                if(ispin==3) {
                    ana9[ispin][0]->DrawNorm(h9, "PE");
                }else{
                    ana9[ispin][0]->DrawNorm(h9, "PE same");
                } 
                
                if(ispin == 3) {
                    h9->SetMarkerStyle(20);
                    h9->SetMarkerColor(kRed);
                    h9->SetLineColor(kRed);
                } else if(ispin == 4) {
                    h9->SetMarkerStyle(24);
                    h9->SetMarkerColor(kOrange+2);
                    h9->SetLineColor(kOrange+2);
                } else if(ispin == 5) {
                    h9->SetMarkerStyle(24);
                    h9->SetMarkerColor(kGreen+2);
                    h9->SetLineColor(kGreen+2);
                } else if(ispin == 6) {
                    h9->SetMarkerStyle(23);
                    h9->SetMarkerColor(kBlue);
                    h9->SetLineColor(kBlue);
                }
            }
            // 凡例を描画 (各角度帯を示す)
            TLegend *leg9 = new TLegend(0.58, 0.22, 0.88, 0.42);
            leg9->SetBorderSize(0);
            leg9->SetFillStyle(0);
            // ダミーのグラフ用に同じスタイルのマーカーを使ってエントリを作成
            TH1D *hleg_r = new TH1D("hleg_r", "", 1, 0, 1); hleg_r->SetMarkerStyle(20); hleg_r->SetMarkerColor(kRed);
            TH1D *hleg_m = new TH1D("hleg_m", "", 1, 0, 1); hleg_m->SetMarkerStyle(24); hleg_m->SetMarkerColor(kOrange+2);
            TH1D *hleg_g = new TH1D("hleg_g", "", 1, 0, 1); hleg_g->SetMarkerStyle(24); hleg_g->SetMarkerColor(kGreen+2);
            TH1D *hleg_b = new TH1D("hleg_b", "", 1, 0, 1); hleg_b->SetMarkerStyle(23); hleg_b->SetMarkerColor(kBlue);
            leg9->AddEntry(hleg_r, "0-45 deg", "lp");
            leg9->AddEntry(hleg_m, "45-90 deg", "lp");
            leg9->AddEntry(hleg_g, "90-135 deg", "lp");
            leg9->AddEntry(hleg_b, "135-180 deg", "lp");
            leg9->Draw();
        }

        // ---------------------------------------------------------
        // c_cf10 : Parallel and Anti-Parallel (Centrality 10-80%)
        //----------------------------------------------------------
        c_cf10 = new TCanvas("c_cf10", "CF10_10-80%", 1300, 900);
        c_cf10->Divide(4, 3);

        CFAnalyzer* ana_Para_10to80[4][4]; // [ipid][ipair]
        CFAnalyzer* ana_Anti_10to80[4][4];
        CFAnalyzer* ana_Para_20to80[4][4]; // [ipid][ipair]
        CFAnalyzer* ana_Anti_20to80[4][4];

        TH1D* hPara_10to80[4] = {nullptr};  TH1D* hAnti_10to80[4] = {nullptr};
        TH1D* hPara_20to80[4] = {nullptr};  TH1D* hAnti_20to80[4] = {nullptr};

        for(int ipid = 0; ipid < 4; ipid++) {
            for(int ipair = 0; ipair < 4; ipair++){

                ana_Para_10to80[ipid][ipair] = new CFAnalyzer(hSame_10to80per[ipid][ipair][1], hMix_10to80per[ipid][ipair][1], Form("Para %s 10-80%% %s", pidName[ipid], pairName[ipair]), pL[1], pA[1]); // [ipid][ipair][ispin]
                ana_Anti_10to80[ipid][ipair] = new CFAnalyzer(hSame_10to80per[ipid][ipair][2], hMix_10to80per[ipid][ipair][2], Form("Anti %s 10-80%% %s", pidName[ipid], pairName[ipair]), pL[1], pA[1]);
                ana_Para_20to80[ipid][ipair] = new CFAnalyzer(hSame_20to80per[ipid][ipair][1], hMix_20to80per[ipid][ipair][1], Form("Para %s 20-80%% %s", pidName[ipid], pairName[ipair]), pL[2], pA[2]); // [ipid][ipair][ispin]
                ana_Anti_20to80[ipid][ipair] = new CFAnalyzer(hSame_20to80per[ipid][ipair][2], hMix_20to80per[ipid][ipair][2], Form("Anti %s 20-80%% %s", pidName[ipid], pairName[ipair]), pL[2], pA[2]);

                ana_Para_10to80[ipid][ipair]->Calculate(0.5, 0.8);
                ana_Anti_10to80[ipid][ipair]->Calculate(0.5, 0.8);   
                ana_Para_20to80[ipid][ipair]->Calculate(0.5, 0.8);
                ana_Anti_20to80[ipid][ipair]->Calculate(0.5, 0.8);  
            }

            
            ana_Para_10to80[ipid][0]->CorrectCF(ana_Para_10to80[ipid][1]->GetNorm(), ana_Para_10to80[ipid][2]->GetNorm(), ana_Para_10to80[ipid][3]->GetNorm());
            ana_Anti_10to80[ipid][0]->CorrectCF(ana_Anti_10to80[ipid][1]->GetNorm(), ana_Anti_10to80[ipid][2]->GetNorm(), ana_Anti_10to80[ipid][3]->GetNorm());
            ana_Para_20to80[ipid][0]->CorrectCF(ana_Para_20to80[ipid][1]->GetNorm(), ana_Para_20to80[ipid][2]->GetNorm(), ana_Para_20to80[ipid][3]->GetNorm());
            ana_Anti_20to80[ipid][0]->CorrectCF(ana_Anti_20to80[ipid][1]->GetNorm(), ana_Anti_20to80[ipid][2]->GetNorm(), ana_Anti_20to80[ipid][3]->GetNorm());

            hPara_10to80[ipid] = ana_Para_10to80[ipid][0]->GetNormCorr();
            hAnti_10to80[ipid] = ana_Anti_10to80[ipid][0]->GetNormCorr();
            hPara_20to80[ipid] = ana_Para_20to80[ipid][0]->GetNormCorr();
            hAnti_20to80[ipid] = ana_Anti_20to80[ipid][0]->GetNormCorr();


            if (!hPara_10to80[ipid] || !hAnti_10to80[ipid] ) {
                std::cout << "⚠️ Skip c_cf10 drawing for " << pidName[ipid] << " due to missing histograms." << std::endl;
                continue;
            }

            hPara_10to80[ipid]->GetXaxis()->SetRangeUser(0.0, 0.6);
            hAnti_10to80[ipid]->GetXaxis()->SetRangeUser(0.0, 0.6);    

            c_cf10->cd(ipid+1);
            DrawSpin(ana5_Para[ipid][0]->GetNormCorr(), ana5_Anti[ipid][0]->GetNormCorr(), Form("%s 0-80%%", pidName[ipid]));
        
            c_cf10->cd(ipid+5);
            DrawSpin(hPara_10to80[ipid], hAnti_10to80[ipid], Form("%s 10-80%%", pidName[ipid]));

            c_cf10->cd(ipid+9);
            DrawSpin(hPara_20to80[ipid], hAnti_20to80[ipid], Form("%s 20-80%%", pidName[ipid]));
        }

    
    TString outFileName = "check_isobar_cutR0p06.root";
    TFile *f_my = new TFile(outFileName, "UPDATE");

    for(int ipid = 0; ipid < 4; ipid++) {
    
        hSave_CF_All[ipid]->Write("", TObject::kOverwrite);
        hSave_CF_Para[ipid]->Write("", TObject::kOverwrite);
        hSave_CF_Anti[ipid]->Write("", TObject::kOverwrite);

        hSave_Same_All[ipid]->Write("", TObject::kOverwrite);
        hSave_Mix_All[ipid]->Write("", TObject::kOverwrite);
        hSave_Same_Para[ipid]->Write("", TObject::kOverwrite);
        hSave_Mix_Para[ipid]->Write("", TObject::kOverwrite);
        hSave_Same_Anti[ipid]->Write("", TObject::kOverwrite);
        hSave_Mix_Anti[ipid]->Write("", TObject::kOverwrite);
    }
    
    f_my->Close();
    std::cout << "✅ Saved histograms to " << outFileName << std::endl;
       


    //=============Save Canvas=============
    /*
    c_cf1[ipid]->SaveAs(Form("CF_Norm_Pair_%s.pdf", pidName[ipid]));
    c_cf2[ipid]->SaveAs(Form("CF_Norm_Spin_%s.pdf", pidName[ipid]));
    c_cf3[ipid]->SaveAs(Form("CF_Normed_Full_%s.pdf", pidName[ipid]));
    c_cf4[ipid]->SaveAs(Form("CF_Spin_Full_run16_%s.pdf", pidName[ipid]));
    c_cf5[ipid]->SaveAs(Form("CF_Spin_Low_run16_%s.pdf", pidName[ipid]));
    c_cf6[ipid]->SaveAs(Form("CF_Corr_Spin_run16_%s.pdf", pidName[ipid]));
    c_cf7[ipid]->SaveAs(Form("CF_Spin_Cent_run16_%s.pdf", pidName[ipid]));
    c_cf8[ipid]->SaveAs(Form("CF_Spin_Cent3_run16_%s.pdf", pidName[ipid]));
    c_cf5->SaveAs("CF_isobar_cut_SB.pdf");
    c_cf6->SaveAs("CF_isobar_cut_CutCorr.pdf");
    c_cf9->SaveAs("CF_isobar_cut_SpinAngle.pdf");*/
    //c_cf10->SaveAs("CF_isobar_10to80_mycut_cutfakeV0_full.pdf");
    } // end of Macro

    /*
     // 1. Parallel (cos theta > 0) の比較
        Double_t p_para    = hPara_LL->Chi2Test(hPara_AA, "WW");
        Double_t chi2_para = hPara_LL->Chi2Test(hPara_AA, "WW CHI2");
        Double_t ndf_para  = hPara_LL->Chi2Test(hPara_AA, "WW CHI2/NDF");
        std::cout << "[Parallel]" << std::endl;
        std::cout << "  p-value  : " << p_para << std::endl;
        std::cout << "  Chi2/ndf : " << chi2_para << " / " << (chi2_para / (ndf_para ? ndf_para : 1.0)) << " = " << ndf_para << std::endl;

        // 2. Anti-Parallel (cos theta < 0) の比較
        Double_t p_anti    = hAnti_LL->Chi2Test(hAnti_AA, "WW");
        Double_t chi2_anti = hAnti_LL->Chi2Test(hAnti_AA, "WW CHI2");
        Double_t ndf_anti  = hAnti_LL->Chi2Test(hAnti_AA, "WW CHI2/NDF");
        std::cout << "[Anti-Parallel]" << std::endl;
        std::cout << "  p-value  : " << p_anti << std::endl;
        std::cout << "  Chi2/ndf : " << chi2_anti << " / " << (chi2_anti / (ndf_anti ? ndf_anti : 1.0)) << " = " << ndf_anti << std::endl;
        std::cout << "=================================================\n" << std::endl;
*/



