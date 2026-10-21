#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"
#include "TLatex.h"
#include "AxisStyle.h"
#include "CFCalculator.h"
#include <iostream>
#include <cmath>

void Overlay_isobar_cut(){

    TFile *f_ref_Para = new TFile("../ana_KFParticle/angel0_90.root", "READ");
    TFile *f_ref_Anti = new TFile("../ana_KFParticle/angel90_180.root", "READ");
    TFile *f_ref_Para_smallbin = new TFile("../ana_KFParticle/smallerbin0_90.root", "READ");
    TFile *f_ref_Anti_smallbin = new TFile("../ana_KFParticle/smallerbin90_180.root", "READ");
/*
    TFile *f_my_S       = new TFile("check_isobar_cutS.root", "READ"); // Δη<0.02, Δφ*<0.02
    TFile *f_my_R0p04   = new TFile("check_isobar_cutR0p04.root", "READ"); // Δη<0.04, Δφ*<0.04
    TFile *f_my_R0p06   = new TFile("check_isobar_cutR0p06.root", "READ"); // Δη<0.06, Δφ*<0.06
   */
    TFile *f_my_phi7eta6   = new TFile("check_isobar_cutphi7eta6_R1p4.root", "READ"); // Δη<0.06, Δφ*<0.07
    TFile *f_my_phi5eta4   = new TFile("check_isobar_cutphi5eta4_R1p4.root", "READ"); // Δη<0.06, Δφ*<0.05
    TFile *f_my_phi3eta2   = new TFile("check_isobar_cutphi3eta2_R1p4.root", "READ"); // Δη<0.06, Δφ*<0.03

    if (!f_ref_Para || f_ref_Para->IsZombie() || !f_ref_Anti || f_ref_Anti_smallbin->IsZombie()) {
        std::cout << "!!! Error: Cannot open one or more input Kehao files !!!" << std::endl;
    }
/*
    if(!f_my_S || f_my_S->IsZombie() || !f_my_R0p04 || f_my_R0p04->IsZombie() || !f_my_R0p06 || f_my_R0p06->IsZombie()) {
        std::cout << "!!! Error: Cannot open my files !!!" << std::endl;
        return;
    }*/

    TH1D *h_ref[2] = {nullptr, nullptr}; // 0: Para, 1: Anti (Qinv* スケール, Max 5.0 GeV/c)
    TH1D *h_ref_smallbin[2] = {nullptr, nullptr};
    /*TH1D *h_my_cutS[2];
    TH1D *h_my_R0p04[2];
    TH1D *h_my_R0p06[2];*/
    TH1D *h_my_phi7eta6[2] = {nullptr, nullptr};
    TH1D *h_my_phi5eta4[2] = {nullptr, nullptr};
    TH1D *h_my_phi3eta2[2] = {nullptr, nullptr};

    h_ref[0] = (TH1D*)f_ref_Para->Get("hCF_cut2_s_cent0");
    h_ref[1] = (TH1D*)f_ref_Anti->Get("hCF_cut2_s_cent0");
    h_ref_smallbin[0] = (TH1D*)f_ref_Para_smallbin->Get("hCF_cut2_s_cent0");
    h_ref_smallbin[1] = (TH1D*)f_ref_Anti_smallbin->Get("hCF_cut2_s_cent0");/*
    h_my_cutS[0]  = (TH1D*)f_my_S->Get("hCF_Para_LL+AA");
    h_my_cutS[1]  = (TH1D*)f_my_S->Get("hCF_Anti_LL+AA");
    h_my_R0p04[0]  = (TH1D*)f_my_R0p04->Get("hCF_Para_LL+AA");
    h_my_R0p04[1]  = (TH1D*)f_my_R0p04->Get("hCF_Anti_LL+AA");
    h_my_R0p06[0]  = (TH1D*)f_my_R0p06->Get("hCF_Para_LL+AA");
    h_my_R0p06[1]  = (TH1D*)f_my_R0p06->Get("hCF_Anti_LL+AA");*/
    h_my_phi7eta6[0]  = (TH1D*)f_my_phi7eta6->Get("hCF_Para_LL+AA");
    h_my_phi7eta6[1]  = (TH1D*)f_my_phi7eta6->Get("hCF_Anti_LL+AA");
    h_my_phi5eta4[0]  = (TH1D*)f_my_phi5eta4->Get("hCF_Para_LL+AA");
    h_my_phi5eta4[1]  = (TH1D*)f_my_phi5eta4->Get("hCF_Anti_LL+AA");
    h_my_phi3eta2[0]  = (TH1D*)f_my_phi3eta2->Get("hCF_Para_LL+AA");
    h_my_phi3eta2[1]  = (TH1D*)f_my_phi3eta2->Get("hCF_Anti_LL+AA");

    for (int j = 0; j < 2; j++) {
        if (!h_my_phi7eta6[j] || !h_my_phi5eta4[j] || !h_my_phi3eta2[j]) {
            std::cout << "!!! Error: My histograms not found for index " << j << " !!!" << std::endl;
            return; // 自分のヒストグラムが無い場合はどうしようもないので終了
        }
        if (!h_ref[j] || !h_ref_smallbin[j]) {
            std::cout << "--- Notice: Kehao histogram not found for index " << j << " (will draw only my hists) ---" << std::endl;
        }
    }

    double scale = 1000.0; 
    for (int j = 0; j < 2; j++) {
        if (h_ref_smallbin[j]) {
            double ref_xmin = h_ref_smallbin[j]->GetXaxis()->GetXmin();
            double ref_xmax = h_ref_smallbin[j]->GetXaxis()->GetXmax();
            // 軸の限界値をスケーリング
            h_ref_smallbin[j]->GetXaxis()->SetLimits(ref_xmin / scale, ref_xmax / scale);
        }
    }


    // Kehao (20: 塗りつぶし丸)
    for (int j = 0; j < 2; j++) {
        if (h_ref_smallbin[j]) {
            h_ref_smallbin[j]->SetMarkerStyle(24);
            h_ref_smallbin[j]->SetMarkerSize(0.8);
            h_ref_smallbin[j]->SetMarkerColor(kMagenta);
            h_ref_smallbin[j]->SetLineColor(kMagenta);
        }

        h_my_phi7eta6[j]->SetMarkerStyle(20);
        h_my_phi7eta6[j]->SetMarkerSize(0.7);
        h_my_phi5eta4[j]->SetMarkerStyle(20);
        h_my_phi5eta4[j]->SetMarkerSize(0.7);
        h_my_phi3eta2[j]->SetMarkerStyle(20);
        h_my_phi3eta2[j]->SetMarkerSize(0.7);

        h_my_phi7eta6[j]->SetMarkerColor(kRed);
        h_my_phi7eta6[j]->SetLineColor(kRed);
        h_my_phi5eta4[j]->SetMarkerColor(kGreen+2);
        h_my_phi5eta4[j]->SetLineColor(kGreen+2);
        h_my_phi3eta2[j]->SetMarkerColor(kBlue);
        h_my_phi3eta2[j]->SetLineColor(kBlue);

    }


    const char* titles[2] = {"Para", "Anti"};
    double y_min[2] = {0.6, 0.8}; 
    double x_max[2] = {2.0, 0.6};

    for(int i=0; i<2; i++){
        TCanvas *c1_ol = new TCanvas(Form("c1_ol[%d]", i), "(0.0, 2.0)GeV/c", 1600, 500);
        c1_ol->Divide(2, 1);

        for (int j = 0; j < 2; j++) {
            c1_ol->cd(j + 1);

            TH1D *h_draw = nullptr;
            if (h_ref_smallbin[j]) {
                h_draw = (TH1D*)h_ref_smallbin[j]->Clone(Form("h_base_clone_%d_%d", i, j));
            } else {
                h_draw = (TH1D*)h_my_phi7eta6[j]->Clone(Form("h_base_clone_%d_%d", i, j));
            }
            h_draw->SetTitle(Form("CF %s; Q^{inv} (GeV/c); C(Q_{inv})", titles[j]));

            h_draw->GetXaxis()->SetRangeUser(0.0, x_max[i]); 
            h_draw->SetMaximum(1.1);
            h_draw->SetMinimum(y_min[j]);

            h_draw->Draw("PE");
            h_my_phi7eta6[j]->Draw("PE SAME");
            h_my_phi5eta4[j]->Draw("PE SAME");
            h_my_phi3eta2[j]->Draw("PE SAME");

            TLegend *leg = new TLegend(0.4, 0.15, 0.9, 0.35);
            leg->SetBorderSize(0);
            leg->AddEntry(h_draw, "Kehao #Delta#eta < 0.02, #Delta#phi < 0.03 (remove pairs)", "pe");
            leg->AddEntry(h_my_phi7eta6[j], "#Delta#eta < 0.02, #Delta#phi < 0.02 (remove tracks)", "pe");
            leg->AddEntry(h_my_phi5eta4[j], "#Delta#eta < 0.04, #Delta#phi < 0.04 (remove tracks) ", "pe");
            leg->AddEntry(h_my_phi3eta2[j], "#Delta#eta < 0.06, #Delta#phi < 0.06 (remove tracks) ", "pe");
            leg->Draw("same");

            TLine *line = new TLine(0.0, 1.0, x_max[i], 1.0);
            line->SetLineStyle(2); 
            line->SetLineColor(kBlack);
            line->SetLineWidth(2);
            line->Draw("same");
        }
    }
    //c1_ol->SaveAs("check_Kehao_70per_cut2.pdf");
}

