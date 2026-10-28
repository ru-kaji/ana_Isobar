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

void Overlay_isobar_Kehao(){

    TFile *f_ref_Para = new TFile("../ana_KFParticle/angel0_90.root", "READ");
    TFile *f_ref_Anti = new TFile("../ana_KFParticle/angel90_180.root", "READ");
    TFile *f_ref_Para_smallbin = new TFile("../ana_KFParticle/smallerbin0_90.root", "READ");
    TFile *f_ref_Anti_smallbin = new TFile("../ana_KFParticle/smallerbin90_180.root", "READ");

    TFile *f_topo_phi3eta2square   = new TFile("check_isobar_KehaoMergeTopocut_cuttracks_R1p4.root", "READ"); // Δη<0.02, Δφ*<0.03　長方形
    TFile *f_my_phi3eta2           = new TFile("check_isobar_cutphi3eta2_cuttracks_R1p4.root", "READ"); // Δη<0.02, Δφ*<0.03 楕円
    TFile *f_my_phi5eta4           = new TFile("check_isobar_cutphi5eta4_cuttracks_R1p4.root", "READ"); // Δη<0.04, Δφ*<0.05 楕円

   

    if (!f_ref_Para || f_ref_Para->IsZombie() || !f_ref_Anti || f_ref_Anti_smallbin->IsZombie()) {
        std::cout << "!!! Error: Cannot open one or more input Kehao files !!!" << std::endl;
        return;
    }
/*
    if(!f_my_S || f_my_S->IsZombie() || !f_my_R0p04 || f_my_R0p04->IsZombie() || !f_my_R0p06 || f_my_R0p06->IsZombie()) {
        std::cout << "!!! Error: Cannot open my files !!!" << std::endl;
        return;
    }*/

    TH1D *h_ref[2]; // 0: Para, 1: Anti (Qinv* スケール, Max 5000)
    TH1D *h_ref_smallbin[2];
    TH1D *h_Kehao_topo[2]; 
    TH1D *h_my_phi3eta2[2]; 
    TH1D *h_my_phi5eta4[2]; 

    h_ref[0] = (TH1D*)f_ref_Para->Get("hCF_cut2_s_cent0");
    h_ref[1] = (TH1D*)f_ref_Anti->Get("hCF_cut2_s_cent0");
    h_ref_smallbin[0] = (TH1D*)f_ref_Para_smallbin->Get("hCF_cut2_s_cent0");
    h_ref_smallbin[1] = (TH1D*)f_ref_Anti_smallbin->Get("hCF_cut2_s_cent0");
    h_Kehao_topo[0]  = (TH1D*)f_topo_phi3eta2square->Get("hCF_Para_LL+AA");
    h_Kehao_topo[1]  = (TH1D*)f_topo_phi3eta2square ->Get("hCF_Anti_LL+AA");
    h_my_phi3eta2[0] = (TH1D*)f_my_phi3eta2->Get("hCF_Para_LL+AA");
    h_my_phi3eta2[1] = (TH1D*)f_my_phi3eta2->Get("hCF_Anti_LL+AA");
    h_my_phi5eta4[0] = (TH1D*)f_my_phi5eta4->Get("hCF_Para_LL+AA");
    h_my_phi5eta4[1] = (TH1D*)f_my_phi5eta4->Get("hCF_Anti_LL+AA");

    for (int j = 0; j < 2; j++) {
        if (!h_ref[j] || !h_ref_smallbin[j] || !h_Kehao_topo[j]) {
            std::cout << "!!! Error: Histogram not found for index " << j << " !!!" << std::endl;
            continue;
        }
    }

    double scale = 1000.0; 
    for (int j = 0; j < 2; j++) {
        double ref_xmin = h_ref_smallbin[j]->GetXaxis()->GetXmin();
        double ref_xmax = h_ref_smallbin[j]->GetXaxis()->GetXmax();
        // 軸の限界値をスケーリング
        h_ref_smallbin[j]->GetXaxis()->SetLimits(ref_xmin / scale, ref_xmax / scale);
    }


    // Kehao (20: 塗りつぶし丸)
    for (int j = 0; j < 2; j++) {
        h_ref_smallbin[j]->SetMarkerStyle(24);
        h_ref_smallbin[j]->SetMarkerSize(0.8);
        h_ref_smallbin[j]->SetMarkerColor(kMagenta);
        h_ref_smallbin[j]->SetLineColor(kMagenta);

        h_Kehao_topo[j]->SetMarkerStyle(20);
        h_Kehao_topo[j]->SetMarkerSize(0.7);
        h_Kehao_topo[j]->SetMarkerColor(kRed);
        h_Kehao_topo[j]->SetLineColor(kRed);

        h_my_phi3eta2[j]->SetMarkerStyle(20);
        h_my_phi3eta2[j]->SetMarkerSize(0.7);
        h_my_phi3eta2[j]->SetMarkerColor(kGreen+2);
        h_my_phi3eta2[j]->SetLineColor(kGreen+2);

        h_my_phi5eta4[j]->SetMarkerStyle(20);
        h_my_phi5eta4[j]->SetMarkerSize(0.7);
        h_my_phi5eta4[j]->SetMarkerColor(kBlue);
        h_my_phi5eta4[j]->SetLineColor(kBlue);
    }


    const char* titles[2] = {"Parallel", "Anti-Parallel"};
    double y_min[2] = {0.6, 0.8}; 
    double x_max[2] = {2.0, 0.6};

    for(int i=0; i<2; i++){
        TCanvas *c1_ol = new TCanvas(Form("c1_ol[%d]", i), "(0.0, 2.0)GeV/c", 1600, 500);
        c1_ol->Divide(2, 1);

        for (int j = 0; j < 2; j++) {
            c1_ol->cd(j + 1);

            TH1D *h_draw = (TH1D*)h_ref_smallbin[j]->Clone(Form("h_ref_clone_%d_%d", i, j));
            h_draw->SetTitle(Form("%s; Q^{inv} (GeV/c); C(Q_{inv})", titles[j]));

            h_draw->GetXaxis()->SetRangeUser(0.0, x_max[i]); 
            h_draw->SetMaximum(1.1);
            h_draw->SetMinimum(y_min[j]);

            h_draw->Draw("PE");
            h_Kehao_topo[j]->Draw("PE SAME");
            h_my_phi3eta2[j]->Draw("PE SAME");
            h_my_phi5eta4[j]->Draw("PE SAME");

            TLegend *leg = new TLegend(0.4, 0.15, 0.8, 0.45);
            leg->SetBorderSize(0);
            leg->AddEntry(h_draw, "Kehao #Delta#eta < 0.02, #Delta#phi^{*} < 0.03 (remove pairs)", "pe");
            leg->AddEntry(h_Kehao_topo[j], "#Delta#eta < 0.02, #Delta#phi^{*} < 0.03 (remove tracks)", "pe");
            leg->AddEntry(h_my_phi3eta2[j], "#frac{#Delta#eta^{2}}{0.02^{2}} + #frac{#Delta#phi^{*2}}{0.03^{2}} < 1.0 (remove tracks)", "pe");
            leg->AddEntry(h_my_phi5eta4[j], "#frac{#Delta#eta^{2}}{0.04^{2}} + #frac{#Delta#phi^{*2}}{0.05^{2}} < 1.0 (remove tracks)", "pe");
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

