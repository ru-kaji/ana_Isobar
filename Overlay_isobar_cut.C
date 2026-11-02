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

    TH1D *h_ref[2]; // 0: Para, 1: Anti (Qinv* スケール, Max 5000)
    TH1D *h_ref_smallbin[2];

    h_ref[0] = (TH1D*)f_ref_Para->Get("hCF_cut2_s_cent0");
    h_ref[1] = (TH1D*)f_ref_Anti->Get("hCF_cut2_s_cent0");
    h_ref_smallbin[0] = (TH1D*)f_ref_Para_smallbin->Get("hCF_cut2_s_cent0");
    h_ref_smallbin[1] = (TH1D*)f_ref_Anti_smallbin->Get("hCF_cut2_s_cent0");

    if (!f_ref_Para || f_ref_Para->IsZombie() || !f_ref_Anti || f_ref_Anti_smallbin->IsZombie()) {
        std::cout << "!!! Error: Cannot open one or more input Kehao files !!!" << std::endl;
        return;
    }

    for (int j = 0; j < 2; j++) {
        if (!h_ref[j] || !h_ref_smallbin[j] ) {
            std::cout << "!!! Error: Kehao Histogram not found for index " << j << " !!!" << std::endl;
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

    for (int j = 0; j < 2; j++) {
        h_ref_smallbin[j]->SetMarkerStyle(24);
        h_ref_smallbin[j]->SetMarkerSize(0.8);
        h_ref_smallbin[j]->SetMarkerColor(kMagenta);
        h_ref_smallbin[j]->SetLineColor(kMagenta);
    }

    // -------------------- f1 ------------------------- 
    TFile *f1    = new TFile("check_isobar_KehaoMergeTopocut_cuttracks_R1p4.root", "READ"); // Δη<0.02, Δφ*<0.02

    TH1D *h1[2] = {nullptr, nullptr};  // 0: Para, 1: Anti (Qinv スケール, Max 0.5)
    h1[0]  = (TH1D*)f1->Get("hCF_Para_LL+AA");
    h1[1]  = (TH1D*)f1->Get("hCF_Anti_LL+AA");

    // -------------------- f2 ------------------------- 
    TFile *f2   = new TFile("check_isobar_KehaoMergeTopocut_cuttracks_R1p4_highpT.root", "READ"); // Δη<0.04, Δφ*<0.04

    TH1D *h2[2] = {nullptr, nullptr};
    h2[0]  = (TH1D*)f2->Get("hCF_Para_LL+AA");
    h2[1]  = (TH1D*)f2->Get("hCF_Anti_LL+AA");

    // -------------------- f3 ------------------------- 
    TFile *f3   = new TFile("check_isobar_cutR0p06.root", "READ"); // Δη<0.06, Δφ*<0.06

    TH1D *h3[2] = {nullptr, nullptr};
    h3[0]  = (TH1D*)f3->Get("hCF_Para_LL+AA");
    h3[1]  = (TH1D*)f3->Get("hCF_Anti_LL+AA");


    for (int j = 0; j < 2; j++) {

        h1[j]->SetMarkerStyle(20);
        h1[j]->SetMarkerSize(0.7);
        h1[j]->SetMarkerColor(kRed);
        h1[j]->SetLineColor(kRed);

        h2[j]->SetMarkerStyle(20);
        h2[j]->SetMarkerSize(0.7);
        h2[j]->SetMarkerColor(kGreen+2);
        h2[j]->SetLineColor(kGreen+2);

        h3[j]->SetMarkerStyle(20);
        h3[j]->SetMarkerSize(0.7);
        h3[j]->SetMarkerColor(kBlue);
        h3[j]->SetLineColor(kBlue);
    }


    const char* titles[2] = {"Para", "Anti"};
    double y_min[2] = {0.6, 0.8}; 
    double x_max[2] = {2.0, 0.6};

    for(int i=0; i<2; i++){
        TCanvas *c1_ol = new TCanvas(Form("c1_ol[%d]", i), "(0.0, 2.0)GeV/c", 1600, 500);
        c1_ol->Divide(2, 1);

        for (int j = 0; j < 2; j++) {
            c1_ol->cd(j + 1);

            TH1D *h_ref = (TH1D*)h_ref_smallbin[j]->Clone(Form("h_ref_clone_%d_%d", i, j));
            h_ref->SetTitle(Form("CF %s; Q^{inv} (GeV/c); C(Q_{inv})", titles[j]));

            h_ref->GetXaxis()->SetRangeUser(0.0, x_max[i]); 
            h_ref->SetMaximum(1.1);
            h_ref->SetMinimum(y_min[j]);

            h_ref->Draw("PE");
            if(h1[j]) h1[j]->Draw("PE SAME");
            if(h2[j]) h2[j]->Draw("PE SAME");
            //if(h3[j]) h3[j]->Draw("PE SAME");

            TLegend *leg = new TLegend(0.4, 0.15, 0.9, 0.35);
            leg->SetBorderSize(0);
            leg->AddEntry(h_ref, "Kehao (0.4< p_{T} <3.0)(remove pairs)", "pe");
            if(h1[j]) leg->AddEntry(h1[j], "Rurika (0.4< p_{T} <3.0)(remove tracks)", "pe");
            if(h2[j]) leg->AddEntry(h2[j], "Rurika (0.4< p_{T})(remove tracks)", "pe");
            //if(h3[j]) leg->AddEntry(h3[j], "#frac{#Delta#eta^{2}}{0.02^{2}} + #frac{#Delta#phi^{*2}}{0.03^{2}} < 1.0 ", "pe");
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

