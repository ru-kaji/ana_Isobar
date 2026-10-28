#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "TBox.h"
#include "AxisStyle.h"
#include <iostream>
using namespace std;

TLegend* SetLegend(TH1* h1, TH1* h2);

TH1* ScaleHist(TH1* h_Same, TH1* h_Mix, double xmin, double xmax){

    int b_start = h_Same->FindBin(xmin);
    int b_end   = h_Same->FindBin(xmax);

	TString newName = TString(h_Same->GetName()) + "_scaled_" + TString(h_Mix->GetName());
	TH1* h_scaled = (TH1*)h_Mix->Clone(newName);

	if (h_Mix->Integral(b_start, b_end) > 0) {
		double scale = h_Same->Integral(b_start, b_end) / h_Mix->Integral(b_start, b_end);
		h_scaled->Scale(scale);
	}

	return h_scaled;
}

TH1* DivideHist(TH1* h_Same, TH1* h_Mix_scaled){

	h_Same->Sumw2();
	h_Mix_scaled->Sumw2();
    
	TString newName = TString(h_Same->GetName()) + "_ratio";
	TH1* h_ratio = (TH1*)h_Same->Clone(newName);

	if (h_Mix_scaled->Integral() > 0) h_ratio->Divide(h_Mix_scaled);

	return h_ratio;
}

double min_Angle = -180.0;
double max_Angle = 180.0;
double min_cosAngle = -1.0;
double max_cosAngle = 1.0;
double min_DAngle = -1.0;
double max_DAngle = 1.0;
double min_Mass = 2.5;
double max_Mass = 3.0;
double min_Qinv = 0.7;
double max_Qinv = 1.0;


void DrawCanvas_isobar(TString fileName){ 

	std::cout << "DEBUG: DrawCanvas.C is loading file: " << fileName << std::endl;

    TFile *fin = new TFile(fileName, "READ");
	if (!fin || fin->IsZombie()) {
        std::cout << "!!! Error: Cannot open file .root !!!" << std::endl;
        return;
    }

	TH1D *hcosAngle_same[3], *hAngle_same[3], *hcosAngle_mix[3], *hAngle_mix[3];
	TH1D *hmass_same[5][3][4][7], *hmass_mix[5][3][4][7];
	TH1D *hqinv_same[5][3][4][7][9], *hqinv_mix[5][3][4][7][9];
	TH1D *hmass_same_low[3][4][7], *hmass_mix_low[3][4][7];
	TH1D *hmass_mix_weigh[3][4], *hqinv_mix_weigh[3][4];

	const char* pidName[3]  = {"LL", "AA", "LA"};
    const char* pidTitle[3] = {"#Lambda#Lambda", "#bar{#Lambda}#bar{#Lambda}", "#Lambda#bar{#Lambda}"};
    
    const char* pairName[4]  = {"SS", "SB", "BS", "BB"};
    const char* pairTitle[4] = {"SS", "SB", "BS", "BB"};
    
    const char* spinName[3]  = {"All", "Para", "Anti"};
	const char* spinName7[7] = {"All", "0-90", "90-180", "0-45", "45-90", "90-135", "135-180"};
    const char* spinTitle[3] = {"All-Spin", "Parallel", "Anti-Parallel"};
	const char* spinTitle7[7] = {"All-Spin", "0-90", "90-180", "0-45", "45-90", "90-135", "135-180"};

	const char* centName[9] = {"0-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "All"};
	const char* centTitle[9] = {"0-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "All"};

	for(int ipid = 0; ipid < 3; ipid++){

		cout << "Get file, ipid=" << ipid << endl;

		hcosAngle_same[ipid] = (TH1D*)fin->Get(Form("hcosAngle_same_%s", pidName[ipid]));
		hAngle_same[ipid]    = (TH1D*)fin->Get(Form("hAngle_same_%s", pidName[ipid]));
		hcosAngle_mix[ipid]  = (TH1D*)fin->Get(Form("hcosAngle_mix_%s", pidName[ipid]));
		hAngle_mix[ipid]     = (TH1D*)fin->Get(Form("hAngle_mix_%s", pidName[ipid]));
		if (!hcosAngle_same[ipid]) std::cout << "!!! Missing: hcosAngle_same_" << pidName[ipid] << " !!!" << std::endl;

		for(int ipair = 0; ipair < 4; ipair++) {
			for(int ispin = 0; ispin < 3; ispin++) {
				for(int icut = 0; icut <5; icut++){
					TString nameMassS = Form("hmass_same_%d_%s_%s_%s", icut, pidName[ipid], pairName[ipair], spinName7[ispin]);
					TString nameQinvS = Form("hqinv_same_%d_%s_%s_%s_%s", icut, pidName[ipid], pairName[ipair], spinName7[ispin], centName[8]);
					TString nameMassM = Form("hmass_mix_%d_%s_%s_%s", icut, pidName[ipid], pairName[ipair], spinName7[ispin]);
					TString nameQinvM = Form("hqinv_mix_%d_%s_%s_%s_%s", icut, pidName[ipid], pairName[ipair], spinName7[ispin], centName[8]);
					hmass_same[icut][ipid][ipair][ispin]      = (TH1D*)fin->Get(nameMassS);
					hqinv_same[icut][ipid][ipair][ispin][8]   = (TH1D*)fin->Get(nameQinvS);
					hmass_mix[icut][ipid][ipair][ispin]       = (TH1D*)fin->Get(nameMassM);
					hqinv_mix[icut][ipid][ipair][ispin][8]    = (TH1D*)fin->Get(nameQinvM);
					if (!hmass_same[icut][ipid][ipair][ispin]) std::cout << "!!! Missing: hmass_same" << pidName[ipid] << " !!!" << std::endl;
					if (!hqinv_mix[icut][ipid][ipair][ispin][8])  std::cout << "!!! Missing: hqinv_mix" << pidName[ipid] << " !!!" << std::endl;
				}
			
				TString nameMassS_low = Form("hmass_same_low_%s_%s_%s", pidName[ipid], pairName[ipair], spinName7[ispin]);
				TString nameMassM_low = Form("hmass_mix_low_%s_%s_%s",  pidName[ipid], pairName[ipair], spinName7[ispin]);
				hmass_same_low[ipid][ipair][ispin] = (TH1D*)fin->Get(nameMassS_low);
				hmass_mix_low[ipid][ipair][ispin]  = (TH1D*)fin->Get(nameMassM_low);
			} // end of spin loop

			if( ipid == 0 || ipid == 1) {
				hmass_mix_weigh[ipid][ipair] = (TH1D*)fin->Get(Form("hmass_mix_weigh_%s_%s", pidName[ipid], pairName[ipair]));
				hqinv_mix_weigh[ipid][ipair] = (TH1D*)fin->Get(Form("hqinv_mix_weigh_%s_%s", pidName[ipid], pairName[ipair]));
				if (!hmass_mix_weigh[ipid][ipair]) std::cout << "!!! Missing: hmass_mix_weigh" << pidName[ipid] << " !!!" << std::endl;
				if (!hqinv_mix_weigh[ipid][ipair])  std::cout << "!!! Missing: hqinv_mix_weigh" << pidName[ipid] << " !!!" << std::endl;
			}
		} // end of pair loop
	} // end of pid loop
		
	TCanvas *c3[3][3], *c4[3], *c5[3], *c6[3], *c7, *c8[3], *c9[3], *c10, *c11, *c_weigh;
	
	gStyle->SetPadBottomMargin(0.14);
	gStyle->SetPadLeftMargin(0.15);
	

	//--------------------------------------------------
	//  c1  : Angle   		 (particle) = (LL, AA, LA) \ (Same/Mix)
	//  c2  : cosAngle		 (particle) = (LL, AA, LA) \ (Same/Mix)
	//  c6  : Mass	 		 (same, mix, same/mix) \ (SS, SB, BB)
	//  c7	: Qinv	 		 (same, mix, same/mix) \ (SS, SB, BB)
	//  c8  : Low Qinv		 (same, mix, same/mix) \ (All, Para, Anti)
	//  c9  : Mass		 	 (same, mix, same/mix) \ (All, Para, Anti)
	//  c10 : Qinv			 (same, mix, same/mix) \ (All, Para, Anti)
	//--------------------------------------------------

		/*
		//--------------------------------------------------
		//  c1: Angle //SSのみ
		//--------------------------------------------------
		TCanvas *c1 = new TCanvas("c1", "Histgram1", 1200, 700);
		c1->Divide(3,2);

		for(int ipid=0; ipid<3; ipid++){ 

			if (hAngle_same[ipid] && hAngle_mix[ipid]) {
				c1->cd(ipid+1);
				hAngle_same[ipid]->SetTitle(Form("Same %s", pidName[ipid]));
				hAngle_same[ipid]->SetXTitle("Angle [degree]");
				hAngle_same[ipid]->Draw();
				SetAxis(hAngle_same[ipid]);

				c1->cd(ipid+4);
			TH1 *hAngle_scaled = ScaleHist(hAngle_same[ipid], hAngle_mix[ipid], min_Angle, max_Angle);
			TH1 *hAngle_ratio  = DivideHist(hAngle_same[ipid], hAngle_scaled);
			hAngle_ratio->SetTitle(Form("(Same/Mix) %s", pidName[ipid]));
			hAngle_ratio->SetXTitle("Angle [degree]");
			hAngle_ratio->Draw();
			SetAxis(hAngle_ratio);
		}
		}

		//--------------------------------------------------
		//  c2: cosAngle
		//--------------------------------------------------
		TCanvas *c2 = new TCanvas("c2", "Histgram2", 1200, 700);
		c2->Divide(3,2);

		for(int ipid=0; ipid<3; ipid++){

			if (hcosAngle_same[ipid] && hcosAngle_mix[ipid]) {

			c2->cd(ipid+1);
			hcosAngle_same[ipid]->SetTitle(Form("Same %s", pidName[ipid]));
			hcosAngle_same[ipid]->SetXTitle("cos#theta_{12}");
			hcosAngle_same[ipid]->Draw();
			SetAxis(hcosAngle_same[ipid]);

			c2->cd(ipid+4);
			TH1 *hcosAngle_scaled = ScaleHist(hcosAngle_same[ipid], hcosAngle_mix[ipid], min_cosAngle, max_cosAngle);
			TH1 *hcosAngle_ratio  = DivideHist(hcosAngle_same[ipid], hcosAngle_scaled);
			hcosAngle_ratio->SetTitle(Form("(Same/Mix) %s", pidName[ipid]));
			hcosAngle_ratio->SetXTitle("cos#theta_{12}");
			hcosAngle_ratio->Draw();
			SetAxis(hcosAngle_ratio);
			}
		}*/

	//===========================================================================

/*
		//--------------------------------------------------
		//  c6: M_inv (pair) (All spin)
		//--------------------------------------------------
		c6[ipid] = new TCanvas(Form("c6_%s", pidName[ipid]), Form("M_inv_%s", pidName[ipid]), 1200, 600);
		c6[ipid]->Divide(3, 2);

		for(int ipair=0; ipair<4; ipair++){ // 0:SS, 1:SB, 2:BS, 3:BB

			if (!hmass_same[4][ipid][ipair][0] || !hmass_mix[4][ipid][ipair][0]) {
				std::cout << "Skip mass plot (c6) for ipid=" << ipid << " ipair=" << ipair << " (missing hist)" << std::endl;
				continue;
			}

			c6[ipid]->cd(ipair+1);
			//hmass_same[ipid][ipair][0]->Draw();
			TH1 *hmass_mix_scaled = ScaleHist(hmass_same[4][ipid][ipair][0], hmass_mix[4][ipid][ipair][0], min_Mass, max_Mass);
			hmass_mix_scaled->SetTitle(Form("Same %s %s", pidName[ipid], pairName[ipair]));
			//hmass_mix_scaled->Draw("same");
			SetRed(hmass_same[4][ipid][ipair][0]);
			SetBlack(hmass_mix_scaled);
			SetAxis(hmass_same[4][ipid][ipair][0]);
			SetLegend(hmass_same[4][ipid][ipair][0], hmass_mix_scaled)->Draw();

			c6[ipid]->cd(ipair+4);
			TH1 *hmass_diff = (TH1*)hmass_same[4][ipid][ipair][0]->Clone(Form("hmass_diff_%s_%s", pidName[ipid], pairName[ipair]));
			hmass_diff->Add(hmass_mix_scaled, -1.0);
			hmass_diff->SetTitle(Form("(Same-Mix) %s %s", pidName[ipid], pairName[ipair]));
			hmass_diff->SetXTitle("M_{inv}");
			hmass_diff->SetYTitle("Counts");
			//hmass_diff->Draw();
			SetAxis(hmass_diff);
			SetRed(hmass_diff);
			
		}
	} // end of pid loop

		//--------------------------------------------------
		//  c7: Qinv (All spin) ← ParaとAntiどうしよう
		//---------------------------------------------------	
		c7 = new TCanvas("c7", "Histogram7", 1000, 800);
		c7->Divide(3, 3);

		for(int ipid=0; ipid<3; ipid++){ // 0:LL, 1:AA, 2:LA
			for(int ipair=0; ipair<3; ipair++){ // 0:SS, 1:SB&BS, 2:BB

				if (!hqinv_same[0][ipid][ipair][0][8] || !hqinv_mix[0][ipid][ipair][0][8]) {
					std::cout << "Skip Qinv plot (c7) for ipid=" << ipid << " ipair=" << ipair << " (missing hist)" << std::endl;
					continue;
				}

				c7->cd(ipair+ipid*3+1);
				hqinv_same[0][ipid][ipair][0][8]->Draw();
				TH1 *hqinv_mix_scaled = ScaleHist(hqinv_same[0][ipid][ipair][0][8], hqinv_mix[0][ipid][ipair][0][8], min_Qinv, max_Qinv);
				hqinv_mix_scaled->SetTitle(Form("%s %s", pidName[ipid], pairName[ipair]));
				hqinv_mix_scaled->Draw("same");
				
				SetAxis(hqinv_same[ipid][ipair][0][8]);
				SetLegend(hqinv_same[ipid][ipair][0][8], hqinv_mix_scaled)->Draw();
			}
		}
		*/
/*
		//---------------------------------------------------
		//  c8: Mass (Low Qinv)
		//---------------------------------------------------	
	for(int ipid=0; ipid<3; ipid++){ // 0:LL, 1:AA, 2:LA
		c8[ipid] = new TCanvas(Form("c8_%s", pidName[ipid]), "Histogram8", 1400, 700);
		c8[ipid]->Divide(3, 2);
	
		for(int ipair=0; ipair<3; ipair++){ // 0:SS, 1:SB&BS, 2:BB

			if (!hmass_same_low[ipid][ipair][0] || !hmass_mix_low[ipid][ipair][0]) {
				std::cout << "Skip low-mass plot (c8) for ipid=" << ipid << " ipair=" << ipair << " (missing hist)" << std::endl;
				continue;
			}

			c8[ipid]->cd(ipair+1);
			hmass_same_low[ipid][ipair][0]->Draw();
			TH1 *hmass_mix_low_scaled = ScaleHist(hmass_same_low[ipid][ipair][0], hmass_mix_low[ipid][ipair][0], 2.4, 2.6);
			hmass_mix_low_scaled->Draw("same");
			SetRed(hmass_same_low[ipid][ipair][0]);
			SetBlack(hmass_mix_low_scaled);
			SetAxis(hmass_same_low[ipid][ipair][0]);
			SetLegend(hmass_same_low[ipid][ipair][0], hmass_mix_low_scaled)->Draw();

			c8[ipid]->cd(ipair+4);
			TH1 *hmass_low_diff = (TH1*)hmass_same_low[ipid][ipair][0]->Clone("hmass_low_diff");
			hmass_low_diff->Add(hmass_mix_low_scaled, -1.0);
			hmass_low_diff->SetTitle(Form("(Same-Mix) %s", pidName[ipid]));
			hmass_low_diff->Draw();
			SetAxis(hmass_low_diff);
			
		}


		//--------------------------------------------------
		//  c9: Mass (Spin) (SS)
		//---------------------------------------------------	
		c9[ipid] = new TCanvas(Form("c9_%s", pidName[ipid]), "Histogram9", 1600, 700);
		c9[ipid]->Divide(4, 2); 

		for(int ispin=0; ispin<3; ispin++){ // 0:All, 1:Para, 2:AntiPara
			if (!hmass_same[ipid][0][ispin] || !hmass_mix[ipid][0][ispin]) {
				std::cout << "Skip mass spin plot (c9) for ipid=" << ipid << " ispin=" << ispin << " (missing hist)" << std::endl;
				continue;
			}
			c9[ipid]->cd(ispin+1);
			hmass_same[ipid][0][ispin]->Draw();
			TH1 *hmass_mix_spin_scaled = ScaleHist(hmass_same[ipid][0][ispin], hmass_mix[ipid][0][ispin], min_Mass, max_Mass);
			hmass_mix_spin_scaled->Draw("same");
			SetRed(hmass_same[ipid][0][ispin]);
			SetBlack(hmass_mix_spin_scaled);
			SetAxis(hmass_same[ipid][0][ispin]);
			SetLegend(hmass_same[ipid][0][ispin], hmass_mix_spin_scaled)->Draw();

			c9[ipid]->cd(ispin+5);
			TH1 *hmass_spin_diff = (TH1*)hmass_same[ipid][0][ispin]->Clone("hmass_spin_diff");
			hmass_spin_diff->Add(hmass_mix_spin_scaled, -1.0);
			hmass_spin_diff->SetTitle(Form("(Same-Mix) %s", spinName[ispin]));
			hmass_spin_diff->Draw();
			SetAxis(hmass_spin_diff); 	
			SetRed(hmass_spin_diff); 		
			
			double xmin = hmass_spin_diff->GetXaxis()->GetXmin();
			double xmax = hmass_spin_diff->GetXaxis()->GetXmax();
			TLine *line = new TLine(xmin, 0.0, xmax, 0.0);
			line->SetLineStyle(2); // dot line
			line->SetLineColor(kBlack);
			line->SetLineWidth(2);
			line->Draw("same");
		}

		c9[ipid]->cd(4);
		hmass_same[4][ipid][0][1]->Draw();
		hmass_same[4][ipid][0][2]->Draw("same");
		SetRed(hmass_same[4][ipid][0][1]);
		SetBlue(hmass_same[4][ipid][0][2]);
		SetAxis(hmass_same[4][ipid][0][1]);
		SetLegend_spin(hmass_same[4][ipid][0][1], hmass_same[4][ipid][0][2])->Draw();

		c9[ipid]->cd(8);
		TH1 *hmass_Para_diff = (TH1*)hmass_same[ipid][0][1]->Clone("hmass_Para_diff");
		TH1 *hmass_mix_Para_scaled = ScaleHist(hmass_same[ipid][0][1], hmass_mix[ipid][0][1], min_Mass, max_Mass);
		hmass_Para_diff->Add(hmass_mix_Para_scaled, -1.0);
		hmass_Para_diff->SetTitle(Form("(Same-Mix) %s", pidName[ipid]));
		hmass_Para_diff->Draw();

		TH1 *hmass_Anti_diff = (TH1*)hmass_same[ipid][0][2]->Clone("hmass_Anti_diff");
		TH1 *hmass_mix_Anti_scaled = ScaleHist(hmass_same[ipid][0][2], hmass_mix[ipid][0][2], min_Mass, max_Mass);
		hmass_Anti_diff->Add(hmass_mix_Anti_scaled, -1.0);
		hmass_Anti_diff->Draw("same");
		SetAxis(hmass_Para_diff);
		SetAxis(hmass_Anti_diff);
		SetRed(hmass_Para_diff);
		SetBlue(hmass_Anti_diff);

		double xmin = hmass_Anti_diff->GetXaxis()->GetXmin();
		double xmax = hmass_Anti_diff->GetXaxis()->GetXmax();
		TLine *line = new TLine(xmin, 0.0, xmax, 0.0);
		line->SetLineStyle(2); // dot line
		line->SetLineColor(kBlack);
		line->SetLineWidth(2);
		line->Draw("same");

	}// end of pid loop
		//--------------------------------------------------
		//  c10: Qinv (Spin) (SS)
		//---------------------------------------------------	
		c10 = new TCanvas("c10", "Histogram10", 1300, 800);
		c10->Divide(4, 3); 

		for(int ipid=0; ipid<3; ipid++){ // 0:LL, 1:AA, 2:LA

			for(int ispin=0; ispin<3; ispin++){ // 0:All, 1:Para, 2:AntiPara
				if (!hqinv_same[ipid][0][ispin][8] || !hqinv_mix[ipid][0][ispin][8]) {
					std::cout << "Skip Qinv spin plot (c10) for ipid=" << ipid << " ispin=" << ispin << " (missing hist)" << std::endl;
					continue;
				}
				c10->cd(ispin+ipid*4+1);
				hqinv_same[ipid][0][ispin][8]->Draw();
				TH1 *hqinv_mix_spin_scaled = ScaleHist(hqinv_same[ipid][0][ispin][8], hqinv_mix[ipid][0][ispin][8], min_Qinv, max_Qinv);
				hqinv_mix_spin_scaled->SetTitle(Form("%s %s", pidName[ipid], spinName[ispin]));
				hqinv_mix_spin_scaled->Draw("same");
				SetRed(hqinv_same[ipid][0][ispin][8]);
				SetBlack(hqinv_mix_spin_scaled);
				SetAxis(hqinv_same[ipid][0][ispin][8]);
				SetLegend(hqinv_same[ipid][0][ispin][8], hqinv_mix_spin_scaled)->Draw();
			}

			c10->cd(ipid*4+4); // Overlay spin
			hqinv_same[ipid][0][1][8]->Draw();
			hqinv_same[ipid][0][2][8]->SetTitle(Form("%s", pidName[ipid]));
			hqinv_same[ipid][0][2][8]->Draw("same");
			SetRed(hqinv_same[ipid][0][1][8]);
			SetBlue(hqinv_same[ipid][0][2][8]);
			SetAxis(hqinv_same[ipid][0][1][8]);
			SetLegend_spin(hqinv_same[ipid][0][1][8], hqinv_same[ipid][0][2][8])->Draw();
		}

		//--------------------------------------------------
		//  c11: Mass (Spin) (SS)
		//---------------------------------------------------	
		c11 = new TCanvas("c11", "Histogram11", 900, 700);
		c11->Divide(2, 2); 

		TH1 *hmass_mix_scaled_LA[3], *hmass_ratio_LA[3];
		TH1 *hmass_mix_scaled_SS[2][3], *hmass_ratio_SS[2][3];
		TH1 *hmass_diff_LA[2][3];
		

		for(int ispin=0; ispin<3; ispin++){
				hmass_mix_scaled_LA[ispin] = ScaleHist(hmass_same[2][0][ispin], hmass_mix[2][0][ispin], 3.5, 4.0);
			hmass_ratio_LA[ispin] = DivideHist(hmass_same[2][0][ispin], hmass_mix_scaled_LA[ispin]);

			for(int ipid=0; ipid<2; ipid++){
				hmass_mix_scaled_SS[ipid][ispin] = ScaleHist(hmass_same[ipid][0][ispin], hmass_mix[ipid][0][ispin], 3.5, 4.0);
				hmass_ratio_SS[ipid][ispin] = DivideHist(hmass_same[ipid][0][ispin], hmass_mix_scaled_SS[ipid][ispin]);
			}
				hmass_diff_LA[0][ispin] = (TH1*)hmass_ratio_SS[0][ispin]->Clone(Form("hmass_ratio_LA_%s",spinName[ispin]));
				hmass_diff_LA[0][ispin]->Add(hmass_ratio_LA[ispin], -1.0);
				hmass_diff_LA[1][ispin] = (TH1*)hmass_ratio_SS[1][ispin]->Clone(Form("hmass_ratio_LA_%s",spinName[ispin]));
				hmass_diff_LA[1][ispin]->Add(hmass_ratio_LA[ispin], -1.0);
			
		}

		for(int ipid=0; ipid<2; ipid++){ 
			
			c11->cd(ipid+1);
			hmass_diff_LA[ipid][0]->SetTitle(Form("%s-LA; M_{inv}(GeV/c); (Same/Mix) (%s-LA)",pidName[ipid], pidName[ipid]));
			hmass_diff_LA[ipid][0]->SetLineColor(kBlack);
			hmass_diff_LA[ipid][0]->Draw();

			double xmin1 = hmass_diff_LA[ipid][0]->GetXaxis()->GetXmin();
			double xmax1 = hmass_diff_LA[ipid][0]->GetXaxis()->GetXmax();
			TLine *line1 = new TLine(xmin1, 0.0, xmax1, 0.0);
			line1->SetLineStyle(2); // dot line
			line1->SetLineColor(kBlack);
			line1->SetLineWidth(2);
			line1->Draw("same");

			c11->cd(ipid+3);
			hmass_diff_LA[ipid][1]->SetTitle(Form("%s-LA, All-spin; M_{inv}(GeV/c); (Same/Mix) (LL-LA)", pidName[ipid]));
			hmass_diff_LA[ipid][1]->SetLineColor(kRed);
			hmass_diff_LA[ipid][1]->Draw();

			hmass_diff_LA[ipid][2]->SetLineColor(kBlue);
			hmass_diff_LA[ipid][2]->Draw("SAME");

			double xmin2 = hmass_diff_LA[ipid][2]->GetXaxis()->GetXmin();
			double xmax2 = hmass_diff_LA[ipid][2]->GetXaxis()->GetXmax();
			TLine *line2 = new TLine(xmin2, 0.0, xmax2, 0.0);
			line2->SetLineStyle(2); // dot line
			line2->SetLineColor(kBlack);
			line2->SetLineWidth(2);
			line2->Draw("same");
		}*/

	// ====================================================================
	//  c_weigh: Same - Mix(weight)
	// ====================================================================	
	c_weigh = new TCanvas("c_weigh", "Weighted M_inv", 1400, 600);
	c_weigh->Divide(3, 2);

	for(int ipid = 0; ipid < 2; ipid++){ // 0:LL, 1:AA

		TH1D *hSame_for_weigh = hmass_same[0][ipid][0][2];
		TH1D *hMix_noweigh = hmass_mix[0][ipid][0][2];
		TH1D *hMix_for_weigh  = hmass_mix_weigh[ipid][0];
		if(!hSame_for_weigh || !hMix_for_weigh) {
			std::cout << "Skip weighted mass plot for ipid=" << ipid << " (missing hist)" << std::endl;
			continue;
		}

		if (hSame_for_weigh->GetSumw2N() == 0) hSame_for_weigh->Sumw2();
		if (hMix_for_weigh->GetSumw2N() == 0)  hMix_for_weigh->Sumw2();
		
		TH1 *hMixWeigh_scaled = ScaleHist(hSame_for_weigh, hMix_for_weigh, 2.24, 2.2799);
		TH1 *hMixnoWeigh_scaled = ScaleHist(hSame_for_weigh, hMix_noweigh, 2.24, 2.2799);
		c_weigh->cd(ipid*3+1);
		hSame_for_weigh->Draw();
		hMixWeigh_scaled->SetLineColor(kRed);
		hMixWeigh_scaled->Draw("same");

		TH1 *hDiff = (TH1*)hSame_for_weigh->Clone(Form("hmass_diff_weigh_%s", pidName[ipid]));
		hDiff->Add(hMixWeigh_scaled, -1.0);

		TH1 *hDiff_noweigh = (TH1*)hSame_for_weigh->Clone(Form("hmass_diff_noweigh_%s", pidName[ipid]));
		hDiff_noweigh->Add(hMixnoWeigh_scaled, -1.0);

		// 5. 描画
		c_weigh->cd(ipid*3+2);
		hDiff_noweigh->SetTitle(Form("(Same - Mix_{weight}) %s", pidName[ipid]));
		hDiff_noweigh->SetTitle(Form("Difference %s", pidName[ipid]));
		hDiff_noweigh->SetXTitle("M_{inv} (GeV/c^{2})");
		hDiff_noweigh->SetYTitle("Counts");
		hDiff_noweigh->SetMarkerStyle(20);
		hDiff_noweigh->SetMarkerSize(0.8);
		hDiff_noweigh->SetLineColor(kBlue);
		hDiff_noweigh->SetMarkerColor(kBlue);
		hDiff_noweigh->GetYaxis()->SetRangeUser(-500, 500); // 最初に描く方にRangeを設定
		hDiff_noweigh->Draw("PE");

		// 2. 後から赤色 (weight) を重ねる
		hDiff->SetMarkerStyle(24);
		hDiff->SetMarkerSize(0.8);
		hDiff->SetMarkerColor(kRed);
		hDiff->SetLineColor(kRed);
		hDiff->Draw("PE same");


		// ゼロのラインを引く
		TLine *line = new TLine(hDiff->GetXaxis()->GetXmin(), 0.0, hDiff->GetXaxis()->GetXmax(), 0.0);
		line->SetLineStyle(2);
		line->SetLineColor(kBlack);
		line->Draw("same");

		TLegend *leg = new TLegend(0.5, 0.2, 0.9, 0.3);
		leg->SetBorderSize(0);
		leg->SetFillStyle(0);
		leg->AddEntry(hDiff, "Same - Mix_{weight}", "l");
		leg->AddEntry(hDiff_noweigh, "Same - Mix_{no weight}", "l");
		leg->Draw();

		c_weigh->cd(ipid*3+3);
		TH1D *h_test_ratio = (TH1D*)hMix_for_weigh->Clone("h_test_ratio");
		h_test_ratio->Divide(hMix_noweigh);
		h_test_ratio->SetTitle("Mix_weight / Mix_noweight; M_{inv}; Ratio");
		h_test_ratio->GetYaxis()->SetRangeUser(0.0, 2.0);
		h_test_ratio->SetMinimum(0.95);
		h_test_ratio->SetMaximum(1.05);
		h_test_ratio->Draw();

	}
	
} // end of Canvas


TLegend* SetLegend(TH1* h1, TH1* h2){

    TLegend *leg = new TLegend(0.4,0.2,0.7,0.4);
    leg->SetBorderSize(0); // 枠線をなくす
    leg->SetFillStyle(0);  // 背景を透明にする
    leg->SetTextFont(43);
    leg->SetTextSize(19);

    leg->AddEntry(h1,"Same","l");
    leg->AddEntry(h2,"Mix(Norm)","l");

    return leg;
}
