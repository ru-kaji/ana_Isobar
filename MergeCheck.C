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

	TString newName = TString(h_Same->GetName()) + "_scaled";
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

double min_DAngle = -1.0;
double max_DAngle = 1.0;


void MergeCheck(TString fileName){ 

	std::cout << "Detadphi.C is loading file: " << fileName << std::endl;

    TFile *fin = new TFile(fileName, "READ");
	if (!fin || fin->IsZombie()) {
        std::cout << "!!! Error: Cannot open file .root !!!" << std::endl;
        return;
    }

	TH1D *hDeta[3][4], *hDphi[3][4], *hDeta_mix[3][4], *hDphi_mix[3][4]; // [pid][dgpair]
	TH2F *hDetaDphi[3][4], *hDetaDphi_mix[3][4];
	TH2D *h2D_Ang[3][4][3], *h2D_Ang_mix[3][4][3]; //[pid][dgpair][ispin]

	const char* dgPairName[4] = {"pp", "pipi", "ppi", "pip"};
	const char* dgPairName_2D[4] = {"p-p", "#pi-#pi", "p-#pi", "#pi-p"};

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

		for(int idgpair = 0; idgpair < 4; idgpair++){
			hDeta[ipid][idgpair] = (TH1D*)fin->Get(Form("h1D_Deta_%s_%s", pidName[ipid], dgPairName[idgpair]));
			hDphi[ipid][idgpair] = (TH1D*)fin->Get(Form("h1D_Dphi_%s_%s", pidName[ipid], dgPairName[idgpair]));
			hDeta_mix[ipid][idgpair] = (TH1D*)fin->Get(Form("h1D_Deta_mix_%s_%s", pidName[ipid], dgPairName[idgpair]));
			hDphi_mix[ipid][idgpair] = (TH1D*)fin->Get(Form("h1D_Dphi_mix_%s_%s", pidName[ipid], dgPairName[idgpair]));
		}

		for(int ispin=0; ispin<3; ispin++){
			h2D_Ang[ipid][0][ispin]     = (TH2D*)fin->Get(Form("hDetaDphiStar1_%s_%s", pidName[ipid], spinName[ispin]));
			h2D_Ang[ipid][1][ispin]     = (TH2D*)fin->Get(Form("hDetaDphiStar2_%s_%s", pidName[ipid], spinName[ispin]));
			h2D_Ang[ipid][2][ispin]     = (TH2D*)fin->Get(Form("hDetaDphiStar3_%s_%s", pidName[ipid], spinName[ispin]));
			h2D_Ang[ipid][3][ispin]     = (TH2D*)fin->Get(Form("hDetaDphiStar4_%s_%s", pidName[ipid], spinName[ispin]));
			if (!h2D_Ang[ipid][3][ispin]) std::cout << "!!! Missing: h2D_Ang" << pidName[ipid] << " !!!" << std::endl;

			h2D_Ang_mix[ipid][0][ispin] = (TH2D*)fin->Get(Form("hDetaDphiStar1_mix_%s_%s", pidName[ipid], spinName[ispin]));
			h2D_Ang_mix[ipid][1][ispin] = (TH2D*)fin->Get(Form("hDetaDphiStar2_mix_%s_%s", pidName[ipid], spinName[ispin]));
			h2D_Ang_mix[ipid][2][ispin] = (TH2D*)fin->Get(Form("hDetaDphiStar3_mix_%s_%s", pidName[ipid], spinName[ispin]));
			h2D_Ang_mix[ipid][3][ispin] = (TH2D*)fin->Get(Form("hDetaDphiStar4_mix_%s_%s", pidName[ipid], spinName[ispin]));
			if (!h2D_Ang_mix[ipid][1][ispin]) std::cout << "!!! Missing: h2D_Ang_mix" << pidName[ipid] << " !!!" << std::endl;

		}
		
	} // end of pid loop
			
		TCanvas *c3[3][3], *c4[3], *c5[3], *c6[3], *c7, *c8[3], *c9[3], *c10, *c11, *c_weigh[2];
	
	gStyle->SetPadBottomMargin(0.14);
	gStyle->SetPadLeftMargin(0.15);
	

	//--------------------------------------------------
	//  c1  : Angle   		 (particle) = (LL, AA, LA) \ (Same/Mix)
	//  c2  : cosAngle		 (particle) = (LL, AA, LA) \ (Same/Mix)
	//  c3	: QA eta vs phi* (same, mix, same/mix),  (LL, AA, LA)
	//  c4  : Δη	         (same & mix, same/mix), (LL, AA, LA)
	//  c5  : ΔΦ*	         (same & mix, same/mix), (LL, AA, LA)
	//  c6  : Mass	 		 (same, mix, same/mix) \ (SS, SB, BB)
	//  c7	: Qinv	 		 (same, mix, same/mix) \ (SS, SB, BB)
	//  c8  : Low Qinv		 (same, mix, same/mix) \ (All, Para, Anti)
	//  c9  : Mass		 	 (same, mix, same/mix) \ (All, Para, Anti)
	//  c10 : Qinv			 (same, mix, same/mix) \ (All, Para, Anti)
	//--------------------------------------------------

	//===========================================================================
	TH2D *h2D_Ang_ratio[3][4][3] = {nullptr}; // [pid][dgpair][ispin]
	TH2D *h2D_Ang_ratio_pid[3][4] = {nullptr};
	TH2D *h2D_Ang_ratio_spin[3][4][3] = {nullptr};

	for(int ipid=0; ipid<3; ipid++){ // 0:LL, 1:AA, 2:LA
		//--------------------------------------------------
		//  c3: QA eta vs phi* (LL, AA, LA)
		//--------------------------------------------------
		for(int ispin=0; ispin<3; ispin++){

			if (!h2D_Ang[ipid][1][ispin] || !h2D_Ang_mix[ipid][1][ispin]) {
				std::cout << "Skip QA plot for ipid=" << ipid << ", ispin=" << ispin << " due to missing histograms." << std::endl;
				continue;
			}
			c3[ipid][ispin] = new TCanvas(Form("c3_%s_%s", pidName[ipid], spinName[ispin]), Form("%s %s", pidName[ipid], spinName[ispin]), 1000, 700);
			c3[ipid][ispin]->Divide(4,3);

			for(int idgpair=0; idgpair<4; idgpair++){ // 0:pp, 1:ππ, 2:pπ, 3:πp

				if (!h2D_Ang[ipid][idgpair][ispin] || !h2D_Ang_mix[ipid][idgpair][ispin]) {
					std::cout << "Skip drawing: ipid=" << ipid << ", idgpair=" << idgpair << std::endl;
					continue;
				}

				c3[ipid][ispin]->cd(idgpair+1);          
				gPad->SetLogz(); 
				h2D_Ang[ipid][idgpair][ispin]->SetStats(0);
				h2D_Ang[ipid][idgpair][ispin]->SetTitle(Form("Same %s %s %s", dgPairName_2D[idgpair], pidTitle[ipid], spinTitle[ispin]));
				h2D_Ang[ipid][idgpair][ispin]->SetXTitle("#Delta#eta");
				h2D_Ang[ipid][idgpair][ispin]->SetYTitle("#Delta#phi*");
				SetAxis2D(h2D_Ang[ipid][idgpair][ispin]);
				h2D_Ang[ipid][idgpair][ispin]->Draw("colz");

				c3[ipid][ispin]->cd(idgpair+5);          
				gPad->SetLogz(); 
				h2D_Ang_mix[ipid][idgpair][ispin]->SetStats(0);
				h2D_Ang_mix[ipid][idgpair][ispin]->SetTitle(Form("Mix %s %s %s", dgPairName_2D[idgpair], pidTitle[ipid], spinTitle[ispin]));
				h2D_Ang_mix[ipid][idgpair][ispin]->SetXTitle("#Delta#eta");
				h2D_Ang_mix[ipid][idgpair][ispin]->SetYTitle("#Delta#phi*");
				SetAxis2D(h2D_Ang_mix[ipid][idgpair][ispin]);
				h2D_Ang_mix[ipid][idgpair][ispin]->Draw("colz");

				c3[ipid][ispin]->cd(idgpair + 9); 
				gPad->SetLogz(0); 
				//Double_t nSame = h2D_Ang[ipid][idgpair][ispin]->GetEntries();
				//Double_t nMix  = h2D_Ang_mix[ipid][idgpair][ispin]->GetEntries();
				TAxis *xAxis = h2D_Ang[ipid][idgpair][ispin]->GetXaxis();
                int bX_neg_min = xAxis->FindBin(-0.5);
                int bX_neg_max = xAxis->FindBin(-0.2); // 凹み・ピークを避ける
                int bX_pos_min = xAxis->FindBin(0.2);  // 凹み・ピークを避ける
                int bX_pos_max = xAxis->FindBin(0.5);

				int bY_min = 1;
                int bY_max = h2D_Ang[ipid][idgpair][ispin]->GetNbinsY();

				Double_t nSame = h2D_Ang[ipid][idgpair][ispin]->Integral(bX_neg_min, bX_neg_max, bY_min, bY_max)
                               + h2D_Ang[ipid][idgpair][ispin]->Integral(bX_pos_min, bX_pos_max, bY_min, bY_max);

                Double_t nMix  = h2D_Ang_mix[ipid][idgpair][ispin]->Integral(bX_neg_min, bX_neg_max, bY_min, bY_max)
                               + h2D_Ang_mix[ipid][idgpair][ispin]->Integral(bX_pos_min, bX_pos_max, bY_min, bY_max);

				h2D_Ang_ratio[ipid][idgpair][ispin] = (TH2D*)h2D_Ang[ipid][idgpair][ispin]->Clone(Form("h2D_Ang%s_%s", pidName[ipid], dgPairName_2D[idgpair]));
				h2D_Ang_ratio[ipid][idgpair][ispin]->Sumw2();
				if(nSame > 0 && nMix > 0) h2D_Ang_ratio[ipid][idgpair][ispin]->Scale(nMix / nSame);
				h2D_Ang_ratio[ipid][idgpair][ispin]->Divide(h2D_Ang_mix[ipid][idgpair][ispin]);
				h2D_Ang_ratio[ipid][idgpair][ispin]->SetStats(0);
				h2D_Ang_ratio[ipid][idgpair][ispin]->SetMinimum(0.05);
				h2D_Ang_ratio[ipid][idgpair][ispin]->SetMaximum(1.2);
				h2D_Ang_ratio[ipid][idgpair][ispin]->SetTitle(Form("(Same/Mix) %s %s %s", dgPairName_2D[idgpair], pidTitle[ipid], spinTitle[ispin]));
				h2D_Ang_ratio[ipid][idgpair][ispin]->SetXTitle("#Delta#eta");
				h2D_Ang_ratio[ipid][idgpair][ispin]->SetYTitle("#Delta#phi*");
				SetAxis2D(h2D_Ang_ratio[ipid][idgpair][ispin]);
				//h2D_Ang_ratio[ipid][idgpair][ispin]->Draw("colz");

				//float min_X=0.0, max_X=0.0, min_Y=0.0, max_Y=0.0;
				float min_X=-0.1, max_X=0.1, min_Y=-0.2, max_Y=0.2;
/*
				if(ipid==0 | ipid==1){
					min_X = -0.05;
					max_X = 0.05;
					min_Y = -0.05;
					max_Y = 0.05;
				} else if(ipid==2){
					min_X = -0.1;
					max_X = 0.1;
					min_Y = -0.2;
					max_Y = 0.2;
				}
*/
				h2D_Ang[ipid][idgpair][ispin]    ->GetXaxis()->SetRangeUser(-0.5, 0.5);
				h2D_Ang[ipid][idgpair][ispin]    ->GetYaxis()->SetRangeUser(-1.5, 1.5);
				h2D_Ang_mix[ipid][idgpair][ispin]->GetXaxis()->SetRangeUser(-0.5, 0.5);
				h2D_Ang_mix[ipid][idgpair][ispin]->GetYaxis()->SetRangeUser(-1.5, 1.5);
				h2D_Ang_ratio[ipid][idgpair][ispin]	 		  ->GetXaxis()->SetRangeUser(min_X, max_X);
				h2D_Ang_ratio[ipid][idgpair][ispin]	 		  ->GetYaxis()->SetRangeUser(min_Y, max_Y);

			}// end of daughter pair loop
		}// end of spin loop
	
		/*

		//--------------------------------------------------
		//  c4: Δη
		//---------------------------------------------------	
		c4[ipid] = new TCanvas(Form("c4_%s", pidName[ipid]), Form("Delta_eta_%s", pidName[ipid]), 1200, 600);
		c4[ipid]->Divide(4, 2);

		for(int idgpair=0; idgpair<4; idgpair++){ // 0:SS, 1:SB&BS, 2:BB

			// 安全対策: ヒストグラムが存在しない場合は描画をスキップ
			if (!hDeta[ipid][idgpair] || !hDeta_mix[ipid][idgpair]) {
				std::cout << "Skip Δη plot for ipid=" << ipid << " idgpair=" << idgpair << " (missing hist)" << std::endl;
				continue;
			}

			c4[ipid]->cd(idgpair+1);
			hDeta[ipid][idgpair]->SetTitle(Form("Same %s %s", dgPairName[idgpair], pidTitle[ipid]));
			hDeta[ipid][idgpair]->SetXTitle("#Delta#eta");
			hDeta[ipid][idgpair]->Draw();
			SetAxis(hDeta[ipid][idgpair]);
			SetRed(hDeta[ipid][idgpair]);

				TH1 *hDeta_mix_scaled = ScaleHist(hDeta[ipid][idgpair], hDeta_mix[ipid][idgpair], min_DAngle, max_DAngle);
				hDeta_mix_scaled->SetXTitle("#Delta#eta");
				hDeta_mix_scaled->SetTitle(Form("Same/Mix %s %s", dgPairName[idgpair], pidTitle[ipid]));
				hDeta_mix_scaled->Draw("same");
				SetAxis(hDeta_mix_scaled);
				SetBlack(hDeta_mix_scaled);
				SetLegend(hDeta[ipid][idgpair], hDeta_mix_scaled)->Draw();

				c4[ipid]->cd(idgpair+5);
				TH1 *hDeta_ratio = (TH1*)hDeta[ipid][idgpair]->Clone(Form("hDeta_ratio_%s_%s", pidName[ipid], dgPairName[idgpair]));
                
				hDeta_ratio->Sumw2();
				hDeta_ratio->Divide(hDeta_mix_scaled);
				hDeta_ratio->SetTitle(Form("(Same/Mix) %s %s", pidName[ipid], dgPairName[idgpair]));
				hDeta_ratio->SetXTitle("#Delta#eta");
				hDeta_ratio->Draw();
				SetAxis(hDeta_ratio);
			
		} // end of daughter pair Loop
		

		//--------------------------------------------------
		//  c5: ΔΦ*
		//---------------------------------------------------	
		c5[ipid] = new TCanvas(Form("c5_%s", pidName[ipid]), Form("Delta_phi*_%s", pidName[ipid]), 1200, 600);
		c5[ipid]->Divide(4, 2);

		for(int idgpair=0; idgpair<4; idgpair++){ // 0:SS, 1:SB&BS, 2:BB

			// 安全対策: ヒストグラムが存在しない場合は描画をスキップ
			if (!hDphi[ipid][idgpair] || !hDphi_mix[ipid][idgpair]) {
				std::cout << "Skip ΔΦ* plot for ipid=" << ipid << " idgpair=" << idgpair << " (missing hist)" << std::endl;
				continue;
			}

			c5[ipid]->cd(idgpair+1);
			hDphi[ipid][idgpair]->SetXTitle("#Delta#phi*");
			hDphi[ipid][idgpair]->Draw();
			SetAxis(hDphi[ipid][idgpair]);
			SetRed(hDphi[ipid][idgpair]);

			TH1 *hDphi_mix_scaled = ScaleHist(hDphi[ipid][idgpair], hDphi_mix[ipid][idgpair], min_DAngle, max_DAngle);
			hDphi_mix_scaled->SetXTitle("#Delta#phi*");
			hDphi_mix_scaled->SetTitle(Form("Same %s %s", dgPairName[idgpair], pidTitle[ipid]));
			hDphi_mix_scaled->Draw("same");
			SetAxis(hDphi_mix_scaled);
			SetBlack(hDphi_mix_scaled);
			SetLegend(hDphi[ipid][idgpair], hDphi_mix_scaled)->Draw();

			c5[ipid]->cd(idgpair+5);
			TH1 *hDphi_ratio = (TH1*)hDphi[ipid][idgpair]->Clone(Form("hDphi_ratio_%s_%s", pidName[ipid], dgPairName[idgpair]));
			hDphi_ratio->Sumw2();
			hDphi_ratio->Divide(hDphi_mix_scaled);
			hDphi_ratio->SetTitle(Form("(Same/Mix) %s %s", pidName[ipid], dgPairName[idgpair]));
			hDphi_ratio->SetXTitle("#Delta#phi*");
			hDphi_ratio->Draw();
			SetAxis(hDphi_ratio);
			
		} // end of daughter pair Loop	
		 */
	} // end of pid loop

	float R_TPC = 1.7;
	TString title_R = TString::Format("R=%.1f", R_TPC);
		//--------------------------------------------------
		//  c3_ratio_pid: Δη vs Δφ* (LL, AA, LA)
		//--------------------------------------------------
		TCanvas *c3_ratio_pid = new TCanvas("c3_ratio_pid", "QA ratio (pid)", 1000, 800);
		TPad *padMain_pid = new TPad("padMain_pid", "padMain_pid", 0.0, 0.0, 1.0, 0.96);
		padMain_pid->Draw();
		padMain_pid->Divide(4, 3);

		c3_ratio_pid->cd(0);
		TLatex *latex = new TLatex();
		latex->SetNDC();
		latex->SetTextSize(0.03);
		latex->SetTextAlign(22);
		latex->DrawLatex(0.5, 0.975, title_R);

		for(int ipid=0; ipid<3; ipid++){
			for(int idgpair=0; idgpair<4; idgpair++){
				padMain_pid->cd(ipid*4 + idgpair + 1);
				gPad->SetLogz(0);
				
				h2D_Ang_ratio_pid[ipid][idgpair] = (TH2D*)h2D_Ang_ratio[ipid][idgpair][0]->Clone(Form("h2D_Ang_ratio_pid_%s_%s", pidName[ipid], dgPairName_2D[idgpair]));
				SetAxis2D(h2D_Ang_ratio_pid[ipid][idgpair]);
				h2D_Ang_ratio_pid[ipid][idgpair]->Draw("colz");
				
			}
		}

		//--------------------------------------------------
		//  c3_ratio_spin: Δη vs Δφ* (ALL, Para, Anti)
		//--------------------------------------------------
		TCanvas *c3_ratio_spin[3];
		for(int ipid=0; ipid<3; ipid++){

			c3_ratio_spin[ipid]= new TCanvas(Form("c3_ratio_spin_%s", pidName[ipid]), Form("QA ratio (spin) %s", pidName[ipid]), 1000, 800);
			TPad *padMain_spin = new TPad(Form("padMain_spin_%s", pidName[ipid]), "padMain_spin", 0.0, 0.0, 1.0, 0.96);
			padMain_spin->Draw();
			padMain_spin->Divide(4, 3);

			c3_ratio_spin[ipid]->cd(0);
			latex->DrawLatex(0.5, 0.975, title_R);

		
			for(int ispin=0; ispin<3; ispin++){
				for(int idgpair=0; idgpair<4; idgpair++){
					padMain_spin->cd(ispin*4 + idgpair + 1);
					gPad->SetLogz(0);
					if(h2D_Ang_ratio[ipid][idgpair][ispin]){
						SetAxis2D(h2D_Ang_ratio[ipid][idgpair][ispin]);
						h2D_Ang_ratio[ipid][idgpair][ispin]->DrawCopy("colz");
					}
				}
			}
		}

	c3_ratio_pid->SaveAs(Form("h2D_isobar_pid_%.1f.pdf", R_TPC));

	for(int ipid=0; ipid<3; ipid++){
		for(int ispin=0; ispin<3; ispin++){
			//c3[ipid][ispin]->SaveAs(Form("h2D_isobar_%s_%s.pdf", pidName[ipid], spinName[ispin]));	
		}
		
		c3_ratio_spin[ipid]->SaveAs(Form("h2D_isobar_spin_%s_%.1f.pdf", pidName[ipid], R_TPC));
		//c4[ipid]->SaveAs(Form("deta_isobar_%s.pdf", pidName[ipid]));
		//c5[ipid]->SaveAs(Form("dphiStar_isobar_%s.pdf", pidName[ipid]));
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
