TF1 *f_integrandF1;
double QS(double q, double r0);
double F1(double z);
double F2(double z);
TComplex ScAmp(double f0, double d0, double k);
TGraph* makeC2LL(int ndat, double qmin, double qmax, double* par, TGraph *gcomp[3]);
const double hbar = 0.1973;

TGraphErrors *gC2star[2];
TGraph       *gLLstar[2];

void SetAxis( TH1F *ax, float xoff=1.5, float yoff=1.6 );
void SetLegend( TLegend *l );

// functions for fit
Int_t FIT_flag;
Int_t fErrFlag;
Double_t par7[7], epar7[7], chi2;
Int_t ndf;
Int_t andf;
TH1F *hact1D;
TH1F *hbak1D;
TH1F *Hcoul = NULL;
TMinuit *gMinu1D;
void LogLikelihood1D(Int_t &npar, Double_t *gin, Double_t &f, Double_t *u, Int_t flag);
Double_t func1D(Double_t x, Double_t fcw, Double_t *par);
Double_t funcLL(Double_t *x, Double_t *par);
Double_t func_noBW(Double_t *x, Double_t *par);
Double_t BW_term(Double_t Q, Double_t K_scale, Double_t E_H, Double_t Gamma_H) {
    return K_scale * (Gamma_H / (2 * TMath::Pi())) / ( (Q - E_H)*(Q - E_H) + Gamma_H*Gamma_H/4.0 );
}
Double_t func_baseline(Double_t *x, Double_t *par);
void Fit1DLog(TH1F *RHist, TH1F* BHist, TH1F *i_Hcoul);
void Minimization();
void FitAngle();


double QS(double q, double r0){
	return -0.5*exp(-r0*r0*q*q);
}

double F1(double z){
	f_integrandF1->SetParameter(0,z);
	return f_integrandF1->Integral(1.0e-05,z);
}

double F2(double z){
	if( z==0 ){
		cout <<"F2(): invalid input!! z="<<z<<" set to z="<<1e-03<< endl;
	   	z = 1e-03;
	}
	return (1.0-exp(-z*z))/z;
}

TComplex ScAmp(double f0, double d0, double k){
	TComplex fkinv( (1.0/f0+0.5*d0*k*k), -k );
	return 1./fkinv;
}

TGraph* makeC2LL(int ndat, double qmin, double qmax, double* par, TGraph *gcomp[3]){
	// par[0] = lambda
	// par[1] = r0
	// par[2] = f0 (= -a0 in traditional convention) 
	// par[3] = d0 = reff
	// par[4] = a_res
	// par[5] = r_res
	// par[6] = N

	TGraph *g = new TGraph(ndat);
	for( int k=0; k<3; k++ ){ 
		gcomp[k] = new TGraph(ndat);
		gcomp[k]->SetLineWidth(3);
		gcomp[k]->SetLineStyle(7);
	}
	gcomp[0]->SetLineColor(kBlue);
	gcomp[1]->SetLineColor(kRed);
	gcomp[2]->SetLineColor(kGreen+2);

	for(int ix=0; ix<ndat; ix++ ){
		float q = (qmax-qmin)/ndat*ix+qmin;
		double termQS = QS(q/hbar,par[1]);

		TComplex fk = ScAmp(par[2],par[3],q/2./hbar);
		double BI1 = 0.25*fk.Rho2()/par[1]/par[1]*( 1.0-0.5*par[3]/par[1]/sqrt(TMath::Pi()) );
		double BI2 = fk.Re()*F1(q*par[1]/hbar)/sqrt(TMath::Pi())/par[1];
		double BI3 = -0.5*fk.Im()*F2(q*par[1]/hbar)/par[1];

		double termRes = par[4]*exp(-par[5]*par[5]*q*q/hbar/hbar);

		double c2 = 1.0 + par[0]*( termQS + BI1 + BI2 + BI3 ) + termRes; 
		//double c2 = 1.0 + termRes; 
		//double c2 = 1.0 + par[0]*termQS; 

		g->SetPoint( ix, q, c2*par[6] );
		gcomp[0]->SetPoint( ix, q, 1.0+par[0]*termQS );
		gcomp[1]->SetPoint( ix, q, 1.0+par[0]*(BI1+BI2+BI3) );
		gcomp[2]->SetPoint( ix, q, 1.0+termRes );
	}

	return g;
}

// fit
//////////////////////////////////////////////////////////////
void LogLikelihood1D(Int_t &npar, Double_t *gin, Double_t &f, Double_t *u, Int_t flag) {

	f    = 0;
	andf = 0;

	Float_t   fcw = 1;
	for(Int_t binx=1; binx<=hact1D->GetNbinsX(); binx++) {
		if( Hcoul!=NULL ) fcw = (Float_t)Hcoul->GetBinContent(binx);
				
		Double_t x[3] = {0};
		x[0]  = hact1D->GetXaxis()->GetBinCenter(binx);
		//Double_t fu = func1D(x,fcw,u);
		Double_t fu = funcLL(x,u);

		Double_t cu1 = hact1D->GetBinContent(binx);
		Double_t eu1 = hact1D->GetBinError  (binx);
		
		Double_t cu2 = hbak1D->GetBinContent(binx);
		Double_t eu2 = hbak1D->GetBinError  (binx);

		Double_t err = eu1*eu1+eu2*eu2*fu*fu;
		if( cu1 > 0 && cu2 > 0 ) {
			//	  fsum = (fu*cu2 - cu1) - err*TMath::Log( (fu*cu2 - cu1 + err)/err );
			Double_t fsum =  cu1*log( (fu*(cu1+cu2))/(cu1*(fu+1)) ) + cu2*log( (cu1+cu2)/(cu2*(fu+1)) ) ;
			f += fsum;
			andf++;
		}
	}
	f = -2*f;
}


Double_t func1D(Double_t x, Double_t fcw, Double_t *par){

	Double_t v = 0;
	if( FIT_flag==1 )  //----- sinyukov's fitting function ------ //
		v = par[0] * ( 1.0 + par[1] * ( ( 1 + exp((-par[2]*par[2]*x*x)/hbar/hbar) ) * fcw - 1 ) );
	else               //----- simplest fit function ------ //
		v = par[0] * ( 1.0 + par[1] * exp((-par[2]*par[2]*x*x)/hbar/hbar) );

	return v;
}

Double_t funcLL(Double_t *x, Double_t *par){

	TComplex fk  = ScAmp(par[2],par[3],x[0]/2./hbar);
	Double_t BI1 = 0.25*fk.Rho2()/par[1]/par[1]*( 1.0-0.5*par[3]/par[1]/sqrt(TMath::Pi()) );
	Double_t BI2 = fk.Re()*F1(x[0]*par[1]/hbar)/sqrt(TMath::Pi())/par[1];
	Double_t BI3 = -0.5*fk.Im()*F2(x[0]*par[1]/hbar)/par[1];
	Double_t termQS  = QS(x[0]/hbar,par[1]);
	Double_t termRes = par[4]*exp(-par[5]*par[5]*x[0]*x[0]/hbar/hbar);

	Double_t val = par[6] * ( 1.0 + par[0]*( termQS + BI1 + BI2 + BI3 ) + termRes ); 

	return val;
}

Double_t func_withBW(Double_t *x, Double_t *par){
    Double_t q = x[0] / hbar;
    // QS項
    Double_t QS_S0 = exp(-par[1]*par[1]*q*q);
    Double_t QS_S1 = -exp(-par[1]*par[1]*q*q);
    // SI項 (S=0のみ)
    TComplex fk  = ScAmp(par[2], par[3], x[0]/2./hbar);
    Double_t BI1 = 0.25*fk.Rho2()/par[1]/par[1]*( 1.0 - 0.5*par[3]/par[1]/sqrt(TMath::Pi()) );
    Double_t BI2 = fk.Re()*F1(x[0]*par[1]/hbar)/sqrt(TMath::Pi())/par[1];
    Double_t BI3 = -0.5*fk.Im()*F2(x[0]*par[1]/hbar)/par[1];
    Double_t C_S0 = 1.0 + QS_S0 + BI1 + BI2 + BI3;
    Double_t C_S1 = 1.0 + QS_S1;
    // ベースライン
    Double_t baseline = (1.0 + par[5]*exp(-par[6]*par[6]*q*q)) * (par[7] + par[8]*x[0]);
    // シグナル項
    Double_t BW = BW_term(x[0], par[9], par[10], par[11]);
    
    return ( 1.0 + par[0] * ( par[4]*C_S0 + (1.0 - par[4])*C_S1 + BW - 1.0 ) ) * baseline;
}

Double_t func_noBW(Double_t *x, Double_t *par){
	// par[0]: lambda_genuine
	// par[1]: r0
	// par[2]: f0 (S=0)
	// par[3]: d0 (S=0)
	// par[4]: a_weight (S=0 の寄与割合)
	// par[5]: a_res
	// par[6]: r_res
	// par[7]: p0
	// par[8]: p1

	Double_t q = x[0] / hbar;
	
	// 量子統計(QS)項
	Double_t QS_S0 = exp(-par[1]*par[1]*q*q);
	Double_t QS_S1 = -exp(-par[1]*par[1]*q*q);

	// 強い相互作用(SI)項 (S=0 のみ)
	TComplex fk  = ScAmp(par[2], par[3], x[0]/2./hbar);
	Double_t BI1 = 0.25*fk.Rho2()/par[1]/par[1]*( 1.0 - 0.5*par[3]/par[1]/sqrt(TMath::Pi()) );
	Double_t BI2 = fk.Re()*F1(x[0]*par[1]/hbar)/sqrt(TMath::Pi())/par[1];
	Double_t BI3 = -0.5*fk.Im()*F2(x[0]*par[1]/hbar)/par[1];

	// スピン状態ごとの相関関数
	Double_t C_S0 = 1.0 + QS_S0 + BI1 + BI2 + BI3;
	Double_t C_S1 = 1.0 + QS_S1;

	// 非フェムトスコピー効果のベースライン
	Double_t termRes = 1.0 + par[5]*exp(-par[6]*par[6]*q*q);
	Double_t baseline = termRes * (par[7] + par[8]*x[0]);

	// 全体の関数 (シグナル項なし)
	Double_t val = ( 1.0 + par[0] * ( par[4]*C_S0 + (1.0 - par[4])*C_S1 - 1.0 ) ) * baseline;

	return val;
}

Double_t func_baseline(Double_t *x, Double_t *par){
    // par[0]: a_res
    // par[1]: r_res
    // par[2]: p0
    // par[3]: p1
    Double_t q = x[0] / hbar;
    Double_t termRes = 1.0 + par[0]*exp(-par[1]*par[1]*q*q);
    return termRes * (par[2] + par[3]*x[0]);
}


void Minimization(){

	gMinu1D = new TMinuit(7);
	gMinu1D->SetMaxIterations(500);
	gMinu1D->SetFCN(LogLikelihood1D);
	gMinu1D->SetPrintLevel(-1);

	Double_t arglist[10];
	Int_t ierflg = 0;

	arglist[0] = 1;
	gMinu1D->mnexcm("SET ERR",    arglist ,1,ierflg);
	gMinu1D->mnexcm("SET ERRDEF", arglist ,1,ierflg);

	//----  Setting for initial values and step sizes for parameters -----
	Double_t vstart[7] = { 0.3, 2.96, -1.0, 5.0, -0.01, 1.0, 1.0 };
	Double_t step  [7] = { 0.1, 0.10,  0.1, 0.1, 0.001, 0.1, 0.1 };
	TString pname  [7] = { "lambda", "r0", "f0", "d0", "a_res", "r_res", "N" };
	for(int i=0; i<7; i++) gMinu1D->mnparm(i, pname[i], vstart[i], step[i], 0, 0, ierflg);
	gMinu1D->FixParameter(1);
	gMinu1D->FixParameter(6);

	//----- ready for minimization step -----
	arglist[0] = 500;
	arglist[1] = 1.;
	gMinu1D->mnexcm("MINI", arglist, 2, ierflg);

	//----- Print results -----
	Double_t amin,edm,errdef;
	Int_t nvpar,nparx,icstat;
	gMinu1D->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
	gMinu1D->mnprin(7,amin);

	for( int k=0; k<7; k++ ) gMinu1D->GetParameter(k,par7[k],epar7[k]);

	chi2 = amin;
	ndf  = andf;
	fErrFlag = ierflg;
}

void Fit1DLog(TH1F* RHist, TH1F* BHist, TH1F *in_Hcoul, float qlow, float qhig){

	if( in_Hcoul!=NULL ){
		Hcoul->Reset();
		Hcoul->Add( (TH1F *)in_Hcoul );
	}
	
	float sum_r = RHist->Integral( RHist->FindBin(qlow), RHist->FindBin(qhig) );
	float sum_m = BHist->Integral( BHist->FindBin(qlow), BHist->FindBin(qhig) );
	float nfl = 1;
	if( sum_m != 0 ) nfl = sum_r / sum_m;

	hact1D = (TH1F*)RHist->Clone(); 
	hbak1D = (TH1F*)BHist->Clone(); 
	hbak1D->Scale(nfl);
	Minimization();
}


void FitAngle() {
    // 1. ファイルとヒストグラムの読み込み
    //TFile *in = TFile::Open("check_isobar_KehaoMergeTopocut_cuttracks_R1p4.root");
	TFile *in_Para = TFile::Open("smallerbin0_90.root");
	TFile *in_Anti = TFile::Open("smallerbin90_180.root");
    /*if(!in || in->IsZombie()) {
        std::cerr << "Error: Could not open root file." << std::endl;
        return;
    }
	if(!in_Para || !i_Anti) {
        std::cerr << "Error: Could not open  Kehao data root file." << std::endl;
        return;
    }
    if(!in || in->IsZombie()) TH1D *hC2_All  = (TH1D*)in->Get("hCF_All_LL");
    if(!in || in->IsZombie()) TH1D *hC2_Anti = (TH1D*)in->Get("hCF_Anti_LL");
    if(!in || in->IsZombie())TH1D *hC2_Para = (TH1D*)in->Get("hCF_Para_LL");*/
	TH1D *hC2_Para = (TH1D*)in_Para->Get("hCF_cut2_s_cent0");
    TH1D *hC2_Anti = (TH1D*)in_Anti->Get("hCF_cut2_s_cent0");

	double scale = 1000.0; 
	double xmin_Para = hC2_Para->GetXaxis()->GetXmin();
	double xmax_Para = hC2_Para->GetXaxis()->GetXmax();
	hC2_Para->GetXaxis()->SetLimits(xmin_Para / scale, xmax_Para / scale);

	double xmin_Anti = hC2_Anti->GetXaxis()->GetXmin();
	double xmax_Anti = hC2_Anti->GetXaxis()->GetXmax();
	hC2_Anti->GetXaxis()->SetLimits(xmin_Anti / scale, xmax_Anti / scale);
    

    TCanvas *cAll = new TCanvas("cAll", "Angle Dependent CF Fits", 1200, 400);
    cAll->Divide(3, 1);
    /*
    // ========================================================
    // [Pad 1] 0-180度: BWなし、a_weight = 0.25 (固定)
    // ========================================================
    cAll->cd(1);
    TF1 *fit_All = new TF1("fit_All", func_withBW, 0.01, 0.6, 12);
    // 初期値設定 (lambda, r0, f0, d0, a_weight, a_res, r_res, p0, p1, K_scale, E_H, Gamma_H)
    Double_t par_All[12] = {0.18, 1.77, 0.60, 9.80, 0.25, -0.044, 0.68, 0.992, 0.00283, 10.75, 0.1, 0.02};
	if(!hC2_All) return;
	fit_All->FixParameter(0, 0.18);        // lambda (固定)
    fit_All->SetParLimits(1, 1.0, 3.0);    // r0
    fit_All->SetParLimits(2, -10.0, 10.0); // f0
    fit_All->SetParLimits(3, 0.0, 20.0);   // d0
    fit_All->FixParameter(4, 0.25);        // a_weight (0-180度は0.25で固定)
    fit_All->SetParLimits(5, -0.2, 0.2);   // a_res
    fit_All->SetParLimits(6, 0.0, 2.0);    // r_res
    fit_All->SetParLimits(7, 0.9, 1.1);    // p0 (ベースライン高さ)
    fit_All->SetParLimits(8, 0.0, 0.01);   // p1 (ベースライン傾き)
    fit_All->SetParLimits(9, 0.0, 20.0);   // K_scale (BW振幅)
    fit_All->SetParLimits(10, 0.08, 0.12); // E_H (Q ~ 100MeV/c付近のピーク位置)
    fit_All->SetParLimits(11, 0.015, 0.05);

    fit_All->SetParameters(par_All);
    fit_All->FixParameter(4, 0.25); // a_weight を 0.25 (1/4) に固定
    fit_All->SetLineColor(kRed);
    hC2_All->SetTitle("All Spin");
    hC2_All->GetXaxis()->SetRangeUser(0.0, 0.6);
    hC2_All->Draw("P");
    hC2_All->Fit(fit_All, "R");

    // ☆ 重み関数の重ね書き (シグナル項なし)
    TF1 *weight_All = new TF1("weight_All", func_noBW, 0.01, 0.6, 9);
    for(int i = 0; i < 9; i++) {
        weight_All->SetParameter(i, fit_All->GetParameter(i));
    }
    weight_All->SetLineColor(kMagenta);
    weight_All->SetLineStyle(2);
    weight_All->Draw("same");

	TLegend *leg1 = new TLegend(0.45, 0.2, 0.85, 0.4);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->AddEntry(hC2_All, "Data", "lp");
    leg1->AddEntry(fit_All, "LL+Baseline+BW", "l");
    leg1->AddEntry(weight_All, "LL+Baseline", "l");
    leg1->Draw();
*/
    // ========================================================
    // [Pad 2] 90-180度: BWなし、a_weight フリー ＋ 重み関数重ね書き
    // ========================================================
    cAll->cd(2);
    TF1 *fit_Anti = new TF1("fit_Anti", func_withBW, 0.01, 0.6, 12);
    Double_t par_Anti[12] = {0.18, 1.64, 0.60, 9.80, 0.32, -0.07, 0.90, 0.992, 0.00247, 5.0, 0.1, 0.02};

	fit_Anti->FixParameter(0, 0.18);       // lambda (固定)
    fit_Anti->SetParLimits(1, 1.0, 3.0);   // r0
    fit_Anti->SetParLimits(2, -10.0, 10.0);// f0
    fit_Anti->SetParLimits(3, 0.0, 20.0);  // d0
    fit_Anti->SetParLimits(4, 0.0, 1.0);   // a_weight (フリー)
    fit_Anti->SetParLimits(5, -0.2, 0.2);  // a_res
    fit_Anti->SetParLimits(6, 0.0, 2.0);   // r_res
    fit_Anti->SetParLimits(7, 0.9, 1.1);   // p0
    fit_Anti->SetParLimits(8, 0.0, 0.01);  // p1
    fit_Anti->SetParLimits(9, 0.0, 10.0);  // K_scale
    fit_Anti->SetParLimits(10, 0.08, 0.12);// E_H
    fit_Anti->SetParLimits(11, 0.015, 0.05);// Gamma_H

    fit_Anti->SetParameters(par_Anti);
    fit_Anti->SetLineColor(kRed);
    hC2_Anti->SetTitle("Anti-Parallel");
    hC2_Anti->GetXaxis()->SetRangeUser(0.0, 0.6);
    hC2_Anti->Draw("P");
    hC2_Anti->Fit(fit_Anti, "R");

    // ☆ 重み関数の重ね書き (シグナル項なし)
    TF1 *weight_Anti = new TF1("weight_Anti", func_noBW, 0.01, 0.6, 9);
    for(int i = 0; i < 9; i++) {
        weight_Anti->SetParameter(i, fit_Anti->GetParameter(i));
    }
    weight_Anti->SetLineColor(kMagenta);
    weight_Anti->SetLineStyle(2);
    weight_Anti->Draw("same");

	// Baseline
	TF1 *draw_Base = new TF1("draw_Base", func_baseline, 0.01, 0.6, 4);
    
    // 全体フィット(fit_Anti)のパラメータ 5,6,7,8 をセット
    draw_Base->SetParameters(
        fit_Anti->GetParameter(5), // a_res
        fit_Anti->GetParameter(6), // r_res
        fit_Anti->GetParameter(7), // p0
        fit_Anti->GetParameter(8)  // p1
    );
    
    draw_Base->SetLineColor(kBlack);       // 好みの色に
    draw_Base->SetLineStyle(7);           // 破線など
    draw_Base->Draw("same");

	TLegend *leg2 = new TLegend(0.45, 0.2, 0.85, 0.4);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->AddEntry(hC2_Anti, "Data", "lp");
    leg2->AddEntry(fit_Anti, "LL+Baseline+BW", "l");
    leg2->AddEntry(weight_Anti, "LL+Baseline", "l");
	leg2->AddEntry(draw_Base, "Baseline", "l");
    leg2->Draw();

    // ========================================================
    // [Pad 3] 0-90度: BWなし、a_weight フリー
    // ========================================================
    cAll->cd(3);
    TF1 *fit_Para = new TF1("fit_Para", func_noBW, 0.01, 0.6, 9);
    // 初期値設定 (lambda, r0, f0, d0, a_weight, a_res, r_res, p0, p1)
    fit_Para->SetParameters(0.18, 1.83, 0.60, 9.80, 0.20, -0.032, 0.54, 0.992, 0.00302);
/*
	fit_Para->FixParameter(0, 0.18);
    fit_Para->SetParLimits(1, 1.0, 3.0);
    fit_Para->SetParLimits(2, -10.0, 10.0);
    fit_Para->SetParLimits(3, 0.0, 20.0);
    fit_Para->SetParLimits(4, 0.0, 1.0);
    fit_Para->SetParLimits(5, -0.2, 0.2);
    fit_Para->SetParLimits(6, 0.0, 2.0);
    fit_Para->SetParLimits(7, 0.9, 1.1);
    fit_Para->SetParLimits(8, 0.0, 0.01);*/

    fit_Para->SetLineColor(kRed);
    hC2_Para->SetTitle("Parallel");
    hC2_Para->GetXaxis()->SetRangeUser(0.0, 0.6);
    hC2_Para->Draw("P");
    hC2_Para->Fit(fit_Para, "R");

	TLegend *leg3 = new TLegend(0.45, 0.2, 0.85, 0.4);
    leg3->SetBorderSize(0);
    leg3->SetFillStyle(0);
    leg3->AddEntry(hC2_Para, "Data", "lp");
    leg3->AddEntry(fit_Para, "LL+Baseline", "l");
    leg3->Draw();

	TFile *fOutWeight = new TFile("weight_func_KehaoTopoMerge.root", "RECREATE");
    // 0.0 ~ 2.0 GeV/c を 2000分割（十分に細かい精度）
    TH1D *hWeight = new TH1D("hWeight", "Weight Function for Mix", 2000, 0.0, 2.0);
    for(int i = 1; i <= 2000; i++){
		
        double q = hWeight->GetBinCenter(i);
        hWeight->SetBinContent(i, weight_Anti->Eval(q));
    }
    hWeight->Write();

	TH1D *hWeight_BW = new TH1D("hWeight_BW", "Weight Function for Mix", 2000, 0.0, 2.0);
    for(int i = 1; i <= 2000; i++){
        double q = hWeight_BW->GetBinCenter(i);
        hWeight_BW->SetBinContent(i, fit_Anti->Eval(q));
    }
    hWeight_BW->Write();
    fOutWeight->Close();

    std::cout << "✅ Saved Weight Function to weight_func.root" << std::endl;
}