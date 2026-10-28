TF1 *f_integrandF1;
double QS(double q, double r0);
double F1(double z);
double F2(double z);
TComplex ScAmp(double f0, double d0, double k);
TGraph* makeC2LL(int ndat, double qmin, double qmax, double* par, TGraph *gcomp[3]);
const double hbar = 0.1973;

TGraphErrors *gC2star[2];
TGraph       *gLLstar[2];
void LoadSTARC2();

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
void FitPML();
void FitAngle();

void fitLL(){

	gStyle->SetOptStat(1111111);

	LoadSTARC2();
	f_integrandF1 = new TF1("integrandF1","exp(x*x-[0]*[0])/[0]",0,1);

	TF1 *fit1D = new TF1("fit1D",funcLL,1e-05,1,7);
	fit1D->SetLineColor(2);

	FitPML();


	// lambda, r0, f0, d0, a_res, r_res, N
	double para    [7] = { 0.18, 2.96, -1.10, 8.52, -0.044, 0.43, 1.006 }; // STAR PRL114.022301.2015
	int ndata = 50;

	TGraph *gcomp[3];
	TGraph *gC2 = makeC2LL(ndata,0.001,1.001,para,gcomp);
	gC2->SetMarkerStyle(24);
	gC2->SetMarkerSize(0.6);
	gC2->SetLineColor(1);
	gC2->SetLineWidth(3);


	TF1 *one = new TF1("one","[0]",0,5);
	one->SetParameter(0,1);
	one->SetLineStyle(2);
	one->SetLineColor(kGray+1);

	TLatex lt;
	lt.SetNDC();
	lt.SetTextSize(13);
	lt.SetTextFont(43);

	TLegend *lg = new TLegend(0.4,0.44,0.7,0.58);
	SetLegend(lg);
	lg->AddEntry( gC2star[0], "#Lambda#Lambda PRL114.022301", "P" );
	lg->AddEntry( gLLstar[1], "STAR Lednicly fit", "L" );

	TLegend *lg2 = new TLegend(0.4,0.22,0.7,0.4);
	SetLegend(lg2);
	lg2->AddEntry( gC2, "Lednicky model", "L" );
	lg2->AddEntry( gcomp[0], "QS", "L" );
	lg2->AddEntry( gcomp[1], "SI", "L" );
	lg2->AddEntry( gcomp[2], "residual", "L" );


	TCanvas *c1 = new TCanvas("c1","c1",300,300);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	TH1F *ax = gPad->DrawFrame( 0, 0.65, 0.7, 1.1);
	SetAxis(ax,1.1,1.1);
	ax->SetXTitle("Q_{inv} [GeV/c] ");
	ax->SetYTitle("C_{2} ");
	one->Draw("same");
	gC2star[0]->Draw("P");
	gLLstar[1]->Draw("L");

	gC2    ->Draw("L");
	gcomp[0]->Draw("L");
	gcomp[1]->Draw("L");
	gcomp[2]->Draw("L");
	lg->Draw();
	lg2->Draw();
	lt.DrawLatex( 0.35, 0.83, "r_{0}=2.96, f_{0}=-1.10, d_{0}=8.52" );


	TGraphErrors *gData = (TGraphErrors*)gC2star[0]->Clone("gData");

	TCanvas *c2 = new TCanvas("c2","c2",300,300);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	ax = gPad->DrawFrame( 0, 0.65, 0.7, 1.1);
	SetAxis(ax,1.1,1.1);
	ax->SetXTitle("Q_{inv} [GeV/c] ");
	ax->SetYTitle("C_{2} ");
	one->Draw("same");
	gData->Draw("P");

	//double para    [7] = { 0.18, 2.96, -1.10, 8.52, -0.044, 0.43, 1.006 }; // STAR PRL114.022301.2015
	// lambda, r0, f0, d0, a_res, r_res, N
	fit1D->SetParNames( "lambda", "r0", "f0", "d0", "a_res", "r_res", "N" );
	fit1D->SetParameters( 0.3, 2.96, -1, 5, -0.01, 1.0, 1.0 );
	fit1D->FixParameter(1,2.96);
	//fit1D->FixParameter(4,-0.04);
	//fit1D->FixParameter(5,0.43);
	//fit1D->FixParameter(4,0);
	//fit1D->FixParameter(5,0.43);
	fit1D->FixParameter(6,1.0);
	//fit1D->SetParLimits(6,0.8,1.2); // norm
	TVirtualFitter::SetMaxIterations(9000);
	gData->Fit(fit1D,"R","",0.01,0.5);

	//lt.DrawLatex( 0.4, 0.56, "fit results (fixed r_{0}=2.96)");
	lt.DrawLatex( 0.4, 0.50, Form("N = %3.2f#pm%3.2f",      fit1D->GetParameter(6),fit1D->GetParError(6)) );
	lt.DrawLatex( 0.4, 0.45, Form("#lambda = %3.2f#pm%3.2f",fit1D->GetParameter(0),fit1D->GetParError(0)) );
	lt.DrawLatex( 0.4, 0.40, Form("f_{0} = %3.2f#pm%3.2f",  fit1D->GetParameter(2),fit1D->GetParError(2)) );
	lt.DrawLatex( 0.4, 0.35, Form("d_{0} = %3.2f#pm%3.2f",  fit1D->GetParameter(3),fit1D->GetParError(3)) );
	lt.DrawLatex( 0.4, 0.30, Form("a_{res} = %3.2f#pm%3.2f",fit1D->GetParameter(4),fit1D->GetParError(4)) );
	lt.DrawLatex( 0.4, 0.25, Form("r_{res} = %3.2f#pm%3.2f",fit1D->GetParameter(5),fit1D->GetParError(5)) );

	TF1 *fitWeight = new TF1("fitWeight", func_noBW, 1e-05, 1, 9);
	
	// Figure 73 の "LL" (ピンク破線) に合わせて色とスタイルを設定
	fitWeight->SetLineColor(kMagenta); 
	fitWeight->SetLineStyle(1); 
	
	// パラメータの設定 (Table 11: 90°-180°領域のフィット結果を参考に初期値をセット)
	// lambda_genuineは論文に合わせて仮置き(例: 0.18)しています。
	fitWeight->SetParameters(0.18, 1.64, 0.60, 9.80, 0.32, -0.07, 0.90, 0.992, 0.00247);

	c2->cd(); // グラフを描画するキャンバスを指定
	fitWeight->Draw("same"); // 既存のグラフの上に重ね書き
	
	// 凡例(テキスト)の追加
	lt.SetTextColor(kMagenta);
	lt.DrawLatex( 0.4, 0.20, "LL + Baseline (No Signal)" );
	lt.SetTextColor(kBlack);

	double parFit[7];
	for( int k=0; k<7; k++ ) parFit[k] = fit1D->GetParameter(k);
	TGraph *gcompFit[3];
	TGraph *gC2Fit = makeC2LL(ndata,0.001,1.001,parFit,gcompFit);
	gcompFit[0]->Draw("L");
	gcompFit[1]->Draw("L");
	gcompFit[2]->Draw("L");

	FitAngle();
}

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

// STAR Lambda-C2
void LoadSTARC2(){

	float q   [25];
	float C2  [25];
	float eC2 [25];
	float seC2[25];
	float LL  [25];
	float LLr [25];
	float eq  [25]={0};

	//Data for figure - 3
	//Q(GeV/c)      C(Q)               Stat_err                Sys_err                LL_noresidual       LL_withresidual
	q[ 0]=0.01;     C2[ 0]=0.8228;     eC2[ 0]=0.0767992 ;     seC2[ 0]=0.118214 ;    LL[ 0]=0.909636 ;   LLr[ 0]=0.835034;
	q[ 1]=0.03;     C2[ 1]=0.8775;     eC2[ 1]=0.0297943 ;     seC2[ 1]=0.116843 ;    LL[ 1]=0.910942 ;   LLr[ 1]=0.8524685;
	q[ 2]=0.05;     C2[ 2]=0.855 ;     eC2[ 2]=0.0174406 ;     seC2[ 2]=0.0479224;    LL[ 2]=0.9134897;   LLr[ 2]=0.879512;
	q[ 3]=0.07;     C2[ 3]=0.9129;     eC2[ 3]=0.0132806 ;     seC2[ 3]=0.0350589;    LL[ 3]=0.9171532;   LLr[ 3]=0.906085;
	q[ 4]=0.09;     C2[ 4]=0.9356;     eC2[ 4]=0.0105106 ;     seC2[ 4]=0.0418513;    LL[ 4]=0.9217586;   LLr[ 4]=0.9256351;
	q[ 5]=0.11;     C2[ 5]=0.9295;     eC2[ 5]=0.00860613;     seC2[ 5]=0.0253097;    LL[ 5]=0.9270960;   LLr[ 5]=0.9367368;
	q[ 6]=0.13;     C2[ 6]=0.9256;     eC2[ 6]=0.00732562;     seC2[ 6]=0.0221106;    LL[ 6]=0.932936 ;   LLr[ 6]=0.9411453;
	q[ 7]=0.15;     C2[ 7]=0.9463;     eC2[ 7]=0.00648223;     seC2[ 7]=0.0221893;    LL[ 7]=0.9390461;   LLr[ 7]=0.94188 ;
	q[ 8]=0.17;     C2[ 8]=0.9505;     eC2[ 8]=0.00572692;     seC2[ 8]=0.0200541;    LL[ 8]=0.9452055;   LLr[ 8]=0.943844;
	q[ 9]=0.19;     C2[ 9]=0.9417;     eC2[ 9]=0.00512124;     seC2[ 9]=0.0223172;    LL[ 9]=0.9512183;   LLr[ 9]=0.951205;
	q[10]=0.21;     C2[10]=0.9618;     eC2[10]=0.00476733;     seC2[10]=0.0264969;    LL[10]=0.956922 ;   LLr[10]=0.959835;
	q[11]=0.23;     C2[11]=0.9661;     eC2[11]=0.00436716;     seC2[11]=0.0273789;    LL[11]=0.96219  ;   LLr[11]=0.965885;
	q[12]=0.25;     C2[12]=0.9803;     eC2[12]=0.00410732;     seC2[12]=0.0131275;    LL[12]=0.96695  ;   LLr[12]=0.96983 ;
	q[13]=0.27;     C2[13]=0.9757;     eC2[13]=0.00379868;     seC2[13]=0.0238427;    LL[13]=0.97115  ;   LLr[13]=0.972706;
	q[14]=0.29;     C2[14]=0.9691;     eC2[14]=0.0035437 ;     seC2[14]=0.0121333;    LL[14]=0.974789 ;   LLr[14]=0.975055;
	q[15]=0.31;     C2[15]=0.9767;     eC2[15]=0.00335767;     seC2[15]=0.0150337;    LL[15]=0.977886 ;   LLr[15]=0.97714 ;
	q[16]=0.33;     C2[16]=0.9826;     eC2[16]=0.00319313;     seC2[16]=0.0154463;    LL[16]=0.980490 ;   LLr[16]=0.979087;
	q[17]=0.35;     C2[17]=0.9776;     eC2[17]=0.00301848;     seC2[17]=0.0106255;    LL[17]=0.9826677;   LLr[17]=0.980939;
	q[18]=0.37;     C2[18]=0.9772;     eC2[18]=0.00289148;     seC2[18]=0.0154458;    LL[18]=0.9844946;   LLr[18]=0.982722;
	q[19]=0.39;     C2[19]=0.9815;     eC2[19]=0.0027688 ;     seC2[19]=0.0127251;    LL[19]=0.986054 ;   LLr[19]=0.984445;
	q[20]=0.41;     C2[20]=0.9874;     eC2[20]=0.00267598;     seC2[20]=0.0107529;    LL[20]=0.9874315;   LLr[20]=0.986108;
	q[21]=0.43;     C2[21]=0.9916;     eC2[21]=0.00257775;     seC2[21]=0.0127466;    LL[21]=0.9887008;   LLr[21]=0.987708;
	q[22]=0.45;     C2[22]=0.9892;     eC2[22]=0.00247993;     seC2[22]=0.0108048;    LL[22]=0.989926 ;   LLr[22]=0.98924 ;
	q[23]=0.47;     C2[23]=0.9914;     eC2[23]=0.00240154;     seC2[23]=0.0110781;    LL[23]=0.9911504;   LLr[23]=0.990710;
	q[24]=0.49;     C2[24]=0.9959;     eC2[24]=0.00234215;     seC2[24]=0.0109378;    LL[24]=0.992395 ;   LLr[24]=0.992104;
	
	gC2star[0] = new TGraphErrors(25,q,C2,eq,eC2);
	gC2star[1] = new TGraphErrors(25,q,C2,eq,seC2);
	gLLstar[0] = new TGraph(25,q,LL);
	gLLstar[1] = new TGraph(25,q,LLr);

	gC2star[0]->SetMarkerStyle(24);
	gLLstar[1]->SetLineColor(kMagenta);
}

void SetAxis( TH1F *ax, float xff, float yoff ){
	ax->SetTitleFont(43,"X");
	ax->SetTitleSize(16,"X");
	ax->SetTitleOffset(xff,"X");
	ax->SetTitleFont(43,"Y");
	ax->SetTitleSize(16,"Y");
	ax->SetTitleOffset(yoff,"Y");
	ax->SetNdivisions(505,"X");
	ax->SetNdivisions(505,"Y");
	ax->SetLabelFont(43,"X");
	ax->SetLabelFont(43,"Y");
	ax->SetLabelSize(12,"X");
	ax->SetLabelSize(12,"Y");
}

void SetLegend( TLegend *l ){

	l->SetFillColor(0);
	l->SetBorderSize(0);
	l->SetTextFont(43);
	l->SetTextSize(13);
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

void FitPML(){

	cout <<"FitPM()............."<< endl;
	//TFile *in = TFile::Open("c2test.root");
	TFile *in = TFile::Open("c2test.root");
	if(!in) return;
	TH1F *hQr = (TH1F*)in->Get("hQr");
	TH1F *hQm = (TH1F*)in->Get("hQm");
	TH1F *hC2 = (TH1F*)in->Get("hC2");

	TH1F *hCoul = NULL;
	Fit1DLog(hQr,hQm,hCoul,0.5,1.0);

	TGraph *gcomp[3];
	TGraph *gC2 = makeC2LL(50,0.001,1.001,par7,gcomp);
	gC2->SetMarkerStyle(24);
	gC2->SetMarkerSize(0.6);
	gC2->SetLineColor(1);
	gC2->SetLineWidth(3);

	TCanvas *cPML = new TCanvas("cPML","cPML",400,400);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	TH1F *ax = gPad->DrawFrame( 0, 0.65, 0.7, 1.1);
	SetAxis(ax,1.1,1.1);
	ax->SetXTitle("Q_{inv} [GeV/c] ");
	ax->SetYTitle("C_{2} ");
	hC2->Draw("P,same");

	gC2->Draw("L");
	gcomp[0]->Draw("L");
	gcomp[1]->Draw("L");
	gcomp[2]->Draw("L");
	
	cout <<"FitPM()............."<< endl;
}

void FitAngle() {
    // 1. ファイルとヒストグラムの読み込み
    TFile *in = TFile::Open("check_isobar_KehaoMergeTopocut_cuttracks_R1p4.root");
    if(!in || in->IsZombie()) {
        std::cerr << "Error: Could not open root file." << std::endl;
        return;
    }
    TH1D *hC2_All  = (TH1D*)in->Get("hCF_All_LL");
    TH1D *hC2_Anti = (TH1D*)in->Get("hCF_Anti_LL");
    TH1D *hC2_Para = (TH1D*)in->Get("hCF_Para_LL");

    TCanvas *cAll = new TCanvas("cAll", "Angle Dependent CF Fits", 1200, 400);
    cAll->Divide(3, 1);
    
    // ========================================================
    // [Pad 1] 0-180度: BWなし、a_weight = 0.25 (固定)
    // ========================================================
    cAll->cd(1);
    TF1 *fit_All = new TF1("fit_All", func_withBW, 0.01, 0.6, 12);
    // 初期値設定 (lambda, r0, f0, d0, a_weight, a_res, r_res, p0, p1, K_scale, E_H, Gamma_H)
    Double_t par_All[12] = {0.18, 1.77, 0.60, 9.80, 0.25, -0.044, 0.68, 0.992, 0.00283, 10.75, 0.1, 0.02};

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