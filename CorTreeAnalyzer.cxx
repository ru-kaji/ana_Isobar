//ISOBAR
// My header files
#include "CorTreeAnalyzer.h"
#include "EventClass.h"
#include "StRefMultCorr/StRefMultCorr.h"

// ROOT header files 
#include "TSystem.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TVector2.h"
#include "TMath.h"
#include "TLorentzVector.h"

// C/C++ heder
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

enum StBeamDirection {east   = 0,
                      yellow = 0,    // yellow beam is going west -> east
                      west   = 1,
                      blue   = 1};   // blue beam is going east -> west

const Float_t massLam    = 1.115683;
const Float_t massLamMin = 1.105; //AuAu: 1.1105
const Float_t massLamMax = 1.1205;

int nPair = 0;

const Float_t meanH[3][3] = {
	{ 1.11578, 1.32199, 1.67263 }, // run11
	{ 1.11606, 1.32225, 1.67277 }, // run14
	{ 1.11610, 1.32228, 1.67283 }  // run16
};

const Float_t sgm_para_Lm[3][3] = {
	{0.00173021, -8.55547e-06, 5.84024e-08}, // run11
	{0.00207058, -4.47636e-06, 2.81876e-08}, // run14
	{0.00212029, -5.08523e-06, 2.65629e-08}, // run16
};

// The constructor. Initialize data members here.
CorTreeAnalyzer::CorTreeAnalyzer() : mEventCounter(0), mFile(0), mRunNumber(0), ANA_MODE(0), mTrkEff(0), mSelfCor(0), mCalibMode(0), mTotalEvents(0)
{ 
	mFileName = "muAnalysis.root";
	gRandom->SetSeed(0);
	mEventCounterRef = 0;
	mColl = 0;
}

CorTreeAnalyzer::CorTreeAnalyzer(const Char_t *outfile, EventClass* in_p_eve) : mEventCounter(0), mFile(0), mRunNumber(0), ANA_MODE(0), mTrkEff(0), mSelfCor(0), mCalibMode(0)
{ 
	mFileName = outfile;
	gRandom->SetSeed(0);
	p_eve = in_p_eve;
	mEventCounterRef = 0;
	mColl = 0;
}

CorTreeAnalyzer::~CorTreeAnalyzer() { /* noop */ }

// Called once at the beginning.
Int_t CorTreeAnalyzer::Init()
{
		// Output histograms
    mFile =  new TFile( mFileName.c_str(), "RECREATE", "Output of CorTreeAnalyzer: 3-PC analysis", 9 );
    mFile->cd();
	cout << "\tHistograms will be stored in file '"
	 <<  mFileName.c_str() << "'" << endl;

	


	// Global variables
	// -------------------------------------------------------

	// Track variables
	// -------------------------------------------------------

	// Event planes and flow vectors
	// -------------------------------------------------------
	
	// Topological cuts
	// -------------------------------------------------------
	// 0-10%
	mDcaP[0][0] = 0.5;  mDcaPi[0][0] = 1.6;  mDcaV0[0][0] = 0.4;  mDL[0][0] = 5.0;  mDcaDh[0][0] = 1.0;
	mDcaP[0][1] = 0.5;  mDcaPi[0][1] = 1.7;  mDcaV0[0][1] = 0.4;  mDL[0][1] = 7.0;  mDcaDh[0][1] = 1.0;
	// 10-20%
	mDcaP[1][0] = 0.3;  mDcaPi[1][0] = 1.6;  mDcaV0[1][0] = 0.4;  mDL[1][0] = 5.0;  mDcaDh[1][0] = 1.0;
	mDcaP[1][1] = 0.4;  mDcaPi[1][1] = 1.7;  mDcaV0[1][1] = 0.4;  mDL[1][1] = 7.0;  mDcaDh[1][1] = 1.0;
	// 20-30%
	mDcaP[2][0] = 0.3;  mDcaPi[2][0] = 1.3;  mDcaV0[2][0] = 0.4;  mDL[2][0] = 4.0;  mDcaDh[2][0] = 1.0;
	mDcaP[2][1] = 0.4;  mDcaPi[2][1] = 1.6;  mDcaV0[2][1] = 0.4;  mDL[2][1] = 5.0;  mDcaDh[2][1] = 1.0;
	// 30-40%
	mDcaP[3][0] = 0.3;  mDcaPi[3][0] = 1.2;  mDcaV0[3][0] = 0.4;  mDL[3][0] = 3.0;  mDcaDh[3][0] = 1.0;
	mDcaP[3][1] = 0.3;  mDcaPi[3][1] = 1.3;  mDcaV0[3][1] = 0.4;  mDL[3][1] = 4.0;  mDcaDh[3][1] = 1.0;
	// 40-50%
	mDcaP[4][0] = 0.2;  mDcaPi[4][0] = 1.0;  mDcaV0[4][0] = 0.4;  mDL[4][0] = 3.0;  mDcaDh[4][0] = 1.0;
	mDcaP[4][1] = 0.3;  mDcaPi[4][1] = 1.2;  mDcaV0[4][1] = 0.4;  mDL[4][1] = 4.0;  mDcaDh[4][1] = 1.0;
	// 50-60%
	mDcaP[5][0] = 0.1;  mDcaPi[5][0] = 0.9;  mDcaV0[5][0] = 0.4;  mDL[5][0] = 3.0;  mDcaDh[5][0] = 1.0;
	mDcaP[5][1] = 0.3;  mDcaPi[5][1] = 1.2;  mDcaV0[5][1] = 0.4;  mDL[5][1] = 4.0;  mDcaDh[5][1] = 1.0;
	// 60-70%
	mDcaP[6][0] = 0.1;  mDcaPi[6][0] = 0.8;  mDcaV0[6][0] = 0.4;  mDL[6][0] = 3.0;  mDcaDh[6][0] = 1.0;
	mDcaP[6][1] = 0.3;  mDcaPi[6][1] = 1.2;  mDcaV0[6][1] = 0.4;  mDL[6][1] = 4.0;  mDcaDh[6][1] = 1.0;
	// 70-80%
	mDcaP[7][0] = 0.1;  mDcaPi[7][0] = 0.7;  mDcaV0[7][0] = 0.4;  mDL[7][0] = 3.0;  mDcaDh[7][0] = 1.0;
	mDcaP[7][1] = 0.3;  mDcaPi[7][1] = 1.2;  mDcaV0[7][1] = 0.4;  mDL[7][1] = 4.0;  mDcaDh[7][1] = 1.0;


	for( int ict=0; ict<8; ict++){
		mDcaP [ict][2] = 0.5;  
		mDcaPi[ict][2] = 1.5;  
		mDcaV0[ict][2] = 0.4;  //元は1.0
		mDL   [ict][2] = 4.0;  
		mDcaDh[ict][2] = 1.0;
	}


	// Ach study
	mAchId = -1;
	for( int ict=0; ict<NCent; ict++ ) mAcWd[ict] = 1.0;
	//SetChAsymWidthRun11();
	SetChAsymWidthRun14();

	string sSw[] = { "OFF", "ON" };
	//string sSwEff[] = { "OFF", "ON (run14/run11)", "ON (run14/run11+embedding-eff.)" };
	string sSwEff[] = { "OFF", "ON (run11 eff)", "ON (run11-eff * relative-eff)" };

	cout << endl << "Set up parameters" << endl;
	cout << "=========================================" << endl; 
	cout << "Q Calib mode = " << mCalibMode << endl;
	cout << "Track efficiency is " << sSwEff[mTrkEff].c_str() << endl;
	cout << "Removing self-correltaion is " << sSw[mSelfCor].c_str() << endl;
	cout << "Systemaic flag = "<< mSysFlag << endl;
	cout << "=========================================" << endl << endl; 

	DefHistograms();
	if( mSelfCor || mCalibMode ){
		if( LoadRunList()==-1 )     return -1;
	}

	TFile *fWeight = TFile::Open("weight_func.root", "READ");
    if(fWeight && !fWeight->IsZombie()){
        hWeightFunc = (TH1D*)fWeight->Get("hWeight");
        if(hWeightFunc)hWeightFunc->SetDirectory(0); // ファイルを閉じてもメモリに保持させる

		hWeightFunc_BW = (TH1D*)fWeight->Get("hWeight_BW");
		if(hWeightFunc_BW)hWeightFunc_BW->SetDirectory(0); // ファイルを閉じてもメモリに保持させる
        fWeight->Close();
        std::cout << "✅ Loaded weight function (hWeight_BW) successfully!" << std::endl;
    } else {
        std::cout << "⚠️ Warning: Could not load weight_func.root! Weight will be 1.0" << std::endl;
        hWeightFunc = nullptr;
        hWeightFunc_BW = nullptr;
    }

	/*if( mSelfCor ){
		if( LoadEPcalibPara()==-1 ) return -1;
	}*/

	if( mTrkEff ){ 
		//if( mTrkEff>1 && SetEfficiencyRun11()==-1 )         return -1;
		//if( mTrkEff>0 && SetRelativeEfficiencyRun14()==-1 ) return -1;
		//if( SetEfficiencyRun11()==-1 )         return -1;

		// isobar
		SetNormEffRun11();
		SetRelativeEffIsobarToRun11();
	
		cout <<"eff(0.6)="<<f_NormEffPtRun11->Eval(0.6) << endl;
		cout <<"eff(0.8)="<<f_NormEffPtRun11->Eval(0.8) << endl;
		cout <<"eff(1.0)="<<f_NormEffPtRun11->Eval(1.0) << endl;
	}

	fRefMultCorr = new StRefMultCorr("refmult","Isobar");

	// for mass window
	//====================================================
	int irun;
	if     ( mRunYear==2011 ) irun = 0;
	else if( mRunYear==2014 ) irun = 1;
	else                      irun = 2;
	
	f_sgmIM[0] = new TF1( "f_sgmIM", "pol2", 0, 100 );
	for( int k=0; k<3; k++ ){
		f_sgmIM[0]->SetParameter(k,sgm_para_Lm[irun][k]);
	}

    return 0;
}

// Called every event after Make(). 
void CorTreeAnalyzer::Clear()
{
	p_eve->InitValue();
}

// Called once at the end.
Int_t CorTreeAnalyzer::End()
{
	// Summarize the run.
    cout << "CorTreeAnalyzer::End()\n";
    cout << "\tProcessed " << mEventCounter << " events." << endl;
	cout << "unmatched events with StRefMult = "<< mEventCounterRef <<", fraction = "<< (float)mEventCounterRef/mEventCounter << endl; 


	// Write histos to file and close it.
    if( mFile ){
      mFile->cd();
		WriteHistograms();
      //mFile->Write();  
    	mFile->Close();
    }

    return 0;
}

bool CorTreeAnalyzer::IsPairClose_pp(float eta1, float phi1, float eta2, float phi2, int ipid) {
    
    float deta = fabs(eta1 - eta2);
	if(deta > 0.05) return false; // my cut

    float dphi = TVector2::Phi_mpi_pi(phi1 - phi2); // -pi ~ pi の範囲で差分を取る

	if(ipid == 2){
		if(deta < 0.025 && fabs(dphi)<0.05 ) return true;

		//if( dphi < -0.1 && dphi > -0.4 ) { //R_TPC=0.7
		if( -0.08 < dphi && dphi < 0.16 ) {   //R_TPC=1.5
			return true; // too close (Merging or Splitting)
		}
	}else if(ipid == 0 || ipid == 1){
		float cut_eta = 0.02;
		float cut_phi = 0.03;

		float r2 = (deta * deta) / (cut_eta * cut_eta) + (dphi * dphi) / (cut_phi * cut_phi);
	
		if (r2 < 1.0) return true; 
	}
	//if(deta < 0.02 && fabs(dphi) < 0.03) return true; // Kehao cut

    return false;
}

bool CorTreeAnalyzer::IsPairClose_pipi(float eta1, float phi1, float eta2, float phi2, int ipid) {
    
    float deta = fabs(eta1 - eta2);
	if(deta >= 0.05) return false;

    float dphi = TVector2::Phi_mpi_pi(phi1 - phi2); // -pi ~ pi の範囲で差分を取る
	
	if(ipid == 0 || ipid == 1){
		float cut_eta = 0.03;
		float cut_phi = 0.03;

		float r2 = (deta * deta) / (cut_eta * cut_eta) + (dphi * dphi) / (cut_phi * cut_phi);
	
		if (r2 < 1.0) return true; 
	}
	
    return false;
}

bool CorTreeAnalyzer::IsPairClose_ppi(float eta1, float phi1, float eta2, float phi2, int ipid) {
    
    float deta = fabs(eta1 - eta2);
	if(deta >= 0.02) return false;

    float dphi = TVector2::Phi_mpi_pi(phi1 - phi2); // -pi ~ pi の範囲で差分を取る

	if(ipid == 2){
		if(0.02<dphi && dphi<0.04) return true;
	}
	
    return false;
}

bool CorTreeAnalyzer::IsPairClose_pip(float eta1, float phi1, float eta2, float phi2, int ipid) {
    
    float deta = fabs(eta1 - eta2);
	if(deta >= 0.02) return false;

    float dphi = TVector2::Phi_mpi_pi(phi1 - phi2); // -pi ~ pi の範囲で差分を取る

	if(ipid == 2){
		if(0.02<dphi && dphi<0.04) return true;
	}
	
    return false;
}

bool CorTreeAnalyzer::IsSplitCheck(float eta1, float phi1, float eta2, float phi2) {
    
    float deta = fabs(eta1 - eta2);
    float dphi = TVector2::Phi_mpi_pi(phi1 - phi2); // -pi ~ pi の範囲で差分を取る

	const float CUT_DETA = 0.02;
	const float CUT_DPHI = 0.02;

	if( deta < CUT_DETA && fabs(dphi) < CUT_DPHI ) {
		return true;
	}
	
    return false;
}

// Phi Normalization →(-pi < phi < pi) 
double CorTreeAnalyzer::NormalizePhi(double phi) {
	if (std::isnan(phi) || std::isinf(phi)) return -999.0;

	// fmodを使って、何周回っていても一気に -2π 〜 2π の範囲に収める
    phi = std::fmod(phi, 2 * M_PI); //fmod : 整数以外にも使える、余りを返す

    if (phi > M_PI) phi -= 2 * M_PI;
    if (phi < -M_PI) phi += 2 * M_PI;
    return phi;
}

const double B_TESLA = -0.5;
const double R_VERTEX = 0.0;
//const double R_TPC_METERS = 1.5; 
const double R_TPC_METERS = 0.7; 


double CorTreeAnalyzer::GetPhiStar(double phi_meas, double pt, double charge){

	// phi_meas : 測定点(R_VERTEX)でのphi
	// R_VERTEX   : phiが測定された半径 (m)
	// R_TPC_METERS : 計算したい目的の半径 (m) ※マージ判定を行う位置
	// 曲率定数 C = 0.3*B*q/2

	// 0.3は光速由来の定数 (0.29979...)
	const double k = 0.299792458 * B_TESLA * charge / 2.0;

	// arcsin(C*R/pt)
	double arg_meas   = (k * R_VERTEX) / pt;		//arg?, meas?
	double arg_target = (k * R_TPC_METERS) / pt;

	// 低pT粒子(Looper)対策 
    // 引数が -1 ~ 1 の範囲を超えると asin は NaN を返すため除外または値を丸める
    if (std::abs(arg_meas) > 1.0 || std::abs(arg_target) > 1.0) {
        // その半径まで粒子が到達しない(回転してしまう)場合
        // ここでは便宜上、元のphiを返すが、解析によっては除外フラグを立てるべき
        return -999.9; 
	}


    // 厳密式: phi_new = phi_meas + asin(meas) - asin(target)
    // 符号の組み合わせは座標系定義によりますが、
    // 「Vertex(R=0)でのphi」に戻す項と、「目的のRへ進める項」の差になります。
    double term1 = phi_meas;
    double term2 = std::asin(arg_meas);
    double term3 = std::asin(arg_target);

    // Trackの曲がる方向と座標系に合わせて符号を調整してください
    // 一般的な定義では以下の形（R_VERTEXからR_TPC_METERSへの射影）
    double phi_star = term1 + term2 - term3;   //term1 + term2 + term3にしていた

    return NormalizePhi(phi_star);
}

double CorTreeAnalyzer::GetAngle(double phi1, double phi2, double theta1, double theta2){

	double p1_x = std::sin(theta1) * std::cos(phi1);
    double p1_y = std::sin(theta1) * std::sin(phi1);
    double p1_z = std::cos(theta1);
    double p2_x = std::sin(theta2) * std::cos(phi2);
    double p2_y = std::sin(theta2) * std::sin(phi2);
    double p2_z = std::cos(theta2);

	double cross_x = p1_y * p2_z - p1_z * p2_y;
	double cross_y = p1_z * p2_x - p1_x * p2_z;
	double cross_z = p1_x * p2_y - p1_y * p2_x;

	double sinTheta = std::sqrt(cross_x*cross_x + cross_y*cross_y + cross_z*cross_z);
	double signed_sin = (cross_z >= 0) ? sinTheta : -sinTheta;
	double cosTheta = (p1_x * p2_x) + (p1_y * p2_y) + (p1_z * p2_z);

	double angle_rad = std::atan2(signed_sin, cosTheta);

    return angle_rad * 180.0 / TMath::Pi(); // -180° ~ 180°
}

double CorTreeAnalyzer::GetcosAngle(double phi1, double phi2, double theta1, double theta2){

	double p1_x = std::sin(theta1) * std::cos(phi1);
    double p1_y = std::sin(theta1) * std::sin(phi1);
    double p1_z = std::cos(theta1);
    double p2_x = std::sin(theta2) * std::cos(phi2);
    double p2_y = std::sin(theta2) * std::sin(phi2);
    double p2_z = std::cos(theta2);

	double cosAngle = (p1_x * p2_x) + (p1_y * p2_y) + (p1_z * p2_z);

	if (cosAngle > 1.0) cosAngle = 1.0;
    if (cosAngle < -1.0) cosAngle = -1.0;

	return cosAngle;
}


// This method is called every event.
Int_t CorTreeAnalyzer::ProcessEvent()
{
    mEventCounter++; 
	if( (mEventCounter%100000)==0 ) cout << "Processed events = " << mEventCounter << "/" << mTotalEvents<< endl;

	if     ( mRunYear==2011 ) mRunFlag = 0;
	else if( mRunYear==2014 ) mRunFlag = 1;
	else  					  mRunFlag = 2;

    if( !p_eve ){
		cout << "CorTreeAnalyzer::ProcessEvent() : p_eve is NULL !!" << endl;
    	return -1;      
    }

	//if( p_eve->getrunid() ){
	//	for( int ich=0; ich<2; ich++ ){
	//		for( int ict=0; ict<8; ict++ ){
	//			f_EffPtIsoRltvRu[ich][ict] = f_EffPtIsoRltvRu[ich][ict];
	//		}
	//	}
	//}

	// event info
    short refMult = p_eve->get_refMult();
    short fQcent  = p_eve->get_centid();  // 16 bins
    float vz      = p_eve->get_vz();
    float fTrgEff = 1.0; //p_eve->get_trgEff();

	if (mEventCounter <= 5) {
    	std::cout << "Event " << mEventCounter << ": fQcent = " << fQcent << ", vz = " << vz << std::endl;
	}

    if( fQcent<0 || fQcent>=16 ) return 0;
    //if( fabs(vz)>6 )            return 0;
    if( vz<-35 || vz>25 )        return 0;

	//int iCent = (15-fQcent) / 2 - 1; // 0=10-20%, 6=70-80%
	//int Cent = (iCent+1.5)*10.0;
	int iCent = (15-fQcent) / 2; 
		// 0=0-10%, 1=10-20%, 2=20-30%, 3=30-40%, 4=40-50%, 5=50-60%, 6=60-70%, 7=70-80%
	int Cent = (iCent+0.5)*10.0;

	// cent9
	//int fCent = iCent + 2; // iCent=0->10-20%
	int fCent = iCent + 1; // iCent=0->0-10%
	if     ( fQcent==15 ) fCent = 0;
	else if( fQcent==14 ) fCent = 1;


	// StRefMultCorr
	fRefMultCorr->init(p_eve->get_runid());
	fRefMultCorr->initEvent(refMult,vz,p_eve->get_zdcRate());
	fTrgEff = fRefMultCorr->getWeight();
	if( fRefMultCorr->getRefMultCorr()>50 ) fTrgEff = 1.0;

	//if( fQcent!=fRefMultCorr->getCentralityBin16() ){ 
	//	//cout <<"refM="<< refMult <<" centid-tree="<< fQcent <<" centid-stref="<< fRefMultCorr->getCentralityBin16() <<" trgW="<< fRefMultCorr->getWeight() << endl;
	//	mEventCounterRef++;
	//} // ~4% difference likely due to float vs double


	hVz ->Fill(vz);
	hRM ->Fill(refMult);
	hRMw->Fill(refMult,fTrgEff);
	hRMc ->Fill(fRefMultCorr->getRefMultCorr());
	hRMcw->Fill(fRefMultCorr->getRefMultCorr(),fTrgEff);
	hCentBin ->Fill(fQcent);
	hCentBinw->Fill(fQcent,fTrgEff);
	hCentBin2 ->Fill(fRefMultCorr->getCentralityBin16());
	hCentBinw2->Fill(fRefMultCorr->getCentralityBin16(),fTrgEff);

	// vz
	//float vzstep = 1.0; // vz bin step
	//int fVz = (int)( (vz+6.0)/vzstep ); 
	float vzstep = 5.0; // vz bin step
	int fVz = (int)( (vz+35.0)/vzstep ); 
	if( fVz<0 || fVz>11 ) return 0;

	// calibrate raw EP
	// ------------------------------------------------------
	float fPsi[2][3];
	size_t runindex;
	if( mCalibMode || mSelfCor ){
		int runid    = p_eve->get_runid();
		auto result = std::find( vrunid.begin(), vrunid.end(), runid );
		if( result==vrunid.end() ){
			cout <<"Failed to find run-index!"<< endl;
			return 0;
		}
		runindex = std::distance( vrunid.begin(), result );

		for( int ih=0; ih<2; ih++ ){

			// backward
			float tpcQxb = p_eve->get_tpcQxB(ih) + p_eve->get_tpcQxB(ih+4); // ch>0 + ch<0
			float tpcQyb = p_eve->get_tpcQyB(ih) + p_eve->get_tpcQyB(ih+4); // ch>0 + ch<0
			// forward
			float tpcQxf = p_eve->get_tpcQxF(ih) + p_eve->get_tpcQxF(ih+4); // ch>0 + ch<0
			float tpcQyf = p_eve->get_tpcQyF(ih) + p_eve->get_tpcQyF(ih+4); // ch>0 + ch<0
			// combined 
			float tpcQxbf = tpcQxb + tpcQxf;
			float tpcQybf = tpcQyb + tpcQyf;

			float tpcQw[3];
			tpcQw[0] = p_eve->get_tpcQwB(0) + p_eve->get_tpcQwB(2); // M_chp + M_chm (forward)
			tpcQw[1] = p_eve->get_tpcQwF(0) + p_eve->get_tpcQwF(2); // (backward)
			tpcQw[2] = tpcQw[0] + tpcQw[1];

			float fPsi_raw[3];
			fPsi_raw[0] = atan2(  tpcQyb, tpcQxb  ) / (ih+1.0);
			fPsi_raw[1] = atan2(  tpcQyf, tpcQxf  ) / (ih+1.0);
			fPsi_raw[2] = atan2( tpcQybf, tpcQxbf ) / (ih+1.0);


			if( mCalibMode ){ 
				hQvRec[ih][0][fQcent][fVz]->Fill( runindex, tpcQxb );
				hQvRec[ih][1][fQcent][fVz]->Fill( runindex, tpcQyb );
				hQvRec[ih][2][fQcent][fVz]->Fill( runindex, tpcQxf );
				hQvRec[ih][3][fQcent][fVz]->Fill( runindex, tpcQyf );
				hQvRec[ih][4][fQcent][fVz]->Fill( runindex, tpcQxbf );
				hQvRec[ih][5][fQcent][fVz]->Fill( runindex, tpcQybf );
			}else{

				// recentering
				doRecenterQv( tpcQxb,ih,0,0,fQcent,fVz,runindex);
				doRecenterQv( tpcQyb,ih,1,0,fQcent,fVz,runindex);
				doRecenterQv( tpcQxf,ih,0,1,fQcent,fVz,runindex);
				doRecenterQv( tpcQyf,ih,1,1,fQcent,fVz,runindex);
				doRecenterQv(tpcQxbf,ih,0,2,fQcent,fVz,runindex);
				doRecenterQv(tpcQybf,ih,1,2,fQcent,fVz,runindex);

				hQvRec[ih][0][fQcent][fVz]->Fill( runindex, tpcQxb );
				hQvRec[ih][1][fQcent][fVz]->Fill( runindex, tpcQyb );
				hQvRec[ih][2][fQcent][fVz]->Fill( runindex, tpcQxf );
				hQvRec[ih][3][fQcent][fVz]->Fill( runindex, tpcQyf );
				hQvRec[ih][4][fQcent][fVz]->Fill( runindex, tpcQxbf );
				hQvRec[ih][5][fQcent][fVz]->Fill( runindex, tpcQybf );


				fPsi[ih][0] = atan2(  tpcQyb, tpcQxb  ) / (ih+1.0);
				fPsi[ih][1] = atan2(  tpcQyf, tpcQxf  ) / (ih+1.0);
				fPsi[ih][2] = atan2( tpcQybf, tpcQxbf ) / (ih+1.0);

				for( int iep=0; iep<3; iep++ ){
					if( tpcQw[iep]>0 ){
						hPsi[ih][iep][fQcent/2][0]->Fill( fPsi_raw[iep] );
						hPsi[ih][iep][fQcent/2][1]->Fill( fPsi[ih][iep] );
					}
				}

				if( tpcQw[0]>0 && tpcQw[1]>0 ){
					hPsiCor[ih][0]->Fill( fCent, cos( (ih+1.0)*(fPsi[ih][0]-fPsi[ih][1]) ) );
					hPsiCor[ih][1]->Fill( fCent, sin( (ih+1.0)*(fPsi[ih][0]-fPsi[ih][1]) ) );
				}
			}
		}// ih

		if( mCalibMode ) return 0;
	}// mCalib or mSelfCor

	// ------------------------------------------------------
	//if( iCent<0 ) return 0; // reject 0-10% centrality


	// get EPs 
	// ------------------------------------------------------
	float zdcPsi[3];
	for( int iep=0; iep<3; iep++ ){
		zdcPsi[iep] = p_eve->get_zdcPsi(iep);
	}
	if( zdcPsi[0]>-100 && zdcPsi[1]>-100 ){
		hZDCPsiCor[0]->Fill( fCent, cos(zdcPsi[0]-zdcPsi[1]));
		hZDCPsiCor[1]->Fill( fCent, sin(zdcPsi[0]-zdcPsi[1]));
		hZDCPsiCor[2]->Fill( fCent, cos(2.0*(zdcPsi[0]-zdcPsi[1])) );
		hZDCPsiCor[3]->Fill( fCent, sin(2.0*(zdcPsi[0]-zdcPsi[1])) );
	}

	float bbcPsiN[2][3];
	for( int isub=0; isub<3; isub++ ){
		bbcPsiN[0][isub] = p_eve->get_bbcPsi1(isub);
		bbcPsiN[1][isub] = p_eve->get_bbcPsi2(isub);
	}
	for( int ih=0; ih<2; ih++ ){
		if( bbcPsiN[ih][0]>-100 && bbcPsiN[ih][1]>-100 ){
			hBBCPsiCor[ih]->Fill( fCent, cos( (ih+1.0)*(bbcPsiN[ih][0]-bbcPsiN[ih][1])) );
		}
	}

	float tpcPsiN[3][6];
	for( int isub=0; isub<6; isub++ ){ // (E,W,E+W) x (|eta|>0.1,|eta|>0.5)
		tpcPsiN[0][isub] = p_eve->get_tpcPsi2(isub);
		tpcPsiN[1][isub] = p_eve->get_tpcPsi3(isub);
		tpcPsiN[2][isub] = p_eve->get_tpcPsi4(isub);
	}
	for( int ih=0; ih<3; ih++ ){
		if( tpcPsiN[ih][0]>-100 && tpcPsiN[ih][1]>-100 ){
			hTPCPsiCor[ih]->Fill( fCent, cos( (ih+2.0)*(tpcPsiN[ih][0]-tpcPsiN[ih][1]) ) );
		}
		if( tpcPsiN[ih][3]>-100 && tpcPsiN[ih][4]>-100 ){
			hTPCPsiCor[ih+3]->Fill( fCent, cos( (ih+2.0)*(tpcPsiN[ih][3]-tpcPsiN[ih][4]) ) );
		}
	}

	float epdPsiN[3][3];
	for( int isub=0; isub<3; isub++ ){ // (E,W,E+W) 
		epdPsiN[0][isub] = p_eve->get_epdPsi2(isub);
		epdPsiN[1][isub] = p_eve->get_epdPsi3(isub);
		epdPsiN[2][isub] = p_eve->get_epdPsi4(isub);
	}
	for( int ih=0; ih<3; ih++ ){
		if( epdPsiN[ih][0]>-100 && epdPsiN[ih][1]>-100 ){
			hEPDPsiCor[ih]->Fill( fCent, cos( (ih+2.0)*(epdPsiN[ih][0]-epdPsiN[ih][1]) ) );
		}
	}
	if( p_eve->get_epdPsi1(2)>-100 && zdcPsi[0]>-100 && zdcPsi[1]>-100 ){
		hEPDZDCPsi1Cor[0]->Fill( fCent, cos(p_eve->get_epdPsi1(2)-zdcPsi[0]) );
		hEPDZDCPsi1Cor[1]->Fill( fCent, cos(p_eve->get_epdPsi1(2)-zdcPsi[1]) );
	}


	if( tpcPsiN[0][2]>-100 ) hCentIso8->Fill((15-fQcent)/2);

	std::vector< pair<float,float> > prtA[2];
	std::vector< pair<float,float> > prtB[2];
	std::vector< pair<float,float> > prtAtpc[2];
	std::vector< pair<float,float> > prtBtpc[2];


	

	// Lambda Track loop
	// ----------------------------------------------------------
    for (int l=0; l<(int)p_eve->get_ntrk(); l++){ 

		//if( !isGoodTrack(p_eve,l) ) continue;
    	
    	float ptv0  = p_eve->get_ptv0(l);
    	float etav0 = p_eve->get_etav0(l);
    	float phiv0 = p_eve->get_phiv0(l);
		//float rapv0 = p_eve->get_rap(l);  // rapidity
		short tofid = p_eve->get_tofid(l);
		short ch    = p_eve->get_charge(l);
		int fCH = ( ch>0 ) ? 0 : 1;
		
		float invM  = p_eve->get_mass(l);
		float dcaA  = p_eve->get_dcaA(l); // proton
		float dcaB  = p_eve->get_dcaB(l); // pion
		float v0dca = p_eve->get_v0dca(l);
		float dgdca = p_eve->get_dgdca(l);
		float dL    = p_eve->get_dL(l);

		// daughter baryon info
    	float ptA    = p_eve->get_ptA(l);
    	float etaA   = p_eve->get_etaA(l);
    	float phiA   = p_eve->get_phiA(l);

		// daughter meson info (trancated)
		float ptB  = p_eve->get_ptB(l);
		float etaB = p_eve->get_etaB(l);
		float phiB = p_eve->get_phiB(l);

		//float yv0 = log( (sqrt(massLam*massLam+ptv0*ptv0*TMath::CosH(etav0)*TMath::CosH(etav0)) + ptv0*TMath::SinH(etav0)) / sqrt(massLam*massLam+ptv0*ptv0) );
		TLorentzVector lvec;
		lvec.SetPtEtaPhiM(ptv0,etav0,phiv0,massLam);

		// proton in Lambda rest frame
		float phiA_RF = p_eve->get_phARF(l);
		float theA_RF = p_eve->get_thARF(l);
		float cosTheStr = cos(theA_RF);

		if( fabs(lvec.Rapidity())>0.8 ) continue;
		//if( ptv0<0.5 || ptv0>6 )      continue;
		if( ptv0<0.6 )      continue;
		if( ptv0>3.0 )      continue; //@
		if( tofid<0 || tofid>3 )      continue;

		// v0 eta and pt bin
		int fEta = (int)( (etav0+1.0)/(2.0/6) );
		int iPt = -1;
		if     ( ptv0<0.8 ) iPt = 0; 
		else if( ptv0<1.1 ) iPt = 1; 
		else if( ptv0<1.4 ) iPt = 2; 
		else if( ptv0<1.8 ) iPt = 3;
		else if( ptv0<2.5 ) iPt = 4; 
		else if( ptv0<3.5 ) iPt = 5; 
		else                iPt = 6;

		if( iPt<0 || iPt>=NptPol ) continue;


		int iEt = -1;
		if( etav0>-0.9 && etav0<=-0.6 ) iEt = 0;
		else if( etav0<=-0.3 ) iEt = 1;
		else if( etav0<= 0.0 ) iEt = 2;
		else if( etav0<= 0.3 ) iEt = 3;
		else if( etav0<= 0.6 ) iEt = 4;
		else if( etav0<= 0.9 ) iEt = 5;



		float fTrkEff = 1.0;
		//if( mTrkEff>0 ) fTrkEff = getTrkEff(ptv0,iCent,fCH,0); 
		if( mTrkEff>0 ) fTrkEff = getTrkEffIso(ptv0,iCent,fCH); 

		// systematic study
		//if( mSysFlag ){ 

		//	for( int ih=0; ih<3; ih++ ){
		//		SysStudyN(ih,iCent,fCH,iPt,tpcPsiN[ih][2],fTrgEff,fTrkEff,invM,cosTheStr,dcaA,dcaB,dgdca,dL,v0dca,phiv0);
		//	}
		//}

		// default cut 
		if( !passTPcut(dcaA,dcaB,v0dca,dL,dgdca,iCent,0,0) ) continue;
		//cout << "before fill hinvM, invM=" << invM << endl;
		//hinvM->Fill(invM);
		//cout << "after fill hinvM, invM=" << invM << endl;

		if( invM>massLamMin && invM<massLamMax ) hAPTvsCent[fCH]->Fill( Cent, ptv0 );


		// for tof matching efficiecny
		if( fabs(etav0)<0.5 ){
			hIMvsPtv0[fCH][iCent][tofid]->Fill( invM, ptv0, fTrgEff );
		}
		
		
		if( invM>1.07 && invM<1.20 ){ 

			//cout << "DEBUG: invM = " << invM << endl;
			// Λ pt cut
			TLorentzVector Lvec;
			Lvec.SetPtEtaPhiM(ptv0, etav0, phiv0, invM);
			float rap = Lvec.Rapidity();
			if( ptv0 < 0.6 ) continue;
			if( ptv0 > 3.0 ) continue; 
			if( fabs(rap) > 0.8 ) continue;

			// only select tof-matched
			if( tofid==0 || tofid==3 ) prtA[fCH].push_back( make_pair(p_eve->get_phiA(l),p_eve->get_etaA(l)) );
			else                       prtAtpc[fCH].push_back( make_pair(p_eve->get_phiA(l),p_eve->get_etaA(l)) );

			if( tofid==0 || tofid==2 ) prtB[fCH].push_back( make_pair(p_eve->get_phiB(l),p_eve->get_etaB(l)) );    // matched
			else                       prtBtpc[fCH].push_back( make_pair(p_eve->get_phiB(l),p_eve->get_etaB(l)) ); // unmatched
		

			// +reserve daughters 
			trkid [fCH].push_back( l );
			mass  [fCH].push_back( p_eve->get_mass(l) );
			dgApt [fCH].push_back( p_eve->get_ptA (l) );
			dgAeta[fCH].push_back( p_eve->get_etaA(l) );
			dgAphi[fCH].push_back( p_eve->get_phiA(l) );
			dgBpt [fCH].push_back( p_eve->get_ptB (l) );
			dgBeta[fCH].push_back( p_eve->get_etaB(l) );
			dgBphi[fCH].push_back( p_eve->get_phiB(l) );
		}

    }//track loop

	if( mSharedTrkFlag  == 2){
		id_shared.clear();
		pid_shared.clear();
		
		// 1. 同種粒子間 (LL, AA) の Sharing Check
		for( int ip=0; ip<2; ip++ ){ 
			for(unsigned int itrk1=0; itrk1<trkid[ip].size(); itrk1++ ){
				for(unsigned int itrk2=itrk1+1; itrk2<trkid[ip].size(); itrk2++ ){
					if( 
						   (fabs(dgApt[ip][itrk1]-dgApt[ip][itrk2]) < 1.0e-5 && fabs(dgAeta[ip][itrk1]-dgAeta[ip][itrk2]) < 1.0e-5 && fabs(dgAphi[ip][itrk1]-dgAphi[ip][itrk2]) < 1.0e-5) // p vs. p
						|| (fabs(dgBpt[ip][itrk1]-dgBpt[ip][itrk2]) < 1.0e-5 && fabs(dgBeta[ip][itrk1]-dgBeta[ip][itrk2]) < 1.0e-5 && fabs(dgBphi[ip][itrk1]-dgBphi[ip][itrk2]) < 1.0e-5) // pi vs. pi
					){
						// 平均質量 (meanH[mRunFlag][0]) から遠い方を fake とする
						int id_fake = ( fabs(mass[ip][itrk1] - meanH[mRunFlag][0]) > fabs(mass[ip][itrk2] - meanH[mRunFlag][0]) ) ? trkid[ip][itrk1] : trkid[ip][itrk2];
						id_shared.push_back( id_fake );
						pid_shared.push_back( ip );
					}
				}
			}
		}

		// 2. 異種粒子間 (LA) の Sharing Check (クロスチェック: p vs pi)
		for(unsigned int itrk1=0; itrk1<trkid[0].size(); itrk1++ ){ // 0: Lambda
			for(unsigned int itrk2=0; itrk2<trkid[1].size(); itrk2++ ){ // 1: AntiLambda
				if( 
					(fabs(dgApt[0][itrk1]-dgBpt[1][itrk2]) < 1.0e-5 && fabs(dgAeta[0][itrk1]-dgBeta[1][itrk2]) < 1.0e-5 && fabs(dgAphi[0][itrk1]-dgBphi[1][itrk2]) < 1.0e-5) 
					|| (fabs(dgBpt[0][itrk1]-dgApt[1][itrk2]) < 1.0e-5 && fabs(dgBeta[0][itrk1]-dgAeta[1][itrk2]) < 1.0e-5 && fabs(dgBphi[0][itrk1]-dgAphi[1][itrk2]) < 1.0e-5) 
				){
					// 質量が本物のΛ(meanH[mRunFlag][0])から遠い方を捨てる
					int id_fake = -1;
					int pid_fake = -1;
					if ( fabs(mass[0][itrk1] - meanH[mRunFlag][0]) > fabs(mass[1][itrk2] - meanH[mRunFlag][0]) ) {
						id_fake = trkid[0][itrk1];
						pid_fake = 0;
					} else {
						id_fake = trkid[1][itrk2];
						pid_fake = 1;
					}
					id_shared.push_back( id_fake );
					pid_shared.push_back( pid_fake );
				}
			}
		}

		// 3. 重複してリストアップされたFake IDを整理
		std::vector<std::pair<int, int>> combined;
		for (size_t i = 0; i < id_shared.size(); ++i) {
			combined.push_back({id_shared[i], pid_shared[i]});
		}
		std::sort(combined.begin(), combined.end());
		combined.erase( std::unique(combined.begin(), combined.end()), combined.end() );

		id_shared.clear();
		pid_shared.clear();
		for (const auto& p : combined) {
			id_shared.push_back(p.first);     
			pid_shared.push_back(p.second);   
		}
	}

	// define daughter baryon in parent rest frame
		vector<float> 			PhiDB[3], ThetaDB[3]; //[pid]
		vector<DhVector> 		currentDhs[3]; //[pid]
		vector<TLorentzVector> 	currentPrts[3]; //[pid]
		vector<int> 			regionid[3]; //[pid]

	for(int pid = 0; pid < 2; pid++){ // 0:Λ, 1:Λbar 
	
		vector<int> 			&prts_ids = trkid[pid];
	
		regionid[pid].assign(prts_ids.size(), -1);
		currentDhs[pid].resize(prts_ids.size());
		currentPrts[pid].resize(prts_ids.size());
		PhiDB[pid].resize(prts_ids.size());
		ThetaDB[pid].resize(prts_ids.size());

		for(size_t i=0; i<prts_ids.size(); i++){
			int id = prts_ids[i];

			// check daughter sharing
			bool isFake = false;
			if(mSharedTrkFlag == 2){
				for(size_t m=0; m<id_shared.size(); m++){
					if(id_shared[m] == id && pid_shared[m] == pid){
						isFake = true;
						break;
					}
				}
			}
			if(isFake) continue;

			// lambda
				float pt   = p_eve->get_ptv0(id);
				float eta  = p_eve->get_etav0(id); 
				float phi  = p_eve->get_phiv0(id);
				float mass = p_eve->get_mass(id);

				hinvM[pid][iCent]->Fill(mass);

				//cout  << "mass=" << mass <<", pt=" << pt << ", eta=" << eta << endl;

				double mass_diff = fabs(mass-meanH[mRunFlag][0]); //pid=0(Lambda)で固定
				double sgm_L = f_sgmIM[0]->Eval(iCent*10+5);

				if( mass_diff < sgm_L*3.0 ){
					
					hpt[pid]->Fill(pt);
					heta[pid]->Fill(eta);
					hphi[pid]->Fill(phi);
					//cout << "before fill hinvM_Sig" << endl;
					hinvM_Sig[pid]->Fill(mass);
					//cout << "after fill hinvM_Sig" << endl;
					
					currentPrts[pid][i].SetPtEtaPhiM(pt, eta, phi, 1.115683);
					regionid[pid][i] = 0;

				}else if( mass_diff > sgm_L*5.0 && mass_diff < sgm_L*8.0 ){

					float invM_BG = mass;
					//cout << "before fill hinvM_BG" << endl;
					hinvM_BG[pid]->Fill(invM_BG);
					//cout << "after fill hinvM_BG" << endl;
					currentPrts[pid][i].SetPtEtaPhiM(pt, eta, phi, invM_BG); //pdgを使うべき？
					regionid[pid][i] = 1;
				}

				if( !isfinite(p_eve->get_phARF(id)) || !isfinite(p_eve->get_thARF(id)) ) continue;
				PhiDB[pid][i] = p_eve->get_phARF(id);
				ThetaDB[pid][i] = p_eve->get_thARF(id);
				
				// daughter
				currentDhs[pid][i] = DhVector(dgAeta[pid][i], dgAphi[pid][i], dgApt[pid][i],
											  dgBeta[pid][i], dgBphi[pid][i], dgBpt[pid][i]);
		
		} // end of Lambda Loop 
	} // end of pid

	float etaA1[3], phiA1[3], ptA1[3];
	float etaB1[3], phiB1[3], ptB1[3];
	float etaA2[3], phiA2[3], ptA2[3];
	float etaB2[3], phiB2[3], ptB2[3];
	double phiStarA1[3], phiStarA2[3], phiStarB1[3], phiStarB2[3];
	float deta_pp[3], deta_pipi[3], deta_ppi[3], deta_pip[3];
	float dphi_pp[3], dphi_pipi[3], dphi_ppi[3], dphi_pip[3];
	
	// ==============================================================================================================
	// Same Event Loop
	// ==============================================================================================================
	
	for(int ipid = 0; ipid < 3; ipid++){ //0: LL, 1: AA, 2: LA
		//if(ipid!=2) continue;
		//cout << "=== Same Event Loop start: ipid = " << ipid << " ===" << endl;

		int pid1 = (ipid == 1) ? 1 : 0; 
		int pid2 = (ipid == 0) ? 0 : 1;

		for(size_t i=0; i<trkid[pid1].size(); i++){
			size_t j_start = (ipid == 2) ? 0 : i+1;
			//cout << "--- i Loop start ---" << endl;
			for(size_t j=j_start; j<trkid[pid2].size(); j++){
			//cout << "--- j Loop start ---" << endl;

			size_t idx1 = i;
			size_t idx2 = j;

			if((ipid ==0 || ipid == 1) && gRandom->Rndm() > 0.5) {
				idx1 = j;
				idx2 = i;
			}

			// region cut
			//---------------------------------------------------------------
			int r1 = regionid[pid1][idx1];
			int r2 = regionid[pid2][idx2];
			
			if(r1==-1||r2==-1){
				//cout << "DEBUG: skip by regionid error, r1=" << r1 << ",r2=" << r2 << endl;
				continue;
			}

			int pairid = -1;
			if(r1==0 && r2	==0) pairid = 0; // SS
			else if(r1==0 && r2==1) pairid = 1; // SB
			else if(r1==1 && r2==0) pairid = 2; // BS
			else if(r1==1 && r2==1) pairid = 3; // BB
			
			// Track Spliting & Merging Cut
			// ---------------------------------------------------------------------------
			// A1 = proton from Λ1,  B1 = pion from Λ1
			etaA1[pid1] = dgAeta[pid1][idx1]; phiA1[pid1] = dgAphi[pid1][idx1]; ptA1[pid1] = dgApt[pid1][idx1];
            etaB1[pid1] = dgBeta[pid1][idx1]; phiB1[pid1] = dgBphi[pid1][idx1]; ptB1[pid1] = dgBpt[pid1][idx1];
            etaA2[pid2] = dgAeta[pid2][idx2]; phiA2[pid2] = dgAphi[pid2][idx2]; ptA2[pid2] = dgApt[pid2][idx2];
            etaB2[pid2] = dgBeta[pid2][idx2]; phiB2[pid2] = dgBphi[pid2][idx2]; ptB2[pid2] = dgBpt[pid2][idx2];

			double chargeA1 = (pid1 == 0) ? 1.0 : -1.0;  // 条件式 ? (True時の値) : (False時の値)
			double chargeB1 = (pid1 == 0) ? -1.0 : 1.0;
			double chargeA2 = (pid2 == 0) ? 1.0 : -1.0;
			double chargeB2 = (pid2 == 0) ? -1.0 : 1.0;

			phiStarA1[pid1] = GetPhiStar(phiA1[pid1], ptA1[pid1], chargeA1);
            phiStarA2[pid2] = GetPhiStar(phiA2[pid2], ptA2[pid2], chargeA2);
            phiStarB1[pid1] = GetPhiStar(phiB1[pid1], ptB1[pid1], chargeB1);
            phiStarB2[pid2] = GetPhiStar(phiB2[pid2], ptB2[pid2], chargeB2);

		
			// Here assuming if -999, we skip the cut or force continue (conservative)
			if(phiStarA1[pid1] < -900 || phiStarA2[pid2] < -900 || phiStarB1[pid1] < -900 || phiStarB2[pid2] < -900){
				//cout << "DEBUG: skip by phiStar error" <<endl;
				continue;
				}

		
			if(mSharedTrkFlag == 1){
				// daughter sharing check
				// LL & AA
				if( (fabs(ptA1[pid1]-ptA2[pid2])<1.0e-5 && fabs(etaA1[pid1]-etaA2[pid2])<1.0e-5 && fabs(phiA1[pid1]-phiA2[pid2])<1.0e-5) || //prton
					(fabs(ptB1[pid1]-ptB2[pid2])<1.0e-5 && fabs(etaB1[pid1]-etaB2[pid2])<1.0e-5 && fabs(phiB1[pid1]-phiB2[pid2])<1.0e-5)  ){//pion
						//cout << "DEBUG: skip by dg shearing check" << endl;
						continue;
					}

				// LA
				if( (fabs(ptA1[pid1]-ptB2[pid2])<1.0e-5 && fabs(etaA1[pid1]-etaB2[pid2])<1.0e-5 && fabs(phiA1[pid1]-phiB2[pid2])<1.0e-5) || 
					(fabs(ptB1[pid1]-ptA2[pid2])<1.0e-5 && fabs(etaB1[pid1]-etaA2[pid2])<1.0e-5 && fabs(phiB1[pid1]-phiA2[pid2])<1.0e-5) ){
					continue;
				}
			}

			float angle_same    = GetAngle(PhiDB[pid1][idx1], PhiDB[pid2][idx2], ThetaDB[pid1][idx1], ThetaDB[pid2][idx2]);
			float cosAngle_same = GetcosAngle(PhiDB[pid1][idx1], PhiDB[pid2][idx2], ThetaDB[pid1][idx1], ThetaDB[pid2][idx2]);

			std::vector<int> spinid;
			spinid.push_back(0); // 0: All

			if( cosAngle_same >= 0.0 ) {
				spinid.push_back(1); // 1: 0-90 (Para)
				if( cosAngle_same >= 1.0/sqrt(2) ) spinid.push_back(3); // 3: 0-45
				else spinid.push_back(4); // 4: 45-90
			} else {
				spinid.push_back(2); // 2: 90-180 (Anti)
				if( cosAngle_same >= -1.0/sqrt(2) ) spinid.push_back(5); // 5: 90-135
				else spinid.push_back(6); // 6: 135-180
			}

			// QA 
			// 1. proton - proton
            deta_pp[ipid] = etaA1[pid1] - etaA2[pid2];
            dphi_pp[ipid] = TVector2::Phi_mpi_pi(phiStarA1[pid1] - phiStarA2[pid2]);
			hDeta[ipid][0]->Fill(deta_pp[ipid]);
			hDeta[ipid][0]->Fill(-deta_pp[ipid]);
			hDphi[ipid][0]->Fill(dphi_pp[ipid]); //phi*
			hDphi[ipid][0]->Fill(-dphi_pp[ipid]);

			hDetaDphi1[ipid]->Fill(deta_pp[ipid], phiA1[pid1]-phiA2[pid2]);
			//hDetaDphi1[ipid]->Fill(-deta_pp[ipid], -phiA1[pid1]+phiA2[pid2]);

            // 2. pion - pion
            deta_pipi[ipid] = etaB1[pid1] - etaB2[pid2];
            dphi_pipi[ipid] = TVector2::Phi_mpi_pi(phiStarB1[pid1] - phiStarB2[pid2]);
			hDeta[ipid][1]->Fill(deta_pipi[ipid]);
			hDeta[ipid][1]->Fill(-deta_pipi[ipid]);
			hDphi[ipid][1]->Fill(dphi_pipi[ipid]);
			hDphi[ipid][1]->Fill(-dphi_pipi[ipid]);

			hDetaDphi2[ipid]->Fill(deta_pipi[ipid], phiB1[pid1]-phiB2[pid2]);
			//hDetaDphi2[ipid]->Fill(-deta_pipi[ipid], -phiB1[pid1]+phiB2[pid2]);

            // 3. proton - pion
            deta_ppi[ipid] = etaA1[pid1] - etaB2[pid2];
            dphi_ppi[ipid] = TVector2::Phi_mpi_pi(phiStarA1[pid1] - phiStarB2[pid2]);
			hDeta[ipid][2]->Fill(deta_ppi[ipid]);
			hDphi[ipid][2]->Fill(dphi_ppi[ipid]);

			hDetaDphi3[ipid]->Fill(deta_ppi[ipid], phiA1[pid1]-phiB2[pid2]);

            // 4. pion - proton
            deta_pip[ipid] = etaB1[pid1] - etaA2[pid2];
            dphi_pip[ipid] = TVector2::Phi_mpi_pi(phiStarB1[pid1] - phiStarA2[pid2]);
			hDeta[ipid][3]->Fill(deta_pip[ipid]);
			hDphi[ipid][3]->Fill(dphi_pip[ipid]);

			hDetaDphi3[ipid]->Fill(-deta_pip[ipid], -phiB1[pid1]+phiA2[pid2]);

			hDetaDphiStar1[ipid][0]->Fill(deta_pp[ipid], dphi_pp[ipid]);
			hDetaDphiStar2[ipid][0]->Fill(deta_pipi[ipid], dphi_pipi[ipid]);
			hDetaDphiStar3[ipid][0]->Fill(deta_ppi[ipid], dphi_ppi[ipid]);
			hDetaDphiStar3[ipid][0]->Fill(deta_pip[ipid], dphi_pip[ipid]);
			hDetaDphiStar4[ipid][0]->Fill(deta_pip[ipid], dphi_pip[ipid]);

			if(ipid==0 || ipid==1){
				hDetaDphiStar1[ipid][0]->Fill(-deta_pp[ipid], -dphi_pp[ipid]);
				hDetaDphiStar2[ipid][0]->Fill(-deta_pipi[ipid], -dphi_pipi[ipid]);
			}

			if(cosAngle_same > 0.0){

				hDetaDphiStar1[ipid][1]->Fill(deta_pp[ipid], dphi_pp[ipid]);
				hDetaDphiStar2[ipid][1]->Fill(deta_pipi[ipid], dphi_pipi[ipid]);
				hDetaDphiStar3[ipid][1]->Fill(deta_ppi[ipid], dphi_ppi[ipid]);
				hDetaDphiStar3[ipid][1]->Fill(deta_pip[ipid], dphi_pip[ipid]);
				hDetaDphiStar4[ipid][1]->Fill(deta_pip[ipid], dphi_pip[ipid]);

				if(ipid==0 || ipid==1){
					hDetaDphiStar1[ipid][1]->Fill(-deta_pp[ipid], -dphi_pp[ipid]);
					hDetaDphiStar2[ipid][1]->Fill(-deta_pipi[ipid], -dphi_pipi[ipid]);
				}

			}else if( cosAngle_same < 0.0 ) {

				hDetaDphiStar1[ipid][2]->Fill(deta_pp[ipid], dphi_pp[ipid]);
				hDetaDphiStar2[ipid][2]->Fill(deta_pipi[ipid], dphi_pipi[ipid]);
				hDetaDphiStar3[ipid][2]->Fill(deta_ppi[ipid], dphi_ppi[ipid]);
				hDetaDphiStar3[ipid][2]->Fill(deta_pip[ipid], dphi_pip[ipid]);
				hDetaDphiStar4[ipid][2]->Fill(deta_pip[ipid], dphi_pip[ipid]);

				if(ipid==0 || ipid==1){
					hDetaDphiStar1[ipid][2]->Fill(-deta_pp[ipid], -dphi_pp[ipid]);
					hDetaDphiStar2[ipid][2]->Fill(-deta_pipi[ipid], -dphi_pipi[ipid]);
				}
			} 

			TLorentzVector v1, v2;
			v1 = currentPrts[pid1][idx1];
			v2 = currentPrts[pid2][idx2];
			//cout << endl;
			//cout << "v1 = " << v1.M() << ", v2 = " << v2.M() << endl;
	
			// 	M, Qinv
			float massL = v1.M();// Lambda mass check
			//hmassL->Fill(massL);
		
			TLorentzVector vPair = v1 + v2; //mass 
			float massLL = vPair.M();

			TLorentzVector vDiff = v1 - v2;
			float qinv_same = sqrt(-vDiff.M2()); // .M2() means p^2 - E^2

			if (pairid >= 0 && pairid <= 3) {
				for (int s : spinid) {
					hmass_same_nocut[ipid][pairid][s]->Fill(massLL);
					hqinv_same_nocut[ipid][pairid][s][iCent]->Fill(qinv_same);
					hqinv_same_nocut[ipid][pairid][s][8]->Fill(qinv_same);
				}
			}

			// check　merging
			if( IsPairClose_pp(etaA1[pid1], phiStarA1[pid1], etaA2[pid2], phiStarA2[pid2], ipid) ) continue;// p - p
				for (int s : spinid) {
					hmass_same_cut[0][ipid][pairid][s]->Fill(massLL);
					hqinv_same_cut[0][ipid][pairid][s][iCent]->Fill(qinv_same);
					hqinv_same_cut[0][ipid][pairid][s][8]->Fill(qinv_same);
				}
			
			if( IsPairClose_pipi(etaB1[pid1], phiStarB1[pid1], etaB2[pid2], phiStarB2[pid2], ipid) ) continue; // pi - pi
				 for (int s : spinid) {
					hmass_same_cut[1][ipid][pairid][s]->Fill(massLL);
					hqinv_same_cut[1][ipid][pairid][s][iCent]->Fill(qinv_same);
					hqinv_same_cut[1][ipid][pairid][s][8]->Fill(qinv_same);
				}
		
			if( IsPairClose_ppi(etaA1[pid1], phiStarA1[pid1], etaB2[pid2], phiStarB2[pid2], ipid) ) continue; // p - pi
				 for (int s : spinid) {
					hmass_same_cut[2][ipid][pairid][s]->Fill(massLL);
					hqinv_same_cut[2][ipid][pairid][s][iCent]->Fill(qinv_same);
					hqinv_same_cut[2][ipid][pairid][s][8]->Fill(qinv_same);
				}
			
			if( IsPairClose_pip(etaB1[pid1], phiStarB1[pid1], etaA2[pid2], phiStarA2[pid2], ipid) ) continue; // pi - p 
				 for (int s : spinid) {
					hmass_same_cut[3][ipid][pairid][s]->Fill(massLL);
					hqinv_same_cut[3][ipid][pairid][s][iCent]->Fill(qinv_same);
					hqinv_same_cut[3][ipid][pairid][s][8]->Fill(qinv_same);
				}
/*
			// check Spliting
			if( ipid == 2){
				if( IsSplitCheck(etaA1[pid1], phiStarA1[pid1], etaB2[pid2], phiStarB2[pid2]) ) continue; // p - pi
				if( IsSplitCheck(etaB1[pid1], phiStarB1[pid1], etaA2[pid2], phiStarA2[pid2]) ) continue; // pi - p

			
			if (pairid >= 0 && pairid <= 3) {
				for (int s : spinid) {
					hmass_same_cut[4][ipid][pairid][s]->Fill(massLL);
					hqinv_same_cut[4][ipid][pairid][s][iCent]->Fill(qinv_same);
					hqinv_same_cut[4][ipid][pairid][s][8]->Fill(qinv_same);
					if( qinv_same < 2.0 ) hmass_same_low[ipid][pairid][s]->Fill(massLL);
				}
			}
			}*/


			if(pairid == 0){ //SS
				hcosAngle_same[ipid]->Fill(cosAngle_same);
				hAngle_same[ipid]->Fill(angle_same);
			}

			nPair++;
			//cout << endl;
			} // end of j loop
		} // end of i loop
	} // end of ipid loop



	// ==============================================================================================================
	// Mixing
	// ==============================================================================================================
	for(size_t ipid=0; ipid<3; ipid++){  // ipid=0:LL, 1:AA, 2:LA

		//cout << "=== Mixing Loop start: ipid = " << ipid << " ===" << endl;

		std::vector<std::pair<int, int>> pid_pairs;
		if(ipid == 0) pid_pairs.push_back({0, 0});
		else if(ipid == 1) pid_pairs.push_back({1, 1});
		else {
			// For mixed-event LA pairs use a single ordering to avoid
			// double-counting/from-order asymmetry (process L->A only).
			pid_pairs.push_back({0, 1}); // LA
		}
		
		for(auto mp : pid_pairs){ // autp: pid_pairsの型・ループ回数を自動で取得

			int pid_curr = mp.first;
			int pid_past = mp.second;
			//cout << "pid_curr = " << pid_curr << ", pid_past = " << pid_past << endl;
		
			vector< vector<TLorentzVector> >     &buf      = mixingBuffer[pid_past][iCent][fVz];
			vector< vector<int> >                &regionBuf = mixingregionBuffer[pid_past][iCent][fVz];
			vector< vector<DhVector> > 			 &dgBuf     = mixingDhBuffer[pid_past][iCent][fVz];
			vector< vector<pair<float, float>> > &spinBuf   = mixingSpinBuffer[pid_past][iCent][fVz];

			/*if(buffer.size() > 0 && currentLambdas.size() > 0) {
				cout << "--- Mixing Start for Event " << mEventCounter << " ---" << endl;
				cout << "Current Lambdas: " << currentLambdas.size() << ", Buffer Events: " << NL << endl;
			}*/

			int mixPair = 0;

			for(size_t ievt=0; ievt<buf.size(); ievt++){

				vector<TLorentzVector>	   &pastPrts     = buf[ievt];
				vector<int> 		  	   &pastregionid = regionBuf[ievt];
				vector<DhVector>    	   &pastDgs      = dgBuf[ievt];
				vector<pair<float, float>> &pastSpins    = spinBuf[ievt];

				if( ievt >= spinBuf.size() ) {
					//cout << "DEBUG: ievt >= spinBuf.size()" << endl;
					continue;
				}

				for(size_t i=0; i<currentPrts[pid_curr].size(); i++){
					for(size_t j=0; j<pastPrts.size();j++){

						//--------------------region cut--------------------
						int r1 = regionid[pid_curr][i];
						int r2 = pastregionid[j];
						
						if(r1==-1||r2==-1) continue;

						int pairid = -1;
						if(r1==0 && r2==0) pairid = 0; // SS
						else if(r1==0 && r2==1) pairid = 1; // SB
						else if(r1==1 && r2==0) pairid = 2; // BS
						else if(r1==1 && r2==1) pairid = 3; // BB
						//cout << "pairid = " << pairid << ", r1 = " << r1 << ", r2 = " << r2 << endl;
		

						//-------------Track Spliting & Merging Cut---------------
						
						// past
						float etaP_past  = pastDgs[j].p_eta;  float etaPi_past = pastDgs[j].pi_eta; float phiP_past  = pastDgs[j].p_phi;	
						float phiPi_past = pastDgs[j].pi_phi; float ptP_past   = pastDgs[j].p_pt;   float ptPi_past  = pastDgs[j].pi_pt;

						// calitulate phi star
						double chargeA_curr = (pid_curr == 0 ) ? 1.0 : -1.0; 
						double chargeB_curr = (pid_curr == 0 ) ? -1.0 : 1.0;
						double chargeA_past = ( pid_past == 0) ? 1.0 : -1.0; 
						double chargeB_past = ( pid_past == 0) ? -1.0 : 1.0;

						double phiStarP_curr  = GetPhiStar(dgAphi[pid_curr][i], dgApt[pid_curr][i], chargeA_curr);
						double phiStarP_past  = GetPhiStar(phiP_past, ptP_past, chargeA_past);
						double phiStarPi_curr = GetPhiStar(dgBphi[pid_curr][i], dgBpt[pid_curr][i], chargeB_curr);
						double phiStarPi_past = GetPhiStar(phiPi_past, ptPi_past, chargeB_past);

						// Here assuming if -999, we skip the cut or force continue (conservative)
						if(phiStarP_curr < -900 || phiStarP_past < -900 || phiStarPi_curr < -900 || phiStarPi_past < -900) continue;

						// daughter sharing check
						if(mSharedTrkFlag == 1){
							// LL & AA
							if( (fabs(dgApt[pid_curr][i]-ptP_past)<1.0e-5 && fabs(dgAeta[pid_curr][i]-etaP_past)<1.0e-5 && fabs(dgAphi[pid_curr][i]-phiP_past)<1.0e-5) || //prton
								(fabs(dgBpt[pid_curr][i]-ptPi_past)<1.0e-5 && fabs(dgBeta[pid_curr][i]-etaPi_past)<1.0e-5 && fabs(dgBphi[pid_curr][i]-phiPi_past)<1.0e-5)  )//pion
								continue;

							// LA
							if( (fabs(dgApt[pid_curr][i]-ptPi_past)<1.0e-5 && fabs(dgAeta[pid_curr][i]-etaPi_past)<1.0e-5 && fabs(dgAphi[pid_curr][i]-phiPi_past)<1.0e-5) || 
								(fabs(dgBpt[pid_curr][i]-ptP_past)<1.0e-5 && fabs(dgBeta[pid_curr][i]-etaP_past)<1.0e-5 && fabs(dgBphi[pid_curr][i]-phiP_past)<1.0e-5) ){
								continue;
							}
						}

						float phi1_mix   = PhiDB[pid_curr][i];
						float theta1_mix = ThetaDB[pid_curr][i];

						float phi2_mix   = pastSpins[j].first;
						float theta2_mix = pastSpins[j].second;

						float angle_mix    = GetAngle(phi1_mix, phi2_mix, theta1_mix, theta2_mix);
						float cosAngle_mix = GetcosAngle(phi1_mix, phi2_mix, theta1_mix, theta2_mix);

						std::vector<int> spinid;
						spinid.push_back(0); // 0: All

						if (cosAngle_mix >= 0.0) {
							spinid.push_back(1); // 1: 0-90 (Para)
							if (cosAngle_mix >= 1.0/sqrt(2)) spinid.push_back(3); // 3: 0-45
							else spinid.push_back(4); // 4: 45-90
						} else {
							spinid.push_back(2); // 2: 90-180 (Anti)
							if (cosAngle_mix >= -1.0/sqrt(2)) spinid.push_back(5); // 5: 90-135
							else spinid.push_back(6); // 6: 135-180
						}

						// QA 
						// 1. proton - proton
						float deta_pp_mix = dgAeta[pid_curr][i] - etaP_past;
						float dphi_pp_mix = TVector2::Phi_mpi_pi(phiStarP_curr - phiStarP_past);
						
						hDeta_mix[ipid][0]->Fill(deta_pp_mix);
						hDeta_mix[ipid][0]->Fill(-deta_pp_mix);
						hDphi_mix[ipid][0]->Fill(dphi_pp_mix);
						hDphi_mix[ipid][0]->Fill(-dphi_pp_mix);
						
						hDetaDphi1_mix[ipid]->Fill(deta_pp_mix, dgAphi[pid_curr][i]-phiP_past);
						hDetaDphi1_mix[ipid]->Fill(-deta_pp_mix, -dgAphi[pid_curr][i]+phiP_past);

						// 2. pion - pion
						float deta_pipi_mix = dgBeta[pid_curr][i] - etaPi_past;
						float dphi_pipi_mix = TVector2::Phi_mpi_pi(phiStarPi_curr - phiStarPi_past);
						hDeta_mix[ipid][1]->Fill(deta_pipi_mix);
						hDeta_mix[ipid][1]->Fill(-deta_pipi_mix);
						hDphi_mix[ipid][1]->Fill(dphi_pipi_mix);
						hDphi_mix[ipid][1]->Fill(-dphi_pipi_mix);
						hDetaDphi2_mix[ipid]->Fill(deta_pipi_mix, dgBphi[pid_curr][i]-phiPi_past);
						hDetaDphi2_mix[ipid]->Fill(-deta_pipi_mix, -dgBphi[pid_curr][i]+phiPi_past);

						// 3. proton - pion
						float deta_ppi_mix = dgAeta[pid_curr][i] - etaPi_past;
						float dphi_ppi_mix = TVector2::Phi_mpi_pi(phiStarP_curr - phiStarPi_past);
						hDeta_mix[ipid][2]->Fill(deta_ppi_mix);
						hDphi_mix[ipid][2]->Fill(dphi_ppi_mix);
						hDetaDphi3_mix[ipid]->Fill(deta_ppi_mix, dgAphi[pid_curr][i]-phiPi_past);

						// 4. pion - proton
						float deta_pip_mix = dgBeta[pid_curr][i] - etaP_past;
						float dphi_pip_mix = TVector2::Phi_mpi_pi(phiStarPi_curr - phiStarP_past);
						hDeta_mix[ipid][3]->Fill(deta_pip_mix);
						hDphi_mix[ipid][3]->Fill(dphi_pip_mix);
						hDetaDphi3_mix[ipid]->Fill(-deta_pip_mix, -dgBphi[pid_curr][i]+phiP_past);

						hDetaDphiStar1_mix[ipid][0]->Fill(deta_pp_mix, dphi_pp_mix);
						hDetaDphiStar2_mix[ipid][0]->Fill(deta_pipi_mix, dphi_pipi_mix);
						hDetaDphiStar3_mix[ipid][0]->Fill(deta_ppi_mix, dphi_ppi_mix);
						hDetaDphiStar3_mix[ipid][0]->Fill(deta_pip_mix, dphi_pip_mix);
						hDetaDphiStar4_mix[ipid][0]->Fill(deta_pip_mix, dphi_pip_mix);
						
						if( ipid==0 || ipid==1 ){
							hDetaDphiStar1_mix[ipid][0]->Fill(-deta_pp_mix, -dphi_pp_mix);
							hDetaDphiStar2_mix[ipid][0]->Fill(-deta_pipi_mix, -dphi_pipi_mix);
						}
						
						if(cosAngle_mix > 0.0){

							hDetaDphiStar1_mix[ipid][1]->Fill(deta_pp_mix, dphi_pp_mix);
							hDetaDphiStar2_mix[ipid][1]->Fill(deta_pipi_mix, dphi_pipi_mix);
							hDetaDphiStar3_mix[ipid][1]->Fill(deta_ppi_mix, dphi_ppi_mix);
							hDetaDphiStar3_mix[ipid][1]->Fill(deta_pip_mix, dphi_pip_mix);
							hDetaDphiStar4_mix[ipid][1]->Fill(deta_pip_mix, dphi_pip_mix);

							if( ipid==0 || ipid==1 ){
								hDetaDphiStar1_mix[ipid][1]->Fill(-deta_pp_mix, -dphi_pp_mix);
								hDetaDphiStar2_mix[ipid][1]->Fill(-deta_pipi_mix, -dphi_pipi_mix);
							}


						}else if( cosAngle_mix < 0.0 ) {
							hDetaDphiStar1_mix[ipid][2]->Fill(deta_pp_mix, dphi_pp_mix);
							hDetaDphiStar2_mix[ipid][2]->Fill(deta_pipi_mix, dphi_pipi_mix);
							hDetaDphiStar3_mix[ipid][2]->Fill(deta_ppi_mix, dphi_ppi_mix);
							hDetaDphiStar3_mix[ipid][2]->Fill(deta_pip_mix, dphi_pip_mix);
							hDetaDphiStar4_mix[ipid][2]->Fill(deta_pip_mix, dphi_pip_mix);

							if( ipid==0 || ipid==1 ){
								hDetaDphiStar1_mix[ipid][2]->Fill(-deta_pp_mix, -dphi_pp_mix);
								hDetaDphiStar2_mix[ipid][2]->Fill(-deta_pipi_mix, -dphi_pipi_mix);
							}
						}

						TLorentzVector v_L1, v_L2;
						v_L1 = currentPrts[pid_curr][i];
						v_L2 = pastPrts[j];

						// M, Qinv
						TLorentzVector vPair_mix = v_L1 + v_L2;
						TLorentzVector vDiff_mix = v_L1 - v_L2;

						float massLL_mix = vPair_mix.M();
						float qinv_mix = sqrt(-vDiff_mix.M2()); 

						if (pairid >= 0 && pairid <= 3) { // nocut
							for (int s : spinid) {
								hmass_mix_nocut[ipid][pairid][s]->Fill(massLL_mix);
								hqinv_mix_nocut[ipid][pairid][s][iCent]->Fill(qinv_mix);
								hqinv_mix_nocut[ipid][pairid][s][8]->Fill(qinv_mix);
							}
						}
						
						// previous pair-close check moved earlier to avoid inconsistent fills between spin bins
						if( IsPairClose_pp(dgAeta[pid_curr][i], phiStarP_curr, etaP_past, phiStarP_past, ipid) ) continue; // p - p
							 if (pairid >= 0 && pairid <= 3) {
								for (int s : spinid) {
									hmass_mix_cut[0][ipid][pairid][s]->Fill(massLL_mix);
									hqinv_mix_cut[0][ipid][pairid][s][iCent]->Fill(qinv_mix);
									hqinv_mix_cut[0][ipid][pairid][s][8]->Fill(qinv_mix);
								}
							}
						
						if( IsPairClose_pipi(dgBeta[pid_curr][i], phiStarPi_curr, etaPi_past, phiStarPi_past, ipid) ) continue; // pi - pi
							 if (pairid >= 0 && pairid <= 3) {
								for (int s : spinid) {
									hmass_mix_cut[1][ipid][pairid][s]->Fill(massLL_mix);
									hqinv_mix_cut[1][ipid][pairid][s][iCent]->Fill(qinv_mix);
									hqinv_mix_cut[1][ipid][pairid][s][8]->Fill(qinv_mix);
								}
							}
						
						if( IsPairClose_ppi(dgAeta[pid_curr][i], phiStarP_curr, etaPi_past, phiStarPi_past, ipid) ) continue; // p - pi
							 if (pairid >= 0 && pairid <= 3) {
								for (int s : spinid) {
									hmass_mix_cut[2][ipid][pairid][s]->Fill(massLL_mix);
									hqinv_mix_cut[2][ipid][pairid][s][iCent]->Fill(qinv_mix);
									hqinv_mix_cut[2][ipid][pairid][s][8]->Fill(qinv_mix);
								}
							}
						if( IsPairClose_pip(dgBeta[pid_curr][i], phiStarPi_curr, etaP_past, phiStarP_past, ipid) ) continue; // pi - p		

						double weight = 1.0;
						double weight_BW = 1.0;
                        if ((ipid == 0 || ipid == 1) && hWeightFunc != nullptr && hWeightFunc_BW != nullptr) {
                            if (qinv_mix > 0.0 && qinv_mix < 2.0) { 
                                int bin = hWeightFunc->FindBin(qinv_mix);
								weight = hWeightFunc->GetBinContent(bin);

								int bin_BW = hWeightFunc_BW->FindBin(qinv_mix);
								weight_BW = hWeightFunc_BW->GetBinContent(bin_BW);  
                            }
                        }

							 if (pairid >= 0 && pairid <= 3) {
								for (int s : spinid) {
									hmass_mix_cut[3][ipid][pairid][s]->Fill(massLL_mix);
									hqinv_mix_cut[3][ipid][pairid][s][iCent]->Fill(qinv_mix);
									hqinv_mix_cut[3][ipid][pairid][s][8]->Fill(qinv_mix);

									if( ipid==0 || ipid==1 ){
										if(s==2){
											hmass_mix_weigh[ipid][pairid]->Fill(massLL_mix, weight);
											hqinv_mix_weigh[ipid][pairid]->Fill(qinv_mix, weight);

											hmass_mix_weigh_BW[ipid][pairid]->Fill(massLL_mix, weight_BW);
											hqinv_mix_weigh_BW[ipid][pairid]->Fill(qinv_mix, weight_BW);
										}
									}
								}
							}
						/*
						if(ipid == 2){
							if( IsSplitCheck(dgAeta[pid_curr][i], phiStarP_curr, etaPi_past, phiStarPi_past)) continue; // p - pi
							if( IsSplitCheck(dgBeta[pid_curr][i], phiStarPi_curr, etaP_past, phiStarP_past) ) continue; // pi - p

						if (pairid >= 0 && pairid <= 3) {
							for (int s : spinid) {
								hmass_mix_cut[4][ipid][pairid][s]->Fill(massLL_mix);
								hqinv_mix_cut[4][ipid][pairid][s][iCent]->Fill(qinv_mix);
								hqinv_mix_cut[4][ipid][pairid][s][8]->Fill(qinv_mix);
								if( qinv_mix<2.0 ) hmass_mix_low[ipid][pairid][s]->Fill(massLL_mix);
							}
						}// end of pairid check
						}*/

						if(pairid == 0){ //SS
							hcosAngle_mix[ipid]->Fill(cosAngle_mix);
							hAngle_mix[ipid]->Fill(angle_mix);
						}
						

						/* //Mixing Check
						cout << "[MixCheck] Pair " << mixPair << ": "
							<< "Curr[" << i << "] (Pt=" << v_L1.Pt() << ") " << "[" << iCent << "][" << fVz << "]" // 現在の粒子のIDとPt
							<< "+ BufEvt[" << ievt << "].Trk[" << j << "] (Pt=" << v_L2.Pt() << ") "  // 過去の粒子のIDとPt
							<< "-> Mass=" << vPair_mix.M() // 計算された質量
							<< endl;*/

						//mixPair ++;
					}//end of jTrack
				}//end of iTrack
			}//end of ievt
		}// end of pid_curr_pairs
	}//end of Mixing


	for(size_t ipid=0; ipid<3; ipid++){
		vector< vector<TLorentzVector> >     &buf       = mixingBuffer[ipid][iCent][fVz];
		vector< vector<int> >                &regionBuf = mixingregionBuffer[ipid][iCent][fVz];
		vector< vector<DhVector> >           &dgBuf     = mixingDhBuffer[ipid][iCent][fVz];
		vector< vector<pair<float, float>> > &spinBuf   = mixingSpinBuffer[ipid][iCent][fVz];
	
		bool ValidPrts = false;
		for(size_t k=0; k<regionid[ipid].size() ; k++){
			if(regionid[ipid][k] != -1){
				ValidPrts = true;
				break;
			}
		}
		
		if(ValidPrts){
			buf.push_back(currentPrts[ipid]);
			regionBuf.push_back(regionid[ipid]);

			vector<pair<float, float>> currentSpins;
			currentSpins.reserve(PhiDB[ipid].size());
			for(size_t k=0; k<PhiDB[ipid].size(); k++){
				currentSpins.push_back( make_pair(PhiDB[ipid][k], ThetaDB[ipid][k]) );
			}
			dgBuf.push_back(currentDhs[ipid]);
			spinBuf.push_back(currentSpins);
		}

		if((int)buf.size()>N_MIX){
			buf.erase(buf.begin());
			regionBuf.erase(regionBuf.begin());
			//cout << "----Erased-----" <<endl;

			dgBuf.erase(dgBuf.begin());
			spinBuf.erase(spinBuf.begin());
		}

	} // end of pid

//===================================
	// reset
	for( int ich=0; ich<2; ich++ ){
		prtA[ich].clear();
		prtB[ich].clear();
		prtAtpc[ich].clear();
		prtBtpc[ich].clear();

		trkid[ich].clear(); //+
		mass[ich].clear();
		dgApt[ich].clear(); dgAeta[ich].clear(); dgAphi[ich].clear();
		dgBpt[ich].clear(); dgBeta[ich].clear(); dgBphi[ich].clear();
	}

    return 0;
}


//______________________________________________________________________________
//void CorTreeAnalyzer::ResetSglTrk( sgl_trk *trk ){
//
//	trk->pt  = 0.0;
//	trk->eta = 0.0;
//	trk->phi = 0.0;
//	trk->ch  = 0;
//}

int CorTreeAnalyzer::GetZvtxBin( float zvtx ){

	int fZvtx = static_cast<int>( ( zvtx + 30.0 ) / 60.0 * NQvZvtx );
	return fZvtx;
}

bool CorTreeAnalyzer::passTPcut(float &dcaA, float &dcaB, float &dcaV0, float &dL, float &dcaAB, int fCent, int fcutID, int id){
	// fcutID : tof hit info. (ignoring for tree analysis)
	// id=0 -> default cut 
	// id=1 -> strict cut for direct profile
	// id=2 -> fixed cut

	bool isV0 = false;
	if( dcaA>mDcaP[fCent][id] && dcaB>mDcaPi[fCent][id] && dcaV0<mDcaV0[fCent][id] && dL>mDL[fCent][id] && dcaAB<mDcaDh[fCent][id] ) isV0 = true;
	return isV0;
}

void CorTreeAnalyzer::SetChAsymWidthRun11(){

	// run11 Ach RMS
	mAcWd[0] = 0.0317639; // temporary set for 0-5%
	mAcWd[1] = 0.0317639; // temporary set for 5-10%
	mAcWd[2] = 0.0317639;
	mAcWd[3] = 0.0373613;
	mAcWd[4] = 0.0456029;
	mAcWd[5] = 0.0579191;
	mAcWd[6] = 0.0765909;
	mAcWd[7] = 0.1054216; 
	mAcWd[8] = 0.1511892; 

	//mAkWd[0] = 0.174804; // temporary set for 0-5%
	//mAkWd[1] = 0.174804; // temporary set for 5-10%
	//mAkWd[2] = 0.174804;
	//mAkWd[3] = 0.208053;
	//mAkWd[4] = 0.256363;
	//mAkWd[5] = 0.324206;
	//mAkWd[6] = 0.411659;
	//mAkWd[7] = 0.505960; 
	//mAkWd[8] = 0.524032; 
}

void CorTreeAnalyzer::DefHistograms(){

	hVz  = new TH1F("hVz","hVz",200,-35,35);
	hRM  = new TH1F("hRM","hRM",300,0,600);
	hRMw = new TH1F("hRMw","hRMw",300,0,600);
	hRMc  = new TH1F("hRMc","hRMc",300,0,600);
	hRMcw = new TH1F("hRMcw","hRMcw",300,0,600);
	hCentBin  = new TH1F("hCentBin","hCentBin",16,0,16);
	hCentBinw = new TH1F("hCentBinw","hCentBinw",16,0,16);
	hCentBin2  = new TH1F("hCentBin2","hCentBin2",16,0,16);
	hCentBinw2 = new TH1F("hCentBinw2","hCentBinw2",16,0,16);
	hCentIso8 = new TH1D("hCentIso8","hCentIso8",8,0,8);

	const char* pidName2[2]  = {"L", "A"};
	const char* pidName[3]  = {"LL", "AA", "LA"};
    const char* pidTitle[3] = {"#Lambda#Lambda", "#bar{#Lambda}#bar{#Lambda}", "#Lambda#bar{#Lambda}"};

	const char* dgPairName[4] = {"pp", "pipi", "ppi", "pip"};
    
    const char* pairName[4]  = {"SS", "SB", "BS", "BB"};
    const char* pairTitle[4] = {"Sig-Sig", "Sig-BG", "BG-Sig", "BG-BG"};
    
    const char* spinName[3]   = {"All", "Para", "Anti"};
	const char* spinName7[7]  = {"All", "0-90", "90-180", "0-45", "45-90", "90-135", "135-180"};
    //const char* spinTitle[3]  = {"All-Spin", "Parallel", "Anti-Parallel",};
    const char* spinTitle7[7] = {"All-Spin", "0-90", "90-180", "0-45", "45-90", "90-135", "135-180"};

	const char* centName[9]  = {"0-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "All"};
	//const char* centTitle[9] = {"0-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "All"};

	int bin_m = 200;
	int bin_q = 250;

	for(int ipid =0; ipid<2; ipid++){

		hpt[ipid] = new TH1F(Form("hpt_%s", pidName2[ipid]), Form("hpt_%s", pidName2[ipid]), 200, 0.0, 10.0);
		heta[ipid] = new TH1F(Form("heta_%s", pidName2[ipid]), Form("heta_%s", pidName2[ipid]), 100, -1.5, 1.5);
		hphi[ipid] = new TH1F(Form("hphi_%s", pidName2[ipid]), Form("hphi_%s", pidName2[ipid]), 100, -3.14159, 3.14159);

		hinvM_BG[ipid] = new TH1D(Form("hinvM_BG_%s", pidName2[ipid]), Form("hinvM_BG_%s", pidName2[ipid]), 200, 1.05, 1.2);
		hinvM_Sig[ipid] = new TH1D(Form("hinvM_Sig_%s", pidName2[ipid]), Form("hinvM_Sig_%s", pidName2[ipid]), 200, 1.05, 1.2);

		for(int icent = 0; icent < 9; icent++) {
			hinvM[ipid][icent] = new TH1D(Form("hinvM_%s_%s", pidName2[ipid], centName[icent]), Form("hinvM_%s_%s", pidName2[ipid], centName[icent]), 200, 1.05, 1.2);
		}
	}

	for(int ipid = 0; ipid < 3; ipid++) {

		hcosAngle_same[ipid] = new TH1D(Form("hcosAngle_same_%s", pidName[ipid]), Form("hcosAngle_same_%s", pidName[ipid]), 100, -1.0, 1.0);
		hAngle_same[ipid] 	 = new TH1D(Form("hAngle_same_%s", pidName[ipid]), Form("hAngle_same_%s", pidName[ipid]), 100, -180.0, 180.0);
		hcosAngle_mix[ipid]  = new TH1D(Form("hcosAngle_mix_%s", pidName[ipid]), Form("hcosAngle_mix_%s", pidName[ipid]), 100, -1.0, 1.0);
		hAngle_mix[ipid] 	 = new TH1D(Form("hAngle_mix_%s", pidName[ipid]), Form("hAngle_mix_%s", pidName[ipid]), 100, -180.0, 180.0);


		for(int idgpair = 0; idgpair < 4; idgpair++){
			hDeta[ipid][idgpair] 	 = new TH1D(Form("h1D_Deta_%s_%s", pidName[ipid], dgPairName[idgpair]), Form("h1D_Deta_%s_%s", pidName[ipid], dgPairName[idgpair]), 100, -1.5, 1.5);
			hDphi[ipid][idgpair] 	 = new TH1D(Form("h1D_Dphi_%s_%s", pidName[ipid], dgPairName[idgpair]), Form("h1D_Dphi_%s_%s", pidName[ipid], dgPairName[idgpair]), 100, -1.5, 1.5);
			
			hDeta_mix[ipid][idgpair] = new TH1D(Form("h1D_Deta_mix_%s_%s", pidName[ipid], dgPairName[idgpair]), Form("h1D_Deta_mix_%s_%s", pidName[ipid], dgPairName[idgpair]), 100, -1.5, 1.5);
			hDphi_mix[ipid][idgpair] = new TH1D(Form("h1D_Dphi_mix_%s_%s", pidName[ipid], dgPairName[idgpair]), Form("h1D_Dphi_mix_%s_%s", pidName[ipid], dgPairName[idgpair]), 100, -1.5, 1.5);
		}

		for(int ispin = 0; ispin < 3; ispin++) {
			
			hDetaDphiStar1[ipid][ispin] = new TH2D(Form("hDetaDphiStar1_%s_%s", pidName[ipid], spinName[ispin]), Form("hDetaDphiStar1_%s_%s_%s", dgPairName[0], pidName[ipid], spinName[ispin]), 100, -0.5, 0.5, 100, -0.5, 0.5); //100, -0.5, 0.5(Kehao)
			hDetaDphiStar2[ipid][ispin] = new TH2D(Form("hDetaDphiStar2_%s_%s", pidName[ipid], spinName[ispin]), Form("hDetaDphiStar2_%s_%s_%s", dgPairName[1], pidName[ipid], spinName[ispin]), 100, -0.5, 0.5, 100, -0.5, 0.5);
			hDetaDphiStar3[ipid][ispin] = new TH2D(Form("hDetaDphiStar3_%s_%s", pidName[ipid], spinName[ispin]), Form("hDetaDphiStar3_%s_%s_%s", dgPairName[2], pidName[ipid], spinName[ispin]), 100, -0.5, 0.5, 100, -0.5, 0.5);
			hDetaDphiStar4[ipid][ispin] = new TH2D(Form("hDetaDphiStar4_%s_%s", pidName[ipid], spinName[ispin]), Form("hDetaDphiStar3_%s_%s_%s", dgPairName[3], pidName[ipid], spinName[ispin]), 100, -0.5, 0.5, 100, -0.5, 0.5);
			
			hDetaDphiStar1_mix[ipid][ispin] = new TH2D(Form("hDetaDphiStar1_mix_%s_%s", pidName[ipid], spinName[ispin]), Form("hDetaDphiStar_mix_%s_%s_%s", dgPairName[0], pidName[ipid], spinName[ispin]), 100, -0.5, 0.5, 100, -0.5, 0.5);
			hDetaDphiStar2_mix[ipid][ispin] = new TH2D(Form("hDetaDphiStar2_mix_%s_%s", pidName[ipid], spinName[ispin]), Form("hDetaDphiStar_mix_%s_%s_%s", dgPairName[1], pidName[ipid], spinName[ispin]), 100, -0.5, 0.5, 100, -0.5, 0.5);
			hDetaDphiStar3_mix[ipid][ispin] = new TH2D(Form("hDetaDphiStar3_mix_%s_%s", pidName[ipid], spinName[ispin]), Form("hDetaDphiStar_mix_%s_%s_%s", dgPairName[2], pidName[ipid], spinName[ispin]), 100, -0.5, 0.5, 100, -0.5, 0.5);
			hDetaDphiStar4_mix[ipid][ispin] = new TH2D(Form("hDetaDphiStar4_mix_%s_%s", pidName[ipid], spinName[ispin]), Form("hDetaDphiStar_mix_%s_%s_%s", dgPairName[3], pidName[ipid], spinName[ispin]), 100, -0.5, 0.5, 100, -0.5, 0.5);
	
		}

		hDetaDphi1[ipid]     = new TH2D(Form("hDetaDphi1_%s", pidName[ipid]), Form("hDetaDphi1_%s", pidName[ipid]), 100, -1.5, 1.5, 100, -1.5, 1.5);
		hDetaDphi2[ipid]     = new TH2D(Form("hDetaDphi2_%s", pidName[ipid]), Form("hDetaDphi2_%s", pidName[ipid]), 100, -1.5, 1.5, 100, -1.5, 1.5);
		hDetaDphi3[ipid]     = new TH2D(Form("hDetaDphi3_%s", pidName[ipid]), Form("hDetaDphi3_%s", pidName[ipid]), 100, -1.5, 1.5, 100, -1.5, 1.5);

		hDetaDphi1_mix[ipid]     = new TH2D(Form("hDetaDphi1_mix_%s", pidName[ipid]), Form("hDetaDphi1_mix_%s", pidName[ipid]), 100, -1.5, 1.5, 100, -1.5, 1.5);
		hDetaDphi2_mix[ipid]     = new TH2D(Form("hDetaDphi2_mix_%s", pidName[ipid]), Form("hDetaDphi2_mix_%s", pidName[ipid]), 100, -1.5, 1.5, 100, -1.5, 1.5);
		hDetaDphi3_mix[ipid]     = new TH2D(Form("hDetaDphi3_mix_%s", pidName[ipid]), Form("hDetaDphi3_mix_%s", pidName[ipid]), 100, -1.5, 1.5, 100, -1.5, 1.5);

		for(int ipair = 0; ipair < 4; ipair++) {
            for(int ispin = 0; ispin < 7; ispin++) {
				for(int icut = 0; icut< 5; icut++){
					TString nameMassS_cut = Form("hmass_same_%d_%s_%s_%s", icut, pidName[ipid], pairName[ipair], spinName7[ispin]);
                	TString titleMassS_cut = Form("Same %d %s %s %s;M_{inv} (GeV/c^{2});Counts", icut, pidTitle[ipid], pairTitle[ipair], spinTitle7[ispin]);
					hmass_same_cut[icut][ipid][ipair][ispin] = new TH1D(nameMassS_cut, titleMassS_cut, 100, 2.23, 2.28);

					TString nameMassM_cut = Form("hmass_mix_%d_%s_%s_%s", icut, pidName[ipid], pairName[ipair], spinName7[ispin]);
                	TString titleMassM_cut = Form("Mix %d %s %s %s;M_{inv} (GeV/c^{2});Counts", icut, pidTitle[ipid], pairTitle[ipair], spinTitle7[ispin]);
					hmass_mix_cut[icut][ipid][ipair][ispin] = new TH1D(nameMassM_cut, titleMassM_cut, 100, 2.23, 2.28);
				}

				// === Same Event ===
				TString nameMassS = Form("hmass_same_%s_%s_%s", pidName[ipid], pairName[ipair], spinName7[ispin]);
                TString titleMassS = Form("Same %s %s %s;M_{inv} (GeV/c^{2});Counts", pidTitle[ipid], pairTitle[ipair], spinTitle7[ispin]);
				hmass_same_nocut[ipid][ipair][ispin] = new TH1D(nameMassS+"nocut", titleMassS+"nocut", bin_m, 1.0, 5.5);

                // === Mixed Event ===
                TString nameMassM = Form("hmass_mix_%s_%s_%s", pidName[ipid], pairName[ipair], spinName7[ispin]);
                TString titleMassM = Form("Mix %s %s %s;M_{inv} (GeV/c^{2});Counts", pidTitle[ipid], pairTitle[ipair], spinTitle7[ispin]);
				hmass_mix_nocut[ipid][ipair][ispin] = new TH1D(nameMassM+"nocut", titleMassM+"nocut", bin_m, 1.0, 5.5);

				//hmass_same_zoom[ipid][ipair][ispin] = new TH1F(nameMassS, titleMassS, 200, 0.1, 0.5);
				//hmass_mix_zoom[ipid][ipair][ispin] = new TH1F(nameMassM, titleMassM, 200, 0.1, 0.5);
				TString nameMassS_low = Form("hmass_same_low_%s_%s_%s", pidName[ipid], pairName[ipair], spinName7[ispin]);
				TString titleMassS_low = Form("Same Low %s %s %s;M_{inv} (GeV/c^{2});Counts", pidTitle[ipid], pairTitle[ipair], spinTitle7[ispin]);
				hmass_same_low[ipid][ipair][ispin] = new TH1D(nameMassS_low, titleMassS_low, bin_m, 2.22, 2.6);

				TString nameMassM_low = Form("hmass_mix_low_%s_%s_%s", pidName[ipid], pairName[ipair], spinName7[ispin]);
				TString titleMassM_low = Form("Mix Low %s %s %s;M_{inv} (GeV/c^{2});Counts", pidTitle[ipid], pairTitle[ipair], spinTitle7[ispin]);
				hmass_mix_low[ipid][ipair][ispin] = new TH1D(nameMassM_low, titleMassM_low, bin_m, 2.22, 2.6);
            
				for(int icent = 0; icent < 9; icent++){
					for(int icut = 0; icut< 5; icut++){
						TString nameQinvS_cut = Form("hqinv_same_%d_%s_%s_%s_%s", icut, pidName[ipid], pairName[ipair], spinName7[ispin], centName[icent]);
						TString titleQinvS_cut = Form("Same %d %s %s %s %s;Q_{inv} (GeV/c);Counts", icut, pidTitle[ipid], pairTitle[ipair], spinTitle7[ispin], centName[icent]);
						hqinv_same_cut[icut][ipid][ipair][ispin][icent] = new TH1D(nameQinvS_cut, titleQinvS_cut, bin_q, 0.0, 4.0);

						TString nameQinvM_cut = Form("hqinv_mix_%d_%s_%s_%s_%s", icut, pidName[ipid], pairName[ipair], spinName7[ispin], centName[icent]);
                		TString titleQinvM_cut = Form("Mix %d %s %s %s %s;Q_{inv} (GeV/c);Counts", icut, pidTitle[ipid], pairTitle[ipair], spinTitle7[ispin], centName[icent]);
                		hqinv_mix_cut[icut][ipid][ipair][ispin][icent] = new TH1D(nameQinvM_cut, titleQinvM_cut, bin_q, 0.0, 4.0);
					}
					
					TString nameQinvS_nocut = Form("hqinv_same_%s_%s_%s_%s", pidName[ipid], pairName[ipair], spinName7[ispin], centName[icent]);
					TString titleQinvS_nocut = Form("Same %s %s %s %s;Q_{inv} (GeV/c);Counts", pidTitle[ipid], pairTitle[ipair], spinTitle7[ispin], centName[icent]);
					hqinv_same_nocut[ipid][ipair][ispin][icent] = new TH1D(nameQinvS_nocut+"nocut", titleQinvS_nocut+"nocut", bin_q, 0.0, 4.0);

					TString nameQinvM_nocut = Form("hqinv_mix_%s_%s_%s_%s", pidName[ipid], pairName[ipair], spinName7[ispin], centName[icent]);
                	TString titleQinvM_nocut = Form("Mix %s %s %s %s;Q_{inv} (GeV/c);Counts", pidTitle[ipid], pairTitle[ipair], spinTitle7[ispin], centName[icent]);
					hqinv_mix_nocut[ipid][ipair][ispin][icent] = new TH1D(nameQinvM_nocut+"nocut", titleQinvM_nocut+"nocut", bin_q, 0.0, 4.0);

				} // icent roop
			}// ispin loop

			hmass_mix_weigh[ipid][ipair] = new TH1D(Form("hmass_mix_weigh_%s_%s", pidName[ipid], pairName[ipair]), Form("Mix (weigh) %s %s;M_{inv} (GeV/c^{2});Counts", pidTitle[ipid], pairTitle[ipair]), 100, 2.23, 2.28);
			hqinv_mix_weigh[ipid][ipair] = new TH1D(Form("hqinv_mix_weigh_%s_%s", pidName[ipid], pairName[ipair]), Form("Mix (weigh) %s %s;Q_{inv} (GeV/c);Counts", pidTitle[ipid], pairTitle[ipair]), bin_q, 0.0, 4.0);

			hmass_mix_weigh_BW[ipid][ipair] = new TH1D(Form("hmass_mix_weigh_BW_%s_%s", pidName[ipid], pairName[ipair]), Form("Mix (weigh, BW) %s %s;M_{inv} (GeV/c^{2});Counts", pidTitle[ipid], pairTitle[ipair]), 100, 2.23, 2.28);
			hqinv_mix_weigh_BW[ipid][ipair] = new TH1D(Form("hqinv_mix_weigh_BW_%s_%s", pidName[ipid], pairName[ipair]), Form("Mix (weigh, BW) %s %s;Q_{inv} (GeV/c);Counts", pidTitle[ipid], pairTitle[ipair]), bin_q, 0.0, 4.0);
		} // ipair loop
	} // ipid loop 

	//Double_t achbins[Nasm+1] = { -4.0, -0.86, -0.26, 0.26, 0.86, 4.0 };
	//hAchBin  = new TProfile( Form("hAchBin"), Form("hAchBin"), Nasm, achbins, -10, 10 );
	//hAchpBin = new TProfile( Form("hAchpBin"), Form("hAchpBin"), Nasm, achbins, -10, 10 );
	//for( int k=0; k<4; k++ ){
	//	if( k<2 ) hAchCent[k] = new TH2F( Form("hAchCent_%d",k), Form("hAchCent_%d",k), 8, 0, 8, 50, -0.5, 0.5 );
	//	else      hAchCent[k] = new TH2F( Form("hAchCent_%d",k), Form("hAchCent_%d",k), 8, 0, 8, 50, -6, 6 );
	//}

	for( int ih=0; ih<2; ih++ ){
		for( int ics=0; ics<2; ics++ ){ 
			hPsiCor[ih][ics] = new TProfile( Form("hPsiCor_%d_%d",ih,ics), Form("hPsiCor_%d_%d",ih,ics), 9, 0, 9 );
			hPsiCor[ih][ics]->Sumw2();
		}
	}

	for( int ih=0; ih<2; ih++ ){
		for( int ics=0; ics<2; ics++ ){
			hZDCPsiCor[ih*2+ics] = new TProfile(Form("hZDCPsiCor_%d",ih*2+ics),Form("hZDCPsiCor_%d",ih*2+ics),9,0,9);
		}
	}
	for( int ih=0; ih<6; ih++ ) hTPCPsiCor[ih] = new TProfile(Form("hTPCPsiCor_%d",ih),Form("hTPCPsiCor_%d",ih),9,0,9);
	for( int ih=0; ih<2; ih++ ) hBBCPsiCor[ih] = new TProfile(Form("hBBCPsiCor_%d",ih),Form("hBBCPsiCor_%d",ih),9,0,9);
	for( int ih=0; ih<3; ih++ )	hEPDPsiCor[ih] = new TProfile(Form("hEPDPsiCor_%d",ih),Form("hEPDPsiCor_%d",ih),9,0,9);
	for( int k=0; k<2; k++ )	hEPDZDCPsi1Cor[k] = new TProfile(Form("hEPDZDCPsi1Cor_%d",k),Form("hEPDZDCPsi1Cor_%d",k),9,0,9);
	
	// centrality dep. of Py
	// InvM
	for( int ich=0; ich<2; ich++ ){
		for( int iep=0; iep<NsubSMD; iep++ ){
			for( int ict=0; ict<NcentPolZ; ict++ ){  // 0-80 with 10% step
				for( int iph=0; iph<NphiPol; iph++ ){
					hInvM_LamCent[ich][iep][ict][iph] = new TH1D( Form("hInvM_LamCent_%d_%d_%d_%d",ich,iep,ict,iph), Form("hInvM_LamCent_%d_%d_%d_%d",ich,iep,ict,iph), 200, 1.07, 1.17 );
					hInvM_LamCent[ich][iep][ict][iph]->Sumw2();
	}	}	}	}

	// IM-method
	for( int ich=0; ich<2; ich++ ){
		for( int iep=0; iep<NsubSMD; iep++ ){
			for( int ict=0; ict<NcentPolZ; ict++ ){
				hInvMvsPolCent  [ich][iep][ict] = new TProfile( Form("hInvMvsPolCent_%d_%d_%d",ich,iep,ict), Form("hInvMvsPolCent_%d_%d_%d",ich,iep,ict), 100, 1.07, 1.17, -10, 10 );
				hInvMvsPolCent  [ich][iep][ict]->Sumw2();
	}	}	}

	// direct profile
	for( int ich=0; ich<2; ich++ ){
		for( int iep=0; iep<NsubSMD; iep++ ){
			hPyCent[ich][iep] = new TProfile( Form("hPyCent_%d_%d",ich,iep), Form("hPyCent_%d_%d",ich,iep), 8, 0, 80, -10, 10 );
			hPyCent[ich][iep]->Sumw2();
		}
	}

	// pT dep. of Py
	for( int ich=0; ich<2; ich++ ){
		for( int ict=0; ict<2; ict++ ){  // 0-20, 20-60
			for( int ipt=0; ipt<NptPol; ipt++ ){
				hIM_Pt[ich][ict][ipt] = new TH1D( Form("hIM_Pt_%d_%d_%d",ich,ict,ipt), Form("hIM_Pt_%d_%d_%d",ich,ict,ipt), 200, 1.07, 1.17 );
				hIM_Pt[ich][ict][ipt]->Sumw2();

				hIMvsPy_Pt[ich][ict][ipt] = new TProfile( Form("hIMvsPy_Pt_%d_%d_%d",ich,ict,ipt), Form("hIMvsPy_Pt_%d_%d_%d",ich,ict,ipt), 200, 1.07, 1.17 );
				hIMvsPy_Pt[ich][ict][ipt]->Sumw2();
			}
			hA0pt[ich][ict] = new TProfile( Form("hA0pt_%d_%d",ich,ict), Form("hA0_pt_%d_%d",ich,ict), NptPol, 0, NptPol );
			hA0pt[ich][ict]->Sumw2();
		}
	}


	// Ach study
	//for( int ich=0; ich<2; ich++ ){
	//	for( int iep=0; iep<NsubSMD; iep++ ){
	//		for( int ias=0; ias<Nasm; ias++ ){
	//			hInvMvsPolAch [ich][iep][ias] = new TProfile( Form("hInvMvsPolAch_%d_%d_%d",ich,iep,ias), Form("hInvMvsPolAch_%d_%d_%d",ich,iep,ias), 100, 1.07, 1.17, -10, 10 );
	//			hInvMvsPolAchp[ich][iep][ias] = new TProfile( Form("hInvMvsPolAchp_%d_%d_%d",ich,iep,ias), Form("hInvMvsPolAchp_%d_%d_%d",ich,iep,ias), 100, 1.07, 1.17, -10, 10 );
	//			hInvMvsPolAch [ich][iep][ias]->Sumw2();
	//			hInvMvsPolAchp[ich][iep][ias]->Sumw2();

	//			hInvMAch [ich][iep][ias] = new TH1F( Form("hInvMAch_%d_%d_%d",ich,iep,ias), Form("hInvMAch_%d_%d_%d",ich,iep,ias), 100, 1.07, 1.17 );
	//			hInvMAchp[ich][iep][ias] = new TH1F( Form("hInvMAchp_%d_%d_%d",ich,iep,ias), Form("hInvMAchp_%d_%d_%d",ich,iep,ias), 100, 1.07, 1.17 );
	//			hInvMAch [ich][iep][ias]->Sumw2();
	//			hInvMAchp[ich][iep][ias]->Sumw2();
	//}	}	}


	// z-component
	for( int ich=0; ich<2; ich++ ){
		for( int iep=0; iep<NpolRP; iep++ ){
			for( int ixyz=0; ixyz<3; ixyz++ ){
				for( int iph=0; iph<NphiRP; iph++ ){
					hInvMvsPolXYZ[ich][iep][ixyz][iph] = new TProfile( Form("hInvMvsPolXYZ_%d_%d_%d_%d",ich,iep,ixyz,iph), Form("hInvMvsPolXYZ_%d_%d_%d_%d",ich,iep,ixyz,iph), 100, 1.07, 1.17, -10, 10 );
					hInvMvsPolXYZ[ich][iep][ixyz][iph]->Sumw2();
	}	}	}	}

	for( int ich=0; ich<2; ich++ ){
		for( int iep=0; iep<NpolRP; iep++ ){
			for( int iph=0; iph<NphiRP; iph++ ){
				hInvMXYZ[ich][iep][iph] = new TH1D( Form("hInvMXYZ_%d_%d_%d",ich,iep,iph), Form("hInvMXYZ_%d_%d_%d",ich,iep,iph), 100, 1.07, 1.17 );
				hInvMXYZ[ich][iep][iph]->Sumw2();
	}	}	}

	// eta 
	for( int ich=0; ich<2; ich++ ){
		for( int iep=0; iep<NpolRP; iep++ ){
			for( int ieta=0; ieta<6; ieta++ ){
				for( int iph=0; iph<NphiRP; iph++ ){
					hIMvsPzEta[ich][iep][ieta][iph] = new TProfile( Form("hIMvsPzEta_%d_%d_%d_%d",ich,iep,ieta,iph), Form("hIMvsPzEta_%d_%d_%d_%d",ich,iep,ieta,iph), 100, 1.07, 1.17, -10, 10 );
					hInvMEta  [ich][iep][ieta][iph] = new TH1D( Form("hInvMEta_%d_%d_%d_%d",ich,iep,ieta,iph), Form("hInvMEta_%d_%d_%d_%d",ich,iep,ieta,iph), 100, 1.07, 1.17 );
					hIMvsPzEta[ich][iep][ieta][iph]->Sumw2();
					hInvMEta  [ich][iep][ieta][iph]->Sumw2();
	}	}	}	}

	// higher order
	for( int ih=0; ih<2; ih++ ){
		for( int ich=0; ich<2; ich++ ){
			for( int irp=0; irp<3; irp++ ){
				for( int ict=0; ict<2; ict++ ){
					for( int iph=0; iph<NphiRP; iph++ ){
						hIMvsPzN[ih][ich][irp][ict][iph] = new TProfile( Form("hIMvsPzN_%d_%d_%d_%d_%d",ih,ich,irp,ict,iph), Form("hIMvsPzN_%d_%d_%d_%d_%d",ih,ich,irp,ict,iph), 100, 1.07, 1.17, -10, 10 );
						hIMPzN  [ih][ich][irp][ict][iph] = new TH1D( Form("hIMPzN_%d_%d_%d_%d_%d",ih,ich,irp,ict,iph), Form("hIMPzN_%d_%d_%d_%d_%d",ih,ich,irp,ict,iph), 100, 1.07, 1.17 );
						hIMvsPzN[ih][ich][irp][ict][iph]->Sumw2();
						hIMPzN  [ih][ich][irp][ict][iph]->Sumw2();
	}	}	}	}	}

	for( int ih=0; ih<2; ih++ ){
		for( int ich=0; ich<2; ich++ ){
			for( int ict=0; ict<NcentPolZ; ict++ ){
				for( int iph=0; iph<NphiRP; iph++ ){
					hIMvsPzNCent[ih][ich][ict][iph] = new TProfile( Form("hIMvsPzNCent_%d_%d_%d_%d",ih,ich,ict,iph), Form("hIMvsPzNCent_%d_%d_%d_%d",ih,ich,ict,iph), 100, 1.07, 1.17, -10, 10 );
					hIMvsPzNCent[ih][ich][ict][iph]->Sumw2();
	}	}	}	}

	for( int ih=0; ih<2; ih++ ){
		for( int ich=0; ich<2; ich++ ){
			for( int ict=0; ict<2; ict++ ){
				for( int ipt=0; ipt<NptPol; ipt++ ){
					for( int iph=0; iph<NphiRP; iph++ ){
						hIMvsPzNPtDph[ih][ich][ict][ipt][iph] = new TProfile( Form("hIMvsPzNPtDph_%d_%d_%d_%d_%d",ih,ich,ict,ipt,iph), Form("hIMvsPzNPtDph_%d_%d_%d_%d_%d",ih,ich,ict,ipt,iph), 100, 1.07, 1.17, -10, 10 );
						hIMvsPzNPtDph[ih][ich][ict][ipt][iph]->Sumw2();
	}	}	}	}	}

	// profile
	for( int ich=0; ich<2; ich++ ){
		for( int irp=0; irp<NpolRP; irp++ ){
			for( int isc=0; isc<2; isc++ ){
				//hPzCt[ich][irp][isc] = new TProfile( Form("hPzCt_%d_%d_%d",ich,irp,isc), Form("hPzCt_%d_%d_%d",ich,irp,isc), 11, 0, 110 );
				//hPzPt[ich][irp][isc] = new TProfile( Form("hPzPt_%d_%d_%d",ich,irp,isc), Form("hPzPt_%d_%d_%d",ich,irp,isc), 18, 0.5, 5.9 );
				//hPzEt[ich][irp][isc] = new TProfile( Form("hPzEt_%d_%d_%d",ich,irp,isc), Form("hPzEt_%d_%d_%d",ich,irp,isc), 6, -1, 1 );
				//hPzCt[ich][irp][isc]->Sumw2();
				//hPzPt[ich][irp][isc]->Sumw2();
				//hPzEt[ich][irp][isc]->Sumw2();

				for( int ict=0; ict<NcentPolZ; ict++ ){ 
					hIMvsPz2Ct[ich][irp][ict][isc] = new TProfile( Form("hIMvsPz2Ct_%d_%d_%d_%d",ich,irp,ict,isc), Form("hIMvsPz2Ct_%d_%d_%d_%d",ich,irp,ict,isc), 100, 1.07, 1.17, -10, 10 );
					hIMvsPz2Ct[ich][irp][ict][isc]->Sumw2();
				}
				for( int ieta=0; ieta<6; ieta++ ){ 
					hIMvsPz2Et[ich][irp][ieta][isc] = new TProfile( Form("hIMvsPz2Et_%d_%d_%d_%d",ich,irp,ieta,isc), Form("hIMvsPz2Et_%d_%d_%d_%d",ich,irp,ieta,isc), 100, 1.07, 1.17, -10, 10 );
					hIMvsPz2Et[ich][irp][ieta][isc]->Sumw2();
				}
	}	}	}

	for( int ich=0; ich<2; ich++ ){
		for( int ipt=0; ipt<NptPol; ipt++ ){
			hIMvsPz2PtEP[ich][ipt] = new TProfile( Form("hIMvsPz2PtEP_%d_%d",ich,ipt), Form("hIMvsPz2Et_%d_%d",ich,ipt), 100, 1.07, 1.17, -10, 10 );
			hIMvsPz2PtEP[ich][ipt]->Sumw2();
		}
	}


	// added on 9/27/2022
	for( int ih=0; ih<3; ih++ ){
		for( int ich=0; ich<2; ich++ ){
			for( int iw=0; iw<2; iw++ ){
				hPzNprof2D[ih][ich][iw] = new TProfile2D( Form("hPzNprof2D_%d_%d_%d",ih,ich,iw), Form("hPzNprof2D_%d_%d_%d",ih,ich,iw), NptPol, 0, NptPol, NcentPolZ , 0, NcentPolZ, -10, 10 );
				hPzNprof2D[ih][ich][iw]->Sumw2();
			}
		}
	}
	for( int ich=0; ich<2; ich++ ){ 
		for( int iw=0; iw<2; iw++ ){
			hCos2prof2D[ich][iw] = new TProfile2D( Form("hCos2prof2D_%d_%d",ich,iw), Form("hCos2prof2D_%d_%d",ich,iw), NptPol, 0, NptPol, NcentPolZ, 0, NcentPolZ, -10, 10 );
			hCos2prof2D[ich][iw]->Sumw2();

			hIMvsPtCent[ich][iw] = new TH3D( Form("hIMvsPtCent_%d_%d",ich,iw), Form("hIMvsPtCent_%d_%d",ich,iw), 100, 1.07, 1.17, NptPol, 0, NptPol, NcentPolZ, 0, NcentPolZ );
			hIMvsPtCent[ich][iw]->Sumw2();
		}
		hIMvsPtCentIso[ich] = new TH3D( Form("hIMvsPtCentIso_%d",ich), Form("hIMvsPtCentIso_%d",ich), 100, 1.07, 1.17, 10, 0.5, 2.5, NcentPolZ, 0, NcentPolZ );
		hIMvsPtCentIso[ich]->Sumw2();
	}



	// off-peak check
	for( int ich=0; ich<2; ich++ ){
		for( int irp=0; irp<NpolRP; irp++ ){
			for( int im=0; im<2; im++ ){
				for( int isc=0; isc<2; isc++ ){
					hBgCt[ich][irp][im][isc] = new TProfile( Form("hBgCt_%d_%d_%d_%d",ich,irp,im,isc), Form("hBgCt_%d_%d_%d_%d",ich,irp,im,isc), 10, 0, 100 );
					hBgEt[ich][irp][im][isc] = new TProfile( Form("hBgEt_%d_%d_%d_%d",ich,irp,im,isc), Form("hBgEt_%d_%d_%d_%d",ich,irp,im,isc), 6, -1, 1 );
					hBgCt[ich][irp][im][isc]->Sumw2();
					hBgEt[ich][irp][im][isc]->Sumw2();
	}	}	}	}
	for( int ih=0; ih<3; ih++ ){
		for( int ich=0; ich<2; ich++ ){
			for( int im=0; im<2; im++ ){
				hBgCtLg[ih][ich][im] = new TProfile( Form("hBgCtLg_%d_%d_%d",ih,ich,im), Form("hBgCtLg_%d_%d_%d",ih,ich,im), 10, 0, 100 );
				hBgCtLg[ih][ich][im]->Sumw2();
			}
		}
	}

	// higher harmonic
	for( int ih=0; ih<2; ih++ ){
		for( int ich=0; ich<2; ich++ ){
			for( int irp=0; irp<3; irp++ ){
				for( int isc=0; isc<2; isc++ ){
					hPzNCt[ih][ich][irp][isc] = new TProfile( Form("hPzNCt_%d_%d_%d_%d",ih,ich,irp,isc), Form("hPzNCt_%d_%d_%d_%d",ih,ich,irp,isc), 11, 0, 110 );
					hPzNCt[ih][ich][irp][isc]->Sumw2();
				}
	}	}	}

	for( int ih=0; ih<2; ih++ ){
		for( int ich=0; ich<2; ich++ ){
			for( int ict=0; ict<NcentPolZ; ict++ ){ 
				hIMvsPzNCt[ih][ich][ict] = new TProfile( Form("hIMvsPzNCt_%d_%d_%d",ih,ich,ict), Form("hIMvsPzNCt_%d_%d_%d",ih,ich,ict), 100, 1.07, 1.17, -10, 10 );
				hIMvsPzNCt[ih][ich][ict]->Sumw2();
	}	}	}



	// acc study
	for( int ich=0; ich<2; ich++ ){
		for( int iep=0; iep<NpolRP; iep++ ){
			for( int iz=0; iz<2; iz++ ){
				hInvMvsPzVz[ich][iep][iz] = new TProfile( Form("hInvMvsPzVz_%d_%d_%d",ich,iep,iz), Form("hInvMvsPzVz_%d_%d_%d",ich,iep,iz), 100, 1.07, 1.17, -10, 10 );
				hInvMvsPzVz[ich][iep][iz]->Sumw2();

				for( int isc=0; isc<2; isc++ ){
					hPzCtVz[ich][iep][iz][isc] = new TProfile( Form("hPzCtVz_%d_%d_%d_%d",ich,iep,iz,isc), Form("hPzCtVz_%d_%d_%d_%d",ich,iep,iz,isc), 11, 0, 110 );
					hPzCtVz[ich][iep][iz][isc]->Sumw2();
				}
	}	}	}


	// EP-method check
	for( int ich=0; ich<2; ich++ ){
		for( int ict=0; ict<NcentPolZ; ict++ ){
			hOffIMvsDphi[ich][ict] = new TH1F( Form("hOffIMvsDphi_%d_%d",ich,ict), Form("hOffIMvsDphi_%d_%d",ich,ict), NphiPol, -TMath::Pi(), TMath::Pi() );
		}
	}

	// v2 test
	for( int ich=0; ich<2; ich++ ){
		for( int ict=0; ict<2; ict++ ){
			hLamV2pt[ich][ict] = new TProfile( Form("hLamV2pt_%d_%d",ich,ict), Form("hLamV2pt_%d_%d",ich,ict), 10, 0.5, 5.5, -10, 10 );
		}
		hLamV2ct[ich] = new TProfile( Form("hLamV2ct_%d",ich), Form("hLamV2ct_%d",ich), 8, 0, 80, -10, 10 );

		for( int ipt=0; ipt<5; ipt++ ){
			//hOnIMvsDphiPt[ich][ipt] = new TH1F( Form("hOnIMvsDphiPt_%d_%d",ich,ipt), Form("hOnIMvsDphiPt_%d_%d",ich,ipt), 10, -TMath::Pi(), TMath::Pi() );
			hOnIMvsDphiPt[ich][ipt] = new TH1D( Form("hOnIMvsDphiPt_%d_%d",ich,ipt), Form("hOnIMvsDphiPt_%d_%d",ich,ipt), 6, -TMath::Pi(), TMath::Pi() );
		}
		for( int ict=0; ict<NcentPolZ; ict++ ){
			//hOnIMvsDphiCt[ich][ict] = new TH1F( Form("hOnIMvsDphiCt_%d_%d",ich,ict), Form("hOnIMvsDphiCt_%d_%d",ich,ict), 10, -TMath::Pi(), TMath::Pi() );
			hOnIMvsDphiCt[ich][ict] = new TH1D( Form("hOnIMvsDphiCt_%d_%d",ich,ict), Form("hOnIMvsDphiCt_%d_%d",ich,ict), 6, -TMath::Pi(), TMath::Pi() );
		}

		hCosLamPPt[ich] = new TProfile( Form("hCosLamPPt_%d",ich), Form("hCosLamPPt_%d",ich), 10, 1, 6, -10, 10 );
		hCosLamPCt[ich] = new TProfile( Form("hCosLamPCt_%d",ich), Form("hCosLamPCt_%d",ich), 8, 0, 80, -10, 10 );
	}

	// sys study
	if( mSysFlag ){
		for( int ich=0; ich<2; ich++ ){
			for( int ict=0; ict<2; ict++ ){
				for( int iph=0; iph<NphiRP; iph++ ){
					hIMPzSys  [ich][ict][iph] = new TH2D( Form("hIMPzSys_%d_%d_%d",ich,ict,iph), Form("hIMPzSys_%d_%d_%d",ich,ict,iph), 100, 1.07, 1.17, 11, 0, 11 );
					hIMPzSys  [ich][ict][iph]->Sumw2();
					for(int ih=0; ih<3; ih++){
						hIMvsPzNSys[ih][ich][ict][iph] = new TProfile2D( Form("hIMvsPzSys_%d_%d_%d_%d",ih,ich,ict,iph), Form("hIMvsPzSys_%d_%d_%d_%d",ih,ich,ict,iph), 100, 1.07, 1.17, 11, 0, 11 );
						hIMvsPzNSys[ih][ich][ict][iph]->Sumw2();
					}
				}
			}
		}

		for( int ich=0; ich<2; ich++ ){
			for( int ict=0; ict<NcentPolZ; ict++ ){
				hIMctSys     [ich][ict] = new TH2D( Form("hIMctSys_%d_%d",ich,ict), Form("hIMctSys_%d_%d",ich,ict), 100, 1.07, 1.17, 11, 0, 11 );
				hIMctSys     [ich][ict]->Sumw2();

				for(int ih=0; ih<3; ih++){
					hIMvsPzNctSys[ih][ich][ict] = new TProfile2D( Form("hIMvsPz2ctSys_%d_%d_%d",ih,ich,ict), Form("hIMvsPz2ctSys_%d_%d_%d",ih,ich,ict), 100, 1.07, 1.17, 11, 0, 11 );
					hIMvsPzNctSys[ih][ich][ict]->Sumw2();
				}
			}
		}

		for( int ich=0; ich<2; ich++ ){
			for( int ict=0; ict<2; ict++ ){
				for( int ipt=0; ipt<NptPol; ipt++ ){
					hIMptSys     [ich][ict][ipt] = new TH2D( Form("hIMptSys_%d_%d_%d",ich,ict,ipt), Form("hIMptSys_%d_%d_%d",ich,ict,ipt), 100, 1.07, 1.17, 11, 0, 11 );
					hIMptSys     [ich][ict][ipt]->Sumw2();

					for(int ih=0; ih<3; ih++){ 
						hIMvsPzNptSys[ih][ich][ict][ipt] = new TProfile2D( Form("hIMvsPz2ptSys_%d_%d_%d_%d",ih,ich,ict,ipt), Form("hIMvsPz2ptSys_%d_%d_%d_%d",ih,ich,ict,ipt), 100, 1.07, 1.17, 11, 0, 11 );
						hIMvsPzNptSys[ih][ich][ict][ipt]->Sumw2();
					}
				}
			}
		}
	}

	for( int ich=0; ich<2; ich++ ){
		for( int ict=0; ict<2; ict++ ){
			hAPTvsPtBin[ich][ict] = new TProfile( Form("hAPTvsPtBin_%d_%d",ich,ict), Form("hAPTvsPtBin_%d_%d",ich,ict), NptPol, 0, NptPol );
			hAPTvsPtBin[ich][ict]->Sumw2();
		}
		hAPTvsCent[ich] = new TProfile( Form("hAPTvsCent_%d",ich), Form("hAPTvsCent_%d",ich), 8, 0, 80 );
	}

	// Random EP  
	for( int ich=0; ich<2; ich++ ){
		for( int ict=0; ict<NcentPolZ; ict++ ){ 
			hIMvsPz2CtRnd[ich][ict] = new TProfile( Form("hIMvsPz2CtRnd_%d_%d",ich,ict), Form("hIMvsPz2CtRnd_%d_%d",ich,ict), 100, 1.07, 1.17, -10, 10 );
			hIMvsPz2CtRnd[ich][ict]->Sumw2();
		}
		for( int ipt=0; ipt<NptPol; ipt++ ){
			hIMvsPz2PtRnd[ich][ipt] = new TProfile( Form("hIMvsPz2PtRnd_%d_%d",ich,ipt), Form("hIMvsPz2PtRnd_%d_%d",ich,ipt), 100, 1.07, 1.17, -10, 10 );
			hIMvsPz2PtRnd[ich][ipt]->Sumw2();
		}
	}


	for( int ich=0; ich<2; ich++ ){
		for( int ict=0; ict<NcentPolZ; ict++ ){
			for( int iph=0; iph<NphiRP; iph++ ){
				hIMcent    [ich][ict][iph] = new TH1D( Form("hIMcent_%d_%d_%d",ich,ict,iph), Form("hIMcent_%d_%d_%d",ich,ict,iph), 100, 1.07, 1.17 );
				hIMvsPzCent[ich][ict][iph] = new TProfile( Form("hIMvsPzCent_%d_%d_%d",ich,ict,iph), Form("hIMvsPzCent_%d_%d_%d",ich,ict,iph), 100, 1.07, 1.17 );
				hIMcent    [ich][ict][iph]->Sumw2();
				hIMvsPzCent[ich][ict][iph]->Sumw2();
			}
		}

		for( int ict=0; ict<2; ict++ ){
			for( int ipt=0; ipt<NptPol; ipt++ ){
				for( int iph=0; iph<NphiRP; iph++ ){
					hIMpt    [ich][ict][ipt][iph] = new TH1D( Form("hIMpt_%d_%d_%d_%d",ich,ict,ipt,iph), Form("hIMpt_%d_%d_%d_%d",ich,ict,ipt,iph), 100, 1.07, 1.17 );
					hIMvsPzPt[ich][ict][ipt][iph] = new TProfile( Form("hIMvsPzPt_%d_%d_%d_%d",ich,ict,ipt,iph), Form("hIMvsPzPt_%d_%d_%d_%d",ich,ict,ipt,iph), 100, 1.07, 1.17 );
					hIMpt    [ich][ict][ipt][iph]->Sumw2();
					hIMvsPzPt[ich][ict][ipt][iph]->Sumw2();
				}
			}
	}	}


	// acceptance correction factor
	for( int ich=0; ich<2; ich++ ){
		hCos2cent[ich] = new TProfile( Form("hCos2cent_%d",ich), Form("hCos2cent_%d",ich), 8, 0, 8 );
		hCos2cent[ich]->Sumw2();

		hCos2eta[ich] = new TProfile( Form("hCos2eta_%d",ich), Form("hCos2eta_%d",ich), 10, -1, 1 );
		hCos2eta[ich]->Sumw2();

		for( int ict=0; ict<2; ict++ ){
			hCos2pt[ich][ict] = new TProfile( Form("hCos2pt_%d_%d",ich,ict), Form("hCos2pt_%d_%d",ich,ict), NptPol, 0, NptPol );
			hCos2pt[ich][ict]->Sumw2();
		}
	}

	// check on Py phi dependence
	for( int ich=0; ich<2; ich++ ){
		for( int irp=0; irp<4; irp++ ){
			for( int iph=0; iph<NepPol; iph++ ){
				hIM_EP  [ich][irp][iph] = new TH1D( Form("hIM_EP_%d_%d_%d",ich,irp,iph), Form("hIM_EP_%d_%d_%d",ich,irp,iph), 200, 1.07, 1.17 );
				hIMvsPy_EP[ich][irp][iph] = new TProfile( Form("hIMvsPy_EP_%d_%d_%d",ich,irp,iph), Form("hIMvsPy_EP_%d_%d_%d",ich,irp,iph), 200, 1.07, 1.17 );
				hIM_EP    [ich][irp][iph]->Sumw2();
				hIMvsPy_EP[ich][irp][iph]->Sumw2();
			}
		}
	}



	// EP calibration
	for( int ih=0; ih<2; ih++ ){
		for( int k=0; k<6; k++ ){
			for( int ict=0; ict<16; ict++ ){
				for( int iz=0; iz<12; iz++ ){
					hQvRec[ih][k][ict][iz] = new TProfile( Form("hQvRec_%d_%d_%d_%d",ih,k,ict,iz), Form("hQvRec_%d_%d_%d_%d",ih,k,ict,iz), 1550, -0.5, 1550-0.5, -9000, 9000, "s" );
				}
		   }
	   }
	}

	for( int ih=0; ih<2; ih++ ){
		for( int isub=0; isub<3; isub++ ){
			for( int ict=0; ict<8; ict++ ){
				for( int icr=0; icr<2; icr++ ){
					hPsi[ih][isub][ict][icr] = new TH1F( Form("hPsi_%d_%d_%d_%d",ih,isub,ict,icr), Form("hPsi_%d_%d_%d_%d",ih,isub,ict,icr), 100, -TMath::Pi()/(ih+1.0), TMath::Pi()/(ih+1.0) );
	}	}	}	}


	// test on self-corr.
	for( int ich=0; ich<2; ich++ ){
		for( int ict=0; ict<NcentPolZ; ict++ ){ 
			for( int ic=0; ic<2; ic++ ){
				hIMvsPz2CtSC[ich][ict][ic] = new TProfile( Form("hIMvsPz2CtSC_%d_%d_%d",ich,ict,ic), Form("hIMvsPz2CtSC_%d_%d_%d",ich,ict,ic), 100, 1.07, 1.17, -10, 10 );
				hIMvsPz2CtSC[ich][ict][ic]->Sumw2();
	}	}	}
	for( int ics=0; ics<2; ics++ ){
		hPsi2CorMod[ics] = new TProfile( Form("hPsi2CorMod_%d",ics), Form("hPsi2CorMod_%d",ics), 9, 0, 9 );
		hPsi2CorMod[ics]->Sumw2();
	}
	for( int ich=0; ich<2; ich++ ){
		for( int ipt=0; ipt<NptPol; ipt++ ){ 
			for( int ic=0; ic<2; ic++ ){
				hIMvsPz2PtSC[ich][ipt][ic] = new TProfile( Form("hIMvsPz2PtSC_%d_%d_%d",ich,ipt,ic), Form("hIMvsPz2PtSC_%d_%d_%d",ich,ipt,ic), 100, 1.07, 1.17, -10, 10 );
				hIMvsPz2PtSC[ich][ipt][ic]->Sumw2();
	}	}	}
	
	for( int ich=0; ich<2; ich++ ){
		for( int iet=0; iet<6; iet++ ){ 
			hIMetSC[ich][iet] = new TH1D( Form("hIMetSC_%d_%d",ich,iet), Form("hIMetSC_%d_%d",ich,iet), 100, 1.07, 1.17 ); 
			for( int ic=0; ic<3; ic++ ){
				hIMvsPz2etSC[ich][iet][ic] = new TProfile( Form("hIMvsPz2etSC_%d_%d_%d",ich,iet,ic), Form("hIMvsPz2etSC_%d_%d_%d",ich,iet,ic), 100, 1.07, 1.17, -10, 10 );
			}
		}
	}



	for( int ich=0; ich<2; ich++ ){
		hIMvsVz_Pz[ich] = new TProfile2D( Form("hIMvsVz_Pz_%d",ich), Form("hIMvsVz_Pz_%d",ich), 100, 1.07, 1.17, 12, -6, 6 );
		hIMvsCut_Pz[ich] = new TProfile2D( Form("hIMvsCut_Pz_%d",ich), Form("hIMvsCut_Pz_%d",ich), 100, 1.07, 1.17, 8, 0, 8 );
	}

	// tof matching
	for( int ich=0; ich<2; ich++ ){
		for( int ict=0; ict<NcentPolZ; ict++ ){
			for( int id=0; id<4; id++ ){
				hIMvsPtv0[ich][ict][id] = new TH2F( Form("hIMvsPtv0_%d_%d_%d",ich,ict,id), Form("hIMvsPtv0_%d_%d_%d",ich,ict,id), 100, 1.07, 1.17, 40, 0.5, 4.5 );
	}	}	}

	// slefcor study
	for( int ich=0; ich<2; ich++ ){
		for( int k=0; k<8; k++ ){
			hPz2SCvsCent[ich][k] = new TProfile( Form("hPz2SCvsCent_%d_%d",ich,k), Form("hPz2SCvsCent_%d_%d",ich,k), 8, 0, 8, -10, 10 );
			hPz2SCvsCent[ich][k]->Sumw2();
	}	}



}

void CorTreeAnalyzer::WriteHistograms(){

	hVz ->Write();
	hRM ->Write();
	hRMw->Write();
	hRMc ->Write();
	hRMcw->Write();
	hCentBin ->Write();
	hCentBinw->Write();
	hCentBin2 ->Write();
	hCentBinw2->Write();
	hCentIso8 ->Write();
	
	for(int ipid=0; ipid < 2; ipid++){
		heta[ipid]->Write();
		hphi[ipid]->Write();
		hpt[ipid]->Write();

		hinvM_BG[ipid]->Write();
		hinvM_Sig[ipid]->Write();

		for(int icent = 0; icent < 9; icent++){
			hinvM[ipid][icent]->Write();
		}
	}

	for(int ipid = 0; ipid < 3; ipid++) {

		hcosAngle_same[ipid]->Write();
		hAngle_same[ipid]->Write();
		hcosAngle_mix[ipid]->Write();
		hAngle_mix[ipid]->Write();

		for(int idgpair = 0; idgpair < 4; idgpair++){
			hDeta[ipid][idgpair]->Write();
			hDphi[ipid][idgpair]->Write();

			hDeta_mix[ipid][idgpair]->Write();
			hDphi_mix[ipid][idgpair]->Write();
		}

		for(int ispin = 0; ispin < 3; ispin++) {
		
			hDetaDphiStar1[ipid][ispin]->Write();
			hDetaDphiStar2[ipid][ispin]->Write();
			hDetaDphiStar3[ipid][ispin]->Write();
			hDetaDphiStar4[ipid][ispin]->Write();

			hDetaDphiStar1_mix[ipid][ispin]->Write();
			hDetaDphiStar2_mix[ipid][ispin]->Write();
			hDetaDphiStar3_mix[ipid][ispin]->Write();
			hDetaDphiStar4_mix[ipid][ispin]->Write();
		}
		
		hDetaDphi1[ipid]->Write();
		hDetaDphi2[ipid]->Write();
		hDetaDphi3[ipid]->Write();

		hDetaDphi1_mix[ipid]->Write();
		hDetaDphi2_mix[ipid]->Write();
		hDetaDphi3_mix[ipid]->Write();

		for(int ipair = 0; ipair < 4; ipair++) {
            for(int ispin = 0; ispin < 7; ispin++) {	// All, 0<θ<90, 90<θ<180, 0<θ<45, 45<θ<90, 90<θ<135, 135<θ<180
				for(int icut=0; icut < 5; icut++){
					hmass_same_cut[icut][ipid][ipair][ispin]->Write();
					hmass_mix_cut[icut][ipid][ipair][ispin]->Write();
				}
				hmass_same_nocut[ipid][ipair][ispin]->Write();
				hmass_mix_nocut[ipid][ipair][ispin]->Write();

				//hqinv_same_zoom[ipid][ipair][ispin]->Write();
				//hqinv_mix_zoom[ipid][ipair][ispin]->Write();

				hmass_same_low[ipid][ipair][ispin]->Write();
				hmass_mix_low[ipid][ipair][ispin]->Write();

				for(int icent = 0; icent < 9; icent++){
					for(int icut=0; icut<5; icut++){
						hqinv_same_cut[icut][ipid][ipair][ispin][icent]->Write();
						hqinv_mix_cut[icut][ipid][ipair][ispin][icent]->Write();
					}
					hqinv_same_nocut[ipid][ipair][ispin][icent]->Write();
					hqinv_mix_nocut[ipid][ipair][ispin][icent]->Write();
					
				} // icent loop
			} // ispin loop

			if( ipid == 0 || ipid == 1 ){
                hmass_mix_weigh[ipid][ipair]->Write();
                hqinv_mix_weigh[ipid][ipair]->Write();
				hmass_mix_weigh_BW[ipid][ipair]->Write();
                hqinv_mix_weigh_BW[ipid][ipair]->Write();
            }
		} // ipair loop
	} // pid loop

	/*// EP calibration
	if( mCalibMode==1 ){ 
		for( int ih=0; ih<2; ih++ ){
			for( int k=0; k<6; k++ ){
				for( int ict=0; ict<16; ict++ ){
					for( int iz=0; iz<12; iz++ ){
						hQvRec[ih][k][ict][iz]->Write();
		}	}	}	}
	}
	if( mCalibMode==1 ) return;

	for( int ih=0; ih<2; ih++ ){
		for( int isub=0; isub<3; isub++ ){
			for( int ict=0; ict<8; ict++ ){
				for( int icr=0; icr<2; icr++ ){
					hPsi[ih][isub][ict][icr]->Write();
	}	}	}	}

	for( int ih=0; ih<2; ih++ ){
		for( int ics=0; ics<2; ics++ ) hPsiCor[ih][ics]->Write();
	}



	for( int ih=0; ih<2; ih++ ){
		for( int ics=0; ics<2; ics++ ){
			hZDCPsiCor[ih*2+ics]->Write();
		}
	}
	for( int ih=0; ih<6; ih++ ) hTPCPsiCor[ih]->Write();
	for( int ih=0; ih<2; ih++ ) hBBCPsiCor[ih]->Write();
	for( int ih=0; ih<3; ih++ ) hEPDPsiCor[ih]->Write();
	for( int k=0; k<2; k++ )	hEPDZDCPsi1Cor[k]->Write();*/


	//hAchBin ->Write();
	//hAchpBin->Write();
	//for( int k=0; k<4; k++ ) hAchCent[k]->Write();


	// centrality dep. of Py
	// InvMass dist. 
	for( int ich=0; ich<2; ich++ ){
		for( int iep=0; iep<NsubSMD; iep++ ){
			for( int ict=0; ict<NcentPolZ; ict++ ){  // 10-80 with 10% step
				for( int iph=0; iph<NphiPol; iph++ ){
					hInvM_LamCent[ich][iep][ict][iph]->Write();
	}	}	}	}

	// IM-method
	for( int ich=0; ich<2; ich++ ){
		for( int iep=0; iep<NsubSMD; iep++ ){
			for( int ict=0; ict<NcentPolZ; ict++ ){
				hInvMvsPolCent  [ich][iep][ict]->Write();
				//hInvMvsPolCentCr[ich][iep][ict]->Write();
	}	}	}

	// direct profile
	for( int ich=0; ich<2; ich++ ){
		for( int iep=0; iep<NsubSMD; iep++ ){
			hPyCent[ich][iep]->Write();
		}
	}

	// pT dep. of Py
	for( int ich=0; ich<2; ich++ ){
		for( int ict=0; ict<2; ict++ ){  // 0-20, 20-60
			for( int ipt=0; ipt<NptPol; ipt++ ){
				hIM_Pt    [ich][ict][ipt]->Write();
				hIMvsPy_Pt[ich][ict][ipt]->Write();
			}
			hA0pt[ich][ict]->Write();
		}
	}

	// ach study
	//for( int ich=0; ich<2; ich++ ){
	//	for( int iep=0; iep<NsubSMD; iep++ ){
	//		for( int ias=0; ias<Nasm; ias++ ){
	//			hInvMvsPolAch [ich][iep][ias]->Write();
	//			hInvMvsPolAchp[ich][iep][ias]->Write();

	//			hInvMAch [ich][iep][ias]->Write();
	//			hInvMAchp[ich][iep][ias]->Write();
	//}	}	}

	// z-component
	//for( int ich=0; ich<2; ich++ ){
	//	for( int iep=0; iep<NpolRP; iep++ ){
	//		for( int ixyz=0; ixyz<3; ixyz++ ){
	//			for( int iph=0; iph<NphiRP; iph++ ){
	//				hInvMvsPolXYZ[ich][iep][ixyz][iph]->Write();
	//}	}	}	}

	//for( int ich=0; ich<2; ich++ ){
	//	for( int iep=0; iep<NpolRP; iep++ ){
	//		for( int iph=0; iph<NphiRP; iph++ ){
	//			hInvMXYZ[ich][iep][iph]->Write();
	//}	}	}

	//// eta 
	for( int ich=0; ich<2; ich++ ){
		for( int iep=0; iep<NpolRP; iep++ ){
			for( int ieta=0; ieta<6; ieta++ ){
				for( int iph=0; iph<NphiRP; iph++ ){
					hIMvsPzEta[ich][iep][ieta][iph]->Write();
					hInvMEta[ich][iep][ieta][iph]->Write();
	}	}	}	}

	//// higher order
	//for( int ih=0; ih<2; ih++ ){
	//	for( int ich=0; ich<2; ich++ ){
	//		for( int irp=0; irp<3; irp++ ){
	//			for( int ict=0; ict<2; ict++ ){
	//				for( int iph=0; iph<NphiRP; iph++ ){
	//					hIMvsPzN[ih][ich][irp][ict][iph]->Write();
	//					hIMPzN  [ih][ich][irp][ict][iph]->Write();
	//}	}	}	}	}
	//for( int ih=0; ih<2; ih++ ){
	//	for( int ich=0; ich<2; ich++ ){
	//		for( int ict=0; ict<NcentPolZ; ict++ ){
	//			for( int iph=0; iph<NphiRP; iph++ ){
	//				hIMvsPzNCent[ih][ich][ict][iph]->Write();
	//}	}	}	}

	//for( int ih=0; ih<2; ih++ ){
	//	for( int ich=0; ich<2; ich++ ){
	//		for( int ict=0; ict<2; ict++ ){
	//			for( int ipt=0; ipt<NptPol; ipt++ ){
	//				for( int iph=0; iph<NphiRP; iph++ ){
	//					hIMvsPzNPtDph[ih][ich][ict][ipt][iph]->Write();
	//}	}	}	}	}

	// profile
	for( int ich=0; ich<2; ich++ ){
		for( int irp=0; irp<NpolRP; irp++ ){
			for( int isc=0; isc<2; isc++ ){
				//hPzCt[ich][irp][isc]->Write();
				//hPzPt[ich][irp][isc]->Write();
				//hPzEt[ich][irp][isc]->Write();
				for( int ict=0; ict<NcentPolZ; ict++ ) hIMvsPz2Ct[ich][irp][ict ][isc]->Write();
				for( int ieta=0; ieta<6; ieta++ )      hIMvsPz2Et[ich][irp][ieta][isc]->Write();
			}
	}	}
	for( int ich=0; ich<2; ich++ ){
		for( int ipt=0; ipt<NptPol; ipt++ ){
			hIMvsPz2PtEP[ich][ipt]->Write();
		}
	}

	// 2D profile
	for( int ih=0; ih<3; ih++ ){
		for( int ich=0; ich<2; ich++ ){
			for( int iw=0; iw<2; iw++ ){
				hPzNprof2D[ih][ich][iw]->Write();
			}
		}
	}
	for( int ich=0; ich<2; ich++ ){ 
		for( int iw=0; iw<2; iw++ ){
			hCos2prof2D[ich][iw]->Write();
			hIMvsPtCent[ich][iw]->Write();
		}
		hIMvsPtCentIso[ich]->Write();
	}



	// off-peak check
	for( int ich=0; ich<2; ich++ ){
		for( int irp=0; irp<NpolRP; irp++ ){
			for( int im=0; im<2; im++ ){
				for( int isc=0; isc<2; isc++ ){
					hBgCt[ich][irp][im][isc]->Write();
					hBgEt[ich][irp][im][isc]->Write();
	}	}	}	}
	for( int ih=0; ih<3; ih++ ){
		for( int ich=0; ich<2; ich++ ){
			for( int im=0; im<2; im++ ){
				hBgCtLg[ih][ich][im]->Write();
			}
		}
	}

	// higher harmonic
	for( int ih=0; ih<2; ih++ ){
		for( int ich=0; ich<2; ich++ ){
			for( int irp=0; irp<3; irp++ ){
				for( int isc=0; isc<2; isc++ ){
					hPzNCt[ih][ich][irp][isc]->Write();
				}
	}	}	}
	for( int ih=0; ih<2; ih++ ){
		for( int ich=0; ich<2; ich++ ){
			for( int ict=0; ict<NcentPolZ; ict++ ) hIMvsPzNCt[ih][ich][ict]->Write();
		}
	}


	// acc study
	for( int ich=0; ich<2; ich++ ){
		for( int iep=0; iep<NpolRP; iep++ ){
			for( int iz=0; iz<2; iz++ ){
				hInvMvsPzVz[ich][iep][iz]->Write();
				for( int isc=0; isc<2; isc++ ){
					hPzCtVz[ich][iep][iz][isc]->Write();
				}
	}	}	}

	// EP-method check
	for( int ich=0; ich<2; ich++ ){
		for( int ict=0; ict<NcentPolZ; ict++ ){
			hOffIMvsDphi[ich][ict]->Write();
		}
	}

	// v2 test
	for( int ich=0; ich<2; ich++ ){
		for( int ict=0; ict<2; ict++ ){
			hLamV2pt[ich][ict]->Write();
		}
		hLamV2ct[ich]->Write();

		for( int ipt=0; ipt<5; ipt++ )         hOnIMvsDphiPt[ich][ipt]->Write();
		for( int ict=0; ict<NcentPolZ; ict++ ) hOnIMvsDphiCt[ich][ict]->Write();
		
		hCosLamPPt[ich]->Write();
		hCosLamPCt[ich]->Write();
	}

	// sys study
	if( mSysFlag ){
		for( int ich=0; ich<2; ich++ ){
			for( int ict=0; ict<2; ict++ ){
				for( int iph=0; iph<NphiRP; iph++ ){
					hIMPzSys  [ich][ict][iph]->Write();
					for(int ih=0; ih<3; ih++) hIMvsPzNSys[ih][ich][ict][iph]->Write();
				}
			}
		}

		for( int ich=0; ich<2; ich++ ){
			for( int ict=0; ict<NcentPolZ; ict++ ){
				hIMctSys     [ich][ict]->Write();
				for(int ih=0; ih<3; ih++) hIMvsPzNctSys[ih][ich][ict]->Write();
			}
		}

		for( int ich=0; ich<2; ich++ ){
			for( int ict=0; ict<2; ict++ ){
				for( int ipt=0; ipt<NptPol; ipt++ ){
					hIMptSys     [ich][ict][ipt]->Write();
					for(int ih=0; ih<3; ih++) hIMvsPzNptSys[ih][ich][ict][ipt]->Write();
				}
			}
		}
	}// mSysFlag

	for( int ich=0; ich<2; ich++ ){ 
		for( int ict=0; ict<2; ict++ ){
			hAPTvsPtBin[ich][ict]->Write();
		}
		hAPTvsCent[ich]->Write();
	}


	// Random EP  
	for( int ich=0; ich<2; ich++ ){
		for( int ict=0; ict<NcentPolZ; ict++ ) hIMvsPz2CtRnd[ich][ict]->Write();
		for( int ipt=0; ipt<NptPol; ipt++ )    hIMvsPz2PtRnd[ich][ipt]->Write();
	}

	// EP-method
	for( int ich=0; ich<2; ich++ ){
		for( int ict=0; ict<NcentPolZ; ict++ ){
			for( int iph=0; iph<NphiRP; iph++ ){
				hIMcent    [ich][ict][iph]->Write();
				hIMvsPzCent[ich][ict][iph]->Write();
			}
		}
		for( int ict=0; ict<2; ict++ ){
			for( int ipt=0; ipt<NptPol; ipt++ ){
				for( int iph=0; iph<NphiRP; iph++ ){
					hIMpt    [ich][ict][ipt][iph]->Write();
					hIMvsPzPt[ich][ict][ipt][iph]->Write();
				}
			}
		}
	}


	// acceptance correction factor
	for( int ich=0; ich<2; ich++ ){
		hCos2cent[ich]->Write();
		hCos2eta [ich]->Write();
		for( int ict=0; ict<2; ict++ ) hCos2pt[ich][ict]->Write();
	}

	// check on Py phi dependence
	for( int ich=0; ich<2; ich++ ){
		for( int irp=0; irp<4; irp++ ){
			for( int iph=0; iph<NepPol; iph++ ){
				hIM_EP    [ich][irp][iph]->Write();
				hIMvsPy_EP[ich][irp][iph]->Write();
			}
		}
	}

	// test on self-corr.
	for( int ich=0; ich<2; ich++ ){
		for( int ict=0; ict<NcentPolZ; ict++ ){ 
			for( int ic=0; ic<2; ic++ ){
				hIMvsPz2CtSC[ich][ict][ic]->Write();
	}	}	}
	for( int ics=0; ics<2; ics++ ) hPsi2CorMod[ics]->Write();
	for( int ich=0; ich<2; ich++ ){
		for( int ipt=0; ipt<NptPol; ipt++ ){ 
			for( int ic=0; ic<2; ic++ ){
				hIMvsPz2PtSC[ich][ipt][ic]->Write();
	}	}	}

	for( int ich=0; ich<2; ich++ ){
		for( int iet=0; iet<6; iet++ ){ 
			hIMetSC[ich][iet]->Write();
			for( int ic=0; ic<3; ic++ ){
				hIMvsPz2etSC[ich][iet][ic]->Write();
			}
		}
	}

	for( int ich=0; ich<2; ich++ ){
		hIMvsVz_Pz[ich]->Write();
		hIMvsCut_Pz[ich]->Write();
	}

	// tof matching
	for( int ich=0; ich<2; ich++ ){
		for( int ict=0; ict<NcentPolZ; ict++ ){
			for( int id=0; id<4; id++ ){
				hIMvsPtv0[ich][ict][id]->Write();
	}	}	}

	// slefcor study
	for( int ich=0; ich<2; ich++ ){
		for( int k=0; k<8; k++ ){
			hPz2SCvsCent[ich][k]->Write();
	}	}
}

void CorTreeAnalyzer::SetChAsymWidthRun14(){

	// run14
	//10-20%  RMS-Ach = 0.0352416  RMS-Ak = 0.188165
	//20-30%  RMS-Ach = 0.0407078  RMS-Ak = 0.222432
	//30-40%  RMS-Ach = 0.0490833  RMS-Ak = 0.272094
	//40-50%  RMS-Ach = 0.0616321  RMS-Ak = 0.339315
	//50-60%  RMS-Ach = 0.0812182  RMS-Ak = 0.426382

	mAcWd[0] = 0.0352416; // temporary set
	mAcWd[1] = 0.0352416; // temporary set
	mAcWd[2] = 0.0352416;
	mAcWd[3] = 0.0407078;
	mAcWd[4] = 0.0490833;
	mAcWd[5] = 0.0616321;
	mAcWd[6] = 0.0812182;
	mAcWd[7] = 0.1116241; // temporary set
	mAcWd[8] = 0.1602741; // temporary set

	//mAkWd[0] = 0.188165; // temporary set
	//mAkWd[1] = 0.188165; // temporary set
	//mAkWd[2] = 0.188165;
	//mAkWd[3] = 0.222432;
	//mAkWd[4] = 0.272094;
	//mAkWd[5] = 0.339315;
	//mAkWd[6] = 0.426382;
	//mAkWd[7] = 0.509916; // temporary set
	//mAkWd[8] = 0.531882; // temporary set

}

void CorTreeAnalyzer::SysStudyN(int fhar, int fct, int fch, int fpt, float &fRP, float &fTrgEff, float &fTrkEff, float &invM, float &fcosTheStr, float &dcaA, float &dcaB, float& dcaDg, float& dcL, float &dcaV0, float &phiv0 ){

	if( fRP<-9000 ) return;

	float dphiRP = phiv0 - fRP;
	float sinN   = sin( (fhar+2.0)*dphiRP );
	float cosN   = cos( (fhar+2.0)*dphiRP );
	dphiRP = atan2( sinN, cosN ) / (fhar+2.0); 
	if( dphiRP<0 ) dphiRP = dphiRP + 2.0*TMath::Pi()/(fhar+2.0);

	float delta_phi = 2.0*TMath::Pi()/(fhar+2.0)/NphiRP;
	int fph = (int)( dphiRP / delta_phi );
	if( fph<0 || fph>=NphiRP ) return;


	if( dcaDg<mDcaDh[fct][0] && dcL>mDL[fct][0] && dcaV0<mDcaV0[fct][0] ){ 
		if( dcaA>mDcaP[fct][0]     && dcaB>mDcaPi[fct][0] )     fillSysHistoN(0,fhar,fct,fch,fpt,fph,invM,fcosTheStr,sinN,fTrgEff,fTrkEff); // 1 default
		if( dcaA>mDcaP[fct][0]     && dcaB>mDcaPi[fct][0]-0.1 ) fillSysHistoN(1,fhar,fct,fch,fpt,fph,invM,fcosTheStr,sinN,fTrgEff,fTrkEff); // 2 
		if( dcaA>mDcaP[fct][0]     && dcaB>mDcaPi[fct][0]+0.1 ) fillSysHistoN(2,fhar,fct,fch,fpt,fph,invM,fcosTheStr,sinN,fTrgEff,fTrkEff); // 3
		if( dcaA>mDcaP[fct][0]-0.1 && dcaB>mDcaPi[fct][0] )     fillSysHistoN(3,fhar,fct,fch,fpt,fph,invM,fcosTheStr,sinN,fTrgEff,fTrkEff); // 4
		if( dcaA>mDcaP[fct][0]+0.1 && dcaB>mDcaPi[fct][0] )     fillSysHistoN(4,fhar,fct,fch,fpt,fph,invM,fcosTheStr,sinN,fTrgEff,fTrkEff); // 5
	}

	if( dcaA>mDcaP[fct][0] && dcaB>mDcaPi[fct][0] ){
		if( dcaDg<mDcaDh[fct][0] && dcL>mDL[fct][0] ){ 
			if( dcaV0<mDcaV0[fct][0]-0.1 ) fillSysHistoN(5,fhar,fct,fch,fpt,fph,invM,fcosTheStr,sinN,fTrgEff,fTrkEff); // 6
			if( dcaV0<mDcaV0[fct][0]+0.1 ) fillSysHistoN(6,fhar,fct,fch,fpt,fph,invM,fcosTheStr,sinN,fTrgEff,fTrkEff); // 7
		}

		if( dcaDg<mDcaDh[fct][0] && dcaV0<mDcaV0[fct][0] ){
			if( dcL>mDL[fct][0]-1.0 ) fillSysHistoN(7,fhar,fct,fch,fpt,fph,invM,fcosTheStr,sinN,fTrgEff,fTrkEff); // 8
			if( dcL>mDL[fct][0]+1.0 ) fillSysHistoN(8,fhar,fct,fch,fpt,fph,invM,fcosTheStr,sinN,fTrgEff,fTrkEff); // 9
		}

		if( dcL>mDL[fct][0] && dcaV0<mDcaV0[fct][0] ){
			if( dcaDg<mDcaDh[fct][0]-0.1 ) fillSysHistoN( 9,fhar,fct,fch,fpt,fph,invM,fcosTheStr,sinN,fTrgEff,fTrkEff); // 10
			if( dcaDg<mDcaDh[fct][0]+0.1 ) fillSysHistoN(10,fhar,fct,fch,fpt,fph,invM,fcosTheStr,sinN,fTrgEff,fTrkEff); // 11
		}
	}

	// off-peak profile with larger BG
	if( invM>1.135 && invM<1.17 ) hBgCtLg[fhar][fch][0]->Fill( (fct+0.5)*10.0, fcosTheStr*sinN, fTrgEff*fTrkEff );
	if( invM>1.085 && invM<1.10 ) hBgCtLg[fhar][fch][1]->Fill( (fct+0.5)*10.0, fcosTheStr*sinN, fTrgEff*fTrkEff );
}

void CorTreeAnalyzer::fillSysHistoN(int sid, int fhar, int fct, int fch, int fpt, int fph, float &invM, float &cosTheStr, float &sinN, float &trgEff, float &trkEff ){

	// fct corresponds to 0=0-10%, 1=10-20%,....
	if( fct<2 ){ // 0-20%
		if( fhar==0 ) hIMptSys[fch][0][fpt]->Fill( invM, sid, trgEff*trkEff );
		hIMvsPzNptSys[fhar][fch][0][fpt]->Fill( invM, sid, cosTheStr*sinN, trgEff*trkEff );

		if( fhar==0 ) hIMPzSys[fch][0][fph]->Fill( invM, sid, trgEff*trkEff );
		hIMvsPzNSys     [fhar][fch][0][fph]->Fill( invM, sid, cosTheStr, trgEff*trkEff );
	}else if( fct>=2 && fct<=5 ){ // 20-60%
		if( fhar==0 ) hIMptSys[fch][1][fpt]->Fill( invM, sid, trgEff*trkEff );
		hIMvsPzNptSys[fhar][fch][1][fpt]->Fill( invM, sid, cosTheStr*sinN, trgEff*trkEff );

		if( fhar==0 ) hIMPzSys[fch][1][fph]->Fill( invM, sid, trgEff*trkEff );
		hIMvsPzNSys     [fhar][fch][1][fph]->Fill( invM, sid, cosTheStr, trgEff*trkEff );
	}

	if( fhar==0 ) hIMctSys[fch][fct]->Fill( invM, sid, trgEff*trkEff );
	hIMvsPzNctSys[fhar][fch][fct]->Fill( invM, sid, cosTheStr*sinN, trgEff*trkEff );
}


int CorTreeAnalyzer::LoadRunList(){
/*
	ifstream in;
	//in.open("/global/homes/t/tniida/work/AnaTree_run14/runlist_ana.txt");
	in.open("/Users/kajihara/ana_isobar");
	if( !in ){
		cout <<"Could not open runlist_wsu_all"<< endl;
		return -1;
	}

	int runid;
	while( !in.eof() ){
		in >> runid;
		if( in.eof() ) break;
		vrunid.push_back(runid);
	}

	//for( int i=0; i<vrunid.size(); i++ ){
	//	cout <<i<<" run="<< vrunid[i] << endl; 
	//}
	cout <<"LoadRunLis(): load run list, nrun="<< vrunid.size() << endl; 
	return 0;*/

	std::ifstream in;
    in.open(mRunListFile.c_str());

    if(!in.is_open()){
        cout <<"Could not open run list file: " << mRunListFile << endl;
        return -1;
    }

    int tmp;
    while(in >> tmp){
        vrunid.push_back(tmp); // AuAuの mRunid を vrunid に変更
    }
    
    cout <<"LoadRunList(): Total number of runs = "<< vrunid.size() << endl;
    return 0;
}
/*
int CorTreeAnalyzer::LoadEPcalibPara(){

	TFile *in = TFile::Open("/Users/niida/work/ana/isb_pol/jobs/epcalib_isobar_v2.root");
	if( !in ){
		cout <<"Could not open EP calib file!!"<< endl;
		return -1;
	}

	for( int ih=0; ih<2; ih++ ){
		for( int isub=0; isub<6; isub++ ){
			for( int ict=0; ict<16; ict++ ){
				for( int iz=0; iz<12; iz++ ){
					hQvPar[ih][isub][ict][iz] = (TProfile*)in->Get( Form("hQvRec_%d_%d_%d_%d",ih,isub,ict,iz) );
					if( !hQvPar[ih][isub][ict][iz] || hQvPar[ih][isub][ict][iz]->IsZombie() ){
						cout <<"LoadEPcalibPara(): Failed to get hQvPara !!"<< endl;
					}
	}	}	}	}

	return 0;
}*/

void CorTreeAnalyzer::doRecenterQv(float &fqv, int fhar, int fxy, int fbf, int fcent, int fz, int index){

	int bin  = hQvPar[fhar][fxy+fbf*2][fcent][fz]->FindBin(index);
	float mn = hQvPar[fhar][fxy+fbf*2][fcent][fz]->GetBinContent(bin);
	float sg = hQvPar[fhar][fxy+fbf*2][fcent][fz]->GetBinError  (bin);
	if( sg!=0 ){
		fqv = (fqv - mn) / sg;
	}
}

float CorTreeAnalyzer::getTrkEff(float& pt, int& fcent, int& fprt, int fid){

	float eff = 1.0;
	float ptbase = 2.4; 
	if( fprt==2 ) return eff;

	// tracking efficiency
	//if( mTrkEff==2 ){
	if( mTrkEff ){
		if( pt<ptbase ) eff = f_EffPt[fid][fprt][fcent]->Eval(pt);
		else            eff = f_EffPt[fid][fprt][fcent]->Eval(ptbase);
		eff /= mEffPtBase[fid][fprt];
		//if( eff<0.01 ) eff = 0.01; 
	}
	
	// relative efficiecny to run11
	//if( mTrkEff>0 ){
	//	if( pt<2.8 ) eff *= f_EffPtRltv[fprt][fcent]->Eval(pt);
	//	else         eff *= f_EffPtRltv[fprt][fcent]->Eval(2.8);
	//}
	if( eff<0.02 ) eff = 0.02; 

	if( eff!=0 ) eff = 1.0 / eff;
	else         eff = 1.0;
	
	return eff;
}

int CorTreeAnalyzer::SetEfficiencyRun11(){

	cout <<"Called SetEfficiencyRun11()"<< endl;

    std::ifstream in1, in2;
	std::string sf1, sf2;
	sf1 = "LamEffPt_run11AuAu200_v3.txt";
	sf2 = "LamBarEffPt_run11AuAu200_v3.txt";
	in1.open( Form("%s/%s",gSystem->Getenv("EFFDATA"),sf1.c_str()) );
	in2.open( Form("%s/%s",gSystem->Getenv("EFFDATA"),sf2.c_str()) );
	if( !in1 || !in2 ){ cout << "Could not open an efficiency file !!" << endl; return -1; }
	else {
		cout << Form("Efficiency correction -----> %s was loaded",sf1.c_str()) << endl;
		cout << Form("Efficiency correction -----> %s was loaded",sf2.c_str()) << endl << endl;
	}

	const int nct = 8;
	float tmp [3];
	float tmp2[3];
	for( int id=0; id<11; id++ ){
		for( int ict=0; ict<nct; ict++ ){

			f_EffPt[id][0][ict] = new TF1( Form("f_EffPt_%d_%d_%d",id,0,ict), "[0]*exp( -pow([1]/x,[2]) )", 0.2, 6 );
			in1 >> tmp[0] >> tmp[1] >> tmp[2];
			for( int k=0; k<3; k++ ) f_EffPt[id][0][ict]->SetParameter( k, tmp[k] );


			f_EffPt[id][1][ict] = new TF1( Form("f_EffPt_%d_%d_%d",id,1,ict), "[0]*exp( -pow([1]/x,[2]) )", 0.2, 6 );
			in2 >> tmp2[0] >> tmp2[1] >> tmp2[2];
			for( int k=0; k<3; k++ ) f_EffPt[id][1][ict]->SetParameter( k, tmp2[k] );

			if( id==0 ){
				for( int k=0; k<3; k++ ) cout << tmp[k] << "\t";
				for( int k=0; k<3; k++ ) cout << tmp2[k] << "\t";
				cout << endl;
			}
	}	}

	for( int id=0; id<11; id++ ){
		for( int ich=0; ich<2; ich++ ){
			mEffPtBase[id][ich] = f_EffPt[id][ich][nct-1]->Eval(2.4);
		}
	}

	// temporarily set 10-20%
	//for( int id=0; id<11; id++ ){
	//	for( int ich=0; ich<2; ich++ ){
	//		f_EffPt[id][ich][0] = new TF1( Form("f_EffPt_%d_%d_%d",id,ich,0), "[0]*exp( -pow([1]/x,[2]) )", 0.2, 6 );
	//		f_EffPt[id][ich][0]->SetParameters( f_EffPt[id][ich][1]->GetParameters() );
	//	}
	//}
	
	in1.close();
	in2.close();

	return 0;
}

int CorTreeAnalyzer::SetRelativeEfficiencyRun14(){

	std::ifstream in1;
	std::string sf1;
	sf1 = "LamEff_AuAu200_Run14toRun11_0-80.txt";
	in1.open( Form("%s/%s",gSystem->Getenv("EFFDATA"),sf1.c_str()) );

	if( !in1 ){ cout << "Could not open the relative efficiency file !!" << endl; return -1; }
	else {
		cout << Form("Efficiency correction -----> %s was loaded",sf1.c_str()) << endl;
	}


	const int nctE = 8;
	float tmp[3];
	float ptbase = 2.8;
	for( int ip=0; ip<2; ip++ ){
		for( int ict=0; ict<nctE; ict++ ){

			f_EffPtRltv[ip][ict] = new TF1( Form("f_EffPtRltv_%d_%d",ip,ict), "[0]*exp( -pow([1]/x,[2]) )", 0.4, 6 );
			in1 >> tmp[0] >> tmp[1] >> tmp[2];
			for( int k=0; k<3; k++ ) f_EffPtRltv[ip][ict]->SetParameter( k, tmp[k] );

			for( int k=0; k<3; k++ ) cout << tmp[k] << "\t"; 
			cout << endl;
		}
	}

	// normalize
	for( int ip=0; ip<2; ip++ ){
		float base = f_EffPtRltv[ip][nctE-1]->Eval(ptbase); // nctE would be needed to modify
		for( int ict=0; ict<nctE; ict++ ){
			f_EffPtRltv[ip][ict]->SetParameter( 0, f_EffPtRltv[ip][ict]->GetParameter(0)/base );
		}
	}

	in1.close();
	return 0;	
}

void CorTreeAnalyzer::SetNormEffRun11(){

	// [p0]*exp(-pow([p1]/x,[p2]))
	float par[3] = {0};
	par[0]=1.2195;
	par[1]=0.819046;
	par[2]=1.75485;

	//f_NormEffPtRun11 = new TF1("f_NormEffPtRun11","[0]*exp(-pow([p1]/x,[p2]))", 0.4, 2.4 );
	f_NormEffPtRun11 = new TF1("f_NormEffPtRun11","[0]*exp(-pow([1]/x,[2]))", 0.4, 2.4 );
	for( int k=0; k<3; k++ ) f_NormEffPtRun11->SetParameter(k,par[k]);

}
void CorTreeAnalyzer::SetRelativeEffIsobarToRun11(){

	// Fit function: ([p0]+[p1]*x+[p2]*pow(x,2)+[p3]*pow(x,3)+[p4]*pow(x,4))
	float par[2][8][5] = {{{0}}};

	// Ru
	if( mColl==0 ){
		par[0][0][0]=-3.56193;	par[0][0][1]=14.3226;	par[0][0][2]=-14.0442;	par[0][0][3]=5.61382;	par[0][0][4]=-0.799469;	
		par[0][1][0]=-2.75704;	par[0][1][1]=10.3279;	par[0][1][2]=-9.45335;	par[0][1][3]=3.5902;	par[0][1][4]=-0.488499;	
		par[0][2][0]=-2.38783;	par[0][2][1]=9.31973;	par[0][2][2]=-8.58978;	par[0][2][3]=3.28486;	par[0][2][4]=-0.449365;	
		par[0][3][0]=-2.29129;	par[0][3][1]=9.4743;	par[0][3][2]=-8.9036;	par[0][3][3]=3.44429;	par[0][3][4]=-0.475463;	
		par[0][4][0]=-2.12918;	par[0][4][1]=8.82221;	par[0][4][2]=-8.20371;	par[0][4][3]=3.15213;	par[0][4][4]=-0.432957;	
		par[0][5][0]=-1.91852;	par[0][5][1]=8.03708;	par[0][5][2]=-7.42528;	par[0][5][3]=2.83738;	par[0][5][4]=-0.386532;	
		par[0][6][0]=-1.83893;	par[0][6][1]=7.87568;	par[0][6][2]=-7.26816;	par[0][6][3]=2.77754;	par[0][6][4]=-0.379597;	
		par[0][7][0]=-2.04249;	par[0][7][1]=8.62924;	par[0][7][2]=-8.07936;	par[0][7][3]=3.13021;	par[0][7][4]=-0.434101;	
		par[1][0][0]=0.124225;	par[1][0][1]=5.33;	par[1][0][2]=-5.94551;	par[1][0][3]=2.43857;	par[1][0][4]=-0.342963;	
		par[1][1][0]=-0.295406;	par[1][1][1]=4.33894;	par[1][1][2]=-4.05975;	par[1][1][3]=1.47459;	par[1][1][4]=-0.184244;	
		par[1][2][0]=0.215774;	par[1][2][1]=2.78315;	par[1][2][2]=-2.53729;	par[1][2][3]=0.84826;	par[1][2][4]=-0.0904312;	
		par[1][3][0]=0.606617;	par[1][3][1]=2.22648;	par[1][3][2]=-2.23356;	par[1][3][3]=0.776409;	par[1][3][4]=-0.0847593;	
		par[1][4][0]=0.69961;	par[1][4][1]=1.6816;	par[1][4][2]=-1.56935;	par[1][4][3]=0.476864;	par[1][4][4]=-0.0386798;	
		par[1][5][0]=0.779724;	par[1][5][1]=1.19374;	par[1][5][2]=-1.04375;	par[1][5][3]=0.257194;	par[1][5][4]=-0.00557906;	
		par[1][6][0]=0.999373;	par[1][6][1]=0.622777;	par[1][6][2]=-0.466128;	par[1][6][3]=0.0103898;	par[1][6][4]=0.0315835;	
		par[1][7][0]=1.01899;	par[1][7][1]=0.552102;	par[1][7][2]=-0.336193;	par[1][7][3]=-0.0956052;	par[1][7][4]=0.0577001;	

	// Zr
	}else{
		par[0][0][0]=-3.61386;	par[0][0][1]=14.5638;	par[0][0][2]=-14.3;	par[0][0][3]=5.71934;	par[0][0][4]=-0.814902;	
		par[0][1][0]=-2.85039;	par[0][1][1]=10.7158;	par[0][1][2]=-9.85994;	par[0][1][3]=3.75827;	par[0][1][4]=-0.513204;	
		par[0][2][0]=-2.50001;	par[0][2][1]=9.80487;	par[0][2][2]=-9.10178;	par[0][2][3]=3.4955;	par[0][2][4]=-0.479984;	
		par[0][3][0]=-2.43089;	par[0][3][1]=10.1205;	par[0][3][2]=-9.6157;	par[0][3][3]=3.74942;	par[0][3][4]=-0.521523;	
		par[0][4][0]=-2.22937;	par[0][4][1]=9.35027;	par[0][4][2]=-8.78883;	par[0][4][3]=3.39746;	par[0][4][4]=-0.468788;	
		par[0][5][0]=-2.02632;	par[0][5][1]=8.55058;	par[0][5][2]=-7.99185;	par[0][5][3]=3.07962;	par[0][5][4]=-0.423002;	
		par[0][6][0]=-1.91605;	par[0][6][1]=8.27081;	par[0][6][2]=-7.71572;	par[0][6][3]=2.97318;	par[0][6][4]=-0.409677;	
		par[0][7][0]=-2.06917;	par[0][7][1]=8.87817;	par[0][7][2]=-8.35293;	par[0][7][3]=3.23312;	par[0][7][4]=-0.446688;	
		par[1][0][0]=0.150405;	par[1][0][1]=5.37909;	par[1][0][2]=-6.02822;	par[1][0][3]=2.4764;	par[1][0][4]=-0.348716;	
		par[1][1][0]=-0.303987;	par[1][1][1]=4.52907;	par[1][1][2]=-4.29311;	par[1][1][3]=1.57678;	par[1][1][4]=-0.199825;	
		par[1][2][0]=0.235543;	par[1][2][1]=2.94121;	par[1][2][2]=-2.75158;	par[1][2][3]=0.942284;	par[1][2][4]=-0.104461;	
		par[1][3][0]=0.66759;	par[1][3][1]=2.34783;	par[1][3][2]=-2.44599;	par[1][3][3]=0.876604;	par[1][3][4]=-0.100276;	
		par[1][4][0]=0.781907;	par[1][4][1]=1.76564;	par[1][4][2]=-1.75295;	par[1][4][3]=0.565837;	par[1][4][4]=-0.0524761;	
		par[1][5][0]=0.819825;	par[1][5][1]=1.36568;	par[1][5][2]=-1.32476;	par[1][5][3]=0.397243;	par[1][5][4]=-0.0288167;	
		par[1][6][0]=1.07579;	par[1][6][1]=0.59124;	par[1][6][2]=-0.486658;	par[1][6][3]=0.0212265;	par[1][6][4]=0.0306087;	
		par[1][7][0]=1.10005;	par[1][7][1]=0.649718;	par[1][7][2]=-0.549483;	par[1][7][3]=0.0155459;	par[1][7][4]=0.0389102;	
	}
	

	for( int ich=0; ich<2; ich++ ){
		for( int ict=0; ict<8; ict++ ){
			f_EffPtIsoRltv[ich][ict] = new TF1( Form("f_EffPtIsoRltv_%d_%d",ich,ict), "pol4", 0.5, 2.2 );
			for( int ip=0; ip<5; ip++ ) f_EffPtIsoRltv[ich][ict]->SetParameter(ip,par[ich][ict][ip]);
		}
	}
}

float CorTreeAnalyzer::getTrkEffIso(float& pt, int& fcent, int& fprt){

	float eff = 1.0;
	float ptbase = 2.4; 
	if( fprt>1 ) return eff;

	// run11 tracking efficiency
	if( mTrkEff ){
		if( pt<ptbase ) eff = f_NormEffPtRun11->Eval(pt);
		else            eff = 1.0;
	}
	if( eff>1.0 ) eff = 1.0;
	
	// relative efficiecny 
	if( mTrkEff>1 ){
		if( pt<2.4 ) eff *= f_EffPtIsoRltv[fprt][fcent]->Eval(pt);
		else         eff *= 1.0;
	}
	if( eff<0.02 ) eff = 0.02; 

	if( eff!=0 ) eff = 1.0 / eff;
	else         eff = 1.0;
	
	return eff;
}

