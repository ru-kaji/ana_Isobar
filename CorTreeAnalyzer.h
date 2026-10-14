/*!
 * \class  CorTreeAnalyzer
 * \brief   
 *
 * -------------------------------------------------------------------------
 *
 * -------------------------------------------------------------------------
 */
// ISOBAR
#ifndef CorTreeAnalyzer_hh     
#define CorTreeAnalyzer_hh

// Include files
//#include "RtypesCore.h"
#include "Rtypes.h"
//#include "StTrackDef.h"
#include "TLorentzVector.h"
#include <string>
#include <vector>
#include <iostream>
using std::string;
using std::vector;

// Forward declarations
class EventClass;
class TFile;
class TH1F;
class TH2F;
class TH1D;
class TH2D;
class TH3D;
class TF1;
class TProfile;
class TProfile2D;
class TVector2;
class StRefMultCorr;

class TGraph;
class TCanvas;
class TLatex;

typedef unsigned int UINT;

class DhVector{
	public:
		float p_eta, p_phi, p_pt; //proton
		float pi_eta, pi_phi, pi_pt; //pion

		DhVector(){}
		DhVector(float ea, float pa, float pta,
				float eb, float pb, float ptb)
				: p_eta(ea), p_phi(pa), p_pt(pta), pi_eta(eb), pi_phi(pb), pi_pt(ptb) {}
};

class CorTreeAnalyzer {
public:

    CorTreeAnalyzer();// constructor
    CorTreeAnalyzer(const Char_t *outfile, EventClass* in_p_eve);// constructor
    ~CorTreeAnalyzer();                                                     // destructor
    
    void Clear();                    // called after every event to cleanup 
    Int_t  Init();                   // called once at the beginning of your job
    Int_t  ProcessEvent();           // invoked for every event
    Int_t  End();                    // called once at the end
	void SetRunNumber(Int_t runnum){ mRunNumber = runnum; }
	Int_t GetRunNumber() const { return mRunNumber; } 
	void SetAnaMode( int i_mode ){ ANA_MODE = i_mode; }
	void SetCollision( int i_col ){ mColl = i_col; };
	void SetSysFlag      (int val){ mSysFlag   = val; }
	void SetTrkEfficiency(int val){ mTrkEff    = val; }
	void SetSelfCorON    (int val){ mSelfCor   = val; }
	void SetCalibMode    (int val){ mCalibMode = val; }

	void SetRunListFile(std::string filename) { mRunListFile = filename; }
	void SetTotalEvents(Long64_t n) { mTotalEvents = n; } // Set the total number of events to process
	void SetSharedTrkFlag(int val){
		if( val>=0 && val<=2 ){ 
			mSharedTrkFlag = val;
			std::cout <<"CorTreeAnalyzer: SetDgSharedFlag("<<val<<")"<< std::endl;
		}else{                  
			mSharedTrkFlag = 0;
			std::cout <<"CorTreeAnalyzer: Invalid input... set as default. SetDgSharedFlag("<<mSharedTrkFlag<<")"<< std::endl;
		}
	}

private:

	EventClass *p_eve;
	StRefMultCorr *fRefMultCorr;
	Long64_t mTotalEvents;

	// data member  
    int        mEventCounter;  //!
    int        mEventCounterRef;  //!
    string     mFileName;      //!
    TFile      *mFile;         //!
	int        mRunNumber;
	int        mRunYear; // 2011, 2014, 2016
	int        mRunFlag;
	int ANA_MODE;
	int mSysFlag;
	int mTrkEff;
		// 0= no trk efficiency
		// 1= include run14/run11 eff.
		// 2= include both run14/run11 and tpc eff.(from embedding)
	int mSelfCor;
		// 0= OFF
		// 1= ON (calculate removing self-cor.)
	int mCalibMode;
		// 0= not calib mode
		// 1= Q calib mode
	int mColl;
		// 0=Ru, 1=Zr
	int mSharedTrkFlag;

	std::string mRunListFile;

	//static const int NQvCent = 8;
	static const int NQvZvtx  = 6; // 10 for run11, 6 for run14
	static const int NsubSMD  = 3;
	//static const int NsubBBC  = 6;
	//static const int NsubTPC  = 12;

	static const int NCent = 9; // for analysis
	// 0-5%, 5-10%, 10-20%, 20-30%, 30-40%, 40-50%, 50-60%, 60-70%, 70-80%

	Int_t  mAchId;
	Int_t  mAchpId;
	double mAcWd[NCent];

	//float mPionDCA;
	//float mProtonDCA;
	//float mDcaDaughters;
	//float mDCAV0;
	//float mDecayLength;

	float mDcaP [8][3];
	float mDcaPi[8][3];
	float mDcaV0[8][3];
	float mDcaDh[8][3];
	float mDL   [8][3];

	//Trks trks;
	//TrksBuff trks_buf[NCentMix][NZvtxMix][NrpMix];
	//std::vector<RPArry> rp_buf[NCentMix][NZvtxMix][NrpMix];

	std::vector<int>   trkid[7];
	std::vector<float> mass  [7];
	std::vector<float> dgApt [7];
	std::vector<float> dgAeta[7];
	std::vector<float> dgAphi[7];
	std::vector<float> dgBpt [7];
	std::vector<float> dgBeta[7];
	std::vector<float> dgBphi[7];

	std::vector<int> id_shared;
	std::vector<int> pid_shared;

	// define physical quantities
	static const int N_event = 10; 
	static const int N_Cent = 16;  
	static const int N_Vz = 12; //Isobar:12, AuAu:10
	static const int N_MIX = 3;

	//int Count[N_Cent][N_Vz];
	std::vector< std::vector<TLorentzVector> > mixingBuffer[3][N_Cent][N_Vz]; // 0: Lambda, 1: Anti-Lambda
	std::vector< std::vector<int> > mixingregionBuffer[3][N_Cent][N_Vz];
	std::vector< std::vector<DhVector> > mixingDhBuffer[3][N_Cent][N_Vz];
	std::vector< std::vector< std::pair<float, float> > > mixingSpinBuffer[3][N_Cent][N_Vz];

	TF1 *f_sgmIM[1];


	// Functions
	bool IsPairClose_pp(float eta1, float phi1, float eta2, float phi2, int ipid);
	bool IsPairClose_pipi(float eta1, float phi1, float eta2, float phi2, int ipid);
	bool IsPairClose_ppi(float eta1, float phi1, float eta2, float phi2, int ipid);
	bool IsPairClose_pip(float eta1, float phi1, float eta2, float phi2, int ipid);
	bool IsSplitCheck(float eta1, float phi1, float eta2, float phi2);

	double NormalizePhi(double phi);
	double GetAngle(double phi1, double phi2, double theta1, double theta2);
	double GetcosAngle(double phi1, double phi2, double theta1, double theta2);
    double GetPhiStar(double phi_meas,double pt, double charge);

	int GetZvtxBin      ( float zvtx    );
	bool passTPcut(float &dcaA, float &dcaB, float &dcaV0, float &dL, float &dcaAB, int fCent, int cutID, int id);
	void DefHistograms();
	void WriteHistograms();
	void SetChAsymWidthRun11();
	void SetChAsymWidthRun14();
	void SysStudyN(int fhar, int fct, int fch, int fpt, float &fRP, float &fTrgEff, float &fTrkEff, float &invM, float &fcosTheStr, float &dcaA, float &dcaB, float& dcaDg, float& dcL, float &dcaV0, float &phiv0 );
	void fillSysHistoN(int sid, int fhar, int fct, int fch, int fpt, int fph, float &invM, float &cosTheStr, float &sin2, float &trgEff, float &trkEff );

	// Pol. analysis
	//static const int NcentPol = 7;//5;
	static const int NetaPol  = 4;
	static const int NptPol   = 7;//5;
	static const int NepPol   = 5;

	static const int NphiPol  = 8;
	static const int NcentPolZ= 8;


	static const int NCent9 = 9; // for analysis
	// 0-5%, 5-10%, 10-20%, 20-30%, 30-40%, 40-50%, 50-60%, 60-70%, 70-80%
	static const int NCent8 = 8; // for analysis
	// 0-10%, 10-20%, 20-30%, 30-40%, 40-50%, 50-60%, 60-70%, 70-80%


	int LoadRunList();
	vector<int> vrunid;
	//int LoadEPcalibPara();
	void doRecenterQv(float &fqv, int fhar, int fxy, int fbf, int fcent, int fz, int index);

	// efficiency correction
	TF1 *f_EffPt [11][2][8];
	TF1 *f_EffPtRltv[2][8];
	Float_t mEffPtBase[11][2];
	float getTrkEff(float& pt, int& fcent, int& fprt, int fid);
	int SetEfficiencyRun11();
	int SetRelativeEfficiencyRun14();
	
	TF1 *f_EffPtIsoRltv [2][8];
	void SetRelativeEffIsobarToRun11();
	TF1 *f_NormEffPtRun11;
	void SetNormEffRun11();
	float getTrkEffIso(float& pt, int& fcent, int& fprt);
	

public:
	// ROOT histograms/profiles
	TH1F *hVz;
	TH1F *hRM;
	TH1F *hRMw;
	TH1F *hRMc;
	TH1F *hRMcw;
	TH1F *hCentBin;
	TH1F *hCentBinw;
	TH1F *hCentBin2;
	TH1F *hCentBinw2;

	TH1D *hcosAngle_same[3];
	TH1D *hcosAngle_mix[3];
	TH1D *hAngle_same[3];
	TH1D *hAngle_mix[3];

	TH1F *heta[2];
	TH1F *hphi[2];
	TH1F *hpt[2];
	TH1D *hinvM[2][9]; //[pid][cent]
	TH1D *hinvM_BG[2];
	TH1D *hinvM_Sig[2];

	// -------------------------------------------------------------
	// [pid][pairid][spin][cent]
	// pid    : 0=LL, 1=LbarLbar, 2=LLbar 
	// pairid : 0=SS, 1=SB&BS, 2=BB
	// spin   : 0=All, 1=[0, 90], 2=[90, 180], 3=[0, 45], 4=[45, 90], 5=[90, 135], 6=[135, 180], 7=No spin info or error
	// cent   : 0-7 for 8 centrality bins, 8 for all centrality
	// -------------------------------------------------------------
	enum { kLL=0, kAA=1, kLA=2 };        // 粒子ペア
	enum { kSS=0, kSB=1, kBB=2 };           // シグナル・背景領域
	enum { kAll=0, kPara=1, kAnti=2 };      // スピン状態

	TH1D *hmass_same_cut[5][3][4][7]; //[icut][ipid][ipair][ispin]
	TH1D *hmass_mix_cut[5][3][4][7];
    TH1D *hqinv_same_cut[5][3][4][7][9]; //[icut][ipid][ipair][ispin][icent]
    TH1D *hqinv_mix_cut[5][3][4][7][9];

	TH1D *hmass_same_nocut[3][4][7];
	TH1D *hmass_mix_nocut[3][4][7];
    TH1D *hqinv_same_nocut[3][4][7][9];
    TH1D *hqinv_mix_nocut[3][4][7][9];

	TH1D *hmass_same_low[3][4][7];
    TH1D *hmass_mix_low[3][4][7];

	TH1D *hWeightFunc;
	TH1D *hmass_mix_weigh[3][4];
	TH1D *hqinv_mix_weigh[3][4];

	// QA
	TH2D *hDetaDphiStar1[3][3];
	TH2D *hDetaDphiStar2[3][3];
	TH2D *hDetaDphiStar3[3][3];
	TH2D *hDetaDphiStar4[3][3];

	TH2D *hDetaDphi1[3];
	TH2D *hDetaDphi2[3];
	TH2D *hDetaDphi3[3];

	// QA mixing
	TH2D *hDetaDphiStar1_mix[3][3];
	TH2D *hDetaDphiStar2_mix[3][3];
	TH2D *hDetaDphiStar3_mix[3][3];
	TH2D *hDetaDphiStar4_mix[3][3];

	TH2D *hDetaDphi1_mix[3];
	TH2D *hDetaDphi2_mix[3];
	TH2D *hDetaDphi3_mix[3];

	TH1D *hDeta[3][4];
	TH1D *hDphi[3][4];

	TH1D *hDeta_mix[3][4];
	TH1D *hDphi_mix[3][4];



	//TH1F *hEtaLam[2][6]; // [ch][vz]

	// Py vs cent
	TH1D *hInvM_LamCent     [2][NsubSMD][NcentPolZ][NphiPol];
	TProfile *hInvMvsPolCent[2][NsubSMD][NcentPolZ];
	TProfile *hPyCent       [2][NsubSMD];

	// Py vs pT
	TH1D     *hIM_Pt    [2][2][NptPol];
	TProfile *hIMvsPy_Pt[2][2][NptPol];
	TProfile *hA0pt     [2][2];

	// Py vs dphi
	TH1D     *hIM_EP    [2][4][NepPol];
	TProfile *hIMvsPy_EP[2][4][NepPol];



	// Py vs Ach
	//static const int Nasm = 5;
	//TH1F *hInvMAch          [2][NsubSMD][Nasm];
	//TH1F *hInvMAchp         [2][NsubSMD][Nasm];
	//TProfile *hInvMvsPolAch [2][NsubSMD][Nasm];
	//TProfile *hInvMvsPolAchp[2][NsubSMD][Nasm];
	//TH2F *hAchCent[4];
	//TProfile *hAchBin;
	//TProfile *hAchpBin;

	// z-component
	static const int NpolRP = 7; // # of RPs
	static const int NphiRP = 6; // # of azimuthal binning
	TH1D *hInvMXYZ         [2][NpolRP][NphiRP];
	TProfile *hInvMvsPolXYZ[2][NpolRP][3][NphiRP]; // [ch][rp][xyz][dphi]

	TH1D *hInvMEta      [2][NpolRP][6][NphiRP];
	TProfile *hIMvsPzEta[2][NpolRP][6][NphiRP];

	// profile
	//TProfile *hPzCt [2][NpolRP][2]; // [ch][rp][sin, cos]
	//TProfile *hPzPt [2][NpolRP][2]; // [ch][rp][sin, cos]
	//TProfile *hPzEt [2][NpolRP][2]; // [ch][rp][sin, cos]
	TProfile *hIMvsPz2Ct[2][NpolRP][NcentPolZ][2];
	TProfile *hIMvsPz2Et[2][NpolRP][6][2];
	TProfile *hIMvsPz2PtEP[2][NptPol];

	// Profile 2D
	TProfile2D *hPzNprof2D[3][2][2];
	TProfile2D *hCos2prof2D  [2][2];
	TH3D *hIMvsPtCent[2][2];

	TH1D *hCentIso8;
	TH3D *hIMvsPtCentIso[2];


	TProfile *hPzNCt[2][2][3][2];  // [har][ch][rp][sin, cos]
	TProfile *hIMvsPzNCt[2][3][NcentPolZ]; // [har][ch][cent]

	// higher order
	TProfile *hIMvsPzN[2][2][NpolRP][2][NphiRP];
	TH1D     *hIMPzN  [2][2][NpolRP][2][NphiRP];

	// acc study
	TProfile *hInvMvsPzVz[2][NpolRP][2];    // [ch][rp][vz]
	TProfile *hPzCtVz    [2][NpolRP][2][2]; // [ch][rp][vz][cs]

	// off-mass test
	TH1F *hOffIMvsDphi[2][NcentPolZ];
	TProfile *hBgEt[2][NpolRP][2][2]; // [ch][rp][mass][sc]
	TProfile *hBgCt[2][NpolRP][2][2]; // [ch][rp][mass][sc]
	TProfile *hBgCtLg[3][2][2]; // [nhar][ch][mass]

	// Lambda v2
	TProfile *hLamV2pt[2][2];
	TProfile *hLamV2ct[2];
	TH1D *hOnIMvsDphiPt[2][5];
	TH1D *hOnIMvsDphiCt[2][NcentPolZ];
	TProfile *hCosLamPPt[2];
	TProfile *hCosLamPCt[2];


	// phi-eta acceptance
	TH2F *hphietamap[2][12];
	TH2F *hphietamapA[2][12];
	TH2F *hphietamapB[2][12];

	TH2F *hphietamapPt[2][4][3];

	TH2F *hphietamapTof[2][4][4];
	TH2F *hphietamapTofA[2][12][2];
	TH2F *hphietamapTofB[2][12][2];

	//sys
	TH2D       *hIMPzSys      [2][2][NphiRP];
	TProfile2D *hIMvsPzNSys[3][2][2][NphiRP];
	// cent
	TH2D       *hIMctSys        [2][NcentPolZ];//     [ch][cent]
	TProfile2D *hIMvsPzNctSys[3][2][NcentPolZ];//[har][ch][cent]
	// pt
	TH2D       *hIMptSys        [2][2][NptPol]; //      [ch][cent][pt]
	TProfile2D *hIMvsPzNptSys[3][2][2][NptPol]; // [har][ch][cent][pt]

	TProfile   *hAPTvsPtBin[2][2]; // [ch][cent]
	TProfile   *hAPTvsCent[2];


	// random EP test
	TProfile *hIMvsPz2CtRnd[2][NcentPolZ];
	TProfile *hIMvsPz2PtRnd[2][NptPol];
	
	// EP-method
	TH1D *hIMcent        [2][NcentPolZ][NphiRP];
	TProfile *hIMvsPzCent[2][NcentPolZ][NphiRP];
	TH1D *hIMpt        [2][2][NptPol][NphiRP];
	TProfile *hIMvsPzPt[2][2][NptPol][NphiRP];

	TProfile *hIMvsPzNCent [2][2][NcentPolZ][NphiRP]; // [har][ch][cent][phi]
	TProfile *hIMvsPzNPtDph[2][2][2][NptPol][NphiRP]; // [har][ch][cent][pt][phi]

	// acceptance factor
	TProfile *hCos2cent[2];
	TProfile *hCos2pt  [2][2];
	TProfile *hCos2eta [2];

	// test on self-corr.
	TProfile *hQvRec[2][6][16][12]; // [har][sub][cent][vz]
	TProfile *hQvPar[2][6][16][12];
	TH1F *hPsi[2][3][8][2];
	TProfile *hPsiCor[2][2]; //[har][c,s]
	TProfile *hPsi2CorMod[2];
	TProfile *hIMvsPz2CtSC[2][NcentPolZ][2];
	TProfile *hIMvsPz2PtSC[2][NptPol][2];

	TProfile2D *hIMvsVz_Pz[2];
	TProfile2D *hIMvsCut_Pz[2];

	// tof matching
	TH2F *hIMvsPtv0[2][NcentPolZ][4];

	TProfile *hZDCPsiCor[4];
	TProfile *hBBCPsiCor[2];
	TProfile *hTPCPsiCor[6];
	TProfile *hEPDPsiCor[3];
	TProfile *hEPDZDCPsi1Cor[2];

	// Pz selfcor study
	TProfile *hPz2SCvsCent[2][8];

	// eta
	TH1D *hIMetSC[2][6];
	TProfile *hIMvsPz2etSC[2][6][3];

    //ClassDef(CorTreeAnalyzer,0)
};
#endif
