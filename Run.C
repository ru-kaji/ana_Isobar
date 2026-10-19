// ISOBAR
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TString.h"
#include "TStopwatch.h"
#include "EventClass.h"
#include "CorTreeAnalyzer.h"
#include <iostream>
#include <fstream>
#include <time.h>
#include <stdlib.h>
using namespace std;

//void readTree(const Char_t *inputFile="12126104.list", const Char_t *outputFile="12126104.root"){
int main(int argc, char* argv[]){

	if( argc!=4 ){ std::cout <<"Usage: ./Run inputlist outputname runid" << std::endl; return 0; }

	int fSysFlag = 1;
		// 0= OFF
		// 1= ON

  int fRemoveSharedTrk = 2;
     // 0= not remove daughter sharing
     // 1= remove shared pairs (fakes remain in single track pool)
     // 2= remove daughter sharing from single track pool

	int fTrkEff = 2;
		// 0= no trk efficiency
		// 1= include run11 eff.
		// 2= run11-eff + relative-eff
		//
		// old
		// 1= include run14/run11 eff.
		// 2= include both run14/run11 and tpc eff.(from embedding)

	int fSelfCor = 0; //キャリブレーションしないからオフにした
		// 0= OFF
		// 1= ON

	int fCalibMode = 0;
		// 0= OFF
		// 1= ON

	TStopwatch sw;
	sw.Start();

	gROOT->ProcessLine("#include <vector>");
	//gSystem->Setenv("BBCALIB","/Users/niida/calib");
	//gSystem->Setenv("EPCALIB","/Users/niida/calib");
	//gSystem->Setenv("EFFDATA","/global/homes/t/tniida/work/AnaTree_run14/eff");
	//gSystem->Setenv("EFFDATA","/Users/niida/work/ana/isb_pol"); ←これはなに？
	gSystem->Load("libStRefMultCorr.so");

	// Setup for input and output files
	//-----------------------------------------------------------------------
	TString inputFile = Form("%s",argv[1]);
	TString outFile   = Form("OutAnaTree_%s.root",argv[2]);

    //Read data list each run
	//-----------------------------------------------------------------------
	ifstream runlist;
	runlist.open( inputFile );
	if( !runlist ){
		cout << "Could not open " << inputFile << endl;
		return -1;
	}

  string test = inputFile.Data();
	int ln = test.length();
	string runnum = test.substr( ln-12, 8 ); // assuming that filename is "runnumber.root"
	cout << "input=" << test << endl;
	cout << "runnum=" << runnum << endl;

	TChain *tch = new TChain("tree");
	string tmp;
	while( runlist >> tmp ){
		cout << "Load " << tmp << endl;
		tch->Add( tmp.c_str() );
		if( runlist.eof() ) break;
	}
	if( tch==NULL || tch->IsZombie() ) cout << "Could not get tree object..." << endl;
	cout << tch->GetEntries() << endl;


	// Analyze tree data
	//-----------------------------------------------------------------------
	//TFile *in = TFile::Open( tmp.c_str() );
	//if( !in ){ cout <<"Could not open tree file!!!"<< endl; return 0; }
	//TTree *tr = (TTree*)in->Get("tree");
	//if( !tr ){ cout <<"Could not get tree..."<< endl; return 0; }
	//cout <<"tree "<< tr->GetEntries() << endl;


	EventClass *p_eve = new EventClass( tch );
	p_eve->InitRead();

	CorTreeAnalyzer *p_ana = new CorTreeAnalyzer( outFile.Data(), p_eve );
	p_ana->SetRunListFile(inputFile.Data());
	p_ana->SetRunNumber( atoi(runnum.c_str()) ); //should be before "Init()"
	p_ana->SetAnaMode(0); // no meaning for now
	p_ana->SetSysFlag(fSysFlag);
  p_ana->SetSharedTrkFlag(fRemoveSharedTrk);
	p_ana->SetTrkEfficiency(fTrkEff);
	p_ana->SetSelfCorON(fSelfCor);
	p_ana->SetCalibMode(fCalibMode);
	p_ana->SetCollision( atoi(argv[3]) );
	int flg = p_ana->Init();
	if( flg==-1 ) return -1;

	cout << "Event loop starts from here" << endl;
	Long64_t nevent = tch->GetEntries();
  p_ana->SetTotalEvents(nevent);
	for( Long64_t ieve=0; ieve<nevent; ieve++ ){

		//if( ieve>200000 ) break;

		tch->GetEntry(ieve);
		p_ana->ProcessEvent();
		p_ana->Clear();
	}
	p_ana->End();
	tch->ResetBranchAddresses();
	//in->Close();
	//-----------------------------------------------------------------------

	sw.Stop();
	cout << endl << "==================================================" << endl;
	sw.Print();
	cout << "==================================================" << endl << endl;

	return 0;
}
