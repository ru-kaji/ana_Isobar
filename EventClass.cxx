#include "EventClass.h"
//#include "StTrackClass.h"
#include "TTree.h"

//ClassImp(EventClass);

EventClass::EventClass(){

	InitValue();
}

EventClass::EventClass( TTree *t ) : ptree(t) {

	InitValue();
}

EventClass::~EventClass(){
}

void EventClass::InitValue(){

	runid   = -1;
	centid  = -1;
	refMult = 0;
	tofMult = 0;
	vx      = -9999.0;
	vy      = -9999.0;
	vz      = -9999.0;
	zdcRate = -9999.0;
	trgEff  = 1.0;
	ach     = -9999.0;
	achp    = -9999.0;

	for( int isub=0; isub<6; isub++ ){
		tpcPsi2[isub] = -9999.0;
		tpcPsi3[isub] = -9999.0;
		tpcPsi4[isub] = -9999.0;
	}
	for( int isub=0; isub<3; isub++ ){
		bbcPsi1[isub] = -9999.0;
		bbcPsi2[isub] = -9999.0;
		zdcPsi [isub] = -9999.0;
		epdPsi1[isub] = -9999.0;
		epdPsi2[isub] = -9999.0;
		epdPsi3[isub] = -9999.0;
		epdPsi4[isub] = -9999.0;
	}
	for( int ih=0; ih<8; ih++ ){
		tpcQxB[ih] = -9999.0;
		tpcQyB[ih] = -9999.0;
		tpcQxF[ih] = -9999.0;
		tpcQyF[ih] = -9999.0;

		tpcQxB02[ih] = -9999.0;
		tpcQyB02[ih] = -9999.0;
		tpcQxF02[ih] = -9999.0;
		tpcQyF02[ih] = -9999.0;
	}
	for( int iw=0; iw<4; iw++ ){
		tpcQwB[iw] = -9999.0;
		tpcQwF[iw] = -9999.0;

		tpcQwB02[iw] = -9999.0;
		tpcQwF02[iw] = -9999.0;
	}
	for( int k=0; k<3; k++ ){
		zsmdX[k] = -9999.0;
		zsmdY[k] = -9999.0;
		if( k<2 ){
			zsmdXw[k] = -9999.0;
			zsmdYw[k] = -9999.0;
		}
	}


	ntrk    = 0;
	crNtrk  = 0;
	for( unsigned int itrk=0; itrk<ntrk_max; itrk++ ){
		ptv0     [itrk] = -9999.0;
		phiv0    [itrk] = -9999.0;
		etav0    [itrk] = -9999.0;
		phARF    [itrk] = -9999.0;
		thARF    [itrk] = -9999.0;

		dL       [itrk] = -9999.0;
		v0dca    [itrk] = -9999.0;
		dgdca    [itrk] = -9999.0;
		mass     [itrk] = -9999.0;
		dcaA     [itrk] = -9999.0;
		dcaB     [itrk] = -9999.0;
		
		tofid   [itrk] = -1;
		ptA     [itrk] = -9999;
		ptB     [itrk] = -9999;
		etaA    [itrk] = -9999;
		etaB    [itrk] = -9999;
		phiA    [itrk] = -9999;
		phiB    [itrk] = -9999;
	}
}

void EventClass::InitWrite(){

	ptree->Branch( "runid",   &runid,   "runid/I"      );
	ptree->Branch( "centid",  &centid,  "centid/B"     );
	ptree->Branch( "refMult", &refMult, "refMult/S"    );
	ptree->Branch( "tofMult", &tofMult, "tofMult/s"    );
	ptree->Branch( "vx",      &vx,      "vx/F"         );
	ptree->Branch( "vy",      &vy,      "vy/F"         );
	ptree->Branch( "vz",      &vz,      "vz/F"         );
	ptree->Branch( "zdcRate", &zdcRate, "zdcRate/F"    );
	ptree->Branch( "trgEff",  &trgEff,  "trgEff/F"     );
	ptree->Branch( "ach",     &ach,     "ach/F"        );
	ptree->Branch( "achp",    &achp,    "achp/F"       );
	ptree->Branch( "ntrk",    &ntrk,    "ntrk/i"       );

	ptree->Branch( "zdcPsi",  zdcPsi,   "zdcPsi[3]/F"  );
	ptree->Branch( "tpcPsi2", tpcPsi2,  "tpcPsi2[6]/F" );
	ptree->Branch( "tpcPsi3", tpcPsi3,  "tpcPsi3[6]/F" );
	ptree->Branch( "tpcPsi4", tpcPsi4,  "tpcPsi4[6]/F" );
	ptree->Branch( "bbcPsi1", bbcPsi1,  "bbcPsi1[3]/F" );
	ptree->Branch( "bbcPsi2", bbcPsi2,  "bbcPsi2[3]/F" );
	ptree->Branch( "epdPsi1", epdPsi1,  "epdPsi1[3]/F" );
	ptree->Branch( "epdPsi2", epdPsi2,  "epdPsi2[3]/F" );
	ptree->Branch( "epdPsi3", epdPsi3,  "epdPsi3[3]/F" );
	ptree->Branch( "epdPsi4", epdPsi4,  "epdPsi4[3]/F" );

	// |eta|>0.1
	ptree->Branch( "tpcQxB",  tpcQxB,   "tpcQxB[8]/F"  );
	ptree->Branch( "tpcQyB",  tpcQyB,   "tpcQyB[8]/F"  );
	ptree->Branch( "tpcQxF",  tpcQxF,   "tpcQxF[8]/F"  );
	ptree->Branch( "tpcQyF",  tpcQyF,   "tpcQyF[8]/F"  );
	ptree->Branch( "tpcQwB",  tpcQwB,   "tpcQwB[4]/F"  );
	ptree->Branch( "tpcQwF",  tpcQwF,   "tpcQwF[4]/F"  );
	// |eta|>0.2
	ptree->Branch( "tpcQxB02",  tpcQxB02,   "tpcQxB02[8]/F"  );
	ptree->Branch( "tpcQyB02",  tpcQyB02,   "tpcQyB02[8]/F"  );
	ptree->Branch( "tpcQxF02",  tpcQxF02,   "tpcQxF02[8]/F"  );
	ptree->Branch( "tpcQyF02",  tpcQyF02,   "tpcQyF02[8]/F"  );
	ptree->Branch( "tpcQwB02",  tpcQwB02,   "tpcQwB02[4]/F"  );
	ptree->Branch( "tpcQwF02",  tpcQwF02,   "tpcQwF02[4]/F"  );

	ptree->Branch( "zsmdX",  zsmdX,   "zsmdX[3]/F"  );
	ptree->Branch( "zsmdY",  zsmdY,   "zsmdY[3]/F"  );
	ptree->Branch( "zsmdXw", zsmdXw,  "zsmdXw[2]/F"  );
	ptree->Branch( "zsmdYw", zsmdYw,  "zsmdYw[2]/F"  );


	ptree->Branch( "ptv0",   ptv0,  "ptv0[ntrk]/F"     );
	ptree->Branch( "phiv0",  phiv0, "phiv0[ntrk]/F"    );
	ptree->Branch( "etav0",  etav0, "etav0[ntrk]/F"    );
	ptree->Branch( "phARF",  phARF, "phARF[ntrk]/F"    );
	ptree->Branch( "thARF",  thARF, "thARF[ntrk]/F"    );

	ptree->Branch( "dL",     dL,    "dL[ntrk]/F"     );
	ptree->Branch( "v0dca",  v0dca, "v0dca[ntrk]/F"  );
	ptree->Branch( "dgdca",  dgdca, "dgdca[ntrk]/F"  );
	ptree->Branch( "mass",   mass,  "mass[ntrk]/F"   );
	ptree->Branch( "dcaA",   dcaA,  "dcaA[ntrk]/F"   );
	ptree->Branch( "dcaB",   dcaB,  "dcaB[ntrk]/F"   );

	ptree->Branch( "tofid",  tofid, "tofid[ntrk]/S" );
	ptree->Branch( "ptA",    ptA,   "ptA[ntrk]/S"   );
	ptree->Branch( "ptB",    ptB,   "ptB[ntrk]/S"   );
	ptree->Branch( "etaA",   etaA,  "etaA[ntrk]/S"  );
	ptree->Branch( "etaB",   etaB,  "etaB[ntrk]/S"  );
	ptree->Branch( "phiA",   phiA,  "phiA[ntrk]/S"  );
	ptree->Branch( "phiB",   phiB,  "phiB[ntrk]/S"  );


}

void EventClass::InitRead(){

	ptree->SetBranchAddress( "runid",   &runid    );
	ptree->SetBranchAddress( "centid",  &centid   );
	ptree->SetBranchAddress( "refMult", &refMult  );
	ptree->SetBranchAddress( "tofMult", &tofMult  );
	ptree->SetBranchAddress( "vx",      &vx       );
	ptree->SetBranchAddress( "vy",      &vy       );
	ptree->SetBranchAddress( "vz",      &vz       );
	ptree->SetBranchAddress( "zdcRate", &zdcRate  );
	ptree->SetBranchAddress( "trgEff",  &trgEff   );
	ptree->SetBranchAddress( "ach",     &ach      );
	ptree->SetBranchAddress( "achp",    &achp     );
	ptree->SetBranchAddress( "ntrk",    &ntrk     );
	ptree->SetBranchAddress( "tpcPsi2", tpcPsi2   );
	ptree->SetBranchAddress( "tpcPsi3", tpcPsi3   );
	ptree->SetBranchAddress( "tpcPsi4", tpcPsi4   );
	ptree->SetBranchAddress( "bbcPsi1", bbcPsi1   );
	ptree->SetBranchAddress( "bbcPsi2", bbcPsi2   );
	ptree->SetBranchAddress( "zdcPsi",  zdcPsi    );
	ptree->SetBranchAddress( "epdPsi1", epdPsi1   );
	ptree->SetBranchAddress( "epdPsi2", epdPsi2   );
	ptree->SetBranchAddress( "epdPsi3", epdPsi3   );
	ptree->SetBranchAddress( "epdPsi4", epdPsi4   );

	// |eta|>0.1
	ptree->SetBranchAddress( "tpcQxB",  tpcQxB   );
	ptree->SetBranchAddress( "tpcQyB",  tpcQyB   );
	ptree->SetBranchAddress( "tpcQxF",  tpcQxF   );
	ptree->SetBranchAddress( "tpcQyF",  tpcQyF   );
	ptree->SetBranchAddress( "tpcQwB",  tpcQwB   );
	ptree->SetBranchAddress( "tpcQwF",  tpcQwF   );
	// |eta|>0.2
	ptree->SetBranchAddress( "tpcQxB02",  tpcQxB02   );
	ptree->SetBranchAddress( "tpcQyB02",  tpcQyB02   );
	ptree->SetBranchAddress( "tpcQxF02",  tpcQxF02   );
	ptree->SetBranchAddress( "tpcQyF02",  tpcQyF02   );
	ptree->SetBranchAddress( "tpcQwB02",  tpcQwB02   );
	ptree->SetBranchAddress( "tpcQwF02",  tpcQwF02   );

	ptree->SetBranchAddress( "zsmdX",  zsmdX   );
	ptree->SetBranchAddress( "zsmdY",  zsmdY   );
	ptree->SetBranchAddress( "zsmdXw", zsmdXw  );
	ptree->SetBranchAddress( "zsmdYw", zsmdYw  );

	ptree->SetBranchAddress( "ptv0",    ptv0      );
	ptree->SetBranchAddress( "phiv0",   phiv0     );
	ptree->SetBranchAddress( "etav0",   etav0     );
	ptree->SetBranchAddress( "phARF",   phARF     );
	ptree->SetBranchAddress( "thARF",   thARF     );
	ptree->SetBranchAddress( "dL",      dL        );
	ptree->SetBranchAddress( "v0dca",   v0dca     );
	ptree->SetBranchAddress( "dgdca",   dgdca     );
	ptree->SetBranchAddress( "mass",    mass      );
	ptree->SetBranchAddress( "dcaA",    dcaA      );
	ptree->SetBranchAddress( "dcaB",    dcaB      );
	ptree->SetBranchAddress( "tofid",   tofid     );
	ptree->SetBranchAddress( "ptA",     ptA       );
	ptree->SetBranchAddress( "ptB",     ptB       );
	ptree->SetBranchAddress( "etaA",    etaA      );
	ptree->SetBranchAddress( "etaB",    etaB      );
	ptree->SetBranchAddress( "phiA",    phiA      );
	ptree->SetBranchAddress( "phiB",    phiB      );

}

int EventClass::set_crNtrk( unsigned int itrk ){ 

	if( itrk>=ntrk_max ) return -1;
	crNtrk = itrk; 
	return 0;
}

