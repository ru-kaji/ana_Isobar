#ifndef _STEVENTCLASS_H_
#define _STEVENTCLASS_H_

#include <vector>
#include <iostream>
#include <cmath>

class TTree;
//class StTrackClass;

class EventClass {

	public:
		EventClass( TTree *t );
		EventClass();
		~EventClass();
		void InitValue();
		void InitWrite();
		void InitRead();


		// Global variables
		//// Get functions
		unsigned int get_ntrk (){ return ntrk;     }
		char  get_centid      (){ return centid;   }
		int   get_runid       (){ return runid;    }
		short get_refMult     (){ return refMult;  }
		float get_vz          (){ return vz;       }
		float get_vx          (){ return vx;       }
		float get_vy          (){ return vy;       }
		float get_zdcRate     (){ return zdcRate;  }
		float get_trgEff      (){ return trgEff;   }
		float get_ach         (){ return ach;      }
		float get_achp        (){ return achp;     }
		unsigned short get_tofMult(){ return tofMult;  }
		
		float get_tpcPsi2 ( int i ){ return tpcPsi2[i]; } // [3] x [sub]
		float get_tpcPsi3 ( int i ){ return tpcPsi3[i]; }
		float get_tpcPsi4 ( int i ){ return tpcPsi4[i]; }
		float get_bbcPsi1 ( int i ){ return bbcPsi1[i]; } // [3]
		float get_bbcPsi2 ( int i ){ return bbcPsi2[i]; } // [3]
		float get_zdcPsi  ( int i ){ return zdcPsi [i]; } // [3]
		float get_epdPsi1 ( int i ){ return epdPsi1[i]; } // [east, west, e+w]
		float get_epdPsi2 ( int i ){ return epdPsi2[i]; } // [east, west, e+w]
		float get_epdPsi3 ( int i ){ return epdPsi3[i]; } // [east, west, e+w]
		float get_epdPsi4 ( int i ){ return epdPsi4[i]; } // [east, west, e+w]

		// |eta|>0.1
		float get_tpcQxB  ( int i ){ return tpcQxB [i]; }
		float get_tpcQyB  ( int i ){ return tpcQyB [i]; }
		float get_tpcQxF  ( int i ){ return tpcQxF [i]; }
		float get_tpcQyF  ( int i ){ return tpcQyF [i]; }
		float get_tpcQwB  ( int i ){ return tpcQwB [i]; }
		float get_tpcQwF  ( int i ){ return tpcQwF [i]; }
		// |eta|>0.2
		float get_tpcQxB02( int i ){ return tpcQxB02[i]; }
		float get_tpcQyB02( int i ){ return tpcQyB02[i]; }
		float get_tpcQxF02( int i ){ return tpcQxF02[i]; }
		float get_tpcQyF02( int i ){ return tpcQyF02[i]; }
		float get_tpcQwB02( int i ){ return tpcQwB02[i]; }
		float get_tpcQwF02( int i ){ return tpcQwF02[i]; }

		float get_zsmdX ( int i ){ return zsmdX [i]; }
		float get_zsmdY ( int i ){ return zsmdY [i]; }
		float get_zsmdXw( int i ){ return zsmdXw[i]; }
		float get_zsmdYw( int i ){ return zsmdYw[i]; }
		
		//// Set functions
		void set_ntrk     ( unsigned int val ){ ntrk    = val; }
		void set_centid   ( char  val )       { centid  = val; }
		void set_runid    ( int   val )       { runid   = val; }
		void set_refMult  ( short val )       { refMult = val; }
		void set_vx       ( float val )       { vx      = val; }
		void set_vy       ( float val )       { vy      = val; }
		void set_vz       ( float val )       { vz      = val; }
		void set_zdcRate  ( float val )       { zdcRate = val; }
		void set_trgEff   ( float val )       { trgEff  = val; }
		void set_ach      ( float val )       { ach     = val; }
		void set_achp     ( float val )       { achp    = val; }
		void set_tofMult  ( unsigned short val ){ tofMult = val; }

		void set_tpcPsi2 ( int i, float val ){ tpcPsi2[i] = val; } // [3] x [sub]
		void set_tpcPsi3 ( int i, float val ){ tpcPsi3[i] = val; }
		void set_tpcPsi4 ( int i, float val ){ tpcPsi4[i] = val; }
		void set_bbcPsi1 ( int i, float val ){ bbcPsi1[i] = val; } // [3]
		void set_bbcPsi2 ( int i, float val ){ bbcPsi2[i] = val; } // [3]
		void set_zdcPsi  ( int i, float val ){ zdcPsi [i] = val; } // [3]
		void set_epdPsi1 ( int i, float val ){ epdPsi1[i] = val; } 
		void set_epdPsi2 ( int i, float val ){ epdPsi2[i] = val; } 
		void set_epdPsi3 ( int i, float val ){ epdPsi3[i] = val; } 
		void set_epdPsi4 ( int i, float val ){ epdPsi4[i] = val; } 

		// |eta|>0.1
		void set_tpcQxB  ( int i, float val ){ tpcQxB [i] = val; }
		void set_tpcQyB  ( int i, float val ){ tpcQyB [i] = val; }
		void set_tpcQxF  ( int i, float val ){ tpcQxF [i] = val; }
		void set_tpcQyF  ( int i, float val ){ tpcQyF [i] = val; }
		void set_tpcQwB  ( int i, float val ){ tpcQwB [i] = val; }
		void set_tpcQwF  ( int i, float val ){ tpcQwF [i] = val; }
		// |eta|>0.2
		void set_tpcQxB02( int i, float val ){ tpcQxB02[i] = val; }
		void set_tpcQyB02( int i, float val ){ tpcQyB02[i] = val; }
		void set_tpcQxF02( int i, float val ){ tpcQxF02[i] = val; }
		void set_tpcQyF02( int i, float val ){ tpcQyF02[i] = val; }
		void set_tpcQwB02( int i, float val ){ tpcQwB02[i] = val; }
		void set_tpcQwF02( int i, float val ){ tpcQwF02[i] = val; }

		void set_zsmdX ( int i, float val ){ zsmdX [i] = val; }
		void set_zsmdY ( int i, float val ){ zsmdY [i] = val; }
		void set_zsmdXw( int i, float val ){ zsmdXw[i] = val; }
		void set_zsmdYw( int i, float val ){ zsmdYw[i] = val; }

		// Single track variables
		int set_crNtrk( unsigned int itrk );

		void set_ptv0     ( float val ){ ptv0    [crNtrk] = val; }
		void set_etav0    ( float val ){ etav0   [crNtrk] = val; }
		void set_phiv0    ( float val ){ phiv0   [crNtrk] = val; }
		void set_phARF    ( float val ){ phARF   [crNtrk] = val; }
		void set_thARF    ( float val ){ thARF   [crNtrk] = val; }

		void set_dL       ( float val ){ dL     [crNtrk] = val; }
		void set_v0dca    ( float val ){ v0dca  [crNtrk] = val; }
		void set_dgdca    ( float val ){ dgdca  [crNtrk] = val; }
		void set_mass     ( float val ){ mass   [crNtrk] = val; }
		void set_dcaA     ( float val ){ dcaA   [crNtrk] = val; }
		void set_dcaB     ( float val ){ dcaB   [crNtrk] = val; }

		void set_tofid( short val ){ tofid[crNtrk] = val; }
		void set_ptA  ( short val ){ ptA  [crNtrk] = val; }
		void set_ptB  ( short val ){ ptB  [crNtrk] = val; }
		void set_etaA ( short val ){ etaA [crNtrk] = val; }
		void set_etaB ( short val ){ etaB [crNtrk] = val; }
		void set_phiA ( short val ){ phiA [crNtrk] = val; }
		void set_phiB ( short val ){ phiB [crNtrk] = val; }


		short get_charge   ( int i ){ return ( (ptv0[i]>0) ? +1 : -1); }
		float get_ptv0     ( int i ){ return fabs(ptv0[i]); }
		float get_phiv0    ( int i ){ return phiv0 [i]; }
		float get_etav0    ( int i ){ return etav0 [i]; }
		float get_phARF    ( int i ){ return phARF [i]; }
		float get_thARF    ( int i ){ return thARF [i]; }

		float get_dL       ( int i ){ return dL    [i]; }
		float get_v0dca    ( int i ){ return v0dca [i]; }
		float get_dgdca    ( int i ){ return dgdca [i]; }
		float get_mass     ( int i ){ return mass  [i]; }
		float get_dcaA     ( int i ){ return dcaA  [i]; }
		float get_dcaB     ( int i ){ return dcaB  [i]; }

		short get_tofid( int i ){ return tofid[i];      }
		float get_ptA  ( int i ){ return ptA[i]/1000.;  }
		float get_ptB  ( int i ){ return ptB[i]/1000.;  }
		float get_etaA ( int i ){ return etaA[i]/1000.; }
		float get_etaB ( int i ){ return etaB[i]/1000.; }
		float get_phiA ( int i ){ return phiA[i]/1000.; }
		float get_phiB ( int i ){ return phiB[i]/1000.; }

	private:
		TTree *ptree;
		int   runid;      // run number
		char  centid;     // 0=75-80%,..., 15=0-5%
		short refMult;    // refMult for run11, grefMult for run14
		unsigned short tofMult;    // tofMult 
		float vz;         // TPC-zvertex satsifying |vz-vzvpd|<3
		float vx;
		float vy;
		float zdcRate;    // zdc coincidence rate
		float trgEff;     // trigger weight
		float ach;        // charge asymmetry
		float achp;       // charge asymmetry excluidng lambda's daughters
		float tpcPsi2[6]; // TPC-Psi2,  3x2 (East, West, East+West) x (|eta|>0.1, |eta|>0.5) 
		float tpcPsi3[6]; // TPC-Psi3 
		float tpcPsi4[6]; // TPC-Psi4 
		float bbcPsi1[3]; // BBC-Psi1: East, West, East+West
		float bbcPsi2[3]; // BBC-Psi2: East, West, East+West
		float zdcPsi [3]; // ZDCSMD-Psi1: East, West, East+West
		float epdPsi1[3]; // EPD-Psi2
		float epdPsi2[3]; // EPD-Psi2
		float epdPsi3[3]; // EPD-Psi2
		float epdPsi4[3]; // EPD-Psi2

		// tpc "raw" flow vectors for -1<eta<-0.1 and 0.1<eta<1 and h+/h-
		float tpcQxB[8]; // backward eta, harmonics = 1-4
		float tpcQyB[8]; // 0-3: 4 harmonics for h+, 4-7: 4 harmonics for h-
		float tpcQxF[8]; // forward eta, harmonics = 1, 2, 3, 4 for h+/h-
		float tpcQyF[8]; 
		float tpcQwB[4]; // 0= multiplicity for h+,  1= sum of pT for h+
		float tpcQwF[4]; // 2= multiplicity for h-,  3= sum of pT for h-
						 // "B"ackward and "F"orward

		// tpc "raw" flow vectors for -1<eta<-0.2 and 0.2<eta<1 and h+/h-
		float tpcQxB02[8]; // backward eta, harmonics = 1-4
		float tpcQyB02[8]; // 0-3: 4 harmonics for h+, 4-7: 4 harmonics for h-
		float tpcQxF02[8]; // forward eta, harmonics = 1, 2, 3, 4 for h+/h-
		float tpcQyF02[8]; 
		float tpcQwB02[4]; // 0= multiplicity for h+,  1= sum of pT for h+
		float tpcQwF02[4]; // 2= multiplicity for h-,  3= sum of pT for h-
						 // "B"ackward and "F"orward

		// zdcsmd mean X and Y after recentering
		float zsmdX [3];
		float zsmdY [3];
		float zsmdXw[2];
		float zsmdYw[2];

		unsigned int ntrk;// number of v0
		unsigned int crNtrk;
		static const unsigned int ntrk_max = 1000; 

		// v0 info
		float ptv0 [ntrk_max]; // pt * charge
		float phiv0[ntrk_max];
		float etav0[ntrk_max];
		float phARF[ntrk_max]; // phi   of daughter p(pbar) in v0 rest frame 
		float thARF[ntrk_max]; // theta of daughter p(pbar) in v0 rest frame

		// for topological cuts
		float dL   [ntrk_max]; // decay length
		float v0dca[ntrk_max]; // v0 dca
		float dgdca[ntrk_max]; // dca between p and pi
		float mass [ntrk_max]; // invarinat mass
		float dcaA [ntrk_max]; // dca of p
		float dcaB [ntrk_max]; // dca of pi

		// daughter info
		short tofid[ntrk_max]; // tof matched info.: 0=both matched, 1= both unmatched, 2=only pi matched, 3=only p matched
		short ptA  [ntrk_max]; // pt*1000 of p with short
		short ptB  [ntrk_max]; // pt*1000 of pion with short
		short etaA [ntrk_max]; 
		short etaB [ntrk_max];
		short phiA [ntrk_max];
		short phiB [ntrk_max];

		//ClassDef(EventClass,1);
};

#endif
