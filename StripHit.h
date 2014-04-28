// Andriy Zatserklyaniy, April 17, 2014

#include <TROOT.h>
#include <TMath.h>

#include <iostream>

using std::cout;     using std::endl;

struct StripHit {
   //
   //          V-board chips
   //    6                       5
   //    11                      0
   //
   // V: 0..383  FPGA#3    383..767
   // T: 767---#11--0 767---#10---0
   // V: 0..383  FPGA#2   384..767 
   // T: 767---#9--0  767---#8----0
   //
   // T: 0---#6---767 0----#7----767
   // V: 767..384  FPGA#1    383..0
   // T: 0----#4--767 0----#5----767
   // V: 767..384  FPGA#0    383..0
   //
   //                ^
   //                |
   //                |  Particle Eye View
   //
   // event 125 good event
   // 125	 iFPGA = 0 trackerChip->address =  9 nstrips[0] = 1 nfirst[0] = 9
   // 125	 iFPGA = 1 trackerChip->address =  9 nstrips[0] = 1 nfirst[0] = 14
   // 125	 iFPGA = 2 trackerChip->address =  2 nstrips[0] = 1 nfirst[0] = 6
   // 125	 iFPGA = 3 trackerChip->address =  3 nstrips[0] = 1 nfirst[0] = 62
   // 125	 iFPGA = 4 trackerChip->address =  8 nstrips[0] = 1 nfirst[0] = 35
   // 125	 iFPGA = 6 trackerChip->address =  8 nstrips[0] = 1 nfirst[0] = 45
   // 125	 iFPGA = 9 trackerChip->address =  2 nstrips[0] = 1 nfirst[0] = 13
   // 125	 iFPGA = 11 trackerChip->address =  3 nstrips[0] = 1 nfirst[0] = 63
   //
   Int_t vstrip[4][2*384];    // particle eye view: layer 0: left part 384..768, right part: 0..383
   Int_t tstrip[4][4*384];
   Int_t* stripFPGA[12];
   StripHit() {
      // link v layer
      stripFPGA[0] = &vstrip[0][0];
      stripFPGA[1] = &vstrip[1][0];
      stripFPGA[2] = &vstrip[2][0];
      stripFPGA[3] = &vstrip[3][0];
      // link t layer
      stripFPGA[4] = &tstrip[0][0];
      stripFPGA[5] = &tstrip[0][2*384];
      stripFPGA[6] = &tstrip[1][0];
      stripFPGA[7] = &tstrip[1][2*384];
      stripFPGA[8] = &tstrip[2][0];
      stripFPGA[9] = &tstrip[2][2*384];
      stripFPGA[10] = &tstrip[3][0];
      stripFPGA[11] = &tstrip[3][2*384];

      clear();
   }
   void clear() {
      for (int iFPGA=0; iFPGA<12; ++iFPGA) for (int istrip=0; istrip<2*384; ++istrip) stripFPGA[iFPGA][istrip] = 0;
   }
   void Analyze() {
      cout<< "v-layer" <<endl;
      for (int ilayer=0; ilayer<4; ++ilayer) {
         for (int istrip=0; istrip<2*384; ++istrip) {
            if (vstrip[ilayer][istrip] > 0) cout<< "vstrip[" << ilayer << "][" << istrip << "] = " << vstrip[ilayer][istrip] <<endl;
         }
      }
      cout<< "t-layer" <<endl;
      for (int ilayer=0; ilayer<4; ++ilayer) {
         for (int istrip=0; istrip<4*384; ++istrip) {
            if (tstrip[ilayer][istrip] > 0) cout<< "tstrip[" << ilayer << "][" << istrip << "] = " << tstrip[ilayer][istrip] <<endl;
         }
      }
   }
};

class GeoHit: public StripHit {
public:
   Double_t vhit[4][2*384];
   Double_t thit[4][4*384];
   GeoHit(): StripHit() {
      clear();
   }
   void clear() {
      StripHit::clear();
      for (int ilayer=0; ilayer<4; ++ilayer) {
         for (int istrip=0; istrip<2*384; ++istrip) {
            vhit[ilayer][istrip] = 0;
            thit[ilayer][istrip] = 0;
         }
         for (int istrip=2*384; istrip<4*384; ++istrip) thit[ilayer][istrip] = 0;
      }
   }
};

struct EventOutput {
   Bool_t good;
   Float_t t[4];
   Float_t v[4];
   Float_t u[4];
   Float_t wepl;
   EventOutput() {clear();}
   void clear() {
      good = kFALSE;
      for (int ilayer=0; ilayer<4; ++ilayer) {
         t[ilayer] = 0;
         v[ilayer] = 0;
         u[ilayer] = 0;
      }
      wepl = 0;
   }
};
