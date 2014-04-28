// Andriy Zatserklyaniy, April 17, 2014

#include "DataFormat.h"
#include "StripHit.h"
//#include "stripPosition.h"

#include <TROOT.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TObject.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2.h>
#include <TF1.h>
#include <TMath.h>
#include <TGraph.h>
#include <TLine.h>
#include <TMarker.h>
#include <TTimer.h>

#include <iostream>
#include <iomanip>
#include <cassert>
#include <sstream>
#include <vector>
#include <fstream>
#include <ctime>
#include <cstdio>

#include <vector>
#include <map>

using std::cout;     using std::endl;

// enum SensorVTType {Unknown, V, T};
// 
// std::ostream& operator <<(std::ostream& os, const SensorVTType& type)
// {
//    os<< (type == Sensor::V? "V": "T");
//    return os;
// }

class Hit1D {
public:
   Double_t x_;
};

class Sensor {
   // sensor with 384 strips
public:
   enum VTType {Unknown, V, T};
   VTType type_;                          // info
   Int_t layer_;                          // info
   Int_t sensor_;                         // info
   Double_t pitch_;
   Double_t strip1_;                      // coordinate of the first strip
   Double_t dir_;                         // +1 or -1 direction of the strip numbers wrt axis
   std::map<Int_t,Int_t> clusters_;       // map<NFirst,NStrips> strip 0..383
   std::vector<Double_t> hits_;           // center of gravity
   std::vector<Int_t> dead_;              // dead channels
   std::vector<Int_t> noisy_;             // noisy channels
public:
   Double_t Coordinate(Double_t strip) const {return strip1_ + dir_*strip*pitch_;}
   Sensor(): type_(Unknown), pitch_(0.228) {}
   Sensor(VTType type): type_(type), pitch_(0.228) {}
   void Set_strip1(Double_t strip1, Double_t dir) {
      strip1_ = strip1;
      dir_ = dir;
   }
   void Set_strip1_dir_info(Double_t strip1, Double_t dir, Int_t layer, Int_t sensor) {
      strip1_ = strip1;
      dir_ = dir;
      layer_ = layer;
      sensor_ = sensor;
   }
   void AddCluster(Int_t chip, Int_t nfirst, Int_t nstrips) {
      Int_t istrip = (chip % 6)*64 + (63 - nfirst);
      clusters_[istrip-nstrips+1] = nstrips;
   }
   void GetHits() {
      if (clusters_.size() == 0) return;
      // cout<< "Sensor::GetHits: clusters_.size() = " << clusters_.size() <<endl;
      std::map<Int_t,Int_t>::const_iterator it = clusters_.begin();
      while (it != clusters_.end()) {
         Int_t nfirst = it->first;
         Int_t nstrips = it->second;
         ++it;
         // if (it != clusters_.end()) cout<< "it->first = " << it->first << " it->second = " << it->second <<endl;
         if (it != clusters_.end() && it->first == nfirst+nstrips) {
            // cout<< "add strips for the neibor cluster" <<endl;
            nstrips += it->second;
            ++it;                      // goto the next cluster
         }
         Double_t hit = Coordinate(nfirst + 0.5*(nstrips-1));
         hits_.push_back(hit);
         cout<< "type: " << (type_ == V? "V": "T") << " layer " << layer_ << " sensor " << sensor_ << " nfirst = " << std::setw(3) << nfirst << " nstrips = " << nstrips << " coordinate = " << Coordinate(nfirst) << " hit = " << hit <<endl;
      }
   }
};

class TSensor: public Sensor {
public:
   TSensor(): Sensor(T) {}
};

class VSensor: public Sensor {
public:
   VSensor(): Sensor(V) {}
};

class PCTSensors {
   //
   // sensor 0: 0..383
   // sensor 1: 384..763
   //
   // V3:      v[3][0]           v[3][1]
   // T3:   t[3][3]  t[3][2]  t[3][1]  t[3][0]
   //
   // V2:      v[2][0]           v[2][1]
   // T2:   t[2][3]  t[2][2]  t[2][1]  t[2][0]
   //
   // T1:   t[1][0]  t[1][1]  t[1][2]  t[1][3]
   // V1:      v[1][1]           v[1][0]
   //
   // T0:   t[0][0]  t[0][1]  t[0][2]  t[0][3]
   // V0:      v[0][1]           v[0][0]
   //
   // t-axis direction:  <--------------------
   //
   // v-axis:  chip6          chip5    |
   //          |              |        |
   //          chip11         chip0    V v-axis direction
   //
public:
   VSensor vSensor[4][2];     // [layer][sensor]
   TSensor tSensor[4][4];     // [layer][sensor]
   void GetHits() {
      for (int ilayer=0; ilayer<4; ++ilayer) for (int isensor=0; isensor<2; ++isensor) vSensor[ilayer][isensor].GetHits();
      for (int ilayer=0; ilayer<4; ++ilayer) for (int isensor=0; isensor<4; ++isensor) tSensor[ilayer][isensor].GetHits();
   }
   PCTSensors() {
      const Double_t pitch = 0.228;
      Double_t vPin[4] = {0.055,0.055,0.055,0.055};
      int vBoard[4] = {6,4,2,3};  	// fpga to V-board translation; this will change if spares are swapped
      Double_t firstStripV[7][2] = {		// Distance from the alignment pin to the first strip
         {-43.7193, -43.716},       // Board V0 doesn't exist
         {-43.7193, -43.716},
         {-43.7193, -43.716},
         {-43.7193, -43.716},
         {-43.7193, -43.716},	      // These numbers are actually from V4
         {-43.7193, -43.716},
         {-43.5855, -43.5865}       // these are actually from T6
      };

      //
      // the firstStripV above is actually the strip #0 in the chip #5 (bottom part of the sensor but top part of the setup)
      // Example: the layer 0 uses vBoard[0], which is #6, correspondent distance to the pin is -43.5855,
      // so the position of the top right strip (chip #5, strip #0) is -43.5855 + 0.055 = -43.530 mm
      //
      // I will use as the first strip the strip #63 in the chip #0 (bottom part of the setup)
      // A coordinate of the my first strip then is -43.5855 + 0.055 + 383*0.228 = +43.793 mm
      //
      // The first strip of the high chip address part (chip #6, strip #63) is at the same position like for the Robert's data:
      // 0.055 - 43.5865 = -43.5305 mm
      //

      Int_t layer;                  // 0..3
      Int_t board;                  // board production number: V-boards: 0..6, T-boards 0..7
      Int_t dir;                    // chip direction wrt corresponding axis
      Int_t sensor;                 // sensor number: V-board: 0..1, T-board: 0..3

      Double_t firstStripVabs[7][2];      // absolute coordinate of the first strip
      sensor = 0;
      for (int iboard=0; iboard<7; ++iboard) firstStripVabs[iboard][sensor] = vPin[0] + firstStripV[iboard][sensor] + 383*pitch;  // chip 0, strip 63
      sensor = 1;
      for (int iboard=0; iboard<7; ++iboard) firstStripVabs[iboard][sensor] = vPin[0] + firstStripV[iboard][sensor];  // chip 6, strip 63
      // for (int iboard=0; iboard<7; ++iboard) {
      //    for (int iside=0; iside<2; ++iside) {
      //       cout<< firstStripVabs[iboard][iside] << "\t ";
      //    }
      //    cout<<endl;
      // }

      layer = 0;
      board = vBoard[layer];  dir = -1; sensor = 0;
      vSensor[layer][sensor].Set_strip1_dir_info(firstStripVabs[board][sensor], dir, layer, sensor);
      board = vBoard[layer];  dir = +1; sensor = 1;
      vSensor[layer][sensor].Set_strip1_dir_info(firstStripVabs[board][sensor], dir, layer, sensor);
      //
      layer = 1;  // identical to the layer 0
      board = vBoard[layer];  dir = -1; sensor = 0;
      vSensor[layer][sensor].Set_strip1_dir_info(firstStripVabs[board][sensor], dir, layer, sensor);
      board = vBoard[layer];  dir = +1; sensor = 1;
      vSensor[layer][sensor].Set_strip1_dir_info(firstStripVabs[board][sensor], dir, layer, sensor);
      //
      layer = 2;
      board = vBoard[layer];  dir = -1; sensor = 0;
      vSensor[layer][sensor].Set_strip1_dir_info(firstStripVabs[board][sensor], dir, layer, sensor);
      board = vBoard[layer];  dir = +1; sensor = 1;
      vSensor[layer][sensor].Set_strip1_dir_info(firstStripVabs[board][sensor], dir, layer, sensor);
      //
      layer = 3;  // identical to the layer 2
      board = vBoard[layer];  dir = -1; sensor = 0;      // like sensor 1 in the layer 0
      vSensor[layer][sensor].Set_strip1_dir_info(firstStripVabs[board][sensor], dir, layer, sensor);
      board = vBoard[layer];  dir = +1; sensor = 1;      // like sensor 0 in the layer 0
      vSensor[layer][sensor].Set_strip1_dir_info(firstStripVabs[board][sensor], dir, layer, sensor);

      Double_t tPin[4] = {215.24, 211.25, -211.25, -215.24};   // T coordinate of alignment pin per layer
      //Double_t tDir[4] = {-1.0, -1.0, 1.0, 1.0};               // T direction of increasing strip number per layer

      Int_t tBoard[4] = {5, 4, 1, 3};     // T-layer to physical T board translation. This will change if spares are swapped in!
      Double_t firstStripT[7][4] = {      // First strip location rel to pin for each sensor on each physical board
         {-999., -999., -999., -999.},		// Board 0 doesn't exist.  We manufactured 6 boards.
         {38.58, 126.85, 215.11, 303.37},
         {38.58, 126.85, 215.11, 303.37},
         {38.58, 126.85, 215.11, 303.37},
         {38.58, 126.85, 215.11, 303.37}, // These numbers actually come from board T4
         {38.62, 126.90, 215.16, 303.41}, // these numbers actually come from board T5
         {38.58, 126.85, 215.11, 303.37} 
      };

      //Double_t firstStripTabs[7][4];      // absolute coordinate of the first strip

      //for (int ilayer=0; ilayer<2; ++ilayer) {
      //   // layers 0,1: direction of the chip address is opposite to the t-axis
      //   for (int iboard=1; iboard<7; ++iboard) {     // ignore board #0
      //      for (int isensor=0; isensor<4; ++isensor) {
      //         firstStripTabs[iboard][isensor] = tPin[ilayer]firstStripT[iboard][isensor];
      //      }
      //   }
      //}

      // tSensor[0][0].Set_strip1_dir_info(215.24-38.58,  -1, 0, 0);
      // tSensor[0][1].Set_strip1_dir_info(215.24-126.85, -1, 0, 1);
      // tSensor[0][2].Set_strip1_dir_info(215.24-215.11, -1, 0, 2);
      // tSensor[0][3].Set_strip1_dir_info(215.24-303.37, -1, 0, 3);

      // tSensor[1][0].Set_strip1_dir_info(211.25-38.58,  -1, 1, 0);
      // tSensor[1][1].Set_strip1_dir_info(211.25-126.85, -1, 1, 1);
      // tSensor[1][2].Set_strip1_dir_info(211.25-215.11, -1, 1, 2);
      // tSensor[1][3].Set_strip1_dir_info(211.25-303.37, -1, 1, 3);

      dir = -1;
      for (int ilayer=0; ilayer<2; ++ilayer) {
         Double_t pin = tPin[ilayer];
         board = tBoard[ilayer];                         // production number of the board for this layer
         for (int isensor=0; isensor<4; ++isensor) {
            tSensor[ilayer][isensor].Set_strip1_dir_info(pin + dir*firstStripT[board][isensor], dir, ilayer, isensor);
         }
      }

      dir = +1;
      for (int ilayer=2; ilayer<4; ++ilayer) {
         Double_t pin = tPin[ilayer];
         board = tBoard[ilayer];                         // production number of the board for this layer
         for (int isensor=0; isensor<4; ++isensor) {
            tSensor[ilayer][isensor].Set_strip1_dir_info(pin + dir*firstStripT[board][isensor], dir, ilayer, isensor);
         }
      }

      // tSensor[2][0].Set_strip1_dir_info(-211.25+38.58,  1, 2, 0);
      // tSensor[2][1].Set_strip1_dir_info(-211.25+126.85, 1, 2, 1);
      // tSensor[2][2].Set_strip1_dir_info(-211.25+215.11, 1, 2, 2);
      // tSensor[2][3].Set_strip1_dir_info(-211.25+303.37, 1, 2, 3);

      // tSensor[3][0].Set_strip1_dir_info(-215.24+38.58,  1, 3, 0);
      // tSensor[3][1].Set_strip1_dir_info(-215.24+126.85, 1, 3, 1);
      // tSensor[3][2].Set_strip1_dir_info(-215.24+215.11, 1, 3, 2);
      // tSensor[3][3].Set_strip1_dir_info(-215.24+303.37, 1, 3, 3);
   }
   void AddCluster(Int_t FPGA, Int_t chip, Int_t nfirst, Int_t nstrips) {
      int sensor;
      int ilayer;
      switch (FPGA) {
         // v-layers
         case 0:
            ilayer = 0;
            sensor = chip < 6? 0: 1;
            vSensor[ilayer][sensor].AddCluster(chip, nfirst, nstrips);
            //cout<< "added chip " << chip << " nfirst " << nfirst << " nstrips " << nstrips << " to vSensor[" << ilayer<< "][" << sensor << "]" <<endl;
            break;
         case 1:
            ilayer = 1;
            sensor = chip < 6? 0: 1;
            vSensor[ilayer][sensor].AddCluster(chip, nfirst, nstrips);
            //cout<< "added chip " << chip << " nfirst " << nfirst << " nstrips " << nstrips << " to vSensor[" << ilayer<< "][" << sensor << "]" <<endl;
            break;
         case 2:
            ilayer = 2;
            sensor = chip < 6? 0: 1;
            vSensor[ilayer][sensor].AddCluster(chip, nfirst, nstrips);
            //cout<< "added chip " << chip << " nfirst " << nfirst << " nstrips " << nstrips << " to vSensor[" << ilayer<< "][" << sensor << "]" <<endl;
            break;
         case 3:
            ilayer = 3;
            sensor = chip < 6? 0: 1;
            vSensor[ilayer][sensor].AddCluster(chip, nfirst, nstrips);
            //cout<< "added chip " << chip << " nfirst " << nfirst << " nstrips " << nstrips << " to vSensor[" << ilayer<< "][" << sensor << "]" <<endl;
            break;

         // t-layers
         case 4:     // elements 0..1
            ilayer = 0;
            sensor = chip < 6? 0: 1;
            tSensor[ilayer][sensor].AddCluster(chip, nfirst, nstrips);
            //cout<< "added chip " << chip << " nfirst " << nfirst << " nstrips " << nstrips << " to tSensor[" << ilayer<< "][" << sensor << "]" <<endl;
            break;
         case 5:     // elements 2..3
            ilayer = 0;
            sensor = chip < 6? 0: 1;
            tSensor[ilayer][2+sensor].AddCluster(chip, nfirst, nstrips);
            //cout<< "added chip " << chip << " nfirst " << nfirst << " nstrips " << nstrips << " to tSensor[" << ilayer<< "][" << sensor << "]" <<endl;
            break;
         case 6:     // elements 0..1
            ilayer = 1;
            sensor = chip < 6? 0: 1;
            tSensor[ilayer][sensor].AddCluster(chip, nfirst, nstrips);
            //cout<< "added chip " << chip << " nfirst " << nfirst << " nstrips " << nstrips << " to tSensor[" << ilayer<< "][" << sensor << "]" <<endl;
            break;
         case 7:     // elements 2..3
            ilayer = 1;
            sensor = chip < 6? 0: 1;
            tSensor[ilayer][2+sensor].AddCluster(chip, nfirst, nstrips);
            //cout<< "added chip " << chip << " nfirst " << nfirst << " nstrips " << nstrips << " to tSensor[" << ilayer<< "][" << sensor << "]" <<endl;
            break;
         case 8:     // elements 0..1
            ilayer = 2;
            sensor = chip < 6? 0: 1;
            tSensor[ilayer][sensor].AddCluster(chip, nfirst, nstrips);
            //cout<< "added chip " << chip << " nfirst " << nfirst << " nstrips " << nstrips << " to tSensor[" << ilayer<< "][" << sensor << "]" <<endl;
            break;
         case 9:     // elements 2..3
            ilayer = 2;
            sensor = chip < 6? 0: 1;
            tSensor[ilayer][2+sensor].AddCluster(chip, nfirst, nstrips);
            //cout<< "added chip " << chip << " nfirst " << nfirst << " nstrips " << nstrips << " to tSensor[" << ilayer<< "][" << sensor << "]" <<endl;
            break;
         case 10:    // elements 0..1
            ilayer = 3;
            sensor = chip < 6? 0: 1;
            tSensor[ilayer][sensor].AddCluster(chip, nfirst, nstrips);
            //cout<< "added chip " << chip << " nfirst " << nfirst << " nstrips " << nstrips << " to tSensor[" << ilayer<< "][" << sensor << "]" <<endl;
            break;
         case 11:    // elements 2..3
            ilayer = 3;
            sensor = chip < 6? 0: 1;
            tSensor[ilayer][2+sensor].AddCluster(chip, nfirst, nstrips);
            //cout<< "added chip " << chip << " nfirst " << nfirst << " nstrips " << nstrips << " to tSensor[" << ilayer<< "][" << sensor << "]" <<endl;
            break;
      }
   }
};

/*
root -l pCTraw_Run_51.out-0.root
.L sensor.C+
sensor(345)
*/
void sensor(Int_t event, TTree* tree=0)
{
   if (!tree) tree = (TTree*) gDirectory->Get("t");
   if (!tree) {
      cout<< "Could not find tree t" <<endl;
      return;
   }

   RunHeader* runHeader = (RunHeader*) tree->GetUserInfo()->First();
   cout<< "runHeader->GetRun() = " << runHeader->GetRun() <<endl;
   time_t start_time = runHeader->GetTime();
   cout<< "run start time: " << std::ctime(&start_time);
   cout<< "program version is " << runHeader->GetVersion() <<endl;
   if (runHeader->GetTimeTag()) cout<< "event time tag was written out" <<endl;
   else cout<< "event time tag was not written out" <<endl;
   cout<<endl;

   PCTSensors pCTSensors;

   PCTEvent* pCTEvent = 0;
   tree->SetBranchAddress("event", &pCTEvent);

   if (tree->LoadTree(event) < 0) {
      cout<< "Could not load event " << event <<endl;
      return;
   }
   tree->GetEntry(event);

   for (int iFPGA=0; iFPGA<12; ++iFPGA)
   {
      const TrackerFPGA& trackerFPGA = pCTEvent->trackerFPGA[iFPGA];
      for (unsigned int ichip=0; ichip<trackerFPGA.numChips; ++ichip)
      {
         TrackerChip* trackerChip = (TrackerChip*) trackerFPGA.chips->At(ichip);
         for (int icluster=0; icluster<trackerChip->numClust; ++icluster) {
            Int_t strip = 64*trackerChip->address + (63-trackerChip->nfirst[icluster]);
            cout<< event << "\t iFPGA = " << std::setw(2) << iFPGA
               << " trackerChip->address = " << std::setw(2) << (unsigned) trackerChip->address
               << " nstrips[" << icluster << "] = " << (unsigned) trackerChip->nstrips[icluster]
               << " nfirst[" << icluster << "] = " << (unsigned) trackerChip->nfirst[icluster]
               << " strip = " << strip <<endl;
            // printf("%-8d iFPGA = %2d trackerChip->address = %2d nstrips[%d] = %d nfirst[%d] = %d strip = %d\n",
            //        event,iFPGA,(unsigned) trackerChip->address,icluster,(unsigned) trackerChip->nstrips[icluster],icluster,(unsigned) trackerChip->nfirst[icluster],strip);

            pCTSensors.AddCluster(iFPGA, trackerChip->address, trackerChip->nfirst[icluster], trackerChip->nstrips[icluster]);
         }
      }
   }
   pCTSensors.GetHits();
}

Int_t display(Int_t event, TTree* tree=0, const char* wname="event_display")
{
   if (!tree) tree = (TTree*) gDirectory->Get("t");
   if (!tree) {
      cout<< "Could not find tree t" <<endl;
      return -1;
   }

   // RunHeader* runHeader = (RunHeader*) tree->GetUserInfo()->First();
   // cout<< "runHeader->GetRun() = " << runHeader->GetRun() <<endl;
   // time_t start_time = runHeader->GetTime();
   // cout<< "run start time: " << std::ctime(&start_time);
   // cout<< "program version is " << runHeader->GetVersion() <<endl;
   // if (runHeader->GetTimeTag()) cout<< "event time tag was written out" <<endl;
   // else cout<< "event time tag was not written out" <<endl;
   // cout<<endl;

   PCTSensors pCTSensors;

   PCTEvent* pCTEvent = 0;
   tree->SetBranchAddress("event", &pCTEvent);

   if (tree->LoadTree(event) < 0) {
      cout<< "Could not load event " << event <<endl;
      return -1;
   }
   tree->GetEntry(event);

   for (int iFPGA=0; iFPGA<12; ++iFPGA)
   {
      const TrackerFPGA& trackerFPGA = pCTEvent->trackerFPGA[iFPGA];
      for (unsigned int ichip=0; ichip<trackerFPGA.numChips; ++ichip)
      {
         TrackerChip* trackerChip = (TrackerChip*) trackerFPGA.chips->At(ichip);
         for (int icluster=0; icluster<trackerChip->numClust; ++icluster) {
            Int_t strip = 64*trackerChip->address + (63-trackerChip->nfirst[icluster]);
            cout<< event << "\t iFPGA = " << std::setw(2) << iFPGA
               << " trackerChip->address = " << std::setw(2) << (unsigned) trackerChip->address
               << " nstrips[" << icluster << "] = " << (unsigned) trackerChip->nstrips[icluster]
               << " nfirst[" << icluster << "] = " << (unsigned) trackerChip->nfirst[icluster]
               << " strip = " << strip <<endl;
            // printf("%-8d iFPGA = %2d trackerChip->address = %2d nstrips[%d] = %d nfirst[%d] = %d strip = %d\n",
            //        event,iFPGA,(unsigned) trackerChip->address,icluster,(unsigned) trackerChip->nstrips[icluster],icluster,(unsigned) trackerChip->nfirst[icluster],strip);

            pCTSensors.AddCluster(iFPGA, trackerChip->address, trackerChip->nfirst[icluster], trackerChip->nstrips[icluster]);
         }
      }
   }
   pCTSensors.GetHits();

   Int_t n_v_hits = 0;
   Int_t n_t_hits = 0;

   // v-hits
   std::vector<Double_t> v_hits[4];       // vhits for all 4 layers
   for (int ilayer=0; ilayer<4; ++ilayer) {
      for (int isensor=0; isensor<2; ++isensor) {
         v_hits[ilayer].insert(v_hits[ilayer].end(), pCTSensors.vSensor[ilayer][isensor].hits_.begin(), pCTSensors.vSensor[ilayer][isensor].hits_.end());
         n_v_hits += pCTSensors.vSensor[ilayer][isensor].hits_.size();
      }
   }
   // t-hits
   std::vector<Double_t> t_hits[4];
   for (int ilayer=0; ilayer<4; ++ilayer) {
      for (int isensor=0; isensor<4; ++isensor) {
         t_hits[ilayer].insert(t_hits[ilayer].end(), pCTSensors.tSensor[ilayer][isensor].hits_.begin(), pCTSensors.tSensor[ilayer][isensor].hits_.end());
         n_t_hits += pCTSensors.tSensor[ilayer][isensor].hits_.size();
      }
   }

   //cout<< "n_v_hits = " << n_v_hits << " n_t_hits = " << n_t_hits <<endl;
   if (n_v_hits == 0 && n_t_hits == 0) return -1;

   // draw the event

   // for (int ilayer=0; ilayer<4; ++ilayer) {
   //    cout<< "v_hits[" << ilayer << "].size() = " << v_hits[ilayer].size() << " t_hits[" << ilayer << "].size() = " << t_hits[ilayer].size() <<endl;
   // }

   TCanvas* can = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(wname);
   if (can) gPad = can;
   else can = new TCanvas("event_display", wname, 700,500);

   can->DrawFrame(-300,-250, 300,250, Form("event %d",event));

   Double_t ut[4] = {-211.80, -161.80, 161.80, 211.80};
   Double_t uv[4] = {-217.70, -167.6,  167.6,  217.70};

   TLine line;
   // forward
   line.DrawLine(-217.70,-220.11, -217.70,215.24);
   line.DrawLine(-211.80,-220.11, -211.80,215.24);
   line.DrawLine(-167.60,-224.11, -167.60,211.24);
   line.DrawLine(-161.80,-224.11, -161.80,211.24);
   // back
   line.DrawLine( 217.70,-220.11,  217.70,215.24);
   line.DrawLine( 211.80,-220.11,  211.80,215.24);
   line.DrawLine( 167.60,-224.11,  167.60,211.24);
   line.DrawLine( 161.80,-224.11,  161.80,211.24);

   TMarker v_marker;
   v_marker.SetMarkerStyle(24);
   v_marker.SetMarkerColor(4);
   TMarker t_marker;
   t_marker.SetMarkerStyle(24);
   t_marker.SetMarkerColor(2);

   for (int ilayer=0; ilayer<4; ++ilayer) {
      for (unsigned ihit=0; ihit<v_hits[ilayer].size(); ++ihit) v_marker.DrawMarker(uv[ilayer], v_hits[ilayer][ihit]);
      for (unsigned ihit=0; ihit<t_hits[ilayer].size(); ++ihit) t_marker.DrawMarker(ut[ilayer], t_hits[ilayer][ihit]);
   }

   TLine v_line;
   v_line.SetLineColor(1);
   TLine t_line;
   t_line.SetLineColor(2);

   if (v_hits[0].size() == 1 && v_hits[3].size() == 1) v_line.DrawLine(uv[0],v_hits[0][0], uv[3],v_hits[3][0]);
   if (t_hits[0].size() == 1 && t_hits[3].size() == 1) t_line.DrawLine(ut[0],t_hits[0][0], ut[3],t_hits[3][0]);

   can->Update();

   tree->ResetBranchAddresses();

   return event;
}

void eloop(Int_t evtNo=0, TTree* tree=0, const char* wname="event_display")
{
   if (tree == 0) {
      tree = (TTree*) gDirectory->Get("t");
      if (!tree) {
         cout<< "Could not find tree \"t\" in the current directory " << gDirectory->GetName() <<endl;
         return;
      }

      RunHeader* runHeader = (RunHeader*) tree->GetUserInfo()->First();
      cout<< "runHeader->GetRun() = " << runHeader->GetRun() <<endl;
      time_t start_time = runHeader->GetTime();
      cout<< "run start time: " << std::ctime(&start_time);
      cout<< "program version is " << runHeader->GetVersion() <<endl;
      if (runHeader->GetTimeTag()) cout<< "event time tag was written out" <<endl;
      else cout<< "event time tag was not written out" <<endl;
      cout<<endl;
   }

   TTimer timer("gSystem->ProcessEvents();",50,kFALSE);  //-- process mouse events every 50 ms

   std::vector<std::string> commands;

   bool first_event = true;

   Int_t event_flag = -1;
   Int_t last_event_plot = -1;

   while (true)
   {
      timer.Start();                                     //-- start processing of mouse events

      std::string line;

      if (first_event) {
         //--evtNo;
         first_event = false;
      }
      else {
         // menu choise for the next event
         cout<< "<CR>: " << evtNo << ", -, Clone, Save, Quit, command: ";
         std::getline(cin, line);
      }

      timer.Stop();                                      //-- disable processing of mouse events

      // Possible inputs
      //
      // 0) <CR>: line.size() == 0:
      //    show next event
      // 1) number:
      //    interpret as event number: show this event
      // 2) not a number:
      //    can be one character or more than one character
      //       2.1) one character, line.size() == 1:
      //            interpret as menu command: C or S or Q
      //       2.2) more than one character, line.size() > 1:
      //            interpret as a ROOT command, try to execute

      // 0) <CR>
      if (line.size() == 0)
      {
         event_flag = display(evtNo,tree,wname);
         if (event_flag >= 0) last_event_plot = event_flag;
         ++evtNo;
         continue;
      }

      // 1) number
      std::stringstream ss(line);
      Int_t number = -1;
      if (ss >> number && ss.eof()) {
         evtNo = number;
         event_flag = display(evtNo,tree,wname);
         if (event_flag >= 0) last_event_plot = event_flag;
         ++evtNo;
         continue;
      }

      // 2) not a number
      if (line.size() == 1)
      {
         // 2.1) input was just one character: interpret as a menu command
         TCanvas* can = 0;
         switch (toupper(line[0])) {
            case '-':
               evtNo -= 2;          // previous event
               event_flag = display(evtNo,tree,wname);
               if (event_flag >= 0) last_event_plot = event_flag;
               ++evtNo;
               break;
            case 'C':
               can = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(wname);
               if (!can) {
                  cout<< "Could not find the canvas " << wname <<endl;
                  break;
               }
               can->DrawClone();
               break;
            case 'S':
               can = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(wname);
               if (!can) {
                  cout<< "Could not find the canvas " << wname <<endl;
                  break;
               }
               can->SaveAs(Form("%s_evt_%d.png", gDirectory->GetName(), last_event_plot));
               break;
            case 'Q':
               return;
            default: cout<< "No such key" <<endl;
         }
         continue;
      }
      else {
         // 2.2) input was more than one character: interpret as a ROOT command
         if (line == std::string(".q")) {
            cout<< "To terminate the ROOT session exit the macro first" <<endl;
            break;
         }
         else {
            if (unsigned(line[0]) == 27 && unsigned(line[1]) == 91) {
               // input was some arrow: ESC + '[' + letter A or B or C or D
               for (unsigned i=0; i<commands.size(); ++i) cout<< commands[i] <<endl;
               continue;
            }

            commands.push_back(line);
            gROOT->ProcessLine(line.c_str());
            continue;
         }
      }
   }
}
