// Andriy Zatserklyaniy, April 17, 2014

#ifndef Sensor_h
#define Sensor_h

#include "Track.h"
#include "Geometry.h"

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

ClassImp(SensorHit);

class Sensor {
   // sensor with 384 strips
public:
   static TClonesArray* poolSensorHit_;         //->
   //-- TClonesArray* poolSensorHit_;         //->
public:
   Bool_t debug_;
   Int_t sensorId_;                       // ex: v-sensor 1 in layer 3: 131, t-sensor 1 in layer 3: 231
   Int_t layer_;                          // info
   Double_t pitch_;
   Double_t u_;                           // coordinate of the u plane
   Double_t strip1_;                      // coordinate of the first strip
   Double_t dir_;                         // +1 or -1 direction of the strip numbers wrt axis
   std::map<Int_t,Int_t> clusters_;       // map<NFirst,NStrips> strip 0..383
   std::vector<SensorHit*> hits_;
   //std::vector<SensorHit*> hits_test_;
   std::vector<Int_t> dead_;              // dead channels
   std::vector<Int_t> noisy_;             // noisy channels
public:
   static void ClearPool() {
      if (poolSensorHit_) {
         // cout<< "ClearPool: poolSensorHit_->GetLast()+1 = " << poolSensorHit_->GetLast()+1 <<endl;
         //-- poolSensorHit_->Clear();
      }
   }
   void clear() {
      hits_.clear();
      clusters_.clear();
   }
   static void CreateHitPool() {
      if (!poolSensorHit_) {
         cout<< "Sensor::CreateHitPool" <<endl;
         poolSensorHit_ = new TClonesArray("SensorHit",1440);
      }
   }
   Sensor(): debug_(kFALSE), sensorId_(0), pitch_(0.228) {
      if (!poolSensorHit_) {
         //cout<< "Sensor::Sensor: create poolSensorHit_" <<endl;
         poolSensorHit_ = new TClonesArray("SensorHit",1440);
      }
   }
   virtual ~Sensor() {}
   void Set_strip1(Double_t strip1, Double_t dir) {
      strip1_ = strip1;
      dir_ = dir;
   }
   virtual void Set_strip1_dir_info(Double_t u, Double_t strip1, Double_t dir, Int_t layer, Int_t sensor) = 0;
   Double_t CenterOfGravity(Double_t strip) const {return strip1_ + dir_*strip*pitch_;}
   Double_t StripPosition(Int_t chip, Int_t nfirst) {
      // input data come directly from the readout:
      // chip:    0..12
      // nfirst:  0..63
      Int_t istrip = (chip % 6)*64 + (63 - nfirst);
      return CenterOfGravity(istrip);
   }
   void AddCluster(Int_t chip, Int_t nfirst, Int_t nstrips) {
      // input data come directly from the readout:
      // chip:    0..12
      // nfirst:  0..63
      // nstrips: 0..32
      Int_t istrip = (chip % 6)*64 + (63 - nfirst);
      clusters_[istrip-nstrips+1] = nstrips;
      // for the test purposes: create SensorHit and add to hits_test_ or modify previous hit if they are the same hit as in 345
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
         Double_t pos = CenterOfGravity(nfirst + 0.5*(nstrips-1));
         // SensorHit* sensorHit = new SensorHit(sensorId_, nfirst, nstrips, u_, pos);
         SensorHit* sensorHit = new ((*poolSensorHit_)[poolSensorHit_->GetLast()+1]) SensorHit(sensorId_, nfirst, nstrips, u_, pos);
         hits_.push_back(sensorHit);
         if (debug_) cout<< "sensorId: " << sensorId_ << " layer " << layer_ << " sensor " << sensorId_ << " nfirst = " << std::setw(3) << nfirst << " nstrips = " << nstrips << " coordinate = " << CenterOfGravity(nfirst) <<endl;
      }
   }
};

class VSensor: public Sensor {
public:
   VSensor(): Sensor() {}
   void Set_strip1_dir_info(Double_t u, Double_t strip1, Double_t dir, Int_t layer, Int_t sensor) {
      sensorId_ = 100 + layer*10 + sensor;
      u_ = u;
      strip1_ = strip1;
      dir_ = dir;
      layer_ = layer;
   }
};

class TSensor: public Sensor {
public:
   TSensor(): Sensor() {}
   void Set_strip1_dir_info(Double_t u, Double_t strip1, Double_t dir, Int_t layer, Int_t sensor) {
      sensorId_ = 200 + layer*10 + sensor;
      u_ = u;
      strip1_ = strip1;
      dir_ = dir;
      layer_ = layer;
   }
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
   const Geometry* geometry_;
   VSensor vSensor[4][2];     // [layer][sensor]
   TSensor tSensor[4][4];     // [layer][sensor]
   // hits
   std::vector<SensorHit*> v_hits[4];       // v_hits for all 4 layers
   std::vector<SensorHit*> t_hits[4];       // t_hits for all 4 layers

   void GetHits() {
      for (int ilayer=0; ilayer<4; ++ilayer) for (int isensor=0; isensor<2; ++isensor) vSensor[ilayer][isensor].GetHits();
      for (int ilayer=0; ilayer<4; ++ilayer) for (int isensor=0; isensor<4; ++isensor) tSensor[ilayer][isensor].GetHits();
      // v-hits
      for (int ilayer=0; ilayer<4; ++ilayer) {
         for (int isensor=0; isensor<2; ++isensor) {
            v_hits[ilayer].insert(v_hits[ilayer].end(), vSensor[ilayer][isensor].hits_.begin(), vSensor[ilayer][isensor].hits_.end());
         }
      }
      // t-hits
      for (int ilayer=0; ilayer<4; ++ilayer) {
         for (int isensor=0; isensor<4; ++isensor) {
            t_hits[ilayer].insert(t_hits[ilayer].end(), tSensor[ilayer][isensor].hits_.begin(), tSensor[ilayer][isensor].hits_.end());
         }
      }
   }
   void clear() {
      for (int ilayer=0; ilayer<4; ++ilayer) {
         v_hits[ilayer].clear();
         t_hits[ilayer].clear();
      }
      for (int ilayer=0; ilayer<4; ++ilayer) {
         for (int isensor=0; isensor<2; ++isensor) {
            vSensor[ilayer][isensor].clear();
            tSensor[ilayer][isensor].clear();
            tSensor[ilayer][isensor+2].clear();
            tSensor[ilayer][isensor+2].clear();
         }
      }
   }
   ~PCTSensors() {
      for (int ilayer=0; ilayer<4; ++ilayer) {
         v_hits[ilayer].clear();
         t_hits[ilayer].clear();
      }
   }
   PCTSensors(const Geometry* geometry): geometry_(geometry)
   {
      //cout<< "PCTSensors::PCTSensors" <<endl;
      //Sensor::CreateHitPool();

      const Double_t pitch = 0.228;
      Int_t layer;                  // 0..3
      Int_t board;                  // board production number: V-boards: 0..6, T-boards 0..7
      Int_t dir;                    // chip direction wrt corresponding axis
      Int_t sensor;                 // sensor number: V-board: 0..1, T-board: 0..3

      //
      //  V-sensors
      //

      /// Double_t vPin[4] = {0.055,0.055,0.055,0.055};
      /// int vBoard[4] = {6,4,2,3};  	      // fpga to V-board translation; this will change if spares are swapped
      /// Double_t firstStripV[7][2] = {		// Distance from the alignment pin to the first strip
      ///    {-43.7193, -43.716},             // Board V0 doesn't exist
      ///    {-43.7193, -43.716},
      ///    {-43.7193, -43.716},
      ///    {-43.7193, -43.716},
      ///    {-43.7193, -43.716},	            // These numbers are actually from V4
      ///    {-43.7193, -43.716},
      ///    {-43.5855, -43.5865}             // these are actually from T6
      /// };

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

      /// Double_t uv[4] = {-217.7, -167.6,  166.8,  216.9};
      /// Double_t ut[4] = {-211.8, -161.8,  161.0,  211.0};
      // //
      // // Before the run 51 we moved the upstream cassette to 2 inch upstream
      // // This is equivalent to increase in abs value of all distances to 25.4 mm
      // //
      // for (int ilayer=0; ilayer<4; ++ilayer) {
      //    Double_t shift = 25.4;  // mm
      //    Double_t sign;
      //    sign = uv[ilayer] < 0? -1.: 1.;
      //    uv[ilayer] += sign*shift;
      //    sign = ut[ilayer] < 0? -1.: 1.;
      //    ut[ilayer] += sign*shift;
      // }

      Double_t firstStripVabs[7][2];      // absolute coordinate of the first strip
      sensor = 0;
      for (int iboard=0; iboard<7; ++iboard) firstStripVabs[iboard][sensor] = geometry_->vPin_[0] + geometry_->firstStripV_[iboard][sensor] + 383*pitch;  // chip 0, strip 63
      sensor = 1;
      for (int iboard=0; iboard<7; ++iboard) firstStripVabs[iboard][sensor] = geometry_->vPin_[0] + geometry_->firstStripV_[iboard][sensor];  // chip 6, strip 63
      // for (int iboard=0; iboard<7; ++iboard) {
      //    for (int iside=0; iside<2; ++iside) {
      //       cout<< firstStripVabs[iboard][iside] << "\t ";
      //    }
      //    cout<<endl;
      // }

      layer = 0;
      board = geometry_->vBoardLayer_[layer];  dir = -1; sensor = 0;   //1*100 + layer*10 + 0;
      vSensor[layer][sensor].Set_strip1_dir_info(geometry_->uv_[layer], firstStripVabs[board][sensor], dir, layer, sensor);
      board = geometry_->vBoardLayer_[layer];  dir = +1; sensor = 1;   //1*100 + layer*10 + 1;
      vSensor[layer][sensor].Set_strip1_dir_info(geometry_->uv_[layer], firstStripVabs[board][sensor], dir, layer, sensor);
      //
      layer = 1;  // identical to the layer 0
      board = geometry_->vBoardLayer_[layer];  dir = -1; sensor = 0;   //1*100 + layer*10 + 0;
      vSensor[layer][sensor].Set_strip1_dir_info(geometry_->uv_[layer], firstStripVabs[board][sensor], dir, layer, sensor);
      board = geometry_->vBoardLayer_[layer];  dir = +1; sensor = 1;   //1*100 + layer*10 + 1;
      vSensor[layer][sensor].Set_strip1_dir_info(geometry_->uv_[layer], firstStripVabs[board][sensor], dir, layer, sensor);
      //
      layer = 2;
      board = geometry_->vBoardLayer_[layer];  dir = -1; sensor = 0;   //1*100 + layer*10 + 0;
      vSensor[layer][sensor].Set_strip1_dir_info(geometry_->uv_[layer], firstStripVabs[board][sensor], dir, layer, sensor);
      board = geometry_->vBoardLayer_[layer];  dir = +1; sensor = 1;   //1*100 + layer*10 + 1;
      vSensor[layer][sensor].Set_strip1_dir_info(geometry_->uv_[layer], firstStripVabs[board][sensor], dir, layer, sensor);
      //
      layer = 3;  // identical to the layer 2
      board = geometry_->vBoardLayer_[layer];  dir = -1; sensor = 0;   //1*100 + layer*10 + 0;
      vSensor[layer][sensor].Set_strip1_dir_info(geometry_->uv_[layer], firstStripVabs[board][sensor], dir, layer, sensor);
      board = geometry_->vBoardLayer_[layer];  dir = +1; sensor = 1;   //1*100 + layer*10 + 1;
      vSensor[layer][sensor].Set_strip1_dir_info(geometry_->uv_[layer], firstStripVabs[board][sensor], dir, layer, sensor);

      //
      //  T-sensors
      //

      /// Double_t tPin[4] = {215.24, 211.25, -211.25, -215.24};   // T coordinate of alignment pin per layer
      /// //Double_t tDir[4] = {-1.0, -1.0, 1.0, 1.0};               // T direction of increasing strip number per layer

      /// Int_t tBoard[4] = {5, 4, 1, 3};     // T-layer to physical T board translation. This will change if spares are swapped in!
      /// Double_t firstStripT[7][4] = {      // First strip location rel to pin for each sensor on each physical board
      ///    {-999., -999., -999., -999.},		// Board 0 doesn't exist.  We manufactured 6 boards.
      ///    {38.58, 126.85, 215.11, 303.37},
      ///    {38.58, 126.85, 215.11, 303.37},
      ///    {38.58, 126.85, 215.11, 303.37},
      ///    {38.58, 126.85, 215.11, 303.37}, // These numbers actually come from board T4
      ///    {38.62, 126.90, 215.16, 303.41}, // these numbers actually come from board T5
      ///    {38.58, 126.85, 215.11, 303.37} 
      /// };

      dir = -1;
      for (int ilayer=0; ilayer<2; ++ilayer) {
         Double_t pin = geometry_->tPin_[ilayer];
         board = geometry_->tBoardLayer_[ilayer];                         // production number of the board for this layer
         for (int isensor=0; isensor<4; ++isensor) {
            //Int_t sensorId = 2*100 + ilayer*10 + isensor;
            tSensor[ilayer][isensor].Set_strip1_dir_info(geometry_->ut_[ilayer], pin + dir*geometry_->firstStripT_[board][isensor], dir, ilayer, isensor);
         }
      }

      dir = +1;
      for (int ilayer=2; ilayer<4; ++ilayer) {
         Double_t pin = geometry_->tPin_[ilayer];
         board = geometry_->tBoardLayer_[ilayer];                         // production number of the board for this layer
         for (int isensor=0; isensor<4; ++isensor) {
            //Int_t sensorId = 2*100 + ilayer*10 + isensor;
            tSensor[ilayer][isensor].Set_strip1_dir_info(geometry_->ut_[ilayer], pin + dir*geometry_->firstStripT_[board][isensor], dir, ilayer, isensor);
         }
      }
   }
   void AddCluster(Int_t FPGA, Int_t chip, Int_t nfirst, Int_t nstrips) {
      int sensor;
      int ilayer;
      switch (FPGA) {
         //
         // v-layers
         //
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

            //
            // t-layers
            //
         case 4:     // elements 0..1
            ilayer = 0;
            sensor = chip < 6? 0: 1;
            tSensor[ilayer][sensor].AddCluster(chip, nfirst, nstrips);
            //cout<< "added chip " << chip << " nfirst " << nfirst << " nstrips " << nstrips << " to tSensor[" << ilayer<< "][" << sensor << "]" <<endl;
            break;
         case 5:     // elements 2..3
            ilayer = 0;
            sensor = chip < 6? 0: 1;                                          // with offset 2
            tSensor[ilayer][2+sensor].AddCluster(chip, nfirst, nstrips);      // NB: 2 + (0 or 1) = (2 or 3)
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

#endif   // Sensor_h
