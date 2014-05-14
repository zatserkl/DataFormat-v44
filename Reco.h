// Andriy Zatserklyaniy, April 17, 2014:  production release

#ifndef Reco_h
#define Reco_h

#include "DataFormat.h"
#include "Track.h"
#include "Sensor.h"
#include "Geometry.h"

#include <iostream>
#include <vector>
#include <map>
#include <list>

using std::cout;    using std::endl;

ClassImp(Hit2D);
ClassImp(SensorHit);
ClassImp(CRay2D);
//ClassImp(CRay);
ClassImp(Track2D);
//ClassImp(Track);
ClassImp(SuperTrack2D);
ClassImp(SuperTrack);

class Reco {
public:
   Bool_t debug_;
   Float_t deltaT_;
   const Geometry* geometry_;
   PCTSensors* pCTSensors_;
   const PCTEvent* pCTEvent_;
   Int_t event_;
   std::list<const Track2D*> tin_;
   std::list<const Track2D*> tout_;
   std::list<const Track2D*> vin_;
   std::list<const Track2D*> vout_;

   std::list<const SuperTrack2D*> vSuperTracks_;
   std::list<const SuperTrack2D*> tSuperTracks_;

   std::list<const SuperTrack*> superTracks_;

   static TClonesArray* poolTrack2D_;                 //->
   static TClonesArray* poolSuperTrack2D_;            //->
   static TClonesArray* poolSuperTrack_;              //->
public:
   Reco(const Geometry* geometry, PCTSensors* pCTSensors, const PCTEvent* pCTEvent, Int_t event, bool debug=kFALSE): debug_(debug)
      , geometry_(geometry), pCTSensors_(pCTSensors), pCTEvent_(pCTEvent)
      , event_(event)
   {
      //cout<< "Reco::Reco" <<endl;

      deltaT_ = pCTEvent_->deltaT;
      //cout<< "Reco::Reco: create pool" <<endl;
      //Sensor::CreateHitPool();
      Sensor::ClearPool();

      if (!poolTrack2D_) poolTrack2D_ = new TClonesArray("Track2D", 1024);
      if (!poolSuperTrack2D_) poolSuperTrack2D_ = new TClonesArray("SuperTrack2D", 128);
      if (!poolSuperTrack_) poolSuperTrack_ = new TClonesArray("SuperTrack", 16);

      //cout<< "Sensor::poolSensorHit_->GetLast()+1 = " << Sensor::poolSensorHit_->GetLast()+1 <<endl;
      //cout<< "poolTrack2D_->GetLast()+1 = " << poolTrack2D_->GetLast()+1 <<endl;
      //cout<< "poolSuperTrack2D_->GetLast()+1 = " << poolSuperTrack2D_->GetLast()+1 <<endl;
      //cout<< "poolSuperTrack_->GetLast()+1 = " << poolSuperTrack_->GetLast()+1 <<endl;

      //----------- pCTSensors_ = new PCTSensors();
      pCTSensors_->clear();

      // get hits for this event
      for (int iFPGA=0; iFPGA<12; ++iFPGA)
      {
         const TrackerFPGA& trackerFPGA = pCTEvent_->trackerFPGA[iFPGA];
         for (unsigned int ichip=0; ichip<trackerFPGA.numChips; ++ichip)
         {
            TrackerChip* trackerChip = (TrackerChip*) trackerFPGA.chips->At(ichip);
            for (int icluster=0; icluster<trackerChip->numClust; ++icluster) {
               if (debug_) {
                  Int_t strip = 64*trackerChip->address + (63-trackerChip->nfirst[icluster]);
                  cout<< event_ << "\t iFPGA = " << std::setw(2) << iFPGA
                  << " trackerChip->address = " << std::setw(2) << (unsigned) trackerChip->address
                  << " nstrips[" << icluster << "] = " << (unsigned) trackerChip->nstrips[icluster]
                  << " nfirst[" << icluster << "] = " << (unsigned) trackerChip->nfirst[icluster]
                  << " strip = " << strip <<endl;
                  // printf("%-8d iFPGA = %2d trackerChip->address = %2d nstrips[%d] = %d nfirst[%d] = %d strip = %d\n",
                  //        event,iFPGA,(unsigned) trackerChip->address,icluster,(unsigned) trackerChip->nstrips[icluster],icluster,(unsigned) trackerChip->nfirst[icluster],strip);
               }
               pCTSensors_->AddCluster(iFPGA, trackerChip->address, trackerChip->nfirst[icluster], trackerChip->nstrips[icluster]);
            }
         }
      }
      pCTSensors_->GetHits();

      Int_t n_v_hits = 0;
      for (unsigned ilayer=0; ilayer<4; ++ilayer) n_v_hits += pCTSensors_->v_hits[ilayer].size();
      Int_t n_t_hits = 0;
      for (unsigned ilayer=0; ilayer<4; ++ilayer) n_t_hits += pCTSensors_->t_hits[ilayer].size();

      if (debug_) cout<< "n_v_hits = " << n_v_hits << " n_t_hits = " << n_t_hits <<endl;
   }
   ~Reco() {
      //cout<< "Reco::~Reco" <<endl;

      //for (std::list<const Track2D*>::const_iterator it=tin_.begin(); it!=tin_.end(); ++it) delete *it;
      //for (std::list<const Track2D*>::const_iterator it=tout_.begin(); it!=tout_.end(); ++it) delete *it;
      //for (std::list<const Track2D*>::const_iterator it=vin_.begin(); it!=vin_.end(); ++it) delete *it;
      //for (std::list<const Track2D*>::const_iterator it=vout_.begin(); it!=vout_.end(); ++it) delete *it;
      if (poolTrack2D_) poolTrack2D_->Clear();
      if (poolTrack2D_) poolSuperTrack2D_->Clear();
      if (poolTrack2D_) poolSuperTrack_->Clear();

      //for (std::list<const SuperTrack2D*>::const_iterator it=vSuperTracks_.begin(); it!=vSuperTracks_.end(); ++it) delete *it;
      //for (std::list<const SuperTrack2D*>::const_iterator it=tSuperTracks_.begin(); it!=tSuperTracks_.end(); ++it) delete *it;

      //for (std::list<const SuperTrack*>::const_iterator it=superTracks_.begin(); it!=superTracks_.end(); ++it) delete *it;

      //------------ delete pCTSensors_;
   }
   void GenerateTracks2D()
   {
      // T-sensors

      for (unsigned ihit1=0; ihit1<pCTSensors_->t_hits[0].size(); ++ihit1) {
         const SensorHit* hit1 = pCTSensors_->t_hits[0][ihit1];
         if (debug_) cout<< "hit1: " << *hit1 <<endl;
         for (unsigned ihit2=0; ihit2<pCTSensors_->t_hits[1].size(); ++ihit2) {
            const SensorHit* hit2 = pCTSensors_->t_hits[1][ihit2];
            if (debug_) cout<< "hit2: " << *hit2 <<endl;
            //-- Track2D* cRay2d = new Track2D(hit1, hit2);
            Track2D* cRay2d = new ((*poolTrack2D_)[poolTrack2D_->GetLast()+1]) Track2D(hit1, hit2);
            tin_.push_back(cRay2d);
         }
      }

      for (unsigned ihit1=0; ihit1<pCTSensors_->t_hits[2].size(); ++ihit1) {
         const SensorHit* hit1 = pCTSensors_->t_hits[2][ihit1];
         if (debug_) cout<< "hit1: " << *hit1 <<endl;
         for (unsigned ihit2=0; ihit2<pCTSensors_->t_hits[3].size(); ++ihit2) {
            const SensorHit* hit2 = pCTSensors_->t_hits[3][ihit2];
            if (debug_) cout<< "hit2: " << *hit2 <<endl;
            //-- Track2D* cRay2d = new Track2D(hit1, hit2);
            Track2D* cRay2d = new ((*poolTrack2D_)[poolTrack2D_->GetLast()+1]) Track2D(hit1, hit2);
            tout_.push_back(cRay2d);
         }
      }

      if (debug_) {
         cout<< "tin_" <<endl;
         for (std::list<const Track2D*>::const_iterator it=tin_.begin(); it!=tin_.end(); ++it) {
            cout<< "it->x_ = " << (*it)->x_ <<endl;
         }
      }

      if (debug_) {
         cout<< "tout_" <<endl;
         for (std::list<const Track2D*>::const_iterator it=tout_.begin(); it!=tout_.end(); ++it) {
            cout<< "it->x_ = " << (*it)->x_ <<endl;
         }
      }

      // V-sensors

      for (unsigned ihit1=0; ihit1<pCTSensors_->v_hits[0].size(); ++ihit1) {
         SensorHit* hit1 = pCTSensors_->v_hits[0][ihit1];
         if (debug_) cout<< "hit1: " << *hit1 <<endl;
         for (unsigned ihit2=0; ihit2<pCTSensors_->v_hits[1].size(); ++ihit2) {
            SensorHit* hit2 = pCTSensors_->v_hits[1][ihit2];
            if (debug_) cout<< "hit2: " << *hit2 <<endl;
            //-- Track2D* cRay2d = new Track2D(hit1, hit2);
            Track2D* cRay2d = new ((*poolTrack2D_)[poolTrack2D_->GetLast()+1]) Track2D(hit1, hit2);
            vin_.push_back(cRay2d);
         }
      }

      for (unsigned ihit1=0; ihit1<pCTSensors_->v_hits[2].size(); ++ihit1) {
         SensorHit* hit1 = pCTSensors_->v_hits[2][ihit1];
         if (debug_) cout<< "hit1: " << *hit1 <<endl;
         for (unsigned ihit2=0; ihit2<pCTSensors_->v_hits[3].size(); ++ihit2) {
            SensorHit* hit2 = pCTSensors_->v_hits[3][ihit2];
            if (debug_) cout<< "hit2: " << *hit2 <<endl;
            //-- Track2D* cRay2d = new Track2D(hit1, hit2);
            Track2D* cRay2d = new ((*poolTrack2D_)[poolTrack2D_->GetLast()+1]) Track2D(hit1, hit2);
            vout_.push_back(cRay2d);
         }
      }

      if (debug_) {
         cout<< "vin_" <<endl;
         for (std::list<const Track2D*>::const_iterator it=vin_.begin(); it!=vin_.end(); ++it) {
            cout<< "it->x_ = " << (*it)->x_ <<endl;
         }
      }

      if (debug_) {
         cout<< "vout_" <<endl;
         for (std::list<const Track2D*>::const_iterator it=vout_.begin(); it!=vout_.end(); ++it) {
            cout<< "it->x_ = " << (*it)->x_ <<endl;
         }
      }
   }
   void GenerateVSuperTracks2D(Double_t rmax=10.)
   {
      // start from V-senosors: as a max we can reconstruct two tracks and only if they are from different V-sensors
      for (std::list<const Track2D*>::const_iterator it=vin_.begin(); it!=vin_.end(); ++it) {
         const Track2D* itrack = *it;
         for (std::list<const Track2D*>::const_iterator ot=vout_.begin(); ot!=vout_.end(); ++ot) {
            const Track2D* otrack = *ot;
            //-- SuperTrack2D* superTrack = new SuperTrack2D(itrack, otrack);
            SuperTrack2D* superTrack = new ((*poolSuperTrack2D_)[poolSuperTrack2D_->GetLast()+1]) SuperTrack2D(itrack, otrack);
            vSuperTracks_.push_back(superTrack);
         }
      }

      // apply filter on the distance between the hits in the plane u = 0

      std::map<Double_t, const SuperTrack2D*> mapCloseTracks;            // for the filter on the distance in the plane u = 0

      for (std::list<const SuperTrack2D*>::const_iterator it=vSuperTracks_.begin(); it!=vSuperTracks_.end(); ++it) {
         const SuperTrack2D* superTrack = *it;
         mapCloseTracks[superTrack->Distance()] = superTrack;
      }
      if (debug_) cout<< "GenerateSuperTracks2D: vSuperTracks_.size() = " << vSuperTracks_.size() <<endl;

      if (debug_) {
         cout<< "resulting distance map from SuperTracks2D for V-board" <<endl;
         for (std::map<Double_t, const SuperTrack2D*>::const_iterator it=mapCloseTracks.begin(); it!=mapCloseTracks.end(); ++it) {
            Double_t distance = it->first;
            cout<< std::distance<std::map<Double_t, const SuperTrack2D*>::const_iterator>(mapCloseTracks.begin(), it) << "\t distance = " << distance <<endl;
         }
      }

      // 1) loop over map to remove the tracks with higher distance
      // 2) check for hits overlap in the (small number) of passed tracks

      if (debug_) cout<< "loop over the map to remove tracks with distance above the r" <<endl;

      //--const-- for (std::map<Double_t, const SuperTrack2D*>::const_iterator it=mapCloseTracks.begin(); it!=mapCloseTracks.end();)
      for (std::map<Double_t, const SuperTrack2D*>::iterator it=mapCloseTracks.begin(); it!=mapCloseTracks.end();)
      {
         if (it->first < rmax) ++it;
         else {
            if (debug_) cout<< "remove track with distance " << it->first <<endl;
            vSuperTracks_.remove(it->second);
            //-- it = mapCloseTracks.erase(it);
            mapCloseTracks.erase(it++);
         }
      }

      if (debug_) cout<< "GenerateSuperTracks2D: vSuperTracks_.size() = " << vSuperTracks_.size() <<endl;

      // make sure that the rest of tracks do not share the same hits

      //--const-- for (std::map<Double_t, const SuperTrack2D*>::const_iterator it=mapCloseTracks.begin(); it!=mapCloseTracks.end(); ++it)
      for (std::map<Double_t, const SuperTrack2D*>::iterator it=mapCloseTracks.begin(); it!=mapCloseTracks.end(); ++it)
      {
         //--const-- std::map<Double_t, const SuperTrack2D*>::const_iterator next = it;
         std::map<Double_t, const SuperTrack2D*>::iterator next = it;
         ++next;
         while (next != mapCloseTracks.end()) {
            if (next->second->SharedHits(it->second)) {
               vSuperTracks_.remove(next->second);
               //-- next = mapCloseTracks.erase(next);
               mapCloseTracks.erase(next++);
            }
            else ++next;
         }
      }
   }
   void GenerateTSuperTracks2D(Double_t rmax=10.)
   {
      for (std::list<const Track2D*>::const_iterator it=tin_.begin(); it!=tin_.end(); ++it) {
         const Track2D* itrack2D = *it;
         for (std::list<const Track2D*>::const_iterator ot=tout_.begin(); ot!=tout_.end(); ++ot) {
            const Track2D* otrack2D = *ot;
            //-- SuperTrack2D* superTrack2D = new SuperTrack2D(itrack2D, otrack2D);
            SuperTrack2D* superTrack2D = new ((*poolSuperTrack2D_)[poolSuperTrack2D_->GetLast()+1]) SuperTrack2D(itrack2D, otrack2D);
            tSuperTracks_.push_back(superTrack2D);
         }
      }

      // apply filter on the distance between the hits in the plane u = 0

      std::map<Double_t, const SuperTrack2D*> mapCloseTracks;            // for the filter on the distance in the plane u = 0

      for (std::list<const SuperTrack2D*>::const_iterator it=tSuperTracks_.begin(); it!=tSuperTracks_.end(); ++it) {
         const SuperTrack2D* superTrack = *it;
         mapCloseTracks[superTrack->Distance()] = superTrack;
      }
      if (debug_) cout<< "GenerateSuperTracks2D: tSuperTracks_.size() = " << tSuperTracks_.size() <<endl;

      if (debug_) {
         cout<< "resulting distance map from SuperTracks2D for T-board" <<endl;
         for (std::map<Double_t, const SuperTrack2D*>::const_iterator it=mapCloseTracks.begin(); it!=mapCloseTracks.end(); ++it) {
            Double_t distance = it->first;
            cout<< std::distance<std::map<Double_t, const SuperTrack2D*>::const_iterator>(mapCloseTracks.begin(), it) << "\t distance = " << distance <<endl;
         }
      }

      // 1) loop over map to remove the tracks with higher distance
      // 2) check for hits overlap in the (small number) of passed tracks

      if (debug_) cout<< "loop over the map to remove tracks with distance above the r" <<endl;

      //--const-- for (std::map<Double_t, const SuperTrack2D*>::const_iterator it=mapCloseTracks.begin(); it!=mapCloseTracks.end();)
      for (std::map<Double_t, const SuperTrack2D*>::iterator it=mapCloseTracks.begin(); it!=mapCloseTracks.end();)
      {
         if (it->first < rmax) ++it;
         else {
            if (debug_) cout<< "remove track with distance " << it->first <<endl;
            tSuperTracks_.remove(it->second);
            //-- it = mapCloseTracks.erase(it);
            mapCloseTracks.erase(it++);
         }
      }

      if (debug_) cout<< "GenerateSuperTracks2D: tSuperTracks_.size() = " << tSuperTracks_.size() <<endl;

      // make sure that the rest of tracks do not share the same hits

      //--const-- for (std::map<Double_t, const SuperTrack2D*>::const_iterator it=mapCloseTracks.begin(); it!=mapCloseTracks.end(); ++it)
      for (std::map<Double_t, const SuperTrack2D*>::iterator it=mapCloseTracks.begin(); it!=mapCloseTracks.end(); ++it)
      {
         //--const-- std::map<Double_t, const SuperTrack2D*>::const_iterator next = it;
         std::map<Double_t, const SuperTrack2D*>::iterator next = it;
         ++next;
         while (next != mapCloseTracks.end()) {
            if (next->second->SharedHits(it->second)) {
               tSuperTracks_.remove(next->second);
               //-- next = mapCloseTracks.erase(next);
               mapCloseTracks.erase(next++);
            }
            else ++next;
         }
      }
   }
   void GenerateSuperTracks2D(Double_t rmax=10.) {
      GenerateVSuperTracks2D(rmax);
      GenerateTSuperTracks2D(rmax);
   }
   void GenerateSuperTracks() {
      if (debug_) cout<< "GenerateSuperTracks: vSuperTracks_.size() = " << vSuperTracks_.size() << " tSuperTracks_.size() = " << tSuperTracks_.size() <<endl;
      // inspect the number of T-super tracks 2D: maximum number is 2
      if (tSuperTracks_.size() > 2) return;

      // inspect the number of V-super tracks 2D: maximum number is 2.
      // In case of 2 tracks they should be from different halves of the V-board (at least in one layer)
      if (vSuperTracks_.size() > 2) return;
      if (vSuperTracks_.size() == 2) {
         const SuperTrack2D* track1 = vSuperTracks_.front();
         const SuperTrack2D* track2 = vSuperTracks_.back();
         bool same = true
            && track1->itrack2D_->hit1_->sensorId_ == track2->itrack2D_->hit1_->sensorId_
            && track1->itrack2D_->hit2_->sensorId_ == track2->itrack2D_->hit2_->sensorId_
            && track1->otrack2D_->hit1_->sensorId_ == track2->otrack2D_->hit1_->sensorId_
            && track1->otrack2D_->hit2_->sensorId_ == track2->otrack2D_->hit2_->sensorId_
         ;
         if (same) {
            if (debug_) cout<< "GenerateSuperTracks: same = " << same <<endl;
            return;
         }
      }

      // at this point we have as max two t-superTracks and no more than v-superTracks
      // in principle, we should combine them to match V-board halves with the t-tracks
      //
      // leave that for future
      //
      if (vSuperTracks_.size() != 1) return;
      if (tSuperTracks_.size() != 1) return;

      // match the halves

      //-- SuperTrack* superTrack = new SuperTrack(vSuperTracks_.front(), tSuperTracks_.front());
      SuperTrack* superTrack = new ((*poolSuperTrack_)[poolSuperTrack_->GetLast()+1]) SuperTrack(vSuperTracks_.front(), tSuperTracks_.front());
      superTracks_.push_back(superTrack);
   }
};

class RecoEvent: public TObject {
public:
   Bool_t ok;                       // error flag
   Float_t deltaT;
   TClonesArray* track;             //-> 
   Int_t nt;                        // the number of super tracks
   Float_t a[5];                    // Energy detector channels
   Float_t ped[5];
   Float_t sample[5][16];           // to plot e.g. channel 1:    r->Draw("sample[1][]:Iteration$","Entry$==0") 
   RecoEvent(): TObject(), ok(kTRUE), nt(0) {
      track = new TClonesArray("SuperTrack");
   }
   ~RecoEvent() {delete track;}
   void clear() {
      deltaT = 0;
      track->Clear();
      nt = 0;
      for (int i=0; i<5; ++i) {
         a[i] = 0;
         ped[i] = 0;
      }
      for (int ichan=0; ichan<5; ++ichan) for (int isample=0; isample<16; ++isample) sample[ichan][isample] = 0;
   }
   Float_t SampleSum(Int_t chan, Int_t nfront, Int_t ntail, Double_t pedestal) const {
      if (chan < 0 || chan > 4) return 0;
      // find a position of the maximum
      Int_t imax = 0;
      // assumes that the number of samples is 16
      Float_t sum = sample[chan][imax];
      for (int isample=0; isample<16; ++isample) {
         sum += sample[chan][isample];
         if (sample[chan][isample] > sample[chan][imax]) {
            imax = isample;
         }
      }
      if (sum == 0) return 0;          // there are no samples in this event
      Int_t n1 = imax - nfront;
      if (n1 < 0) n1 = 0;
      Int_t n2 = imax + ntail;
      if (n2 > 15) n2 = 15;

      sum = 0;
      for (int isample=n1; isample<=n2; ++isample) sum += sample[chan][isample];
      sum -= (n2 - n1 + 1)*pedestal;

      return sum;
   }
   void Extract(const Reco& reco) {
      clear();
      deltaT = reco.deltaT_;
      for (std::list<const SuperTrack*>::const_iterator it=reco.superTracks_.begin(); it!=reco.superTracks_.end(); ++it) {
         const SuperTrack* superTrack = *it;
         new ((*track)[track->GetLast()+1]) SuperTrack(*superTrack);
         ++nt;
      }
      // Energy detector
      //a[0] = reco.pCTEvent_->energyBoard[0].pulse[0];
      //a[1] = reco.pCTEvent_->energyBoard[0].pulse[1];
      //a[2] = reco.pCTEvent_->energyBoard[0].pulse[2];
      //a[3] = reco.pCTEvent_->energyBoard[1].pulse[0];
      //a[4] = reco.pCTEvent_->energyBoard[1].pulse[1];

      /////////////////////////////
      int brd=0;
      ///// //--orig int enrgTag0= thisEvent->Event->Board[0].enrgTag;
      ///// //--orig int enrgTag1= thisEvent->Event->Board[1].enrgTag;
      ///// //		if (enrgTag0 != enrgTag1) cout << "enrg tag mismatch: " << enrgTag0 << " vs " << enrgTag1;
      if (reco.pCTEvent_->energyBoard[brd].numChan>0 || reco.pCTEvent_->energyBoard[brd].numSamples>0) {
        /////   //			if (thisEvent->Event->Board[brd].enrgTag != evtNum % 4) cout << "tag mismatch, energy tag=" << thisEvent->Event->Board[brd].enrgTag << "\n";
        if (reco.pCTEvent_->energyBoard[brd].reduced) {
          a[0] = reco.pCTEvent_->energyBoard[brd].pulse[0];  ped[0] = reco.pCTEvent_->energyBoard[brd].pedestal[0];
          a[1] = reco.pCTEvent_->energyBoard[brd].pulse[1];	ped[1] = reco.pCTEvent_->energyBoard[brd].pedestal[1];
          a[2] =	reco.pCTEvent_->energyBoard[brd].pulse[2];	ped[2] = reco.pCTEvent_->energyBoard[brd].pedestal[2];
        } else {
          a[0] = 0.; a[1] = 0.; a[2] = 0.;
          /// enrgSamp *thisSamp = thisEvent->Event->Board[brd].firstSample;
          //EnergySample* energySample = (EnergySample*) reco.pCTEvent_->energyBoard[brd].samples->At(0);
          //if (thisSamp != 0)
          if (reco.pCTEvent_->energyBoard[brd].samples->GetLast()+1 > 0)
          {
            EnergySample* energySample0 = (EnergySample*) reco.pCTEvent_->energyBoard[brd].samples->At(0);
            sample[0][0] = energySample0->pulse[0];   // assign the first sample of the RecoEvent::sample for the first energy board
            sample[1][0] = energySample0->pulse[1];
            sample[2][0] = energySample0->pulse[2];
            int ped0 = energySample0->pulse[0]; ped[0] = ped0;
            int ped1 = energySample0->pulse[1]; ped[1] = ped1;
            int ped2 = energySample0->pulse[2]; ped[2] = ped2;
            //while (thisSamp != 0)
            for (int isample=1; isample<reco.pCTEvent_->energyBoard[brd].samples->GetLast()+1; ++isample)
            {
              EnergySample* energySample = (EnergySample*) reco.pCTEvent_->energyBoard[brd].samples->At(isample);
              sample[0][isample] = energySample->pulse[0];  // assign the rest of 16 samples of the RecoEvent::sample for the first energy board
              sample[1][isample] = energySample->pulse[1];
              sample[2][isample] = energySample->pulse[2];
              a[0] = a[0] + energySample->pulse[0] - ped0;
              a[1] = a[1] + energySample->pulse[1] - ped1;
              a[2] = a[2] + energySample->pulse[2] - ped2;
            }
          }
        }
      } else {
        a[0] = 0; ped[0] = 0;
        a[1] = 0; ped[1] = 0;
        a[2] = 0; ped[2] = 0;
      }
      brd=1;
      if (reco.pCTEvent_->energyBoard[brd].numChan>0 || reco.pCTEvent_->energyBoard[brd].numSamples>0) {
        /////   //			if (thisEvent->Event->Board[brd].enrgTag != evtNum % 4) cout << "tag mismatch, energy tag=" << thisEvent->Event->Board[brd].enrgTag << "\n";
        if (reco.pCTEvent_->energyBoard[brd].reduced) {
          a[3] = reco.pCTEvent_->energyBoard[brd].pulse[0];  ped[3] = reco.pCTEvent_->energyBoard[brd].pedestal[0];
          a[4] = reco.pCTEvent_->energyBoard[brd].pulse[1];	ped[4] = reco.pCTEvent_->energyBoard[brd].pedestal[1];
          //--no such channel-- PhCh5 =	reco.pCTEvent_->energyBoard[brd].pulse[2];
        } else {
          a[3] = 0.; a[4] = 0.; //--no such channel-- PhCh5 = 0.;
          int ped3 = 0; int ped4 = 0; int ped5 = 0;
          /// enrgSamp *thisSamp = thisEvent->Event->Board[brd].firstSample;
          //EnergySample* energySample = (EnergySample*) reco.pCTEvent_->energyBoard[brd].samples->At(0);
          //if (thisSamp != 0)
          if (reco.pCTEvent_->energyBoard[brd].samples->GetLast()+1 > 0)
          {
            EnergySample* energySample0 = (EnergySample*) reco.pCTEvent_->energyBoard[brd].samples->At(0);
            sample[3][0] = energySample0->pulse[0];      // assign the first sample of the RecoEvent::sample for the second energy board
            sample[4][0] = energySample0->pulse[1];
            ped3 = energySample0->pulse[0]; ped[3] = ped3;
            ped4 = energySample0->pulse[1]; ped[4] = ped4;
            ped5 = energySample0->pulse[2];
            //while (thisSamp != 0)
            for (int isample=1; isample<reco.pCTEvent_->energyBoard[brd].samples->GetLast()+1; ++isample)
            {
              EnergySample* energySample = (EnergySample*) reco.pCTEvent_->energyBoard[brd].samples->At(isample);
              sample[3][isample] = energySample->pulse[0];  // assign the rest of 16 samples of the RecoEvent::sample for the second energy board
              sample[4][isample] = energySample->pulse[1];
              a[3] = a[3] + energySample->pulse[0] - ped3;
              a[4] = a[4] + energySample->pulse[1] - ped4;
              //--no such channel-- PhCh5 = PhCh5 + energySample->pulse[2] - ped5;
            }
          }
        }
      } else {
        a[3] = 0; ped[3] = 0;
        a[4] = 0; ped[4] = 0;
        //--no such channel-- PhCh5 = 0;
      }
      /////////////////////////////
   }

   ClassDef(RecoEvent, 5);
};

#ifdef __MAKECINT__
#pragma link C++ class RecoEvent;
#endif

#endif  // Reco_h
