// Andriy Zatserklyaniy, April 17, 2014

#ifndef Track_h
#define Track_h

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include <TEfficiency.h>
#include <TClonesArray.h>
#include <TPaletteAxis.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <map>
#include <iomanip>
#include <cstdio>
#include <algorithm>

using std::cout;     using std::endl;

//-------------------------------- CRay.h begin ----------------------------------

class Hit2D: public TObject {
public:
   Double_t u_;
   Double_t pos_;
   Hit2D(): TObject(), u_(0), pos_(0) {}
   Hit2D(Double_t u, Double_t pos): TObject(), u_(u), pos_(pos) {}
   Hit2D(const Hit2D& hit): TObject(hit), u_(hit.u_), pos_(hit.pos_) {}

   ClassDef(Hit2D, 1);
};

#ifdef __MAKECINT__
#pragma link C++ class Hit2D;
#endif

class CRay2D: public TObject {
friend std::ostream& operator<<(std::ostream&, const CRay2D&);
public:
   Double_t x_;               // 2D radiant: u_ = 0
   Double_t cx_, cu_;	        // direction cosines
public:
   void clear() {
      x_ = 0;
      cx_ = cu_ = 0;
   }
   CRay2D(): TObject() {
      clear();
   }
   CRay2D(const CRay2D& ray2d): TObject(ray2d) {
      x_ = ray2d.x_;
      cx_ = ray2d.cx_;
      cu_ = ray2d.cu_;
   }
   CRay2D(Double_t x1, Double_t u1, Double_t x2, Double_t u2): TObject() {
      Double_t dx = x2 - x1;
      Double_t du = u2 - u1;
      Double_t alpha = TMath::ATan2(dx, du);
      cx_ = TMath::Sin(alpha);                         // cos(pi/2 - alpha) = sin(alpha)
      cu_ = TMath::Cos(alpha);
      // propagate to the plane u = 0
      Double_t p = -u1 / cu_;               // (0 - hit1.u_)/cu_
      x_ = x1 + p*cx_;
   }
   CRay2D(const Hit2D* hit1, const Hit2D* hit2): TObject() {
      Double_t dx = hit2->pos_ - hit1->pos_;
      Double_t du = hit2->u_ - hit1->u_;
      Double_t alpha = TMath::ATan2(dx, du);
      cx_ = TMath::Sin(alpha);                         // cos(pi/2 - alpha) = sin(alpha)
      cu_ = TMath::Cos(alpha);
      // propagate to the plane u = 0
      Double_t p = -hit1->u_ / cu_;          // (0 - hit1->u_)/cu_
      x_ = hit1->pos_ + p*cx_;
   }
   Double_t Theta() const {return TMath::ACos(cu_);}
   static Double_t Angle(const CRay2D* track1, const CRay2D* track2) {     // standalone function
      Double_t scalar = track1->cx_*track2->cx_ + track1->cu_*track2->cu_;
      Double_t theta = TMath::Pi()/2.;
      if (TMath::Abs(scalar) < 1.) theta = TMath::ACos(scalar);
      return theta;
   }
   Double_t at(Double_t uplane) const {
      const Double_t eps = 1e-7;
      Double_t x_proj = 1./eps;
      if (TMath::Abs(cu_) > eps) {
      	 Double_t p = uplane/cu_;        // (uplane - 0)/cu_
      	 x_proj = x_ + cx_*p;
      }
      // cout<< "CRay2D::at: x_proj = " << x_proj <<endl;
      return x_proj;
   }

   ClassDef(CRay2D, 1);
};

#ifdef __MAKECINT__
#pragma link C++ class CRay2D;
#endif

std::ostream& operator << (std::ostream& os, const CRay2D& ray) {
   os << "x_ = " << ray.x_ << " cx_ = " << ray.cx_ << " cu_ = " << ray.cu_;
   return os;
}

//-------------------------------- CRay.h end ----------------------------------

class SensorHit: public Hit2D {
   friend std::ostream& operator<<(std::ostream&, const SensorHit&);
public:
   Int_t sensorId_;
   Int_t nfirst_;                         // the first strip in the cluster (0..383)
   Int_t nstrips_;                        // the number of strips in the cluster
public:
   SensorHit(): Hit2D()
      , sensorId_(-1)
      , nfirst_(-1)
      , nstrips_(0)
   {}
   SensorHit(Int_t sensorId, Int_t nfirst, Int_t nstrips, Double_t u, Double_t pos): Hit2D(u,pos)
      , sensorId_(sensorId)
      , nfirst_(nfirst)
      , nstrips_(nstrips)
   {}
   SensorHit(const SensorHit& hit): Hit2D(hit)
                                    , sensorId_(hit.sensorId_)
                                    , nfirst_(hit.nfirst_)
                                    , nstrips_(hit.nstrips_)
   {}
   ClassDef(SensorHit, 5);
};

std::ostream& operator << (std::ostream& os, const SensorHit& hit) {
   return os
     << "sensorId " << std::setw(4) << hit.sensorId_
     << " nfirst " << std::setw(2) << hit.nfirst_
     << " nstrips " << hit.nstrips_
     << " u " << hit.u_
     << " pos " << hit.pos_;
}

#ifdef __MAKECINT__
#pragma link C++ class SensorHit;
#endif

class Track2D: public CRay2D {
friend std::ostream& operator<<(std::ostream&, const Track2D&);
  // CRay2D with 2 hits
public:
   const SensorHit* hit1_;
   const SensorHit* hit2_;
public:
   Track2D(): CRay2D() {}
   Track2D(const Track2D& track2D): CRay2D(track2D) {
      //-- hit1_ = new SensorHit(*track2D.hit1_);
      //-- hit2_ = new SensorHit(*track2D.hit2_);
      // assign hit pointers from the track2D provided
      hit1_ = track2D.hit1_;
      hit2_ = track2D.hit2_;
   }
   Track2D(const SensorHit* hit1, const SensorHit* hit2): CRay2D(hit1,hit2) {
      hit1_ = hit1;
      hit2_ = hit2;
   }
   Track2D(Double_t x1, Double_t u1, Double_t x2, Double_t u2): CRay2D(x1,u1, x2,u2) {}
   ~Track2D() {
      // do NOT delete the hits
      //delete hit1_;
      //delete hit2_;
   }
   bool hitSensor(Int_t sensorId) const {
      return hit1_->sensorId_ == sensorId || hit2_->sensorId_ == sensorId;
   }

   ClassDef(Track2D, 3)
};

#ifdef __MAKECINT__
#pragma link C++ class Track2D;
#endif

std::ostream& operator << (std::ostream& os, const Track2D& track2D) {
   os << "x_ = " << track2D.x_ << " cx_ = " << track2D.cx_ << " cu_ = " << track2D.cu_ << " hit0 " << *track2D.hit1_ << " hit1 " << *track2D.hit2_;
   return os;
}

class SuperTrack2D: public TObject {
public:
   const Track2D* itrack2D_;
   const Track2D* otrack2D_;
public:
   SuperTrack2D(): TObject(), itrack2D_(0), otrack2D_(0) {}
   SuperTrack2D(const SuperTrack2D& superTrack2D): TObject(superTrack2D) {
      itrack2D_ = new Track2D(*superTrack2D.itrack2D_);
      otrack2D_ = new Track2D(*superTrack2D.otrack2D_);
   }
   SuperTrack2D(const Track2D* itrack2D, const Track2D* otrack2D): TObject(), itrack2D_(itrack2D), otrack2D_(otrack2D) {}
   ~SuperTrack2D() {
      //delete itrack2D_;
      //delete otrack2D_;
   }
   Double_t Angle() const {
      Double_t scalar = itrack2D_->cx_*otrack2D_->cx_ + itrack2D_->cu_*otrack2D_->cu_;
      Double_t theta = TMath::Pi()/2.;
      if (TMath::Abs(scalar) < 1.) theta = TMath::ACos(scalar);
      return theta;
   }
   Double_t Distance() const {
      // Distance between the hits in the plane u = 0. Just distance between the radiants.
      return TMath::Abs(itrack2D_->x_ - otrack2D_->x_);
   }
   Double_t at(Double_t u) const {
      if (u > 0) return otrack2D_->at(u);
      else return itrack2D_->at(u);
   }
   bool iHitSensor(Int_t sensorId) const {
      return itrack2D_->hitSensor(sensorId);
   }
   bool oHitSensor(Int_t sensorId) const {
      return otrack2D_->hitSensor(sensorId);
   }
   // bool HitSensor(Int_t sensorId) const {                      // NB: do I really need this time-consuming method?
   //    return iHitSensor(sensorId) || oHitSensor(sensorId);
   // }
   Bool_t SharedHits(const SuperTrack2D* superTrack2D) const {
      if (itrack2D_->hit1_ == superTrack2D->itrack2D_->hit1_) return kTRUE;
      if (itrack2D_->hit2_ == superTrack2D->itrack2D_->hit2_) return kTRUE;
      if (itrack2D_->hit1_ == superTrack2D->itrack2D_->hit1_) return kTRUE;
      if (itrack2D_->hit2_ == superTrack2D->itrack2D_->hit2_) return kTRUE;
      if (otrack2D_->hit1_ == superTrack2D->otrack2D_->hit1_) return kTRUE;
      if (otrack2D_->hit2_ == superTrack2D->otrack2D_->hit2_) return kTRUE;
      if (otrack2D_->hit1_ == superTrack2D->otrack2D_->hit1_) return kTRUE;
      if (otrack2D_->hit2_ == superTrack2D->otrack2D_->hit2_) return kTRUE;
      return kFALSE;
   }
   bool hitSensor(Int_t sensorId) const {
      //return itrack2D_->hitSensor(sensorId) || itrack2D_->hitSensor(sensorId);
      return iHitSensor(sensorId) || oHitSensor(sensorId);
   }

   ClassDef(SuperTrack2D, 1);
};

#ifdef __MAKECINT__
#pragma link C++ class SuperTrack2D;
#endif

class SuperTrack: public TObject {
public:
   const SuperTrack2D* vTrack_;
   const SuperTrack2D* tTrack_;
   Double_t angle;      // angle between the input and output 3D tracks
   Double_t itheta;
   Double_t otheta;
   // handy vars
   Double_t vcal;       // v at the calorimeter front surface
   Double_t tcal;       // t at the calorimeter front surface
public:
   SuperTrack(): TObject(), vTrack_(0), tTrack_(0), angle(0), itheta(0), otheta(0), vcal(0), tcal(0) {}
   SuperTrack(const SuperTrack& superTrack): TObject(superTrack) {
      vTrack_ = new SuperTrack2D(*superTrack.vTrack_);
      tTrack_ = new SuperTrack2D(*superTrack.tTrack_);
      angle = superTrack.angle;
      itheta = superTrack.itheta;
      otheta = superTrack.otheta;
      at(266.9, vcal, tcal);           // 266.9 is approx u-coordinate of the front surface of the calorimeter
   }
   SuperTrack(const SuperTrack2D* vTrack, const SuperTrack2D* tTrack): TObject(), vTrack_(vTrack), tTrack_(tTrack) {
      angle = Angle();
   }
   ~SuperTrack() {
      //delete vTrack_;
      //delete tTrack_;
   }
   void at(Double_t u, Double_t& v, Double_t& t) const {
      v = vTrack_->at(u);
      t = tTrack_->at(u);
   }
   Double_t V(Double_t u=266.9) const {
     return vTrack_->at(u); 
   }
   Double_t T(Double_t u=266.9) const {
     return tTrack_->at(u); 
   }
   Double_t Angle() {
      Double_t v_xu = 0;         // ratio of the 2D directional cosines (will be used later)
      Double_t t_xu = 0;

      // input ray
      v_xu = vTrack_->itrack2D_->cx_/vTrack_->itrack2D_->cu_;
      t_xu = tTrack_->itrack2D_->cx_/tTrack_->itrack2D_->cu_;
      // 3D directional cosine for u
      Double_t cu_i = 1./TMath::Sqrt(1. + (v_xu*v_xu + t_xu*t_xu));
      itheta = 0;
      if (cu_i < 1.) itheta = TMath::ACos(cu_i);                      // theta of the input track
      // 3D directional cosines for v and t
      Double_t cv_i = v_xu * cu_i;
      Double_t ct_i = t_xu * cu_i;

      // output ray
      v_xu = vTrack_->otrack2D_->cx_/vTrack_->otrack2D_->cu_;
      t_xu = tTrack_->otrack2D_->cx_/tTrack_->otrack2D_->cu_;
      // 3D directional cosine for u
      Double_t cu_o = 1./TMath::Sqrt(1. + (v_xu*v_xu + t_xu*t_xu));
      otheta = 0;
      if (cu_o < 1.) otheta = TMath::ACos(cu_o);                      // theta of the output track
      // 3D directional cosines for v and t
      Double_t cv_o = v_xu * cu_o;
      Double_t ct_o = t_xu * cu_o;

      Double_t scalar = cv_i*cv_o + ct_i*ct_o + cu_i*cu_o;
      Double_t theta = TMath::Pi()/2.;
      if (TMath::Abs(scalar) < 1.) theta = TMath::ACos(scalar);
      return theta;
   }
   bool hitSensor(Int_t sensorId) const {
      return vTrack_->hitSensor(sensorId) || tTrack_->hitSensor(sensorId);
   }

   ClassDef(SuperTrack, 4);
};

#ifdef __MAKECINT__
#pragma link C++ class SuperTrack;
#endif

#endif	// #ifdef Track_h
