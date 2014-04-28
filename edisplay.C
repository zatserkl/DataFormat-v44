// Andriy Zatserklyaniy, April 17, 2014

#include "DataFormat.h"
#include "StripHit.h"
#include "stripPosition.h"

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
#include <map>
#include <vector>
#include <fstream>
#include <ctime>

using std::cout;     using std::endl;

void occupancy(Int_t event1=0, Int_t event2=-1, TTree* tree=0)
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

   PCTEvent* pCTEvent = 0;
   tree->SetBranchAddress("event", &pCTEvent);

   GeoHit geoHit;

   TH1I* hFPGA[12];
   for (int iFPGA=0; iFPGA<12; ++iFPGA) {
      hFPGA[iFPGA] = new TH1I(Form("hFPGA#%d",iFPGA), Form("FPGA %d",iFPGA), 2*384, 0, 2*384);
   }
   
   // TH2F* h2vt = new TH2F("h2vt", "v vs t for the first cassette", 1500,0,1500, 400,0,400);

   TH1F* heffV[4];
   for (int ilayer=0; ilayer<4; ++ilayer) heffV[ilayer] = new TH1F(Form("heffV%d",ilayer), Form("missing hit in V layer %d",ilayer), 384, 0, 384);
   TH1F* heff4V[4];
   for (int ilayer=0; ilayer<4; ++ilayer) heff4V[ilayer] = new TH1F(Form("heff4V%d",ilayer), Form("4 hits for V layer %d",ilayer), 384, 0, 384);
   //TH1F* heffratioV[4];
   //for (int ilayer=0; ilayer<4; ++ilayer) heffratioV[ilayer] = new TH1F(Form("heffratioV%d",ilayer), Form("efficiency ratio for V layer %d",ilayer), 384, 0, 384);
   TGraph* g3hits = new TGraph(3);
   g3hits->SetMarkerStyle(20);
   Int_t nhitsV = 0;
   Int_t nhitT = 0;
   std::vector<Int_t> hitV[4];
   std::vector<Int_t> hitT[4];

   Double_t ut[4] = {-211.80, -161.80, 161.80, 211.80};
   Double_t uv[4] = {-217.70, -167.6,  167.6,  217.70};

   if (event2 < event1) event2 = tree->GetEntries() - 1;

   for (int ientry=event1; ientry<=event2; ientry++)
   {
      if (tree->LoadTree(ientry) < 0) {
         cout<< "Could not load event " << ientry <<endl;
         break;
      }
      tree->GetEntry(ientry);
      if (ientry < 10 || (ientry < 10000 && ientry%1000 == 0) || ientry%100000 == 0) cout<< "processing entry " << ientry <<endl;

      geoHit.clear();

      for (int ilayer=0; ilayer<4; ++ilayer) {
         hitV[ilayer].clear();
         hitT[ilayer].clear();
      }

      for (int iFPGA=0; iFPGA<12; ++iFPGA)
      {
         const TrackerFPGA& trackerFPGA = pCTEvent->trackerFPGA[iFPGA];
         for (unsigned int ichip=0; ichip<trackerFPGA.numChips; ++ichip)
         {
            TrackerChip* trackerChip = (TrackerChip*) trackerFPGA.chips->At(ichip);
            for (int icluster=0; icluster<trackerChip->numClust; ++icluster) {
               Int_t strip = 64*trackerChip->address + (63-trackerChip->nfirst[icluster]);
               //if (iFPGA == 4 || iFPGA == 5) cout<< "strip (FPGA 4 or 5) = " << strip <<endl;
               switch (iFPGA) {
                  case 0:  hitV[0].push_back(strip); break;
                  case 1:  hitV[1].push_back(strip); break;
                  case 2:  hitV[2].push_back(strip); break;
                  case 3:  hitV[3].push_back(strip); break;
                  case 4:  hitT[0].push_back(strip); break;
                  case 5:  hitT[0].push_back(strip+768); break;
                  case 6:  hitT[1].push_back(strip); break;
                  case 7:  hitT[1].push_back(strip+768); break;
                  case 8:  hitT[2].push_back(strip); break;
                  case 9:  hitT[2].push_back(strip+768); break;
                  case 10: hitT[3].push_back(strip); break;
                  case 11: hitT[3].push_back(strip+768); break;
               }
               // cout<< ientry << "\t iFPGA = " << std::setw(2) << iFPGA << " trackerChip->address = " << std::setw(2) << (unsigned) trackerChip->address << " nstrips[" << icluster << "] = " << (unsigned) trackerChip->nstrips[icluster] << " nfirst[" << icluster << "] = " << (unsigned) trackerChip->nfirst[icluster] << " strip = " << strip <<endl;
               geoHit.stripFPGA[iFPGA][strip] = trackerChip->nstrips[icluster];

               for (int icurr=0; icurr<trackerChip->nstrips[icluster]; ++icurr) {
                  assert(strip - icurr >= 0);
                  hFPGA[iFPGA]->Fill(strip-icurr);
               }
            }
         }
      }

      // efficiency
      //for (unsigned ilayer=0; ilayer<4; ++ilayer) {
      //   cout<< ilayer << "\t hitV[ilayer].size() = " << hitV[ilayer].size() << " hitT[ilayer].size() = " << hitT[ilayer].size() <<endl;
      //}
      if (hitV[1].size() == 1 && hitV[2].size() == 1 && hitV[3].size() == 1
          && hitT[0].size() == 1 && hitT[1].size() == 1 && hitT[2].size() == 1 && hitT[3].size() == 1) {
         //-- if (hitT[0][0] >= 768) continue;
         //-- if (hitT[0][0] < 768) continue;
         ++nhitsV;
         // find positions
         for (int ilayer=1; ilayer<4; ++ilayer) {
            Double_t vstrip;
            vstrip = (hitV[ilayer][0] < 384)? hitV[ilayer][0]: 768 - hitV[ilayer][0];
            g3hits->SetPoint(ilayer-1, uv[ilayer], vstrip);
         }
         g3hits->Fit("pol1", "Q", "goff");
         // new TCanvas;
         // g3hits->Draw("ap");
         // g3hits->Fit("pol1");
         TF1* pol1 = g3hits->GetFunction("pol1");
         Double_t vfit = pol1->Eval(uv[0]);
         //cout<< "vfit = " << vfit <<endl;
         Double_t vstrip0 = (hitV[0][0] < 384)? hitV[0][0]: 768 - hitV[0][0];
         if (hitV[0].size() == 0) {
            heffV[0]->Fill(vfit);
         }
         else {
            heff4V[0]->Fill(vstrip0);
         }
      }
   }  // loop over events

   for (int iFPGA=0; iFPGA<12; ++iFPGA) {
      new TCanvas;
      hFPGA[iFPGA]->Draw();
   }

   cout<< "nhitsV = " << nhitsV <<endl;

   new TCanvas;
   heffV[0]->Draw();

   new TCanvas;
   heff4V[0]->Draw();

   TH1F* heffratioV = (TH1F*) heff4V[0]->Clone();
   heffratioV->SetNameTitle("heffratioV", "Efficiency for layer 0");
   heffratioV->Add(heffV[0], -1.);
   heffratioV->Divide(heff4V[0]);
   new TCanvas;
   heffratioV->Draw();

   tree->ResetBranchAddresses();
}

TTree* getTree()
{
   TTree* t = (TTree*) gDirectory->Get("t");
   return t;
}

class GeoSensor {
public:
   // alignment
   struct Geo {
      Double_t halfsize;
      Double_t tpin[2];
      Double_t firstStrip[4];    // position of the first strip of the t-sensors relative the pin
      Geo() {
         halfsize = 43.662;
         tpin[0] = 215.24;
         tpin[1] = 211.25;
         firstStrip[0] = 38.58;
         firstStrip[1] = 126.85;
         firstStrip[2] = 215.11;
         firstStrip[3] = 303.37;
      }
   } geo;
   // hits
   typedef std::map<Int_t,Int_t> SensorMap;
   SensorMap sensorMap[24];
   std::vector<Double_t> tHits[4];
   std::vector<Double_t> vHitsR[4];
   std::vector<Double_t> vHitsL[4];
};

class SSDSensorGeometryT {
public:
   Double_t pin;        // coordinate of the pin position for layers 0,1: t, for layers 2,3: -t
   Double_t dstrip1;    // distance to pin
   Int_t* strips;       // pointer to the array of strips
};

Int_t edisplay(Int_t event, TTree* tree=0, const char* wname="event_display", EventOutput* eventOutput=0)
{
   if (!tree) tree = (TTree*) gDirectory->Get("t");
   if (!tree) {
      cout<< "Could not find tree t" <<endl;
      return -1;
   }

   // constants
   const Double_t pitch = 0.228;
   const Double_t tgap = 0.702;                          // dead space between the sensors, Hartmut data
   //-- const Double_t vgap = 0.542;                          // dead space between the sensors, Hartmut data

   // const Double_t toffset[4] = {0, 0, 0, 0};
   const Double_t toffset[4] = {0.0, 4.0, 4.0, 0.0};
   // const Double_t toffset[4] = {-2.0, -4.0, 4.0, 2.0};
   Double_t ut[4] = {-211.80, -161.80, 161.80, 211.80};
   Double_t uv[4] = {-217.70, -167.6,  167.6,  217.70};

   // const Double_t voffset = 0;

   TCanvas* can = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(wname);
   if (can) gPad = can;
   else can = new TCanvas("event_display", wname, 700,500);

   ///   RunHeader* runHeader = (RunHeader*) tree->GetUserInfo()->First();
   ///   cout<< "runHeader->GetRun() = " << runHeader->GetRun() <<endl;
   ///   time_t start_time = runHeader->GetTime();
   ///   cout<< "run start time: " << std::ctime(&start_time);
   ///   cout<< "program version is " << runHeader->GetVersion() <<endl;
   ///   if (runHeader->GetTimeTag()) cout<< "event time tag was written out" <<endl;
   ///   else cout<< "event time tag was not written out" <<endl;
   ///   cout<<endl;

   PCTEvent* pCTEvent = 0;
   tree->SetBranchAddress("event", &pCTEvent);

   GeoHit geoHit;

   if (tree->LoadTree(event) < 0) {
      cout<< "Could not load event " << event <<endl;
      return -1;
   }
   tree->GetEntry(event);

   geoHit.clear();

   typedef std::map<Int_t,Int_t> SensorMap;
   SensorMap sensor[24];                     // 24 sensors total in all 12 FPGAs

   Int_t nhits_total = 0;

   for (int iFPGA=0; iFPGA<12; ++iFPGA)
   {
      const TrackerFPGA& trackerFPGA = pCTEvent->trackerFPGA[iFPGA];
      for (unsigned int ichip=0; ichip<trackerFPGA.numChips; ++ichip)
      {
         TrackerChip* trackerChip = (TrackerChip*) trackerFPGA.chips->At(ichip);
         for (int icluster=0; icluster<trackerChip->numClust; ++icluster) {
            Int_t strip = 64*trackerChip->address + (63-trackerChip->nfirst[icluster]);
            cout<< event << "\t iFPGA = " << std::setw(2) << iFPGA << " trackerChip->address = " << std::setw(2) << (unsigned) trackerChip->address << " nstrips[" << icluster << "] = " << (unsigned) trackerChip->nstrips[icluster] << " nfirst[" << icluster << "] = " << (unsigned) trackerChip->nfirst[icluster] << " strip = " << strip <<endl;
            //--readout-- geoHit.stripFPGA[iFPGA][strip] = trackerChip->nstrips[icluster];
            geoHit.stripFPGA[iFPGA][strip - trackerChip->nstrips[icluster] + 1] = trackerChip->nstrips[icluster];
            ++nhits_total;

            // sensor map
            Int_t nsensor = iFPGA*2 + strip/384;
            sensor[nsensor][strip - trackerChip->nstrips[icluster] + 1] = trackerChip->nstrips[icluster];
         }
      }
   }

   cout<< "sensor maps" <<endl;
   for (int isensor=0; isensor<24; ++isensor) {
      SensorMap& sensorMap = sensor[isensor];
      for (std::map<Int_t,Int_t>::const_iterator it=sensorMap.begin(); it!=sensorMap.end(); ++it) {
         cout<< "sensor" << std::setw(3) << isensor << " it->first = " << it->first << " it->second = " << it->second <<endl;
      }
   }

   // combine clusters
   SensorMap sensorClusters[24];                      // for all 24 sensors in all 12 FPGAs
   for (int isensor=0; isensor<24; ++isensor) {
      SensorMap& sensorMap = sensor[isensor];
      std::map<Int_t,Int_t>::const_iterator it=sensorMap.begin();
      while (it!=sensorMap.end()) {
         Int_t strip = it->first;
         Int_t nstrips = it->second;
         cout<< "sensor" << std::setw(3) << isensor << " it->first = " << it->first << " it->second = " << it->second <<endl;
         sensorClusters[isensor][strip] = nstrips;
         // look at the next cluster
         ++it;
         if (it != sensorMap.end()) {
            if (it->first == strip + nstrips) {
               sensorClusters[isensor][strip] += it->second;
               ++it;
            }
         }
      }
   }

   cout<< "sensorClusters maps" <<endl;
   for (int isensor=0; isensor<24; ++isensor) {
      SensorMap& sensorMap = sensorClusters[isensor];
      for (std::map<Int_t,Int_t>::const_iterator it=sensorMap.begin(); it!=sensorMap.end(); ++it) {
         cout<< "sensor" << std::setw(3) << isensor << " it->first = " << it->first << " it->second = " << it->second <<endl;
      }
   }

   if (nhits_total == 0) return -1;

   // Display the event

   can->DrawFrame(-300,-250, 300,250, Form("event %d   tgap=%0.1f toffset: %0.2f %0.2f %0.2f %0.2f",event,tgap,toffset[0],toffset[1],toffset[2],toffset[3]));
   TLine* line = new TLine();

   // forward
   line->DrawLine(-217.70,-220.11, -217.70,215.24);
   line->DrawLine(-211.80,-220.11, -211.80,215.24);
   line->DrawLine(-167.60,-224.11, -167.60,211.24);
   line->DrawLine(-161.80,-224.11, -161.80,211.24);
   // back
   line->DrawLine( 217.70,-220.11,  217.70,215.24);
   line->DrawLine( 211.80,-220.11,  211.80,215.24);
   line->DrawLine( 167.60,-224.11,  167.60,211.24);
   line->DrawLine( 161.80,-224.11,  161.80,211.24);

   TMarker* tmarker = new TMarker();
   tmarker->SetMarkerStyle(24);
   tmarker->SetMarkerColor(2);
   TLine* tline = new TLine();
   tline->SetLineColor(2);

   TMarker* vmarker = new TMarker();
   vmarker->SetMarkerStyle(24);
   vmarker->SetMarkerColor(4);

   TMarker* v0marker = new TMarker();
   v0marker->SetMarkerStyle(24);
   v0marker->SetMarkerColor(4);
   // TLine* v0line = new TLine();
   // v0line->SetLineColor(4);

   TMarker* v1marker = new TMarker();
   v1marker->SetMarkerStyle(24);
   v1marker->SetMarkerColor(8);
   // TLine* v1line = new TLine();
   // v1line->SetLineColor(8);

   TLine* vline = new TLine();
   vline->SetLineColor(1);

   Double_t marker_ut[4], marker_t[4];
   for (int ilayer=0; ilayer<4; ++ilayer) marker_ut[ilayer] = 0;     // hit will have non-zero coordinate
   Double_t marker_uv0[4], marker_v0[4];
   for (int ilayer=0; ilayer<4; ++ilayer) marker_uv0[ilayer] = 0;    // hit will have non-zero coordinate
   Double_t marker_uv1[4], marker_v1[4];
   for (int ilayer=0; ilayer<4; ++ilayer) marker_uv1[ilayer] = 0;    // hit will have non-zero coordinate

   Int_t nhits_t[4];    // the number of hits in the layer
   Int_t nhits_v[4];
   for (int ilayer=0; ilayer<4; ++ilayer) {
      nhits_t[ilayer] = 0;
      nhits_v[ilayer] = 0;
   }

	Double_t tpin[4] = {215.24, 211.25, -211.25, -215.24};  //T coordinate of alignment pin per layer
	Double_t tfirst[4] = {38.58, 126.85, 215.11, 303.37};   //These numbers actually come from board T4

   for (int ilayer=0; ilayer<4; ++ilayer) {
      // t-layer
      for (int istrip=0; istrip<4*384; ++istrip) {
         Int_t nstrips = geoHit.tstrip[ilayer][istrip];
         if (nstrips > 0) {
            ++nhits_t[ilayer];
            for (int icurr=0; icurr<nstrips; ++icurr) {
               Int_t strip = istrip + icurr;
               assert(strip >= 0);
               Int_t ngaps = strip/384;
               Double_t treadout = strip*pitch + ngaps*tgap - (768*pitch + 1.5*tgap);
               //Double_t t = ilayer < 2? -treadout: treadout;
               Double_t t = treadout + toffset[ilayer];
               if (ilayer < 2) t = -t;
               tmarker->DrawMarker(ut[ilayer], t);
               marker_ut[ilayer] = ut[ilayer];
               marker_t[ilayer] = t;
               cout<< "t  layer: ut[" << ilayer << "] = " << ut[ilayer] << " istrip = " << istrip << " strip = " << strip << " t = " << t <<endl;
            }
         }
      }

      // v-layer right: lower chip addresses for layers 0,1 higher chip addresses for layers 2,3
      // v-layer left:  higher chip addresses for layers 0,1 lower chip addresses for layers 2,3

      // v-layer lower chip addresses
      for (int istrip=0; istrip<384; ++istrip) {
         Int_t nstrips = geoHit.vstrip[ilayer][istrip];
         if (nstrips > 0) {
            ++nhits_v[ilayer];
            // cout<< "layer " << ilayer << " incremented nhits_v 0..383" <<endl;
            for (int icurr=0; icurr<nstrips; ++icurr) {
               Int_t strip = istrip + icurr;
               assert(strip >= 0);
               Double_t v = strip*pitch - (384/2)*pitch;
               v *= -1;                      // invert direction for chips 0..5

               if (ilayer < 2) {
                  marker_uv0[ilayer] = uv[ilayer];
                  marker_v0[ilayer] = v;
                  vmarker->SetMarkerColor(4);                  // right (negative t) part: blue
                  vmarker->DrawMarker(uv[ilayer], v);
                  cout<< "v0 layer: uv[" << ilayer << "] = " << uv[ilayer] << " istrip = " << istrip << " strip = " << strip << " v = " << v <<endl;
               }
               else {
                  marker_uv1[ilayer] = uv[ilayer];
                  marker_v1[ilayer] = v;
                  vmarker->SetMarkerColor(8);                  // left (positive t) part: green
                  vmarker->DrawMarker(uv[ilayer], v);
                  cout<< "v1 layer: uv[" << ilayer << "] = " << uv[ilayer] << " istrip = " << istrip << " strip = " << strip << " v = " << v <<endl;
               }
            }
         }
      }
      // v-layer higher chip addresses
      for (int istrip=384; istrip<2*384; ++istrip) {
         Int_t nstrips = geoHit.vstrip[ilayer][istrip];
         if (nstrips > 0) {
            ++nhits_v[ilayer];
            // cout<< "layer " << ilayer << " incremented nhits_v 384..767" <<endl;
            for (int icurr=0; icurr<nstrips; ++icurr) {
               Int_t strip = istrip + icurr - 384;
               assert(strip >= 0);
               //Double_t v = (strip-384)*pitch - (384/2)*pitch;
               Double_t v = strip*pitch - (384/2)*pitch;
               if (ilayer < 2) {
                  marker_uv1[ilayer] = uv[ilayer];
                  marker_v1[ilayer] = v;
                  vmarker->SetMarkerColor(8);            // left (positive t) part: green
                  vmarker->DrawMarker(uv[ilayer], v);
                  cout<< "v1 layer: uv[" << ilayer << "] = " << uv[ilayer] << " istrip = " << istrip << " strip = " << strip << " v = " << v <<endl;
               }
               else {
                  marker_uv0[ilayer] = uv[ilayer];
                  marker_v0[ilayer] = v;
                  vmarker->SetMarkerColor(4);            // right (negative t) part: blue
                  vmarker->DrawMarker(uv[ilayer], v);
                  cout<< "v0 layer: uv[" << ilayer << "] = " << uv[ilayer] << " istrip = " << istrip << " strip = " << strip << " v = " << v <<endl;
               }
            }
         }
      }
   }

   // for (int ilayer=0; ilayer<4; ++ilayer) {
   //    cout<< "nhits_t[" << ilayer << "] = " << nhits_t[ilayer] << " nhits_v[" << ilayer << "] = " << nhits_v[ilayer] <<endl;
   // }
   bool good_t = nhits_t[0] == 1 && nhits_t[1] == 1 && nhits_t[2] == 1 && nhits_t[3] == 1;
   bool good_v = nhits_v[0] == 1 && nhits_v[1] == 1 && nhits_v[2] == 1 && nhits_v[3] == 1;

   if (good_t) tline->DrawLine(marker_ut[0], marker_t[0], marker_ut[3], marker_t[3]);
   if (good_v) {
      Double_t vup = marker_uv0[0]? marker_v0[0]: marker_v1[0];
      Double_t vdown = marker_uv0[3]? marker_v0[3]: marker_v1[3];
      vline->DrawLine(uv[0], vup, uv[3], vdown);
      //cout<< "vup = " << vup << " vdown = " << vdown <<endl;
   }

   if (eventOutput)
   {
      eventOutput->clear();
      if (good_t && good_v) {
         eventOutput->good = kTRUE;
         for (int ilayer=0; ilayer<4; ++ilayer) {
            eventOutput->u[ilayer] = marker_ut[ilayer];
            eventOutput->t[ilayer] = marker_t[ilayer];
            eventOutput->v[ilayer] = marker_uv0[ilayer]? marker_v0[ilayer]: marker_v1[ilayer];
         }
         eventOutput->wepl = 10.;
      }
   }

   can->Update();

   tree->ResetBranchAddresses();

   delete tmarker;
   delete v0marker;
   delete v1marker;
   delete tline;
   //delete v0line;
   //delete v1line;
   delete vline;

   return event;
}

Int_t edisplay_input(Int_t event, TTree* tree=0, const char* wname="event_display", EventOutput* eventOutput=0)
{
   if (!tree) tree = (TTree*) gDirectory->Get("t");
   if (!tree) {
      cout<< "Could not find tree t" <<endl;
      return -1;
   }

   // constants
   const Double_t pitch = 0.228;
   const Double_t tgap = 0.702;                          // dead space between the sensors, Hartmut data
   //-- const Double_t vgap = 0.542;                          // dead space between the sensors, Hartmut data

   // const Double_t toffset[4] = {0, 0, 0, 0};
   const Double_t toffset[4] = {0.0, 4.0, 4.0, 0.0};
   // const Double_t toffset[4] = {-2.0, -4.0, 4.0, 2.0};
   Double_t ut[4] = {-211.80, -161.80, 161.80, 211.80};
   Double_t uv[4] = {-217.70, -167.6,  167.6,  217.70};

   // const Double_t voffset = 0;

   TCanvas* can = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(wname);
   if (can) gPad = can;
   else can = new TCanvas("event_display", wname, 700,500);

   ///   RunHeader* runHeader = (RunHeader*) tree->GetUserInfo()->First();
   ///   cout<< "runHeader->GetRun() = " << runHeader->GetRun() <<endl;
   ///   time_t start_time = runHeader->GetTime();
   ///   cout<< "run start time: " << std::ctime(&start_time);
   ///   cout<< "program version is " << runHeader->GetVersion() <<endl;
   ///   if (runHeader->GetTimeTag()) cout<< "event time tag was written out" <<endl;
   ///   else cout<< "event time tag was not written out" <<endl;
   ///   cout<<endl;

   PCTEvent* pCTEvent = 0;
   tree->SetBranchAddress("event", &pCTEvent);

   GeoHit geoHit;

   if (tree->LoadTree(event) < 0) {
      cout<< "Could not load event " << event <<endl;
      return -1;
   }
   tree->GetEntry(event);

   geoHit.clear();

   typedef std::map<Int_t,Int_t> SensorMap;
   SensorMap sensor[24];                     // 24 sensors total in all 12 FPGAs

   Int_t nhits_total = 0;

   for (int iFPGA=0; iFPGA<12; ++iFPGA)
   {
      const TrackerFPGA& trackerFPGA = pCTEvent->trackerFPGA[iFPGA];
      for (unsigned int ichip=0; ichip<trackerFPGA.numChips; ++ichip)
      {
         TrackerChip* trackerChip = (TrackerChip*) trackerFPGA.chips->At(ichip);
         for (int icluster=0; icluster<trackerChip->numClust; ++icluster) {
            Int_t strip = 64*trackerChip->address + (63-trackerChip->nfirst[icluster]);
            cout<< event << "\t iFPGA = " << std::setw(2) << iFPGA << " trackerChip->address = " << std::setw(2) << (unsigned) trackerChip->address << " nstrips[" << icluster << "] = " << (unsigned) trackerChip->nstrips[icluster] << " nfirst[" << icluster << "] = " << (unsigned) trackerChip->nfirst[icluster] << " strip = " << strip <<endl;
            //--readout-- geoHit.stripFPGA[iFPGA][strip] = trackerChip->nstrips[icluster];
            geoHit.stripFPGA[iFPGA][strip - trackerChip->nstrips[icluster] + 1] = trackerChip->nstrips[icluster];
            ++nhits_total;

            // sensor map
            Int_t nsensor = iFPGA*2 + strip/384;
            sensor[nsensor][strip - trackerChip->nstrips[icluster] + 1] = trackerChip->nstrips[icluster];
         }
      }
   }

   cout<< "sensor maps" <<endl;
   for (int isensor=0; isensor<24; ++isensor) {
      SensorMap& sensorMap = sensor[isensor];
      for (std::map<Int_t,Int_t>::const_iterator it=sensorMap.begin(); it!=sensorMap.end(); ++it) {
         cout<< "sensor" << std::setw(3) << isensor << " it->first = " << it->first << " it->second = " << it->second <<endl;
      }
   }

   // combine clusters
   SensorMap sensorClusters[24];                      // for all 24 sensors in all 12 FPGAs
   for (int isensor=0; isensor<24; ++isensor) {
      SensorMap& sensorMap = sensor[isensor];
      std::map<Int_t,Int_t>::const_iterator it=sensorMap.begin();
      while (it!=sensorMap.end()) {
         Int_t strip = it->first;
         Int_t nstrips = it->second;
         cout<< "sensor" << std::setw(3) << isensor << " it->first = " << it->first << " it->second = " << it->second <<endl;
         sensorClusters[isensor][strip] = nstrips;
         // look at the next cluster
         ++it;
         if (it != sensorMap.end()) {
            if (it->first == strip + nstrips) {
               sensorClusters[isensor][strip] += it->second;
               ++it;
            }
         }
      }
   }

   cout<< "sensorClusters maps" <<endl;
   for (int isensor=0; isensor<24; ++isensor) {
      SensorMap& sensorMap = sensorClusters[isensor];
      for (std::map<Int_t,Int_t>::const_iterator it=sensorMap.begin(); it!=sensorMap.end(); ++it) {
         cout<< "sensor" << std::setw(3) << isensor << " it->first = " << it->first << " it->second = " << it->second <<endl;
      }
   }

   if (nhits_total == 0) return -1;

   // Display the event

   can->DrawFrame(-300,-250, 300,250, Form("event %d   tgap=%0.1f toffset: %0.2f %0.2f %0.2f %0.2f",event,tgap,toffset[0],toffset[1],toffset[2],toffset[3]));
   TLine* line = new TLine();

   // forward
   line->DrawLine(-217.70,-220.11, -217.70,215.24);
   line->DrawLine(-211.80,-220.11, -211.80,215.24);
   line->DrawLine(-167.60,-224.11, -167.60,211.24);
   line->DrawLine(-161.80,-224.11, -161.80,211.24);
   // back
   line->DrawLine( 217.70,-220.11,  217.70,215.24);
   line->DrawLine( 211.80,-220.11,  211.80,215.24);
   line->DrawLine( 167.60,-224.11,  167.60,211.24);
   line->DrawLine( 161.80,-224.11,  161.80,211.24);

   TMarker* tmarker = new TMarker();
   tmarker->SetMarkerStyle(24);
   tmarker->SetMarkerColor(2);
   TLine* tline = new TLine();
   tline->SetLineColor(2);

   TMarker* vmarker = new TMarker();
   vmarker->SetMarkerStyle(24);
   vmarker->SetMarkerColor(4);

   TMarker* v0marker = new TMarker();
   v0marker->SetMarkerStyle(24);
   v0marker->SetMarkerColor(4);
   // TLine* v0line = new TLine();
   // v0line->SetLineColor(4);

   TMarker* v1marker = new TMarker();
   v1marker->SetMarkerStyle(24);
   v1marker->SetMarkerColor(8);
   // TLine* v1line = new TLine();
   // v1line->SetLineColor(8);

   TLine* vline = new TLine();
   vline->SetLineColor(1);

   Double_t marker_ut[4], marker_t[4];
   for (int ilayer=0; ilayer<4; ++ilayer) marker_ut[ilayer] = 0;     // hit will have non-zero coordinate
   Double_t marker_uv0[4], marker_v0[4];
   for (int ilayer=0; ilayer<4; ++ilayer) marker_uv0[ilayer] = 0;    // hit will have non-zero coordinate
   Double_t marker_uv1[4], marker_v1[4];
   for (int ilayer=0; ilayer<4; ++ilayer) marker_uv1[ilayer] = 0;    // hit will have non-zero coordinate

   Int_t nhits_t[4];    // the number of hits in the layer
   Int_t nhits_v[4];
   for (int ilayer=0; ilayer<4; ++ilayer) {
      nhits_t[ilayer] = 0;
      nhits_v[ilayer] = 0;
   }

	Double_t tpin[4] = {215.24, 211.25, -211.25, -215.24};  //T coordinate of alignment pin per layer
	Double_t tfirst[4] = {38.58, 126.85, 215.11, 303.37};   //These numbers actually come from board T4

   for (int ilayer=0; ilayer<4; ++ilayer) {
      // t-layer
      for (int istrip=0; istrip<4*384; ++istrip) {
         Int_t nstrips = geoHit.tstrip[ilayer][istrip];
         if (nstrips > 0) {
            ++nhits_t[ilayer];
            for (int icurr=0; icurr<nstrips; ++icurr) {
               Int_t strip = istrip + icurr;
               assert(strip >= 0);
               Int_t ngaps = strip/384;
               Double_t treadout = strip*pitch + ngaps*tgap - (768*pitch + 1.5*tgap);
               //Double_t t = ilayer < 2? -treadout: treadout;
               Double_t t = treadout + toffset[ilayer];
               if (ilayer < 2) t = -t;
               tmarker->DrawMarker(ut[ilayer], t);
               marker_ut[ilayer] = ut[ilayer];
               marker_t[ilayer] = t;
               cout<< "t  layer: ut[" << ilayer << "] = " << ut[ilayer] << " istrip = " << istrip << " strip = " << strip << " t = " << t <<endl;
            }
         }
      }

      // v-layer right: lower chip addresses for layers 0,1 higher chip addresses for layers 2,3
      // v-layer left:  higher chip addresses for layers 0,1 lower chip addresses for layers 2,3

      // v-layer lower chip addresses
      for (int istrip=0; istrip<384; ++istrip) {
         Int_t nstrips = geoHit.vstrip[ilayer][istrip];
         if (nstrips > 0) {
            ++nhits_v[ilayer];
            // cout<< "layer " << ilayer << " incremented nhits_v 0..383" <<endl;
            for (int icurr=0; icurr<nstrips; ++icurr) {
               Int_t strip = istrip + icurr;
               assert(strip >= 0);
               Double_t v = strip*pitch - (384/2)*pitch;
               v *= -1;                      // invert direction for chips 0..5

               if (ilayer < 2) {
                  marker_uv0[ilayer] = uv[ilayer];
                  marker_v0[ilayer] = v;
                  vmarker->SetMarkerColor(4);                  // right (negative t) part: blue
                  vmarker->DrawMarker(uv[ilayer], v);
                  cout<< "v0 layer: uv[" << ilayer << "] = " << uv[ilayer] << " istrip = " << istrip << " strip = " << strip << " v = " << v <<endl;
               }
               else {
                  marker_uv1[ilayer] = uv[ilayer];
                  marker_v1[ilayer] = v;
                  vmarker->SetMarkerColor(8);                  // left (positive t) part: green
                  vmarker->DrawMarker(uv[ilayer], v);
                  cout<< "v1 layer: uv[" << ilayer << "] = " << uv[ilayer] << " istrip = " << istrip << " strip = " << strip << " v = " << v <<endl;
               }
            }
         }
      }
      // v-layer higher chip addresses
      for (int istrip=384; istrip<2*384; ++istrip) {
         Int_t nstrips = geoHit.vstrip[ilayer][istrip];
         if (nstrips > 0) {
            ++nhits_v[ilayer];
            // cout<< "layer " << ilayer << " incremented nhits_v 384..767" <<endl;
            for (int icurr=0; icurr<nstrips; ++icurr) {
               Int_t strip = istrip + icurr - 384;
               assert(strip >= 0);
               //Double_t v = (strip-384)*pitch - (384/2)*pitch;
               Double_t v = strip*pitch - (384/2)*pitch;
               if (ilayer < 2) {
                  marker_uv1[ilayer] = uv[ilayer];
                  marker_v1[ilayer] = v;
                  vmarker->SetMarkerColor(8);            // left (positive t) part: green
                  vmarker->DrawMarker(uv[ilayer], v);
                  cout<< "v1 layer: uv[" << ilayer << "] = " << uv[ilayer] << " istrip = " << istrip << " strip = " << strip << " v = " << v <<endl;
               }
               else {
                  marker_uv0[ilayer] = uv[ilayer];
                  marker_v0[ilayer] = v;
                  vmarker->SetMarkerColor(4);            // right (negative t) part: blue
                  vmarker->DrawMarker(uv[ilayer], v);
                  cout<< "v0 layer: uv[" << ilayer << "] = " << uv[ilayer] << " istrip = " << istrip << " strip = " << strip << " v = " << v <<endl;
               }
            }
         }
      }
   }

   // for (int ilayer=0; ilayer<4; ++ilayer) {
   //    cout<< "nhits_t[" << ilayer << "] = " << nhits_t[ilayer] << " nhits_v[" << ilayer << "] = " << nhits_v[ilayer] <<endl;
   // }
   bool good_t = nhits_t[0] == 1 && nhits_t[1] == 1 && nhits_t[2] == 1 && nhits_t[3] == 1;
   bool good_v = nhits_v[0] == 1 && nhits_v[1] == 1 && nhits_v[2] == 1 && nhits_v[3] == 1;

   if (good_t) tline->DrawLine(marker_ut[0], marker_t[0], marker_ut[3], marker_t[3]);
   if (good_v) {
      Double_t vup = marker_uv0[0]? marker_v0[0]: marker_v1[0];
      Double_t vdown = marker_uv0[3]? marker_v0[3]: marker_v1[3];
      vline->DrawLine(uv[0], vup, uv[3], vdown);
      //cout<< "vup = " << vup << " vdown = " << vdown <<endl;
   }

   if (eventOutput)
   {
      eventOutput->clear();
      if (good_t && good_v) {
         eventOutput->good = kTRUE;
         for (int ilayer=0; ilayer<4; ++ilayer) {
            eventOutput->u[ilayer] = marker_ut[ilayer];
            eventOutput->t[ilayer] = marker_t[ilayer];
            eventOutput->v[ilayer] = marker_uv0[ilayer]? marker_v0[ilayer]: marker_v1[ilayer];
         }
         eventOutput->wepl = 10.;
      }
   }

   can->Update();

   tree->ResetBranchAddresses();

   delete tmarker;
   delete v0marker;
   delete v1marker;
   delete tline;
   //delete v0line;
   //delete v1line;
   delete vline;

   return event;
}

void processEvents(const char* input_filename="pCTraw_Run_51.out-0.root")
{
   TFile* input_file = new TFile(input_filename);
   if (!input_file) {
      cout<< "Could not open file " << input_filename <<endl;
      return;
   }

   TTree* tree = (TTree*) input_file->Get("t");
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

   EventOutput eventOutput;

   std::vector<Int_t> events;
   events.push_back(4);
   events.push_back(241);
   events.push_back(345);
   events.push_back(575);
   events.push_back(589);
   events.push_back(621);
   events.push_back(634);
   events.push_back(673);
   events.push_back(683);
   events.push_back(864);

   const char* tname[4] = {"temporary_file_t0file.dat", "temporary_file_t1file.dat", "temporary_file_t2file.dat", "temporary_file_t3file.dat"};
   const char* vname[4] = {"temporary_file_v0file.dat", "temporary_file_v1file.dat", "temporary_file_v2file.dat", "temporary_file_v3file.dat"};
   const char* uname[4] = {"temporary_file_u0file.dat", "temporary_file_u1file.dat", "temporary_file_u2file.dat", "temporary_file_u3file.dat"};
   const char* wname = "temporary_file_wfile.dat";

   std::fstream tfile[4];
   std::fstream vfile[4];
   std::fstream ufile[4];
   std::fstream wfile;
   for (int ilayer=0; ilayer<4; ++ilayer) {
      //cout<< "open file " << tname[ilayer] <<endl;
      tfile[ilayer].open(tname[ilayer], std::ios::binary | std::ios::out);
      tfile[ilayer].close();
      tfile[ilayer].open(tname[ilayer], std::ios::binary | std::ios::in | std::ios::out);
      //if (tfile[ilayer].is_open()) cout<< "successfully opened file " << tname[ilayer] <<endl;
      //else cout<< "Could not open file " << tname[ilayer] <<endl;

      vfile[ilayer].open(vname[ilayer], std::ios::binary | std::ios::out);
      vfile[ilayer].close();
      vfile[ilayer].open(vname[ilayer], std::ios::binary | std::ios::in | std::ios::out);

      ufile[ilayer].open(uname[ilayer], std::ios::binary | std::ios::out);
      ufile[ilayer].close();
      ufile[ilayer].open(uname[ilayer], std::ios::binary | std::ios::in | std::ios::out);
   }
   wfile.open(wname, std::ios::binary | std::ios::out);
   wfile.close();
   wfile.open(wname, std::ios::binary | std::ios::in | std::ios::out);

   Int_t ngood = 0;
   for (unsigned ievent=0; ievent<events.size(); ++ievent) {
      edisplay_input(events[ievent], tree, "event_display", &eventOutput);
      if (eventOutput.good) {
         ++ngood;
         cout<< "--> event " << events[ievent] << ": ";
         for (int ilayer=0; ilayer<4; ilayer++) cout<< eventOutput.t[ilayer] << " ";
         for (int ilayer=0; ilayer<4; ilayer++) cout<< eventOutput.v[ilayer] << " ";
         for (int ilayer=0; ilayer<4; ilayer++) cout<< eventOutput.u[ilayer] << " ";
         cout<< eventOutput.wepl <<endl;

         for (int ilayer=0; ilayer<4; ilayer++) tfile[ilayer].write((const char*) &eventOutput.t[ilayer], sizeof(Float_t));
         for (int ilayer=0; ilayer<4; ilayer++) vfile[ilayer].write((const char*) &eventOutput.v[ilayer], sizeof(Float_t));
         for (int ilayer=0; ilayer<4; ilayer++) ufile[ilayer].write((const char*) &eventOutput.u[ilayer], sizeof(Float_t));
         wfile.write((const char*) &eventOutput.wepl, sizeof(Float_t));
      }
   }
   cout<<endl<< ngood << " good events out of " << events.size() <<endl;

   cout<< "look at the output file" <<endl;

   for (int ilayer=0; ilayer<4; ilayer++) {
      tfile[ilayer].seekg(0, tfile[ilayer].beg);
      vfile[ilayer].seekg(0, vfile[ilayer].beg);
      ufile[ilayer].seekg(0, ufile[ilayer].beg);
   }
   wfile.seekg(0, wfile.beg);

   // output file
   const char* uniname = Form("%s.bin",input_filename);
   std::ofstream unifile(uniname, std::ios::binary);

   char magicNumber[4];
   magicNumber[0] = 'P';
   magicNumber[1] = 'C';
   magicNumber[2] = 'T';
   magicNumber[3] = 'D';
   unifile.write(magicNumber, 4);

   Int_t versionNumberIdentifier = 0;
   unifile.write((const char*) &versionNumberIdentifier, sizeof(Int_t));

   Int_t numberEvents = ngood;
   unifile.write((const char*) &numberEvents, sizeof(Int_t));
   
   Float_t projectionAngle = 0;
   unifile.write((const char*) &projectionAngle, sizeof(Float_t));

   Float_t beamEnergy = 200;
   unifile.write((const char*) &beamEnergy, sizeof(Float_t));

   Int_t acquisitionDate = start_time;
   unifile.write((const char*) &acquisitionDate, sizeof(Int_t));

   time_t timer = std::time(NULL);
   Int_t preprocessingDate = timer;
   unifile.write((const char*) &preprocessingDate, sizeof(Int_t));

   Int_t variableStringSize = 0;       // will be used for each of the variable length string

   std::string phantomName = "Very nice phantom";
   variableStringSize = phantomName.size() + 1;
   unifile.write((const char*) &variableStringSize, sizeof(Int_t));
   unifile.write(phantomName.c_str(), variableStringSize);

   std::string dataSource = "Data Source";
   variableStringSize = dataSource.size() + 1;
   unifile.write((const char*) &variableStringSize, sizeof(Int_t));
   unifile.write(dataSource.c_str(), variableStringSize);

   std::string preparedBy = "Tia";
   variableStringSize = preparedBy.size() + 1;
   unifile.write((const char*) &variableStringSize, sizeof(Int_t));
   unifile.write(preparedBy.c_str(), variableStringSize);

   for (int ilayer=0; ilayer<4; ++ilayer) {
      std::istreambuf_iterator<char> begin_source(tfile[ilayer]);
      std::istreambuf_iterator<char> end_source;
      std::ostreambuf_iterator<char> begin_dest(unifile);
      std::copy(begin_source, end_source, begin_dest);
   }

   for (int ilayer=0; ilayer<4; ++ilayer) {
      std::istreambuf_iterator<char> begin_source(vfile[ilayer]);
      std::istreambuf_iterator<char> end_source;
      std::ostreambuf_iterator<char> begin_dest(unifile);
      std::copy(begin_source, end_source, begin_dest);
   }

   for (int ilayer=0; ilayer<4; ++ilayer) {
      std::istreambuf_iterator<char> begin_source(ufile[ilayer]);
      std::istreambuf_iterator<char> end_source;
      std::ostreambuf_iterator<char> begin_dest(unifile);
      std::copy(begin_source, end_source, begin_dest);
   }
   std::istreambuf_iterator<char> begin_source(wfile);
   std::istreambuf_iterator<char> end_source;
   std::ostreambuf_iterator<char> begin_dest(unifile);
   std::copy(std::istreambuf_iterator<char>(wfile), std::istreambuf_iterator<char>(), std::ostreambuf_iterator<char>(unifile));

   cout<< "Wrote " << ngood << " events into output binary file " << uniname <<endl;
}

void events_tv(Int_t event1=0, Int_t event2=-1, TTree* tree=0)
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

   PCTEvent* pCTEvent = 0;
   tree->SetBranchAddress("event", &pCTEvent);

   GeoHit geoHit;

   if (event2 < event1) event2 = tree->GetEntries() - 1;

   Int_t ngood = 0;

   for (int ientry=event1; ientry<=event2; ++ientry) {
      if (tree->LoadTree(ientry) < 0) break;
      tree->GetEntry(ientry);

      if (event2-event1 < 100 || pCTEvent->evt%10000 == 0) cout<< "--> pCTEvent->evt = " << pCTEvent->evt << "\t event = " << ientry <<endl;

      geoHit.clear();

      Int_t nhits_total = 0;

      for (int iFPGA=0; iFPGA<12; ++iFPGA)
      {
         const TrackerFPGA& trackerFPGA = pCTEvent->trackerFPGA[iFPGA];
         for (unsigned int ichip=0; ichip<trackerFPGA.numChips; ++ichip)
         {
            TrackerChip* trackerChip = (TrackerChip*) trackerFPGA.chips->At(ichip);
            for (int icluster=0; icluster<trackerChip->numClust; ++icluster) {
               // cout<< ientry << "\t iFPGA = " << iFPGA << " trackerChip->address = " << std::setw(2) << (unsigned) trackerChip->address << " nstrips[" << icluster << "] = " << (unsigned) trackerChip->nstrips[icluster] << " nfirst[" << icluster << "] = " << (unsigned) trackerChip->nfirst[icluster] <<endl;
               *(geoHit.stripFPGA[iFPGA] + 64*trackerChip->address + (63 - trackerChip->nfirst[icluster])) = trackerChip->nstrips[icluster];
               ++nhits_total;
            }
         }
      }

      if (nhits_total == 0) continue;

      Int_t t_hit[4];
      Int_t v_hit[4];
      for (int ilayer=0; ilayer<4; ilayer++) {
         t_hit[ilayer] = 0;
         v_hit[ilayer] = 0;
      }

      for (int ilayer=0; ilayer<4; ++ilayer) {
         // t-layer
         for (int istrip=0; istrip<4*384; ++istrip) {
            Int_t nstrips = geoHit.tstrip[ilayer][istrip];
            if (nstrips > 0) for (int icurr=0; icurr<nstrips; ++icurr) {
               ++t_hit[ilayer];
            }
         }

         // v-layer right (lower chip addresses)
         for (int istrip=0; istrip<384; ++istrip) {
            Int_t nstrips = geoHit.vstrip[ilayer][istrip];
            if (nstrips > 0) for (int icurr=0; icurr<nstrips; ++icurr) {
               ++v_hit[ilayer];
            }
         }
         // v-layer left (higher chip addresses)
         for (int istrip=384; istrip<2*384; ++istrip) {
            Int_t nstrips = geoHit.vstrip[ilayer][istrip];
            if (nstrips > 0) for (int icurr=0; icurr<nstrips; ++icurr) {
               ++v_hit[ilayer];
            }
         }
      }

      bool good = true;
      for (int ilayer=0; ilayer<4; ++ilayer) {
         if (t_hit[ilayer] == 0) {
            good = false;
            break;
         }
         if (v_hit[ilayer] == 0) {
            good = false;
            break;
         }
      }

      if (good) {
         ngood++;
         cout<< ngood << "\t " << ientry <<endl;
      }

   }  // loop over entries

   tree->ResetBranchAddresses();
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
         event_flag = edisplay(evtNo,tree,wname);
         if (event_flag >= 0) last_event_plot = event_flag;
         ++evtNo;
         continue;
      }

      // 1) number
      std::stringstream ss(line);
      Int_t number = -1;
      if (ss >> number && ss.eof()) {
         evtNo = number;
         event_flag = edisplay(evtNo,tree,wname);
         if (event_flag >= 0) last_event_plot = event_flag;
         ++evtNo;
         continue;
      }

      // 2) not a number
      if (line.size() == 1)
      {
         // 2.1) input was just one character: interpret as a menu command
         switch (toupper(line[0])) {
            case '-':
               evtNo -= 2;          // previous event
               event_flag = edisplay(evtNo,tree,wname);
               if (event_flag >= 0) last_event_plot = event_flag;
               ++evtNo;
               break;
            case 'C':
               gPad->DrawClone();
               break;
            case 'S':
               if (getTree()) gPad->SaveAs(Form("%s_evt_%d.png", getTree()->GetCurrentFile()->GetName(), last_event_plot));
               else gPad->SaveAs(Form("evt_%d.png", evtNo));
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
