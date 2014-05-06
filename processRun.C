// Andriy Zatserklyaniy, April 17, 2014

#include "Geometry.h"
#include "Reco.h"

#include <TCanvas.h>
#include <TH2.h>
#include <TFile.h>

#include <string>

//-- ClassImp(RecoEvent);

//
// static pointer to the geometry
//
Geometry* Geometry::geometry_ = 0;
//
// static hit pool
//
TClonesArray* Sensor::poolSensorHit_ = 0;
TClonesArray* Reco::poolTrack2D_ = 0;
TClonesArray* Reco::poolSuperTrack2D_ = 0;
TClonesArray* Reco::poolSuperTrack_ = 0;

void reco(const char* ifname, Int_t event1=0, Int_t event2=-1, const char* dbname="rundb-Mar2014.dat")
{
   Bool_t debug = kFALSE;
   // Bool_t debug = kTRUE;

   TFile* ifile = new TFile(ifname);
   if (!ifile) {
      cout<< "Could not open file " << ifname <<endl;
      return;
   }

   TTree* tree = (TTree*) gDirectory->Get("t");
   if (!tree) {
      cout<< "Could not find tree t" <<endl;
      return;
   }

   if (event2 < event1) event2 = tree->GetEntries()-1;

   const PCTEvent* pCTEvent = 0;
   tree->SetBranchAddress("event", &pCTEvent);

   RecoEvent* recoEvent = new RecoEvent();

   std::string ofname;
   if (event1 != 0 || event2 != tree->GetEntries()-1) Form("%s-%d-%d.reco.root",ifname,event1,event2);
   else ofname = ofname = Form("%s.reco.root",ifname);

   TFile* ofile = new TFile(ofname.c_str(), "recreate");
   TTree* otree = new TTree("r", Form("Reconstruction of %s",ifname));
   otree->SetMarkerColor(2);
   otree->Branch("revent", "RecoEvent", &recoEvent);

   RunHeader* runHeader = (RunHeader*) tree->GetUserInfo()->First();
   cout<< "runHeader->GetRun() = " << runHeader->GetRun() <<endl;
   time_t start_time = runHeader->GetTime();
   cout<< "run start time: " << std::ctime(&start_time);
   cout<< "program version is " << runHeader->GetVersion() <<endl;
   if (runHeader->GetTimeTag()) cout<< "event time tag was written out" <<endl;
   else cout<< "event time tag was not written out" <<endl;
   cout<<endl;

   //
   // read the database
   //
   Geometry* geometry = new Geometry(runHeader->GetRun(), dbname);

   PCTSensors* pCTSensors = new PCTSensors(geometry);

   Int_t tnchan = 400, vnchan = 10;
   Double_t tlow = -200., tup = 200., vlow = -50, vup = 50;
   TH2F* hcal[5];
   for (int i=0; i<5; ++i) hcal[i] = new TH2F(Form("hcal%d",i), Form("Cal response for channel %d",i), tnchan,tlow,tup, vnchan,vlow,vup);
   TH2F* hcal_i[5];
   for (int i=0; i<5; ++i) hcal_i[i] = new TH2F(Form("hcal%d_i",i), Form("Nevents for Cal response for channel %d",i), tnchan,tlow,tup, vnchan,vlow,vup);

   if (event2 < event1) event2 = tree->GetEntries()-1;

   for (int ientry=event1; ientry<=event2; ++ientry)
   {
      if (tree->LoadTree(ientry) < 0) {
         cout<< "Could not load event " << ientry <<endl;
         return;
      }
      tree->GetEntry(ientry);

      if (false
          || (ientry < 10)
          || (ientry < 10000 && ientry%1000 == 0)
          || (ientry%100000 == 0)
      ) cout<< "---> processing entry " << ientry <<endl; 

      Reco reco(geometry, pCTSensors, pCTEvent, ientry, debug);
      //cout<< "------------------------ generate 2D tracks ---" <<endl;
      reco.GenerateTracks2D();
      //cout<< "------------------------ generate 2D SuperTracks ---" <<endl;
      reco.GenerateSuperTracks2D();
      //cout<< "------------------------ generate super tracks ---" <<endl;
      reco.GenerateSuperTracks();

      // cout<< "recoEvent: superTracks_.size() = " << reco.superTracks_.size() <<endl;

      // for (std::list<const SuperTrack*>::const_iterator it=reco.superTracks_.begin(); it!=reco.superTracks_.end(); ++it) {
      //    cout<< std::distance<std::list<const SuperTrack*>::const_iterator>(reco.superTracks_.begin(), it) << " angle = " << (*it)->angle <<endl;
      //    // const SuperTrack* superTrack = *it;
      //    // cout<< std::distance<std::list<const SuperTrack*>::const_iterator>(reco.superTracks_.begin(), it) << " angle = " << superTrack->angle <<endl;
      // }

      recoEvent->Extract(reco);

      // for (int itrack=0; itrack<recoEvent->track->GetLast()+1; ++itrack) {
      //    const SuperTrack* superTrack = (const SuperTrack*) recoEvent->track->At(itrack);
      //    cout<< itrack << " angle = " << superTrack->angle <<endl;
      //    Double_t v = 0, t = 0;
      //    superTrack->at(266.9, v, t);
      //    cout<< "v = " << v << " t = " << t <<endl;
      // }
      otree->Fill();

      for (int itrack=0; itrack<recoEvent->track->GetLast()+1; ++itrack) {
         const SuperTrack* superTrack = (const SuperTrack*) recoEvent->track->At(itrack);
         for (int ical=0; ical<5; ++ical) {
            Float_t adc = recoEvent->a[ical];
            if (adc < 100) adc = 0;
            // hcal[ical]->Fill(superTrack->tcal, superTrack->vcal, adc);
            // hcal_i[ical]->Fill(superTrack->tcal, superTrack->vcal);
            hcal[ical]->Fill(superTrack->T(-128.), superTrack->V(-128.), adc);
            hcal_i[ical]->Fill(superTrack->T(-128.), superTrack->V(-128.));
         }
      }

      // //----- debug -----
      // for (int ihit=0; ihit<Sensor::poolSensorHit_->GetLast()+1; ++ihit) {
      //    const SensorHit* hit = (const SensorHit*) Sensor::poolSensorHit_->At(ihit);
      //    cout<< ihit << "\t" << *hit <<endl;
      // }
   }

   //new TCanvas;
   //hcal[0]->DrawCopy("colz");

   for (int ical=0; ical<5; ++ical) hcal[ical]->Divide(hcal_i[ical]);

   new TCanvas;
   hcal[0]->DrawCopy("colz");

   cout<< "write " << otree->GetEntries() << " events into output file " << otree->GetCurrentFile()->GetName() <<endl;
   ofile->Write();
}

void recoEvent(Int_t event=100128, TTree* tree=0, bool debug=true, const char* dbname="rundb-Mar2014.dat")
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

   //
   // read the database
   //
   Geometry* geometry = new Geometry(runHeader->GetRun(), dbname);

   const PCTEvent* pCTEvent = 0;
   tree->SetBranchAddress("event", &pCTEvent);

   RecoEvent* recoEvent = new RecoEvent();

   PCTSensors* pCTSensors = new PCTSensors(geometry);

   if (tree->LoadTree(event) < 0) {
      cout<< "Could not load event " << event <<endl;
      return;
   }
   tree->GetEntry(event);

   Reco reco(geometry, pCTSensors, pCTEvent, event, debug);
   cout<< "------------------------ generate 2D tracks ---" <<endl;
   reco.GenerateTracks2D();
   cout<< "------------------------ generate 2D SuperTracks ---" <<endl;
   reco.GenerateSuperTracks2D();
   cout<< "------------------------ generate super tracks ---" <<endl;
   reco.GenerateSuperTracks();

   // cout<< "------------------------ apply filter for distance at u=0 ---" <<endl;
   // Double_t r = 10.; // mm
   // reco.FilterDistanceU0(r);

   // cout<< "recoEvent: superTracks_.size() = " << reco.superTracks_.size() <<endl;

   for (std::list<const SuperTrack*>::const_iterator it=reco.superTracks_.begin(); it!=reco.superTracks_.end(); ++it) {
      cout<< std::distance<std::list<const SuperTrack*>::const_iterator>(reco.superTracks_.begin(), it) << " angle = " << (*it)->angle <<endl;
      // const SuperTrack* superTrack = *it;
      // cout<< std::distance<std::list<const SuperTrack*>::const_iterator>(reco.superTracks_.begin(), it) << " angle = " << superTrack->angle <<endl;
   }

   recoEvent->Extract(reco);

   for (int itrack=0; itrack<recoEvent->track->GetLast()+1; ++itrack) {
      const SuperTrack* superTrack = (const SuperTrack*) recoEvent->track->At(itrack);
      cout<< itrack << " angle = " << superTrack->angle <<endl;
      Double_t v = 0, t = 0;
      superTrack->at(266.9, v, t);
      cout<< "v = " << v << " t = " << t <<endl;
   }
   for (int ichan=0; ichan<5; ++ichan) cout<< "a[" << ichan << "] = " << recoEvent->a[ichan] << " "; cout<<endl;
}

Int_t display(Int_t event, TTree* tree=0, const char* dbname="rundb-Mar2014.dat", const char* wname="event_display")
{
   if (!tree) tree = (TTree*) gDirectory->Get("t");
   if (!tree) {
      cout<< "Could not find tree t" <<endl;
      return -1;
   }

   RunHeader* runHeader = (RunHeader*) tree->GetUserInfo()->First();
   cout<< "runHeader->GetRun() = " << runHeader->GetRun() <<endl;
   time_t start_time = runHeader->GetTime();
   cout<< "run start time: " << std::ctime(&start_time);
   cout<< "program version is " << runHeader->GetVersion() <<endl;
   if (runHeader->GetTimeTag()) cout<< "event time tag was written out" <<endl;
   else cout<< "event time tag was not written out" <<endl;
   cout<<endl;

   //
   // read the database
   //
   Geometry* geometry = new Geometry(runHeader->GetRun(), dbname);

   PCTSensors pCTSensors(geometry);

   const PCTEvent* pCTEvent = 0;
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
   //std::vector<Double_t> v_hits[4];       // vhits for all 4 layers
   for (int ilayer=0; ilayer<4; ++ilayer) {
      for (int isensor=0; isensor<2; ++isensor) {
         //v_hits[ilayer].insert(v_hits[ilayer].end(), pCTSensors.vSensor[ilayer][isensor].hits_.begin(), pCTSensors.vSensor[ilayer][isensor].hits_.end());
         n_v_hits += pCTSensors.vSensor[ilayer][isensor].hits_.size();
      }
   }
   // t-hits
   //std::vector<Double_t> t_hits[4];
   for (int ilayer=0; ilayer<4; ++ilayer) {
      for (int isensor=0; isensor<4; ++isensor) {
         //t_hits[ilayer].insert(t_hits[ilayer].end(), pCTSensors.tSensor[ilayer][isensor].hits_.begin(), pCTSensors.tSensor[ilayer][isensor].hits_.end());
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

   /// Double_t uv[4] = {-217.7, -167.6,  166.8,  216.9};
   /// Double_t ut[4] = {-211.8, -161.8,  161.0,  211.0};
   const Double_t* uv = geometry->uv_;
   const Double_t* ut = geometry->uv_;

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

   TLine line;
   // forward
   // line.DrawLine(-217.70,-220.11, -217.70,215.24);
   // line.DrawLine(-211.80,-220.11, -211.80,215.24);
   // line.DrawLine(-167.60,-224.11, -167.60,211.24);
   // line.DrawLine(-161.80,-224.11, -161.80,211.24);
   // // back
   // line.DrawLine( 217.70,-220.11,  217.70,215.24);
   // line.DrawLine( 211.80,-220.11,  211.80,215.24);
   // line.DrawLine( 167.60,-224.11,  167.60,211.24);
   // line.DrawLine( 161.80,-224.11,  161.80,211.24);

   // forward
   line.DrawLine(uv[0],-220.11, uv[0],215.24);
   line.DrawLine(ut[0],-220.11, ut[0],215.24);
   line.DrawLine(ut[1],-224.11, ut[1],211.24);
   line.DrawLine(uv[1],-224.11, uv[1],211.24);
   // back
   line.DrawLine(ut[2],-220.11, ut[2],215.24);
   line.DrawLine(uv[2],-220.11, uv[2],215.24);
   line.DrawLine(ut[3],-224.11, ut[3],211.24);
   line.DrawLine(uv[3],-224.11, uv[3],211.24);

   TMarker v_marker;
   v_marker.SetMarkerStyle(24);
   v_marker.SetMarkerColor(4);
   TMarker t_marker;
   t_marker.SetMarkerStyle(24);
   t_marker.SetMarkerColor(2);

   for (int ilayer=0; ilayer<4; ++ilayer) {
      for (unsigned ihit=0; ihit<pCTSensors.v_hits[ilayer].size(); ++ihit) v_marker.DrawMarker(uv[ilayer], pCTSensors.v_hits[ilayer][ihit]->pos_);
      for (unsigned ihit=0; ihit<pCTSensors.t_hits[ilayer].size(); ++ihit) t_marker.DrawMarker(ut[ilayer], pCTSensors.t_hits[ilayer][ihit]->pos_);
   }

   TLine v_line;
   v_line.SetLineColor(1);
   TLine t_line;
   t_line.SetLineColor(2);

   if (pCTSensors.v_hits[0].size() == 1 && pCTSensors.v_hits[3].size() == 1) v_line.DrawLine(uv[0],pCTSensors.v_hits[0][0]->pos_, uv[3],pCTSensors.v_hits[3][0]->pos_);
   if (pCTSensors.t_hits[0].size() == 1 && pCTSensors.t_hits[3].size() == 1) t_line.DrawLine(ut[0],pCTSensors.t_hits[0][0]->pos_, ut[3],pCTSensors.t_hits[3][0]->pos_);

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
