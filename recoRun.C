// Andriy Zatserklyaniy, April 17, 2014

#include "Geometry.h"
#include "Reco.h"
//--------------------- #include "Reco-gap.h"

// WEPL
#include "CalTV.h"
#include "Wepl.h"
#include <TH2.h>
#include <TF2.h>
#include <TProfile2D.h>

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

void recoRun(const char* ifname
             , Int_t event1=0, Int_t event2=-1
             , const char* dbname="rundb-May2015.dat"
             , const char* tv_calib_fname="TVcalib.txt"
             , const char* wet_calib_fname="wet5calibExp.txt"
             , bool debug=false
             )
{
   bool do_wepl = true;
   if (!tv_calib_fname || !tv_calib_fname[0]) do_wepl = false;
   if (!wet_calib_fname || !wet_calib_fname[0]) do_wepl = false;

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

   // const RecoEvent* recoEvent = 0;                 // pointer to the reco event buffer
   RecoEvent* recoEvent = 0;                          // pointer to the reco event buffer

   std::string ofname;
   if (event1 != 0 || event2 != tree->GetEntries()-1) ofname = Form("%s-%d-%d.reco.root",ifname,event1,event2);
   else ofname = Form("%s.reco.root",ifname);

   //cout<< "trying to open output file " << ofname.str() <<endl;
   TFile* ofile = new TFile(ofname.c_str(), "recreate");
   TTree* otree = new TTree("r", Form("Reconstruction of %s",ifname));
   otree->SetMarkerColor(602);
   otree->Branch("revent", "RecoEvent", &recoEvent);

   RunHeader* runHeader = new RunHeader(*((RunHeader*) tree->GetUserInfo()->First()));  // use copy constructor
   cout<< "runHeader->GetRun() = " << runHeader->GetRun() <<endl;
   time_t start_time = runHeader->GetTime();
   cout<< "run start time: " << std::ctime(&start_time);
   cout<< "program version is " << runHeader->GetVersion() <<endl;
   if (runHeader->GetTimeTag()) cout<< "event time tag was written out" <<endl;
   else cout<< "event time tag was not written out" <<endl;
   cout<< "runHeader->GetAngle() = " << runHeader->GetAngle() <<endl;
   cout<<endl;

   //
   // read the database
   //
   Geometry* geometry = new Geometry(runHeader->GetRun(), dbname);

   //--new-- // add the angle and the phantom name to the runHeader
   //--new-- runHeader->SetAngle(geometry->angle_);
   //--new-- runHeader->SetPhantom(geometry->phantom_.Data());

   // pass the runHeader from the input to the output tree
   otree->GetUserInfo()->Add(runHeader);

   PCTSensors* pCTSensors = new PCTSensors(geometry);

   //-- Geant4 setup
   //-- Double_t RbeamSpot = 10.;                                            // Geant4
   //-- BeamSpot beamSpotIn(0,0,-3500, RbeamSpot);                           // Geant4
   //-- BeamSpot beamSpotOut(0,0,-3500, RbeamSpot);                          // Geant4

   // Double_t RbeamSpot = 0.5*(14.16 + 17.08);                         // Sep2014 beam test
   // BeamSpot beamSpotIn(35.52,-1.07,-3500, RbeamSpot);                // Sep2014 beam test
   // BeamSpot beamSpotOut(13.76,-6.96,-3500, 140.);                    // Sep2014 beam test

   Double_t RbeamSpot = 10.;                                       // May2015 beam test
   BeamSpot beamSpotIn(-10.5,9.7,-2100, RbeamSpot);                // May2015 beam test
   BeamSpot beamSpotOut(-5.0,13.0,-2100, RbeamSpot);               // May2015 beam test

   Int_t tnchan = 400, vnchan = 20;
   Double_t tlow = -200., tup = 200., vlow = -50, vup = 50;
   TH2F* hcal_a[5];
   for (int i=0; i<5; ++i) hcal_a[i] = new TH2F(Form("hcal%d_a",i), Form("Cal response for channel %d",i), tnchan,tlow,tup, vnchan,vlow,vup);
   TH2F* hcal_i[5];
   for (int i=0; i<5; ++i) hcal_i[i] = new TH2F(Form("hcal%d_i",i), Form("Nevents for Cal response for channel %d",i), tnchan,tlow,tup, vnchan,vlow,vup);

   //
   // WEPL reconstruction
   //

   // Book histo's 
   TProfile2D* hcal[5];
   for (int i=0; i<5; ++i) hcal[i] = new TProfile2D(Form("hcal%d",i), Form("Stage %d response, MeV, t-v corrected",i),60,-180,180,18,-45,45,0,80);
   TH1F* hcalr[5];
   for (int i=0; i<5; ++i) hcalr[i] = new TH1F(Form("hcalr%d",i), Form("Stage %d row",i),400, 0,8000);  	
   TH1F* hcalc[5];
   for (int i=0; i<5; ++i) hcalc[i] = new TH1F(Form("hcalc%d",i), Form("Stage %d tv-corrected",i),400, 0,100);   
   TH1F *h7 = new TH1F("h7","Reconstructed WET",300,-10.5,289.5);
   // h7->SetStats(kFALSE);
   TProfile2D* hwepl= new TProfile2D("hwepl","Reconstructed WEPL",180,-6.35*30,6.35*30,45,-45,45,-10,180);	

   // ADC Pedestals

   // Double_t ped[5] = {9.645, -20.484, -201.987, 62.966, -7.747};     // Celeste data
   // Jul2014-- Double_t ped[5] = {121.3, -71.5, -1137, 346.2, -49.};     // New pedestals (x6, reduced data)
   // Double_t ped[5] = {431,-130,-20,224,60};     // Sept. 2014 pedestals (x6, reduced data)
   Double_t ped[5] = {549,92,204,575,385};     // New pedestals, May 2015 (x6, reduced data)

   //   Prepare stuff for TV correction and convertion ADC->energy(MeV)

   Double_t ucal = 216.9 + 40;            // approx position for the calorimeter entrance  
   Float_t par[5];
   Float_t adc;
   Float_t Estage[5];

   CalF* f2cal[5]; 	                     // Calibration functions for 5 stages 

   if (do_wepl) {
      // open TV-correction parameters file
      std::ifstream TVcalfile(tv_calib_fname);
      if (!TVcalfile) {
         cout<< "Could not open TV-correction file " << tv_calib_fname <<endl;
         return;
      }

      for (int i=0; i<5; ++i) {            // read fit parameters:
         TVcalfile >> par[0] >> par[1] >> par[2] >> par[3] >> par[4] ;
         f2cal[i]=new CalF(i,par);        // initialize t-v calibration  function  
      }
      TVcalfile.close(); 
   }

   // Prepare stuff for WEPL calibration 	

   Float_t wepl_par[45];
   Float_t Wet;

   if (do_wepl) {
      // open text file with WEPL calibration data (parameters of 9 pol4 curves)   
      std::ifstream WEPLcalfile(wet_calib_fname); 
      if (!WEPLcalfile) {
         cout<< "Could not find file " << wet_calib_fname <<endl;
         return;
      }

      for (int i=0; i<45; ++i) WEPLcalfile >> wepl_par[i];
      WEPLcalfile.close(); 
   }

   Wepl* WEPL = 0;
   if (do_wepl) {
      WEPL = new Wepl(wepl_par);                // initialize WEPL calibr. 
      WEPL->SetEthresholds(1,.99,1,1,1);        // Set all stage thresholds to 1 MeV
   }

   //----------------- end of WEPL reconstruction stuff -----------------

   Int_t nevents_with_tracks = 0;

   if (event2 < event1) event2 = tree->GetEntries()-1;

   for (int ientry=event1; ientry<=event2; ++ientry)
   {
      if (tree->LoadTree(ientry) < 0) {
         cout<< "Could not load event " << ientry <<endl;
         break;
      }
      tree->GetEntry(ientry);

      if (false
          || (ientry < 10)
          || (ientry < 10000 && ientry%1000 == 0)
          || (ientry%100000 == 0)
      ) cout<< "---> processing entry " << ientry <<endl; 

      //--new-- Reco reco(geometry, &beamSpot, pCTSensors, pCTEvent, ientry, debug);
      Reco reco(geometry, pCTSensors, &beamSpotIn, &beamSpotOut, pCTEvent, ientry, debug);
      reco.Tracking();

      if (debug) cout<< "call recoEvent->Extract(reco)" <<endl;
      recoEvent->Extract(reco);
      if (recoEvent->track->GetLast()+1 > 0) ++nevents_with_tracks;

      recoEvent->wepl = -2000;
      if (do_wepl) {
         //
         //    calculate WEPL and assign to the recoEvent->wepl
         //
         if (debug) cout<< "calculate WEPL: recoEvent->track->GetLast()+1 = " << recoEvent->track->GetLast()+1 <<endl;
         if (recoEvent->track->GetLast()+1 > 0)
         {
            const SuperTrack* superTrack = (const SuperTrack*) recoEvent->track->At(0);   // take the first reconstructed track

            for (int ical=0; ical<5; ++ical) {
               adc = recoEvent->a[ical];
               adc = adc - ped[ical]; 
               //	if (adc < 100) adc = 0;
               hcalr[ical]->Fill(adc);  

               // Apply TV correction and convert ADC channel to E in MeV
               if (debug) cout<< "Apply TV correction and convert ADC channel to E in MeV" <<endl;
               Estage[ical]=f2cal[ical]->CalTVf(superTrack->T(ucal),superTrack->V(ucal),adc);   

               // Fill TV-corr control histo's
               hcalc[ical]->Fill(Estage[ical]);          
               hcal[ical]->Fill(superTrack->T(ucal), superTrack->V(ucal),Estage[ical]);
            }

            // Get Wet=WEPL from Estage
            if (debug) cout<< "Get Wet=WEPL from Estage" <<endl;
            Wet=WEPL->EtoWEPL(Estage);
            if (debug) cout<< "Wet = " << Wet <<endl;

            //-- if(Wet>999. || Wet<-999.) continue;
            //--new-- recoEvent->SetWEPL(Wet);
            if (debug) cout<< "assign Wet to recoEvent->wepl" <<endl;
            if(Wet>-999. && Wet<999.) recoEvent->wepl = Wet;
            if (debug) cout<< "   done" <<endl;

            // Fill WEPL-calib control histo's
            h7->Fill(Wet);  
            hwepl->Fill(superTrack->T(0), superTrack->V(0),Wet);
            if (debug) cout<< "finished with WEPL" <<endl;
         }
      }

      if (debug) cout<< "otree->Fill()" <<endl;
      otree->Fill();
      if (debug) cout<< "   done otree->Fill()" <<endl;

      for (int itrack=0; itrack<recoEvent->track->GetLast()+1; ++itrack) {
         const SuperTrack* superTrack = (const SuperTrack*) recoEvent->track->At(itrack);
         for (int ical=0; ical<5; ++ical) {
            Float_t adc = recoEvent->a[ical];
            if (adc < 100) adc = 0;
            hcal_a[ical]->Fill(superTrack->T(-128.), superTrack->V(-128.), adc);
            hcal_i[ical]->Fill(superTrack->T(-128.), superTrack->V(-128.));
         }
      }
   }

   for (int ical=0; ical<5; ++ical) hcal_a[ical]->Divide(hcal_i[ical]);

   new TCanvas;
   hcal_a[0]->DrawCopy("colz");

   cout<< "Reconstructed tracks in " << nevents_with_tracks << " events" <<endl;

   if (otree->GetCurrentFile()) cout<< "write " << otree->GetEntries() << " events into output file " << otree->GetCurrentFile()->GetName() <<endl;
   else cout<< "write " << otree->GetEntries() << " events into output file " << ofile->GetName() <<endl;
   ofile->Write();
}

void recoEvent(Int_t event, TTree* tree=0
               , bool debug=true
               , const char* dbname="rundb-May2015.dat"
              )
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
   cout<< "runHeader->GetAngle() = " << runHeader->GetAngle() <<endl;
   cout<<endl;

   //
   // read the database
   //
   Geometry* geometry = new Geometry(runHeader->GetRun(), dbname);

   const PCTEvent* pCTEvent = 0;
   tree->SetBranchAddress("event", &pCTEvent);

   RecoEvent* recoEvent = new RecoEvent();

   //--new-- PCTSensors* pCTSensors = new PCTSensors(geometry, debug);
   PCTSensors* pCTSensors = new PCTSensors(geometry);

   Double_t RbeamSpot = 0.5*(14.16 + 17.08);                         // Sep2014 beam test
   BeamSpot beamSpotIn(35.52,-1.07,-3500, RbeamSpot);                // Sep2014 beam test
   BeamSpot beamSpotOut(13.76,-6.96,-3500, 140.);                    // Sep2014 beam test

   if (tree->LoadTree(event) < 0) {
      cout<< "Could not load event " << event <<endl;
      return;
   }
   tree->GetEntry(event);

   Reco reco(geometry, pCTSensors, &beamSpotIn, &beamSpotOut, pCTEvent, event, debug);
   reco.Tracking();

   cout<< "--> call Extract" <<endl;
   recoEvent->Extract(reco);
   cout<< "--> recoEvent->track->GetLast()+1 = " << recoEvent->track->GetLast()+1 <<endl;

   for (int itrack=0; itrack<recoEvent->track->GetLast()+1; ++itrack) {
      const SuperTrack* superTrack = (const SuperTrack*) recoEvent->track->At(itrack);
      cout<< itrack << " angle = " << superTrack->angle <<endl;
      Double_t v = 0, t = 0;
      superTrack->at(266.9, v, t);
      cout<< "v = " << v << " t = " << t <<endl;
      //--new-- // cout<< *superTrack <<endl;
      //--new-- cout<< itrack << "\tangle = " << superTrack->angle << " distance = " << superTrack->Distance() <<endl;
      //--new-- cout<< "\tsuperTrack->itrack_->vTrack_->hit1_ = " << *superTrack->itrack_->vTrack_->hit1_ << "\tsuperTrack->itrack_->vTrack_->hit2_ = " << *superTrack->itrack_->vTrack_->hit2_ <<endl;
      //--new-- cout<< "\tsuperTrack->itrack_->tTrack_->hit1_ = " << *superTrack->itrack_->tTrack_->hit1_ << "\tsuperTrack->itrack_->tTrack_->hit2_ = " << *superTrack->itrack_->tTrack_->hit2_ <<endl;
      //--new-- cout<< "\tsuperTrack->otrack_->vTrack_->hit1_ = " << *superTrack->otrack_->vTrack_->hit1_ << "\tsuperTrack->otrack_->vTrack_->hit2_ = " << *superTrack->otrack_->vTrack_->hit2_ <<endl;
      //--new-- cout<< "\tsuperTrack->otrack_->tTrack_->hit1_ = " << *superTrack->otrack_->tTrack_->hit1_ << "\tsuperTrack->otrack_->tTrack_->hit2_ = " << *superTrack->otrack_->tTrack_->hit2_ <<endl;
   }
   for (int ichan=0; ichan<5; ++ichan) cout<< "a[" << ichan << "] = " << recoEvent->a[ichan] << " "; cout<<endl;

   // cout<< "All hits from the Sensor::poolSensorHit_" <<endl;
   // for (int ihit=0; ihit<Sensor::poolSensorHit_->GetLast()+1; ++ihit) {
   //    const SensorHit* hit = (const SensorHit*) Sensor::poolSensorHit_->At(ihit);
   //    cout<< ihit << "\t" << *hit <<endl;
   // }
}

//
// NB: display does not use Reco
//

Int_t display(Int_t event, TTree* tree=0, const char* dbname="rundb-May2015.dat", const char* wname="event_display")
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
   cout<< "runHeader->GetAngle() = " << runHeader->GetAngle() <<endl;
   cout<<endl;

   //
   // read the database
   //
   Geometry* geometry = new Geometry(runHeader->GetRun(), dbname);

   //BeamSpot beamSpot(0,0,-3500);    // the display does not use Reco (and the BeamSpot)

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

void eloop(Int_t evtNo=0, TTree* tree=0, const char* dbname="rundb-May2015.dat", const char* wname="event_display")
{
   if (tree == 0) {
      tree = (TTree*) gDirectory->Get("t");
      if (!tree) {
         cout<< "Could not find tree \"t\" in the current directory " << gDirectory->GetName() <<endl;
         return;
      }

      // RunHeader* runHeader = (RunHeader*) tree->GetUserInfo()->First();
      // cout<< "runHeader->GetRun() = " << runHeader->GetRun() <<endl;
      // time_t start_time = runHeader->GetTime();
      // cout<< "run start time: " << std::ctime(&start_time);
      // cout<< "program version is " << runHeader->GetVersion() <<endl;
      // if (runHeader->GetTimeTag()) cout<< "event time tag was written out" <<endl;
      // else cout<< "event time tag was not written out" <<endl;
      // cout<< "runHeader.GetAngle() = " << runHeader.GetAngle() <<endl;
      // cout<<endl;
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
         event_flag = display(evtNo,tree,dbname,wname);
         if (event_flag >= 0) last_event_plot = event_flag;
         ++evtNo;
         continue;
      }

      // 1) number
      std::stringstream ss(line);
      Int_t number = -1;
      if (ss >> number && ss.eof()) {
         evtNo = number;
         event_flag = display(evtNo,tree,dbname,wname);
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
               event_flag = display(evtNo,tree,dbname,wname);
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
