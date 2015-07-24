//////////////////////////////////////////////////
//
//  recoRun.cpp:
//
//  Created: 07/07/2015, Andriy Zatserklyaniy <zatserkl@fnal.gov>
//
//  See $ROOTSYS/main/src/rmain.cxx
//  
/////////////////////////////////////////////////

/*
Linux:
rootcint -f Reco_dict.cxx -c Reco.h Reco_linkdef.h
g++ -Wall `$ROOTSYS/bin/root-config --cflags --glibs` recoRun.cpp -o recoRun Reco_dict.cxx

Mac OS X:
rootcint -f Reco_dict.cxx -c Reco.h Reco_linkdef.h
clang++ -Wall `$ROOTSYS/bin/root-config --cflags --glibs` recoRun.cpp -o recoRun Reco_dict.cxx
*/

#include "Geometry.h"
#include "Reco.h"

// WEPL
#include "CalTV.h"
#include "Wepl.h"
#include <TH2.h>
#include <TF2.h>
#include <TProfile2D.h>

#include <TCanvas.h>
#include <TH2.h>
#include <TFile.h>
#include <TRint.h>

#include <string>

/// ClassImp(Hit2D);
/// ClassImp(SensorHit);
/// ClassImp(CRay2D);
/// //ClassImp(CRay);
/// ClassImp(Track2D);
/// //ClassImp(Track);
/// ClassImp(SuperTrack2D);
/// ClassImp(SuperTrack);
/// ClassImp(RecoEvent);

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

int main(int argc, char *argv[])
{
    int n = 1;

    const char* ifname = 0;
    if (argc > n) ifname = argv[n];
    else {
        cout<< "***Warning: Expected filename for the data" <<endl;
        cout<< "Usage: " << argv[0] << " ifname event1=0 event2=-1 dbname=rundb-May2015.dat tv_calib_fname=TVcalib.txt wet_calib_fname=wet5calibExp.txt debug=0" <<endl;
        return 0;
    }
    ++n;

    Int_t event1 = 0;
    if (argc > n) event1 = atoi(argv[n]);
    else {
        cout<< "***Warning: Expected the first event to process" <<endl;
        cout<< "Usage: " << argv[0] << " ifname event1=0 event2=-1 dbname=rundb-May2015.dat tv_calib_fname=TVcalib.txt wet_calib_fname=wet5calibExp.txt debug=0" <<endl;
        return 0;
    }
    ++n;

    Int_t event2 = 0;
    if (argc > n) event2 = atoi(argv[n]);
    else {
        cout<< "***Warning: Expected the last event to process" <<endl;
        cout<< "Usage: " << argv[0] << " ifname event1=0 event2=-1 dbname=rundb-May2015.dat tv_calib_fname=TVcalib.txt wet_calib_fname=wet5calibExp.txt debug=0" <<endl;
        return 0;
    }
    ++n;

    const char* dbname = "rundb-May2015.dat";
    if (argc > n) dbname = argv[n];
    else {
        cout<< "***Warning: Expected filename for rundb" <<endl;
        cout<< "Usage: " << argv[0] << " ifname event1=0 event2=-1 dbname=rundb-May2015.dat tv_calib_fname=TVcalib.txt wet_calib_fname=wet5calibExp.txt debug=0" <<endl;
        return 0;
    }
    ++n;

    const char* tv_calib_fname = "TVcalib.txt";
    if (argc > n) tv_calib_fname = argv[n];
    else {
        cout<< "***Warning: Expected filename for TV calibration coefficients" <<endl;
        cout<< "Usage: " << argv[0] << " ifname event1=0 event2=-1 dbname=rundb-May2015.dat tv_calib_fname=TVcalib.txt wet_calib_fname=wet5calibExp.txt debug=0" <<endl;
        return 0;
    }
    ++n;

    const char* wet_calib_fname = "wet5calibExp.txt";
    if (argc > n) wet_calib_fname = argv[n];
    else {
        cout<< "***Warning: Expected filename for WET calibration coefficients" <<endl;
        cout<< "Usage: " << argv[0] << " ifname event1=0 event2=-1 dbname=rundb-May2015.dat tv_calib_fname=TVcalib.txt wet_calib_fname=wet5calibExp.txt debug=0" <<endl;
        return 0;
    }
    ++n;

    bool debug = false;
    if (argc > n) debug = atoi(argv[n]);
    ++n;

    recoRun(ifname
            , event1, event2
            , dbname
            , tv_calib_fname
            , wet_calib_fname
            , debug
           );

    //
    //  uncomment to run live ROOT session
    //
    //TRint* theApp = new TRint("Rint", &argc, argv, 0, 0, 1); // do not show splash screen
    //theApp->Run();
    //delete theApp;

    return 0;
}
