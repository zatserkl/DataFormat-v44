/*
Linux:
rootcint -f Reco_dict.cxx -c Reco.h Reco_linkdef.h
g++ -Wall `$ROOTSYS/bin/root-config --cflags --glibs` writeReco.cpp -o writeReco Reco_dict.cxx

Mac OS X:
rootcint -f Reco_dict.cxx -c Reco.h Reco_linkdef.h
clang++ -Wall `$ROOTSYS/bin/root-config --cflags --glibs` writeReco.cpp -o writeReco Reco_dict.cxx
*/

// T-V correction, ADC-MeV convertion and WEPL reconstruction implemented  
//  into Reco.C by Andriy Zatserklyaniy,  and renamed weplReco
// by V.A.Bashkirov (on May 9, 2014, modified May 16, 2014)
// In Root, assuming Andriy's classes are loaded, type .L weplReco.C+  and
// run CalibTVEW("<file name>.reco.root") to execute and see control histo's;
// e.q.  CalibTVEW("pCTraw_Run_705.out.root.reco.root") 

#include "EventOutput.h"

#include "Reco.h"
#include "DataFormat.h"
#include "CalTV.h"
#include "Wepl.h"

#include <TFile.h>

#include <fstream>
#include <string>
#include <sstream>
#include <cstdio>

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

struct OutputFiles {
   Int_t run;
   Double_t angle;
   Int_t nevents;
   std::stringstream vname[4];
   std::stringstream tname[4];
   std::stringstream uname[4];
   std::stringstream wname;

   std::fstream* vfile[4];
   std::fstream* tfile[4];
   std::fstream* ufile[4];
   std::fstream* wfile;

   OutputFiles(): run(0), angle(0), nevents(0) {
   }
   ~OutputFiles() {
      close();
      for (int ilayer=0; ilayer<4; ++ilayer) {
         remove(vname[ilayer].str().c_str());
         remove(tname[ilayer].str().c_str());
         remove(uname[ilayer].str().c_str());
      }
      std::remove(wname.str().c_str());
   }
   void open() {
      //
      //    open file in out mode, close and open again in out/in mode
      //
      for (int ilayer=0; ilayer<4; ++ilayer) {
         vname[ilayer] << "temporary_file_v" << ilayer << "file-" << std::setfill('0') << std::setw(3) << angle << ".dat";
         tname[ilayer] << "temporary_file_t" << ilayer << "file-" << std::setfill('0') << std::setw(3) << angle << ".dat";
         uname[ilayer] << "temporary_file_u" << ilayer << "file-" << std::setfill('0') << std::setw(3) << angle << ".dat";

         vfile[ilayer] = new fstream(vname[ilayer].str().c_str(), std::ios::binary | std::ios::out);
         vfile[ilayer]->close();
         vfile[ilayer]->open(vname[ilayer].str().c_str(), std::ios::binary | std::ios::in | std::ios::out);

         tfile[ilayer] = new fstream(tname[ilayer].str().c_str(), std::ios::binary | std::ios::out);
         tfile[ilayer]->close();
         tfile[ilayer]->open(tname[ilayer].str().c_str(), std::ios::binary | std::ios::in | std::ios::out);

         ufile[ilayer] = new fstream(uname[ilayer].str().c_str(), std::ios::binary | std::ios::out);
         ufile[ilayer]->close();
         ufile[ilayer]->open(uname[ilayer].str().c_str(), std::ios::binary | std::ios::in | std::ios::out);
      }
      wname << "temporary_file_wfile-" << std::setfill('0') << std::setw(3) << angle << ".dat";
      wfile = new fstream(wname.str().c_str(), std::ios::binary | std::ios::out);
      wfile->close();
      wfile->open(wname.str().c_str(), std::ios::binary | std::ios::in | std::ios::out);
   }
   void close() {
      for (int ilayer=0; ilayer<4; ++ilayer) {
         vfile[ilayer]->close();
         tfile[ilayer]->close();
         ufile[ilayer]->close();
      }
      wfile->close();
   }
   void writeOutputFile(const char* uniname)
   {
       OutputFiles& outputFiles = *this;

       Int_t ngood = outputFiles.nevents;

       //
       // combine output files
       //

       cout<<endl<< ngood << " good events" <<endl;

       cout<< "look at the output file" <<endl;

       for (int ilayer=0; ilayer<4; ilayer++) {
           outputFiles.tfile[ilayer]->seekg(0, outputFiles.tfile[ilayer]->beg);
           outputFiles.vfile[ilayer]->seekg(0, outputFiles.vfile[ilayer]->beg);
           outputFiles.ufile[ilayer]->seekg(0, outputFiles.ufile[ilayer]->beg);
       }
       outputFiles.wfile->seekg(0, outputFiles.wfile->beg);

       // output file
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

       Float_t projectionAngle = outputFiles.angle;                      // take angle from outputFiles
       unifile.write((const char*) &projectionAngle, sizeof(Float_t));

       Float_t beamEnergy = 200;
       unifile.write((const char*) &beamEnergy, sizeof(Float_t));

       //-- Int_t acquisitionDate = start_time;              //--TODO
       Int_t acquisitionDate = 0;
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

       cout<< "combine the data files" <<endl;

       // write T first
       for (int ilayer=0; ilayer<4; ++ilayer) {
           std::istreambuf_iterator<char> begin_source(*outputFiles.tfile[ilayer]);
           std::istreambuf_iterator<char> end_source;
           std::ostreambuf_iterator<char> begin_dest(unifile);
           std::copy(begin_source, end_source, begin_dest);
       }

       // then V
       for (int ilayer=0; ilayer<4; ++ilayer) {
           std::istreambuf_iterator<char> begin_source(*outputFiles.vfile[ilayer]);
           std::istreambuf_iterator<char> end_source;
           std::ostreambuf_iterator<char> begin_dest(unifile);
           std::copy(begin_source, end_source, begin_dest);
       }

       // followed by U
       for (int ilayer=0; ilayer<4; ++ilayer) {
           std::istreambuf_iterator<char> begin_source(*outputFiles.ufile[ilayer]);
           std::istreambuf_iterator<char> end_source;
           std::ostreambuf_iterator<char> begin_dest(unifile);
           std::copy(begin_source, end_source, begin_dest);
       }
       // and WEPL
       std::istreambuf_iterator<char> begin_source(*outputFiles.wfile);
       std::istreambuf_iterator<char> end_source;
       std::ostreambuf_iterator<char> begin_dest(unifile);
       std::copy(std::istreambuf_iterator<char>(*outputFiles.wfile), std::istreambuf_iterator<char>(), std::ostreambuf_iterator<char>(unifile));

       //cout<< "OutputFiles: Wrote " << ngood << " events into output binary file " << uniname <<endl;
       outputFiles.close();
   }
};  // struct OutputFiles

Int_t writeReco(Int_t event1=0, Int_t event2=-1
               , TTree* tree=0
               , const char* dbname="rundb-May2015.dat"
               , const char* tv_calib_fname="TVcalib.txt"
               , const char* wet_calib_fname="wet5calibExp.txt"
               , OutputFiles* outputFiles=0
               , bool debug=false
               )
{
   bool do_wepl = true;
   if (!tv_calib_fname || !tv_calib_fname[0]) do_wepl = false;
   if (!wet_calib_fname || !wet_calib_fname[0]) do_wepl = false;
   if (!do_wepl) {
       cout<< "There are no TV-correction and/or WET calibration file provided. Nothing to do" <<endl;
       return 0;
   }

   if (!tree) tree = (TTree*) gDirectory->Get("t");                 // DataFormat tree
   if (!tree) {
      cout<< "Could not find tree t" <<endl;
      return 0;
   }

   const PCTEvent* pCTEvent = 0;
   tree->SetBranchAddress("event", &pCTEvent);

   RecoEvent* recoEvent = new RecoEvent();

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

   //
   // use geometry to fill up outputFiles data and to get u-positions of the layers
   //
   outputFiles->run = runHeader->GetRun();
   //--from-rundb-- outputFiles->angle = geometry->angle_;
   outputFiles->angle = runHeader->GetAngle();                  // take angle from the RunHeader instead of rundb
   //
   //   open temporary output files. The angle will be embedded into the file name
   //
   outputFiles->open();

   Double_t uv[4];
   Double_t ut[4];
   for (int ilayer=0; ilayer<4; ++ilayer) {
      uv[ilayer] = geometry->uv_[ilayer];
      ut[ilayer] = geometry->ut_[ilayer];
      if (debug) cout<< "uv[" << ilayer << "] = " << uv[ilayer] << " ut[" << ilayer << "] = " << ut[ilayer] <<endl;
   }

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

   //
   // WEPL reconstruction
   //

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
         return 0;
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
         return 0;
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
   cout<< "event1 = " << event1 << " event2 = " << event2 <<endl;

   EventOutput eventOutput;

   //
   //   event loop over DataFormat tree
   //
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

      //-- recoEvent->clear();   // probably we don't need it, the recoEvent->Extract(reco) will do the job

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

            //
            // angle cut
            //
            if (superTrack->angle > 0.070) continue;

            for (int ical=0; ical<5; ++ical) {
               adc = recoEvent->a[ical];
               adc = adc - ped[ical]; 
               //	if (adc < 100) adc = 0;

               // Apply TV correction and convert ADC channel to E in MeV
               if (debug) cout<< "Apply TV correction and convert ADC channel to E in MeV" <<endl;
               Estage[ical]=f2cal[ical]->CalTVf(superTrack->T(ucal),superTrack->V(ucal),adc);   
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

            //
            //  writeReco stuff
            //
            ++outputFiles->nevents;
            eventOutput.good = kTRUE;
            for (int ilayer=0; ilayer<4; ++ilayer) {
                eventOutput.u[ilayer] = ut[ilayer];
                eventOutput.v[ilayer] = superTrack->V(ut[ilayer]);
                eventOutput.t[ilayer] = superTrack->T(ut[ilayer]);
            }
            eventOutput.wepl = Wet;

            for (int ilayer=0; ilayer<4; ilayer++) outputFiles->tfile[ilayer]->write((const char*) &eventOutput.t[ilayer], sizeof(Float_t));
            for (int ilayer=0; ilayer<4; ilayer++) outputFiles->vfile[ilayer]->write((const char*) &eventOutput.v[ilayer], sizeof(Float_t));
            for (int ilayer=0; ilayer<4; ilayer++) outputFiles->ufile[ilayer]->write((const char*) &eventOutput.u[ilayer], sizeof(Float_t));
            outputFiles->wfile->write((const char*) &eventOutput.wepl, sizeof(Float_t));
         }
      }
   }

   cout<< "Wrote " << outputFiles->nevents << " events into temporary files" <<endl;
   return outputFiles->nevents;
}



Int_t writeReco(const char* ifname
               , Int_t event1=0, Int_t event2=-1
               , const char* dbname="rundb-May2015.dat"
               , const char* tv_calib_fname="TVcalib.txt"
               , const char* wet_calib_fname="wet5calibExp.txt"
               , bool debug=false
              )
{
   TFile* ifile = new TFile(ifname);
   if (!ifile) {
      cout<< "Could not open file " << ifname <<endl;
      return 0;
   }

   TTree* tree = (TTree*) gDirectory->Get("t");                 // DataFormat tree
   if (!tree) {
      cout<< "Could not find tree DataFormat tree t" <<endl;
      return 0;
   }

   OutputFiles outputFiles;

   Int_t ntot = writeReco(event1, event2
             , tree
             , dbname
             , tv_calib_fname
             , wet_calib_fname
             , &outputFiles
             , debug
            );

   // output file
   // const char* uniname = Form("%s-%0.0f.bin",ifname,outputFiles.angle);    // angle from outputFiles
   const char* uniname = Form("%s.bin",ifname);

   //
   // combine output files
   //

   outputFiles.writeOutputFile(uniname);
   cout<< "Wrote " << outputFiles.nevents << " good events out of " << ntot << " into output binary file " << uniname <<endl;
   return outputFiles.nevents;
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

    Int_t event2 = -1;
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

    writeReco(ifname
            , event1, event2
            , dbname
            , tv_calib_fname 
            , wet_calib_fname 
            , debug
           );

    return 0;
}
