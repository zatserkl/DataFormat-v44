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
#include <TCanvas.h>
#include <TH2.h>
#include <TF2.h>
#include <TProfile2D.h>
#include <TFile.h>
#include <TStyle.h>

#include <fstream>
#include <string>
#include <sstream>
#include <cstdio>

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
      wname << "temporary_file_wfile.dat";
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
};

void writeReco(Int_t event1=0, Int_t event2=-1, TTree* tree=0, OutputFiles* outputFiles=0, const char* dbname="rundb-Sep2014.dat", bool plot=false)
{
   Bool_t debug = kFALSE;
   if (debug) cout<< "debug in on" <<endl;

   if (!tree) tree = (TTree*) gDirectory->Get("r");
   if (!tree) {
      cout<< "Could not find tree r" <<endl;
      return;
   }

   // if (tree->GetCurrentFile()) cout<< "tree->GetCurrentFile()->GetName() = " << tree->GetCurrentFile()->GetName() <<endl;
   // else {
   //    cout<< "the current name is not available" <<endl;
   //    return;
   // }

   // use runHeader to get the run number for the database search
   RunHeader* runHeader = (RunHeader*) tree->GetUserInfo()->First();
   cout<< "runHeader->GetRun() = " << runHeader->GetRun() <<endl;
   time_t start_time = runHeader->GetTime();
   cout<< "run start time: " << std::ctime(&start_time);
   cout<< "program version is " << runHeader->GetVersion() <<endl;
   if (runHeader->GetTimeTag()) cout<< "event time tag was written out" <<endl;
   else cout<< "event time tag was not written out" <<endl;
   cout<<endl;

   // Geometry

   Geometry geometry(runHeader->GetRun(), dbname);   // get the run number from the RunHeader
   outputFiles->run = runHeader->GetRun();
   outputFiles->angle = geometry.angle_;
   outputFiles->open();

   Double_t uv[4];
   Double_t ut[4];
   for (int ilayer=0; ilayer<4; ++ilayer) {
      uv[ilayer] = geometry.uv_[ilayer];
      ut[ilayer] = geometry.ut_[ilayer];
      if (debug) cout<< "uv[" << ilayer << "] = " << uv[ilayer] << " ut[" << ilayer << "] = " << ut[ilayer] <<endl;
   }

   const RecoEvent* recoEvent = 0;
   tree->SetBranchAddress("revent", &recoEvent);
   // Book histo's 
   TProfile2D* hcal[5];
   for (int i=0; i<5; ++i) 
      hcal[i] = new TProfile2D(Form("hcal%d",i), Form("Stage %d response, MeV, t-v corrected",i),60,-180,180,18,-45,45,0,80);
   TH1F* hcalr[5];
   for (int i=0; i<5; ++i) hcalr[i] = new TH1F(Form("hcalr%d",i), Form("Stage %d row",i),400, 0,8000);  	
   TH1F* hcalc[5];
   for (int i=0; i<5; ++i) hcalc[i] = new TH1F(Form("hcalc%d",i), Form("Stage %d tv-corrected",i),400, 0,100);   
   TH1F *h7 = new TH1F("h7","Reconstructed WET",300,-10.5,289.5);
   // h7->SetStats(kFALSE);
   TProfile2D* hwepl= new TProfile2D("hwepl","Reconstructed WEPL",180,-6.35*30,6.35*30,45,-45,45,-10,180);	

   // ADC Pedestals
   //  Double_t ped[5] = {9.645, -20.484, -201.987, 62.966, -7.747};     // Celeste data
   Double_t ped[5] = {121.3, -71.5, -1137, 346.2, -49.};     // New pedestals (x6, reduced data)
   //   Prepare stuff for TV correction and convertion ADC->energy(MeV)
   Double_t ucal = 216.9 + 40; // approx position for the calorimeter entrance  
   Float_t par[5]; Float_t adc; Float_t Estage[5];
   CalF* f2cal[5]; 	//                   Declare calibration functions for 5 stages 
   std::fstream TVcalfile("TVcalib.txt", std::ios_base::in); //text file with TV-corr data   

   for (int i=0; i<5; ++i) {            // read fit parameters:
      TVcalfile >> par[0] >> par[1] >> par[2] >> par[3] >> par[4] ;
      f2cal[i]=new CalF(i,par);  // initialize t-v calibration  function  
   }
   TVcalfile.close(); 
   // Prepare stuff for WEPL calibration 	
   Float_t wepl_par[45];Float_t Wet;
   //open text file with WEPL calibration data (parameters of 9 pol4 curves)   
   std::fstream WEPLcalfile("wet5calibExp.txt", std::ios_base::in); 
   for (int i=0; i<45; ++i) WEPLcalfile >> wepl_par[i];
   WEPLcalfile.close(); 
   Wepl * WEPL=new Wepl(wepl_par); //initialize WEPL calibr. 
   WEPL->SetEthresholds(1,.99,1,1,1); // Set all stage thresholds to 1 MeV
   // Start loop over events

   if (event2 < event1) event2 = tree->GetEntries()-1;
   cout<< "Events to process: " << event2 - event1 + 1 << " out of total events " << tree->GetEntries() <<endl;

   EventOutput eventOutput;

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
         ) cout<< "---> processing entry " << ientry << "\t\toutputFiles->nevents = " << outputFiles->nevents <<endl; 

      for (int itrack=0; itrack<recoEvent->track->GetLast()+1; ++itrack) {
         const SuperTrack* superTrack = (const SuperTrack*) recoEvent->track->At(itrack);
         //	if(recoEvent->deltaT<14.) continue; 

         if (superTrack->angle > 0.07) continue;      // apply cut on the angle between the tracks in the supertrack
         for (int ical=0; ical<5; ++ical) {
            //            adc = recoEvent->SampleSum(ical, 1, 4, ped[ical]);
            adc = recoEvent->a[ical];
            adc=adc-ped[ical]; 
            hcalr[ical]->Fill(adc);  
            //	if (adc < 100) adc = 0;
            // Apply TV correction and convert ADC channel to E in MeV
            Estage[ical]=f2cal[ical]->CalTVf(superTrack->T(ucal),superTrack->V(ucal),adc);   
            // Fill TV-corr control histo's
            hcalc[ical]->Fill(Estage[ical]);          
            hcal[ical]->Fill(superTrack->T(ucal), superTrack->V(ucal),Estage[ical]);
         }
         // Get Wet=WEPL from Estage
         Wet=WEPL->EtoWEPL(Estage);
         if(Wet>999. || Wet<-999.) continue;
         // Fill WEPL-calib control histo's
         h7->Fill(Wet);  
         hwepl->Fill(superTrack->T(0), superTrack->V(0),Wet);

         // write output files

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

   if (plot) {
     // Plot control histo's
     gStyle->SetOptFit();
     for (int ical=0; ical<5; ++ical){
       new TCanvas;
       hcal[ical]->Draw("lego2");               
       new TCanvas;
       hcalc[ical]->Fit("gaus");;                
       new TCanvas;
       hcalr[ical]->Fit("gaus");;               
     } 
     new TCanvas;
     hwepl->Draw("lego2");               
     new TCanvas;
     h7->Draw();               
   }

   cout<< "Wrote " << outputFiles->nevents << " events into temporary files" <<endl;
}

void writeReco(const char* ifname, Int_t event1=0, Int_t event2=-1, const char* dbname="rundb-Sep2014.dat")
{
   Bool_t debug = kFALSE;
   if (debug) cout<< "debug in on" <<endl;

   TFile* ifile = new TFile(ifname);
   if (!ifile) {
      cout<< "Could not open file " << ifname <<endl;
      return;
   }

   TTree* tree = (TTree*) gDirectory->Get("r");
   if (!tree) {
      cout<< "Could not find tree r" <<endl;
      return;
   }

   OutputFiles outputFiles;

   writeReco(event1, event2, tree, &outputFiles, dbname);

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
   // const char* uniname = Form("%s-%0.0f.bin",ifname,outputFiles.angle);    // angle from outputFiles
   const char* uniname = Form("%s.bin",ifname);
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

   cout<< "Wrote " << ngood << " events into output binary file " << uniname <<endl;

   //outputFiles.close();
}
