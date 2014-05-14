// Andriy Zatserklyaniy, April 17, 2014

#include "Reco.h"

#include <TCanvas.h>
#include <TH2.h>
#include <TFile.h>

void plotReco(Int_t event1=0, Int_t event2=-1, TTree* tree=0)
{
   Bool_t debug = kFALSE;
   if (debug) cout<< "debug in on" <<endl;

   if (!tree) tree = (TTree*) gDirectory->Get("r");
   if (!tree) {
      cout<< "Could not find tree r" <<endl;
      return;
   }

   Double_t ped[5] = {9.645, -20.484, -201.987, 62.966, -7.747};     // Celeste data
   Double_t ucal = 216.9 + 40;                                       // approx position for the calorimeter entrance for runs >= 51

   const RecoEvent* recoEvent = 0;
   tree->SetBranchAddress("revent", &recoEvent);

   TH2F* hcal[5];
   // for (int i=0; i<5; ++i) hcal[i] = new TH2F(Form("hcal%d",i), Form("Cal response for channel %d",i), 200,-200,200, 50,-50,50);
   for (int i=0; i<5; ++i) hcal[i] = new TH2F(Form("hcal%d",i), Form("Cal response for channel %d",i), 400,-200,200, 10,-50,50);
   TH2F* hcal_i[5];
   // for (int i=0; i<5; ++i) hcal_i[i] = new TH2F(Form("hcal%d_i",i), Form("Nevents for Cal response for channel %d",i), 200,-200,200, 50,-50,50);
   for (int i=0; i<5; ++i) hcal_i[i] = new TH2F(Form("hcal%d_i",i), Form("Nevents for Cal response for channel %d",i), 400,-200,200, 10,-50,50);

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

      for (int itrack=0; itrack<recoEvent->track->GetLast()+1; ++itrack) {
         const SuperTrack* superTrack = (const SuperTrack*) recoEvent->track->At(itrack);

         if (superTrack->angle > 0.07) continue;      // apply cut on the angle between the tracks in the supertrack
         for (int ical=0; ical<5; ++ical) {
            //--old-- Float_t adc = recoEvent->a[ical];
            Float_t adc = recoEvent->SampleSum(ical, 1, 2, ped[ical]);
            if (adc < 100) adc = 0;
            // hcal[ical]->Fill(superTrack->tcal, superTrack->vcal, adc);
            // hcal_i[ical]->Fill(superTrack->tcal, superTrack->vcal);
            hcal[ical]->Fill(superTrack->T(ucal), superTrack->V(ucal), adc);
            hcal_i[ical]->Fill(superTrack->T(ucal), superTrack->V(ucal));
         }
      }
   }

   //new TCanvas;
   //hcal[0]->DrawCopy("colz");

   for (int ical=0; ical<5; ++ical) hcal[ical]->Divide(hcal_i[ical]);

   new TCanvas;
   hcal[0]->Draw("colz");
}

void plotReco(const char* ifname, Int_t event1=0, Int_t event2=-1)
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

   plotReco(event1, event2, tree);
}
