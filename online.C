// Andriy Zatserklyaniy, April 17, 2014

#include "DataFormat.h"

#include <TROOT.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TObject.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2.h>
#include <TF1.h>
#include <TGraph.h>

#include <iostream>
#include <cassert>

using std::cout;     using std::endl;

class HTree: public TH1 {
private:
   Int_t Fill(Double_t, Double_t) {return 0;}
   Int_t Fill(const char*, Double_t) {return 0;}
public:
   TTree* tree;
   Float_t x;
public:
   HTree(const char* name="ht", const char* title="ht", const char* varname="x"): TH1() {
      if (name && *name) SetName(name);
      if (title && *title) SetTitle(title);
      tree = new TTree(Form("%s_tree",GetName()),Form("%s tree",GetName()));
      tree->Branch(varname, &x, Form("%s/F",GetName()));
      tree->SetDirectory(0);
      tree->SetEstimate(-1);           // use all entries to estimate the range

      gDirectory->Append(this);
   }
   virtual ~HTree() {
      tree->SetDirectory(0);
      delete tree;
   }
   void SetBranchName(const char* varname) {
      TBranch* branch = (TBranch*) tree->GetListOfBranches()->First();
      branch->SetName(varname);
   }
   Int_t Fill(Double_t value) {
      x = value;
      tree->Fill();
      return tree->GetEntries();
   }
   virtual void Draw(Option_t* option="") {
      if (!tree) return;
      if (tree->GetEntries() == 0) return;
      SetEstimate();
      TBranch* branch = (TBranch*) tree->GetListOfBranches()->First();
      if (!branch) return;
      const char* bname = branch->GetName();
      tree->Draw(bname, option);
      if (*GetTitle()) tree->GetHistogram()->SetTitle(GetTitle());
      if (*GetXaxis()->GetTitle()) tree->GetHistogram()->GetXaxis()->SetTitle(GetXaxis()->GetTitle());
      if (*GetYaxis()->GetTitle()) tree->GetHistogram()->GetYaxis()->SetTitle(GetYaxis()->GetTitle());
   }
   void SetLineColor(Color_t color) {tree->SetLineColor(color);}
   void SetFillColor(Color_t color) {tree->SetFillColor(color);}
   void SetFillStyle(Style_t style) {tree->SetFillStyle(style);}
   void SetEstimate() {tree->SetEstimate(tree->GetEntries());}
   ClassDef(HTree, 1);
};
ClassImp(HTree);

void EnergyDet(Int_t nfront=1, Int_t ntail=4, TTree* tree=0, Int_t event1=0, Int_t event2=-1)
{
   //-- Double_t ped[5] = {0, 0, 0, 0, 0};
   Double_t ped[5] = {
      9.645,
      -20.484,
      -201.987,
      62.966,
      -7.747
   };

   if (!tree) tree = (TTree*) gDirectory->Get("t");
   if (!tree) {
      cout<< "Could not find tree t" <<endl;
      return;
   }

   RunHeader* runHeader = (RunHeader*) tree->GetUserInfo()->First();
   cout<< "runHeader->GetRun() = " << runHeader->GetRun() <<endl;
   time_t start_time = runHeader->GetTime();
   //--g++ cout<< "run start time: " << std::ctime(&start_time);
   cout<< "run start time: " << ctime(&start_time);
   cout<< "program version is " << runHeader->GetVersion() <<endl;
   if (runHeader->GetTimeTag()) cout<< "event time tag was written out" <<endl;
   else cout<< "event time tag was not written out" <<endl;
   cout<<endl;

   PCTEvent* pCTEvent = 0;
   tree->SetBranchAddress("event", &pCTEvent);

   // HTree* ht_ch[5];
   // for (int channel=0; channel<5; ++channel) {
   //    ht_ch[channel] = new HTree(Form("ht_ch%d",channel), Form("Samples for channel %d",channel), "adc");
   // }

   TH1F* ht_ch[5];
   for (int channel=0; channel<5; ++channel) {
      ht_ch[channel] = new TH1F(Form("ht_ch%d",channel), Form("channel %d %s",channel,tree->GetTitle()), 1000,-100,10000);
   }

   if (event2 < event1) event2 = tree->GetEntries() - 1;

   cout<< "event1 = " << event1 << " event2 = " << event2 <<endl;

   for (int ientry=event1; ientry<=event2; ++ientry) {
      //pCTEvent->Clear();
      if (tree->LoadTree(ientry) < 0) {
         cout<< "LoadTree failed for entry " << ientry <<endl;
         break;
      }
      tree->GetEntry(ientry);

      if (event2-event1 < 100 || pCTEvent->evt%100000 == 0) cout<< "--> pCTEvent->evt = " << pCTEvent->evt << "\t event = " << ientry <<endl;

      for (int channel=0; channel<5; ++channel)
      {
         Int_t board = channel / 3;
         Int_t chan  = channel % 3;

         //--old-- if (pCTEvent->energyBoard[0].reduced == 0)
         if (pCTEvent->energyBoard[board].reduced == 0)           // should not make any difference
         {
            /// // look at the energy samples
            /// for (int iboard=0; iboard<2; ++iboard) {
            ///    // cout<< "pCTEvent->energyBoard[" << iboard << "].numSamples = " << pCTEvent->energyBoard[iboard].numSamples <<endl;
            ///    for (unsigned int isample=0; isample<pCTEvent->energyBoard[iboard].numSamples; ++isample) {
            ///       EnergySample* energySample = (EnergySample*) pCTEvent->energyBoard[iboard].samples->At(isample);
            ///       for (int ipulse=0; ipulse<3; ++ipulse) {
            ///          cout<< "pulse[" << ipulse << "] = " << energySample->pulse[ipulse] <<endl;
            ///       }
            ///    }
            /// }

            Int_t nsamples = pCTEvent->energyBoard[board].numSamples;
            //cout<< "nsamples = " << nsamples <<endl;

            Double_t x[16];         // actually 16 is the maximum possible number of samples
            Double_t y[16];
            for (int isample=0; isample<nsamples; ++isample) {
               x[isample] = isample;
               const EnergySample* sample = (const EnergySample*) pCTEvent->energyBoard[board].samples->At(isample);
               y[isample] = sample->pulse[chan];
            }
            // find a position of the maximum
            Int_t imax = 0;
            // assumes that the number of samples is 16
            Float_t sum = y[imax];
            for (int isample=0; isample<16; ++isample) {
               sum += y[isample];
               if (y[isample] > y[imax]) {
                  imax = isample;
               }
            }
            //-----------------------if (sum == 0) return;          // there are no samples in this event
            Int_t n1 = imax - nfront;
            if (n1 < 0) n1 = 0;
            Int_t n2 = imax + ntail;
            if (n2 > 15) n2 = 15;

            sum = 0;
            for (int isample=n1; isample<=n2; ++isample) sum += y[isample];
            sum -= (n2 - n1 + 1)*ped[channel];
            ht_ch[channel]->Fill(sum);

            //cout<< "channel = " << channel << " sum = " << sum <<endl;
         }
         else {
            ht_ch[channel]->Fill(pCTEvent->energyBoard[board].pulse[chan]);
         }
      }
   }

   cout<< "plot channels" <<endl;

   for (int channel=0; channel<5; ++channel) {
      new TCanvas;
      ht_ch[channel]->Draw();
      // ht_ch[channel]->tree->Draw("adc","adc>0&&adc<25000");
   }
}

void DataFormat(TTree* tree, Int_t event1=0, Int_t event2=-1)
{
   RunHeader* runHeader = (RunHeader*) tree->GetUserInfo()->First();
   cout<< "runHeader->GetRun() = " << runHeader->GetRun() <<endl;
   time_t start_time = runHeader->GetTime();
   //--g++ cout<< "run start time: " << std::ctime(&start_time);
   cout<< "run start time: " << ctime(&start_time);
   cout<< "program version is " << runHeader->GetVersion() <<endl;
   if (runHeader->GetTimeTag()) cout<< "event time tag was written out" <<endl;
   else cout<< "event time tag was not written out" <<endl;
   cout<<endl;

   PCTEvent* pCTEvent = 0;
   tree->SetBranchAddress("event", &pCTEvent);

   //TTree* otree = new TTree("tt", "Tracker tree");
   //otree->SetMarkerColor(2);
   //Tree::book(otree);

   // v-board
   HTree* ht_fpga0_address = new HTree("ht_fpga0_address", "Chip address for V-board, FPGA0", "address");
   HTree* ht_fpga0_nstrips = new HTree("ht_fpga0_nstrips", "The number of strips for V-board, FPGA#0", "nstrips");
   TH1F* h_v_fpga0_first = new TH1F("h_v_fpga0_first", "The first strip for V-board FPGA#0", 384, 0, 384);
   TH1F* h_v_fpga0_center = new TH1F("h_v_fpga0_center", "Cluster center of gravity for V-board FPGA#0", 384, 0, 384);

   // t-board
   HTree* ht_fpga5_address = new HTree("ht_fpga5_address", "Chip address for T-board, FPGA#5", "address");
   HTree* ht_fpga5_nstrips = new HTree("ht_fpga5_nstrips", "The number of strips for T-board, FPGA#5", "nstrips");
   TH1F* h_t_fpga5_first = new TH1F("h_t_fpga5_first", "The first strip for T-board FPGA#5", 768, 0, 768);
   TH1F* h_t_fpga5_center = new TH1F("h_t_fpga5_center", "Cluster center of gravity for T-board FPGA#5", 768, 0, 768);
   //
   HTree* ht_fpga7_address = new HTree("ht_fpga7_address", "Chip address for T-board, FPGA#7", "address");
   HTree* ht_fpga7_nstrips = new HTree("ht_fpga7_nstrips", "The number of strips for T-board, FPGA#7", "nstrips");
   TH1F* h_t_fpga7_first = new TH1F("h_t_fpga7_first", "The first strip for T-board FPGA#7", 768, 0, 768);
   TH1F* h_t_fpga7_center = new TH1F("h_t_fpga7_center", "Cluster center of gravity for T-board FPGA#7", 768, 0, 768);

   // energy detector
   HTree* ht_ch0 = new HTree("ht_ch0", "Energy detector channel 0", "pulse0");
   HTree* ht_ch1 = new HTree("ht_ch1", "Energy detector channel 1", "pulse1");
   HTree* ht_ch2 = new HTree("ht_ch2", "Energy detector channel 2", "pulse2");
   HTree* ht_ch3 = new HTree("ht_ch3", "Energy detector channel 3", "pulse3");
   HTree* ht_ch4 = new HTree("ht_ch4", "Energy detector channel 4", "pulse4");

   if (event2 < event1) event2 = tree->GetEntries() - 1;

   cout<< "event1 = " << event1 << " event2 = " << event2 <<endl;

   for (int ientry=event1; ientry<=event2; ++ientry) {
      //pCTEvent->Clear();
      if (tree->LoadTree(ientry) < 0) break;
      tree->GetEntry(ientry);

      if (event2-event1 < 100 || pCTEvent->evt%10000 == 0) cout<< "--> pCTEvent->evt = " << pCTEvent->evt << "\t event = " << ientry <<endl;

      //Tree::clear();
      cout<< "print out strips" <<endl;
      for (int iFPGA=0; iFPGA<12; ++iFPGA)
      {
         //Tree::fpga = iFPGA;
         TrackerFPGA& trackerFPGA = pCTEvent->trackerFPGA[iFPGA];
         if (trackerFPGA.numChips > 0) {
            for (unsigned ichip=0; ichip<trackerFPGA.numChips; ++ichip) {
               TrackerChip* trackerChip = (TrackerChip*) trackerFPGA.chips->At(ichip);
               cout<< "  -- trackerChip->address = " << (unsigned) trackerChip->address <<endl;
               //Tree::chip_address = trackerChip->address;
               for (int icluster=0; icluster<trackerChip->numClust; ++icluster) {
                  //Tree::nfirst = trackerChip->nfirst[icluster];
                  //Tree::nstrips = trackerChip->nstrips[icluster];
                  cout<< "  -- iFPGA = " << iFPGA << " nstrips[" << icluster << "] = " << (unsigned) trackerChip->nstrips[icluster] << " nfirst[" << icluster << "] = " << (unsigned) trackerChip->nfirst[icluster] <<endl;
               }
            }
         }
         //otree->Fill();
      }

      for (int iFPGA=0; iFPGA<12; ++iFPGA)
      {
         TrackerFPGA& trackerFPGA = pCTEvent->trackerFPGA[iFPGA];
         //cout<< iFPGA << "\ttrackerFPGA.numChips = " << trackerFPGA.numChips << " trackerFPGA.chips->GetLast()+1 = " << trackerFPGA.chips->GetLast()+1 <<endl;
         if (iFPGA == 0)
         {
            // upstream V-board
            if (trackerFPGA.numChips > 0) {
               for (unsigned int ichip=0; ichip<trackerFPGA.numChips; ++ichip)
               {
                  TrackerChip* trackerChip = (TrackerChip*) trackerFPGA.chips->At(ichip);
                  unsigned nfirst_abs[10];                  // absolute strip number 0..383    <-- ignore gap for now
                  unsigned nstrips[10];
                  unsigned ncluster = 0;
                  ht_fpga0_address->Fill(trackerChip->address);
                  for (int icluster=0; icluster<trackerChip->numClust; ++icluster) {
                     // cout<< "  -- iFPGA = " << iFPGA << " trackerChip->address = " << trackerChip->address << " nstrips[" << icluster << "] = " << trackerChip->nstrips[icluster] << " nfirst[" << icluster << "] = " << trackerChip->nfirst[icluster] <<endl;
                     ht_fpga0_nstrips->Fill(trackerChip->nstrips[icluster]);

                     assert(ncluster < 10);
                     nstrips[ncluster] = trackerChip->nstrips[icluster];
                     nfirst_abs[ncluster] = (5 - trackerChip->address)*64 + trackerChip->nfirst[icluster];  // V-board chip address 0..5
                     h_v_fpga0_first->Fill(nfirst_abs[ncluster]);
                     ++ncluster;
                  }

                  // done with this chip. Compute the center of gravity.
                  unsigned icurr = 0;
                  while (icurr < ncluster) {
                     unsigned first_strip = nfirst_abs[icurr];
                     unsigned last_strip = first_strip + nstrips[icurr];
                     unsigned next_cluster = icurr + 1;
                     // look at the next cluster if any
                     if (next_cluster < ncluster) {
                        if (first_strip <= last_strip) {
                           // merge clusters: use the last strip of the next cluster as the upper cluster boundary
                           last_strip = first_strip + nstrips[next_cluster];
                           // next cluster has been merged with icurr, increament the icurr to avoid processing of the merged cluster
                           ++icurr;
                        }
                     }
                     Double_t center = (first_strip + last_strip)/2.;
                     h_v_fpga0_center->Fill(center);

                     ++icurr;
                  }
               }
            }
         }

         if (iFPGA == 5)
         {
            // upstream T-board
            if (trackerFPGA.numChips > 0) {
               for (int ichip=0; ichip<trackerFPGA.numChips; ++ichip)
               {
                  TrackerChip* trackerChip = (TrackerChip*) trackerFPGA.chips->At(ichip);
                  unsigned nfirst_abs[10];                  // absolute strip number 0..383    <-- ignore gap for now
                  unsigned nstrips[10];
                  unsigned ncluster = 0;
                  ht_fpga5_address->Fill(trackerChip->address);
                  for (int icluster=0; icluster<trackerChip->numClust; ++icluster) {
                     // cout<< "  -- iFPGA = " << iFPGA << " trackerChip->address = " << trackerChip->address << " nstrips[" << icluster << "] = " << trackerChip->nstrips[icluster] << " nfirst[" << icluster << "] = " << trackerChip->nfirst[icluster] <<endl;
                     ht_fpga5_nstrips->Fill(trackerChip->nstrips[icluster]);

                     assert(ncluster < 10);
                     nstrips[ncluster] = trackerChip->nstrips[icluster];
                     nfirst_abs[ncluster] = (11 - trackerChip->address)*64 + trackerChip->nfirst[icluster];    // T-board chip address 0..11
                     h_t_fpga5_first->Fill(nfirst_abs[ncluster]);
                     ++ncluster;
                  }

                  // done with this chip. Compute the center of gravity.
                  unsigned icurr = 0;
                  while (icurr < ncluster) {
                     unsigned first_strip = nfirst_abs[icurr];
                     unsigned last_strip = first_strip + nstrips[icurr];
                     unsigned next_cluster = icurr + 1;
                     // look at the next cluster if any
                     if (next_cluster < ncluster) {
                        if (first_strip <= last_strip) {
                           // merge clusters: use the last strip of the next cluster as the upper cluster boundary
                           last_strip = first_strip + nstrips[next_cluster];
                           // next cluster has been merged with icurr, increament the icurr to avoid processing of the merged cluster
                           ++icurr;
                        }
                     }
                     Double_t center = (first_strip + last_strip)/2.;
                     h_t_fpga5_center->Fill(center);

                     ++icurr;
                  }
               }
            }
         }

         if (iFPGA == 7)
         {
            // upstream T-board
            if (trackerFPGA.numChips > 0) {
               for (int ichip=0; ichip<trackerFPGA.numChips; ++ichip)
               {
                  TrackerChip* trackerChip = (TrackerChip*) trackerFPGA.chips->At(ichip);
                  unsigned nfirst_abs[10];                  // absolute strip number 0..383    <-- ignore gap for now
                  unsigned nstrips[10];
                  unsigned ncluster = 0;
                  ht_fpga7_address->Fill(trackerChip->address);
                  for (int icluster=0; icluster<trackerChip->numClust; ++icluster) {
                     // cout<< "  -- iFPGA = " << iFPGA << " trackerChip->address = " << trackerChip->address << " nstrips[" << icluster << "] = " << trackerChip->nstrips[icluster] << " nfirst[" << icluster << "] = " << trackerChip->nfirst[icluster] <<endl;
                     ht_fpga7_nstrips->Fill(trackerChip->nstrips[icluster]);

                     assert(ncluster < 10);
                     nstrips[ncluster] = trackerChip->nstrips[icluster];
                     nfirst_abs[ncluster] = (11 - trackerChip->address)*64 + trackerChip->nfirst[icluster];    // T-board chip address 0..11
                     h_t_fpga7_first->Fill(nfirst_abs[ncluster]);
                     ++ncluster;
                  }

                  // done with this chip. Compute the center of gravity.
                  unsigned icurr = 0;
                  while (icurr < ncluster) {
                     unsigned first_strip = nfirst_abs[icurr];
                     unsigned last_strip = first_strip + nstrips[icurr];
                     unsigned next_cluster = icurr + 1;
                     // look at the next cluster if any
                     if (next_cluster < ncluster) {
                        if (first_strip <= last_strip) {
                           // merge clusters: use the last strip of the next cluster as the upper cluster boundary
                           last_strip = first_strip + nstrips[next_cluster];
                           // next cluster has been merged with icurr, increament the icurr to avoid processing of the merged cluster
                           ++icurr;
                        }
                     }
                     Double_t center = (first_strip + last_strip)/2.;
                     h_t_fpga7_center->Fill(center);

                     ++icurr;
                  }
               }
            }
         }
      }

      // for (int iboard=0; iboard<2; ++iboard) {
      //    EnergyBoard* energyBoard = (EnergyBoard*) pCTEvent->energyBoard[iboard];
      // }
      ht_ch0->Fill(pCTEvent->energyBoard[0].pulse[0]);
      ht_ch1->Fill(pCTEvent->energyBoard[0].pulse[1]);
      ht_ch2->Fill(pCTEvent->energyBoard[0].pulse[2]);
      ht_ch3->Fill(pCTEvent->energyBoard[1].pulse[0]);
      ht_ch4->Fill(pCTEvent->energyBoard[1].pulse[1]);

      // cout<< "pCTEvent->energyBoard[0].reduced = " << pCTEvent->energyBoard[0].reduced <<endl;
      if (pCTEvent->energyBoard[0].reduced == 0)
      {
         // look at the energy samples
         for (int iboard=0; iboard<2; ++iboard) {
            // cout<< "pCTEvent->energyBoard[" << iboard << "].numSamples = " << pCTEvent->energyBoard[iboard].numSamples <<endl;
            for (unsigned int isample=0; isample<pCTEvent->energyBoard[iboard].numSamples; ++isample) {
               EnergySample* energySample = (EnergySample*) pCTEvent->energyBoard[iboard].samples->At(isample);
               for (int ipulse=0; ipulse<3; ++ipulse) {
                  cout<< "pulse[" << ipulse << "] = " << energySample->pulse[ipulse] <<endl;
               }
            }
         }
      }
   }  // loop over entries

   new TCanvas;
   ht_fpga0_address->Draw();

   new TCanvas;
   h_v_fpga0_first->Draw();

   new TCanvas;
   h_v_fpga0_center->Draw();

   new TCanvas;
   ht_fpga5_address->Draw();

   new TCanvas;
   h_t_fpga5_first->Draw();

   new TCanvas;
   h_t_fpga5_center->Draw();

   new TCanvas;
   ht_fpga7_address->Draw();

   new TCanvas;
   h_t_fpga7_first->Draw();

   new TCanvas;
   h_t_fpga7_center->Draw();

   tree->ResetBranchAddresses();
}

void DataFormat(const char* ifname, Int_t event1=0, Int_t event2=-1)
{
   TFile* ifile = new TFile(ifname);
   if (!ifile) {
      cout<< "Could not open file " << ifname <<endl;
      return;
   }
   TTree* tree = (TTree*) ifile->Get("t");
   if (!tree) {
      cout<< "Could not find tree \"t\" in the file " << ifname <<endl;
      return;
   }

   DataFormat(tree, event1, event2);
}

void plotSamples(Int_t channel, Int_t event, TTree* tree=0, TCanvas* can=0)
{
   if (!tree) tree = (TTree*) gDirectory->Get("t");
   if (!tree) {
      cout<< "Could not find tree \"t\"" <<endl;
      return;
   }

   PCTEvent* pCTEvent = 0;
   tree->SetBranchAddress("event", &pCTEvent);

   Int_t nread = tree->GetEntry(event);
   if (nread == 0) {
      cout<< "Could not load entry " << event <<endl;
      return;
   }

   Int_t board = channel / 3;
   Int_t chan  = channel % 3;

   Int_t nsamples = pCTEvent->energyBoard[board].numSamples;
   cout<< "nsamples = " << nsamples <<endl;

   Double_t x[16];         // actually 16 is the maximum possible number of samples
   Double_t y[16];
   for (int isample=0; isample<nsamples; ++isample) {
      x[isample] = isample;
      const EnergySample* sample = (const EnergySample*) pCTEvent->energyBoard[board].samples->At(isample);
      y[isample] = sample->pulse[chan];
   }

   TGraph* gsample = new TGraph(nsamples, x, y);
   gsample->SetNameTitle(Form("gsample_ch%d_evt%d",channel,event), Form("Energy samples for channel %d, event %d;sample",channel,event));
   gsample->SetMarkerStyle(20);
   gsample->SetMarkerColor(2);
   
   if (can) can->cd();
   else new TCanvas;
   gsample->Draw("ap");
}

PCTEvent* GetPCTEvent(TTree* tree=0)
{
   if (!tree) tree = (TTree*) gDirectory->Get("t");
   if (!tree) {
      cout<< "Could not find tree \"t\"" <<endl;
      return 0;
   }

   PCTEvent* pCTEvent = new PCTEvent;
   tree->SetBranchAddress("event", &pCTEvent);

   return pCTEvent;
}

struct StripHit {
   //
   //          V-board chips
   //    6                       5
   //    11                      0
   //
   // V: 0..383  FPGA#3    383..767
   // T: 767---#11--0 767---#10---0
   // V: 0..383  FPGA#2   384..767 
   // T: 767---#9--0  767---#8----0
   //
   // T: 0---#6---767 0----#7----767
   // V: 767..384  FPGA#1    383..0
   // T: 0----#4--767 0----#5----767
   // V: 767..384  FPGA#0    383..0
   //
   //                ^
   //                |
   //                |  Particle Eye View
   //
   // event 125 good event
   // 125	 iFPGA = 0 trackerChip->address =  9 nstrips[0] = 1 nfirst[0] = 9
   // 125	 iFPGA = 1 trackerChip->address =  9 nstrips[0] = 1 nfirst[0] = 14
   // 125	 iFPGA = 2 trackerChip->address =  2 nstrips[0] = 1 nfirst[0] = 6
   // 125	 iFPGA = 3 trackerChip->address =  3 nstrips[0] = 1 nfirst[0] = 62
   // 125	 iFPGA = 4 trackerChip->address =  8 nstrips[0] = 1 nfirst[0] = 35
   // 125	 iFPGA = 6 trackerChip->address =  8 nstrips[0] = 1 nfirst[0] = 45
   // 125	 iFPGA = 9 trackerChip->address =  2 nstrips[0] = 1 nfirst[0] = 13
   // 125	 iFPGA = 11 trackerChip->address =  3 nstrips[0] = 1 nfirst[0] = 63
   //
   Int_t vstrip[4][2*384];    // particle eye view: layer 0: left part 384..768, right part: 0..383
   Int_t tstrip[4][4*384];
   Int_t* stripFPGA[12];
   StripHit() {
      // link v layer
      stripFPGA[0] = &vstrip[0][0];
      stripFPGA[1] = &vstrip[1][0];
      stripFPGA[2] = &vstrip[2][0];
      stripFPGA[3] = &vstrip[3][0];
      // link t layer
      stripFPGA[4] = &tstrip[0][0];
      stripFPGA[5] = &tstrip[0][2*384];
      stripFPGA[6] = &tstrip[1][0];
      stripFPGA[7] = &tstrip[1][2*384];
      stripFPGA[8] = &tstrip[2][0];
      stripFPGA[9] = &tstrip[2][2*384];
      stripFPGA[10] = &tstrip[3][0];
      stripFPGA[11] = &tstrip[3][2*384];

      clear();
   }
   void clear() {
      for (int iFPGA=0; iFPGA<12; ++iFPGA) for (int istrip=0; istrip<2*384; ++istrip) stripFPGA[iFPGA][istrip] = 0;
   }
   void Analyze() {
      cout<< "v-layer" <<endl;
      for (int ilayer=0; ilayer<4; ++ilayer) {
         for (int istrip=0; istrip<2*384; ++istrip) {
            if (vstrip[ilayer][istrip] > 0) cout<< "vstrip[" << ilayer << "][" << istrip << "] = " << vstrip[ilayer][istrip] <<endl;
         }
      }
      cout<< "t-layer" <<endl;
      for (int ilayer=0; ilayer<4; ++ilayer) {
         for (int istrip=0; istrip<4*384; ++istrip) {
            if (tstrip[ilayer][istrip] > 0) cout<< "tstrip[" << ilayer << "][" << istrip << "] = " << tstrip[ilayer][istrip] <<endl;
         }
      }
   }
};

class GeoHit: public StripHit {
public:
   Double_t vhit[4][2*384];
   Double_t thit[4][4*384];
   GeoHit(): StripHit() {
      clear();
   }
   void clear() {
      StripHit::clear();
      for (int ilayer=0; ilayer<4; ++ilayer) {
         for (int istrip=0; istrip<2*384; ++istrip) {
            vhit[ilayer][istrip] = 0;
            thit[ilayer][istrip] = 0;
         }
         for (int istrip=2*384; istrip<4*384; ++istrip) thit[ilayer][istrip] = 0;
      }
   }
};

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
   cout<< "run start time: " << ctime(&start_time);
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
   // Int_t nhitT = 0;
   std::vector<Int_t> hitV[4];
   std::vector<Int_t> hitT[4];

   // Double_t ut[4] = {-211.80, -161.80, 161.80, 211.80};
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

   // for (int iFPGA=0; iFPGA<12; ++iFPGA) {
   //    new TCanvas;
   //    hFPGA[iFPGA]->Draw();
   // }

   TCanvas* vcanvas = new TCanvas("vcanvas", "v FPGAs", 700,250);
   vcanvas->Divide(4,1);
   Int_t i_vcanvas = 0;
   for (int iFPGA=0; iFPGA<4; ++iFPGA) {
      vcanvas->cd(++i_vcanvas);
      hFPGA[iFPGA]->Draw();
   }
   vcanvas->cd();

   TCanvas* tcanvas = new TCanvas("tcanvas", "t FPGAs", 700,500);
   tcanvas->Divide(4,2);
   Int_t i_tcanvas = 0;
   for (int iFPGA=4; iFPGA<12; ++iFPGA) {
      tcanvas->cd(++i_tcanvas);
      hFPGA[iFPGA]->Draw();
   }
   tcanvas->cd();

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
   heffratioV->Fit("pol0", "", "", 40,360);

   tree->ResetBranchAddresses();
}
