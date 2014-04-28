// Andriy Zatserklyaniy, April 17, 2014

#include "DataFormat.h"

#include <TROOT.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TObject.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2.h>
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
   ClassDef(HTree, 1);
};
ClassImp(HTree);

// ClassImp(TrackerChip);
// ClassImp(TrackerFPGA);
// ClassImp(EnergySample);
// ClassImp(EnergyBoard);
// ClassImp(PCTEvent);
// 
// // 1440 = room for 10 elements for each of 12 chip in each of 12 tracker FPGA
// UInt_t TrackerChipBuffer::nfirst[1440];
// UInt_t TrackerChipBuffer::nstrips[1440];
// UInt_t TrackerChipBuffer::ncurrent = 0;

//namespace Tree
//{
//   Int_t fpga;
//   Int_t chip_address;
//   Int_t nfirst;
//   Int_t nstrips;
//   Double_t center;
//
//   void clear() {
//      fpga = -1;
//      chip_address = -1;
//      nfirst = -1;
//      nstrips = 0;
//      center = -100;
//   }
//
//   void book(TTree* tree) {
//      tree->Branch("fpga",    &fpga,   "fpga/I");
//      tree->Branch("chip_address",    &chip_address,   "chip_address/I");
//      tree->Branch("nfirst",  &nfirst, "nfirst /I");
//      tree->Branch("nstrips",  &nstrips, "nstrips /I");
//      tree->Branch("center",    &center,   "center/D");
//   }
//
//   void connect(TTree* tree) {
//      tree->SetBranchAddress("fpga",    &fpga);
//      tree->SetBranchAddress("chip_address",    &chip_address);
//      tree->SetBranchAddress("nfirst",  &nfirst);
//      tree->SetBranchAddress("nstrips",  &nstrips);
//      tree->SetBranchAddress("center",    &center);
//   }
//}  // namespace Tree;

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
