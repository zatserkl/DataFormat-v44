/*
   Windows. Make sure to setup VC environment with 
   vcvars32
   
   rootcint -f DataFormat_dict.cxx -c DataFormat.h DataFormat_linkdef.h
   cl -Z7 -MDd -GR -EHsc -MT -Ox pctroot.cpp DataFormat_dict.cxx -I C:\root\include /link -LIBPATH:C:\root\lib libCore.lib libCint.lib libGui.lib libRIO.lib libTree.lib libHist.lib libGraf.lib libMinuit.lib

   Mac OSX:
   rootcint -f DataFormat_dict.cxx -c DataFormat.h DataFormat_linkdef.h
   clang++ -Wall `$ROOTSYS/bin/root-config --cflags --glibs` pctroot.cpp -o pctroot DataFormat_dict.cxx

   Scientific Linux 6.5
   rootcint -f DataFormat_dict.cxx -c DataFormat.h DataFormat_linkdef.h
   g++ -Wall `$ROOTSYS/bin/root-config --cflags --glibs` pctroot.cpp -o pctroot DataFormat_dict.cxx

   for debugging
   clang++ -g -Wall `$ROOTSYS/bin/root-config --cflags --glibs` pctroot.cpp -o pctroot DataFormat_dict.cxx
 */

// Analyze a pCT unformatted raw data file.
// R. Johnson   June 8, 2013
// July 21, 2013:  revised for new data format changes and to work for all 14 FPGAs
// July 22, 2013:  fixed up so that it works when the system does not have all boards installed
// Andriy's modifications October 10, 2013 to do away with the binary string data
// November 21, 2013, modify to work with reduced events that have no pedestals stored
// Andriy Zatserklyaniy, April 17, 2014: production release

#include "DataFormat.h"
#include "BitBuffer.h"

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TObject.h"

//#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <bitset>
#include <vector>
#include <cassert>

using std::cout;      using std::endl;

ClassImp(TrackerChip);
ClassImp(TrackerFPGA);
ClassImp(EnergySample);
ClassImp(EnergyBoard);
ClassImp(PCTEvent);
ClassImp(RunHeader);

//Method to parse the binary string data and load a single event into the event data structure

int loadEvent(BitBuffer& bitBuffer, PCTEvent* pCTEvent)
{
  // cout<< "pCT_Event::loadEvent" <<'\n';

  // read 36-bit event time tag (in clock cycles)
  std::bitset<36> eventTimeTag;
  for (int i=35; i>=0; --i) eventTimeTag.set(i,bitBuffer.getbit());
  pCTEvent->timeTag = eventTimeTag.to_ulong();
  // cout<< "pCTEvent->timeTag = " << pCTEvent->timeTag <<endl;

  // read 12-bit time since the previous trigger
  std::bitset<12> DeltaT;
  for (int i=11; i>=0; --i) DeltaT.set(i,bitBuffer.getbit());
  pCTEvent->deltaT = DeltaT.to_ulong();

  // Parse the event header

  // read start bits 10
  Bool_t digit1 = bitBuffer.getbit();
  Bool_t digit2 = bitBuffer.getbit();
  if (digit1 != 1 || digit2 != 0) {
    cout<< "***Read error: " <<endl;
    if (digit1 != 1) cout<< "digit1 != 1" <<endl;
    if (digit2 != 0) cout<< "digit2 != 0" <<endl;

    bitBuffer.PrintBuffer();
    throw BitBufferException(&bitBuffer, "*** Error read start bit pattern 10 ***");
  }

  pCTEvent->addressError = bitBuffer.getbit();
  pCTEvent->tagError = bitBuffer.getbit();
  pCTEvent->CRCError = bitBuffer.getbit();
  bitBuffer.getbit(); // pCTEvent->TBDError = bitBuffer.getbit();

  std::bitset<18> evtNumber;
  for (int i=17; i>=0; --i) evtNumber.set(i,bitBuffer.getbit());
  pCTEvent->evt = evtNumber.to_ulong();	

  //
  // Loop over tracker layers
  //
  for (int iFPGA=0; iFPGA<12; iFPGA++)
  {
    //... 4-bit FPGA address
    pCTEvent->trackerFPGA[iFPGA].layerFPGA = bitBuffer.fillpart(4);

    //... 3-bit event tag
    pCTEvent->trackerFPGA[iFPGA].evtTag = bitBuffer.fillpart(3);
    /// cout << "Layer " << iFPGA << " evtTag=" << Event->Lyr[iFPGA].evtTag << "\n";

    //... 1-bit error flag (set in case of a trigger tag mismatch)
    Int_t tagError = bitBuffer.getbit();
    pCTEvent->trackerFPGA[iFPGA].tagError = tagError;

    //... 4-bit the number of chips reporting cluster data
    pCTEvent->trackerFPGA[iFPGA].numChips = bitBuffer.fillpart(4);

    //--debug-- cout<< "trackerFPGA[" << iFPGA << "].layerFPGA = " << pCTEvent->trackerFPGA[iFPGA].layerFPGA
    //--debug--   << " the number of chips with clusters trackerFPGA[" << iFPGA << "].numChips = " << pCTEvent->trackerFPGA[iFPGA].numChips
    //--debug-- <<endl;

    // Loop over hit chips
    for (int ichip=0; ichip<pCTEvent->trackerFPGA[iFPGA].numChips; ichip++)
    {
      // fired chips in readout order

      TrackerChip* trackerChip = (TrackerChip*) pCTEvent->trackerFPGA[iFPGA].chips->ConstructedAt(ichip);

      //... chip overflow bit
      trackerChip->overflow = bitBuffer.getbit();

      //... unused bit, set to 0 to differentiate this from a TOT packet
      bitBuffer.getbit();

      //... 4-bit number of clusters
      trackerChip->numClust = bitBuffer.fillpart(4);

      //... chip error bit
      trackerChip->error = bitBuffer.getbit();

      //... chip parity error bit
      trackerChip->errorParity = bitBuffer.getbit();

      //... 4-bit chip address
      trackerChip->address = bitBuffer.fillpart(4);

      //--debug-- cout<<"\t"<< ichip << "\t trackerChip->address = " << trackerChip->address << " numClust = " << trackerChip->numClust <<endl;

      // Loop over clusters
      for (int icluster=0; icluster<trackerChip->numClust; icluster++)
      {
        //... 6-bit number of strips minus 1 (0..10)
        trackerChip->nstrips[icluster] = bitBuffer.fillpart(6) + 1;    // assigned nstrips[icluster]

        //... 6-bit first strip (0..63)
        trackerChip->nfirst[icluster] = bitBuffer.fillpart(6);         // assigned nfirst[icluster]

        //--debug-- cout<<"\t\t"<< icluster << "\t trackerChip->nstrips[" << icluster << "] = " << (unsigned) trackerChip->nstrips[icluster] << " trackerChip->nfirst[icluster] = " << (unsigned) trackerChip->nfirst[icluster] <<endl;
      }
    }
  }  // End loop over tracker boards

  /// cout<< "print out strips" <<endl;
  /// for (int iFPGA=0; iFPGA<12; ++iFPGA) {
  ///   TrackerFPGA& trackerFPGA = pCTEvent->trackerFPGA[iFPGA];
  ///   if (trackerFPGA.numChips > 0) {
  ///     for (int ichip=0; ichip<trackerFPGA.numChips; ++ichip) {
  ///       TrackerChip* trackerChip = (TrackerChip*) trackerFPGA.chips->At(ichip);
  ///       cout<< "  -- trackerChip->address = " << (unsigned) trackerChip->address <<endl;
  ///       for (int icluster=0; icluster<trackerChip->numClust; ++icluster) {
  ///         cout<< "  -- iFPGA = " << iFPGA << " nstrips[" << icluster << "] = " << (unsigned) trackerChip->nstrips[icluster] << " nfirst[" << icluster << "] = " << (unsigned) trackerChip->nfirst[icluster] <<endl;
  ///       }
  ///     }
  ///   }
  /// }

  //
  // Loop over energy digitizer boards
  //
  for (int iboard=0; iboard<2; iboard++)
  {
    std::bitset<5> trgCntE;
    for (int i=4; i>=0; --i) trgCntE.set(i,bitBuffer.getbit());

    pCTEvent->energyBoard[iboard].reduced = bitBuffer.getbit();    // reduced: the last bit of the 6-bit Start Pattern

    if (pCTEvent->energyBoard[iboard].reduced)
    {
      pCTEvent->energyBoard[iboard].pedsOut = trgCntE.test(4);   // pedsOut: bit from the 6-bit Start Pattern
      trgCntE.set(4,0);                                           //--??? why we need to do that? -- to make next to_ulong enrgTag positive?
    }
    else pCTEvent->energyBoard[iboard].pedsOut = 0;

    pCTEvent->energyBoard[iboard].energyTag = trgCntE.to_ulong();

    pCTEvent->energyBoard[iboard].error = bitBuffer.getbit();

    //... some bits from Nsamp will be used for Nch below
    std::bitset<5> Nsamp;
    for (int i=4; i>=0; --i) Nsamp.set(i,bitBuffer.getbit());

    if (!Nsamp.any()) {
      pCTEvent->energyBoard[iboard].numChan = 0;
      pCTEvent->energyBoard[iboard].numSamples = 0;
      continue;   // No energy board was attached, so we see only an empty header
    }

    if (pCTEvent->energyBoard[iboard].reduced)
    {
      // set two bits of Nch from Nsamp
      std::bitset<2> Nch;
      for (int i=1; i>=0; --i) Nch.set(i,Nsamp[i+3]);          // sequence of this bit chunk: NCH,OTR2,OTR1,OTR0. Here NCH uses 2 bits
      pCTEvent->energyBoard[iboard].numChan = Nch.to_ulong();

      pCTEvent->energyBoard[iboard].numSamples = 0;
      pCTEvent->energyBoard[iboard].OTR[2] = Nsamp[2];
      pCTEvent->energyBoard[iboard].OTR[1] = Nsamp[1];
      pCTEvent->energyBoard[iboard].OTR[0] = Nsamp[0];

      //--old-- if (pCTEvent->energyBoard[iboard].numChan==2)  //deltaT is included for a board with reduced data only if there are just 2 channels
      //--old-- {
      //--old--   std::bitset<12> DeltaT;
      //--old--   for (int i=11; i>=0; --i) DeltaT.set(i,bitBuffer.getbit());
      //--old--   pCTEvent->deltaT = DeltaT.to_ulong();
      //--old-- }

      for (int ichan=0; ichan<pCTEvent->energyBoard[iboard].numChan; ichan++)
      {
        if (pCTEvent->energyBoard[iboard].pedsOut)
        {
          std::bitset<8> Ped;
          for (int i=7; i>=0; --i) Ped.set(i,bitBuffer.getbit());

          if (Ped.test(7)) {
            Ped.flip();
            pCTEvent->energyBoard[iboard].pedestal[ichan] = (-1) * (Ped.to_ulong() + 1);   //Negative number
          }
          else pCTEvent->energyBoard[iboard].pedestal[ichan] = Ped.to_ulong();
        }
        else pCTEvent->energyBoard[iboard].pedestal[ichan] = 0;

        std::bitset<16> Ph;
        for (int i=15; i>=0; --i) Ph.set(i,bitBuffer.getbit());

        if (Ph.test(15)) {                  //--AZ: to make positive: flip 0 <--> 1 and add 1
          Ph.flip();
          pCTEvent->energyBoard[iboard].pulse[ichan] = (-1) * (Ph.to_ulong() + 1);  //Negative number
        }
        else pCTEvent->energyBoard[iboard].pulse[ichan] = Ph.to_ulong();
      }

      if (pCTEvent->energyBoard[iboard].numChan == 2 && !pCTEvent->energyBoard[iboard].pedsOut) {
        for (int ibit=0; ibit<4; ++ibit) bitBuffer.getbit();  //Skip 4 unused bits
      }
    }
    else {
      // not reduced: pCTEvent->energyBoard[iboard].reduced is false here
      pCTEvent->energyBoard[iboard].numSamples = Nsamp.to_ulong();
      pCTEvent->energyBoard[iboard].numChan = 0;

      for (int i=0; i<3; i++) {
        pCTEvent->energyBoard[iboard].pedestal[i] = 0;
        pCTEvent->energyBoard[iboard].pulse[i] = 0;
        pCTEvent->energyBoard[iboard].OTR[i] = 0;
      }

      //--old-- std::bitset<12> DeltaT;
      //--old-- for (int i=11; i>=0; --i) DeltaT.set(i,bitBuffer.getbit());
      //--old-- pCTEvent->deltaT = DeltaT.to_ulong();

      for (int isample=0; isample<pCTEvent->energyBoard[iboard].numSamples; isample++)
      {
        EnergySample* energySample = (EnergySample*) pCTEvent->energyBoard[iboard].samples->ConstructedAt(isample);

        // read but ignore 3-bit sample count (this wraps around if the number of samples is greater than 7)
        //
        // std::bitset<3> Scnt;
        // for (int i=2; i>=0; --i) {
        //   Scnt.set(i,bitBuffer.getbit());
        // }
        bitBuffer.fillpart(3);

        // read 3 15-bit digitizations (each is 14 bits plus MSB=overflow indicator)
        for (int ich=0; ich<3; ich++)
        {
          // 1-bit overflow indicator
          energySample->OTR[ich] = bitBuffer.getbit();

          std::bitset<14> Samp;
          for (int i=13; i>=0; --i) Samp.set(i,bitBuffer.getbit());

          if (Samp.test(13)) {
            Samp.flip();
            energySample->pulse[ich] = (-1) * (Samp.to_ulong() + 1);   //Negative number
          }
          else energySample->pulse[ich] = Samp.to_ulong();
        }
      }
    }
  }  //End of loop over energy detector boards
  return 0;
}

//Simple analysis program to read in and print pCT events and construct an ntuple.
//R. Johnson  June 8, 2013
//July 21, 2013:  revised for new data format changes and to work for all 14 FPGAs

int main (int argc, char* argv[])
{
  for (int i=0; i<argc; i++) printf("%s ",argv[i]);
  printf("\n");
  if (argc<2) {
    cout<< "Usage: " << argv[0] << "filename [event1] [event2]" <<endl;
    return 0;
  }

  // read arguments 2 and 3 before the 1
  int event1 = 0;
  int event2 = -1;
  if (argc > 2) event1 = atoi(argv[2]);
  if (argc > 3) event2 = atoi(argv[3]);
  cout<< "event1 = " << event1 << " event2 = " << event2 <<endl;

  std::string inputFileName = argv[1];

  // BitBuffer bitBuffer(inputFileName.c_str(), 32);
  // BitBuffer bitBuffer(inputFileName.c_str(), 128);
  BitBuffer bitBuffer(inputFileName.c_str(), 1024);
  // BitBuffer bitBuffer(inputFileName.c_str(), 1048576);

  // Define Root stuff for the output ntuple
  Float_t PhCh0, PhCh1, PhCh2, PhCh3, PhCh4, PhCh5;
  Float_t PdCh0, PdCh1, PdCh2, PdCh3, PdCh4;            //--orig , PdCh5;
  Int_t numLyrs, numClust, evtNum, CRCerr, NhitLyr[12], numChipLyr[12];
  Float_t stripLyr[12], strip2Lyr[12], strip3Lyr[12], deltaT;

  TFile *F = new TFile(Form("%s-pCT_ntuple.root",inputFileName.c_str()),"recreate");
  TTree *T = new TTree("T","pct");
  T->Branch("PdCh0",&PdCh0,"PdCh0/F");				//Digitizer channel 0 pedestal
  T->Branch("PdCh1",&PdCh1,"PdCh1/F");				//Digitizer channel 1 pedestal
  T->Branch("PdCh2",&PdCh2,"PdCh2/F");				//Digitizer channel 2 pedestal
  T->Branch("PdCh3",&PdCh3,"PdCh3/F");				//Digitizer channel 3 pedestal
  T->Branch("PdCh4",&PdCh4,"PdCh4/F");				//Digitizer channel 4 pedestal
  T->Branch("PhCh0",&PhCh0,"PhCh0/F");				//Digitizer channel 0 pulse height
  T->Branch("PhCh1",&PhCh1,"PhCh1/F");				//Digitizer channel 1 pulse height
  T->Branch("PhCh2",&PhCh2,"PhCh2/F");				//Digitizer channel 2 pulse height
  T->Branch("PhCh3",&PhCh3,"PhCh3/F");				//Digitizer channel 3 pulse height
  T->Branch("PhCh4",&PhCh4,"PhCh4/F");				//Digitizer channel 4 pulse height
  T->Branch("PhCh5",&PhCh5,"PhCh5/F");				//Digitizer channel 5 pulse height (not normally present in pCT)
  T->Branch("numLyrs",&numLyrs,"numLyrs/I");			//Number of tracker layers hit
  T->Branch("numChipLyr0",&numChipLyr[0],"numChipLyr0/I");   //Number of chips in layer 0
  T->Branch("numChipLyr5",&numChipLyr[5],"numChipLyr5/I");   //Number of chips in layer 5
  T->Branch("numClust",&numClust,"numClust/I");		//Number of clusters total
  T->Branch("evtNum",&evtNum,"evtNum/I");				//Event number
  T->Branch("CRCerr",&CRCerr,"CRCerr/I");				//CRC error flag
  T->Branch("NhitLyr0",&NhitLyr[0],"NhitLyr0/I");
  T->Branch("NhitLyr1",&NhitLyr[1],"NhitLyr1/I");
  T->Branch("NhitLyr5",&NhitLyr[5],"NhitLyr5/I");
  T->Branch("stripLyr0",&stripLyr[0],"stripLyr0/F");	//Strip hit in fpga 0
  T->Branch("stripLyr1",&stripLyr[1],"stripLyr1/F");	//Strip hit in fpga 1
  T->Branch("stripLyr5",&stripLyr[5],"stripLyr5/F");	//Strip hit in fpga 5
  T->Branch("stripLyr7",&stripLyr[7],"stripLyr7/F");	//Strip hit in fpga 7
  T->Branch("strip2Lyr0",&strip2Lyr[0],"strip2Lyr0/F");	//Strip hit in fpga 0
  T->Branch("strip2Lyr1",&strip2Lyr[1],"strip2Lyr1/F");	//Strip hit in fpga 1
  T->Branch("strip2Lyr5",&strip2Lyr[5],"strip2Lyr5/F");	//Strip hit in fpga 5
  T->Branch("strip2Lyr7",&strip2Lyr[7],"strip2Lyr7/F");	//Strip hit in fpga 7
  T->Branch("strip3Lyr0",&strip3Lyr[0],"strip3Lyr0/F");	//Strip hit in fpga 0
  T->Branch("strip3Lyr1",&strip3Lyr[1],"strip3Lyr1/F");	//Strip hit in fpga 1
  T->Branch("strip3Lyr5",&strip3Lyr[5],"strip3Lyr5/F");	//Strip hit in fpga 5
  T->Branch("strip3Lyr7",&strip3Lyr[7],"strip3Lyr7/F");	//Strip hit in fpga 7
  T->Branch("deltaT",&deltaT,"deltaT/F");				//Time since previous trigger

  //--
  //--  ROOT class tree
  //--
  TFile* ofile = 0;
  if (event1 == 0 && event2 == -1) ofile = new TFile(Form("%s.root", inputFileName.c_str()), "recreate");
  else if (event2 == -1) ofile = new TFile(Form("%s-%d.root", inputFileName.c_str(),event1), "recreate");
  else ofile = new TFile(Form("%s-%d-%d.root", inputFileName.c_str(),event1,event2), "recreate");

  PCTEvent* pCTEvent = new PCTEvent();

  TTree* tree = new TTree("t", "pCT readout tree");
  tree->SetMarkerColor(2);
  tree->Branch("event", "PCTEvent", &pCTEvent);

  // Read run header

  const unsigned char pattern_run_header_identifier[3] = {0xD2, 0x55, 0x4E};     // 1R U N

  // run header identifier

  unsigned char buf_run_header_identifier[4];
  for (int ibyte=0; ibyte<3; ++ibyte) buf_run_header_identifier[ibyte] = bitBuffer.fillbyte();

  bool good_run_header = true
    && buf_run_header_identifier[0] == pattern_run_header_identifier[0]
    && buf_run_header_identifier[1] == pattern_run_header_identifier[1]
    && buf_run_header_identifier[2] == pattern_run_header_identifier[2]
  ;

  int ntry = 0;
  while (!good_run_header) {
    ++ntry;
    // cout<< "ntry = " << ntry <<endl;
    if (ntry > 10) break;
    buf_run_header_identifier[0] = buf_run_header_identifier[1];
    buf_run_header_identifier[1] = buf_run_header_identifier[2];
    buf_run_header_identifier[2] = bitBuffer.fillbyte();

    good_run_header = true
      && buf_run_header_identifier[0] == pattern_run_header_identifier[0]
      && buf_run_header_identifier[1] == pattern_run_header_identifier[1]
      && buf_run_header_identifier[2] == pattern_run_header_identifier[2]
    ;
  }

  if (!good_run_header) {
    cout<< "This is not a pCT binary data file" <<endl;
    return 0;
  }

  // read the run header

  RunHeader runHeader;
  char buf_run[3];
  for (int ibyte=0; ibyte<3; ++ibyte) buf_run[ibyte] = bitBuffer.fillbyte();
  char buf_time[4];
  for (int ibyte=0; ibyte<4; ++ibyte) buf_time[ibyte] = bitBuffer.fillbyte();
  char buf_timeTag = bitBuffer.fillbyte();
  char buf_version = bitBuffer.fillbyte();
  runHeader.Fill(buf_run, buf_time, buf_timeTag, buf_version);
  cout<< "runHeader: Run = " << runHeader.GetRun() << " timeTag = " << runHeader.GetTimeTag() << " version = " << runHeader.GetVersion() <<endl;
  time_t start_time = runHeader.GetTime();
  //--g++ cout<< "run start time: " << std::ctime(&start_time);
  cout<< "run start time: " << ctime(&start_time);

  tree->GetUserInfo()->Add(&runHeader);

  // Loop over events in the input stream

  int evtCnt = 0;
  int event = -1;                         // AZ event counter

  const unsigned char pattern[3] = {0xF0, 0x43, 0x54};                // 1p C T
  unsigned char buf[3];

  try                                     // to catch EOF exception
  {
    while (bitBuffer.Inside())
    {
      // search for beginning of event
      // cout<< "search for beginning of the event" <<'\n';

      buf[0] = bitBuffer.fillbyte();
      buf[1] = bitBuffer.fillbyte();
      buf[2] = bitBuffer.fillbyte();
      bool found_event = true
        && buf[0] == pattern[0]
        && buf[1] == pattern[1]
        && buf[2] == pattern[2]
      ;

      //-- if (!found_event) {
      //--   cout<< "found_event = false. Current buf: ";
      //--   for (int i=0; i<3; ++i) cout<< std::hex << std::hex << std::setfill('0') << std::setw(2) << (unsigned) buf[i] << std::dec << " ";
      //--   cout<<endl;
      //-- }

      Int_t nzeros = 0;
      while (!found_event)
      {
        //-- buf[0] = buf[1];
        //-- buf[1] = buf[2];
        //-- buf[2] = bitBuffer.fillbyte();

        buf[0] = (buf[0] << 4) | (buf[1] >> 4);
        buf[1] = (buf[1] << 4) | (buf[2] >> 4);
        buf[2] = (buf[2] << 4) | bitBuffer.fillpart(4);
        //-- cout<< "search for event header: " <<  bitBuffer <<endl;

        //-- cout<< "current buf: ";
        //-- for (int i=0; i<3; ++i) cout<< std::hex << std::hex << std::setfill('0') << std::setw(2) << (unsigned) buf[i] << std::dec << " ";
        //-- cout<<endl;

        found_event = true
          && buf[0] == pattern[0]
          && buf[1] == pattern[1]
          && buf[2] == pattern[2]
        ;

        if (buf[2] == 0) nzeros++;
        if (nzeros > 3) {
          // probably we see the end of file filled by zeros
          cout<< "found 4 consecutive zero bytes: probably the rest of file is filled with zeros" <<endl;
          break;
        }
      }

      if (!found_event) {
        // cout<< "Failed to find the next event" <<'\n';
        break;
      }
      // else cout<< "found_event the next event" <<'\n';

      ++event;

      if (event < event1) continue;
      if (event2 >= event1 && event > event2) break;

      if (event < 10 || (event < 100000 && event%10000 == 0) || (event%100000 == 0)) cout<< "------------------ event = " << event << " --------------------" <<'\n';
      //-- cout<< "\n------------------ event = " << event << " --------------------" <<'\n';

      pCTEvent->Clear("C");
      loadEvent(bitBuffer, pCTEvent);

      //--debug-- cout<< "---------- before tree->Fill -----------" <<endl;
      //--debug-- cout<< "print out strips" <<endl;
      //--debug-- for (int iFPGA=0; iFPGA<12; ++iFPGA) {
      //--debug--   TrackerFPGA& trackerFPGA = pCTEvent->trackerFPGA[iFPGA];
      //--debug--   if (trackerFPGA.numChips > 0) {
      //--debug--     for (int ichip=0; ichip<trackerFPGA.numChips; ++ichip) {
      //--debug--       TrackerChip* trackerChip = (TrackerChip*) trackerFPGA.chips->At(ichip);
      //--debug--       cout<< "  -- trackerChip->address = " << (unsigned) trackerChip->address <<endl;
      //--debug--       for (int icluster=0; icluster<trackerChip->numClust; ++icluster) {
      //--debug--         cout<< "  -- iFPGA = " << iFPGA << " nstrips[" << icluster << "] = " << (unsigned) trackerChip->nstrips[icluster] << " nfirst[" << icluster << "] = " << (unsigned) trackerChip->nfirst[icluster] <<endl;
      //--debug--       }
      //--debug--     }
      //--debug--   }
      //--debug-- }
      //--debug-- cout<< "---------- tree->Fill -----------" <<endl;
      tree->Fill();

      //      Event analysis:  fill up a Root ntuple
      ///// evtNum = thisEvent->Event->evtNum;
      evtNum = pCTEvent->evt;
      ///// if (evtCnt % 1000 == 0) cout << "Event count " << evtCnt << " Event number " << evtNum << "\n";
      ///// CRCerr = thisEvent->Event->CRCerr;
      CRCerr = pCTEvent->CRCError;
      ///// deltaT = thisEvent->Event->DeltaT;
      deltaT = pCTEvent->deltaT;
      int brd=0;
      ///// //--orig int enrgTag0= thisEvent->Event->Board[0].enrgTag;
      ///// //--orig int enrgTag1= thisEvent->Event->Board[1].enrgTag;
      ///// //		if (enrgTag0 != enrgTag1) cout << "enrg tag mismatch: " << enrgTag0 << " vs " << enrgTag1;
      if (pCTEvent->energyBoard[brd].numChan>0 || pCTEvent->energyBoard[brd].numSamples>0) {
        /////   //			if (thisEvent->Event->Board[brd].enrgTag != evtNum % 4) cout << "tag mismatch, energy tag=" << thisEvent->Event->Board[brd].enrgTag << "\n";
        if (pCTEvent->energyBoard[brd].reduced) {
          PhCh0 = pCTEvent->energyBoard[brd].pulse[0];  PdCh0 = pCTEvent->energyBoard[brd].pedestal[0];
          PhCh1 = pCTEvent->energyBoard[brd].pulse[1];	PdCh1 = pCTEvent->energyBoard[brd].pedestal[1];
          PhCh2 =	pCTEvent->energyBoard[brd].pulse[2];	PdCh2 = pCTEvent->energyBoard[brd].pedestal[2];
        } else {
          PhCh0 = 0.; PhCh1 = 0.; PhCh2 = 0.;
          /// enrgSamp *thisSamp = thisEvent->Event->Board[brd].firstSample;
          //EnergySample* energySample = (EnergySample*) pCTEvent->energyBoard[brd].samples->At(0);
          //if (thisSamp != 0)
          if (pCTEvent->energyBoard[brd].samples->GetLast()+1 > 0)
          {
            EnergySample* energySample0 = (EnergySample*) pCTEvent->energyBoard[brd].samples->At(0);
            int ped0 = energySample0->pulse[0]; PdCh0 = ped0;
            int ped1 = energySample0->pulse[1]; PdCh1 = ped1;
            int ped2 = energySample0->pulse[2]; PdCh2 = ped2;
            //while (thisSamp != 0)
            for (int isample=1; isample<pCTEvent->energyBoard[brd].samples->GetLast()+1; ++isample)
            {
              EnergySample* energySample = (EnergySample*) pCTEvent->energyBoard[brd].samples->At(isample);
              PhCh0 = PhCh0 + energySample->pulse[0] - ped0;
              PhCh1 = PhCh1 + energySample->pulse[1] - ped1;
              PhCh2 = PhCh2 + energySample->pulse[2] - ped2;
            }
          }
        }
      } else {
        PhCh0 = 0; PdCh0 = 0;
        PhCh1 = 0; PdCh1 = 0;
        PhCh2 = 0; PdCh2 = 0;
      }
      brd=1;
      if (pCTEvent->energyBoard[brd].numChan>0 || pCTEvent->energyBoard[brd].numSamples>0) {
        /////   //			if (thisEvent->Event->Board[brd].enrgTag != evtNum % 4) cout << "tag mismatch, energy tag=" << thisEvent->Event->Board[brd].enrgTag << "\n";
        if (pCTEvent->energyBoard[brd].reduced) {
          PhCh3 = pCTEvent->energyBoard[brd].pulse[0];  PdCh3 = pCTEvent->energyBoard[brd].pedestal[0];
          PhCh4 = pCTEvent->energyBoard[brd].pulse[1];	PdCh4 = pCTEvent->energyBoard[brd].pedestal[1];
          PhCh5 =	pCTEvent->energyBoard[brd].pulse[2];
        } else {
          PhCh3 = 0.; PhCh4 = 0.; PhCh5 = 0.;
          int ped3 = 0; int ped4 = 0; int ped5 = 0;
          /// enrgSamp *thisSamp = thisEvent->Event->Board[brd].firstSample;
          //EnergySample* energySample = (EnergySample*) pCTEvent->energyBoard[brd].samples->At(0);
          //if (thisSamp != 0)
          if (pCTEvent->energyBoard[brd].samples->GetLast()+1 > 0)
          {
            EnergySample* energySample0 = (EnergySample*) pCTEvent->energyBoard[brd].samples->At(0);
            ped3 = energySample0->pulse[0]; PdCh3 = ped3;
            ped4 = energySample0->pulse[1]; PdCh4 = ped4;
            ped5 = energySample0->pulse[2];
            //while (thisSamp != 0)
            for (int isample=1; isample<pCTEvent->energyBoard[brd].samples->GetLast()+1; ++isample)
            {
              EnergySample* energySample = (EnergySample*) pCTEvent->energyBoard[brd].samples->At(isample);
              PhCh3 = PhCh3 + energySample->pulse[0] - ped3;
              PhCh4 = PhCh4 + energySample->pulse[1] - ped4;
              PhCh5 = PhCh5 + energySample->pulse[2] - ped5;
            }
          }
        }
      } else {
        PhCh3 = 0; PdCh3 = 0;
        PhCh4 = 0; PdCh4 = 0;
        PhCh5 = 0;
      }

      numLyrs = 0;
      numClust = 0;
      for (int lyr=0; lyr<12; lyr++) {
        NhitLyr[lyr] = 0;
        stripLyr[lyr] = -999.;
        strip2Lyr[lyr] = -999.;
        numChipLyr[lyr] = pCTEvent->trackerFPGA[lyr].numChips;
        if (pCTEvent->trackerFPGA[lyr].numChips > 0) {
          //--orig int lyrTag= thisEvent->Event->Lyr[lyr].evtTag % 4;
          //				if (lyrTag != enrgTag0) cout << "Tkr tag mismatch, layer " << lyr << " tag=" << lyrTag << "\n";
          //				if (thisEvent->Event->Lyr[lyr].evtTag != evtNum % 8) cout << "Tkr tag mismatch, evtTag=" << thisEvent->Event->Lyr[lyr].evtTag << "\n";
          numLyrs++;
          //tkrChip *thisChip = thisEvent->Event->Lyr[lyr].firstChip;
          int stripOne1 = -999;
          int stripOneLen;
          int stripOneChip;
          int stripTwo1 = -999;
          int stripTwoLen = -999;
          int stripTwoChip = -999;
          //while (thisChip != 0)
          for (int ichip=0; ichip<pCTEvent->trackerFPGA[lyr].numChips; ++ichip)
          {
            TrackerChip* trackerChip = (TrackerChip*) pCTEvent->trackerFPGA[lyr].chips->At(ichip);
            numClust = numClust + trackerChip->numClust;
            if (numClust > 0) {
              //tkrClust *thisClust = thisChip->firstClust;
              //while (thisClust != 0)
              for (int icluster=0; icluster<numClust; ++icluster)
              {
                NhitLyr[lyr]++;
                if (NhitLyr[lyr]==1) {
                  stripLyr[lyr] = 64.0*(trackerChip->address) + 63.0-(trackerChip->nfirst[icluster] + (trackerChip->nstrips[icluster]-1.0)/2.);
                  stripOne1 = trackerChip->nfirst[icluster];
                  stripOneLen = trackerChip->nstrips[icluster];
                  stripOneChip = trackerChip->address;
                } 
                if (NhitLyr[lyr]==2) {
                  if (trackerChip->nfirst[icluster] + trackerChip->nstrips[icluster]-1 == 63 && stripOne1 == 0 && trackerChip->address == stripOneChip + 1) {   //Cluster continued from previous chip 
                    stripLyr[lyr] = 64.0*(trackerChip->address) + 63.0-(trackerChip->nfirst[icluster] + (trackerChip->nstrips[icluster]+stripOneLen-1.0)/2.);
                    //									cout << "1: Chip " << stripOneChip << " length= " << stripOneLen << " strip " << stripOne1 << "\n";
                  }
                  else {
                    strip2Lyr[lyr] = 64.0*(trackerChip->address) + 63.0-(trackerChip->nfirst[icluster] + (trackerChip->nstrips[icluster]-1.0)/2.);
                    stripTwo1 = trackerChip->nfirst[icluster];
                    stripTwoLen = trackerChip->nstrips[icluster];
                    stripTwoChip = trackerChip->address;
                  }
                }
                if (NhitLyr[lyr]==3) {
                  if (trackerChip->nfirst[icluster] + trackerChip->nstrips[icluster] - 1 == 63 && stripTwo1 == 0 && trackerChip->address == stripTwoChip + 1) {   //Cluster continued from previous chip 
                    strip2Lyr[lyr] = 64.0*(trackerChip->address) + 63.0-(trackerChip->nfirst[icluster] + (trackerChip->nstrips[icluster]+stripTwoLen-1.0)/2.);
                    //									cout << "2: Chip " << stripTwoChip << " length= " << stripTwoLen << " strip " << stripTwo1 << "\n";
                  }
                  else {
                    strip3Lyr[lyr] = 64.0*(trackerChip->address) + 63.0-(trackerChip->nfirst[icluster] + (trackerChip->nstrips[icluster]-1.0)/2.);
                  }
                }
                //							printf("stripLyr=%f  chip=%d   First strip=%2d,  Number of strips=%2d \n",stripLyr[lyr],thisChip->Address,thisClust->firstStrip,thisClust->numStrips);
                //--thisClust = thisClust->nextClust;
              }
            }
            //--thisChip = thisChip->nextChip;
          }
        }
      }

      T->Fill();
      evtCnt++;
    } // while block
  } // try block
  catch (const BitBufferException& bitBufferException) {
    cout<< bitBufferException.what() <<endl;
    bitBufferException.Print();
  }

  cout << "\n" << evtCnt << " events processed.\n"; 
  cout << "Wrote " << T->GetEntries() << " events to the ntuple to the file " << F->GetName() << "\n";
  // T->Print();
  F->Write();
  F->Close();

  cout<< "Writing " << tree->GetEntries() << " events into output file " << ofile->GetName() <<endl;
  ofile->Write();
  ofile->Close();

  delete pCTEvent;    pCTEvent = 0;
  return 0;
}
