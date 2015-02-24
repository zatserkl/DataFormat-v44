/*
   Windows. Make sure to setup VC environment with 
   vcvars32
   
   rootcint -f DataFormat_dict.cxx -c DataFormat.h DataFormat_linkdef.h
   cl -Z7 -MDd -GR -EHsc -MT -Ox pctroot.cpp DataFormat_dict.cxx -I C:\root\include /link -LIBPATH:C:\root\lib libCore.lib libCint.lib libGui.lib libRIO.lib libTree.lib libHist.lib libGraf.lib libMinuit.lib

   Mac OSX:
   ROOT 5:
   rootcint -f DataFormat_dict.cxx -c DataFormat.h DataFormat_linkdef.h
   clang++ -Wall `$ROOTSYS/bin/root-config --cflags --glibs` pctroot.cpp -o pctroot DataFormat_dict.cxx
   ROOT 6:
   rootcling -f DataFormat_dict.cxx -c DataFormat.h
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
  //--std-- cout<< "loadEvent: eventTimeTag[0] = " << eventTimeTag[0] << " loadEvent: eventTimeTag[35] = " << eventTimeTag[35] << " pCTEvent->timeTag = " << pCTEvent->timeTag <<endl;

  //--timePix-- std::bitset<1> startBit;
  //--timePix-- startBit.set(0,bitBuffer.getbit());
  //--timePix-- pCTEvent->startBit = startBit.to_ulong();
  //--timePix-- std::bitset<35> eventTimeTag;
  //--timePix-- for (int i=34; i>=0; --i) eventTimeTag.set(i,bitBuffer.getbit());
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
      Int_t numClust_max = trackerChip->numClust;
      if (numClust_max > 10) {
        cout<< "--> reduce the numClust from " << trackerChip->numClust << " to 10" <<endl;
        numClust_max = 10;
      }
      for (int icluster=0; icluster<numClust_max; icluster++)
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

  BitBuffer bitBuffer(inputFileName.c_str(), 1024);

  //--
  //--  ROOT class tree
  //--
  TFile* ofile = 0;
  if (event1 == 0 && event2 == -1) ofile = new TFile(Form("%s.root", inputFileName.c_str()), "recreate");
  else if (event2 == -1) ofile = new TFile(Form("%s-%d.root", inputFileName.c_str(),event1), "recreate");
  else ofile = new TFile(Form("%s-%d-%d.root", inputFileName.c_str(),event1,event2), "recreate");

  PCTEvent* pCTEvent = new PCTEvent();

  TTree* tree = new TTree("t", inputFileName.c_str());
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

  RunHeader runHeader;    // create RunHeader object to store in the output tree the info from the run header
  char buf_run[3];
  for (int ibyte=0; ibyte<3; ++ibyte) buf_run[ibyte] = bitBuffer.fillbyte();
  char buf_time[4];
  for (int ibyte=0; ibyte<4; ++ibyte) buf_time[ibyte] = bitBuffer.fillbyte();
  char buf_timeTag = bitBuffer.fillbyte();
  char buf_version = bitBuffer.fillbyte();

  int angle10 = 0;
  UChar_t version = buf_version;
  if (version > 61)                 // version 61 was used for the Mar2014 beam test
  {
    // the 12-bit stage angle, in tens of degree was introduced from May 2014
    std::bitset<12> angle_msb;
    for (int i=11; i>=0; --i) angle_msb.set(i,bitBuffer.getbit());
    angle10 = angle_msb.to_ulong();
  }

  runHeader.Fill(buf_run, buf_time, buf_timeTag, buf_version, angle10);
  cout<< "runHeader: Run = " << runHeader.GetRun() << " timeTag = " << runHeader.GetTimeTag() << " version = " << runHeader.GetVersion() <<endl;
  time_t start_time = runHeader.GetTime();
  //--g++ cout<< "run start time: " << std::ctime(&start_time);
  cout<< "run start time: " << ctime(&start_time);

  cout<< "runHeader.GetAngle() = " << runHeader.GetAngle() <<endl;

  tree->GetUserInfo()->Add(&runHeader);

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

        const Int_t nzeros_max = 3;
        if (nzeros > nzeros_max) {
          // probably we see the end of file filled by zeros
          cout<< "found " << nzeros_max+1 << " consecutive zero bytes: probably the rest of file is filled with zeros" <<endl;
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

    } // while block
  } // try block
  catch (const BitBufferException& bitBufferException) {
    cout<< bitBufferException.what() <<endl;
    bitBufferException.Print();
  }

  cout<< "Writing " << tree->GetEntries() << " events into output file " << ofile->GetName() <<endl;
  ofile->Write();
  ofile->Close();

  delete pCTEvent;    pCTEvent = 0;
  return 0;
}
