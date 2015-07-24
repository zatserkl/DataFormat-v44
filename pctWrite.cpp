/*
Linux:
rootcint -f DataFormatReco_dict.cxx -c DataFormat.h Reco.h DataFormatReco_linkdef.h
g++ -Wall `$ROOTSYS/bin/root-config --cflags --glibs` pctWrite.cpp -o pctWrite DataFormatReco_dict.cxx

Mac OS X:
rootcint -f DataFormatReco_dict.cxx -c DataFormat.h Reco.h DataFormatReco_linkdef.h
clang++ -Wall `$ROOTSYS/bin/root-config --cflags --glibs` pctWrite.cpp -o pctWrite DataFormatReco_dict.cxx
*/

// T-V correction, ADC-MeV convertion and WEPL reconstruction implemented  
//  into Reco.C by Andriy Zatserklyaniy,  and renamed weplReco
// by V.A.Bashkirov (on May 9, 2014, modified May 16, 2014)
// In Root, assuming Andriy's classes are loaded, type .L weplReco.C+  and
// run CalibTVEW("<file name>.reco.root") to execute and see control histo's;
// e.q.  CalibTVEW("pCTraw_Run_705.out.root.reco.root") 

#include "DataFormat.h"
#include "BitBuffer.h"

#include "EventOutput.h"

#include "Reco.h"
#include "DataFormat.h"
#include "CalTV.h"
#include "Wepl.h"

#include <TFile.h>

//#include "stdafx.h"       // for Windows
#include <iostream>
#include <iomanip>
#include <bitset>
#include <vector>
#include <cassert>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdio>

// ClassImp(TrackerChip);
// ClassImp(TrackerFPGA);
// ClassImp(EnergySample);
// ClassImp(EnergyBoard);
// ClassImp(PCTEvent);
// ClassImp(RunHeader);

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
         vname[ilayer] << "temporary_file_v" << ilayer << "file-" << run << "-" << std::setfill('0') << std::setw(3) << angle << ".dat";
         tname[ilayer] << "temporary_file_t" << ilayer << "file-" << run << "-" << std::setfill('0') << std::setw(3) << angle << ".dat";
         uname[ilayer] << "temporary_file_u" << ilayer << "file-" << run << "-" << std::setfill('0') << std::setw(3) << angle << ".dat";

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

       //cout<< "OutputFiles: Wrote " << ngood << " good events into output binary file " << uniname <<endl;
       outputFiles.close();
   }
};  // struct OutputFiles

Int_t pctWrite(std::string inputFileName=""
              , Int_t event1=0, Int_t event2=-1
              , const char* dbname="rundb-May2015.dat"
              , const char* tv_calib_fname="TVcalib.txt"
              , const char* wet_calib_fname="wet5calibExp.txt"
              , OutputFiles* outputFiles=0
              , bool debug=false
             )
{
    //
    // TODO: beam spot position and pedestals are hardcoded into the code
    //

    bool do_wepl = true;
    if (!tv_calib_fname || !tv_calib_fname[0]) do_wepl = false;
    if (!wet_calib_fname || !wet_calib_fname[0]) do_wepl = false;
    if (!do_wepl) {
        cout<< "There are no TV-correction and/or WET calibration file provided. Nothing to do" <<endl;
        return 0;
    }

    //
    //   create BitBuffer object
    //
    BitBuffer bitBuffer(inputFileName.c_str(), 1024);

    //
    //   read the run header and create a RunHeader object
    //

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

    cout<< "runHeader.GetRun() = " << runHeader.GetRun() <<endl;
    time_t start_time = runHeader.GetTime();
    cout<< "run start time: " << std::ctime(&start_time);
    cout<< "program version is " << runHeader.GetVersion() <<endl;
    if (runHeader.GetTimeTag()) cout<< "event time tag was written out" <<endl;
    else cout<< "event time tag was not written out" <<endl;
    cout<< "runHeader.GetAngle() = " << runHeader.GetAngle() <<endl;
    cout<<endl;

    //
    // Use the run number from the runHeader to read the database
    //
    Geometry* geometry = new Geometry(runHeader.GetRun(), dbname);

    //
    // use angle from the geometry (== from rundb) or from the runHeader
    //
    //--from-rundb-- outputFiles->angle = geometry->angle_;     // take angle from the rundb
    outputFiles->angle = runHeader.GetAngle();                  // take angle from the RunHeader instead of rundb
    outputFiles->run = runHeader.GetRun();                      // take run from the RunHeader
    //
    //   open temporary output files. The angle will be embedded into the file name
    //
    outputFiles->open();

    PCTSensors* pCTSensors = new PCTSensors(geometry);

    //
    //  TODO: hardcoded beam spot
    //

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
    //----------------- prepare stuff for WEPL reconstruction -------------------
    //

    //
    //  TODO: hardcoded ADC Pedestals
    //

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

    EventOutput eventOutput;

    //
    //   DataFormat event object
    //
    PCTEvent* pCTEvent = new PCTEvent();
    //
    //   Reco event object
    //
    RecoEvent* recoEvent = new RecoEvent();

    //
    //   event loop over DataFormat tree
    //

    const unsigned char pattern[3] = {0xF0, 0x43, 0x54};                // 1p C T
    unsigned char buf[3];

    Int_t event = -1;
    try                                     // to catch EOF exception
    {
        while (bitBuffer.Inside())
        {
            //
            // search for beginning of event
            //
            if (debug) cout<< "search for beginning of the event" <<'\n';

            buf[0] = bitBuffer.fillbyte();
            buf[1] = bitBuffer.fillbyte();
            buf[2] = bitBuffer.fillbyte();
            bool found_event = true && buf[0] == pattern[0] && buf[1] == pattern[1] && buf[2] == pattern[2];

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

                found_event = true && buf[0] == pattern[0] && buf[1] == pattern[1] && buf[2] == pattern[2];

                if (buf[2] == 0) nzeros++;

                const Int_t nzeros_max = 3;
                if (nzeros > nzeros_max) {
                    // probably we see the end of file filled by zeros
                    cout<< "found " << nzeros_max+1 << " consecutive zero bytes: probably the rest of file is filled with zeros" <<endl;
                    break;
                }
            }

            if (!found_event) {
                if (debug) cout<< "Failed to find the next event " << event+1 <<'\n';
                break;
            }
            else if (debug) cout<< "found_event the next event " << event+1 <<'\n';

            ++event;

            if (event < event1) continue;
            if (event2 >= event1 && event > event2) break;

            if (event < 10 || (event < 100000 && event%10000 == 0) || (event%100000 == 0)) cout<< "--> processing event = " << event <<'\n';

            //
            //  read this event from the binary file
            //
            pCTEvent->Clear("C");
            loadEvent(bitBuffer, pCTEvent);

            //
            //  reconstruct the event
            //
            Reco reco(geometry, pCTSensors, &beamSpotIn, &beamSpotOut, pCTEvent, event, debug);
            reco.Tracking();

            if (debug) cout<< "call recoEvent->Extract(reco)" <<endl;
            recoEvent->Extract(reco);

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
                    //  pctWrite stuff: write event into temporary files
                    //
                    ++outputFiles->nevents;
                    eventOutput.good = kTRUE;
                    for (int ilayer=0; ilayer<4; ++ilayer) {
                        eventOutput.u[ilayer] = geometry->ut_[ilayer];
                        eventOutput.v[ilayer] = superTrack->V(eventOutput.u[ilayer]);
                        eventOutput.t[ilayer] = superTrack->T(eventOutput.u[ilayer]);
                    }
                    eventOutput.wepl = Wet;

                    for (int ilayer=0; ilayer<4; ilayer++) outputFiles->tfile[ilayer]->write((const char*) &eventOutput.t[ilayer], sizeof(Float_t));
                    for (int ilayer=0; ilayer<4; ilayer++) outputFiles->vfile[ilayer]->write((const char*) &eventOutput.v[ilayer], sizeof(Float_t));
                    for (int ilayer=0; ilayer<4; ilayer++) outputFiles->ufile[ilayer]->write((const char*) &eventOutput.u[ilayer], sizeof(Float_t));
                    outputFiles->wfile->write((const char*) &eventOutput.wepl, sizeof(Float_t));
                }
            }
        }
    }   // try block
    catch (const BitBufferException& bitBufferException) {
        cout<< bitBufferException.what() <<endl;
        bitBufferException.Print();
    }

    cout<< "Processed " << event-event1+1 << " events from the input file " << inputFileName <<endl;
    cout<< "Wrote " << outputFiles->nevents << " events into temporary files" <<endl;
    return event-event1+1;
}

Int_t pctWrite(const char* ifname
              , Int_t event1=0, Int_t event2=-1
              , const char* dbname="rundb-May2015.dat"
              , const char* tv_calib_fname="TVcalib.txt"
              , const char* wet_calib_fname="wet5calibExp.txt"
              , bool debug=false
             )
{
   OutputFiles outputFiles;

   Int_t ntot = pctWrite(ifname
             , event1, event2
             , dbname
             , tv_calib_fname
             , wet_calib_fname
             , &outputFiles
             , debug
            );

   // output file
   const char* uniname = Form("%s.bin",ifname);
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

    pctWrite(ifname
            , event1, event2
            , dbname
            , tv_calib_fname 
            , wet_calib_fname 
            , debug
           );

    return 0;
}
