// Andriy Zatserklyaniy, April 17, 2014

#ifndef DataFormat_h
#define DataFormat_h

#include <TROOT.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TObject.h>
#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <iomanip>

using std::cout;     using std::endl;

/*
        Example of use with TClonesArray::ConstructedAt

int loadEvent(BitBuffer& bitBuffer, PCTEvent* pCTEvent)
{
  ...
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

    // Loop over hit chips
    for (int ichip=0; ichip<pCTEvent->trackerFPGA[iFPGA].numChips; ichip++)
    {
      // fired chips in readout order

      TrackerChip* trackerChip = (TrackerChip*) pCTEvent->trackerFPGA[iFPGA].chips->ConstructedAt(ichip);

      //... chip overflow bit
      trackerChip->overflow = bitBuffer.getbit();
      ...
    }
    ...
  }
  ...
}
*/

class TrackerChip: public TObject {
public:
  UChar_t address;         //[0,11,4]
  Bool_t overflow;         //[0,1,1]
  Bool_t error;            //[0,1,1]
  Bool_t errorParity;      //[0,1,1]
  Int_t numClust;          //[0,10,6]
  UChar_t* nfirst;         //[numClust][0,63,6]
  UChar_t* nstrips;        //[numClust][0,63,6]
  TrackerChip(): TObject()
     , address(0)
     , overflow(false)
     , error(false)
     , errorParity(false)
     , numClust(0)
   {
      // cout<< "constructor TrackerChip" <<endl;
      nfirst = new UChar_t[10];   // the max number of clusters is 32 but truncated to max 10 clusters
      nstrips = new UChar_t[10];
   }
  ~TrackerChip() {
     //cout<< "destructor ~TrackerChip" <<endl;
     delete[] nfirst;
     delete[] nstrips;
  }
   void Clear(Option_t*) {
      //cout<< "TrackerChip::Clear" <<endl;
      address = 0;
      overflow = false;
      error = false;
      errorParity = false;
      numClust = 0;
   }
   ClassDef(TrackerChip, 17);
};

class TrackerFPGA: public TObject {
public:
  UChar_t layerFPGA;         //[0,12,4]      // 4-bit tracker FPGA address
  Bool_t evtTag;            //[0,1,1]
  Bool_t tagError;          //[0,1,1]
  UChar_t numChips;          //[0,11,6]
  TClonesArray* chips;      //->    maximum is 12 chips
  TrackerFPGA(): TObject()
    , layerFPGA(0)
    , evtTag(0)
    , tagError(0)
    , numChips(0)
  {
    chips = new TClonesArray("TrackerChip",12);
    Clear(0);
  }
  ~TrackerFPGA() {
    //cout<< "destructor ~TrackerFPGA" <<endl;
    delete chips;    chips = 0;
  }
  void Clear(Option_t*) {
     //cout<< "TrackerFPGA::Clear" <<endl;
     layerFPGA = 0;
     evtTag = 0;
     tagError = 0;
     numChips = 0;
     chips->Clear("C");    // NB: Clear("C") to call TrackerChip::Clear to delete cluster data (basically, to set numClust = 0)
  }
  void Print(Option_t* option="") const {
     if (!option) return;
     for (int ichip=0; ichip<chips->GetLast()+1; ++ichip) {
        const TrackerChip* chip = (const TrackerChip*) chips->At(ichip);
        for (int icluster=0; icluster<chip->numClust; ++icluster) {
           bool errors = chip->overflow || chip->error || chip->errorParity;
           cout<< "FPGA #" << (unsigned) layerFPGA
              << " chip " << (unsigned) chip->address
              << " nfirst[" << icluster << "] = " << (unsigned) chip->nfirst[icluster]
              << " nstrips[" << icluster << "] = " << (unsigned) chip->nstrips[icluster]
              << " errors = " << errors
           <<endl;
        }
     }
  }
  ClassDef(TrackerFPGA, 7);
};

class EnergySample: public TObject {
public:
  Bool_t OTR[3];                 //[0,1,1] out-of-range
  Short_t pulse[3];                //[0,0,16]
public:
  EnergySample(): TObject() {
     //cout<< "constructor EnergySample" <<endl;
     for (int i=0; i<3; ++i) {
        OTR[i] = kFALSE;
        pulse[i] = 0;
     }
  }
  ClassDef(EnergySample, 8);
};

class EnergyBoard: public TObject {
public:
   UChar_t energyTag;           //[0,7,3]
   Bool_t reduced;            //[0,1,1]
   Bool_t pedsOut;            //[0,1,1]
   Bool_t error;              //[0,1,1]
   UChar_t numChan;             //[0,3,2]
   Bool_t OTR[3];             //[0,1,1]
   Short_t pedestal[3];         //[0,0,8]
   Short_t pulse[3];            //[0,0,16]
   UChar_t numSamples;          //[0,0,5]
   TClonesArray* samples;       //->
public:
   EnergyBoard(): TObject()
      , energyTag(0)
      , reduced(false)
      , error(false)
   {
      samples = new TClonesArray("EnergySample", 16);
      Clear(0);
   }
   ~EnergyBoard() {
      //cout<< "destructor ~EnergyBoard" <<endl;
      delete samples;
   }
   void Clear(Option_t*) {
      //cout<< "EnergyBoard::Clear" <<endl;
      energyTag = 0;
      reduced = false;
      error = false;
      numChan = 0;
      numSamples = 0;
      samples->Clear();    // NB: just Clear(), not Clear("C"), don't need to call EnergySample::Clear()
      for (int i=0; i<3; ++i) {
         OTR[i] = 0;
         pedestal[i] = 0;
         pulse[i] = 0;
      }
   }
   ClassDef(EnergyBoard, 9);
};

class PCTEvent: public TObject {
public:
   Long64_t timeTag;                   // 36-bit event time tag
   UInt_t evt;
   Bool_t addressError;                //[0,1,1]
   Bool_t tagError;                    //[0,1,1]
   Bool_t CRCError;                    //[0,1,1]
   Short_t deltaT;                     //[0,0,16]
   TrackerFPGA trackerFPGA[12];
   EnergyBoard energyBoard[2];
public:
   PCTEvent(): TObject() {}
   void Clear(Option_t*) {
      evt = -1;
      addressError = false;
      tagError = false;
      CRCError = false;
      deltaT = -999;
      for (int iFPGA=0; iFPGA<12; ++iFPGA) trackerFPGA[iFPGA].Clear(0);
      for (int iboard=0; iboard<2; ++iboard) energyBoard[iboard].Clear(0);
   }
   const char* PrintTime() const {
      time_t timeTag_long = timeTag;
      const char* timestr = ctime(&timeTag_long);
      cout<< timestr;
      return timestr;
   }
   ClassDef(PCTEvent, 14);
};

class RunHeader: public TObject {
   // added from version 44
   Int_t run;                          // 24-bit
   Int_t time;                         // 32-bit start time in seconds
   Bool_t timeTag;                     // 8-bit status bit, set if event time tags are written out
   UChar_t version;                    // 8-bit program version number, e.g. 44
public:
   RunHeader(): TObject() {}
   void Fill(char buf_run[3]
             , char buf_time[4]
             , char buf_timeTag
             , char buf_version
   ) {
      union {
         Char_t buf[4];
         Int_t run;
      } run_union;
      run_union.buf[0] = buf_run[2];
      run_union.buf[1] = buf_run[1];
      run_union.buf[2] = buf_run[0];
      run_union.buf[3] = 0;
      run = run_union.run;

      union {
         Char_t buf[4];
         Float_t time;                  // 32 bit run start time, floating point single precision, in seconds
      } time_union;
      time_union.buf[0] = buf_time[3];
      time_union.buf[1] = buf_time[2];
      time_union.buf[2] = buf_time[1];
      time_union.buf[3] = buf_time[0];
      time = (Int_t) time_union.time;   // convert from float to Int_t

      timeTag = buf_timeTag;
      version = buf_version;
   }
   Int_t GetRun() const {return run;}
   Int_t GetTime() const {return time;}
   Bool_t GetTimeTag() const {return timeTag;}
   const char* PrintTime() const {
      time_t time_long = time;
      const char* timestr = ctime(&time_long);
      cout<< timestr;
      return timestr;
   }

   Int_t GetVersion() const {return version;}

   ClassDef(RunHeader, 3);
};

#endif  // DataFormat_h
