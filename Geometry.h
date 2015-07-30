// Andriy Zatserklyaniy, May 5, 2014

#ifndef Geometry_h
#define Geometry_h

#include <TROOT.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using std::cout;     using std::endl;

class Geometry {
public:
   static Geometry* geometry_;
public:
   const Double_t pitch_;                   // = 0.228 mm
   const Double_t offset_rear_telescope_;   // 8 mm in positive t-direction

   // V-board constants
   Double_t vPin_[4];               // coordinate of the pin for the board in each layer
   Int_t vBoardLayer_[4];           // Id of the board for each layer
   Double_t firstStripV_[10][2];    // distance from alignment pin to the first strip for each sensor for all board Ids
   Double_t vdir_[8];               // direction (+1 or -1) of the sensor's abs strip number (same as the sensor chip address) wrt v-axis
   Double_t uv_[4];                 // u-coordinate of each v-layer

   // T-board constants
   Double_t tPin_[4];               // coordinate of the pin for the board in each layer
   Int_t tBoardLayer_[4];           // Id of the board for each layer
   Double_t firstStripT_[10][4];    // distance from alignment pin to the first strip for each sensor for all board Ids
   Double_t tdir_[16];              // direction (+1 or -1) of the sensor's abs strip number (same as the sensor chip address) wrt t-axis
   Double_t ut_[4];                 // u-coordinate of each t-layer
   Double_t angle_;
public:
   Geometry(Int_t the_run, const char* dbname, Double_t offset_rear_telescope=8.): pitch_(0.228)
                                                                                   , offset_rear_telescope_(offset_rear_telescope)
                                                                                   , angle_(-1)
   {
      geometry_ = this;

      // read the database
      //-- assert(readdb(the_run, dbname));
      //--beamtest-- if (!readdb(the_run, dbname)) readdb(0,dbname);     // if run was not found in the database use settings for the run 0
      if (!readdb(the_run, dbname)) {
        cout<< "Run " << the_run << " was not found in the database" <<endl;
        exit(0);
      }

      // set parameters which are constants

      //
      // V-boards
      //

      Double_t vPin[4] = {0.055,0.055,0.055,0.055};
      //int vBoardLayer_[4] = {6,4,2,3};  	      // fpga to V-board translation; this will change if spares are swapped
      //--March-- Double_t firstStripV[7][2] = {		// Distance from the alignment pin to the first strip
      //--March--    {-43.7193, -43.716},             // Board V0 doesn't exist
      //--March--    {-43.7193, -43.716},
      //--March--    {-43.7193, -43.716},
      //--March--    {-43.7193, -43.716},
      //--March--    {-43.7193, -43.716},	            // These numbers are actually from V4
      //--March--    {-43.7193, -43.716},
      //--March--    {-43.5855, -43.5865}             // these are actually from T6
      //--March-- };
      // Double_t firstStripV[7][2] = {		// Distance from the alignment pin to the first strip
      //    {-999., -999.},                  // Board V0 doesn't exist, we manufactured 6 boards
      //    {-43.727, -43.687},              // Board V1
      //    {-43.682, -43.681},              // Board V2
      //    {-43.702, -43.713},              // Board V3
      //    {-43.7193, -43.716},	          // Board V4
      //    {-43.7193, -43.716},             // No survey yet for Board V5
      //    {-43.686, -43.687}               // Board V6, including a typo correction June 9 that changed the numbers by 0.1 mm
      // };

      // Robert Johnson, 10/2/2014
      Double_t firstStripV[7][2] = {		// Distance from the alignment pin to the first strip
         {-999., -999.},                  // Board V0 doesn't exist, we manufactured 6 boards
         {-43.727, -43.687},              // Board V1
         {-43.7033, -43.7185},              // Board V2
         {-43.8453, -43.7855},              // Board V3
         {-43.5806, -43.6415},	          // Board V4
         {-43.7193, -43.716},             // No survey yet for Board V5
         {-43.6603, -43.6515}               // Board V6, including a typo correction June 9 that changed the numbers by 0.1 mm
      };

      for (int ilayer=0; ilayer<4; ++ilayer) vPin_[ilayer] = vPin[ilayer];
      // clear the 10x2 array
      for (int i=0; i<10; ++ i) {
         for (int isensor=0; isensor<2; ++isensor) {
            firstStripV_[i][isensor] = 0;
         }
      }
      // assign the first 7 measured values
      for (int i=0; i<7; ++i) {
         for (int isensor=0; isensor<2; ++isensor) {
            firstStripV_[i][isensor] = firstStripV[i][isensor];
         }
      }
      // directions for the sensors
      for (int ilayer=0; ilayer<4; ++ilayer) {
         // vdir_[0] = -1;
         // vdir_[1] = +1;
         vdir_[2*ilayer]   = -1;
         vdir_[2*ilayer+1] = +1;
      }

      //
      // T-boards
      //
      //-- Double_t tPin[4] = {215.24, 211.25, -211.25, -215.24};   // T coordinate of alignment pin per layer
      Double_t tPin[4] = {215.3125, 211.3595, -211.1630, -215.1940};   // T coordinate of alignment pin per layer //NB: Robert Johnson 10/2/2014
      //Double_t tDir[4] = {-1.0, -1.0, 1.0, 1.0};               // T direction of increasing strip number per layer

      //Int_t tBoard[4] = {5, 4, 1, 3};     // T-layer to physical T board translation. This will change if spares are swapped in!
      //--March-- Double_t firstStripT[7][4] = {      // First strip location rel to pin for each sensor on each physical board
      //--March--    {-999., -999., -999., -999.},		// Board 0 doesn't exist.  We manufactured 6 boards.
      //--March--    {38.58, 126.85, 215.11, 303.37},
      //--March--    {38.58, 126.85, 215.11, 303.37},
      //--March--    {38.58, 126.85, 215.11, 303.37},
      //--March--    {38.58, 126.85, 215.11, 303.37}, // These numbers actually come from board T4
      //--March--    {38.62, 126.90, 215.16, 303.41}, // these numbers actually come from board T5
      //--March--    {38.58, 126.85, 215.11, 303.37} 
      //--March-- };
      Double_t firstStripT[7][4] = {      // First strip location rel to pin for each sensor on each physical board
         {-999., -999., -999., -999.},		// Board 0 doesn't exist.  We manufactured 6 boards.
         {38.60, 126.87, 215.15, 303.42},   // Board T1
         {38.48, 126.76, 215.04, 303.32},   // Board T2
         {38.69, 126.95, 215.23, 303.57},   // Board T3
         {38.58, 126.85, 215.11, 303.37},   // Board T4
         {38.62, 126.90, 215.16, 303.41},   // Board T5
         {38.58, 126.85, 215.11, 303.37}    // Board T6, no survey has been done
      };

      for (int ilayer=0; ilayer<4; ++ilayer) tPin_[ilayer] = tPin[ilayer];
      // clear the 10x4 array
      for (int i=0; i<10; ++ i) {
         for (int isensor=0; isensor<4; ++isensor) {
            firstStripT_[i][isensor] = 0;
         }
      }
      // assign the first 7 measured values
      for (int i=0; i<7; ++i) {
         for (int isensor=0; isensor<4; ++isensor) {
            firstStripT_[i][isensor] = firstStripT[i][isensor];
         }
      }
      // directions for the sensors
      for (int ilayer=0; ilayer<2; ++ilayer) {
         for (int isensor=0; isensor<4; ++isensor) {
            tdir_[isensor] = -1;                         // sensors in the first two layers are in the opposite direction
         }
      }
      for (int ilayer=2; ilayer<4; ++ilayer) {
         for (int isensor=0; isensor<4; ++isensor) {
            tdir_[isensor] = +1;                         // sensors in the last two layers are in the same direction
         }
      }
   }
   bool readdb(Int_t the_run, const char* dbname) {
      //cout<< "readdb: the_run = " << the_run << " dbname = " << dbname <<endl;

      Int_t run = 0;
      for (int ilayer=0; ilayer<4; ++ilayer) {
         uv_[ilayer] = 0;
         ut_[ilayer] = 0;
         vBoardLayer_[ilayer] = -1;
         tBoardLayer_[ilayer] = -1;
      }
      std::string phantom;
      Int_t nbricks = -1;

      const std::string run_str = "run";
      const std::string uvt_str = "uvt";
      const std::string v_board_str = "v-board";
      const std::string t_board_str = "t-board";
      const std::string phantom_str = "phantom";
      const std::string none_str = "none";
      const std::string phantom_pyramid_str = "pyramid";
      const std::string nbricks_str = "nbricks";
      const std::string phantom_water_str = "water";
      const std::string phantom_wire_str = "wire";
      const std::string phantom_catphan_str = "catphan";
      const std::string phantom_sensitom_str = "sensitom";
      const std::string phantom_CIRSFordInf_str = "CIRSFordInf";
      const std::string phantom_CIRSFordSup_str = "CIRSFordSup";

      std::ifstream file(dbname);
      if (!file) {
         cout<< "Could not open file " << dbname <<endl;
         return false;
      }

      cout<< "processing file " << dbname <<endl;

      bool found = false;
      std::string line;
      std::string word;
      while (std::getline(file, line))
      {
         //cout<< "line = " << line <<endl;
         if (line.size() == 0) continue;

         std::istringstream ss(line);          // create a new instance of the stringstream for each iteration to avoid problem with eof

         // read the first word
         ss >> word;

         if (word == run_str) {     // data line should start from the word "run", other lines will be ignored
            ss >> run;
            if (run == the_run) {
               found = true;
               break;
            }
            else continue;
         }
         else continue;
      }

      if (!found) {
         cout<< "Run " << the_run << " was not found in the database" <<endl;
         return false;
      }

      // found line for the_run
      cout<< line <<endl;

      std::istringstream ss(line);
      ss >> word >> run;              // remove from the stream already read items

      // u coordinates
      ss >> word;
      if (word != uvt_str) {
         cout<< "read error while expected " << uvt_str <<endl;
         return false;
      }
      for (int ilayer=0; ilayer<4 && ss >> uv_[ilayer] >> ut_[ilayer]; ++ilayer) {
         cout<< ilayer << "\tuv_[" << ilayer << "] = " << uv_[ilayer] << "\t ut_[" << ilayer << "] = " << ut_[ilayer] <<endl;
      }
      if (ut_[3] == 0) {
         cout<< "error in reading uvt values" <<endl;
         return false;
      }

      // v-boards
      ss >> word;
      if (word != v_board_str) {
         cout<< "read error while expected " << v_board_str <<endl;
         return false;
      }
      cout<< "tBoard: ";
      for (int ilayer=0; ilayer<4 && ss >> vBoardLayer_[ilayer]; ++ilayer) {
         cout<< vBoardLayer_[ilayer] << " ";
      }
      cout<<endl;
      if (vBoardLayer_[3] == -1) {
         cout<< "error in reading " << v_board_str << " values" <<endl;
         return false;
      }

      // t-boards
      ss >> word;
      if (word != t_board_str) {
         cout<< "read error while expected " << t_board_str <<endl;
         return false;
      }
      cout<< "tBoard: ";
      for (int ilayer=0; ilayer<4 && ss >> tBoardLayer_[ilayer]; ++ilayer) {
         cout<< tBoardLayer_[ilayer] << " ";
      }
      cout<<endl;
      if (tBoardLayer_[3] == -1) {
         cout<< "error in reading " << t_board_str << " values" <<endl;
         return false;
      }

      // phantom
      ss >> word;
      if (word != phantom_str) {
         cout<< "read error while expected \"phantom\"" <<endl;
         return false;
      }
      ss >> phantom;
      cout<< "--> phantom = " << phantom <<endl;
      if (phantom == phantom_pyramid_str) {
          // read the word "nbricks"
          ss >> word;
          if (word == nbricks_str) {
              nbricks = -1;
              ss >> nbricks;
              if (nbricks == -1) {
                  cout<< "read error while expected integer number of bricks" <<endl;
                  return false;
              }
              cout<< "phantom: " << phantom << " nbricks = " << nbricks <<endl;
          }
          else {
              cout<< "expected word \"nbricks\"" <<endl;
              return false;
          }
      }
      else if (phantom == phantom_water_str) {
         ss >> word >> angle_;
         cout<< "phantom == " << phantom_water_str << " angle_ = " << angle_ <<endl;
      }
      else if (phantom == phantom_wire_str) {
         cout<< "phantom == " << phantom_wire_str <<endl;
      }
      else if (phantom == phantom_catphan_str) {
         ss >> word >> angle_;
         cout<< "phantom == " << phantom_catphan_str << " angle_ = " << angle_ <<endl;
      }
      else if (phantom == phantom_sensitom_str) {
         ss >> word >> angle_;
         cout<< "phantom == " << phantom_sensitom_str << " angle_ = " << angle_ <<endl;
      }
      else if (phantom == phantom_CIRSFordInf_str) {
         ss >> word >> angle_;
         cout<< "phantom == " << phantom_CIRSFordInf_str << " angle_ = " << angle_ <<endl;
      }
      else if (phantom == phantom_CIRSFordSup_str) {
         ss >> word >> angle_;
         cout<< "phantom == " << phantom_CIRSFordSup_str << " angle_ = " << angle_ <<endl;
      }
      else if (phantom == none_str) {
         cout<< "phantom == " << none_str <<endl;
      }
      else {
         cout<< "read error of " << phantom << " while expected a name of the phantom " <<endl;
         return false;
      }

      return true;
   }
};

#endif  // Geometry_h
