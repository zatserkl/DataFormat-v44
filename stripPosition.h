// Robert Johnson, March 12, 2014

inline float stripPositionV(int fpga, int chip, int strip) {

   // Returns the V coordinate of the center of a strip in a tracker V board.
   // The origin is taken to be roughly the center of the SSD, which should also line
   // up with the center of the energy detector.
   // (In principle we could also input information about the T coordinate in order
   // to know which of the two sensors on a given side was hit, in case that the
   // alignment of the sensors relative to each other is significantly non-zero.)

   int VBoard[12] = {6,4,2,3};  	//fpga to V-board translation; this will change if spares are swapped
   float Vpin[4] = {0.055,0.055,0.055,0.055};  //V coordinate of alignment pin per fpga
   float firstStrip[7][2] = {		//Distance from the alignment pin to the first strip
      {-43.7193, -43.716},    //Board V0 doesn't exist
      {-43.7193, -43.716},
      {-43.7193, -43.716},
      {-43.7193, -43.716},
      {-43.7193, -43.716},	  //These numbers are actually from V4
      {-43.7193, -43.716},
      {-43.5855, -43.5865}   //these are actually from T6
   };

   int gStrip, side;
   if (chip < 6) {
      gStrip = (5 - chip)*64 + strip;
      side = 0;
   } else {
      side = 1;
      gStrip = (chip - 6)*64 + (63 - strip);
   }
   int brd = VBoard[fpga];

   return firstStrip[brd][side] + Vpin[fpga] + gStrip*(0.228);
}

inline float stripPositionT(int fpga, int chip, int strip) {
   // fpga is the readout address, which for T boards runs from 4 to 11.
   // chip is the chip number from 0 to 11 that is read by the given fpga.
   // strip is the internal strip number in the chip, from 0 to 63.

   // The layer ordering is, going along the beam direction:
   // 0= V   0
   // 1= T   4 5       T-layer 0
   // 3= V   1
   // 4= T   6 7       T-layer 1
   // 5= T   8 9       T-layer 2
   // 6= V   2
   // 7= T   10 11     T-layer 3
   // 8= V   3
   // Here I defined an ad-hoc numbering for the T layers themselves

   int Tlayer[12] = {-1, -1, -1, -1, 0, 0, 1, 1, 2, 2, 3, 3};  //fpga to T-layer translation
   float Tpin[4] = {215.24, 211.25, -211.25, -215.24};  //T coordinate of alignment pin per layer
   float Tdir[4] = {-1.0, -1.0, 1.0, 1.0};   //T direction of increasing strip number per layer

   int Tboard[4] = {5, 4, 1, 3};  //T-layer to physical T board translation. This will change if spares are swapped in!
   float firstStrip[7][4] = {    //First strip location rel to pin for each sensor on each physical board
      {-999., -999., -999., -999.},	   //Board 0 doesn't exist.  We manufactured 6 boards.
      {38.58, 126.85, 215.11, 303.37}, // #1
      {38.58, 126.85, 215.11, 303.37}, // #2
      {38.58, 126.85, 215.11, 303.37}, // #3
      {38.58, 126.85, 215.11, 303.37}, // #4   //These numbers actually come from board T4
      {38.62, 126.90, 215.16, 303.41}, // #5   //these numbers actually come from board T5
      {38.58, 126.85, 215.11, 303.37}  // #6 
   };

   int sensor = 2*(fpga % 2) + chip/6;
   int chipInSensor = chip % 6;

   int gStrip = 64*chipInSensor + (63 - strip);
   int lyr= Tlayer[fpga];
   int brd = Tboard[lyr];
   return Tpin[lyr] + Tdir[lyr]*(firstStrip[brd][sensor]) + gStrip*0.228;
}
