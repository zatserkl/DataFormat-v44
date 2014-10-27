#include <Rtypes.h>

struct EventOutput {
   Bool_t good;
   Float_t t[4];
   Float_t v[4];
   Float_t u[4];
   Float_t wepl;
   EventOutput() {clear();}
   void clear() {
      good = kFALSE;
      for (int ilayer=0; ilayer<4; ++ilayer) {
         t[ilayer] = 0;
         v[ilayer] = 0;
         u[ilayer] = 0;
      }
      wepl = 0;
   }
};
