#ifndef __ShowerPlot__
#define __ShowerPlot__
#define NShowerTrack 4
#include "ReadTrack.h"
#include "CorsikaEvent.h"
#include "TCanvas.h"
#include "TLegend.h"
class ShowerPlot{
   public:
   static int jdebug;
   static char tpname[NShowerTrack][20];
   char primary[20];
   char filename[NShowerTrack][100];
   ReadTrack* ptrk[NShowerTrack];
   CorsikaIO* corsika;
   CorsikaEvent* pevt;
   TObjArray* plot;

   public:
   void Init(const char* priname);
   void Release();
   ShowerPlot(const char* priname) {Init(priname);}
   ~ShowerPlot() {Release();}
   bool Add(const char* filename,int type);
   bool Read();
   TCanvas* Draw(int ViewOpt=0);
};

#endif
