#ifndef __Cloud__
#define __Cloud__
#include "TH2Poly.h"
#include "TGraph.h"
#include <vector>
using std::vector;

#define Cntheta 40
#define Cnphi 158

class WFTelescopeArray;
class Cloud {
   public:
   static int drawcircle;
   static double Cbintheta[Cntheta+1];
   static int Cnbinphi[Cntheta];
   static double Cbinphi[Cntheta][Cnphi+1];
   TH2Poly* cloudmap;
   int time;
   double temp;
   vector<TGraph*> graphlist;

   public:
   static void SetBins();
   static void Convert(int index,double &xx,double &yy);
   static int FindBinIndex(double xx,double yy);
   static TGraph* TelView(WFTelescopeArray* pct,int iTel);

   void Init();
   void Reset();
   void Clear();
   Cloud() { Init(); }
   ~Cloud() { Clear(); }
   void ReadCloudMap(char* filename);
   bool ReadTemp(char* filename=0);
   void AveTemp(double &avetemp,double &mintemp,TGraph* gr);
   void Draw(WFTelescopeArray* pct,char* opt=(char*)"colz");
};
#endif
