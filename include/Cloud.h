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
   static bool DoTempCorr;
   static bool drawmoon;
   static int drawcircle;
   static double Cbintheta[Cntheta+1];
   static int Cnbinphi[Cntheta];
   static double Cbinphi[Cntheta][Cnphi+1];
   TH2Poly* cloudmap;
   int mapnbins;
   int time;
   double temp0;
   double temp;
   double humi;
   vector<TGraph*> graphlist;

   public:
   static void SetBins();
   static bool Convert(int index,double &xx,double &yy,double xyboun[][2]);
   static int FindBinIndex(double xx,double yy);
   static TGraph* TelView(WFTelescopeArray* pct,int iTel);
   static void LoadTelSetting(char* filename);

   void Init();
   void Reset();
   void Clear();
   Cloud() { Init(); }
   ~Cloud() { Clear(); }
   int GetNbins() {return mapnbins;}
   void ReadCloudMap(char* filename);
   bool ReadTemp(char* filename=0);
   void AveTemp(double &avetemp,double &mintemp,double &rmstemp,TGraph* gr);
   void Draw(WFTelescopeArray* pct,char* opt=(char*)"colz");
   double GetCorrected(double input);
   double GetTemperature(int itemp);
   double GetHumidity();
   ///get the IB temperature at coor. (xx,yy)
   double GetIBTemp(int ibin);
   double GetIBTemp(double xx,double yy);
   double GetTelAveIBTemp(int iTel);
   double GetTelMinIBTemp(int iTel);
   double GetTelRmsIBTemp(int iTel);
   double GetAveIBTemp(double theta);
   double GetMinIBTemp(double theta);
   double GetRmsIBTemp(double theta);
   double GetAveIBTemp(double theta1,double theta2);
   double GetMinIBTemp(double theta1,double theta2);
   double GetRmsIBTemp(double theta1,double theta2);
   double GetAveIBTemp(TGraph* gr);
   double GetMinIBTemp(TGraph* gr);
   double GetRmsIBTemp(TGraph* gr);
};
#endif
