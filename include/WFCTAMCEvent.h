#ifndef __WFCTAMCEvent__
#define __WFCTAMCEvent__
#include "TSelector.h"
#include "WFTelescope.h"
#include "WFCamera.h"
#include "common.h"
#include "TH1D.h"
class WFCTAMCEvent
{
   public:
   static bool RecordRayTrace;
   static double fAmpLow;
   static double fAmpHig;
   public:
   Double_t Ngen;
   Int_t Timegen;
   //vector<float> Coogen[3];
   //vector<float> Dirgen[3];
   //vector<float> Wavegen;
   Float_t TubeSignal[NCTMax][NSIPM]; //!
   Float_t eTubeSignal[NCTMax][NSIPM];

   int NArrival[NCTMax]; //!
   double ArrivalTimeMin[NCTMax]; //!
   double ArrivalTimeMax[NCTMax]; //!
   bool OverFlow[NCTMax]; //!
   long int ArrivalTime[NCTMax][MaxTimeBin]; //!
   double ArrivalCount[NCTMax][NSIPM][MaxTimeBin]; //!
   double ArrivalCountE[NCTMax][NSIPM][MaxTimeBin]; //!

   Int_t TubeTrigger[NCTMax][NSIPM]; //!
   Int_t TelTrigger[NCTMax]; //!

   vector<int> RayTrace; //!
   static TH1D* hRayTrace; //! the histogram to store ray tracing
   static TH1D* hWeightRayTrace; //! the histogram to store weighted ray tracing

   public:
   void Init(int size=0);
   WFCTAMCEvent() {Init();}
   ~WFCTAMCEvent() {;}
   void Reset();
   void Copy(WFTelescopeArray* pct=0);
   void GetTubeTrigger();
   void GetTelescopeTrigger(WFTelescopeArray* pct=0);

   //ClassDef(WFCTAMCEvent,1);
};
#endif
