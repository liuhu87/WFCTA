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
   Float_t TubeSignal[NCTMax][NSIPM];
   Float_t eTubeSignal[NCTMax][NSIPM];
   Float_t ArrivalTimeMin[NCTMax][NSIPM];
   Float_t ArrivalTimeMax[NCTMax][NSIPM];
   Float_t ArrivalAccTime[NCTMax][NSIPM];
   Int_t NArrival[NCTMax][NSIPM];
   Int_t TubeTrigger[NCTMax][NSIPM];
   Int_t TelTrigger[NCTMax];

   vector<int> RayTrace; //!
   static TH1D* hRayTrace; //! the histogram to store ray tracing

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
