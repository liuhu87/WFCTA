#ifndef __WFCTAMCEvent__
#define __WFCTAMCEvent__
#include "TSelector.h"
#include "WFTelescope.h"
#include "WFCamera.h"
#include "common.h"
class WFCTAMCEvent
{
   public:
   Int_t iuse;
   long int Ngen;
   Int_t Timegen;
   //vector<float> Coogen[3];
   //vector<float> Dirgen[3];
   //vector<float> Wavegen;
   vector<Int_t> RayTrace;
   Float_t TubeSignal[NCTMax][NSIPM];
   Int_t TubeTrigger[NCTMax][NSIPM];
   Int_t TelTrigger[NCTMax];

   public:
   void Init(int size=0);
   WFCTAMCEvent() {Init();}
   ~WFCTAMCEvent() {;}
   void Reset();
   void Copy(WFTelescopeArray* pct=0);
   void GetTubeTrigger();
   void GetTelescopeTrigger(WFTelescopeArray* pct=0);
};
#endif
