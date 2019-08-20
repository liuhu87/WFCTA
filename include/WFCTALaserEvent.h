#ifndef __WFCTALaserEvent__
#define __WFCTALaserEvent__
#include "common.h"
#include "TH1D.h"
class WFCTALaserEvent{
   public:
   static bool Recordweight;
   int Time;
   float LaserCoo[3];
   float LaserDir[2];
   float wavelength;
   float Intensity;
   float Frequency;
   float pulsetime;

   vector<double> weight; //!
   static TH1D* hweight; //! the histogram to store the weight used

   public:
   void Init();
   void Reset();
   WFCTALaserEvent() {Init();}
   ~WFCTALaserEvent() {;}

   //ClassDef(WFCTALaserEvent,1);
};
#endif
