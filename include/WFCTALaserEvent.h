#ifndef __WFCTALaserEvent__
#define __WFCTALaserEvent__
class WFCTALaserEvent{
   public:
   int Time;
   float LaserCoo[3];
   float LaserDir[2];
   float wavelength;
   float Intensity;
   float Frequency;
   float pulsetime;
   public:
   void Init();
   void Reset();
   WFCTALaserEvent() {Init();}
   ~WFCTALaserEvent() {;}

   //ClassDef(WFCTALaserEvent,1);
};
#endif
