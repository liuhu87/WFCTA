#ifndef __WFCTALaserEvent__
#define __WFCTALaserEvent__
class WFCTALaserEvent{
   public:
   int Time;
   float Frequency;
   float flux;
   public:
   void Init();
   void Reset();
   WFCTALaserEvent() {Init();}
   ~WFCTALaserEvent() {;}
};
#endif
