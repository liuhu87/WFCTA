#ifndef __WFCTALedEvent__
#define __WFCTALedEvent__
class WFCTALedEvent{
   public:
   int Time;
   int Frequency;
   bool DoorOpen;

   public:
   void Init();
   void Reset();
   WFCTALedEvent() {Init();}
   ~WFCTALedEvent() {;}
};
#endif
