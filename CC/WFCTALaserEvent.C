#include "WFCTALaserEvent.h"
void WFCTALaserEvent::Init(){
   Reset();
}
void WFCTALaserEvent::Reset(){
   Time=0;
   for(int ii=0;ii<3;ii++) LaserCoo[ii]=0;
   for(int ii=0;ii<2;ii++) LaserDir[ii]=0;
   wavelength=0;
   Intensity=0;
   Frequency=0;
   pulsetime=0;
}

