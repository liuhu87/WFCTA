#include "WFCTALaserEvent.h"
bool WFCTALaserEvent::Recordweight=false;
TH1D* WFCTALaserEvent::hweight=new TH1D("Laser_weight",";log10(weight);Entries",100,-20,20);
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
   weight.clear();
   hweight->Reset();
}

