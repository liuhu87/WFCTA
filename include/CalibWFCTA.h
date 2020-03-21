#ifndef __CalibWFCTA__
#define __CalibWFCTA__
#define CalibMaxTel 20

#include "TGraphErrors.h"
//#include "camera.h"

class CalibWFCTA {
   private:
   static CalibWFCTA* _Head;

   public:
   static int UseSiPMCalibVer;
   TGraphErrors* sipmcalib_norm[CalibMaxTel];
   TGraphErrors* sipmcalib_slope[CalibMaxTel];
   double unif_factor[CalibMaxTel][1024];
   double deltag_20[CalibMaxTel][1024];

   void Init();
   void Release();
   int LoadCalibSiPM(int version=UseSiPMCalibVer,char* dirname=0);
   CalibWFCTA() {Init(); LoadCalibSiPM();}
   CalibWFCTA(char* dirname) {Init(); LoadCalibSiPM(-1,dirname);}
   ~CalibWFCTA() {Release();}
   static CalibWFCTA* GetHead(char* dirname=0);
   double DoCalibSiPM(int iTel,int isipm,double input,double temperature,int Time,int calibtype=0x7);
};

#endif
