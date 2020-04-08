#ifndef __CalibWFCTA__
#define __CalibWFCTA__
#define CalibMaxTel 20

#include "TGraphErrors.h"
//#include "camera.h"

class CalibWFCTA {
   private:
   static CalibWFCTA* _Head;
   static char DirName[3][200];

   public:
   static bool ForceCorr;
   static int UseSiPMCalibVer;
   static int jdebug;

   //Version 1
   ///info for the tree
   TFile* fin;
   TTree* tree;
   long int entries;
   long int curentry;

   ///variables in the tree
   ///Tel No.
   Short_t curTel;
   ///rabbit Time
   Long64_t rabbitTime;
   ///BaseH_RMS
   Float_t mBaseH_RMS[1024];
   ///BaseL_RMS
   Float_t mBaseL_RMS[1024];
   ///mBaseH
   Float_t mBaseH[1024];
   ///mBaseL
   Float_t mBaseL[1024];
   ///AdcH
   Float_t mAdcH[1024];
   ///AdcL
   Float_t mAdcL[1024];
   ///AdcL_RMS
   Float_t mAdcL_RMS[1024];
   ///AdcH_RMS
   Float_t mAdcH_RMS[1024];
   ///Low_Flag
   Short_t Low_Flag[1024];
   ///HGain
   Double_t H_Gain_Factor[1024];
   ///LGain
   Double_t L_Gain_Factor[1024];
   ///Gain
   Double_t Gain_Factor[1024];

   //Version 2
   TGraphErrors* sipmcalib_norm[CalibMaxTel];
   TGraphErrors* sipmcalib_slope[CalibMaxTel];
   //Version 3
   double unif_factor[CalibMaxTel][1024];
   double deltag_20[CalibMaxTel][1024];

   void Init();
   void Reset();
   void Release();
   static char* GetDirName(int version=UseSiPMCalibVer) {return DirName[version-1];}
   static void SetDirName(int version,char* dirname);
   int LoadCalibSiPM(int version=-1,char* dirname=0);
   CalibWFCTA() {Init(); LoadCalibSiPM();}
   CalibWFCTA(char* dirname) {Init(); LoadCalibSiPM(UseSiPMCalibVer,dirname);}
   ~CalibWFCTA() {Release();}
   static CalibWFCTA* GetHead();
   void SetBranchAddress();
   long int LoadDay(int Day,int iTel);
   int IsTimeEqual(int time1,int time2);
   int LoadTime(int time);
   int LoadFromrabbitTime(double time,int iTel);
   int LoadEntry(long int entry);
   void Dump(int isipm=16*32);
   double DoCalibSiPM(int iTel,int isipm,double input,double temperature,double Time,int calibtype=0x7,bool IsLED=false);
};

#endif
