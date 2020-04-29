#ifndef __CalibWFCTA__
#define __CalibWFCTA__

#include "TGraphErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "common.h"
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
   TTree* tree2;
   long int entries;
   long int curentry;

   ///variables in the tree
   ///HL Gain Ratio
   Double_t Ratio_HL[MAXPMT];
   ///Dead Channel
   Short_t Dead_Channel_Flag[MAXPMT];
   ///Tel No.
   Short_t curTel;
   ///rabbit Time
   Long64_t rabbitTime;
   ///BaseH_RMS
   Float_t mBaseH_RMS[MAXPMT];
   ///BaseL_RMS
   Float_t mBaseL_RMS[MAXPMT];
   ///mBaseH
   Float_t mBaseH[MAXPMT];
   ///mBaseL
   Float_t mBaseL[MAXPMT];
   ///AdcH
   Float_t mAdcH[MAXPMT];
   ///AdcL
   Float_t mAdcL[MAXPMT];
   ///AdcL_RMS
   Float_t mAdcL_RMS[MAXPMT];
   ///AdcH_RMS
   Float_t mAdcH_RMS[MAXPMT];
   ///Low_Flag
   Short_t Low_Flag[MAXPMT];
   ///HGain
   Double_t H_Gain_Factor[MAXPMT];
   ///LGain
   Double_t L_Gain_Factor[MAXPMT];
   ///Gain
   Double_t Gain_Factor[MAXPMT];

   //Version 2
   TGraphErrors* sipmcalib_norm[NCTMax];
   TGraphErrors* sipmcalib_slope[NCTMax];
   //Version 3
   double unif_factor[NCTMax][MAXPMT];
   double deltag_20[NCTMax][MAXPMT];

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
   static int IsTimeEqual(int time1,int time2);
   int LoadTime(int time);
   int LoadFromrabbitTime(double time,int iTel);
   int LoadEntry(long int entry);
   void Dump(int isipm=16*32);
   double DoCalibSiPM(int iTel,int isipm,double input,double temperature,double Time,int calibtype=0x7,int type=0);
};

#endif
