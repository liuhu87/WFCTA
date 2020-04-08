#ifndef __CalibLED__
#define __CalibLED__

#include "common.h"
#include "TFile.h"
#include "TTree.h"

class WFCTAEvent;
class CalibLED{
   private:
   static CalibLED* _Head;
   static char DirName[200];
   public:
   static bool ForceCorr;
   static int TimeDelay;
   static int jdebug;
   ///Tel No.
   int curTel;
   ///info for the tree
   TFile* fin;
   TTree* tree;
   long int entries;
   long int curentry;
   ///variables in the tree
   ///MJD Time
   double T1MJD;
   ///LED Drive Temp Correction Factor
   double LED_DRF;
   ///Current of the high voltage
   double T1HVCurrent;
   ///High Voltage
   double T1HV;
   ///temperature of the drive board of led which point to sipm
   double T1DledDrDStem;

   public:
   void Init();
   void Reset();
   void Clear();
   CalibLED() {Init();}
   ~CalibLED() {Clear();}
   static char* GetDirName() {return DirName;}
   static void SetDirName(char* dirname);
   static CalibLED* GetHead(char* dirname=0);
   void SetBranchAddress();
   long int LoadDay(int Day,int iTel);
   int LoadTime(int time);
   int LoadEntry(long int entry);
   int LoadFromrabbitTime(double time,int iTel);
   void Dump();
   double DoLedDriveTempCorr(double input,int isipm,double time=0,int iTel=0);
   double DoLedDriveTempCorr(double input,int isipm,WFCTAEvent* pev=0);
};
#endif
