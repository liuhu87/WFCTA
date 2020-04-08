#ifndef __WFCTASCDB__
#define __WFCTASCDB__

#include "common.h"
#include "TFile.h"
#include "TTree.h"

class WFCTAEvent;
class WFCTASCDB{
   private:
   static WFCTASCDB* _Head;
   static char DirName[200];
   public:
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
   ///High Voltage
   double T1HV;
   ///Current of the high voltage
   double T1HVCurrent;
   ///Nagative Voltage
   double T1NV;
   ///Current for Nagative Voltage
   double T1NVCurrent;
   ///Positive Voltage 1
   double T1P1V;
   ///Current for Positive Voltage 1
   double T1P1VCurrent;
   ///Positive Voltage 2
   double T1P2V;
   ///Current for Positive Voltage 2
   double T1P2VCurrent;
   ///Positive Voltage 3
   double T1P3V;
   ///Current for Positive Voltage 3
   double T1P3VCurrent;
   ///control temperature of the led which point to sipm
   double T1Dledtem;
   ///temperature of the led which point to sipm
   double T1DledDStem;
   ///control temperature of the drive board of led which point to sipm
   double T1DledDrtem;
   ///temperature of the drive board of led which point to sipm
   double T1DledDrDStem;
   ///control temperature of the led which point to mirror
   double T1Rledtem;
   ///temperature of the led which point to mirror
   double T1RledDStem;
   ///control temperature of the drive board of led which point to mirrot
   double T1RledDrtem;
   ///temperature of the drive board of led which point to mirror
   double T1RledDrDStem;
   ///temperature of wind input side of the sipm
   double T1CaImtem;
   ///temperature of wind output side of the sipm
   double T1CaExtem;
   ///Angle of Telescope Door 1 in degree
   double T1Door1Deg;
   ///Angle of Telescope Door 2 in degree
   double T1Door2Deg;
   ///Humidity
   double T1HM;
   ///temperature of front door
   double T1TemF;
   ///temperature of middle door
   double T1TemM;
   ///temperature of back door
   double T1TemB;
   ///X inclination angle 
   double T1IncXDeg;
   ///Y inclination angle 
   double T1IncYDeg;
   ///Elevation Angle
   double T1EleDeg;
   ///pulse width of led which point to sipm
   int T1DledPW;
   ///pulse frequency of led which point to sipm
   int T1DledPF;
   ///pulse width of led which point to mirror
   int T1RledPW;
   ///pulse frequency of led which point to mirror
   int T1RledPF;

   public:
   void Init();
   void Reset();
   void Clear();
   WFCTASCDB() {Init();}
   ~WFCTASCDB() {Clear();}
   static char* GetDirName() {return DirName;}
   static void SetDirName(char* dirname);
   static WFCTASCDB* GetHead(char* dirname=0);
   void SetBranchAddress();
   long int LoadDay(int Day,int iTel);
   int LoadTime(int time);
   int LoadEntry(long int entry);
   int LoadFromrabbitTime(double time,int iTel);
   void Dump();
   bool IsDLed(WFCTAEvent* pev=0);
   bool IsRLed(WFCTAEvent* pev=0);
   bool DoorOpened(WFCTAEvent* pev=0);
};
#endif
