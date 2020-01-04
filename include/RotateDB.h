#ifndef __RotateDB__
#define __RotateDB__
#define NSWITH 19
#define NVALUE 12

#include "common.h"

class RotateDB {
   private:
   static RotateDB* _Head;
   public:
   static int jdebug;
   static int timedelay;
   static int ntotmin;
   static int nsidemin;
   char buff[500];
   //Li
   int Li;
   //time information
   int time;
   //position in the log file
   long int currpos;
   //gps information
   double latitude;
   double longitude;
   double altitude;
   //swith information
   bool allswith[NSWITH];
   //other variables
   double varinfo[NVALUE];

   public:
   void Init();
   void Reset();
   void Release();
   RotateDB() {Init();}
   RotateDB(RotateDB* pr_in);
   ~RotateDB() {Release();}
   void Copy(RotateDB* pr_in);
   static RotateDB* GetHead();
   static bool LocateFirst(ifstream* fin);
   long int LoadData(int time_in,int Li_in,int pLi=0,int ptime=-1,long int cpos=-1);
   int ProcessTime();
   void ProcessAll();
   void DumpInfo();
   bool GetLaserSwith();
   bool GetDoorSwith();
   double GetTemperature(int itemp);
   double GetHumidity();
   double GetInclination(bool IsX);
   double GetDoorAngle();
   double GetHeight();
   void GetAngles(double angle[3]);
   double GetElevation();
   double GetAzimuth();
   double GetAng();
   int IsFineAngle(double ele_in,double azi_in,int iTel,int &index);
   int GetEleAzi(int time_in,int Li_in,int iTel=-1);
};

#endif
