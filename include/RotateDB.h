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
   ~RotateDB() {Release();}
   static RotateDB* GetHead();
   static bool LocateFirst(ifstream* fin);
   bool LoadData(int time_in,int Li_in);
   int ProcessTime();
   void ProcessAll();
   void DumpInfo();
   bool GetEleAzi(int time_in,int Li_in);
};

#endif
