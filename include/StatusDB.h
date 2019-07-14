#ifndef __Status_DB__
#define __Status_DB__
#include <stdint.h>
#include <stdio.h>
#include "WFCTADecode.h"
using namespace std;

#define BUF_LEN 41552
#define STATUS_BUF_LEN 41552

class StatusDB {
   private:
   static StatusDB* _Head;
   static char DBPath[200];

   static uint8_t *buf;
   static unsigned long buflength;

   static char FILENAME[300];
   static FILE* fp;
   static vector<int> containtime;
   static int readtime;

   static vector<int> failtime;

   static WFCTADecode *wfctaDecode;
   public:
   static int jdebug;
   static int timemargin;

   static int fpgaVersion[10];
   static long clb_initial_Time;
   static double clb_initial_time;
   static int fired_tube;
   static long status_readback_Time;
   static double status_readback_time;
   static short single_thresh[1024];
   static short record_thresh[1024];
   static long single_count[1024];
   static float DbTemp[1024];
   static long single_time[1024];
   static float HV[1024];
   static float PreTemp[1024];
   static float BigResistence[1024];
   static float SmallResistence[1024];
   static long ClbTime[1024];
   static float ClbTemp[1024];


   static void Init();
   static void Reset();
   static void Release();
   StatusDB() {Init();}
   ~StatusDB(){Release();}
   static void SetDirectory(const char* dirname);
   static StatusDB* GetHead();
   static FILE* LocateFile(int Time,double time=0);
   static bool LocateBlk(int Time,double time=0);
   static void ResetBuffer();
   static void Fill();
   static bool Locate(int Time,double time=0);

   static int   GetVersion(int Time,int i);
   static long  GetClbInitTime(int Time);
   static double GetClbInittime(int Time);
   static int   GetFiredTube(int Time);
   static long  GetReadbackTime(int Time);
   static double GetReadbacktime(int Time);
   static short GetSingleThrd(int Time,int i);
   static short GetRecordThrd(int Time,int i);
   static long  GetSingleCount(int Time,int i);
   static float GetDbTemp(int Time,int i);
   static long  GetSingleTime(int Time,int i);
   static float GetHV(int Time,int i);
   static float GetPreTemp(int Time,int i);
   static float GetBigR(int Time,int i);
   static float GetSmallR(int Time,int i);
   static long  GetClbTime(int Time,int i);
   static float GetClbTemp(int Time,int i);
};

#endif
