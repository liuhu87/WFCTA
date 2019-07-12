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
   static int buflength;
   static char FILENAME[300];
   static FILE* fp;
   static int readtime;
   static vector<int> failtime;

   static WFCTADecode *wfctaDecode;
   public:
   static int timemargin;

  //int fpgaVersion[10];
  //long clb_initial_Time;
  //double clb_initial_time;
  //int fired_tube;
  //long status_readback_Time;
  //double status_readback_time;
  ////int sipm[1024];  for(int i=0;i<1024;i++) {sipm[i]=i;}
  //short single_thresh[1024];
  //short record_thresh[1024];
  //long single_count[1024];
  //float DbTemp[1024];
  //long single_time[1024];
  //float HV[1024];
  //float PreTemp[1024];
  //float BigResistence[1024];
  //float SmallResistence[1024];
  //long ClbTime[1024];
  //float ClbTemp[1024];


   static void Init();
   static void Reset();
   static void Release();
   StatusDB() {Init();}
   ~StatusDB(){Release();}
   static void SetDirectory(const char* dirname);
   static StatusDB* GetHead();
   static FILE* LocateFile(int Time,double time=0);
   static bool LocateBlk(int Time,double time=0);
   //static bool Fill(int Time,double time=0);
};

#endif
