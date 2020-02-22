#ifndef __Status_DB__
#define __Status_DB__
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include "WFCTADecode.h"
#include "TFile.h"
#include "TTree.h"
using namespace std;

#define BUF_LEN 41552
#define STATUS_BUF_LEN 41552
#define MAXTel 25

class StatusDB {
   private:
   static StatusDB* _Head;
   static char DBPath0[2][200];
   static char DBPath[200];

   //for raw data
   ///the filename of the raw data file
   static char FILENAME[300];
   static FILE* fp;
   static uint8_t *buf;
   static unsigned long buflength;

   //for decode data
   static TFile* fstatus;
   static TTree* tree;

   //to locate the filename
   static int currentday;
   static int nfiles[MAXTel];
   static int telindex[MAXTel];
   static int filetime[MAXTel][1000];
   static int fileindex[MAXTel][1000];
   static char filename[MAXTel][1000][200];
   //to read data containt
   static int currentfile;
   static int currenttime;
   static long int entryno;
   static vector<int> failday;
   static vector<int> failtime;

   static WFCTADecode *wfctaDecode;

   public:
   static bool UseDecodeData;
   static int jdebug;
   static int timemargin;
   static int MaxFileTime;

   static Short_t iTel;
   static int fpgaVersion[10];
   static int f9mode;
   static int f9pattern;
   static int DbVersion[2][89];
   static int ClbVersion[2][89];
   static Long64_t clb_initial_Time;
   static double clb_initial_time;
   static int fired_tube;
   static Long64_t status_readback_Time;
   static double status_readback_time;
   static int sipm[1024];
   static int mask[1024];
   static short single_thresh[1024];
   static short record_thresh[1024];
   static Long64_t single_count[1024];
   static Long64_t single_time[1024];
   static float DbTemp[1024];
   static float HV[1024];
   static float PreTemp[1024];
   static float BigResistence[1024];
   static float SmallResistence[1024];
   static Long64_t ClbTime[1024];
   static float ClbTemp[1024];

   static TBranch* b_iTel;
   static TBranch* b_fpgaVersion;
   static TBranch* b_f9mode;
   static TBranch* b_f9pattern;
   static TBranch* b_DbVersion;
   static TBranch* b_ClbVersion;
   static TBranch* b_clb_initial_Time;
   static TBranch* b_clb_initial_time;
   static TBranch* b_fired_tube;
   static TBranch* b_status_readback_Time;
   static TBranch* b_status_readback_time;
   static TBranch* b_sipm;
   static TBranch* b_mask;
   static TBranch* b_single_thresh;
   static TBranch* b_record_thresh;
   static TBranch* b_single_count;
   static TBranch* b_single_time;
   static TBranch* b_DbTemp;
   static TBranch* b_HV;
   static TBranch* b_PreTemp;
   static TBranch* b_BigResistence;
   static TBranch* b_SmallResistence;
   static TBranch* b_ClbTime;
   static TBranch* b_ClbTemp;

   static void Init();
   static void Reset();
   static void Release();
   StatusDB() {Init(); SetDirectory();}
   ~StatusDB(){Release();}
   static void SetDirectory(const char* dirname=0);
   static StatusDB* GetHead();
   static int LocateTel(int iTel0);
   static int LoadFile(int whichday);
   static int LocateFile(int iTel0,int Time,double time=0);
   static bool LocateBlk(int Time,double time=0);
   static void ResetBuffer();
   static void Fill();
   static long int LocateEntry(int Time);
   static long int Locate(int iTel0,int Time,double time=0);

   static int   GetVersion(int iTel0,int Time,int i);
   static long  GetClbInitTime(int iTel0,int Time);
   static double GetClbInittime(int iTel0,int Time);
   static int   GetFiredTube(int iTel0,int Time);
   static long  GetReadbackTime(int iTel0,int Time);
   static double GetReadbacktime(int iTel0,int Time);
   static short GetSingleThrd(int iTel0,int Time,int i);
   static short GetRecordThrd(int iTel0,int Time,int i);
   static long  GetSingleCount(int iTel0,int Time,int i);
   static float GetDbTemp(int iTel0,int Time,int i);
   static long  GetSingleTime(int iTel0,int Time,int i);
   static float GetHV(int iTel0,int Time,int i);
   static float GetPreTemp(int iTel0,int Time,int i);
   static float GetBigR(int iTel0,int Time,int i);
   static float GetSmallR(int iTel0,int Time,int i);
   static long  GetClbTime(int iTel0,int Time,int i);
   static float GetClbTemp(int iTel0,int Time,int i);
};

#endif
