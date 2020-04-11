#ifndef __Status_DB__
#define __Status_DB__
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include "WFCTADecode.h"
#include "TFile.h"
#include "TTree.h"
#include "common.h"
using namespace std;

#define BUF_LEN 41552
#define STATUS_BUF_LEN 41552

class StatusDB {
   private:
   static StatusDB* _Head;
   static char DBPath0[2][200];
   static char DBPath[200];

   //for raw data
   ///the filename of the raw data file
   char FILENAME[300];
   FILE* fp;
   uint8_t *buf;
   unsigned long buflength;

   //for decode data
   TFile* fstatus;
   TTree* tree;

   //to locate the filename
   int currentday;
   int nfiles[NCTMax];
   int telindex[NCTMax];
   int filetime[NCTMax][1000];
   int fileindex[NCTMax][1000];
   char filename[NCTMax][1000][200];
   //to read data containt
   int currentfile;
   int currenttime;
   long int entryno;
   vector<int> failday;
   vector<int> failtime;

   WFCTADecode *wfctaDecode;

   public:
   static bool UseDecodeData;
   static int jdebug;
   static int timemargin;
   static int MaxFileTime;

   Short_t iTel;
   int fpgaVersion[10];
   int f9mode;
   int f9pattern;
   int DbVersion[2][89];
   int ClbVersion[2][89];
   Long64_t clb_initial_Time;
   double clb_initial_time;
   int fired_tube;
   Long64_t status_readback_Time;
   double status_readback_time;
   int sipm[MAXPMT];
   int mask[MAXPMT];
   short single_thresh[MAXPMT];
   short record_thresh[MAXPMT];
   Long64_t single_count[MAXPMT];
   Long64_t single_time[MAXPMT];
   float DbTemp[MAXPMT];
   float HV[MAXPMT];
   float PreTemp[MAXPMT];
   float BigResistence[MAXPMT];
   float SmallResistence[MAXPMT];
   Long64_t ClbTime[MAXPMT];
   float ClbTemp[MAXPMT];

   void Init();
   void Reset();
   void Release();
   StatusDB() {Init(); SetDirectory();}
   ~StatusDB(){Release();}
   static void SetDirectory(const char* dirname=0);
   static StatusDB* GetHead();
   int LocateTel(int iTel0);
   int LoadFile(int whichday,char* filename_in=0);
   int LocateFile(int iTel0,int Time,char* filename_in=0,double time=0);
   bool LocateBlk(int Time,double time=0);
   void ResetBuffer();
   void Fill();
   long int LocateEntry(int Time);
   long int Locate(int iTel0,int Time,char* filename_in=0,double time=0);

   int   GetVersion(int iTel0,int Time,int i,char* filename_in=0);
   long  GetClbInitTime(int iTel0,int Time,char* filename_in=0);
   double GetClbInittime(int iTel0,int Time,char* filename_in=0);
   int   GetFiredTube(int iTel0,int Time,char* filename_in=0);
   long  GetReadbackTime(int iTel0,int Time,char* filename_in=0);
   double GetReadbacktime(int iTel0,int Time,char* filename_in=0);
   short GetSingleThrd(int iTel0,int Time,int i,char* filename_in=0);
   short GetRecordThrd(int iTel0,int Time,int i,char* filename_in=0);
   long  GetSingleCount(int iTel0,int Time,int i,char* filename_in=0);
   float GetDbTemp(int iTel0,int Time,int i,char* filename_in=0);
   long  GetSingleTime(int iTel0,int Time,int i,char* filename_in=0);
   float GetHV(int iTel0,int Time,int i,char* filename_in=0);
   float GetPreTemp(int iTel0,int Time,int i,char* filename_in=0);
   float GetBigR(int iTel0,int Time,int i,char* filename_in=0);
   float GetSmallR(int iTel0,int Time,int i,char* filename_in=0);
   long  GetClbTime(int iTel0,int Time,int i,char* filename_in=0);
   float GetClbTemp(int iTel0,int Time,int i,char* filename_in=0);
};

#endif
