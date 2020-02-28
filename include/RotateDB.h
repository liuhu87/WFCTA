#ifndef __RotateDB__
#define __RotateDB__
#define NSWITH 19
#define NVALUE 12

#include "common.h"

class WFCTAEvent;
class RotateDB {
   private:
   static RotateDB* _Head;
   public:
   static int jdebug;
   static int ntotmin;
   static int nsidemin;
   static int nrot;
   static int rotindex[10];
   static int timedelay[10];
   static double rottime[10];
   static int ntel;
   static int telindex[20];
   static double aglmargin;
   static double phimargin;
   static double ccmargin;

   //log version: v1(save the status for every second), v2(save the change of status)
   int version;
   char buff[500];
   char buff2[2][500];
   //Li
   int Li;
   int Li2;
   //time information
   int time;
   int time2[2];
   //position in the log file
   long int currpos;
   long int currpos2;
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
   void CleanBuffer(int ibuf=-1);
   void Release();
   RotateDB() {Init();}
   RotateDB(RotateDB* pr_in);
   ~RotateDB() {Release();}
   void Copy(RotateDB* pr_in);
   static RotateDB* GetHead();
   static bool LocateFirst(ifstream* fin);
   void CleanLog(int ilog);
   long int LoadData(int time_in,int Li_in,int pLi=0,int ptime=-1,long int cpos=-1);
   int ProcessTime();
   void ProcessAll();
   bool IsLogFine(char* buffer,bool IsRotate);
   bool IsLogFine(bool IsRotate=true);
   int ReadData(ifstream* fin,bool godown,bool IsRotate=true);
   int ReadData2(ifstream* fin,bool godown,bool IsRotate=true,int Li_in=2);
   long int LoadData2(int time_in,int Li_in,int pLi=0,int ptime1=0,int ptime2=0,long int cpos=-1);
   int ProcessTime2(int itime);
   void ProcessAll2();
   int ProcessEnv(int time_in,int Li_in);
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
   int GetLi();
   static int GetLi(int Li_in);
   static int GetLi(double rabbittime);
   static int GetTi(int Tindex);
   static double GetMinDistEleAzi(double ele_in,double azi_in,int irot,int itel,double &minele,double &minazi,int &index);
   static int IsFineAngle(double ele_in,double azi_in,int Li_in,int iTel);
   int GetEleAzi1(int time_in,int Li_in,int iTel=-1);
   int GetEleAzi2(int time_in,int Li_in,int iTel=-1);
   int GetEleAzi(int time_in,int Li_in,int iTel=-1);
   int GetEleAzi(WFCTAEvent* pev);
   static void GetMinDistFit(WFCTAEvent* pev,int EleAziIndex,int Li_in,double &minphi,double &mincc);
   static bool IsFineImage(WFCTAEvent* pev,int EleAziIndex,int Li_in=0);
   int LaserIsFine(WFCTAEvent* pev);

   //Get some temperature
   bool GetEnv(int time_in,int Li_in,double *temp);
};

#endif
