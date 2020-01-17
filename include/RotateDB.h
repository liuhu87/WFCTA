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
   int GetLi();
   static int GetLi(int Li_in);
   static int GetLi(double rabbittime);
   static int GetTi(int Tindex);
   static double GetMinDistEleAzi(double ele_in,double azi_in,int irot,int itel,double &minele,double &minazi,int &index);
   static int IsFineAngle(double ele_in,double azi_in,int Li_in,int iTel);
   int GetEleAzi(int time_in,int Li_in,int iTel=-1);
   int GetEleAzi(WFCTAEvent* pev);
   static void GetMinDistFit(WFCTAEvent* pev,int EleAziIndex,int Li_in,double &minphi,double &mincc);
   static bool IsFineImage(WFCTAEvent* pev,int EleAziIndex,int Li_in=0);
   int LaserIsFine(WFCTAEvent* pev);
};

#endif
