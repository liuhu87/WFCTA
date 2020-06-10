#ifndef __RotateDB__
#define __RotateDB__
#define NSWITH 19
#define NVALUE 12

#include "common.h"
#include <fstream>

class WFCTAEvent;
class RotateDB {
   private:
   static RotateDB* _Head;
   static char RotateLogDir[200];
   public:
   static int jdebug;
   static bool UseGPSTime;
   static bool UseGPSInfo;
   static bool UseEnvInfo;
   static bool UsePltInfo;
   static int ntotmin;
   static int nsidemin;
   static int nrot;
   static int rotindex[NRotMax];
   static int timedelay[NRotMax];
   static double rottime[NRotMax];
   static int ntel;
   static int telindex[NCTMax];
   static int TargetIndex;
   static double aglmargin;
   static double phimargin;
   static double ccmargin;
   static double phismargin;
   static double ccsmargin;
   static bool SearchAllAngle;

   //log version: v1(save the status for every second), v2(save the change of status)
   ifstream fin_log1[NRotMax];
   ifstream fin_log2[NRotMax];
   char filename_log1[NRotMax][300];
   char filename_log2[NRotMax][300];
   long int filesize1[NRotMax];
   long int filesize2[NRotMax];

   //about other infos
   ifstream fin_gps[NRotMax];
   ifstream fin_env[NRotMax];
   char filename_gps[NRotMax][300];
   char filename_env[NRotMax][300];

   char rotateinfo[500];
   char rotateinfo2[2][500];

   char pltinfo[2][500];
   char gpsinfo[2][500];
   char envinfo[2][500];

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
   void CleanBuffer(int ibuf=-1);
   void Release();
   RotateDB() {Init();}
   RotateDB(RotateDB* pr_in);
   ~RotateDB() {Release();}
   void Copy(RotateDB* pr_in);
   static char* GetDirName();
   static bool SetDirName(char* dirname);
   static RotateDB* GetHead();
   static RotateDB* GetHead(char* dirname);
   static bool LocateFirst(ifstream* fin);
   static bool LocateEnd(ifstream* fin);
   static bool GetFirstLastLine(ifstream* fin,char* firstline,char* lastline);
   static bool ReadData(ifstream* fin,int godown,char* content);
   static bool ReadAData(ifstream* fin,int godown,char* content,int Type=1,bool cur_godown=true);
   static int GetLogTime(int time_in,int Li_in);
   static bool IsLogFine(char* buffer,int RotateType);
   bool IsLogFine(int RotateType);
   bool FillInfo(int type,char* content1,char* content2);
   int ReadData(int Li_in,int godown=true,char* content=0);
   int ReadAData(int Li_in,int godown,char* content=0,int Type=1,bool cur_godown=true);
   bool LoadData(int time_in,int Li_in);
   bool LoadAData(int time_in,int Li_in,int Type=1);
   bool LoadRotate(int time_in,int Li_in);
   bool LoadGPS(int time_in,int Li_in);
   bool LoadEnv(int time_in,int Li_in);
   bool LoadPlt(int time_in,int Li_in);

   //long int LoadData(int time_in,int Li_in,int pLi=0,int ptime=-1,long int cpos=-1);
   //int ReadData(ifstream* fin,bool godown,int RotateType=1);
   //int ReadData2(ifstream* fin,bool godown,bool IsRotate=true,int Li_in=2);
   //long int LoadData2(int time_in,int Li_in,int pLi=0,int ptime1=0,int ptime2=0,long int cpos=-1);

   static int ProcessTime(const char* content);
   int ProcessTime();
   bool ProcessAll(int time_in,int Li_in,bool update=true);
   static int ProcessTime2(char* content);
   int ProcessTime2(int Type,int itime);
   bool ProcessAll2(int time_in,int Li_in,bool update=true);
   bool ProcessRotate();
   bool ProcessGPS();
   bool ProcessEnv();
   bool ProcessPlant();

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
   static bool GetEleAzi(int index,double &elevation,double &azimuth);
   static double GetMinDistEleAzi(double ele_in,double azi_in,int irot,int itel,double &minele,double &minazi,int &index);
   static int IsFineAngle(double ele_in,double azi_in,int Li_in,int iTel);
   static int SearchVersion(int time_in,int Li_in);
   int GetEleAzi1(int time_in,int Li_in,int iTel=-1);
   int GetEleAzi2(int time_in,int Li_in,int iTel=-1);
   int GetEleAzi(int time_in,int Li_in,int iTel=-1);
   int GetEleAzi(WFCTAEvent* pev);
   static void GetMinDistFit(WFCTAEvent* pev,double ele_in,double azi_in,int Li_in,double &minphi,double &mincc,double &minphi_sigma,double &mincc_sigma);
   static bool IsFineImage(WFCTAEvent* pev,double ele_in,double azi_in);
   static bool IsFineImage(WFCTAEvent* pev,int EleAziIndex);
   static double GetEref(int Li_in,int iTel,int type,int iangle);
   int LaserIsFine(WFCTAEvent* pev);

   //Get some temperature
   bool GetEnv(int time_in,int Li_in,double *temp);
};

#endif
