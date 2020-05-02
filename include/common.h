#ifndef __COMMON__
#define __COMMON__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
//#include "/cvmfs/lhaaso.ihep.ac.cn/anysw/slc5_ia64_gcc73/external/root/6.14.00/etc/cling/lib/clang/5.0.0/include/bits/stat.h"
#include <dirent.h>
#include <time.h>
#include "TH1F.h"
using namespace std;
using std::vector;

#define PI 3.14159265358979312
#define RADDEG 180/PI
const double hplank=6.62607015e-34; //in J*s
const double hplank_gev=4.1356676969e-24; //in GeV*s
const double vlight=2.998e10;
const int MJD19700101=40587;
const int TAI2UTC=37;
extern int Nuse;

#define NCTMax 20
#define NRotMax 5
#define MAXPMT 1024

#define InfPNumber 1.0e20
#define InfNNumber (-1.0e20)

  ///convert time
class CommonTools {
   public:
   ///cerenkov photon arrival information
   static TH1F* HArrival[1024];
   ///time bin width,in unit of nano second
   static double timebinunit[2];
   ///check wheather it is laser event
   static bool IsLaser;

   public:
   static bool Is366(int year);
   static int Convert(double time);
   static double InvConvert(int time);
   static double ConvertMJD(double mjdtime);
   static int ConvertMJD2Time(double mjdtime);
   static int ConvertMJD2time(double mjdtime);
   static double InvConvertMJD(int time);
   static int TimeFlag(double time,int type);
   static int TimeFlag(int time,int type);
   static int GetTelDay(int time);
   static int GetFirstLastLine(const char* filename,char* firstline,char * lastline);
   static void GetFileType(char* Type,const char* filename);
   static int GetTimeFromFileName(const char* filename);
   static int GetTimeFromFileName(const char* filename,int start,int length);
   static bool GetStatusFile(char* statusfile,char* eventfile);
   static int GetBins(int start,int end,double step,double bins[100000]);
   static int GetTelIndex(const char* filename,int start,int length);
   static int GetTelIndex(const char* filename);
   static int get_file_size_time(const char* filename);
   static void getFiles(string path,vector<string>& files);

   static void InitHArrival();
   static void ResetHArrival();

   static double ProcessAngle(double angle,bool IsDegree=false);
   static bool CombineAngleRange(double range1[2],double range2[2],double combrange[2],bool IsDegree=false);
};

/// Max number of Telescopes
#define NCTMax 20

/// Max number of time bins
#define MaxTimeBin 150

#endif
