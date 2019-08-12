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
#include <dirent.h>
#include <time.h>
using namespace std;

#define PI 3.14159265358979312
#define RADDEG 180/PI
const double hplank=6.62607015e-34; //in J*s
const double hplank_gev=4.1356676969e-24; //in GeV*s
const double vlight=2.998e10;
  ///convert time
class CommonTools {
   public:
   static bool Is366(int year);
   static int Convert(double time);
   static double InvConvert(int time);
   static int TimeFlag(double time,int type);
   static int TimeFlag(int time,int type);
   static bool GetFirstLastLine(const char* filename,char* firstline,char * lastline);
   static int GetTimeFromFileName(const char* filename,int start,int length);
   static int GetTelIndex(const char* filename,int start,int length);
   static int get_file_size_time(const char* filename);
   static void getFiles(string path,vector<string>& files);
};

/// Max number of Telescopes
#define NCTMax 20
#endif
