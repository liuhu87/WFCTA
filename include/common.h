#ifndef __COMMON__
#define __COMMON__
#define PI 3.14159265358979312
#define RADDEG 180/PI
const double hplank=6.62607015e-34;
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
   static int GetTimeFromFileName(char* filename,int start,int length);
};

/// Max number of Telescopes
#define NCTMax 20
#endif
