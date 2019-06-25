#ifndef __COMMON__
#define __COMMON__
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
#endif
