#ifndef __CloudDB__
#define __CloudDB__
#include "Cloud.h"
#include "common.h"
class CloudDB{
   private:
   static CloudDB* _Head;
   public:
   static int jdebug;
   Cloud cloud;
   long int ctime;
   void Init();
   void Clear() { ctime=0;}
   CloudDB() {Init();}
   ~CloudDB() {Clear();}
   static CloudDB* GetHead();
   bool LoadIBFile(int Time);
   double GetTemperature(int Time,int itemp);
   double GetHumidity(int Time);
   double GetIBTemp(int Time,int ibin);
   double GetIBTemp(int Time,double xx,double yy);
   double GetTelAveIBTemp(int Time,int iTel);
   double GetTelMinIBTemp(int Time,int iTel);
   double GetTelRmsIBTemp(int Time,int iTel);
   double GetAveIBTemp(int Time,double theta);
   double GetMinIBTemp(int Time,double theta);
   double GetRmsIBTemp(int Time,double theta);
   double GetAveIBTemp(int Time,double theta1,double theta2);
   double GetMinIBTemp(int Time,double theta1,double theta2);
   double GetRmsIBTemp(int Time,double theta1,double theta2);
   double GetAveIBTemp(int Time,TGraph* gr);
   double GetMinIBTemp(int Time,TGraph* gr);
   double GetRmsIBTemp(int Time,TGraph* gr);
};
#endif
