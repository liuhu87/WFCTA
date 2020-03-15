#include "CloudDB.h"
CloudDB* CloudDB::_Head=0;
int CloudDB::jdebug=0;
void CloudDB::Init(){
   Cloud::SetBins();
   cloud.Reset();
   ctime=0;
}
CloudDB* CloudDB::GetHead(){
   if(!_Head) _Head=new CloudDB();
   return _Head;
}
bool CloudDB::LoadIBFile(int Time){
   double time0=ctime*100.;
   int Time0=CommonTools::Convert(time0);
   int time1=Time0-5*60;
   int time2=Time0+5*60;
   if(Time>=time1&&Time<time2){ //still in current IB file
      return true;
   }
   else{
      bool exist=false;
      char filename1[400]="";
      int year,month,day,hour,min;
      int year1,month1,day1;
      int year2,month2,day2;
      for(int ii=-4;ii<=5;ii++){
         int Timei=Time+ii*60;
         year=(CommonTools::TimeFlag(Timei,1)%100)+2000;
         month=CommonTools::TimeFlag(Timei,2);
         day=CommonTools::TimeFlag(Timei,3);
         hour=CommonTools::TimeFlag(Timei,4);
         min=CommonTools::TimeFlag(Timei,5);
         bool before2005=(hour*100+min)<2005;
         int Time1=Timei+24*3600;
         int Time2=Timei+48*3600;
         int xyear=(CommonTools::TimeFlag(Time1,1)%100)+2000;
         int xmonth=CommonTools::TimeFlag(Time1,2);
         int xday=CommonTools::TimeFlag(Time1,3);
         int xyear2=(CommonTools::TimeFlag(Time2,1)%100)+2000;
         int xmonth2=CommonTools::TimeFlag(Time2,2);
         int xday2=CommonTools::TimeFlag(Time2,3);

         year1=before2005?year:xyear;
         month1=before2005?month:xmonth;
         day1=before2005?day:xday;
         year2=before2005?xyear:xyear2;
         month2=before2005?xmonth:xmonth2;
         day2=before2005?xday:xday2;

         char filenamei[400]="";
         strcpy(filenamei,Form("/eos/lhaaso/raw/wfctalaser/IBTemp/Cloudmapdata/%4d/%02d/%4d_%02d%02d/12345_cloud_ir_%4d%02d%02d%02d%02d.dat",year1,month1,year1,month1,day1,year,month,day,hour,min));
         FILE* fin=fopen(filenamei,"r");
         if(jdebug>2) printf("CloudDB::LoadIBFile: read loop ii=%d Time=%d filename=%s FILE=%p\n",ii,Time,filenamei,fin);
         if(fin){
            fclose(fin);
            exist=true;
            strcpy(filename1,filenamei);
            break;
         }
      }
      if(!exist) return false;
      else{
         ctime=(year%100)*((long int)100000000)+month*1000000+day*10000+hour*100+min;
         cloud.ReadCloudMap(filename1);
         char filename2[400]="";
         strcpy(filename2,Form("/eos/lhaaso/raw/wfctalaser/IBTemp/TempHumdata/%4d/%02d/temp/12345_cloud_temp_%4d%02d%02d.txt",year2,month2,year2,month2,day2));
         bool readed=cloud.ReadTemp(filename2);
         if(jdebug>1) printf("CloudDB::LoadIBFile: read humi: Time=%d readed=%d humi=%.2lf(temp=%.2lf,%.2lf) filename=%s\n",Time,readed,cloud.humi,cloud.temp0,cloud.temp,filename2);
         return true;
      }
   }
}

double CloudDB::GetTemperature(int Time,int itemp){
   double temp=1000;
   if(!LoadIBFile(Time)) return temp;
   return cloud.GetTemperature(itemp);
}
double CloudDB::GetHumidity(int Time){
   double humi=-1;
   if(!LoadIBFile(Time)) return humi;
   return cloud.GetHumidity();
}
double CloudDB::GetIBTemp(int Time,int ibin){
   double temp=1000;
   if(!LoadIBFile(Time)) return temp;
   return cloud.GetIBTemp(ibin);
}
double CloudDB::GetIBTemp(int Time,double xx,double yy){
   double temp=1000;
   if(!LoadIBFile(Time)) return temp;
   return cloud.GetIBTemp(xx,yy);
}
double CloudDB::GetTelAveIBTemp(int Time,int iTel){
   double temp=1000;
   if(!LoadIBFile(Time)) return temp;
   return cloud.GetTelAveIBTemp(iTel);
}
double CloudDB::GetTelMinIBTemp(int Time,int iTel){
   double temp=1000;
   if(!LoadIBFile(Time)) return temp;
   return cloud.GetTelMinIBTemp(iTel);
}
double CloudDB::GetTelRmsIBTemp(int Time,int iTel){
   double temp=-1;
   if(!LoadIBFile(Time)) return -1;
   return cloud.GetTelRmsIBTemp(iTel);
}
double CloudDB::GetAveIBTemp(int Time,double theta){
   double temp=1000;
   if(!LoadIBFile(Time)) return temp;
   return cloud.GetAveIBTemp(theta);
}
double CloudDB::GetMinIBTemp(int Time,double theta){
   double temp=1000;
   if(!LoadIBFile(Time)) return temp;
   return cloud.GetMinIBTemp(theta);
}
double CloudDB::GetRmsIBTemp(int Time,double theta){
   double temp=-1;
   if(!LoadIBFile(Time)) return temp;
   return cloud.GetRmsIBTemp(theta);
}
double CloudDB::GetAveIBTemp(int Time,double theta1,double theta2){
   double temp=1000;
   if(!LoadIBFile(Time)) return temp;
   return cloud.GetAveIBTemp(theta1,theta2);
}
double CloudDB::GetMinIBTemp(int Time,double theta1,double theta2){
   double temp=1000;
   if(!LoadIBFile(Time)) return temp;
   return cloud.GetMinIBTemp(theta1,theta2);
}
double CloudDB::GetRmsIBTemp(int Time,double theta1,double theta2){
   double temp=-1;
   if(!LoadIBFile(Time)) return temp;
   return cloud.GetRmsIBTemp(theta1,theta2);
}
double CloudDB::GetAveIBTemp(int Time,TGraph* gr){
   double temp=1000;
   if(!LoadIBFile(Time)) return temp;
   return cloud.GetAveIBTemp(gr);
}
double CloudDB::GetMinIBTemp(int Time,TGraph* gr){
   double temp=1000;
   if(!LoadIBFile(Time)) return temp;
   return cloud.GetMinIBTemp(gr);
}
double CloudDB::GetRmsIBTemp(int Time,TGraph* gr){
   double temp=1000;
   if(!LoadIBFile(Time)) return temp;
   return cloud.GetRmsIBTemp(gr);
}

