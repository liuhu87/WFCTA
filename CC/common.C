#include "common.h"
#include <fstream>
#include "string.h"
#include "stdlib.h"
using namespace std;
///check wheather the year is leap year
bool CommonTools::Is366(int year){
   if((year%4==0)&&(year%100!=0)) return true;
   else if((year%100==0)&&(year%400==0)) return true;
   return false;
}
///convert the time to second since 1970-01-01 08:00:00
int CommonTools::Convert(double time){
   int year=int(time/10000000000);
   if(year>=100) year=year%100;
   double time1=(time-year*10000000000);
   int month=int(time1/100000000);
   int time2=int(time1-month*100000000);
   int day=(time2%100000000)/1000000;
   int h=(time2%1000000)/10000;
   int min=(time2%10000)/100;
   int sec=time2%100;
   int result=0;
   for(int ii=0;ii<year;ii++){
      result+=(Is366(2000+ii)?366:365)*24*3600;
   }
   for(int ii=1;ii<month;ii++){
      if((ii==1||ii==3||ii==5||ii==7)||(ii==8||ii==10||ii==12)) result+=31*24*3600;
      else if(ii==2) result+=(Is366(2000+ii)?29:28)*24*3600;
      else result+=30*24*3600;
   }
   for(int ii=1;ii<day;ii++){
      result+=24*3600;
   }
   for(int ii=0;ii<h;ii++){
      result+=3600;
   }
   for(int ii=0;ii<min;ii++){
      result+=60;
   }
   result+=sec;
   return result+946656000; //946656000 is the number for "2000-01-01 00:00:00"
}
//inverse of convert
double CommonTools::InvConvert(int time){
   int time0=time-946656000;
   if(time0<0) return 0;
   int year,month,day,h,min,sec;
   year=0;
   while(true){
      int result=(Is366(2000+year)?366:365)*24*3600;
      if(time0<result) break;
      else{time0-=result; year++;}
   }
   int result=0;
   for(int ii=1;ii<=12;ii++){
      int result0=result;
      if((ii==1||ii==3||ii==5||ii==7)||(ii==8||ii==10||ii==12)) result+=31*24*3600;
      else if(ii==2) result+=(Is366(2000+ii)?29:28)*24*3600;
      else result+=30*24*3600;
      if(time0<result) {month=ii; time0-=result0; break;}
   }
   day=time0/(24*3600);
   time0-=day*(24*3600);
   day+=1;
   h=time0/3600;
   time0-=h*3600;
   min=time0/60;
   time0-=min*60;
   sec=time0;

   double res=0;
   res+=year*10000000000;
   res+=month*100000000;
   res+=day*1000000;
   res+=h*10000;
   res+=min*100;
   res+=sec;
   return res;
}
///return the year/month/day/hour/minute/second
int CommonTools::TimeFlag(double time,int type){
   int year=int(time/10000000000);
   if(year>=100) year=year%100;
   double time1=(time-year*10000000000);
   int month=int(time1/100000000);
   int time2=int(time1-month*100000000);
   int day=(time2%100000000)/1000000;
   int h=(time2%1000000)/10000;
   int min=(time2%10000)/100;
   int sec=time2%100;
   if(type==1) return year;
   else if(type==2) return month;
   else if(type==3) return day;
   else if(type==4) return h;
   else if(type==5) return min;
   else if(type==6) return sec;
   else return -1;
}
int CommonTools::TimeFlag(int time,int type){
   double time0=InvConvert(time);
   return TimeFlag(time0,type);
}
bool CommonTools::GetFirstLastLine(const char* filename,char* firstline,char * lastline){
   const int maxlen=300;
   char buff[maxlen];
   char timebuff[13];
   long pos0,pos1;
   ifstream fbuff(filename,std::ios::in);
   fbuff.seekg(ios::beg);
   pos0=fbuff.tellg();
   fbuff.getline(buff,maxlen);
   pos1=fbuff.tellg();
   bool res=(!fbuff.eof());
   int length=strlen(buff);
   strcpy(firstline,buff);
   fbuff.seekg(pos0-pos1,ios::end);
   fbuff.getline(buff,maxlen);
   strcpy(lastline,buff);
   return res;
}
int CommonTools::GetTimeFromFileName(const char* filename,int start,int length){
   int strlength=strlen(filename);
   if(start<0) return 0;
   else if(start+length>strlength) return 0;
   else{
      char timebuff[13];
      int np=0;
      for(int ii=start;ii<start+length;ii++){
         timebuff[np++]=filename[ii];
      }
      for(int ii=np;ii<12;ii++){
         timebuff[ii]='0';
         np++;
      }
      timebuff[np++]='\0';
      //printf("name=%s  timebuff=%s(np=%d)\n",filename,timebuff,np);
      int time=Convert(atol(timebuff));
      return time;
   }
}
int CommonTools::GetTelIndex(const char* filename,int start,int length){
   int strlength=strlen(filename);
   if(start<0) return 0;
   else if(start+length>strlength) return 0;
   else{
      char timebuff[100];
      int np=0;
      for(int ii=start;ii<start+length;ii++){
         timebuff[np++]=filename[ii];
      }
      timebuff[np++]='\0';
      int itel=atoi(timebuff);
      int time=GetTimeFromFileName(filename,46,12);
      if(time<Convert(190517160000)) itel+=3;
      return itel;
   }
}

int CommonTools::get_file_size_time(const char* filename){
   struct stat statbuf;
   if(stat(filename,&statbuf) == -1) return -1;
   if(S_ISDIR(statbuf.st_mode)) return 1;
   if(S_ISREG(statbuf.st_mode)) return 0;
}
void CommonTools::getFiles(string path,vector<string>& files){
   files.clear();
   string end(".dat");

   DIR* dirp;
   struct dirent *direntp;
   int status;
   char buf[300];
   if(((status = get_file_size_time(path.data())) == 0) || (status == -1) ) return;
   if( (dirp=opendir(path.data())) == NULL) return;
   while((direntp=readdir(dirp))!=NULL){
      sprintf(buf,"%s/%s",path.data(),direntp->d_name);
      //printf("read %s\n",buf);
      if(get_file_size_time(buf)==-1) break;
      string filename(buf);
      if( filename.length()>=end.length() && end == filename.substr(filename.length()-end.length(),end.length())){
         files.push_back(filename);
      }
      else continue;
   }
   closedir(dirp);
}
