#include "common.h"
#include <fstream>
#include "string.h"
#include "stdlib.h"
using namespace std;

int Nuse=0;
TH1F* CommonTools::HArrival[1024]={0};
double CommonTools::timebinunit[2]={0.3,3000};
bool CommonTools::IsLaser=false;
///check wheather the year is leap year
bool CommonTools::Is366(int year){
   if(((year%4)==0)&&((year%100)!=0)) return true;
   else if(((year%100)==0)&&((year%400)==0)) return true;
   return false;
}
///convert the time to second since 1970-01-01 08:00:00
int CommonTools::Convert(double time){
   int year=int(time*1.e-10);
   double time1=(time-year*10000000000);
   year=(year%100);
   int month=int(time1/100000000);
   int time2=int(time1-month*100000000);
   int day=(time2%100000000)/1000000;
   int h=(time2%1000000)/10000;
   int min=(time2%10000)/100;
   int sec=(time2%100);
   //printf("time=%.0lf %d-%d-%d %d:%d:%d\n",time,year,month,day,h,min,sec);
   int result=0;
   for(int ii=0;ii<year;ii++){
      result+=(Is366(2000+ii)?366:365)*24*3600;
      //printf("add year: ii=%d day=%d\n",ii,(Is366(2000+ii)?366:365));
   }
   //printf("year,res=%d\n",result);
   for(int ii=1;ii<month;ii++){
      if((ii==1||ii==3||ii==5||ii==7)||(ii==8||ii==10||ii==12)) result+=31*24*3600;
      else if(ii==2) result+=(Is366(2000+year)?29:28)*24*3600;
      else result+=30*24*3600;
   }
   //printf("month,res=%d\n",result);
   for(int ii=1;ii<day;ii++){
      result+=24*3600;
   }
   //printf("day,res=%d\n",result);
   for(int ii=0;ii<h;ii++){
      result+=3600;
   }
   //printf("hour,res=%d\n",result);
   for(int ii=0;ii<min;ii++){
      result+=60;
   }
   //printf("minute,res=%d\n",result);
   result+=sec;
   //printf("second,res=%d\n",result);
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
      else if(ii==2) result+=(Is366(2000+year)?29:28)*24*3600;
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
   res+=(year+2000.)*10000000000.;
   res+=month*100000000;
   res+=day*1000000;
   res+=h*10000;
   res+=min*100;
   res+=sec;
   return res;
}

double CommonTools::ConvertMJD(double mjdtime){
  double time = (( mjdtime - MJD19700101 ) * 86400 + TAI2UTC);
  return time;
}
int CommonTools::ConvertMJD2Time(double mjdtime){
  return ((int)ConvertMJD(mjdtime));
}
int CommonTools::ConvertMJD2time(double mjdtime){
  double time = ConvertMJD(mjdtime);
  int Time=((int)time);
  return (int)((time-Time)/(2.0e-8));
}
double CommonTools::InvConvertMJD(int time){
   double mjd=MJD19700101 + (time - TAI2UTC)/86400.;
   return mjd;
}
///return the year/month/day/hour/minute/second
int CommonTools::TimeFlag(double time,int type){
   int year=int(time/10000000000);
   double time1=(time-year*10000000000);
   year=(year%100)+2000;
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
int CommonTools::GetTelDay(int time){
   int hour=TimeFlag(time,4);
   int min=TimeFlag(time,5);
   int sec=TimeFlag(time,6);
   int year=(TimeFlag(time,1)%100)+2000;
   int month=TimeFlag(time,2);
   int day=TimeFlag(time,3);
   int pyear=(TimeFlag(time-24*3600,1)%100)+2000;
   int pmonth=TimeFlag(time-24*3600,2);
   int pday=TimeFlag(time-24*3600,3);
   if(hour*10000+min*100+sec<120000) return pyear*10000+pmonth*100+pday;
   else return year*10000+month*100+day;
}
int CommonTools::GetFirstLastLine(const char* filename,char* firstline,char * lastline){
   if(!filename) return 0;
   const int maxlen=300;
   char buff[maxlen];
   char buff2[maxlen];
   /*char timebuff[13];
   long pos0,pos1,pos2;
   ifstream fbuff(filename,std::ios::in);
   fbuff.seekg(ios::beg);
   pos0=fbuff.tellg();
   fbuff.getline(buff,maxlen);
   pos1=fbuff.tellg();
   int length=strlen(buff);
   strcpy(firstline,buff);
   fbuff.seekg(3*(pos0-pos1)/2,ios::end);
   fbuff.getline(buff,maxlen);
   fbuff.getline(buff,maxlen);
   pos2=fbuff.tellg();
   strcpy(lastline,buff);
   return (pos1>pos0)?(pos2-pos0)/(pos1-pos0):0;*/

   int count=0;
   ifstream fbuff(filename,std::ios::in);
   if(!fbuff.is_open()) return 0;
   do{
      if(count==1) strcpy(firstline,buff);
      for(int ii=0;ii<maxlen;ii++){
         buff2[ii]=buff[ii];
         if(buff[ii]=='\0'||buff[ii]=='\n') break;
      }
      fbuff.getline(buff,maxlen);
      count++;
   }
   while(fbuff.good());
   fbuff.close();
   if(count==1) return 0;
   else strcpy(lastline,buff2);
   return count-1;
}
void CommonTools::GetFileType(char* Type,const char* filename){
   if(!Type) return;
   int length=strlen(filename);
   int iloc=-1;
   for(int ii=length-1;ii>=0;ii--){
      if(filename[ii]=='.'){ iloc=ii; break;}
   }
   if(length<2||iloc<0) {Type[0]='\0'; return;}
   for(int ii=0;ii<length;ii++){
      if(ii>iloc) Type[ii-iloc-1]=filename[ii];
   }
   Type[length-1-iloc]='\0';
}
int CommonTools::GetTimeFromFileName(const char* filename){
   int strlength=strlen(filename);
   int ncount=0,nchar=0;
   char buff0[300];
   bool exist=false;
   for(int ii=strlength;ii>=0;ii--){
      //printf("ii=%d(%d) ss=%c\n",ii,strlength,filename[ii]);
      if(filename[ii]=='.'||filename[ii]=='_'){
         ncount++;
         if(ncount>1) buff0[nchar++]='\0';
         bool istime=true;
         int nnum=0;
         for(int jj=nchar-2;jj>=0;jj--){
            int index=(buff0[jj]-'0');
            //printf("jj=%d char=%c index=%d\n",jj,buff0[jj],index);
            if(index<0||index>9) {istime=false; break;}
            nnum++;
         }
         istime=istime&&(nnum>=12);
         //printf("ncount=%d buff0=%s istime=%d(%d)\n",ncount,buff0,istime,nnum);
         if(istime){
            exist=true;
            break;
         }
         nchar=0; continue;
      }
      buff0[nchar++]=filename[ii];
   }
   if(!exist) return 0;
   char timebuff[300];
   for(int ii=nchar-2;ii>=0;ii--) timebuff[(nchar-2)-ii]=buff0[ii];
   for(int ii=nchar-1;ii<=14;ii++){
      timebuff[ii]=(ii<14)?'0':'\0';
   }
   return Convert(atol(timebuff));
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
bool CommonTools::GetStatusFile(char* statusfile,char* eventfile){
   if(!eventfile) return false;
   char* pos=strstr(eventfile,".event.root");
   if(!pos) return false;
   int nchar=(pos-eventfile);
   if(nchar<=0) return false;
   for(int ii=0;ii<nchar;ii++){
      statusfile[ii]=eventfile[ii];
   }
   statusfile[nchar++]='.';
   statusfile[nchar++]='s';
   statusfile[nchar++]='t';
   statusfile[nchar++]='a';
   statusfile[nchar++]='t';
   statusfile[nchar++]='u';
   statusfile[nchar++]='s';
   statusfile[nchar++]='.';
   statusfile[nchar++]='r';
   statusfile[nchar++]='o';
   statusfile[nchar++]='o';
   statusfile[nchar++]='t';
   statusfile[nchar++]='\0';
   return true;
}
int CommonTools::GetBins(int start,int end,double step,double bins[100000]){
   const int maxbins=100000;
   int h1=8*3600,h2=17*3600+30*60;
   int nbin=0;
   int hour1=TimeFlag(start,4)*3600+TimeFlag(start,5)*60+TimeFlag(start,6);
   int hour2=TimeFlag(end,4)*3600+TimeFlag(end,5)*60+TimeFlag(end,6);
   bins[nbin++]=start;
   while(nbin<100000){
      double newbin=bins[nbin-1]+step;
      int hour=TimeFlag((int)(newbin+0.5),4)*3600+TimeFlag((int)(newbin+0.5),5)*60+TimeFlag((int)(newbin+0.5),6);
      if(hour>=h1&&hour<=h2){
         newbin+=(h2-hour);
      }
      bins[nbin++]=(newbin>=end)?end:newbin;
      //int bintime=(int)(bins[nbin-1]+0.5);
      //printf("nbin=%d bins=%lf %04d-%02d-%02d %02d:%02d:%02d\n",nbin-1,bins[nbin-1],CommonTools::TimeFlag(bintime,1),CommonTools::TimeFlag(bintime,2),CommonTools::TimeFlag(bintime,3),CommonTools::TimeFlag(bintime,4),CommonTools::TimeFlag(bintime,5),CommonTools::TimeFlag(bintime,6));
      if(newbin>=end) break;
      if(nbin>=maxbins) break;
   }
   return nbin-1;
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
int CommonTools::GetTelIndex(const char* filename){
   int strlength=strlen(filename);
   int ncount=0,nchar=0;
   char buff0[300];
   int itel=-1;
   bool IsOld=false;
   for(int ii=0;ii<strlength;ii++){
      //printf("ii=%d(%d) ss=%c\n",ii,strlength,filename[ii]);
      if(filename[ii]=='.'||filename[ii]=='_'){
         ncount++;
         buff0[nchar++]='\0';
         //check the strings
         char telname[10]="";
         int ichar=0;
         if(strstr(buff0,"WFCTA")){
            for(int i0=5;i0<nchar;i0++){
               int iindex=(buff0[i0]-'0');
               if((iindex>=0&&iindex<=9)||buff0[i0]=='\0') telname[ichar++]=buff0[i0];
            }
            itel=atoi(telname);
            break;
         }
         if(strstr(buff0,"FULL")){
            itel=100;
         }
         if(strstr(buff0,"Channel")){
            IsOld=true;
         }
         if(IsOld){
            for(int i0=2;i0<nchar;i0++) telname[ichar++]=buff0[i0];
            if(ichar<2) itel=-1;
            else{
               itel=atoi(telname)+100;
               double time=GetTimeFromFileName(filename);
               if(time<Convert(190517160000.)) itel+=3;
            }
            break;
         }
         nchar=0; continue;
      }
      buff0[nchar++]=filename[ii];
   }
   return itel;
}

int CommonTools::get_file_size_time(const char* filename){
   struct stat statbuf;
   if(stat(filename,&statbuf) == -1) return -1;
   if(S_ISDIR(statbuf.st_mode)) return 1;
   if(S_ISREG(statbuf.st_mode)) return 0;
   return -1;
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
      //if( filename.length()>=end.length() && end == filename.substr(filename.length()-end.length(),end.length())){
         files.push_back(filename);
      //}
      //else continue;
   }
   closedir(dirp);
}

void CommonTools::InitHArrival(){
   for(int ii=0;ii<1024;ii++){
      HArrival[ii]=new TH1F(Form("PMT_%d",ii),"",(MaxTimeBin*2),-MaxTimeBin,MaxTimeBin);
   }
}
void CommonTools::ResetHArrival(){
   for(int ii=0;ii<1024;ii++){
      if(HArrival[ii]) HArrival[ii]->Reset();
   }
}

double CommonTools::ProcessAngle(double angle,bool IsDegree){
   double scale=180./PI;
   if(IsDegree) angle=angle/scale;
   double res;
   if(angle>=0&&angle<(2*PI)){
      res=angle;
   }
   else if(angle<0){
      int nn=int((angle)/(-2*PI)-1.0e-8);
      res=angle+(nn+1)*(2*PI);
   }
   else{
      int nn=int((angle)/(2*PI));
      res=angle-nn*(2*PI);
   }
   return IsDegree?(res*scale):res;
}
