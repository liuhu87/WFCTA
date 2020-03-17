#include "RotateDB.h"
#include "stdlib.h"
#include "WFCTAEvent.h"
#include <fstream>
#include <string>
#include <stdio.h>
RotateDB* RotateDB::_Head=0;
int RotateDB::jdebug=0;
int RotateDB::ntotmin=5;
int RotateDB::nsidemin=2;
int RotateDB::nrot=2;
int RotateDB::rotindex[10]={2,3,0,0,0,0,0,0,0,0};
int RotateDB::timedelay[10]={35,37,0,0,0,0,0,0,0,0};
double RotateDB::rottime[10]={990016000,990845000,0,0,0,0,0,0,0,0};
int RotateDB::ntel=6;
int RotateDB::telindex[20]={1,2,3,4,5,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
double RotateDB::aglmargin=0.02;
double RotateDB::phimargin=10.;
double RotateDB::ccmargin=1.;
RotateDB::RotateDB(RotateDB* pr_in){
   Copy(pr_in);
}
void RotateDB::Init(){
   version=0;
   for(int ii=0;ii<500;ii++){
      buff[ii]='\0';
      buff2[0][ii]='\0';
      buff2[1][ii]='\0';
   }
   Li=0;
   Li2=0;
   time=0;
   time2[0]=0;
   time2[1]=0;
   currpos=-1;
   currpos2=-1;
   Reset();
}
void RotateDB::Reset(){
   latitude=0;
   longitude=0;
   altitude=0;
   for(int ii=0;ii<NSWITH;ii++) allswith[ii]=false;
   for(int ii=0;ii<NVALUE;ii++) varinfo[ii]=0;
}
void RotateDB::CleanLog(int ilog){
   if(ilog<0||(ilog&0x1)>0){
      for(int ii=0;ii<500;ii++){
         buff[ii]='\0';
      }
      Li=0;
      time=0;
      currpos=-1;
   }
   if(ilog<0||(ilog&0x2)>0){
      for(int ii=0;ii<500;ii++){
         buff2[0][ii]='\0';
      }
      Li2=0;
      time2[0]=0;
      time2[1]=0;
      currpos2=-1;
   }
}
void RotateDB::Release(){
   Init();
}
void RotateDB::Copy(RotateDB* pr_in){
   version=pr_in->version;
   for(int ii=0;ii<500;ii++){
      buff[ii]=pr_in->buff[ii];
      buff2[0][ii]=pr_in->buff2[0][ii];
      buff2[1][ii]=pr_in->buff2[1][ii];
   }
   Li=pr_in->Li;
   Li2=pr_in->Li2;
   time=pr_in->time;
   time2[0]=pr_in->time2[0];
   time2[1]=pr_in->time2[1];
   currpos=pr_in->currpos;
   latitude=pr_in->latitude;
   longitude=pr_in->longitude;
   altitude=pr_in->altitude;
   for(int ii=0;ii<NSWITH;ii++) allswith[ii]=pr_in->allswith[ii];
   for(int ii=0;ii<NVALUE;ii++) varinfo[ii]=pr_in->varinfo[ii];
}
RotateDB* RotateDB::GetHead(){
   if(!_Head) _Head=new RotateDB();
   return _Head;
}
bool RotateDB::LocateFirst(ifstream* fin){
   if(!fin) return false;
   if(!fin->good()) return false;
   long int pos0=fin->tellg();
   char current='0';
   while(!(current=='\n')){
      long int posi=fin->tellg();
      if(posi==0) return true;
      fin->seekg(-1,ios::cur);
      current=fin->get();
      fin->seekg(-1,ios::cur);
   }
   if(current=='\n'){
      fin->seekg(1,ios::cur);
      return true;
   }
   else return false;
}
long int RotateDB::LoadData(int time_in,int Li_in,int pLi,int ptime,long int cpos){
   int irot=GetLi(Li_in);
   if(irot<0) return -5;
   int pirot=GetLi(pLi);
   int time_new=time_in-timedelay[irot];
   int year=CommonTools::TimeFlag(time_new,1);
   year=2000+(year%100);
   int month=CommonTools::TimeFlag(time_new,2);
   int day=CommonTools::TimeFlag(time_new,3);
   int hour=CommonTools::TimeFlag(time_new,4);
   int min=CommonTools::TimeFlag(time_new,5);
   int sec=CommonTools::TimeFlag(time_new,6);
   char filename[100]="/scratchfs/lhaaso/hliu/rotate_log/";
   char* namebuff=Form("L%d/%d-%02d-%02d.txt.utf8",Li_in,year,month,day);
   strcat(filename,namebuff);
   if(jdebug>0) printf("RotateDB::LoadData: filename=%s\n",filename);
   bool sameday=false;
   if(ptime>0&&pLi>0&&pirot>=0){
      int ptime_new=ptime-timedelay[pirot];
      int pyear=CommonTools::TimeFlag(ptime_new,1);
      pyear=2000+(pyear%100);
      int pmonth=CommonTools::TimeFlag(ptime_new,2);
      int pday=CommonTools::TimeFlag(ptime_new,3);
      sameday=(Li_in==pLi)&&((year==pyear)&&(month==pmonth)&&(day==pday));
   }

   //read the file content
   ifstream fin;
   fin.open(filename,std::ios::in);
   if(!fin.is_open()){
      //Reset();
      return -1;
   }
   fin.seekg(0,std::ios::beg);
   long int poss=fin.tellg();
   fin.getline(buff,500);
   long int pos1=fin.tellg();
   fin.seekg(0,std::ios::end);
   long int pose=fin.tellg();
   long int linelength=pos1-poss;
   if(jdebug>0) printf("RotateDB::LoadData: file opened. pos={%ld,%ld,%ld},linelength=%ld\n",poss,pos1,pose,linelength);

   //locate the initial position
   bool located=false;
   if(sameday&&cpos>=0){
      fin.seekg(cpos,std::ios::beg);
      if(((long int)fin.tellg())<0){
         fin.close();
         fin.open(filename,std::ios::in);
      }
      else located=true;
   }
   if(!located){
      int time_first=abs(ProcessTime());
      if(jdebug>2) printf("RotateDB::LoadData: time_new=%d time_first=%d buff=%s\n",time_new,time_first,buff);
      if(time_first>time_new){
         //Reset();
         fin.close();
         return -3;
      }
      fin.seekg(-10,std::ios::end);
      int time_last;
      if(!LocateFirst(&fin)){
         //Reset();
         fin.close();
         return -2;
      }
      else{
         fin.getline(buff,500);
         time_last=abs(ProcessTime());
         if(jdebug>2) printf("RotateDB::LoadData: time_new=%d time_last=%d buff=%s\n",time_new,time_last,buff);
         if(time_last<time_new){
            //Reset();
            fin.close();
            return -3;
         }
      }

      if(hour<12){
         long int posi=(time_new-time_first)*linelength;
         if(jdebug>2) printf("RotateDB::LoadData: poss=%ld pose=%ld posi=%d deltaT1=%d\n",poss,pose,posi,(hour*3600+min*60+sec-time_first));
         if(poss+posi>=pose){
            fin.seekg(-10,std::ios::end);
         }
         else{
            fin.seekg(posi,std::ios::beg);
         }
      }
      else{
         long int posi=pose-(time_last-time_new)*linelength;
         if(jdebug>2) printf("RotateDB::LoadData: poss=%ld pose=%ld posi=%d deltaT2=%d\n",poss,pose,posi,(time_last-time_new));
         if(posi<poss){
            fin.seekg(0,std::ios::beg);
         }
         else{
            fin.seekg(-(time_last-time_new)*linelength,std::ios::end);
         }
      }
   }

   LocateFirst(&fin);
   if(jdebug>0) printf("RotateDB::LoadData: initial position located. pos=%ld(%d)\n",fin.tellg(),located);
   fin.getline(buff,500);
   long int posc=fin.tellg();
   int timei=abs(ProcessTime());
   if(jdebug>0) printf("RotateDB::LoadData: time_in=%d time_loop0=%d buff_loop0=%s\n",time_new,timei,buff);
   bool godown=timei<time_new;
   int nloop=0;
   while(timei!=time_new){
      if(godown){
         fin.getline(buff,500);
         timei=abs(ProcessTime());
         if(!(fin.good())) break;
         if(timei>=time_new) break;
      }
      else{
         long int posi=posc-(long int)(1.5*linelength);
         if(posi>=poss){
            fin.seekg(-(long int)(1.5*linelength),std::ios::cur);
            LocateFirst(&fin);
            fin.getline(buff,500);
            timei=abs(ProcessTime());
            if(timei<=time_new) break;
         }
         else break;
      }
      if(jdebug>1) printf("RotateDB::LoadData: loop=%d(%d) time={%d,%d} buff=%s\n",nloop,godown,timei,time_new,buff);
      nloop++;
   }

   if(jdebug>0) printf("RotateDB::LoadData: time_in=%d time_final=%d\n",time_new,timei);
   if(timei==time_new){
      fin.seekg(-(long int)(0.5*linelength),std::ios::cur);
      LocateFirst(&fin);
      long int pos0=fin.tellg();
      version=1;
      Li=Li_in;
      time=time_in;
      currpos=pos0;
      fin.close();
      if(jdebug>0) printf("RotateDB::LoadData: loaded, Li=%d time=%d cpos=%ld\n",Li,time,currpos);
      return currpos;
   }
   else{
      //Reset();
      fin.close();
      return -4;
   }
}
int RotateDB::ProcessTime(){
   int size=strlen(buff);
   char buff0[10];

   int count;
   int nchar1,nchar2,nchar3;
   //Get GPS Time information
   char hour[10],min[10],sec[10];
   for(int i1=0;i1<10;i1++){
      hour[i1]='\0';
      min[i1]='\0';
      sec[i1]='\0';
   }
   int index_time1=0;
   int index_time2=0;
   count=0;
   nchar1=0;nchar2=0;nchar3=0;
   for(int ii=size-1;ii>=0;ii--){
      if(buff[ii]==':'){
         count++;
         if(count==1) index_time2=ii+3;
         else if(count==2) index_time1=ii-4;
         continue;
      }
      if(count==1){
         min[nchar1++]=buff[ii];
      }
      else if(count==2){
         if(buff[ii]==' ') break;
         hour[nchar2++]=buff[ii];
      }
   }
   //reverse
   for(int ii=nchar1-1;ii>=0;ii--) buff0[nchar1-1-ii]=min[ii];
   for(int ii=0;ii<nchar1;ii++) min[ii]=buff0[ii];
   for(int ii=nchar2-1;ii>=0;ii--) buff0[nchar2-1-ii]=hour[ii];
   for(int ii=0;ii<nchar2;ii++) hour[ii]=buff0[ii];
   if(index_time2<=size) {sec[0]=buff[index_time2-2]; sec[1]=buff[index_time2-1];}

   char year[10],month[10],day[10];
   for(int i1=0;i1<10;i1++){
      year[i1]='\0';
      month[i1]='\0';
      day[i1]='\0';
   }
   count=0;
   nchar1=0;nchar2=0;nchar3=0;
   for(int ii=index_time1;ii>=0;ii--){
      if(buff[ii]=='-'){
         count++;
         continue;
      }
      if(count==0){
         day[nchar1++]=buff[ii];
      }
      if(count==1){
         month[nchar2++]=buff[ii];
      }
      else if(count==2){
         if(nchar3>=4) break;
         year[nchar3++]=buff[ii];
      }
   }
   //reverse
   for(int ii=nchar1-1;ii>=0;ii--) buff0[nchar1-1-ii]=day[ii];
   for(int ii=0;ii<nchar1;ii++) day[ii]=buff0[ii];
   for(int ii=nchar2-1;ii>=0;ii--) buff0[nchar2-1-ii]=month[ii];
   for(int ii=0;ii<nchar2;ii++) month[ii]=buff0[ii];
   for(int ii=nchar3-1;ii>=0;ii--) buff0[nchar3-1-ii]=year[ii];
   for(int ii=0;ii<nchar3;ii++) year[ii]=buff0[ii];

   double time0=atoi(year)*10000000000+atoi(month)*100000000+atoi(day)*1000000+atoi(hour)*10000+atoi(min)*100+atoi(sec);
   int result=CommonTools::Convert(time0);
   return (index_time1>125)?result:(-result);
}
void RotateDB::ProcessAll(){
   if(time<=0) return;
   int size=strlen(buff);
   if(jdebug>0) printf("RotateDB::ProcessAll: %s(size=%d)\n",buff,size);
   char buff0[10];

   int count;
   int nchar1,nchar2,nchar3;
   //Get Time information
   int index_time1=0;
   int index_time2=0;
   char hour[10],min[10],sec[10];
   for(int i1=0;i1<10;i1++){
      hour[i1]='\0';
      min[i1]='\0';
      sec[i1]='\0';
   }
   count=0;
   nchar1=0;nchar2=0;nchar3=0;
   for(int ii=size-1;ii>=0;ii--){
      if(buff[ii]==':'){
         count++;
         if(count==1) index_time2=ii+3;
         else if(count==2) index_time1=ii-4;
         continue;
      }
      if(count==1){
         min[nchar1++]=buff[ii];
      }
      else if(count==2){
         if(buff[ii]==' ') break;
         hour[nchar2++]=buff[ii];
      }
   }
   //reverse
   for(int ii=nchar1-1;ii>=0;ii--) buff0[nchar1-1-ii]=min[ii];
   for(int ii=0;ii<nchar1;ii++) min[ii]=buff0[ii];
   for(int ii=nchar2-1;ii>=0;ii--) buff0[nchar2-1-ii]=hour[ii];
   for(int ii=0;ii<nchar2;ii++) hour[ii]=buff0[ii];
   if(index_time2<=size) {sec[0]=buff[index_time2-2]; sec[1]=buff[index_time2-1];}

   int index_day1=0;
   int index_day2=0;
   char year[10],month[10],day[10];
   for(int i1=0;i1<10;i1++){
      year[i1]='\0';
      month[i1]='\0';
      day[i1]='\0';
   }
   count=0;
   nchar1=0;nchar2=0;nchar3=0;
   for(int ii=index_time1;ii>=0;ii--){
      if(buff[ii]=='-'){
         count++;
         if(count==1) index_day2=ii+3;
         else if(count==2) index_day1=((ii-5)>=0)?(ii-5):(ii-4);
         continue;
      }
      if(count==0){
         day[nchar1++]=buff[ii];
      }
      if(count==1){
         month[nchar2++]=buff[ii];
      }
      else if(count==2){
         if(nchar3>=4) break;
         year[nchar3++]=buff[ii];
      }
   }
   //reverse
   for(int ii=nchar1-1;ii>=0;ii--) buff0[nchar1-1-ii]=day[ii];
   for(int ii=0;ii<nchar1;ii++) day[ii]=buff0[ii];
   for(int ii=nchar2-1;ii>=0;ii--) buff0[nchar2-1-ii]=month[ii];
   for(int ii=0;ii<nchar2;ii++) month[ii]=buff0[ii];
   for(int ii=nchar3-1;ii>=0;ii--) buff0[nchar3-1-ii]=year[ii];
   for(int ii=0;ii<nchar3;ii++) year[ii]=buff0[ii];
   if(jdebug>1) printf("RotateDB::ProcessAll: year=%s month=%s day=%s hh=%s min=%s sec=%s\n",year,month,day,hour,min,sec);
   double time0=atoi(year)*10000000000+atoi(month)*100000000+atoi(day)*1000000+atoi(hour)*10000+atoi(min)*100+atoi(sec);
   time=CommonTools::Convert(time0);

   //start and end index for other information
   int index_start=index_time2;
   int index_end=size-1;

   //Get GPS information
   char lat[3][10],lon[3][10],alt[10];
   for(int i1=0;i1<3;i1++){
      for(int i2=0;i2<10;i2++){
         lat[i1][i2]='\0';
         lon[i1][i2]='\0';
         if(i1==0) alt[i2]='\0';
      }
   }
   if(index_time2>100){
      //Get Computer Time Information
      for(int ii=100;ii>=0;ii--){
         if(buff[ii]==':'){
            index_start=ii+3;
            index_end=index_day1;
            break;
         }
      }

      count=0;
      nchar1=0;nchar2=0;nchar3=0;
      int nchar4=0,nchar5=0,nchar6=0,nchar7=0;
      for(int ii=index_time2;ii<size;ii++){
         //printf("nline=%d char=%c\n",nline,buff[ii]);
         if(buff[ii]=='\0') break;
         if(buff[ii]=='\n') break;
         char buff2[4];
         for(int i2=0;i2<4;i2++) buff2[i2]='\0';
         for(int i2=ii;i2<(ii+((count==0||count==3)?2:3));i2++){
            if(i2<size) buff2[i2-ii]=buff[i2];
            //printf("i2=%d,ii=%d,buff2=%c\n",i2,ii,buff2[i2-ii]);
         }
         //printf("ii=%d, c=%c, %s\n",ii,buff[ii],buff2);
         if((strstr(buff2,"°")||strstr(buff2,"′")) || (strstr(buff2,"″")||strstr(buff2," "))){
            count++;
            if(strstr(buff2,"°")) ii=ii+1;
            else if(strstr(buff2,"′")||strstr(buff2,"″")) ii=ii+2;
            continue;
         }
         if(count==0){
            lon[0][nchar1++]=buff[ii];
         }
         else if(count==1){
            lon[1][nchar2++]=buff[ii];
         }
         else if(count==2){
            lon[2][nchar3++]=buff[ii];
         }
         else if(count==3){
            lat[0][nchar4++]=buff[ii];
         }
         else if(count==4){
            lat[1][nchar5++]=buff[ii];
         }
         else if(count==5){
            lat[2][nchar6++]=buff[ii];
         }
         else if(count==6){
            alt[nchar7++]=buff[ii];
         }
         else break;
      }
      if(jdebug>1) printf("RotateDB::ProcessAll: lat={%s,%s,%s} long={%s,%s,%s} alt=%s\n",lat[0],lat[1],lat[2],lon[0],lon[1],lon[2],alt);
      latitude=atoi(lat[0])+atoi(lat[1])/100.+atoi(lat[2])/10000.;
      longitude=atoi(lon[0])+atoi(lon[1])/100.+atoi(lon[2])/10000.;
      altitude=atof(alt);
   }

   //Get all kinds of swith information
   int index_start2=0;
   for(int ii=0;ii<NSWITH;ii++) allswith[ii]=false;
   int nswith=0;
   for(int ii=index_start;ii<=index_end;ii=ii+3){
      char buff3[4];
      for(int i2=0;i2<4;i2++) buff3[i2]='\0';
      for(int i2=0;i2<3;i2++){
         if(i2+ii<size) buff3[i2]=buff[i2+ii];
      }
      if(strstr(buff3,"开")){
         allswith[nswith]=true;
         index_start2=ii+3;
      }
      if(strstr(buff3,"关")){
         allswith[nswith]=false;
         index_start2=ii+3;
      }
      nswith++;
      if(nswith>=NSWITH) break;
   }
   //printf("swith=");
   //for(int ii=0;ii<19;ii++){
   //   printf("%d ",allswith[ii]);
   //}
   //printf("\n");

   //Get other variables information
   char vars[NVALUE][10];
   for(int i1=0;i1<NVALUE;i1++){
      for(int i2=0;i2<10;i2++){
         vars[i1][i2]='\0';
      }
   }
   count=0;
   int nchar[NVALUE];
   for(int ii=0;ii<NVALUE;ii++) nchar[ii]=0;
   for(int ii=index_start2;ii<=index_end;ii++){
      if(buff[ii]=='\0') break;
      if(buff[ii]=='\n') break;
      if(buff[ii]=='.'){
         if(count<5){
            vars[count][nchar[count]]=buff[ii];
            nchar[count]=nchar[count]+1;
            vars[count][nchar[count]]=buff[ii+1];
            nchar[count]=nchar[count]+1;
            ii=ii+1;
         }
         else{
            vars[count][nchar[count]]=buff[ii];
            nchar[count]=nchar[count]+1;
            vars[count][nchar[count]]=buff[ii+1];
            nchar[count]=nchar[count]+1;
            vars[count][nchar[count]]=buff[ii+2];
            nchar[count]=nchar[count]+1;
            ii=ii+2;
         }
         count++;
         if(count==NVALUE) break;
         continue;
      }
      vars[count][nchar[count]]=buff[ii];
      nchar[count]=nchar[count]+1;
   }
   //printf("vars:");
   //for(int ii=0;ii<NVALUE;ii++) printf("%s,",vars[ii]);
   //printf("\n\n");
   for(int ii=0;ii<NVALUE;ii++) varinfo[ii]=atof(vars[ii]);
   if(jdebug>1) printf("RotateDB::ProcessAll: Ele=%.2lf Ai=%.2lf\n",varinfo[9],varinfo[10]);

   //some changes due to the change of slow control software
   if(time<1573541100){
      varinfo[7]*=(-1);
      varinfo[9]*=(-1);
   }
   else{
      varinfo[10]*=(-1);
   }
   if(time<1573541855) varinfo[8]+=511.03;

   return;
}
bool RotateDB::IsLogFine(char* buffer,int RotateType){
   if(!buffer) return false;
   char tokens[5][20]={"","","","",""};
   sscanf(buffer,"%s %s %s %s %s",tokens[0],tokens[1],tokens[2],tokens[3],tokens[4]);
   /////////////////
   //RotateType=0: keyword:Plant
   //RotateType=1: keyword:Rotate
   //RotateType=2: keyword:GPS
   //RotateType=2: keyword:Environment
   /////////////////
   if(RotateType==0) return (strcmp(tokens[0],"Plant")==0);
   else if(RotateType==1) return (strcmp(tokens[0],"Rotate")==0);
   else if(RotateType==2) return (strcmp(tokens[0],"GPS")==0);
   else if(RotateType==3) return (strcmp(tokens[0],"Environment")==0);
   else return false;
}
bool RotateDB::IsLogFine(bool IsRotate){
   char tokens1[5][20]={"","","","",""};
   sscanf(buff2[0],"%s %s %s %s %s",tokens1[0],tokens1[1],tokens1[2],tokens1[3],tokens1[4]);
   char tokens2[5][20]={"","","","",""};
   sscanf(buff2[1],"%s %s %s %s %s",tokens2[0],tokens2[1],tokens2[2],tokens2[3],tokens2[4]);
   bool result;
   if(IsRotate) result=( (strcmp(tokens1[0],"Rotate")==0&&strcmp(tokens2[0],"Rotate")==0) && (strcmp(tokens1[4],"O")==0&&strcmp(tokens2[4],"O")==0) );
   else result=true;
   if(jdebug>2) printf("RotateDB::IsLogFine: token1=%s token2=%s res=%d\n",tokens1[4],tokens2[4]);
   return result;
}
int RotateDB::ReadData(ifstream* fin,bool godown,int RotateType){
   int result=0;
   if(!fin) return result;
   if(RotateType<0||RotateType>3) return result;
   if(!LocateFirst(fin)) return result;
   long int pos0=fin->tellg();
   if(pos0<0) return result;
   long int posres=pos0;
   fin->seekg(0,std::ios::beg);
   long int poss=fin->tellg();
   fin->seekg(0,std::ios::end);
   long int pose=fin->tellg();
   char buffer[200];
   //char target[20];
   /////////////////
   //RotateType=0: keyword:Plant
   //RotateType=1: keyword:Rotate
   //RotateType=2: keyword:GPS
   //RotateType=2: keyword:Environment
   /////////////////
   //char target0[4][20]={"Plant","Rotate","GPS","Environment"};
   //strcpy(target,target0[RotateType]);
   int targetlength[4]={45,64,52,56};
   long int posp,posc;
   long int linelength=targetlength[RotateType];
   if(godown){
      fin->seekg(pos0,std::ios::beg);
      posp=fin->tellg();
      posres=posp;
      if(posp<0||posp>=pose){
         fin->seekg(posres,std::ios::beg);
         if(jdebug>6) printf("RotateDB::ReadData: return0=%d\n",result);
         return result;
      }
      fin->getline(buff2[0],500);
      posc=fin->tellg();
      strcpy(buffer,buff2[0]);
      if(jdebug>6) printf("RotateDB::ReadData: item0 godown=%d pos={%ld,%ld} posse={%ld,%ld} buffer=%s\n",godown,posp,posc,poss,pose,buff2[0]);
      while(!IsLogFine(buffer,RotateType)){
         posp=fin->tellg();
         if(posp<0||posp>=pose) break;
         fin->getline(buff2[0],500);
         strcpy(buffer,buff2[0]);
         posc=fin->tellg();
         if(jdebug>7) printf("RotateDB::ReadData: item0i godown=%d pos={%ld,%ld} posse={%ld,%ld} buffer=%s\n",godown,posp,posc,poss,pose,buff2[0]);
      }
      if(IsLogFine(buffer,RotateType)){
         linelength=posc-posp;
         result=(result|0x1);
      }
      posres=posc;
      if(posc<0||posc>=pose){
         fin->seekg(posres,std::ios::beg);
         if(jdebug>6) printf("RotateDB::ReadData: return1=%d\n",result);
         return result;
      }
      fin->getline(buff2[1],500);
      strcpy(buffer,buff2[1]);
      posc=fin->tellg();
      if(jdebug>6) printf("RotateDB::ReadData: item1 godown=%d pos={%ld,%ld} posse={%ld,%ld} buffer=%s\n",godown,posp,posc,poss,pose,buff2[1]);
      while(!IsLogFine(buffer,RotateType)){
         posp=fin->tellg();
         if(posp<0||posp>=pose) break;
         fin->getline(buff2[1],500);
         strcpy(buffer,buff2[1]);
         posc=fin->tellg();
         if(jdebug>7) printf("RotateDB::ReadData: item1i godown=%d pos={%ld,%ld} posse={%ld,%ld} buffer=%s\n",godown,posp,posc,poss,pose,buff2[1]);
      }
      if(IsLogFine(buffer,RotateType)){
         result=(result|0x2);
      }
      else{
         if(result|0x1) posres=posc;
      }
      fin->seekg(posres,std::ios::beg);
      if(jdebug>6) printf("RotateDB::ReadData: return2=%d\n",result);
      return result;
   }
   else{
      fin->seekg(pos0,std::ios::beg);
      posp=fin->tellg();
      posres=posp;
      if(posp-0.5*linelength<poss){
         fin->seekg(posres,std::ios::beg);
         if(jdebug>6) printf("RotateDB::ReadData: return3=%d\n",result);
         return result;
      }
      fin->seekg(posp-0.5*linelength,std::ios::beg);
      LocateFirst(fin);
      posp=fin->tellg();
      fin->getline(buff2[1],500);
      posc=fin->tellg();
      strcpy(buffer,buff2[1]);
      if(jdebug>6) printf("RotateDB::ReadData: item0 godown=%d pos={%ld,%ld} posse={%ld,%ld} buffer=%s\n",godown,posp,posc,poss,pose,buff2[1]);
      while(!IsLogFine(buffer,RotateType)){
         if(posp<0||posp-0.5*linelength<poss) break;
         fin->seekg(posp-0.5*linelength,std::ios::beg);
         LocateFirst(fin);
         posp=fin->tellg();
         fin->getline(buff2[1],500);
         strcpy(buffer,buff2[1]);
         posc=fin->tellg();
         if(jdebug>7) printf("RotateDB::ReadData: item0i godown=%d pos={%ld,%ld} posse={%ld,%ld} buffer=%s\n",godown,posp,posc,poss,pose,buff2[1]);
      }
      if(IsLogFine(buffer,RotateType)){
         linelength=posc-posp;
         result=(result|0x1);
      }
      posres=posp;
      if(posp<0||posp-0.5*linelength<poss){
         fin->seekg(posres,std::ios::beg);
         if(jdebug>6) printf("RotateDB::ReadData: return4=%d\n",result);
         return result;
      }
      fin->seekg(posp-0.5*linelength,std::ios::beg);
      LocateFirst(fin);
      posp=fin->tellg();
      fin->getline(buff2[0],500);
      strcpy(buffer,buff2[0]);
      if(jdebug>6) printf("RotateDB::ReadData: item1 godown=%d pos={%ld,%ld} posse={%ld,%ld} buffer=%s\n",godown,posp,posc,poss,pose,buff2[0]);
      while(!IsLogFine(buffer,RotateType)){
         if(posp<0||posp-0.5*linelength<poss) break;
         fin->seekg(posp-0.5*linelength,std::ios::beg);
         LocateFirst(fin);
         posp=fin->tellg();
         fin->getline(buff2[0],500);
         strcpy(buffer,buff2[0]);
         if(jdebug>7) printf("RotateDB::ReadData: item1i godown=%d pos={%ld,%ld} posse={%ld,%ld} buffer=%s\n",godown,posp,posc,poss,pose,buff2[0]);
      }
      if(IsLogFine(buffer,RotateType)){
         result=(result|0x2);
      }
      else{
         if(result|0x1) posres=posp;
      }
      fin->seekg(posres,std::ios::beg);
      if(jdebug>6) printf("RotateDB::ReadData: return5=%d\n",result);
      return result;
   }
   return result;
}
int RotateDB::ReadData2(ifstream* fin,bool godown,bool IsRotate,int Li_in){
   if(!fin) return 0;
   char buffer0[2][500]={"",""};
   strcpy(buffer0[0],buff2[0]);
   strcpy(buffer0[1],buff2[1]);
   int RotateType=IsRotate?1:2;
   if(jdebug>5) printf("RotateDB::ReadData2: ReadData_before pos=%ld\n",fin->tellg());
   int res=ReadData(fin,godown,RotateType);
   if(jdebug>5) printf("RotateDB::ReadData2: ReadData readres=%d pos=%ld\n",res,fin->tellg());
   if(res<=0){
      strcpy(buff2[0],buffer0[0]);
      strcpy(buff2[1],buffer0[1]);
      if(jdebug>5) printf("RotateDB::ReadData2: return1=%d\n",res);
      return res;
   }
   else if((res&0x1)==0){ //should not exist
      strcpy(buff2[0],buffer0[0]);
      strcpy(buff2[1],buffer0[1]);
      if(jdebug>5) printf("RotateDB::ReadData2: return2=%d\n",0);
      return 0;
   }
   else if((res&0x2)==0){ //should read the first or last from another file
      if(godown){
         int timei=abs(ProcessTime2(0));
         int nyear=CommonTools::TimeFlag(timei+24*3600,1);
         nyear=2000+(nyear%100);
         int nmonth=CommonTools::TimeFlag(timei+24*3600,2);
         int nday=CommonTools::TimeFlag(timei+24*3600,3);
         char nfilename[150]="/scratchfs/lhaaso/hliu/rotate_log/";
         strcat(nfilename,Form("L%d/rotate/%d/%02d/log_%02d%02d.txt",Li_in,nyear,nmonth,nmonth,nday));
         ifstream fin2;
         fin2.open(nfilename,std::ios::in);
         char buffer1[500];
         strcpy(buffer1,buff2[0]);
         int res2=ReadData(&fin2,godown,RotateType);
         if(jdebug>6) printf("RotateDB::ReadData2: read the log for next day. readres=%d godown=%d filename=%s\n",res2,godown,nfilename);
         fin2.close();
         if((res2&0x1)==0){
            strcpy(buff2[0],buffer0[0]);
            strcpy(buff2[1],buffer0[1]);
            if(jdebug>5) printf("RotateDB::ReadData2: return3=%d\n",0);
            return 0;
         }
         else{
            int timei2=abs(ProcessTime2(0));
            //if(abs(timei2-timei)>3600)
            if(timei2<timei){
               strcpy(buff2[0],buffer0[0]);
               strcpy(buff2[1],buffer0[1]);
               if(jdebug>5) printf("RotateDB::ReadData2: return4=%d\n",0);
               return 0;
            }
            else{
               strcpy(buff2[1],buff2[0]);
               strcpy(buff2[0],buffer1);
               if(jdebug>5) printf("RotateDB::ReadData2: return5=%d\n",(res|0x2));
               return (res|0x2);
            }
         }
      }
      else{
         int timei=abs(ProcessTime2(1));
         int pyear=CommonTools::TimeFlag(timei-24*3600,1);
         pyear=2000+(pyear%100);
         int pmonth=CommonTools::TimeFlag(timei-24*3600,2);
         int pday=CommonTools::TimeFlag(timei-24*3600,3);
         char pfilename[150]="/scratchfs/lhaaso/hliu/rotate_log/";
         strcat(pfilename,Form("L%d/rotate/%d/%02d/log_%02d%02d.txt",Li_in,pyear,pmonth,pmonth,pday));
         ifstream fin2;
         fin2.open(pfilename,std::ios::in);
         fin2.seekg(0,std::ios::end);
         char buffer1[500];
         strcpy(buffer1,buff2[0]);
         int res2=ReadData(&fin2,godown,RotateType);
         if(jdebug>6) printf("RotateDB::ReadData2: read the log for previous day. readres=%d godown=%d filename=%s\n",res2,godown,pfilename);
         fin2.close();
         if((res2&0x2)==0){
            strcpy(buff2[0],buffer0[0]);
            strcpy(buff2[1],buffer0[1]);
            if(jdebug>5) printf("RotateDB::ReadData2: return6=%d\n",0);
            return 0;
         }
         else{
            int timei2=abs(ProcessTime2(1));
            //if(abs(timei2-timei)>3600)
            if(timei2>timei){
               strcpy(buff2[0],buffer0[0]);
               strcpy(buff2[1],buffer0[1]);
               if(jdebug>5) printf("RotateDB::ReadData2: return7=%d\n",0);
               return 0;
            }
            else{
               strcpy(buff2[0],buff2[1]);
               strcpy(buff2[1],buffer1);
               if(jdebug>5) printf("RotateDB::ReadData2: return8=%d\n",(res|0x2));
               return (res|0x2);
            }
         }
      }
      if(jdebug>5) printf("RotateDB::ReadData2: return9=%d\n",(res|0x2));
      return (res|0x2);
   }
   else{
      if(jdebug>5) printf("RotateDB::ReadData2: return10=%d\n",res);
      return res;
   }
}
long int RotateDB::LoadData2(int time_in,int Li_in,int pLi,int ptime1,int ptime2,long int cpos){
   int ptime[2]={ptime1,ptime2};
   int irot=GetLi(Li_in);
   if(irot<0) return -5;
   if(Li_in==pLi&&(ptime[0]>=0&&ptime[1]>=0)&&(time_in>=ptime[0]&&time_in<=ptime[1])){
      if(cpos>=0) return cpos;
   }
   int pirot=GetLi(pLi);
   int time_new=time_in-timedelay[irot];
   int year=CommonTools::TimeFlag(time_new,1);
   year=2000+(year%100);
   int month=CommonTools::TimeFlag(time_new,2);
   int day=CommonTools::TimeFlag(time_new,3);
   int hour=CommonTools::TimeFlag(time_new,4);
   int min=CommonTools::TimeFlag(time_new,5);
   int sec=CommonTools::TimeFlag(time_new,6);

   char filename[150]="/scratchfs/lhaaso/hliu/rotate_log/";
   char* namebuff=Form("L%d/rotate/%d/%02d/log_%02d%02d.txt",Li_in,year,month,month,day);
   strcat(filename,namebuff);
   if(jdebug>0) printf("RotateDB::LoadData2: filename=%s\n",filename);
   bool sameday=false;
   if((ptime[0]>0&&ptime[1]>0)&&pLi>0&&pirot>=0){
      int ptime_new=(ptime[0]+ptime[1])/2-timedelay[pirot];
      int pyear=CommonTools::TimeFlag(ptime_new,1);
      pyear=2000+(pyear%100);
      int pmonth=CommonTools::TimeFlag(ptime_new,2);
      int pday=CommonTools::TimeFlag(ptime_new,3);
      sameday=(Li_in==pLi)&&((year==pyear)&&(month==pmonth)&&(day==pday));
   }

   //read the file content
   char tokens[20];
   ifstream fin;
   fin.open(filename,std::ios::in);
   if(!fin.is_open()){
      //Reset();
      return -1;
   }
   fin.seekg(0,std::ios::beg);
   long int poss=fin.tellg();
   fin.seekg(0,std::ios::end);
   long int pose=fin.tellg();

   //locate the initial position
   bool located=false;
   if(sameday&&cpos>=0){
      fin.seekg(cpos,std::ios::beg);
      if(((long int)fin.tellg())<0){
         fin.close();
         fin.open(filename,std::ios::in);
      }
      located=true;
   }
   if(!located){
      fin.seekg(0,std::ios::beg);
   }

   LocateFirst(&fin);
   long int pos0=fin.tellg();
   long int posp,posc;
   if(jdebug>0) printf("RotateDB::LoadData2: initial position located. pos=%ld(%d)\n",fin.tellg(),located);
   bool godown=true;
   posp=fin.tellg();
   int firstres=ReadData2(&fin,godown,true,Li_in);
   posc=fin.tellg();
   if(firstres<=0){
      godown=false;
      fin.seekg(pos0,std::ios::beg);
      posp=fin.tellg();
      firstres=ReadData2(&fin,godown,true,Li_in);
      posc=fin.tellg();
   }
   if(firstres<=0) return -2;
   int timei0=abs(ProcessTime2(0));
   int timei1=abs(ProcessTime2(1));
   if(timei1<timei0||timei0<1000000000) return -2;
   godown=timei0<time_new;
   if(jdebug>0) printf("RotateDB::LoadData2: time_in=%d time_loop0={%d,%d} buff_loop0={%s,%s}\n",time_new,timei0,timei1,buff2[0],buff2[1]);
   bool loaded=true;
   int nloop=0;
   while(time_new<timei0||time_new>timei1){
      posp=fin.tellg();
      int fillres=ReadData2(&fin,godown,true,Li_in);
      posc=fin.tellg();
      if((fillres&0x3)<=0){
         if(jdebug>1) printf("RotateDB::LoadData2: loop=%d(%d) time=%d buff={%s,%s}\n",nloop,godown,time_new,buff2[0],buff2[1]);
         loaded=false; break;
      }
      else{
         timei0=abs(ProcessTime2(0));
         timei1=abs(ProcessTime2(1));
         if(jdebug>1) printf("RotateDB::LoadData2: loop=%d(%d) time={%d,%d,%d} buff={%s,%s}\n",nloop,godown,timei0,timei1,time_new,buff2[0],buff2[1]);
         if((time_new>=timei0&&time_new<=timei1)&&(!IsLogFine(true))) {loaded=false; break;}
      }
      nloop++;
   }

   if(jdebug>0) printf("RotateDB::LoadData2: time_in=%d time_final={%d,%d}\n",time_new,timei0,timei1);
   if(loaded&&(time_new>=timei0&&time_new<=timei1)){
      version=2;
      Li2=Li_in;
      time2[0]=timei0;
      time2[1]=timei1;
      currpos2=(posp+posc)/2;
      fin.close();
      if(jdebug>0) printf("RotateDB::LoadData2: loaded, Li=%d time_in=%d time={%d,%d} cpos=%ld posse={%ld,%ld} buff={%s,%s}\n",time_new,Li2,time2[0],time2[1],currpos2,poss,pose,buff2[0],buff2[1]);
      return currpos2;
   }
   else{
      //Reset();
      fin.close();
      return -4;
   }
}
int RotateDB::ProcessTime2(int itime){
   if(itime<0||itime>1) return 0;

   char tokens[14][20];
   char date0[20];
   char time0[20];
   for(int ii=0;ii<14;ii++){
      for(int jj=0;jj<20;jj++){
         tokens[ii][jj]='\0';
         if(ii==0){
            date0[jj]='\0';
            time0[jj]='\0';
         }
      }
   }
   sscanf(buff2[itime], "%s %s %s %s %s %s %s %s %s %s %s %s %s %s",tokens[0],tokens[1],tokens[2],tokens[3],tokens[4],tokens[5],tokens[6],tokens[7],tokens[8],tokens[9],tokens[10],tokens[11],tokens[12],tokens[13]);

   //printf("tokens=%s\n",tokens[0]);
   if(strcmp(tokens[0],"Rotate")==0){
      strcpy(date0,tokens[1]);
      strcpy(time0,tokens[2]);
   }
   else if(strcmp(tokens[0],"GPS")==0){
      strcpy(date0,tokens[2]);
      strcpy(time0,tokens[3]);
   }
   else return 0;

   int size=strlen(date0);
   int count;
   int nchar1,nchar2,nchar3;
   //Get Date information
   char year[10],month[10],day[10];
   for(int i1=0;i1<10;i1++){
      year[i1]='\0';
      month[i1]='\0';
      day[i1]='\0';
   }
   count=0;
   nchar1=0;nchar2=0;nchar3=0;
   for(int ii=0;ii<size;ii++){
      if(date0[ii]=='-'){
         count++;
         continue;
      }
      if(count==0){
         year[nchar1++]=date0[ii];
      }
      else if(count==1){
         month[nchar2++]=date0[ii];
      }
      else if(count==2){
         if(date0[ii]=='\0'||date0[ii]=='\n') break;
         day[nchar3++]=date0[ii];
      }
   }
   //printf("date0=%s\n",date0);
   //printf("date=%s %s %s\n",year,month,day);
   char hour[10],min[10],sec[10];
   for(int i1=0;i1<10;i1++){
      hour[i1]='\0';
      min[i1]='\0';
      sec[i1]='\0';
   }
   size=strlen(time0);
   count=0;
   nchar1=0;nchar2=0;nchar3=0;
   for(int ii=1;ii<size;ii++){
      if(time0[ii]==':'){
         count++;
         continue;
      }
      if(count==0){
         hour[nchar1++]=time0[ii];
      }
      else if(count==1){
         min[nchar2++]=time0[ii];
      }
      else if(count==2){
         if(time0[ii]=='\0'||time0[ii]=='\n') break;
         sec[nchar3++]=time0[ii];
      }
   }
   //printf("time0=%s\n",time0);
   //printf("time=%s %s %s\n",hour,min,sec);

   double time00=atoi(year)*10000000000+atoi(month)*100000000+atoi(day)*1000000+atoi(hour)*10000+atoi(min)*100+atoi(sec);
   int result=CommonTools::Convert(time00);
   return (time0[0]=='+')?result:(-result);
}
void RotateDB::ProcessAll2(){
   if(time2<=0) return;
   int size=strlen(buff2[0]);
   if(jdebug>0) printf("RotateDB::ProcessAll2: %s(size=%d)\n",buff2[0],size);
   char tokens[14][20];
   sscanf(buff2[0],"%s %s %s %s %s %s %s %s %s %s %s %s %s %s",tokens[0],tokens[1],tokens[2],tokens[3],tokens[4],tokens[5],tokens[6],tokens[7],tokens[8],tokens[9],tokens[10],tokens[11],tokens[12],tokens[13]);
   int count;
   int nchar1,nchar2,nchar3;
   if(strcmp(tokens[0],"Rotate")==0){
      //get time information
      char date0[3][10];
      size=strlen(tokens[1]);
      count=0;
      nchar1=0,nchar2=0,nchar3=0;
      for(int ii=0;ii<3;ii++){
         for(int jj=0;jj<10;jj++){
            date0[ii][jj]='\0';
         }
      }
      for(int ii=0;ii<size;ii++){
         if(tokens[1][ii]=='-'){
            count++;
            continue;
         }
         if(count==0){
            date0[0][nchar1++]=tokens[1][ii];
         }
         else if(count==1){
            date0[1][nchar2++]=tokens[1][ii];
         }
         else if(count==2){
            if(tokens[1][ii]=='\0'||tokens[1][ii]=='\n') break;
            date0[2][nchar3++]=tokens[1][ii];
         }
      }
      double time0=atoi(date0[0])*10000000000+atoi(date0[1])*100000000+atoi(date0[2])*1000000;
      size=strlen(tokens[2]);
      count=0;
      nchar1=0,nchar2=0,nchar3=0;
      for(int ii=0;ii<3;ii++){
         for(int jj=0;jj<10;jj++){
            date0[ii][jj]='\0';
         }
      }
      for(int ii=1;ii<size;ii++){
         if(tokens[2][ii]=='-'){
            count++;
            continue;
         }
         if(count==0){
            date0[0][nchar1++]=tokens[2][ii];
         }
         else if(count==1){
            date0[1][nchar2++]=tokens[2][ii];
         }
         else if(count==2){
            if(tokens[2][ii]=='\0'||tokens[2][ii]=='\n') break;
            date0[2][nchar3++]=tokens[2][ii];
         }
      }
      time0+=atoi(date0[0])*10000+atoi(date0[1])*100+atoi(date0[2]);
      time=CommonTools::Convert(time0);
      //get swith information
      allswith[4]=(strcmp(tokens[3],"O")==0); //rotate
      allswith[0]=(strcmp(tokens[4],"O")==0); //laser
      allswith[1]=(strcmp(tokens[5],"O")==0); //electric heat
      allswith[2]=(strcmp(tokens[6],"O")==0); //electric air heat
      allswith[3]=(strcmp(tokens[7],"O")==0); //power
      allswith[5]=(strcmp(tokens[8],"O")==0); //door
      for(int ii=6;ii<NSWITH;ii++) allswith[ii]=true;
      //get value information
      for(int ii=0;ii<NVALUE;ii++) varinfo[ii]=0;
      varinfo[8]=atof(tokens[9]); //height
      varinfo[11]=atof(tokens[10]); //angle
      varinfo[9]=atof(tokens[11]); //elevation
      varinfo[10]=atof(tokens[12]); //azimuth
      //some changes due to the change of slow control software
      if(time<1573541100){
         varinfo[7]*=(-1);
         varinfo[9]*=(-1);
      }
      else{
         varinfo[10]*=(-1);
      }
      if(time<1573541855) varinfo[8]+=511.03;
   }
   else if(strcmp(tokens[0],"GPS")==0){
      size=strlen(tokens[4]);
      char coo3[3][10];
      count=0;
      nchar1=0,nchar2=0,nchar3=0;
      for(int ii=0;ii<3;ii++){
         for(int jj=0;jj<10;jj++){
            coo3[ii][jj]='\0';
         }
      }
      for(int ii=0;ii<size;ii++){
         if(tokens[4][ii]=='-'){
            count++;
            continue;
         }
         if(count==0){
            coo3[0][nchar1++]=tokens[4][ii];
         }
         else if(count==1){
            coo3[1][nchar2++]=tokens[4][ii];
         }
         else if(count==2){
            if(tokens[4][ii]=='\0'||tokens[4][ii]=='\n') break;
            coo3[2][nchar3++]=tokens[4][ii];
         }
      }
      latitude=atoi(coo3[0])+atoi(coo3[1])/100.+atoi(coo3[2])/10000.;
      size=strlen(tokens[5]);
      count=0;
      nchar1=0,nchar2=0,nchar3=0;
      for(int ii=0;ii<3;ii++){
         for(int jj=0;jj<10;jj++){
            coo3[ii][jj]='\0';
         }
      }
      for(int ii=0;ii<size;ii++){
         if(tokens[5][ii]=='-'){
            count++;
            continue;
         }
         if(count==0){
            coo3[0][nchar1++]=tokens[5][ii];
         }
         else if(count==1){
            coo3[1][nchar2++]=tokens[5][ii];
         }
         else if(count==2){
            if(tokens[5][ii]=='\0'||tokens[5][ii]=='\n') break;
            coo3[2][nchar3++]=tokens[5][ii];
         }
      }
      longitude=atoi(coo3[0])+atoi(coo3[1])/100.+atoi(coo3[2])/10000.;
      altitude=atof(tokens[6]);
   }
   return;
}
int RotateDB::ProcessEnv(int time_in,int Li_in){
   int irot=GetLi(Li_in);
   if(irot<0) return -4;
   int time_new=time_in-timedelay[irot];
   int year=CommonTools::TimeFlag(time_new,1);
   year=2000+(year%100);
   int month=CommonTools::TimeFlag(time_new,2);
   int day=CommonTools::TimeFlag(time_new,3);
   char filename[100]="/scratchfs/lhaaso/hliu/rotate_log/";
   char* namebuff=Form("L%d/rotate/%04d/%02d/var_%02d%02d.txt",Li_in,year,month,month,day);
   strcat(filename,namebuff);
   if(jdebug>0) printf("RotateDB::ProcessEnv: filename=%s\n",filename);
   ifstream fin;
   fin.open(filename,std::ios::in);
   bool isopen=fin.is_open();
   if(!isopen) return -3;
   bool located=false;
   double vars12[2][6];
   for(int ii=0;ii<6;ii++){
      vars12[0][ii]=-100;
      vars12[1][ii]=-100;
   }
   int timemargin=5*60;
   while(!fin.eof()){
      char linebuff[300];
      fin.getline(linebuff,300);
      int size=strlen(linebuff);
      if(size<10) break;
      char tokens[8][20];
      sscanf(linebuff,"%s %s %s %s %s %s %s %s",tokens[0],tokens[1],tokens[2],tokens[3],tokens[4],tokens[5],tokens[6],tokens[7],tokens[8]);
      if(strcmp(tokens[0],"Environment")==0){
         //time
         size=strlen(tokens[2]);
         int count,nchar1,nchar2,nchar3;
         char hour[10];
         char min[10];
         char sec[10];
         for(int ii=0;ii<10;ii++){
            hour[ii]='\0';
            min[ii]='\0';
            sec[ii]='\0';
         }
         count=0;
         nchar1=0;nchar2=0;nchar3=0;
         for(int ii=1;ii<size;ii++){
            if(tokens[2][ii]==':'){
               count++;
               continue;
            }
            if(count==0){
               hour[nchar1++]=tokens[2][ii];
            }
            else if(count==1){
               min[nchar2++]=tokens[2][ii];
            }
            else if(count==2){
               if(tokens[2][ii]=='\0'||tokens[2][ii]=='\n') break;
               sec[nchar3++]=tokens[2][ii];
            }
         }
         double time00=year*10000000000.+month*100000000.+day*1000000+atoi(hour)*10000+atoi(min)*100+atoi(sec);
         int timei=CommonTools::Convert(time00);
         if(!located){
            if(timei>=time_new-timemargin&&timei<time_new+timemargin){
               version=3;
               located=true;
               vars12[0][0]=timei*1.;
               for(int ii=0;ii<5;ii++) vars12[0][ii+1]=atof(tokens[ii+3]);
            }
         }
         else{
            if(fabs(timei-vars12[0][0])>15*60) break;
            vars12[1][0]=timei*1.;
            for(int ii=0;ii<5;ii++) vars12[1][ii+1]=atof(tokens[ii+3]);
            break;
         }
      }
   }
   if(located){
      if(vars12[1][0]<0){
         for(int ii=0;ii<5;ii++) varinfo[ii]=vars12[0][ii+1];
      }
      else{
         for(int ii=0;ii<5;ii++){
            double ratio=(time_new-vars12[0][0])/(vars12[1][0]-vars12[0][0]);
            varinfo[ii]=vars12[0][ii+1]*(1-ratio)+vars12[1][ii+1]*ratio;
         }
      }
   }
   return located?1:-1;
}
void RotateDB::DumpInfo(){
   printf("RotateDB::DumpInfo: buff=%s\n",buff);
   int year=CommonTools::TimeFlag(time,1);
   int month=CommonTools::TimeFlag(time,2);
   int day=CommonTools::TimeFlag(time,3);
   int hour=CommonTools::TimeFlag(time,4);
   int min=CommonTools::TimeFlag(time,5);
   int sec=CommonTools::TimeFlag(time,6);
   printf("time=%d %4d-%02d-%02d %02d:%02d:%02d\n",time,year,month,day,hour,min,sec);
   if(time<=0) return;
   printf("latitude=%.4lf longitude=%.4lf altitude=%.2lf\n",latitude,longitude,altitude);
   printf("swith:");
   for(int ii=0;ii<NSWITH;ii++) printf("%d ",allswith[ii]);
   printf("\nvars:");
   for(int ii=0;ii<NVALUE;ii++) printf("%.2lf,",varinfo[ii]);
   printf("\n\n");
}
bool RotateDB::GetLaserSwith(){
   return allswith[0];
}
bool RotateDB::GetDoorSwith(){
   return allswith[6];
}
double RotateDB::GetTemperature(int itemp){
   if(itemp<0||itemp>3) return -100;
   return varinfo[itemp];
}
double RotateDB::GetHumidity(){
   return varinfo[4];
}
double RotateDB::GetInclination(bool IsX){
   return varinfo[IsX?5:6];
}
double RotateDB::GetDoorAngle(){
   return varinfo[7];
}
double RotateDB::GetHeight(){
   return varinfo[8];
}
void RotateDB::GetAngles(double angle[3]){
   angle[0]=varinfo[9];
   angle[1]=varinfo[10];
   angle[2]=varinfo[11];
}
double RotateDB::GetElevation(){
   double angle[3];
   GetAngles(angle);
   return angle[0];
}
double RotateDB::GetAzimuth(){
   double angle[3];
   GetAngles(angle);
   return angle[1];
}
double RotateDB::GetAng(){
   double angle[3];
   GetAngles(angle);
   return angle[2];
}

int RotateDB::GetLi(){
   int Li_in=Li;
   return GetLi(Li_in);
}
int RotateDB::GetLi(int Li_in){
   int result=-1;
   for(int irot=0;irot<nrot;irot++){
      if(Li_in==rotindex[irot]) {result=irot; break;}
   }
   return result;
}
int RotateDB::GetLi(double rabbittime){
   int result=-1;
   for(int irot=0;irot<nrot;irot++){
      if(fabs(rabbittime*20-rottime[irot])<1600*20) {result=irot; break;}
   }
   return result;
}
int RotateDB::GetTi(int Tindex){
   int result=-1;
   for(int itel=0;itel<ntel;itel++){
      if(Tindex==telindex[itel]) {result=itel; break;}
   }
   return result;
}
double RotateDB::GetMinDistEleAzi(double ele_in,double azi_in,int irot,int itel,double &minele,double &minazi,int &index){
   index=-1;
   minele=1000;
   minazi=1000;
   double result=TMath::Max(minele,minazi);
   if(irot<0||irot>=2) return result;
   if(itel<0||itel>=6) return result;

   if(rotindex[irot]==2){
      const int nangle1=4;
      double AngleList1[nangle1][2];
      AngleList1[0][0]=40; AngleList1[0][1]=42;
      AngleList1[1][0]=40; AngleList1[1][1]=29;
      AngleList1[2][0]=40; AngleList1[2][1]=19;
      AngleList1[3][0]=40; AngleList1[3][1]=4;
      for(int ii=0;ii<nangle1;ii++){
         bool isfine=false;
         if(ii==0&&(telindex[itel]==5||telindex[itel]==6)) isfine=true;
         if(ii==1&&(telindex[itel]==4||telindex[itel]==5)) isfine=true;
         if(ii==2&&(telindex[itel]==3||telindex[itel]==4)) isfine=true;
         if(ii==3&&(telindex[itel]==1||telindex[itel]==2||telindex[itel]==3)) isfine=true;
         if(!isfine) continue;
         double dist1=fabs(ele_in-AngleList1[ii][0]);
         double dist2=fabs(azi_in-AngleList1[ii][1]);
         double dist=TMath::Max(dist1,dist2);
         if(dist<result){
            minele=dist1;
            minazi=dist2;
            result=dist;
            index=1*10+ii;
         }
      }
   }
   const int nangle2=7;
   double elelist[nangle2]={10,20,30,40,50,55,60};
   double azilist[2][6][nangle2]={{{23.7,19,15,10,1,-5,-1000},
                                   {23.7,17,11,0,-12,-18,-1000},
                                   {24.5,22,19,15,10,3,-1000},
                                   {27,26,26,26,25,24,-1000},
                                   {30.5,33,35,38,42,44,-1000},
                                   {32,38,43,47,54,61,-1000}},
                                  {{-71,-67,-63,-57,-52,-42,-37},
                                   {-74,-72,-72,-69,-67,-66,-65},
                                   {-69,-64,-57,-47,-37,-27,-17},
                                   {-67,-62,-54,-47,-32,-17,-2},
                                   {-68,-63,-57,-48,-32,-13,3},
                                   {-71,-67,-62,-57,-52,-42,-37}}
                                 };
   for(int ii=0;ii<nangle2;ii++){
      double dist1=fabs(ele_in-elelist[ii]);
      double dist2=fabs(azi_in-azilist[irot][itel][ii]);
      double dist=TMath::Max(dist1,dist2);
      if(dist<result){
         minele=dist1;
         minazi=dist2;
         result=dist;
         index=2*10+ii;
      }
   }
   return result;
}
int RotateDB::IsFineAngle(double ele_in,double azi_in,int Li_in,int iTel){
   int irot=GetLi(Li_in);
   int itel=GetTi(iTel);
   int index;
   double mindistele,mindistazi;
   double mindist=GetMinDistEleAzi(ele_in,azi_in,irot,itel,mindistele,mindistazi,index);
   if(mindist<aglmargin) return index;
   else return -1;
}
int RotateDB::GetEleAzi1(int time_in,int Li_in,int iTel){
   if(jdebug>0) printf("RotateDB::GetEleAzi1: timein=%d Lin=%d iTel=%d\n",time_in,Li_in,iTel);
   if(LoadData(time_in,Li_in,Li,time,currpos)<0) return -1;
   ProcessAll();
   RotateDB rotbuff(this);

   double ele0=GetElevation();
   double azi0=GetAzimuth();
   if(jdebug>0) printf("RotateDB::GetEleAzi1: findtimelog ele=%lf azi=%lf\n",ele0,azi0);
   int retval=-1;
   if(iTel>0){
      retval=IsFineAngle(ele0,azi0,Li_in,iTel);
      if(retval<=0){
         if(jdebug>0) printf("RotateDB::GetEleAzi1: IsFineAngle=%d\n",retval);
         return -2;
      }
   }
   else{
      printf("RotateDB::GetEleAzi1: wrong iTel=%d\n",iTel);
      return -3;
   }

   int pLi=rotbuff.Li;
   int ctime=rotbuff.time;
   long int cpos=rotbuff.currpos;

   int ncount1=0,ncount2=0;
   for(int itime=-1;itime>=-(60*100);itime--){
      int timei=time_in+itime;
      if(LoadData(timei,Li_in,pLi,ctime,cpos)<0) break;
      ProcessAll();
      pLi=Li;
      ctime=time;
      cpos=currpos;
      double elei=GetElevation();
      double azii=GetAzimuth();
      if(fabs(ele0-elei)>aglmargin) break;
      if(fabs(azi0-azii)>aglmargin) break;
      ncount1++;
   }
   pLi=rotbuff.Li;
   ctime=rotbuff.time;
   cpos=rotbuff.currpos;
   for(int itime=1;itime<=(60*100);itime++){
      int timei=time_in+itime;
      if(LoadData(timei,Li_in,pLi,ctime,cpos)<0) break;
      ProcessAll();
      pLi=Li;
      ctime=time;
      cpos=currpos;
      double elei=GetElevation();
      double azii=GetAzimuth();
      if(fabs(ele0-elei)>aglmargin) break;
      if(fabs(azi0-azii)>aglmargin) break;
      ncount2++;
   }

   //LoadData(time_in,Li_in);
   //ProcessAll();
   this->Copy(&rotbuff);

   //double type0[4][2]={{40,42},{40,29},{40,19},{40,4}};
   //bool Istype0=false;
   //for(int ii=0;ii<4;ii++){
   //   if(fabs(ele0-type0[ii][0])<aglmargin&&fabs(azi0-type0[ii][1])<aglmargin) Istype0=true;
   //}

   int ncount=ncount1+ncount2+1;
   //int nside1=((ncount%2)==0)?((ncount/2)-13):(((ncount-1)/2)-12);
   //int nside2=((ncount%2)==0)?((ncount/2)-13):(((ncount-1)/2)-13);
   //if(jdebug>0) printf("RotateDB::GetEleAzi: ncount={%d,%d} nside={%d,%d}\n",ncount1,ncount2,nside1,nside2);

   if(ncount<ntotmin) return -4;
   if(ncount1<=nsidemin||ncount2<=nsidemin) return -5;
   else return (Li_in*100+retval);
}
int RotateDB::GetEleAzi2(int time_in,int Li_in,int iTel){
   if(jdebug>0) printf("RotateDB::GetEleAzi2: timein=%d Lin=%d iTel=%d\n",time_in,Li_in,iTel);
   int irot=GetLi(Li_in);
   if(irot<0) return -1;
   if(LoadData2(time_in,Li_in,Li2,time2[0],time2[1],currpos2)<0) return -1;
   ProcessAll2();

   double ele0=GetElevation();
   double azi0=GetAzimuth();
   if(jdebug>0) printf("RotateDB::GetEleAzi2: findtimelog ele=%lf azi=%lf\n",ele0,azi0);
   int retval=-1;
   if(iTel>0){
      retval=IsFineAngle(ele0,azi0,Li_in,iTel);
      if(retval<=0){
         if(jdebug>0) printf("RotateDB::GetEleAzi2: IsFineAngle=%d\n",retval);
         return -2;
      }
   }
   else{
      printf("RotateDB::GetEleAzi2: wrong iTel=%d\n",iTel);
      return -3;
   }

   int timei0=abs(ProcessTime2(0));
   int timei1=abs(ProcessTime2(1));
   int ncount1=(time_in-timedelay[irot])-timei0;
   int ncount2=timei1-(time_in-timedelay[irot]);
   int ncount=ncount1+ncount2+1;

   if(ncount<ntotmin) return -4;
   if(ncount1<=nsidemin||ncount2<=nsidemin) return -5;
   else return (Li_in*100+retval);
}
int RotateDB::GetEleAzi(int time_in,int Li_in,int iTel){
   int irot=GetLi(Li_in);
   if(irot<0) return -10;
   int time_new=time_in-timedelay[irot];
   int year=CommonTools::TimeFlag(time_new,1);
   year=2000+(year%100);
   int month=CommonTools::TimeFlag(time_new,2);
   int day=CommonTools::TimeFlag(time_new,3);

   char name1[300]="";
   char name2[300]="";
   char filename[100]="/scratchfs/lhaaso/hliu/rotate_log/";
   char* namebuff1=Form("L%d/%d-%02d-%02d.txt.utf8",Li_in,year,month,day);
   char* namebuff2=Form("L%d/rotate/%d/%02d/log_%02d%02d.txt",Li_in,year,month,month,day);
   strcat(name1,filename);
   strcat(name1,namebuff1);
   strcat(name2,filename);
   strcat(name2,namebuff2);

   ifstream fin1;
   fin1.open(name1,std::ios::in);
   bool file1=fin1.is_open();
   ifstream fin2;
   fin2.open(name2,std::ios::in);
   bool file2=fin2.is_open();

   if(file1){
      if(!file2) return GetEleAzi1(time_in,Li_in,iTel);
      else{
         int res1=GetEleAzi1(time_in,Li_in,iTel);
         int res2=GetEleAzi2(time_in,Li_in,iTel);
         if(res1==res2) return res1;
         else if(res1<0) return res2;
         else if(res2<0) return res1;
         else return -11;
      }
   }
   else{
      if(file2) return GetEleAzi2(time_in,Li_in,iTel);
      else return -12;
   }
}
int RotateDB::GetEleAzi(WFCTAEvent* pev){
   if(!pev) return -6;
   int time_in=pev->rabbitTime;
   int irot=GetLi((double)pev->rabbittime);
   if(irot<0) return -7;
   int Li_in=rotindex[irot];
   return GetEleAzi(time_in,Li_in,pev->iTel);
}
void RotateDB::GetMinDistFit(WFCTAEvent* pev,int EleAziIndex,int Li_in,double &minphi,double &mincc){
   minphi=1000;
   mincc=1000;
   if(!pev) return;
   if(!pev->minimizer) return;
   double phi=pev->minimizer->X()[3]/PI*180;
   double cc=pev->minimizer->X()[2]/PI*180;

   int irot=Li_in>0?GetLi(Li_in):GetLi((double)pev->rabbittime);
   if(irot<0||irot>=2) return;
   int itel=GetTi(pev->iTel);
   if(itel<0||itel>=6) return;
   
   if(EleAziIndex<=0) return;
   int itype=(EleAziIndex%100)/10;
   int index=EleAziIndex%10;

   if(itype==1&&rotindex[irot]==2){
      const int nangle1=4;
      if(index<0||index>=nangle1) return;
      double phi1[6][nangle1];
      double cc1[6][nangle1];
      for(int ii=0;ii<nangle1;ii++){
         for(int jj=0;jj<6;jj++){
            phi1[jj][ii]=-10000;
            cc1[jj][ii]=-10000;
            if(ii==0&&telindex[jj]==6){
               phi1[jj][ii]=134.93;
               cc1[jj][ii]=-5.55;
            }
            if(ii==0&&telindex[jj]==5){
               phi1[jj][ii]=110.89;
               cc1[jj][ii]=6.30;
            }
            if(ii==1&&telindex[jj]==5){
               phi1[jj][ii]=100.29;
               cc1[jj][ii]=-8.20;
            }
            if(ii==1&&telindex[jj]==4){
               phi1[jj][ii]=82.93;
               cc1[jj][ii]=5.92;
            }
            if(ii==2&&telindex[jj]==4){
               phi1[jj][ii]=76.82;
               cc1[jj][ii]=-5.72;
            }
            if(ii==2&&telindex[jj]==3){
               phi1[jj][ii]=54.42;
               cc1[jj][ii]=6.61;
            }
            if(ii==3&&telindex[jj]==3){
               phi1[jj][ii]=49.25;
               cc1[jj][ii]=-8.80;
            }
            if(ii==3&&telindex[jj]==2){
               phi1[jj][ii]=3.00;
               cc1[jj][ii]=3.03;
            }
            if(ii==3&&telindex[jj]==1){
               phi1[jj][ii]=26.04;
               cc1[jj][ii]=-1.72;
            }
         }
      }

      double dist1=TMath::Min(fabs(phi-phi1[itel][index]),fabs(phi+180-phi1[itel][index]));
      dist1=TMath::Min(dist1,fabs(phi-180-phi1[itel][index]));
      double dist2=fabs(cc-cc1[itel][index]);
      if(dist1<minphi) minphi=dist1;
      if(dist2<mincc) mincc=dist2;
   }
   else if(itype==2){
      const int nangle2=7;
      if(index<0||index>=nangle2) return;
      double phi2[2][6][nangle2]={{{27.97,27.40,27.08,26.44,26.77,26.64,-10000},
                                  {3.00,3.17,3.16,3.13,3.16,3.07,-10000},
                                  {49.29,50.43,51.88,51.91,51.83,50.56,-10000},
                                  {75.00,78.46,79.76,80.80,80.78,80.40,-10000},
                                  {111.32,111.31,110.69,110.68,110.56,110.94,-10000},
                                  {135.70,137.83,137.27,136.46,136.46,136.62,-10000}
                                 },
                                 {{134.08,134.34,133.91,134.05,133.69,133.19,-10000},
                                  {106.63,105.91,104.02,105.04,104.61,-1000,-10000},
                                  {156.15,157.01,157.54,158.09,157.82,157.85,-10000},
                                  {0.00,0.09,0.15,179.82,0.10,179.63,-10000},
                                  {24.68,24.19,24.13,23.95,23.50,22.61,-10000},
                                  {48.75,48.35,48.50,48.62,49.40,48.29,-10000}
                                 }
                                };
      double cc2[2][6][nangle2]={{{-1.99,0.11,2.81,4.13,3.34,3.15,-10000},
                                  {3.88,0.90,2.01,0.70,-0.20,1.23,-10000},
                                  {-2.45,0.94,2.00,2.04,2.37,0.17,-10000},
                                  {-1.11,-0.65,1.30,2.34,2.12,1.82,-10000},
                                  {0.95,2.31,1.66,1.99,2.23,1.86,-10000},
                                  {-2.52,2.98,2.37,-0.09,-0.54,0.59,-10000}
                                 },
                                 {{-0.75,-0.14,-0.84,-0.46,-2.24,0.56,-10000},
                                  {-1.10,0.98,-1.79,-0.01,-0.37,-1000,-10000},
                                  {-0.99,-2.01,-0.69,0.77,-0.81,-0.50,-10000},
                                  {-2.71,1.45,0.19,-2.36,1.25,-0.03,-10000},
                                  {-1.01,1.37,1.88,0.99,-0.90,-3.19,-10000},
                                  {1.16,-0.21,-1.11,-0.39,1.36,-1.30,-10000}
                                 }
                                };
      double dist1=TMath::Min(fabs(phi-phi2[irot][itel][index]),fabs(phi+180-phi2[irot][itel][index]));
      dist1=TMath::Min(dist1,fabs(phi-180-phi2[irot][itel][index]));
      double dist2=fabs(cc-cc2[irot][itel][index]);
      if(dist1<minphi) minphi=dist1;
      if(dist2<mincc) mincc=dist2;
   }
}
bool RotateDB::IsFineImage(WFCTAEvent* pev,int EleAziIndex,int Li_in){
   double minphi,mincc;
   GetMinDistFit(pev,EleAziIndex,Li_in,minphi,mincc);
   return (minphi<phimargin&&mincc<ccmargin);
}
double RotateDB::GetEref(int Li_in,int iTel,int type,int iangle){
   int irot=GetLi(Li_in);
   if(irot<0||irot>=nrot) return -1;
   int itel=GetTi(iTel);
   if(itel<0||itel>=ntel) return -1;
   if(type==1){
      const int nangle1=4;
      double npe[6][nangle1];
      for(int ii=0;ii<6;ii++){
         for(int jj=0;jj<nangle1;jj++){
            npe[ii][jj]=1.31e6;
         }
      }
      if(iangle<0||iangle>=nangle1) return -1;
      else return npe[itel][iangle];
   }
   else if(type==2){
      const int nangle2=6;
      double npe2[6][nangle2];
      for(int ii=0;ii<6;ii++){
         for(int jj=0;jj<nangle2;jj++){
            npe2[ii][jj]=1.31e6;
         }
      }
      npe2[0][0]=2.34e6;
      npe2[0][1]=1.31e6;
      npe2[0][2]=9.53e5;
      npe2[0][3]=7.74e5;
      npe2[0][4]=6.89e5;
      npe2[0][5]=6.1e5;
      if(iangle<0||iangle>=nangle2) return -1;
      else return npe2[itel][iangle];
   }
   else return -1;
}
int RotateDB::LaserIsFine(WFCTAEvent* pev){
   int EleAziIndex=GetEleAzi(pev);
   if(jdebug>0) printf("RotateDB::LaserIsFine: EleAziIndex=%d\n",EleAziIndex);
   if(EleAziIndex<=0) return EleAziIndex;

   bool fitted=pev->DoFit(0,3);
   if(jdebug>1) printf("RotateDB::LaserIsFine: Fit=%d\n",fitted);
   bool image=IsFineImage(pev,EleAziIndex);
   if(jdebug>1) printf("RotateDB::LaserIsFine: FineImage=%d\n",image);

   return image?EleAziIndex:-10;
}
bool RotateDB::GetEnv(int time_in,int Li_in,double *temp){
   if(!temp) return false;
   int irot=GetLi(Li_in);
   if(irot<0) return false;
   int time_new=time_in;
   int year=CommonTools::TimeFlag(time_new,1);
   year=2000+(year%100);
   int month=CommonTools::TimeFlag(time_new,2);
   int day=CommonTools::TimeFlag(time_new,3);

   char name1[300]="";
   char name2[300]="";
   char filename[100]="/scratchfs/lhaaso/hliu/rotate_log/";
   char* namebuff1=Form("L%d/%d-%02d-%02d.txt.utf8",Li_in,year,month,day);
   char* namebuff2=Form("L%d/rotate/%d/%02d/var_%02d%02d.txt",Li_in,year,month,month,day);
   strcat(name1,filename);
   strcat(name1,namebuff1);
   strcat(name2,filename);
   strcat(name2,namebuff2);

   ifstream fin1;
   fin1.open(name1,std::ios::in);
   bool file1=fin1.is_open();
   if(file1) fin1.close();
   ifstream fin2;
   fin2.open(name2,std::ios::in);
   bool file2=fin2.is_open();
   if(file2) fin2.close();

   if(file1){
      if(LoadData(time_in,Li_in,Li,time,currpos)<0){
         if(!file2) return false;
         else{
            if(ProcessEnv(time_in,Li_in)<=0) return false;
            else{
               for(int ii=0;ii<4;ii++) temp[ii]=GetTemperature(ii);
               temp[4]=GetHumidity();
               return true;
            }
         }
      }
      ProcessAll();
      for(int ii=0;ii<4;ii++) temp[ii]=GetTemperature(ii);
      temp[4]=GetHumidity();
      return true;
   }
   else{
      if(file2){
         if(ProcessEnv(time_in,Li_in)<=0) return false;
         else{
            for(int ii=0;ii<4;ii++) temp[ii]=GetTemperature(ii);
            temp[4]=GetHumidity();
            return true;
         }
      }
      else return false;
   }
}
