#include "RotateDB.h"
#include "stdlib.h"
#include <fstream>
#include <string>
RotateDB* RotateDB::_Head=0;
int RotateDB::jdebug=0;
int RotateDB::timedelay=35;
int RotateDB::ntotmin=10;
int RotateDB::nsidemin=3;
RotateDB::RotateDB(RotateDB* pr_in){
   Copy(pr_in);
}
void RotateDB::Init(){
   Li=0;
   time=0;
   currpos=-1;
   Reset();
}
void RotateDB::Reset(){
   for(int ii=0;ii<500;ii++){
      buff[ii]='\0';
   }
   latitude=0;
   longitude=0;
   altitude=0;
   for(int ii=0;ii<NSWITH;ii++) allswith[ii]=false;
   for(int ii=0;ii<NVALUE;ii++) varinfo[ii]=0;
}
void RotateDB::Release(){
   Init();
}
void RotateDB::Copy(RotateDB* pr_in){
   for(int ii=0;ii<500;ii++){
      buff[ii]=pr_in->buff[ii];
   }
   Li=pr_in->Li;
   time=pr_in->time;
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
   int time_new=time_in-timedelay;
   int year=CommonTools::TimeFlag(time_new,1);
   year=2000+(year%100);
   int month=CommonTools::TimeFlag(time_new,2);
   int day=CommonTools::TimeFlag(time_new,3);
   int hour=CommonTools::TimeFlag(time_new,4);
   int min=CommonTools::TimeFlag(time_new,5);
   int sec=CommonTools::TimeFlag(time_new,6);
   char filename[100]="/eos/user/h/hliu/rotate_log/";
   char* namebuff=Form("L%d/%d-%02d-%02d.txt.utf8",Li_in,year,month,day);
   strcat(filename,namebuff);
   if(jdebug>0) printf("RotateDB::LoadData: filename=%s\n",filename);
   bool sameday=false;
   if(ptime>0&&pLi>0){
      int ptime_new=ptime-timedelay;
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
      Reset();
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
      if(hour<12){
         long int posi=(hour*3600+min*60+sec)*linelength;
         if(poss+posi>=pose){
            fin.seekg(-10,std::ios::end);
         }
         else{
            fin.seekg((hour*3600+min*60+sec)*linelength,std::ios::beg);
         }
      }
      else{
         fin.seekg(-10,std::ios::end);
         if(!LocateFirst(&fin)){
            Reset();
            fin.close();
            return -2;
         }
         fin.getline(buff,500);
         int time_last=ProcessTime();
         if(time_last<time_new){
            Reset();
            fin.close();
            return -3;
         }
         long int posi=pose-(time_last-time_new)*linelength;
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
   int timei=ProcessTime();
   if(jdebug>0) printf("RotateDB::LoadData: time_in=%d time_loop0=%d buff_loop0=%s\n",time_new,timei,buff);
   bool godown=timei<time_new;
   int nloop=0;
   while(timei!=time_new){
      if(godown){
         fin.getline(buff,500);
         timei=ProcessTime();
         if(!(fin.good())) break;
         if(timei>=time_new) break;
      }
      else{
         long int posi=posc-(long int)(1.5*linelength);
         if(posi>=poss){
            fin.seekg(-(long int)(1.5*linelength),std::ios::cur);
            LocateFirst(&fin);
            fin.getline(buff,500);
            timei=ProcessTime();
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
      Li=Li_in;
      time=time_in;
      currpos=pos0;
      fin.close();
      if(jdebug>0) printf("RotateDB::LoadData: loaded, Li=%d time=%d cpos=%ld\n",Li,time,currpos);
      return currpos;
   }
   else{
      Reset();
      fin.close();
      return -4;
   }
}
int RotateDB::ProcessTime(){
   int size=strlen(buff);
   char buff0[10];

   int count;
   int nchar1,nchar2,nchar3;
   //Get Time information
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
   return result;
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

int RotateDB::IsFineAngle(double ele_in,double azi_in,int iTel,int &index){
   const int nrot=2;
   int rotindex[nrot]={2,3};
   int irot=-1;
   for(int ii=0;ii<nrot;ii++){
      if(Li==rotindex[ii]) {irot=ii; break;}
   }
   if(irot<0) return 0;
   const int ntel=6;
   int telindex[ntel]={1,2,3,4,5,6};
   int itel=-1;
   for(int ii=0;ii<ntel;ii++){
      if(iTel==telindex[ii]) {itel=ii; break;}
   }
   if(itel<0) return 0;

   double margin=0.02;

   const int nangle1=4;
   double AngleList1[nangle1][2];
   AngleList1[0][0]=40; AngleList1[0][1]=42;
   AngleList1[1][0]=40; AngleList1[1][1]=29;
   AngleList1[2][0]=40; AngleList1[2][1]=19;
   AngleList1[3][0]=40; AngleList1[3][1]=4;
   for(int ii=0;ii<nangle1;ii++){
      if(fabs(ele_in-AngleList1[ii][0])<margin&&fabs(azi_in-AngleList1[ii][1])<margin){
         if(Li!=2) return 0;
         bool isfine=false;
         if(ii==0&&(telindex[itel]==5||telindex[itel]==6)) isfine=true;
         if(ii==1&&(telindex[itel]==4||telindex[itel]==5)) isfine=true;
         if(ii==2&&(telindex[itel]==3||telindex[itel]==4)) isfine=true;
         if(ii==3&&(telindex[itel]==1||telindex[itel]==2||telindex[itel]==3)) isfine=true;
         if(isfine) {index=ii; return 1;}
         else {index=-1; return 0;}
      }
   }

   const int nangle2=7;
   double elelist[nangle2]={10,20,30,40,50,55,60};
   double azilist[nrot][ntel][nangle2]={{{23.7,19,15,10,1,-5,-1000},
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
      if(fabs(ele_in-elelist[ii])<margin&&fabs(azi_in-azilist[irot][itel][ii])<margin){
         index=ii;
         return 2;
      }
   }
   return 0;
}
int RotateDB::GetEleAzi(int time_in,int Li_in,int iTel){
   if(jdebug>0) printf("RotateDB::GetEleAzi: timein=%d Lin=%d iTel=%d\n",time_in,Li_in,iTel);
   if(LoadData(time_in,Li_in,Li,time,currpos)<0) return -1;
   ProcessAll();
   RotateDB rotbuff(this);

   double ele0=GetElevation();
   double azi0=GetAzimuth();
   if(jdebug>0) printf("RotateDB::GetEleAzi: findtimelog ele=%lf azi=%lf\n",ele0,azi0);
   int index=-1;
   int retval=-1;
   if(iTel>0){
      retval=IsFineAngle(ele0,azi0,iTel,index);
      if(retval<=0){
         if(jdebug>0) printf("RotateDB::GetEleAzi: IsFineAngle=%d Index=%d\n",retval,index);
         return -2;
      }
   }
   if(jdebug>0) printf("RotateDB::GetEleAzi: IsFineAngle=%d Index=%d\n",retval,index);

   int pLi=rotbuff.Li;
   int ctime=rotbuff.time;
   long int cpos=rotbuff.currpos;

   double margin=0.02;
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
      if(fabs(ele0-elei)>margin) break;
      if(fabs(azi0-azii)>margin) break;
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
      if(fabs(ele0-elei)>margin) break;
      if(fabs(azi0-azii)>margin) break;
      ncount2++;
   }

   //LoadData(time_in,Li_in);
   //ProcessAll();
   this->Copy(&rotbuff);

   //double type0[4][2]={{40,42},{40,29},{40,19},{40,4}};
   //bool Istype0=false;
   //for(int ii=0;ii<4;ii++){
   //   if(fabs(ele0-type0[ii][0])<margin&&fabs(azi0-type0[ii][1])<margin) Istype0=true;
   //}

   int ncount=ncount1+ncount2+1;
   //int nside1=((ncount%2)==0)?((ncount/2)-13):(((ncount-1)/2)-12);
   //int nside2=((ncount%2)==0)?((ncount/2)-13):(((ncount-1)/2)-13);
   //if(jdebug>0) printf("RotateDB::GetEleAzi: ncount={%d,%d} nside={%d,%d}\n",ncount1,ncount2,nside1,nside2);

   if(ncount<ntotmin) return -3;
   if(ncount1<=nsidemin||ncount2<=nsidemin) return -4;
   else{
      return retval*10+index;
   }
}
