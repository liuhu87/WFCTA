#include "RotateDB.h"
#include "stdlib.h"
#include "WFCTAEvent.h"
#include "TelGeoFit.h"
#include <fstream>
#include <string>
#include <stdio.h>
RotateDB* RotateDB::_Head=0;
char RotateDB::RotateLogDir[200]="/eos/lhaaso/raw/wfctalaser";
int RotateDB::jdebug=0;
bool RotateDB::UseGPSTime=true;
bool RotateDB::UseGPSInfo=false;
bool RotateDB::UseEnvInfo=false;
bool RotateDB::UsePltInfo=false;
int RotateDB::ntotmin=5;
int RotateDB::nsidemin=2;
int RotateDB::nrot=2;
int RotateDB::rotindex[NRotMax]={2,3,0,0,0};
int RotateDB::TargetIndex=-1;
int RotateDB::timedelay[NRotMax]={37,37,0,0,0};
double RotateDB::rottime[NRotMax]={990016000,990845000,0,0,0};
int RotateDB::ntel=6;
int RotateDB::telindex[20]={1,2,3,4,5,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
double RotateDB::aglmargin=0.02;
double RotateDB::phimargin=10.;
double RotateDB::ccmargin=1.;
double RotateDB::phismargin=4.;
double RotateDB::ccsmargin=4.;
bool RotateDB::SearchAllAngle=false;
RotateDB::RotateDB(RotateDB* pr_in){
   Copy(pr_in);
}
void RotateDB::Init(){
   for(int ii=0;ii<500;ii++){
      rotateinfo[ii]='\0';
      rotateinfo2[0][ii]='\0';
      rotateinfo2[1][ii]='\0';
      pltinfo[0][ii]='\0';
      pltinfo[1][ii]='\0';
      gpsinfo[0][ii]='\0';
      gpsinfo[1][ii]='\0';
      envinfo[0][ii]='\0';
      envinfo[1][ii]='\0';
   }
   Li=0;
   time=0;

   for(int ii=0;ii<NRotMax;ii++){
      filesize1[ii]=0;
      filesize2[ii]=0;
   }
   Reset();
}
void RotateDB::Reset(){
   latitude=0;
   longitude=0;
   altitude=0;
   for(int ii=0;ii<NSWITH;ii++) allswith[ii]=false;
   for(int ii=0;ii<NVALUE;ii++) varinfo[ii]=0;
}
void RotateDB::Release(){
   for(int ii=0;ii<NRotMax;ii++){
      if(fin_log1[ii].is_open()) fin_log1[ii].close();
      if(fin_log2[ii].is_open()) fin_log2[ii].close();
      if(fin_gps[ii].is_open()) fin_gps[ii].close();
      if(fin_env[ii].is_open()) fin_env[ii].close();
   }
}
void RotateDB::Copy(RotateDB* pr_in){
   //info
   strcpy(rotateinfo,pr_in->rotateinfo);
   strcpy(rotateinfo2[0],pr_in->rotateinfo2[0]);
   strcpy(rotateinfo2[1],pr_in->rotateinfo2[1]);
   strcpy(pltinfo[0],pr_in->pltinfo[0]);
   strcpy(pltinfo[1],pr_in->pltinfo[1]);
   strcpy(gpsinfo[0],pr_in->gpsinfo[0]);
   strcpy(gpsinfo[1],pr_in->gpsinfo[1]);
   strcpy(envinfo[0],pr_in->envinfo[0]);
   strcpy(envinfo[1],pr_in->envinfo[1]);
   //filename
   for(int ii=0;ii<NRotMax;ii++){
      strcpy(filename_log1[ii],pr_in->filename_log1[ii]);
      strcpy(filename_log2[ii],pr_in->filename_log2[ii]);
      strcpy(filename_gps[ii],pr_in->filename_gps[ii]);
      strcpy(filename_env[ii],pr_in->filename_env[ii]);
      filesize1[ii]=pr_in->filesize1[ii];
      filesize2[ii]=pr_in->filesize2[ii];
      //processing the files
      if(pr_in->fin_log1[ii].is_open()){
         if(fin_log1[ii].is_open()) fin_log1[ii].close();
         fin_log1[ii].open(filename_log1[ii],std::ios::in);
         fin_log1[ii].seekg((long int)pr_in->fin_log1[ii].tellg(),std::ios::beg);
      }
      if(pr_in->fin_log2[ii].is_open()){
         if(fin_log2[ii].is_open()) fin_log2[ii].close();
         fin_log2[ii].open(filename_log2[ii],std::ios::in);
         fin_log2[ii].seekg((long int)pr_in->fin_log2[ii].tellg(),std::ios::beg);
      }
      if(pr_in->fin_gps[ii].is_open()){
         if(fin_gps[ii].is_open()) fin_gps[ii].close();
         fin_gps[ii].open(filename_gps[ii],std::ios::in);
         fin_gps[ii].seekg((long int)pr_in->fin_gps[ii].tellg(),std::ios::beg);
      }
      if(pr_in->fin_env[ii].is_open()){
         if(fin_env[ii].is_open()) fin_env[ii].close();
         fin_env[ii].open(filename_log2[ii],std::ios::in);
         fin_env[ii].seekg((long int)pr_in->fin_env[ii].tellg(),std::ios::beg);
      }
   }

   Li=pr_in->Li;
   time=pr_in->time;
   latitude=pr_in->latitude;
   longitude=pr_in->longitude;
   altitude=pr_in->altitude;
   for(int ii=0;ii<NSWITH;ii++) allswith[ii]=pr_in->allswith[ii];
   for(int ii=0;ii<NVALUE;ii++) varinfo[ii]=pr_in->varinfo[ii];
}
char* RotateDB::GetDirName(){
   return RotateLogDir;
}
bool RotateDB::SetDirName(char* dirname){
   if(!dirname) return false;
   else{
      if(strcpy(RotateLogDir,dirname)) return true;
      else return false;
   }
}
RotateDB* RotateDB::GetHead(){
   if(!_Head) _Head=new RotateDB();
   return _Head;
}
RotateDB* RotateDB::GetHead(char* dirname){
   SetDirName(dirname);
   return GetHead();
}
bool RotateDB::LocateFirst(ifstream* fin){
   if(!fin) return false;
   if(!fin->is_open()) return false;
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
bool RotateDB::LocateEnd(ifstream* fin){
   if(!fin) return false;
   if(!fin->is_open()) return false;
   if(!fin->good()) return false;
   long int pos0=fin->tellg();
   fin->seekg(0,std::ios::end);
   long int pose=fin->tellg();
   fin->seekg(pos0,std::ios::beg);
   char current='0';
   while(!(current=='\n')){
      long int posi=fin->tellg();
      if(posi<0) return false;
      if(posi>=pose) return true;
      current=fin->get();
   }
   if(current=='\n'){
      return true;
   }
   else return false;
}
bool RotateDB::GetFirstLastLine(ifstream* fin,char* firstline,char* lastline){
   if(!fin) return false;
   if(!fin->good()) return false;
   if(!firstline) return false;
   if(!lastline) return false;
   //int size1=sizeof(firstline)/sizeof(firstline[0]);
   //int size2=sizeof(lastline)/sizeof(lastline[0]);
   //if(size1<250||size2<250){
   //   cerr<<"RotateDB::GetFirstLastLine: "<<"size too small. size1="<<size1<<" size2="<<size2<<endl;
   //   printf("firstline=%s,lastline=%s\n",firstline,lastline);
   //   return false;
   //}
   long int cpos=fin->tellg();
   fin->seekg(0,std::ios::beg);
   fin->getline(firstline,250);
   if(!fin->good()){
      fin->seekg(cpos,std::ios::beg);
      return false;
   }
   fin->seekg(-5,std::ios::end);
   if(!LocateFirst(fin)){
      fin->seekg(cpos,std::ios::beg);
      return false;
   }
   fin->getline(lastline,250);
   if(!fin->good()){
      fin->seekg(cpos,std::ios::beg);
      return false;
   }
   fin->seekg(cpos,std::ios::beg);
   return true;
}
bool RotateDB::ReadData(ifstream* fin,int godown,char* content){
   if(!fin) return false;
   if(!content) return false;
   if(!fin->is_open()) return false;
   if(fin->tellg()<0) return false;
   if(godown>0&&fin->eof()) return false;
   if((godown<0)&&fin->tellg()<=0) return false;
   long int cpos=fin->tellg();
   bool exist=false;
   if(godown>0){
      //printf("pos_before=%ld\n",fin->tellg());
      fin->seekg(0,std::ios::end);
      long int pose=fin->tellg();
      fin->seekg(cpos,std::ios::beg);
      if(pose-cpos>2){
         fin->getline(content,300);
         exist=true;
      }
      //printf("pos_after=%ld\n",fin->tellg());
      return exist;
   }
   else{
      //printf("pos_before=%ld\n",fin->tellg());
      if(cpos-2>=0){
         fin->seekg(-2,std::ios::cur);
         LocateFirst(fin);
         if(godown==0){
            fin->getline(content,300);
            return true;
         }
         else{
            cpos=fin->tellg();
            if(cpos-2>=0){
               fin->seekg(-2,std::ios::cur);
               LocateFirst(fin);
               fin->getline(content,300);
               return true;
            }
            else{
               fin->seekg(0,std::ios::beg);
               return false;
            }
         }
      }
      else{
         fin->seekg(0,std::ios::beg);
         if(godown==0){
            fin->getline(content,300);
            return true;
         }
         else return false;
      }
   }
}
bool RotateDB::ReadAData(ifstream* fin,int godown,char* content,int Type,bool cur_godown){
   if(!content) return false;
   bool IsCType=true;
   int ctime=0;
   do{
      bool doread=ReadData(fin,godown,content);
      //if(jdebug>0&&fin) printf("RotateDB::ReadData2: pos=%ld godown=%d content=%s\n",(long int)fin->tellg(),godown,content);
      if(!doread) break;
      else{
         IsCType=IsLogFine(content,Type);
         if(IsCType) ctime=abs(ProcessTime2(content));
      }
      if(godown==0&&(!IsCType)) godown=(cur_godown?1:-1);
   }
   while(!IsCType);
   return (ctime>1300000000);
}
int RotateDB::GetLogTime(int time_in,int Li_in){
   int irot=GetLi(Li_in);
   if(irot<0) return time_in;
   if(UseGPSTime){
      if(Li_in!=3) return time_in-timedelay[irot];
      else{
         if(time_in>1580800045&&time_in<1587295800){
            return time_in-(133+8); //L3 using version2, the log need to be corrected
         }
         else return time_in-timedelay[irot];
      }
   }
   else time_in;
}
bool RotateDB::IsLogFine(char* buffer,int RotateType){
   if(!buffer) return false;
   char tokens[5][20]={"","","","",""};
   sscanf(buffer,"%s %s %s %s %s",tokens[0],tokens[1],tokens[2],tokens[3],tokens[4]);
   /////////////////
   //RotateType=0: keyword:Plant
   //RotateType=1: keyword:Rotate
   //RotateType=2: keyword:GPS
   //RotateType=3: keyword:Environment
   /////////////////
   if(RotateType==0) return (strcmp(tokens[0],"Plant")==0);
   else if(RotateType==1) return (strcmp(tokens[0],"Rotate")==0);
   else if(RotateType==2) return (strcmp(tokens[0],"GPS")==0);
   else if(RotateType==3) return (strcmp(tokens[0],"Environment")==0);
   else return false;
}
bool RotateDB::IsLogFine(int RotateType){
   char tokens1[5][20]={"","","","",""};
   char tokens2[5][20]={"","","","",""};
   if(RotateType==0){
      sscanf(pltinfo[0],"%s %s %s %s %s",tokens1[0],tokens1[1],tokens1[2],tokens1[3],tokens1[4]);
      sscanf(pltinfo[1],"%s %s %s %s %s",tokens2[0],tokens2[1],tokens2[2],tokens2[3],tokens2[4]);
   }
   else if(RotateType==1){
      sscanf(rotateinfo2[0],"%s %s %s %s %s",tokens1[0],tokens1[1],tokens1[2],tokens1[3],tokens1[4]);
      sscanf(rotateinfo2[1],"%s %s %s %s %s",tokens2[0],tokens2[1],tokens2[2],tokens2[3],tokens2[4]);
   }
   else if(RotateType==2){
      sscanf(gpsinfo[0],"%s %s %s %s %s",tokens1[0],tokens1[1],tokens1[2],tokens1[3],tokens1[4]);
      sscanf(gpsinfo[1],"%s %s %s %s %s",tokens2[0],tokens2[1],tokens2[2],tokens2[3],tokens2[4]);
   }
   else if(RotateType==3){
      sscanf(envinfo[0],"%s %s %s %s %s",tokens1[0],tokens1[1],tokens1[2],tokens1[3],tokens1[4]);
      sscanf(envinfo[1],"%s %s %s %s %s",tokens2[0],tokens2[1],tokens2[2],tokens2[3],tokens2[4]);
   }
   else return false;
   bool result;
   //if(IsRotate) result=( (strcmp(tokens1[0],"Rotate")==0&&strcmp(tokens2[0],"Rotate")==0) && (strcmp(tokens1[4],"O")==0&&strcmp(tokens2[4],"O")==0) );
   if(RotateType==1) result=( (strcmp(tokens1[0],"Rotate")==0&&strcmp(tokens2[0],"Rotate")==0) );
   else result=true;
   if(jdebug>2) printf("RotateDB::IsLogFine: token1=%s token2=%s res=%d\n",tokens1[4],tokens2[4],result);
   return result;
}
bool RotateDB::FillInfo(int type,char* content1,char* content2){
   if(!content1) return false;
   if(!content2) return false;
   if(type==0){
      strcpy(pltinfo[0],content1);
      strcpy(pltinfo[1],content2);
      return true;
   }
   else if(type==1){
      strcpy(rotateinfo2[0],content1);
      strcpy(rotateinfo2[1],content2);
      return true;
   }
   else if(type==2){
      strcpy(gpsinfo[0],content1);
      strcpy(gpsinfo[1],content2);
      return true;
   }
   else if(type==3){
      strcpy(envinfo[0],content1);
      strcpy(envinfo[1],content2);
      return true;
   }
   else return false;
}
int RotateDB::ReadData(int Li_in,int godown,char* content){
   int irot=GetLi(Li_in);
   if(irot<0) return 0;
   char content0[500];
   if(ReadData(&(fin_log1[irot]),godown,content0)){
      int ctime=ProcessTime(content0);
      if(abs(ctime)>1300000000){
         if(content) strcpy(content,content0);
         return ctime;
      }
      else return 0;
   }
   else return 0;
}
int RotateDB::ReadAData(int Li_in,int godown,char* content,int Type,bool cur_godown){
   int irot=GetLi(Li_in);
   if(irot<0) return 0;
   char content0[500];
   ifstream* fbuff=0;
   if(Type==1) fbuff=&(fin_log2[irot]);
   else if(Type==2) fbuff=&(fin_gps[irot]);
   else if(Type==3) fbuff=&(fin_env[irot]);
   if(ReadAData(fbuff,godown,(char*)content0,Type,cur_godown)){
      int ctime=ProcessTime2(content0);
      if(abs(ctime)>1300000000){
         if(content) strcpy(content,content0);
         return ctime;
      }
      else return 0;
   }
   else return 0;
}
bool RotateDB::LoadData(int time_in,int Li_in){
   int irot=GetLi(Li_in);
   if(irot<0) return false;
   int time_new=GetLogTime(time_in,rotindex[irot]);
   int year=2000+(CommonTools::TimeFlag(time_new,1)%100);
   int month=CommonTools::TimeFlag(time_new,2);
   int day=CommonTools::TimeFlag(time_new,3);
   int hour=CommonTools::TimeFlag(time_new,4);
   int min=CommonTools::TimeFlag(time_new,5);
   int sec=CommonTools::TimeFlag(time_new,6);
   char filename[200]="";
   strcpy(filename,RotateLogDir);
   char* namebuff=Form("/Laser%02d/L%d_SCdata/L%d/%d-%02d-%02d.txt.utf8",Li_in,Li_in,Li_in,year,month,day);
   strcat(filename,namebuff);
   if(strcmp(filename_log1[irot],filename)!=0){
      if(fin_log1[irot].is_open()) fin_log1[irot].close();
      strcpy(filename_log1[irot],filename);
      fin_log1[irot].open(filename_log1[irot],std::ios::in);
      fin_log1[irot].seekg(0,std::ios::end);
      filesize1[irot]=(long int)fin_log1[irot].tellg();
      fin_log1[irot].seekg(0,std::ios::beg);
      if(jdebug>0) printf("RotateDB::LoadData1: time_in=%d open new file %s filesize=%ld\n",time_new,filename_log1[irot],filesize1[irot]);
   }
   //read current, previous and next data
   int time_read=-1;
   char content[500];
   time_read=abs(ReadData(Li_in,0,content));
   if(jdebug>1) printf("RotateDB::LoadData1: time_in=%d read current item %s\n",time_new,content);
   long int pos_cur=fin_log1[irot].tellg();
   if(time_new==time_read){
      strcpy(rotateinfo,content);
      return true;
   }

   time_read=abs(ReadData(Li_in,1,content));
   if(jdebug>1) printf("RotateDB::LoadData1: time_in=%d read next item %s\n",time_new,content);
   long int pos_nex=fin_log1[irot].tellg();
   if(time_new==time_read){
      strcpy(rotateinfo,content);
      fin_log1[irot].seekg(pos_cur,std::ios::beg);
      return true;
   }
   int time_cur=time_read;

   fin_log1[irot].seekg(pos_cur,std::ios::beg);
   time_read=abs(ReadData(Li_in,-1,content));
   if(jdebug>1) printf("RotateDB::LoadData1: time_in=%d read previous item %s\n",time_new,content);
   long int pos_pre=fin_log1[irot].tellg();
   if(time_new==time_read){
      strcpy(rotateinfo,content);
      fin_log1[irot].seekg(pos_cur,std::ios::beg);
      return true;
   }

   //get first and last line
   char firstline[500];
   char lastline[500];
   GetFirstLastLine(&(fin_log1[irot]),firstline,lastline);
   int time0=abs(ProcessTime(firstline));
   int time1=abs(ProcessTime(lastline));
   if(time_new<time0||time_new>time1) return false;
   long int pos0,pos1;
   if(time_new<time_cur){
      pos0=0;
      pos1=pos_cur;
      time1=time_cur;
   }
   else{
      pos0=pos_cur;
      time0=time_cur;
      fin_log1[irot].seekg(-5,std::ios::end);
      LocateFirst(&(fin_log1[irot]));
      pos1=fin_log1[irot].tellg();
   }
   if(jdebug>1) printf("RotateDB::LoadData1: time_in=%d time_first=%d(%ld) time_last=%d(%ld)\n",time_new,time0,pos0,time1,pos1);
   bool exist=false;
   int nloop=0;
   while(nloop<1000){
      long int posi=(pos0+pos1)/2;
      fin_log1[irot].seekg(posi,std::ios::beg);
      LocateFirst(&(fin_log1[irot]));
      long int posi2=fin_log1[irot].tellg();
      if(posi2==pos0||posi2==pos1) break;
      fin_log1[irot].getline(content,500);
      int timei=abs(ProcessTime(content));
      if(jdebug>2) printf("RotateDB::LoadData1: time_in=%d nloop=%d timei=%d pos=(%ld,%ld,%ld)\n",time_new,nloop,timei,pos0,posi2,pos1);
      if(timei==time_new){
         strcpy(rotateinfo,content);
         exist=true;
         break;
      }
      else if(timei<time_new){
         pos0=posi2;
         time0=timei;
      }
      else{
         pos1=posi2;
         time1=timei;
      }
      nloop++;
   }
   if(jdebug>1) printf("RotateDB::LoadData1: time_in=%d exist=%d pos=%ld\n",time_new,exist,fin_log1[irot].tellg());
   return exist;
}
bool RotateDB::LoadAData(int time_in,int Li_in,int Type){
   int irot=GetLi(Li_in);
   if(irot<0) return false;
   int time_new=GetLogTime(time_in,rotindex[irot]);
   int year=2000+(CommonTools::TimeFlag(time_new,1)%100);
   int month=CommonTools::TimeFlag(time_new,2);
   int day=CommonTools::TimeFlag(time_new,3);
   int hour=CommonTools::TimeFlag(time_new,4);
   int min=CommonTools::TimeFlag(time_new,5);
   int sec=CommonTools::TimeFlag(time_new,6);
   char filename[200]="";
   strcpy(filename,RotateLogDir);
   ifstream* filein=0;
   if(Type==0||Type==2){
      strcat(filename,Form("/Laser%02d/L%d_SCdata/%04d/info_%02d.txt",Li_in,Li_in,year,month));
      filein=&(fin_gps[irot]);
      if(strcmp(filename_gps[irot],filename)!=0){
         if(filein->is_open()) filein->close();
         strcpy(filename_gps[irot],filename);
         filein->open(filename_gps[irot],std::ios::in);
         if(jdebug>0) printf("RotateDB::LoadAData:Type=%d time_in=%d open new file %s filesize=xx\n",Type,time_new,filename_gps[irot]);
      }
   }
   if(Type==1){
      strcat(filename,Form("/Laser%02d/L%d_SCdata/%04d/%02d/log_%02d%02d.txt",Li_in,Li_in,year,month,month,day));
      filein=&(fin_log2[irot]);
      if(strcmp(filename_log2[irot],filename)!=0){
         if(filein->is_open()) filein->close();
         strcpy(filename_log2[irot],filename);
         filein->open(filename_log2[irot],std::ios::in);
         filein->seekg(0,std::ios::end);
         filesize2[irot]=(long int)filein->tellg();
         filein->seekg(0,std::ios::beg);
         if(jdebug>0) printf("RotateDB::LoadAData:Type=%d time_in=%d open new file %s filesize=%ld\n",Type,time_new,filename_gps[irot],filesize2[irot]);
      }
   }
   if(Type==3){
      strcat(filename,Form("/Laser%02d/L%d_SCdata/%04d/%02d/var_%02d%02d.txt",Li_in,Li_in,year,month,month,day));
      filein=&(fin_env[irot]);
      if(strcmp(filename_env[irot],filename)!=0){
         if(filein->is_open()) filein->close();
         strcpy(filename_env[irot],filename);
         filein->open(filename_env[irot],std::ios::in);
         if(jdebug>0) printf("RotateDB::LoadAData:Type=%d time_in=%d open new file %s filesize=xx\n",Type,time_new,filename_env[irot]);
      }
   }
   if(!filein) return false;
   //read current, previous and next data
   int time_pre=0;
   int time_cur=0;
   int time_nex=0;
   char content_pre[500];
   char content_cur[500];
   char content_nex[500];
   time_cur=abs(ReadAData(Li_in,0,(char*)content_cur,Type));
   long int pos_cur=filein->tellg();
   if(jdebug>1) printf("RotateDB::LoadAData:Type=%d Li=%d time_in=%d read current item time_cur=%d,%s\n",Type,Li_in,time_new,time_cur,content_cur);

   time_nex=abs(ReadAData(Li_in,1,(char*)content_nex,Type));
   long int pos_nex=filein->tellg();
   if(jdebug>1) printf("RotateDB::LoadAData:Type=%d Li=%d time_in=%d read next item time_nex=%d,%s\n",Type,Li_in,time_new,time_nex,content_nex);

   filein->seekg(pos_cur,std::ios::beg);
   time_pre=abs(ReadAData(Li_in,-1,content_pre,Type));
   long int pos_pre=filein->tellg();
   if(jdebug>1) printf("RotateDB::LoadAData:Type=%d Li=%d time_in=%d read previous item time_pre=%d,%s\n",Type,Li_in,time_new,time_pre,content_pre);

   //filein->seekg(pos_cur,std::ios::beg);
   if(time_new>=time_pre&&time_new<time_cur){
      FillInfo(Type,content_pre,content_cur);
      filein->seekg(pos_pre,std::ios::beg);
      if(jdebug>1) printf("RotateDB::LoadAData:Type=%d Li=%d time_in=%d exist=1 pos=%ld\n",Type,Li_in,time_new,filein->tellg());
      return true;
   }
   else if(time_new>=time_cur&&time_new<time_nex){
      FillInfo(Type,content_cur,content_nex);
      filein->seekg(pos_cur,std::ios::beg);
      if(jdebug>1) printf("RotateDB::LoadAData:Type=%d Li=%d time_in=%d exist=1 pos=%ld\n",Type,Li_in,time_new,filein->tellg());
      return true;
   }

   //get first and last line
   long int pos0,pos1;
   int time0,time1;
   char firstline[500];
   char lastline[500];
   filein->seekg(0,std::ios::beg);
   time0=abs(ReadAData(Li_in,1,firstline,Type));
   filein->seekg(-2,std::ios::cur);
   LocateFirst(filein);
   pos0=filein->tellg();
   filein->seekg(0,std::ios::end);
   time1=abs(ReadAData(Li_in,0,lastline,Type,false));
   filein->seekg(-2,std::ios::cur);
   LocateFirst(filein);
   pos1=filein->tellg();
   if(jdebug>1) printf("RotateDB::LoadAData:Type=%d time_in=%d time_first=%d pos_first=%ld, time_last=%d pos_last=%ld\n",Type,time_new,time0,pos0,time1,pos1);
   if(time0<1300000000||time1<1300000000) return false;

   //out of file range
   if(time_new<time0){
      //read previous file
      ifstream fin3;
      int pnum=1;
      do{
         char pfilename[200]="";
         strcpy(pfilename,RotateLogDir);
         int pyear,pmonth,pday;
         if(Type==0||Type==2){ //one month
            pyear=year;
            pmonth=month-1*pnum;
            if(pmonth<=0){
              pyear--;
              pmonth=(pmonth+12)%12;
            }
            strcat(pfilename,Form("/Laser%02d/L%d_SCdata/%04d/info_%02d.txt",Li_in,Li_in,pyear,pmonth));
         }
         else{
            int time_new_pre=time_new-24*3600*pnum;
            pyear=2000+(CommonTools::TimeFlag(time_new_pre,1)%100);
            pmonth=CommonTools::TimeFlag(time_new_pre,2);
            pday=CommonTools::TimeFlag(time_new_pre,3);
            strcat(pfilename,Form("/Laser%02d/L%d_SCdata/%04d/%02d/%s_%02d%02d.txt",Li_in,Li_in,pyear,pmonth,Type==1?"log":"var",pmonth,pday));
         }
         fin3.open(pfilename,std::ios::in);
         pnum++;
         if(pnum>10) break;
      }
      while(!fin3.is_open());
      if(!fin3.is_open()) return false;
      fin3.seekg(0,std::ios::end);
      char content[500];
      int ptime_last=abs(ReadAData(&fin3,0,content,Type,false));
      fin3.close();
      if(time_new<ptime_last) return false;
      FillInfo(Type,content,firstline);
      filein->seekg(0,std::ios::beg);
      if(jdebug>1) printf("RotateDB::LoadAData:Type=%d Li=%d time_in=%d, exist in previous file, pos=%ld\n",Type,Li_in,time_new,filein->tellg());
      return true;
   }
   else if(time_new>=time1){
      //read next file
      ifstream fin3;
      int nnum=1;
      do{
         char nfilename[200]="";
         strcpy(nfilename,RotateLogDir);
         int nyear,nmonth,nday;
         if(Type==0||Type==2){ //one month
            nyear=year;
            nmonth=month+1*nnum;
            if(nmonth>12){
              nyear++;
              nmonth=(nmonth-12)%12;
            }
            strcat(nfilename,Form("/Laser%02d/L%d_SCdata/%04d/info_%02d.txt",Li_in,Li_in,nyear,nmonth));
         }
         else{
            int time_new_nex=time_new+24*3600*nnum;
            int nyear=2000+(CommonTools::TimeFlag(time_new_nex,1)%100);
            int nmonth=CommonTools::TimeFlag(time_new_nex,2);
            int nday=CommonTools::TimeFlag(time_new_nex,3);
            strcat(nfilename,Form("/Laser%02d/L%d_SCdata/%04d/%02d/%s_%02d%02d.txt",Li_in,Li_in,nyear,nmonth,Type==1?"log":"var",nmonth,nday));
         }
         fin3.open(nfilename,std::ios::in);
         nnum++;
         if(nnum>10) break;
      }
      while(!fin3.is_open());
      if(!fin3.is_open()) return false;
      fin3.seekg(0,std::ios::beg);
      char content[500];
      int ntime_first=abs(ReadAData(&fin3,0,content,Type,true));
      fin3.close();
      if(time_new>=ntime_first) return false;
      FillInfo(Type,lastline,content);
      filein->seekg(0,std::ios::end);
      if(jdebug>1) printf("RotateDB::LoadAData:Type=%d Li=%d time_in=%d, exist in next file, pos=%ld\n",Type,Li_in,time_new,filein->tellg());
      return true;
   }

   if(time_new<time_pre){
      pos1=pos_cur;
      time1=time_nex;
   }
   else{
      pos0=pos_cur;
      time0=time_nex;
   }
   if(jdebug>1) printf("RotateDB::LoadAData:Type=%d Li=%d time_in=%d time_first=%d(%ld) time_last=%d(%ld)\n",Type,Li_in,time_new,time0,pos0,time1,pos1);
   bool exist=false;
   int nloop=0;
   bool searchdown=true;
   while(nloop<1000){
      long int posi=(pos0+pos1)/2;
      filein->seekg(posi,std::ios::beg);
      LocateEnd(filein);
      int timei=abs(ReadAData(Li_in,0,content_cur,Type,searchdown));
      long int posi1=filein->tellg();
      filein->seekg(-2,std::ios::cur);
      LocateFirst(filein);
      long int posi2=filein->tellg();
      if(jdebug>2) printf("RotateDB::LoadAData:Type=%d Li=%d time_in=%d nloop=%d timei=(%d,%d,%d) pos=(%ld,%ld,%ld,%ld)\n",Type,Li_in,time_new,nloop,time0,timei,time1,pos0,posi2,posi1,pos1);
      if(posi2==pos0){
         FillInfo(Type,firstline,lastline);
         exist=true;
         filein->seekg(posi1,std::ios::beg);
         break;
      }
      else if(posi2==pos1){
         if(!searchdown){
         FillInfo(Type,firstline,lastline);
         exist=true;
         filein->seekg(posi1,std::ios::beg);
         break;
         }
         else{
            searchdown=false;
            continue;
         }
      }
      if(timei<time0||timei>=time1) break;
      else if(timei<=time_new){
         pos0=posi2;
         time0=timei;
         strcpy(firstline,content_cur);
      }
      else{
         pos1=posi2;
         time1=timei;
         strcpy(lastline,content_cur);
      }
      nloop++;
   }
   if(jdebug>1) printf("RotateDB::LoadAData:Type=%d Li=%d time_in=%d exist=%d pos=%ld\n",Type,Li_in,time_new,exist,(long int)filein->tellg());
   return exist;
}
bool RotateDB::LoadRotate(int time_in,int Li_in){
   return LoadAData(time_in,Li_in,1);
}
bool RotateDB::LoadGPS(int time_in,int Li_in){
   return LoadAData(time_in,Li_in,2);
}
bool RotateDB::LoadEnv(int time_in,int Li_in){
   return LoadAData(time_in,Li_in,3);
}
bool RotateDB::LoadPlt(int time_in,int Li_in){
   return LoadAData(time_in,Li_in,0);
}
//long int RotateDB::LoadData(int time_in,int Li_in,int pLi,int ptime,long int cpos){
//   int irot=GetLi(Li_in);
//   if(irot<0) return -5;
//   int pirot=GetLi(pLi);
//   int time_new=GetLogTime(time_in,rotindex[irot]);
//   int year=2000+(CommonTools::TimeFlag(time_new,1)%100);
//   int month=CommonTools::TimeFlag(time_new,2);
//   int day=CommonTools::TimeFlag(time_new,3);
//   int hour=CommonTools::TimeFlag(time_new,4);
//   int min=CommonTools::TimeFlag(time_new,5);
//   int sec=CommonTools::TimeFlag(time_new,6);
//   char filename[200]="";
//   strcpy(filename,RotateLogDir);
//   char* namebuff=Form("/Laser%02d/L%d_SCdata/L%d/%d-%02d-%02d.txt.utf8",Li_in,Li_in,Li_in,year,month,day);
//   strcat(filename,namebuff);
//   if(jdebug>0) printf("RotateDB::LoadData: filename=%s\n",filename);
//   bool sameday=false;
//   if(ptime>0&&pLi>0&&pirot>=0){
//      int ptime_new=GetLogTime(ptime,rotindex[pirot]);
//      int pyear=2000+(CommonTools::TimeFlag(ptime_new,1)%100);
//      int pmonth=CommonTools::TimeFlag(ptime_new,2);
//      int pday=CommonTools::TimeFlag(ptime_new,3);
//      sameday=(Li_in==pLi)&&((year==pyear)&&(month==pmonth)&&(day==pday));
//   }
//
//   //read the file content
//   ifstream fin;
//   fin.open(filename,std::ios::in);
//   if(!fin.is_open()){
//      //Reset();
//      return -1;
//   }
//   fin.seekg(0,std::ios::beg);
//   long int poss=fin.tellg();
//   fin.getline(buff,500);
//   long int pos1=fin.tellg();
//   fin.seekg(0,std::ios::end);
//   long int pose=fin.tellg();
//   long int linelength=pos1-poss;
//   if(jdebug>0) printf("RotateDB::LoadData: file opened. pos={%ld,%ld,%ld},linelength=%ld\n",poss,pos1,pose,linelength);
//   fin.seekg(0,std::ios::beg);
//
//   //locate the initial position
//   bool located=false;
//   if(sameday&&cpos>=0){
//      fin.seekg(cpos,std::ios::beg);
//      if(((long int)fin.tellg())<0){
//         fin.close();
//         fin.open(filename,std::ios::in);
//      }
//      else located=true;
//   }
//   if(!located){
//      int time_first=abs(ProcessTime());
//      if(jdebug>2) printf("RotateDB::LoadData: time_new=%d time_first=%d buff=%s\n",time_new,time_first,buff);
//      if(time_first>time_new){
//         //Reset();
//         fin.close();
//         return -3;
//      }
//      fin.seekg(-10,std::ios::end);
//      int time_last;
//      if(!LocateFirst(&fin)){
//         //Reset();
//         fin.close();
//         return -2;
//      }
//      else{
//         fin.getline(buff,500);
//         time_last=abs(ProcessTime());
//         if(jdebug>2) printf("RotateDB::LoadData: time_new=%d time_last=%d buff=%s\n",time_new,time_last,buff);
//         if(time_last<time_new){
//            //Reset();
//            fin.close();
//            return -3;
//         }
//      }
//
//      if(hour<12){
//         long int posi=(time_new-time_first)*linelength;
//         if(jdebug>2) printf("RotateDB::LoadData: poss=%ld pose=%ld posi=%d deltaT1=%d\n",poss,pose,posi,(hour*3600+min*60+sec-time_first));
//         if(poss+posi>=pose){
//            fin.seekg(-10,std::ios::end);
//         }
//         else{
//            fin.seekg(posi,std::ios::beg);
//         }
//      }
//      else{
//         long int posi=pose-(time_last-time_new)*linelength;
//         if(jdebug>2) printf("RotateDB::LoadData: poss=%ld pose=%ld posi=%d deltaT2=%d\n",poss,pose,posi,(time_last-time_new));
//         if(posi<poss){
//            fin.seekg(0,std::ios::beg);
//         }
//         else{
//            fin.seekg(-(time_last-time_new)*linelength,std::ios::end);
//         }
//      }
//   }
//
//   LocateFirst(&fin);
//   if(jdebug>0) printf("RotateDB::LoadData: initial position located. pos=%ld(%d)\n",fin.tellg(),located);
//   fin.getline(buff,500);
//   long int posc=fin.tellg();
//   int timei=abs(ProcessTime());
//   if(jdebug>0) printf("RotateDB::LoadData: time_in=%d time_loop0=%d buff_loop0=%s\n",time_new,timei,buff);
//   bool godown=timei<time_new;
//   int nloop=0;
//   while(timei!=time_new){
//      if(godown){
//         fin.getline(buff,500);
//         timei=abs(ProcessTime());
//         if(!(fin.good())) break;
//         if(timei>=time_new) break;
//      }
//      else{
//         long int posi=posc-(long int)(1.5*linelength);
//         if(posi>=poss){
//            fin.seekg(-(long int)(1.5*linelength),std::ios::cur);
//            LocateFirst(&fin);
//            fin.getline(buff,500);
//            timei=abs(ProcessTime());
//            if(timei<=time_new) break;
//         }
//         else break;
//      }
//      if(jdebug>1) printf("RotateDB::LoadData: loop=%d(%d) time={%d,%d} buff=%s\n",nloop,godown,timei,time_new,buff);
//      nloop++;
//   }
//
//   if(jdebug>0) printf("RotateDB::LoadData: time_in=%d time_final=%d\n",time_new,timei);
//   if(timei==time_new){
//      fin.seekg(-(long int)(0.5*linelength),std::ios::cur);
//      LocateFirst(&fin);
//      long int pos0=fin.tellg();
//      version=1;
//      Li=Li_in;
//      time=time_in;
//      currpos=pos0;
//      fin.close();
//      if(jdebug>0) printf("RotateDB::LoadData: loaded, Li=%d time=%d cpos=%ld\n",Li,time,currpos);
//      return currpos;
//   }
//   else{
//      //Reset();
//      fin.close();
//      return -4;
//   }
//}
int RotateDB::ProcessTime(const char* content){
   if(!content) return 0;
   int size=strlen(content);

   int count;
   int nchar1,nchar2,nchar3;
   char hour[10],min[10],sec[10];
   for(int i1=0;i1<10;i1++){
      hour[i1]='\0';
      min[i1]='\0';
      sec[i1]='\0';
   }
   char year[10],month[10],day[10];
   for(int i1=0;i1<10;i1++){
      year[i1]='\0';
      month[i1]='\0';
      day[i1]='\0';
   }

   //Get Start Location
   int index0=-1;
   for(int ii=size-1;ii>=0;ii--){
      if(content[ii]==':'){
         index0=ii;
         break;
      }
   }
   int index_start=index0-16;
   int index_end=index0+2;
   if(!UseGPSTime){
      index_start=0;
      index_end=18;
   }
   if(index_start<0) return 0;
   bool IsGPSTime=(index_start>100);

   int index_time1=0;
   count=0;
   nchar1=0;nchar2=0;nchar3=0;
   for(int ii=index_start;ii<=TMath::Min(size-1,index_end);ii++){
      if(content[ii]=='-'||content[ii]==' '){
         count++;
         if(content[ii]==' ') {index_time1=ii+1; break;}
         continue;
      }
      if(count==0){
         year[nchar1++]=content[ii];
      }
      else if(count==1){
         month[nchar2++]=content[ii];
      }
      else if(count==2){
         day[nchar3++]=content[ii];
      }
   }
   count=0;
   nchar1=0;nchar2=0;nchar3=0;
   for(int ii=index_time1;ii<=TMath::Min(size-1,index_end);ii++){
      //char buff_string[4]; //the string contains two chars
      //for(int i2=0;i2<4;i2++) buff_string[i2]='\0';
      //for(int i2=0;i2<3;i2++){
      //   if(i2+ii<TMath::Min(35,size)) buff_string[i2]=content[i2+ii];
      //}
      //if(content[ii]==':' || (strstr(buff_string,"开")||strstr(buff_string,"关")) ){
      //   count++;
      //   if(strstr(buff_string,"开")||strstr(buff_string,"关")) {break;}
      //   continue;
      //}
      if(content[ii]==':'){
         count++;
         continue;
      }
      if(count==0){
         hour[nchar1++]=content[ii];
      }
      else if(count==1){
         min[nchar2++]=content[ii];
      }
      else if(count==2){
         sec[nchar3++]=content[ii];
      }
   }
   //printf("%s-%s-%s %s:%s:%s\n",year,month,day,hour,min,sec);

   double time0=atoi(year)*10000000000+atoi(month)*100000000+atoi(day)*1000000+atoi(hour)*10000+atoi(min)*100+atoi(sec);
   int result=CommonTools::Convert(time0);
   return IsGPSTime?result:(-result);
}
int RotateDB::ProcessTime(){
   return ProcessTime(rotateinfo);
}
bool RotateDB::ProcessAll(int time_in,int Li_in,bool update){
   int irot=GetLi(Li_in);
   if(irot<0&&irot>=nrot) return false;
   int time_new=GetLogTime(time_in,Li_in);
   if(!update){
      if(abs(time)==time_new&&Li==Li_in){
         if(jdebug>1) printf("RotateDB::ProcessAll: Info of time_in=%d Li_in=%d already exist\n",time_in,Li_in);
         return true;
      }
      else{
         if(jdebug>1) printf("RotateDB::ProcessAll: Info of time_in=%d Li_in=%d need to update\n",time_in,Li_in);
         return false;
      }
   }

   int time_buff=ProcessTime();
   if(abs(time_buff)<1300000000){
      cerr<<"RotateDB::ProcessAll: No Data Loaded for time="<<time_in<<" and Li="<<Li_in<<endl;
      return false;
   }
   time=time_buff;
   Li=Li_in;

   int size=strlen(rotateinfo);
   if(jdebug>0) printf("RotateDB::ProcessAll: time_in=%d Li_in=%d %s(size=%d)\n",time_in,Li_in,rotateinfo,size);

   int count;
   int nchar1,nchar2,nchar3;
   //Get Time information
   int index_time1=0;
   int index_time2=0;
   count=0;
   nchar1=0;nchar2=0;nchar3=0;
   for(int ii=size-1;ii>=0;ii--){
      if(rotateinfo[ii]==':'){
         count++;
         if(count==1) index_time2=ii+3;
         else if(count==2) index_time1=ii-4;
         continue;
      }
   }

   int index_day1=0;
   int index_day2=0;
   count=0;
   for(int ii=index_time1;ii>=0;ii--){
      if(rotateinfo[ii]=='-'){
         count++;
         if(count==1) index_day2=ii+3;
         else if(count==2) index_day1=((ii-5)>=0)?(ii-5):(ii-4);
         continue;
      }
   }

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
         if(rotateinfo[ii]==':'){
            index_start=ii+3;
            index_end=index_day1;
            break;
         }
      }

      count=0;
      nchar1=0;nchar2=0;nchar3=0;
      int nchar4=0,nchar5=0,nchar6=0,nchar7=0;
      for(int ii=index_time2;ii<size;ii++){
         //printf("nline=%d char=%c\n",nline,rotateinfo[ii]);
         if(rotateinfo[ii]=='\0') break;
         if(rotateinfo[ii]=='\n') break;
         char buff2[4];
         for(int i2=0;i2<4;i2++) buff2[i2]='\0';
         for(int i2=ii;i2<(ii+((count==0||count==3)?2:3));i2++){
            if(i2<size) buff2[i2-ii]=rotateinfo[i2];
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
            lon[0][nchar1++]=rotateinfo[ii];
         }
         else if(count==1){
            lon[1][nchar2++]=rotateinfo[ii];
         }
         else if(count==2){
            lon[2][nchar3++]=rotateinfo[ii];
         }
         else if(count==3){
            lat[0][nchar4++]=rotateinfo[ii];
         }
         else if(count==4){
            lat[1][nchar5++]=rotateinfo[ii];
         }
         else if(count==5){
            lat[2][nchar6++]=rotateinfo[ii];
         }
         else if(count==6){
            alt[nchar7++]=rotateinfo[ii];
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
         if(i2+ii<size) buff3[i2]=rotateinfo[i2+ii];
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
      if(rotateinfo[ii]=='\0') break;
      if(rotateinfo[ii]=='\n') break;
      if(rotateinfo[ii]=='.'){
         if(count<5){
            vars[count][nchar[count]]=rotateinfo[ii];
            nchar[count]=nchar[count]+1;
            vars[count][nchar[count]]=rotateinfo[ii+1];
            nchar[count]=nchar[count]+1;
            ii=ii+1;
         }
         else{
            vars[count][nchar[count]]=rotateinfo[ii];
            nchar[count]=nchar[count]+1;
            vars[count][nchar[count]]=rotateinfo[ii+1];
            nchar[count]=nchar[count]+1;
            vars[count][nchar[count]]=rotateinfo[ii+2];
            nchar[count]=nchar[count]+1;
            ii=ii+2;
         }
         count++;
         if(count==NVALUE) break;
         continue;
      }
      vars[count][nchar[count]]=rotateinfo[ii];
      nchar[count]=nchar[count]+1;
   }
   //printf("vars:");
   //for(int ii=0;ii<NVALUE;ii++) printf("%s,",vars[ii]);
   //printf("\n\n");
   for(int ii=0;ii<NVALUE;ii++) varinfo[ii]=atof(vars[ii]);
   if(jdebug>1) printf("RotateDB::ProcessAll: Ele=%.2lf Ai=%.2lf\n",varinfo[9],varinfo[10]);

   //some changes due to the change of slow control software
   if(abs(time)<1573541100){
      varinfo[7]*=(-1);
      varinfo[9]*=(-1);
   }
   else{
      varinfo[10]*=(-1);
   }
   if(abs(time)<1573541855) varinfo[8]+=511.03;

   return true;
}
//int RotateDB::ReadData(ifstream* fin,bool godown,int RotateType){
//   int result=0;
//   if(!fin) return result;
//   if(RotateType<0||RotateType>3) return result;
//   if(!LocateFirst(fin)) return result;
//   long int pos0=fin->tellg();
//   if(pos0<0) return result;
//   long int posres=pos0;
//   fin->seekg(0,std::ios::beg);
//   long int poss=fin->tellg();
//   fin->seekg(0,std::ios::end);
//   long int pose=fin->tellg();
//   char buffer[200];
//   //char target[20];
//   /////////////////
//   //RotateType=0: keyword:Plant
//   //RotateType=1: keyword:Rotate
//   //RotateType=2: keyword:GPS
//   //RotateType=2: keyword:Environment
//   /////////////////
//   //char target0[4][20]={"Plant","Rotate","GPS","Environment"};
//   //strcpy(target,target0[RotateType]);
//   int targetlength[4]={45,64,52,56};
//   long int posp,posc;
//   long int linelength=targetlength[RotateType];
//   if(godown){
//      fin->seekg(pos0,std::ios::beg);
//      posp=fin->tellg();
//      posres=posp;
//      if(posp<0||posp>=pose){
//         fin->seekg(posres,std::ios::beg);
//         if(jdebug>6) printf("RotateDB::ReadData: return0=%d\n",result);
//         return result;
//      }
//      fin->getline(buff2[0],500);
//      posc=fin->tellg();
//      strcpy(buffer,buff2[0]);
//      if(jdebug>6) printf("RotateDB::ReadData: item0 godown=%d pos={%ld,%ld} posse={%ld,%ld} buffer=%s\n",godown,posp,posc,poss,pose,buff2[0]);
//      while(!IsLogFine(buffer,RotateType)){
//         posp=fin->tellg();
//         if(posp<0||posp>=pose) break;
//         fin->getline(buff2[0],500);
//         strcpy(buffer,buff2[0]);
//         posc=fin->tellg();
//         if(jdebug>7) printf("RotateDB::ReadData: item0i godown=%d pos={%ld,%ld} posse={%ld,%ld} buffer=%s\n",godown,posp,posc,poss,pose,buff2[0]);
//      }
//      if(IsLogFine(buffer,RotateType)){
//         linelength=posc-posp;
//         result=(result|0x1);
//      }
//      posres=posc;
//      if(posc<0||posc>=pose){
//         fin->seekg(posres,std::ios::beg);
//         if(jdebug>6) printf("RotateDB::ReadData: return1=%d\n",result);
//         return result;
//      }
//      fin->getline(buff2[1],500);
//      strcpy(buffer,buff2[1]);
//      posc=fin->tellg();
//      if(jdebug>6) printf("RotateDB::ReadData: item1 godown=%d pos={%ld,%ld} posse={%ld,%ld} buffer=%s\n",godown,posp,posc,poss,pose,buff2[1]);
//      while(!IsLogFine(buffer,RotateType)){
//         posp=fin->tellg();
//         if(posp<0||posp>=pose) break;
//         fin->getline(buff2[1],500);
//         strcpy(buffer,buff2[1]);
//         posc=fin->tellg();
//         if(jdebug>7) printf("RotateDB::ReadData: item1i godown=%d pos={%ld,%ld} posse={%ld,%ld} buffer=%s\n",godown,posp,posc,poss,pose,buff2[1]);
//      }
//      if(IsLogFine(buffer,RotateType)){
//         result=(result|0x2);
//      }
//      else{
//         if(result|0x1) posres=posc;
//      }
//      fin->seekg(posres,std::ios::beg);
//      if(jdebug>6) printf("RotateDB::ReadData: return2=%d\n",result);
//      return result;
//   }
//   else{
//      fin->seekg(pos0,std::ios::beg);
//      posp=fin->tellg();
//      posres=posp;
//      if(posp-0.5*linelength<poss){
//         fin->seekg(posres,std::ios::beg);
//         if(jdebug>6) printf("RotateDB::ReadData: return3=%d\n",result);
//         return result;
//      }
//      fin->seekg(posp-0.5*linelength,std::ios::beg);
//      LocateFirst(fin);
//      posp=fin->tellg();
//      fin->getline(buff2[1],500);
//      posc=fin->tellg();
//      strcpy(buffer,buff2[1]);
//      if(jdebug>6) printf("RotateDB::ReadData: item0 godown=%d pos={%ld,%ld} posse={%ld,%ld} buffer=%s\n",godown,posp,posc,poss,pose,buff2[1]);
//      while(!IsLogFine(buffer,RotateType)){
//         if(posp<0||posp-0.5*linelength<poss) break;
//         fin->seekg(posp-0.5*linelength,std::ios::beg);
//         LocateFirst(fin);
//         posp=fin->tellg();
//         fin->getline(buff2[1],500);
//         strcpy(buffer,buff2[1]);
//         posc=fin->tellg();
//         if(jdebug>7) printf("RotateDB::ReadData: item0i godown=%d pos={%ld,%ld} posse={%ld,%ld} buffer=%s\n",godown,posp,posc,poss,pose,buff2[1]);
//      }
//      if(IsLogFine(buffer,RotateType)){
//         linelength=posc-posp;
//         result=(result|0x1);
//      }
//      posres=posp;
//      if(posp<0||posp-0.5*linelength<poss){
//         fin->seekg(posres,std::ios::beg);
//         if(jdebug>6) printf("RotateDB::ReadData: return4=%d\n",result);
//         return result;
//      }
//      fin->seekg(posp-0.5*linelength,std::ios::beg);
//      LocateFirst(fin);
//      posp=fin->tellg();
//      fin->getline(buff2[0],500);
//      strcpy(buffer,buff2[0]);
//      if(jdebug>6) printf("RotateDB::ReadData: item1 godown=%d pos={%ld,%ld} posse={%ld,%ld} buffer=%s\n",godown,posp,posc,poss,pose,buff2[0]);
//      while(!IsLogFine(buffer,RotateType)){
//         if(posp<0||posp-0.5*linelength<poss) break;
//         fin->seekg(posp-0.5*linelength,std::ios::beg);
//         LocateFirst(fin);
//         posp=fin->tellg();
//         fin->getline(buff2[0],500);
//         strcpy(buffer,buff2[0]);
//         if(jdebug>7) printf("RotateDB::ReadData: item1i godown=%d pos={%ld,%ld} posse={%ld,%ld} buffer=%s\n",godown,posp,posc,poss,pose,buff2[0]);
//      }
//      if(IsLogFine(buffer,RotateType)){
//         result=(result|0x2);
//      }
//      else{
//         if(result|0x1) posres=posp;
//      }
//      fin->seekg(posres,std::ios::beg);
//      if(jdebug>6) printf("RotateDB::ReadData: return5=%d\n",result);
//      return result;
//   }
//   return result;
//}
//int RotateDB::ReadData2(ifstream* fin,bool godown,bool IsRotate,int Li_in){
//   if(!fin) return 0;
//   char buffer0[2][500]={"",""};
//   strcpy(buffer0[0],buff2[0]);
//   strcpy(buffer0[1],buff2[1]);
//   int RotateType=IsRotate?1:2;
//   if(jdebug>5) printf("RotateDB::ReadData2: ReadData_before pos=%ld\n",fin->tellg());
//   int res=ReadData(fin,godown,RotateType);
//   if(jdebug>5) printf("RotateDB::ReadData2: ReadData readres=%d pos=%ld\n",res,fin->tellg());
//   if(res<=0){
//      strcpy(buff2[0],buffer0[0]);
//      strcpy(buff2[1],buffer0[1]);
//      if(jdebug>5) printf("RotateDB::ReadData2: return1=%d\n",res);
//      return res;
//   }
//   else if((res&0x1)==0){ //should not exist
//      strcpy(buff2[0],buffer0[0]);
//      strcpy(buff2[1],buffer0[1]);
//      if(jdebug>5) printf("RotateDB::ReadData2: return2=%d\n",0);
//      return 0;
//   }
//   else if((res&0x2)==0){ //should read the first or last from another file
//      if(godown){
//         int timei=abs(ProcessTime2(0));
//         int nyear=CommonTools::TimeFlag(timei+24*3600,1);
//         nyear=2000+(nyear%100);
//         int nmonth=CommonTools::TimeFlag(timei+24*3600,2);
//         int nday=CommonTools::TimeFlag(timei+24*3600,3);
//         char nfilename[200]="";
//         strcpy(nfilename,RotateLogDir);
//         strcat(nfilename,Form("/Laser%02d/L%d_SCdata/%d/%02d/log_%02d%02d.txt",Li_in,Li_in,nyear,nmonth,nmonth,nday));
//         ifstream fin2;
//         fin2.open(nfilename,std::ios::in);
//         char buffer1[500];
//         strcpy(buffer1,buff2[0]);
//         int res2=ReadData(&fin2,godown,RotateType);
//         if(jdebug>6) printf("RotateDB::ReadData2: read the log for next day. readres=%d godown=%d filename=%s\n",res2,godown,nfilename);
//         fin2.close();
//         if((res2&0x1)==0){
//            strcpy(buff2[0],buffer0[0]);
//            strcpy(buff2[1],buffer0[1]);
//            if(jdebug>5) printf("RotateDB::ReadData2: return3=%d\n",0);
//            return 0;
//         }
//         else{
//            int timei2=abs(ProcessTime2(0));
//            //if(abs(timei2-timei)>3600)
//            if(timei2<timei){
//               strcpy(buff2[0],buffer0[0]);
//               strcpy(buff2[1],buffer0[1]);
//               if(jdebug>5) printf("RotateDB::ReadData2: return4=%d\n",0);
//               return 0;
//            }
//            else{
//               strcpy(buff2[1],buff2[0]);
//               strcpy(buff2[0],buffer1);
//               if(jdebug>5) printf("RotateDB::ReadData2: return5=%d\n",(res|0x2));
//               return (res|0x2);
//            }
//         }
//      }
//      else{
//         int timei=abs(ProcessTime2(1));
//         int pyear=CommonTools::TimeFlag(timei-24*3600,1);
//         pyear=2000+(pyear%100);
//         int pmonth=CommonTools::TimeFlag(timei-24*3600,2);
//         int pday=CommonTools::TimeFlag(timei-24*3600,3);
//         char pfilename[200]="";
//         strcpy(pfilename,RotateLogDir);
//         strcat(pfilename,Form("Laser%02d/L%d_SCdata/%d/%02d/log_%02d%02d.txt",Li_in,Li_in,pyear,pmonth,pmonth,pday));
//         ifstream fin2;
//         fin2.open(pfilename,std::ios::in);
//         fin2.seekg(0,std::ios::end);
//         char buffer1[500];
//         strcpy(buffer1,buff2[0]);
//         int res2=ReadData(&fin2,godown,RotateType);
//         if(jdebug>6) printf("RotateDB::ReadData2: read the log for previous day. readres=%d godown=%d filename=%s\n",res2,godown,pfilename);
//         fin2.close();
//         if((res2&0x2)==0){
//            strcpy(buff2[0],buffer0[0]);
//            strcpy(buff2[1],buffer0[1]);
//            if(jdebug>5) printf("RotateDB::ReadData2: return6=%d\n",0);
//            return 0;
//         }
//         else{
//            int timei2=abs(ProcessTime2(1));
//            //if(abs(timei2-timei)>3600)
//            if(timei2>timei){
//               strcpy(buff2[0],buffer0[0]);
//               strcpy(buff2[1],buffer0[1]);
//               if(jdebug>5) printf("RotateDB::ReadData2: return7=%d\n",0);
//               return 0;
//            }
//            else{
//               strcpy(buff2[0],buff2[1]);
//               strcpy(buff2[1],buffer1);
//               if(jdebug>5) printf("RotateDB::ReadData2: return8=%d\n",(res|0x2));
//               return (res|0x2);
//            }
//         }
//      }
//      if(jdebug>5) printf("RotateDB::ReadData2: return9=%d\n",(res|0x2));
//      return (res|0x2);
//   }
//   else{
//      if(jdebug>5) printf("RotateDB::ReadData2: return10=%d\n",res);
//      return res;
//   }
//}
//long int RotateDB::LoadData2(int time_in,int Li_in,int pLi,int ptime1,int ptime2,long int cpos){
//   int ptime[2]={ptime1,ptime2};
//   int irot=GetLi(Li_in);
//   if(irot<0) return -5;
//   if(Li_in==pLi&&(ptime[0]>=0&&ptime[1]>=0)&&(time_in>=ptime[0]&&time_in<=ptime[1])){
//      if(cpos>=0) return cpos;
//   }
//   int pirot=GetLi(pLi);
//   int time_new=GetLogTime(time_in,rotindex[irot]);
//   int year=CommonTools::TimeFlag(time_new,1);
//   year=2000+(year%100);
//   int month=CommonTools::TimeFlag(time_new,2);
//   int day=CommonTools::TimeFlag(time_new,3);
//   int hour=CommonTools::TimeFlag(time_new,4);
//   int min=CommonTools::TimeFlag(time_new,5);
//   int sec=CommonTools::TimeFlag(time_new,6);
//
//   char filename[200]="";
//   strcpy(filename,RotateLogDir);
//   char* namebuff=Form("/Laser%02d/L%d_SCdata/%d/%02d/log_%02d%02d.txt",Li_in,Li_in,year,month,month,day);
//   strcat(filename,namebuff);
//   if(jdebug>0) printf("RotateDB::LoadData2: filename=%s\n",filename);
//   bool sameday=false;
//   if((ptime[0]>0&&ptime[1]>0)&&pLi>0&&pirot>=0){
//      int ptime_new=GetLogTime((ptime[0]+ptime[1])/2,rotindex[pirot]);
//      int pyear=CommonTools::TimeFlag(ptime_new,1);
//      pyear=2000+(pyear%100);
//      int pmonth=CommonTools::TimeFlag(ptime_new,2);
//      int pday=CommonTools::TimeFlag(ptime_new,3);
//      sameday=(Li_in==pLi)&&((year==pyear)&&(month==pmonth)&&(day==pday));
//   }
//
//   //read the file content
//   char tokens[20];
//   ifstream fin;
//   fin.open(filename,std::ios::in);
//   if(!fin.is_open()){
//      //Reset();
//      return -1;
//   }
//   fin.seekg(0,std::ios::beg);
//   long int poss=fin.tellg();
//   fin.seekg(0,std::ios::end);
//   long int pose=fin.tellg();
//
//   //locate the initial position
//   bool located=false;
//   if(sameday&&cpos>=0){
//      fin.seekg(cpos,std::ios::beg);
//      if(((long int)fin.tellg())<0){
//         fin.close();
//         fin.open(filename,std::ios::in);
//      }
//      located=true;
//   }
//   if(!located){
//      fin.seekg(0,std::ios::beg);
//   }
//
//   LocateFirst(&fin);
//   long int pos0=fin.tellg();
//   long int posp,posc;
//   if(jdebug>0) printf("RotateDB::LoadData2: initial position located. pos=%ld(%d)\n",fin.tellg(),located);
//   bool godown=true;
//   posp=fin.tellg();
//   int firstres=ReadData2(&fin,godown,true,Li_in);
//   posc=fin.tellg();
//   if(firstres<=0){
//      godown=false;
//      fin.seekg(pos0,std::ios::beg);
//      posp=fin.tellg();
//      firstres=ReadData2(&fin,godown,true,Li_in);
//      posc=fin.tellg();
//   }
//   if(firstres<=0) return -2;
//   int timei0=abs(ProcessTime2(0));
//   int timei1=abs(ProcessTime2(1));
//   if(timei1<timei0||timei0<1000000000) return -2;
//   godown=timei0<time_new;
//   if(jdebug>0) printf("RotateDB::LoadData2: time_in=%d time_loop0={%d,%d} buff_loop0={%s,%s}\n",time_new,timei0,timei1,buff2[0],buff2[1]);
//   bool loaded=true;
//   int nloop=0;
//   while(time_new<timei0||time_new>timei1){
//      posp=fin.tellg();
//      int fillres=ReadData2(&fin,godown,true,Li_in);
//      posc=fin.tellg();
//      if((fillres&0x3)<=0){
//         if(jdebug>1) printf("RotateDB::LoadData2: loop=%d(%d) time=%d buff={%s,%s}\n",nloop,godown,time_new,buff2[0],buff2[1]);
//         loaded=false; break;
//      }
//      else{
//         timei0=abs(ProcessTime2(0));
//         timei1=abs(ProcessTime2(1));
//         if(jdebug>1) printf("RotateDB::LoadData2: loop=%d(%d) time={%d,%d,%d} buff={%s,%s}\n",nloop,godown,timei0,timei1,time_new,buff2[0],buff2[1]);
//         if((time_new>=timei0&&time_new<=timei1)&&(!IsLogFine(true))) {loaded=false; break;}
//      }
//      nloop++;
//   }
//
//   if(jdebug>0) printf("RotateDB::LoadData2: time_in=%d time_final={%d,%d}\n",time_new,timei0,timei1);
//   if(loaded&&(time_new>=timei0&&time_new<=timei1)){
//      version=2;
//      Li2=Li_in;
//      time2[0]=timei0;
//      time2[1]=timei1;
//      currpos2=(posp+posc)/2;
//      fin.close();
//      if(jdebug>0) printf("RotateDB::LoadData2: loaded, Li=%d time_in=%d time={%d,%d} cpos=%ld posse={%ld,%ld} buff={%s,%s}\n",time_new,Li2,time2[0],time2[1],currpos2,poss,pose,buff2[0],buff2[1]);
//      return currpos2;
//   }
//   else{
//      //Reset();
//      fin.close();
//      return -4;
//   }
//}
int RotateDB::ProcessTime2(char* content){
   if(!content) return 0;

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
   sscanf(content, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s",tokens[0],tokens[1],tokens[2],tokens[3],tokens[4],tokens[5],tokens[6],tokens[7],tokens[8],tokens[9],tokens[10],tokens[11],tokens[12],tokens[13]);

   //printf("tokens=%s\n",tokens[0]);
   if(strcmp(tokens[0],"Rotate")==0){
      strcpy(date0,tokens[1]);
      strcpy(time0,tokens[2]);
   }
   else if(strcmp(tokens[0],"Environment")==0){
      strcpy(date0,tokens[1]);
      strcpy(time0,tokens[2]);
   }
   else if(strcmp(tokens[0],"GPS")==0){
      strcpy(date0,tokens[2]);
      strcpy(time0,tokens[3]);
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
int RotateDB::ProcessTime2(int Type,int itime){
   if(itime<0||itime>1) return 0;
   if(Type==0) return ProcessTime2(pltinfo[itime]);
   else if(Type==1) return ProcessTime2(rotateinfo2[itime]);
   else if(Type==2) return ProcessTime2(gpsinfo[itime]);
   else if(Type==3) return ProcessTime2(envinfo[itime]);
   else return 0;
}
bool RotateDB::ProcessAll2(int time_in,int Li_in,bool update){
   int irot=GetLi(Li_in);
   if(irot<0&&irot>=nrot) return false;
   int time_new=GetLogTime(time_in,Li_in);
   if(!update){
      if(abs(time)==time_new&&Li==Li_in){
         if(jdebug>1) printf("RotateDB::ProcessAll2: Info of time_in=%d Li_in=%d already exist\n",time_in,Li_in);
         return true;
      }
      else{
         if(jdebug>1) printf("RotateDB::ProcessAll2: Info of time_in=%d Li_in=%d need to update\n",time_in,Li_in);
         return false;
      }
   }

   int ninfo=0;
   for(int type=0;type<4;type++){
      int time0=abs(ProcessTime2(type,0));
      int time1=abs(ProcessTime2(type,1));
      if(time0<1300000000||time1<1300000000) continue;
      if(time0>time1) continue;
      if(time_new<time0||time_new>=time1) continue;
      ninfo++;
   }
   if(ninfo<=0){
      cerr<<"RotateDB::ProcessAll2: No Data Loaded for time="<<time_in<<" and Li="<<Li_in<<endl;
      return false;
   }
   time=time_new;
   Li=Li_in;
   ProcessRotate();
   ProcessGPS();
   ProcessEnv();
   ProcessPlant();
   return true;
}
bool RotateDB::ProcessRotate(){
   int time0=ProcessTime2(1,0);
   int time1=ProcessTime2(1,1);
   if(time0<1300000000||time1<1300000000) return false;
   if(time0>time1) return false;
   if(time<time0||time>=time1) return false;

   int size=strlen(rotateinfo2[0]);
   if(jdebug>0) printf("RotateDB::ProcessRotate: %s(size=%d)\n",rotateinfo2[0],size);
   char tokens[14][20];
   sscanf(rotateinfo2[0],"%s %s %s %s %s %s %s %s %s %s %s %s %s %s",tokens[0],tokens[1],tokens[2],tokens[3],tokens[4],tokens[5],tokens[6],tokens[7],tokens[8],tokens[9],tokens[10],tokens[11],tokens[12],tokens[13]);
   if(strcmp(tokens[0],"Rotate")==0){
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
      if(abs(time)<1573541100){
         varinfo[7]*=(-1);
         varinfo[9]*=(-1);
      }
      else{
         varinfo[10]*=(-1);
      }
      if(abs(time)<1573541855) varinfo[8]+=511.03;
   }
   return true;
}
bool RotateDB::ProcessGPS(){
   int time0=ProcessTime2(2,0);
   int time1=ProcessTime2(2,1);
   if(time0<1300000000||time1<1300000000) return false;
   if(time0>time1) return false;
   if(time<time0||time>=time1) return false;

   int size=strlen(gpsinfo[0]);
   if(jdebug>0) printf("RotateDB::ProcessGPS: %s(size=%d)\n",gpsinfo[0],size);
   char tokens[14][20];
   sscanf(gpsinfo[0],"%s %s %s %s %s %s %s %s %s %s %s %s %s %s",tokens[0],tokens[1],tokens[2],tokens[3],tokens[4],tokens[5],tokens[6],tokens[7],tokens[8],tokens[9],tokens[10],tokens[11],tokens[12],tokens[13]);
   int count;
   int nchar1,nchar2,nchar3;
   if(strcmp(tokens[0],"GPS")==0){
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
   return true;
}
bool RotateDB::ProcessEnv(){
   int time0=ProcessTime2(3,0);
   int time1=ProcessTime2(3,1);
   if(time0<1300000000||time1<1300000000) return false;
   if(time0>time1) return false;
   if(time<time0||time>=time1) return false;

   int size=strlen(envinfo[0]);
   if(jdebug>0) printf("RotateDB::ProcessEnv: %s(size=%d)\n",envinfo[0],size);
   char tokens[2][3][20];
   double vars12[2][6];
   sscanf(envinfo[1],"%s %s %s %lf %lf %lf %lf %lf",tokens[1][0],tokens[1][1],tokens[1][2],&vars12[1][0],&vars12[1][1],&vars12[1][2],&vars12[1][3],&vars12[1][4]);
   sscanf(envinfo[0],"%s %s %s %lf %lf %lf %lf %lf",tokens[0][0],tokens[0][1],tokens[0][2],&vars12[0][0],&vars12[0][1],&vars12[0][2],&vars12[0][3],&vars12[0][4]);
   if(strcmp(tokens[0][0],"Environment")==0&&strcmp(tokens[1][0],"Environment")==0){
      for(int ii=0;ii<5;ii++){
         double ratio=(time-time0)/(time1-time0);
         varinfo[ii]=vars12[0][ii]*(1-ratio)+vars12[1][ii]*ratio;
      }
   }
   return true;
}
bool RotateDB::ProcessPlant(){
   int time0=ProcessTime2(0,0);
   int time1=ProcessTime2(0,1);
   if(time0<1300000000||time1<1300000000) return false;
   if(time0>time1) return false;
   if(time<time0||time>=time1) return false;

   int size=strlen(pltinfo[0]);
   if(jdebug>0) printf("RotateDB::ProcessGPS: %s(size=%d)\n",pltinfo[0],size);
   return true;
}
int RotateDB::ProcessEnv(int time_in,int Li_in){
   int irot=GetLi(Li_in);
   if(irot<0) return -4;
   int time_new=GetLogTime(time_in,rotindex[irot]);
   int year=CommonTools::TimeFlag(time_new,1);
   year=2000+(year%100);
   int month=CommonTools::TimeFlag(time_new,2);
   int day=CommonTools::TimeFlag(time_new,3);
   char filename[200]="";
   strcpy(filename,RotateLogDir);
   char* namebuff=Form("/Laser%02d/L%d_SCdata/%04d/%02d/var_%02d%02d.txt",Li_in,Li_in,year,month,month,day);
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
   printf("RotateDB::DumpInfo: rotateinfo=%s\n",rotateinfo);
   printf("rotateinfo2[0]=%s\n",rotateinfo2[0]);
   printf("rotateinfo2[1]=%s\n",rotateinfo2[1]);
   printf("gpsinfo[0]=%s\n",gpsinfo[0]);
   printf("gpsinfo[1]=%s\n",gpsinfo[1]);
   printf("envinfo[0]=%s\n",envinfo[0]);
   printf("envinfo[1]=%s\n",envinfo[1]);
   printf("pltinfo[0]=%s\n",pltinfo[0]);
   printf("pltinfo[1]=%s\n",pltinfo[1]);
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

bool RotateDB::GetEleAzi(int index,double &elevation,double &azimuth){
   if(abs(index)<100) return false;
   index=abs(index);
   int iTel=index/1000;
   int itel=-1;
   for(int ii=0;ii<ntel;ii++){
      if(telindex[ii]==iTel){itel=ii; break;}
   }
   if(itel<0) return false;
   int index0=(index%1000);
   int iRot=index0/100;
   int irot=GetLi(iRot);
   if(irot<0) return false;
   int itype=(index0%100)/10;
   int iangle=(index0%10);
   if(itype==1){
      const int nangle1=4;
      double AngleList1[nangle1][2];
      AngleList1[0][0]=40; AngleList1[0][1]=42;
      AngleList1[1][0]=40; AngleList1[1][1]=29;
      AngleList1[2][0]=40; AngleList1[2][1]=19;
      AngleList1[3][0]=40; AngleList1[3][1]=4;
      if(iangle>3) return false;
      elevation=AngleList1[iangle][0];
      azimuth=AngleList1[iangle][1];
      return true;
   }
   else if(itype==2){
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
      if(irot>1) return false;
      if(itel>5) return false;
      if(iangle>6) return false;
      elevation=elelist[iangle];
      azimuth=azilist[irot][itel][iangle];
      return fabs(azimuth)<500?true:false;
   }
   else return false;
}
double RotateDB::GetMinDistEleAzi(double ele_in,double azi_in,int irot,int itel,double &minele,double &minazi,int &index){
   index=-1;
   minele=1000;
   minazi=1000;
   double result=TMath::Max(minele,minazi);
   if(irot<0||irot>=2) return result;
   bool wrongtel=(itel<0||itel>=6);
   if(itel<0||itel>=6) return result;

   const int nangle1=4;
   double AngleList1[nangle1][2];
   AngleList1[0][0]=40; AngleList1[0][1]=42;
   AngleList1[1][0]=40; AngleList1[1][1]=29;
   AngleList1[2][0]=40; AngleList1[2][1]=19;
   AngleList1[3][0]=40; AngleList1[3][1]=4;

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

   const int nangle3=7;
   double elelist3=20.;
   double azilist3[nangle3]={19.,19.1,19.2,19.5,18.9,18.8,18.5};

   if(TargetIndex>=0){
      int TelTarget=(TargetIndex/1000);
      int iTelTarget=GetTi(TelTarget);
      int TargetIndex0=(TargetIndex%1000);
      if(irot!=GetLi(TargetIndex0/100)) return result;
      int itype=((TargetIndex0%100)/10);
      int iangle=(TargetIndex0%10);
      if(itype==1){
         if(wrongtel) return result;
         for(int ii=0;ii<nangle1;ii++){
            if(ii!=iangle) continue;
            bool isfine=false;
            if(ii==0&&(telindex[itel]==5||telindex[itel]==6)) isfine=true;
            if(ii==1&&(telindex[itel]==4||telindex[itel]==5)) isfine=true;
            if(ii==2&&(telindex[itel]==3||telindex[itel]==4)) isfine=true;
            if(ii==3&&(telindex[itel]==1||telindex[itel]==2||telindex[itel]==3)) isfine=true;
            if(!isfine) continue;
            if(TelTarget>0&&TelTarget!=telindex[itel]) continue;
            double dist1=fabs(ele_in-AngleList1[ii][0]);
            double dist2=fabs(azi_in-AngleList1[ii][1]);
            double dist=TMath::Max(dist1,dist2);
            if(dist<result){
               minele=dist1;
               minazi=dist2;
               result=dist;
               index=1*10+ii+telindex[itel]*1000;
            }
         }
      }
      else if(itype==2){
         for(int ii=0;ii<nangle2;ii++){
            if(ii!=iangle) continue;
            double dist1,dist2,dist;
            dist1=fabs(ele_in-elelist[ii]);
            int whichtel=0;
            if(TelTarget<=0||iTelTarget>=6){
               dist2=1000;
               for(int jj=0;jj<ntel;jj++){
                  double dist20=fabs(azi_in-azilist[irot][jj][ii]);
                  if(dist20<dist2||(dist20==dist2&&(jj==itel&&(!wrongtel)))){
                     dist2=dist20;
                     whichtel=telindex[jj];
                  }
               }
            }
            else{
               dist2=fabs(azi_in-azilist[irot][iTelTarget][ii]);
               whichtel=TelTarget;
            }
            dist=TMath::Max(dist1,dist2);
            if(dist<result){
               minele=dist1;
               minazi=dist2;
               result=dist;
               index=2*10+ii+whichtel*1000;
            }
         }
      }
      return result;
   }

   if(rotindex[irot]==2){
      for(int ii=0;ii<nangle1;ii++){
         bool isfine=false;
         if(!wrongtel){
         if(ii==0&&(telindex[itel]==5||telindex[itel]==6)) isfine=true;
         if(ii==1&&(telindex[itel]==4||telindex[itel]==5)) isfine=true;
         if(ii==2&&(telindex[itel]==3||telindex[itel]==4)) isfine=true;
         if(ii==3&&(telindex[itel]==1||telindex[itel]==2||telindex[itel]==3)) isfine=true;
         }
         if(!isfine) continue;
         double dist1=fabs(ele_in-AngleList1[ii][0]);
         double dist2=fabs(azi_in-AngleList1[ii][1]);
         double dist=TMath::Max(dist1,dist2);
         if(dist<result){
            minele=dist1;
            minazi=dist2;
            result=dist;
            index=1*10+ii+telindex[itel];
         }
      }
   }

   for(int ii=0;ii<nangle2;ii++){
      double dist1,dist2,dist;
      dist1=fabs(ele_in-elelist[ii]);
      int whichtel=0;
      if(SearchAllAngle||wrongtel){
         dist2=1000;
         for(int jj=0;jj<ntel;jj++){
            double dist20=fabs(azi_in-azilist[irot][jj][ii]);
            if(dist20<dist2||(dist20==dist2&&(jj==itel&&(!wrongtel)))){
               dist2=dist20;
               whichtel=telindex[jj];
            }
         }
      }
      else{
         dist2=fabs(azi_in-azilist[irot][itel][ii]);
         whichtel=telindex[itel];
      }
      dist=TMath::Max(dist1,dist2);
      if(dist<result){
         minele=dist1;
         minazi=dist2;
         result=dist;
         index=2*10+ii+whichtel*1000;
      }
   }

   //for(int ii=0;ii<nangle3&&(rotindex[irot]==2&&telindex[itel]==1);ii++){
   //   double dist1=fabs(ele_in-elelist3);
   //   double dist2=fabs(azi_in-azilist3[ii]);
   //   double dist=TMath::Max(dist1,dist2);
   //   if(dist<result){
   //      minele=dist1;
   //      minazi=dist2;
   //      result=dist;
   //      index=3*10+ii;
   //   }
   //}

   return result;
}
int RotateDB::IsFineAngle(double ele_in,double azi_in,int Li_in,int iTel){
   int irot=GetLi(Li_in);
   int itel=GetTi(iTel);
   int index;
   double mindistele,mindistazi;
   double mindist=GetMinDistEleAzi(ele_in,azi_in,irot,itel,mindistele,mindistazi,index);
   if(jdebug>3) printf("RotateDB::IsFineAngle: ele_in=%.2lf azi_in=%.2lf Li_in=%d iTel_in=%d, index=%d mindist_ele=%.2lf mindist_azi=%.2lf\n",ele_in,azi_in,Li_in,iTel,index,mindistele,mindistazi);
   if(mindist<aglmargin) return index;
   else return -1;
}
int RotateDB::SearchVersion(int time_in,int Li_in){
   int irot=GetLi(Li_in);
   if(irot<0) return -3;
   int time_new=GetLogTime(time_in,rotindex[irot]);
   int year=2000+(CommonTools::TimeFlag(time_new,1)%100);
   int month=CommonTools::TimeFlag(time_new,2);
   int day=CommonTools::TimeFlag(time_new,3);

   int ptime=time_new-24*3600;
   int ntime=time_new+24*3600;
   int pyear=2000+(CommonTools::TimeFlag(ptime,1)%100);
   int pmonth=CommonTools::TimeFlag(ptime,2);
   int pday=CommonTools::TimeFlag(ptime,3);
   int nyear=2000+(CommonTools::TimeFlag(ntime,1)%100);
   int nmonth=CommonTools::TimeFlag(ntime,2);
   int nday=CommonTools::TimeFlag(ntime,3);

   char name1[300]="";
   char name2[300]="";
   char pname2[300]="";
   char nname2[300]="";
   char* namebuff1=Form("/Laser%02d/L%d_SCdata/L%d/%d-%02d-%02d.txt.utf8",Li_in,Li_in,Li_in,year,month,day);
   char* namebuff2=Form("/Laser%02d/L%d_SCdata/%d/%02d/log_%02d%02d.txt",Li_in,Li_in,year,month,month,day);
   strcpy(name1,RotateLogDir);
   strcat(name1,namebuff1);
   strcpy(name2,RotateLogDir);
   strcat(name2,namebuff2);
   strcpy(pname2,RotateLogDir);
   strcat(pname2,Form("/Laser%02d/L%d_SCdata/%d/%02d/log_%02d%02d.txt",Li_in,Li_in,pyear,pmonth,pmonth,pday));
   strcpy(nname2,RotateLogDir);
   strcat(nname2,Form("/Laser%02d/L%d_SCdata/%d/%02d/log_%02d%02d.txt",Li_in,Li_in,nyear,nmonth,nmonth,nday));

   ifstream fin1;
   ifstream fin2;
   bool exist1=access(name1,4)>=0;
   bool exist2=access(name2,4)>=0;
   if(exist1&&exist2){
      fin1.open(name1,std::ios::in);
      char firstline1[300];
      char lastline1[300];
      bool isver1=GetFirstLastLine(&fin1,firstline1,lastline1);
      fin1.close();
      fin2.open(name2,std::ios::in);
      char firstline2[300];
      char lastline2[300];
      bool isver2=GetFirstLastLine(&fin2,firstline2,lastline2);
      fin2.close();
      int time1[2]={abs(ProcessTime(firstline1)),abs(ProcessTime(lastline1))};
      int time2[2]={abs(ProcessTime2(firstline2)),abs(ProcessTime2(lastline2))};
      if(time1[0]>1300000000&&isver1){
         if(time_new>=time1[0]&&time_new<=time1[1]) return 1;
      }
      if(time2[0]>1300000000&&isver2){
         if(time_new>=time2[0]&&time_new<=time2[1]) return 2;
         else if(time_new<time2[0]){
            ifstream fin3;
            fin3.open(pname2,std::ios::in);
            isver2=GetFirstLastLine(&fin3,firstline1,lastline1);
            if(fin3.is_open()) fin3.close();
            if(isver2){
               int time_last=abs(ProcessTime2(lastline1));
               if(time_new>=time_last) return 2;
            }
         }
         else{
            ifstream fin3;
            fin3.open(nname2,std::ios::in);
            isver2=GetFirstLastLine(&fin3,firstline1,lastline1);
            if(fin3.is_open()) fin3.close();
            if(isver2){
               int time_first=abs(ProcessTime2(firstline1));
               if(time_new<=time_first) return 2;
            }
         }
      }
      return -1;
   }
   else if(exist1){
      fin1.open(name1,std::ios::in);
      char firstline1[300];
      char lastline1[300];
      bool isver1=GetFirstLastLine(&fin1,(char*)firstline1,(char*)lastline1);
      fin1.close();
      int time1[2]={abs(ProcessTime(firstline1)),abs(ProcessTime(lastline1))};
      if(time1[0]>1300000000&&isver1){
         if(time_new>=time1[0]&&time_new<=time1[1]) return 1;
      }
      return -1;
   }
   else if(exist2){
      fin2.open(name2,std::ios::in);
      char firstline2[300];
      char lastline2[300];
      bool isver2=GetFirstLastLine(&fin2,firstline2,lastline2);
      fin2.close();
      int time2[2]={abs(ProcessTime2(firstline2)),abs(ProcessTime2(lastline2))};
      //printf("RotateDB::SearchVersion: time_new=%d time_first=%d time_last=%d\n",time_new,time2[0],time2[1]);
      if(time2[0]>1300000000&&isver2){
         char firstline1[300];
         char lastline1[300];
         if(time_new>=time2[0]&&time_new<=time2[1]) return 2;
         else if(time_new<time2[0]){
            ifstream fin3;
            fin3.open(pname2,std::ios::in);
            isver2=GetFirstLastLine(&fin3,firstline1,lastline1);
            if(fin3.is_open()) fin3.close();
            if(isver2){
               int time_last=abs(ProcessTime2(lastline1));
               //printf("RotateDB::SearchVersion: time_new=%d time_first_prefile=%d time_last=%d\n",time_new,time_last,time2[0]);
               if(time_new>=time_last) return 2;
            }
         }
         else{
            ifstream fin3;
            fin3.open(nname2,std::ios::in);
            isver2=GetFirstLastLine(&fin3,firstline1,lastline1);
            if(fin3.is_open()) fin3.close();
            if(isver2){
               int time_first=abs(ProcessTime2(firstline1));
               //printf("RotateDB::SearchVersion: time_new=%d time_first=%d time_last_nextfile=%d\n",time_new,time2[1],time_first);
               if(time_new<=time_first) return 2;
            }
         }
      }
      return -1;
   }
   else return -2; //both not exist
}
int RotateDB::GetEleAzi1(int time_in,int Li_in,int iTel){
   if(jdebug>0) printf("RotateDB::GetEleAzi1: timein=%d Lin=%d iTel=%d\n",time_in,Li_in,iTel);
   if(!LoadData(time_in,Li_in)) return -1;
   if(!ProcessAll(time_in,Li_in)) return -2;
   int irot=GetLi(Li_in);
   RotateDB rotbuff(this);

   double ele0=GetElevation();
   double azi0=GetAzimuth();
   if(jdebug>0) printf("RotateDB::GetEleAzi1: findtimelog time_in=%d ele=%lf azi=%lf\n",time_in,ele0,azi0);
   int retval=IsFineAngle(ele0,azi0,Li_in,iTel);
   if(retval<=0){
      if(jdebug>0) printf("RotateDB::GetEleAzi1: IsFineAngle=%d\n",retval);
      return -4;
   }

   int ncount1=0,ncount2=0;
   for(int itime=-1;itime>=-(60*100);itime--){
      int timei=time_in+itime;
      if(!LoadData(timei,Li_in)) break;
      if(!ProcessAll(timei,Li_in)) break;
      double elei=GetElevation();
      double azii=GetAzimuth();
      if(fabs(ele0-elei)>aglmargin) break;
      if(fabs(azi0-azii)>aglmargin) break;
      ncount1++;
   }
   fin_log1[irot].seekg((long int)rotbuff.fin_log1[irot].tellg(),std::ios::beg);
   for(int itime=1;itime<=(60*100);itime++){
      int timei=time_in+itime;
      if(!LoadData(timei,Li_in)) break;
      if(!ProcessAll(timei,Li_in)) break;
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

   if(ncount<ntotmin) return -(Li_in*100+retval);
   if(ncount1<=nsidemin||ncount2<=nsidemin) return -(Li_in*100+retval);
   else return (Li_in*100+retval);
}
int RotateDB::GetEleAzi2(int time_in,int Li_in,int iTel){
   if(jdebug>0) printf("RotateDB::GetEleAzi2: timein=%d Lin=%d iTel=%d\n",time_in,Li_in,iTel);
   int irot=GetLi(Li_in);
   if(irot<0) return -1;
   if(!LoadAData(time_in,Li_in,1)) return -2;
   if(UseGPSInfo){
      if(!LoadAData(time_in,Li_in,2)) return -2;
   }
   if(UseEnvInfo){
      if(!LoadAData(time_in,Li_in,3)) return -2;
   }
   if(UsePltInfo){
      if(!LoadAData(time_in,Li_in,0)) return -2;
   }
   if(!IsLogFine(1)) return -3;
   if(!ProcessAll2(time_in,Li_in)) return -3;

   double ele0=GetElevation();
   double azi0=GetAzimuth();
   if(jdebug>0) printf("RotateDB::GetEleAzi2: findtimelog ele=%lf azi=%lf\n",ele0,azi0);
   int retval=IsFineAngle(ele0,azi0,Li_in,iTel);
   if(retval<=0){
      if(jdebug>0) printf("RotateDB::GetEleAzi2: IsFineAngle=%d\n",retval);
      return -4;
   }

   int timei0=abs(ProcessTime2(1,0));
   int timei1=abs(ProcessTime2(1,1));
   int ncount1=GetLogTime(time_in,rotindex[irot])-timei0;
   int ncount2=timei1-GetLogTime(time_in,rotindex[irot]);
   int ncount=ncount1+ncount2+1;

   if(ncount<ntotmin) return -(Li_in*100+retval);
   if(ncount1<=nsidemin||ncount2<=nsidemin) return -(Li_in*100+retval);
   else return (Li_in*100+retval);
}
int RotateDB::GetEleAzi(int time_in,int Li_in,int iTel){
   int getversion=SearchVersion(time_in,Li_in);
   if(getversion==1) return GetEleAzi1(time_in,Li_in,iTel);
   else if(getversion==2) return GetEleAzi2(time_in,Li_in,iTel);
   else return -10;
}
int RotateDB::GetEleAzi(WFCTAEvent* pev){
   if(!pev) return -6;
   int time_in=pev->rabbitTime;
   int irot=GetLi((double)pev->rabbittime);
   if(irot<0) return -7;
   int Li_in=rotindex[irot];
   return GetEleAzi(time_in,Li_in,pev->iTel);
}
void RotateDB::GetMinDistFit(WFCTAEvent* pev,double ele_in,double azi_in,int Li_in,double &minphi,double &mincc,double &minphi_sigma,double &mincc_sigma){
   minphi=1000;
   mincc=1000;
   minphi_sigma=1000;
   mincc_sigma=1000;
   if(!pev) return;
   if(!pev->minimizer) return;
   double phi=pev->minimizer->X()[3]/PI*180;
   double ephi=pev->minimizer->Errors()[3]/PI*180;
   double cc=pev->minimizer->X()[2]/PI*180;
   double ecc=pev->minimizer->Errors()[2]/PI*180;

   int irot=Li_in>0?GetLi(Li_in):GetLi((double)pev->rabbittime);
   if(irot<0||irot>=2) return;
   int itel=GetTi(pev->iTel);
   if(itel<0||itel>=6) return;
   
   double cc_model,phi_model;
   double rotdir[2];
   TelGeoFit::CalDir_out(ele_in/180*PI,azi_in/180*PI,rotindex[irot],rotdir[0],rotdir[1]);
   rotdir[0]=PI/2-rotdir[0]; //from elevation to zenith
   double rotcoo[3];
   for(int ii=0;ii<3;ii++) rotcoo[ii]=TelGeoFit::GetRotPos(ii,rotindex[irot]);

   WFTelescopeArray* pa=WFTelescopeArray::GetHead();
   int itel0=pa->GetTelescope(pev->iTel);
   if(itel0<0) return;
   WFTelescope* pt=pa->pct[itel0];
   if(!pt) return;
   double zenith=pt->TelZ_;
   double azimuth=pt->TelA_;
   double telpos[3]={pt->Telx_,pt->Tely_,pt->Telz_};

   TelGeoFit::GetCCphi(zenith,azimuth,telpos,rotcoo,rotdir,cc_model,phi_model);
   cc_model*=180/PI;
   phi_model*=180/PI;

   double dist10=TMath::Min(fabs(phi-phi_model),fabs(phi+180-phi_model));
   dist10=TMath::Min(dist10,fabs(phi-180-phi_model));
   double dist20=fabs(cc-cc_model);
   if(dist10<minphi) minphi=dist10;
   if(dist20<mincc) mincc=dist20;
   minphi_sigma=minphi/ephi;
   mincc_sigma=mincc/ecc;
   return;

   //int TelTarget=EleAziIndex/1000;
   //int EleAziIndex0=(EleAziIndex%1000);
   //int itype=(EleAziIndex0%100)/10;
   //int index=EleAziIndex0%10;

   //if(itype==1&&rotindex[irot]==2){
   //   const int nangle1=4;
   //   if(index<0||index>=nangle1) return;
   //   double phi1[6][nangle1];
   //   double cc1[6][nangle1];
   //   for(int ii=0;ii<nangle1;ii++){
   //      for(int jj=0;jj<6;jj++){
   //         phi1[jj][ii]=-10000;
   //         cc1[jj][ii]=-10000;
   //         if(ii==0&&telindex[jj]==6){
   //            phi1[jj][ii]=134.93;
   //            cc1[jj][ii]=-5.55;
   //         }
   //         if(ii==0&&telindex[jj]==5){
   //            phi1[jj][ii]=110.89;
   //            cc1[jj][ii]=6.30;
   //         }
   //         if(ii==1&&telindex[jj]==5){
   //            phi1[jj][ii]=100.29;
   //            cc1[jj][ii]=-8.20;
   //         }
   //         if(ii==1&&telindex[jj]==4){
   //            phi1[jj][ii]=82.93;
   //            cc1[jj][ii]=5.92;
   //         }
   //         if(ii==2&&telindex[jj]==4){
   //            phi1[jj][ii]=76.82;
   //            cc1[jj][ii]=-5.72;
   //         }
   //         if(ii==2&&telindex[jj]==3){
   //            phi1[jj][ii]=54.42;
   //            cc1[jj][ii]=6.61;
   //         }
   //         if(ii==3&&telindex[jj]==3){
   //            phi1[jj][ii]=49.25;
   //            cc1[jj][ii]=-8.80;
   //         }
   //         if(ii==3&&telindex[jj]==2){
   //            phi1[jj][ii]=3.00;
   //            cc1[jj][ii]=3.03;
   //         }
   //         if(ii==3&&telindex[jj]==1){
   //            phi1[jj][ii]=26.04;
   //            cc1[jj][ii]=-1.72;
   //         }
   //      }
   //   }

   //   double phi0=180-phi1[itel][index];
   //   double cc0=-cc1[itel][index];
   //   double dist1=TMath::Min(fabs(phi-phi0),fabs(phi+180-phi0));
   //   dist1=TMath::Min(dist1,fabs(phi-180-phi0));
   //   double dist2=fabs(cc-cc0);
   //   if(dist1<minphi) minphi=dist1;
   //   if(dist2<mincc) mincc=dist2;
   //}
   //else if(itype==2){
   //   const int nangle2=7;
   //   if(index<0||index>=nangle2) return;
   //   double phi2[2][6][nangle2]={{{27.97,27.40,27.08,26.44,26.77,26.64,-10000},
   //                               {3.00,3.17,3.16,3.13,3.16,3.07,-10000},
   //                               {49.29,50.43,51.88,51.91,51.83,50.56,-10000},
   //                               {75.00,78.46,79.76,80.80,80.78,80.40,-10000},
   //                               {111.32,111.31,110.69,110.68,110.56,110.94,-10000},
   //                               {135.70,137.83,137.27,136.46,136.46,136.62,-10000}
   //                              },
   //                              {{134.08,134.34,133.91,134.05,133.69,133.19,-10000},
   //                               {106.63,105.91,104.02,105.04,104.61,-1000,-10000},
   //                               {156.15,157.01,157.54,158.09,157.82,157.85,-10000},
   //                               {0.00,0.09,0.15,179.82,0.10,179.63,-10000},
   //                               {24.68,24.19,24.13,23.95,23.50,22.61,-10000},
   //                               {48.75,48.35,48.50,48.62,49.40,48.29,-10000}
   //                              }
   //                             };
   //   double cc2[2][6][nangle2]={{{-1.99,0.11,2.81,4.13,3.34,3.15,-10000},
   //                               {3.88,0.90,2.01,0.70,-0.20,1.23,-10000},
   //                               {-2.45,0.94,2.00,2.04,2.37,0.17,-10000},
   //                               {-1.11,-0.65,1.30,2.34,2.12,1.82,-10000},
   //                               {0.95,2.31,1.66,1.99,2.23,1.86,-10000},
   //                               {-2.52,2.98,2.37,-0.09,-0.54,0.59,-10000}
   //                              },
   //                              {{-0.75,-0.14,-0.84,-0.46,-2.24,0.56,-10000},
   //                               {-1.10,0.98,-1.79,-0.01,-0.37,-1000,-10000},
   //                               {-0.99,-2.01,-0.69,0.77,-0.81,-0.50,-10000},
   //                               {-2.71,1.45,0.19,-2.36,1.25,-0.03,-10000},
   //                               {-1.01,1.37,1.88,0.99,-0.90,-3.19,-10000},
   //                               {1.16,-0.21,-1.11,-0.39,1.36,-1.30,-10000}
   //                              }
   //                             };
   //   double phi0=(180-phi2[irot][itel][index]);
   //   double cc0=-cc2[irot][itel][index];
   //   double dist1=TMath::Min(fabs(phi-phi0),fabs(phi+180-phi0));
   //   dist1=TMath::Min(dist1,fabs(phi-180-phi0));
   //   double dist2=fabs(cc-cc0);
   //   if(dist1<minphi) minphi=dist1;
   //   if(dist2<mincc) mincc=dist2;
   //}
}
bool RotateDB::IsFineImage(WFCTAEvent* pev,double ele_in,double azi_in){
   int irot=GetLi((double)pev->rabbittime);
   if(irot<0) return false; 
   double minphi,mincc,minphi_sigma,mincc_sigma;
   GetMinDistFit(pev,ele_in,azi_in,rotindex[irot],minphi,mincc,minphi_sigma,mincc_sigma);
   return (minphi<phimargin&&mincc<ccmargin)&&(minphi_sigma<phismargin&&mincc_sigma<ccsmargin);
}
bool RotateDB::IsFineImage(WFCTAEvent* pev,int EleAziIndex){
   double ele_in,azi_in;
   if(!GetEleAzi(EleAziIndex,ele_in,azi_in)) return false;
   return IsFineImage(pev,ele_in,azi_in);
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
      double npe2[2][6][nangle2];
      for(int ii=0;ii<2;ii++){
         for(int jj=0;jj<6;jj++){
            for(int kk=0;kk<nangle2;kk++){
               if(ii==0) npe2[ii][jj][kk]=1.31e6;
               else if(ii==1) npe2[ii][jj][kk]=1.0e6;
               else npe2[ii][jj][kk]=1.0e6;
            }
         }
      }
      npe2[0][0][0]=2.34e6;
      npe2[0][0][1]=1.31e6;
      npe2[0][0][2]=9.53e5;
      npe2[0][0][3]=7.74e5;
      npe2[0][0][4]=6.89e5;
      npe2[0][0][5]=6.1e5;
      npe2[1][4][0]=3.7e5;
      npe2[1][4][1]=3.1e4;
      npe2[1][4][2]=1.0e4;
      npe2[1][4][3]=1.0e4;
      npe2[1][4][4]=1.0e4;
      npe2[1][4][5]=1.0e4;
      npe2[1][5][0]=1.27e5;
      npe2[1][5][1]=2.7e5;
      npe2[1][5][2]=6.4e4;
      npe2[1][5][3]=1.0e4;
      npe2[1][5][4]=1.0e4;
      npe2[1][5][5]=1.0e4;
      if(iangle<0||iangle>=nangle2) return -1;
      else return npe2[irot][itel][iangle];
   }
   else return -1;
}
int RotateDB::LaserIsFine(WFCTAEvent* pev){
   int EleAziIndex=GetEleAzi(pev);
   if(jdebug>0) printf("RotateDB::LaserIsFine: EleAziIndex=%d\n",EleAziIndex);
   if(EleAziIndex<=0) return EleAziIndex;
   double ele_in=GetElevation();
   double azi_in=GetAzimuth();

   bool fitted=pev->DoFit(0,12);
   if(jdebug>1) printf("RotateDB::LaserIsFine: Fit=%d\n",fitted);
   bool image=IsFineImage(pev,ele_in,azi_in);
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
   char* namebuff1=Form("/Laser%02d/L%d_SCdata/L%d/%d-%02d-%02d.txt.utf8",Li_in,Li_in,Li_in,year,month,day);
   char* namebuff2=Form("/Laser%02d/L%d_SCdata/%d/%02d/log_%02d%02d.txt",Li_in,Li_in,year,month,month,day);
   strcpy(name1,RotateLogDir);
   strcat(name1,namebuff1);
   strcpy(name2,RotateLogDir);
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
      if(LoadData(time_in,Li_in)<0){
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
      ProcessAll(time_in,Li_in);
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
