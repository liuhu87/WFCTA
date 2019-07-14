#include "StatusDB.h"
#include "TString.h"
#include "common.h"

StatusDB* StatusDB::_Head=0;
char StatusDB::DBPath[200]="/eos/lhaaso/raw/wfcta";

uint8_t* StatusDB::buf=0;
unsigned long StatusDB::buflength=0;
char StatusDB::FILENAME[300]="";
FILE* StatusDB::fp=0;
vector<int> StatusDB::containtime;
int StatusDB::readtime=0;

vector<int> StatusDB::failtime;

WFCTADecode* StatusDB::wfctaDecode=0;

int StatusDB::jdebug=0;
int StatusDB::timemargin=15;

//status variables
int StatusDB::fpgaVersion[10];
long StatusDB::clb_initial_Time;
double StatusDB::clb_initial_time;
int StatusDB::fired_tube;
long StatusDB::status_readback_Time;
double StatusDB::status_readback_time;
short StatusDB::single_thresh[1024];
short StatusDB::record_thresh[1024];
long StatusDB::single_count[1024];
float StatusDB::DbTemp[1024];
long StatusDB::single_time[1024];
float StatusDB::HV[1024];
float StatusDB::PreTemp[1024];
float StatusDB::BigResistence[1024];
float StatusDB::SmallResistence[1024];
long StatusDB::ClbTime[1024];
float StatusDB::ClbTemp[1024];

void StatusDB::Init(){
   if(!buf){
      buf=new uint8_t[BUF_LEN];
      Reset();
   }
   if(!wfctaDecode) wfctaDecode=new WFCTADecode();
}
void StatusDB::Reset(){
   if(fp) fseek(fp,0,0);
   readtime=0;
   containtime.clear();
   failtime.clear();
   for(int ii=0;ii<BUF_LEN;ii++) buf[ii]=0;
   buflength=0;

   for(int ii=0;ii<10;ii++) fpgaVersion[ii]=0;
   clb_initial_Time=0;
   clb_initial_time=0;
   fired_tube=0;
   status_readback_Time=0;
   status_readback_time=0;
   for(int ii=0;ii<1024;ii++){
   single_thresh[ii]=-1000;
   record_thresh[ii]=-1000;
   single_count[ii]=-1000;
   DbTemp[ii]=-1000;
   single_time[ii]=-1000;
   HV[ii]=-1000;
   PreTemp[ii]=-1000;
   BigResistence[ii]=-1000;
   SmallResistence[ii]=-1000;
   ClbTime[ii]=-1000;
   ClbTemp[ii]=-1000;
   }
}
void StatusDB::Release(){
   if(fp) fclose(fp);
   if(buf) delete buf;
   if(wfctaDecode) delete wfctaDecode;
   _Head=0;
}
void StatusDB::SetDirectory(const char* dirname){
   strcpy(DBPath,dirname);
}
StatusDB* StatusDB::GetHead(){
   if(!_Head) _Head=new StatusDB();
   return _Head;
}
FILE* StatusDB::LocateFile(int iTel,int Time,double time){
   for(int ii=0;ii<containtime.size();ii++){
      if(Time==containtime.at(ii)){
         if(!fp) {printf("StatusDB::LocateFile: containtime error! time=%d\n",Time); return 0;}
         int itel=CommonTools::GetTelIndex(FILENAME,42,1);
         if(itel==iTel) return fp;
         else break;
      }
   }
   double time0=CommonTools::InvConvert(Time);
   if(time0<=0) return 0;
   int year=CommonTools::TimeFlag(time0,1)+2000;
   int month=CommonTools::TimeFlag(time0,2);
   int day=CommonTools::TimeFlag(time0,3);
   int hour=CommonTools::TimeFlag(time0,4);
   int minute=CommonTools::TimeFlag(time0,5);
   int second=CommonTools::TimeFlag(time0,6);
   char dirname[300]="";
   strcat(dirname,DBPath);
   strcat(dirname,Form("/%d",year));
   strcat(dirname,Form("/%02d%02d",month,day));
   string dirpath(dirname);
   vector<string> namebuff;
   CommonTools::getFiles(dirpath,namebuff);
   int itar=-1;
   int cfiletime=0;
   for(int ii=0;ii<namebuff.size();ii++){
      int itel=CommonTools::GetTelIndex(namebuff.at(ii).data(),42,1);
      if(itel!=iTel) continue;
      int filetime=CommonTools::GetTimeFromFileName(namebuff.at(ii).data(),46,12);
      //printf("index=%d %s %s time=%d\n",ii,dirpath.data(),namebuff.at(ii).data(),filetime);
      if(Time>=filetime){
         if(filetime>cfiletime){
            cfiletime=filetime;
            itar=ii;
         }
      }
   }
   if(abs(cfiletime-Time)<600){ //one file contains at most 1 hour data
      string cfilename(FILENAME);
      if(cfilename==namebuff.at(itar)){ //in the current file
         if(jdebug>0) printf("StatusDB::LocateFile: find current file: %d %s time=(%d,%d) Tel=%d\n",itar,namebuff.at(itar).data(),Time,cfiletime,iTel);
         containtime.push_back(Time);
         return fp;
      }
      else{ //in a new file
         if(jdebug>0) printf("StatusDB::LocateFile: find new file: %d %s time=(%d,%d) Tel=%d\n",itar,namebuff.at(itar).data(),Time,cfiletime,iTel);
         if(fp) fclose(fp);
         strcpy(FILENAME,namebuff.at(itar).data());
         fp=fopen(FILENAME,"rb");
         Reset();
         containtime.push_back(Time);
         return fp;
      }
   }
   else{ //there is no file contains this Time 
      if(jdebug>0) printf("StatusDB::LocateFile: Couldn't find file for time=%d Tel=%d\n",Time,iTel);
      return 0;
   }
}
bool StatusDB::LocateBlk(int Time,double time){
   if(!fp) return false;
   for(int ii=0;ii<failtime.size();ii++){
      if(jdebug>1) printf("StatusDB::LocateBlk: Time=%d in the failtime list\n",Time);
      if(Time==failtime.at(ii)) return false;
   }

   double readmargin=64;

   //backup the read pointer and buffer
   long int clocation=ftell(fp);
   uint8_t buffer[BUF_LEN];
   for(int ii=0;ii<buflength;ii++) buffer[ii]=buf[ii];

   //search the position to left  or right
   int64_t packSize;
   size_t size_of_read;
   if(abs(Time-readtime)<=timemargin){ //just read the data from buf
      if(jdebug>1) printf("StatusDB::LocateBlk: Time=%d in the current buffer(readtime=%d)\n",Time,readtime);
      return true;
   }
   else if(Time>readtime){ //search from right
      fseek(fp,buflength,1);
      while((size_of_read=fread((uint8_t *)buf,1,STATUS_BUF_LEN,fp))!=0){
         bool find=wfctaDecode->StatusPackCheck(buf,size_of_read,9)>=0;
         uint64_t Time0=0;
         double time0=0;
         if(find){
            packSize = wfctaDecode->PackSize();
            Time0 = wfctaDecode->GetStatusReadbackTime(buf,packSize);
            time0 = wfctaDecode->GetStatusReadbacktime(buf,packSize);
         }
         if(abs(Time-Time0)<=timemargin){ //find it
            if(jdebug>1) printf("StatusDB::LocateBlk: find Time=%d in the right(readtime=%d readtime_new=%d)\n",Time,readtime,Time0);
            readtime=Time0;
            buflength=(unsigned long)size_of_read;
            fseek(fp,-buflength,1);
            ResetBuffer();
            Fill();
            return true;
         }
         else if(Time<Time0){ //there is no coresponding time
            if(jdebug>1) printf("StatusDB::LocateBlk: couldn't find Time=%d in the right(readtime=%d readtime_new=%d)\n",Time,readtime,Time0);
            failtime.push_back(Time);
            fseek(fp,clocation,0);
            for(int ii=0;ii<buflength;ii++) buf[ii]=buffer[ii];
            return false;
         }
         fseek(fp,-readmargin,1);
      }
   }
   else{ //search from left
      if(clocation<=63){
         failtime.push_back(Time);
         return false;
      }
      bool dopre=fseek(fp,-STATUS_BUF_LEN,1)>=0;
      long int offset=-STATUS_BUF_LEN;
      if(!dopre){
         offset=-clocation;
         fseek(fp,offset,1);
      }
      while((size_of_read=fread((uint8_t *)buf,1,-offset,fp))!=0){
         bool find=wfctaDecode->StatusPackCheck(buf,-offset,9)>=0;
         uint64_t Time0=0;
         double time0=0;
         if(find){
            packSize = wfctaDecode->PackSize();
            Time0 = wfctaDecode->GetStatusReadbackTime(buf,packSize);
            time0 = wfctaDecode->GetStatusReadbacktime(buf,packSize);
         }
         if(abs(Time-Time0)<=timemargin){ //find it
            if(jdebug>1) printf("StatusDB::LocateBlk: find Time=%d in the left(readtime=%d readtime_new=%d)\n",Time,readtime,Time0);
            readtime=Time0;
            buflength=(long int)size_of_read;
            fseek(fp,-buflength,1);
            ResetBuffer();
            Fill();
            return true;
         }
         else if(Time>Time0&&Time0>0){  //there is no coresponding time
            if(jdebug>1) printf("StatusDB::LocateBlk: couldn't find Time=%d in the left(readtime=%d readtime_new=%d)\n",Time,readtime,Time0);
            failtime.push_back(Time);
            fseek(fp,clocation,0);
            for(int ii=0;ii<BUF_LEN;ii++) buf[ii]=buffer[ii];
            return false;
         }
         fseek(fp,-((unsigned long)size_of_read),1);
         fseek(fp,readmargin,1);
         bool dopre=fseek(fp,-STATUS_BUF_LEN,1)>=0;
         offset=-STATUS_BUF_LEN;
         if(!dopre){
            offset=-ftell(fp);
            fseek(fp,offset,1);
         }
      }
   }
}
void StatusDB::ResetBuffer(){
   int find9=wfctaDecode->StatusPackCheck(buf,buflength,9)+63;
   int find1=wfctaDecode->StatusPackCheck(buf,buflength,1);
   if(find1>=0){ //there are all the status packets
      buflength=(find9-find1+1);
      fseek(fp,find1,1);
      size_t size_of_read=fread((uint8_t *)buf,1,buflength,fp);
      fseek(fp,-((long int)size_of_read),1);
      if(jdebug>2) printf("StatusDB::ResetBuffer: there is F1 and F9 in the buffer\n");
   }
   else{ //no F1 status packets
      //backup the read pointer and buffer
      long int clocation=ftell(fp);
      uint8_t buffer[BUF_LEN];
      for(int ii=0;ii<=find9;ii++) buffer[ii]=buf[ii];
      //search STATUS_BUF_LEN until F1 status packet
      fseek(fp,find9+1,1);
      bool dopre=fseek(fp,-STATUS_BUF_LEN,1)>=0;
      long int offset=-STATUS_BUF_LEN;
      if(!dopre){
         offset=-ftell(fp);
         fseek(fp,offset,1);
      }
      size_t size_of_read=fread((uint8_t *)buf,1,-offset,fp);
      find9=wfctaDecode->StatusPackCheck(buf,-offset,9)+63; //F9 packet should exist
      if(find9+1!=-offset) printf("StatusDB::ResetBuffer: Error, there is no F9 packet in the buffer\n");
      find1=wfctaDecode->StatusPackCheck(buf,-offset,1);
      if(find1>=0){
         buflength=(find9-find1+1);
         fseek(fp,offset+find1,1);
         if(jdebug>2) printf("StatusDB::ResetBuffer: Find F1 packet by search to left %ld bytes\n",-offset-find1);
      }
      else{
         //keep remaining packet
         buflength=-offset;
         fseek(fp,offset,1);
         if(jdebug>2) printf("StatusDB::ResetBuffer: No F1 packet after searching to left %ld bytes, keep only the remaining packet\n",-offset);
      }
      size_of_read=fread((uint8_t *)buf,1,buflength,fp);
      fseek(fp,-((long int)size_of_read),1);
   }
}
void StatusDB::Fill(){
   if(buflength<64) return;
   long int offset=0;
   long int packSize;
   int status_pack_marker;
   int nfilled=0;
   while(true){
      status_pack_marker=wfctaDecode->StatusPackCheck(buf+offset,buflength-offset);
      if(status_pack_marker>0){
         packSize = wfctaDecode->PackSize();
         if(status_pack_marker>=1&&status_pack_marker<=8){
            if(jdebug>5) printf("StatusDB::Fill:packetpos=%5ld, Fill F%d packet\n",offset+packSize,status_pack_marker);
            nfilled++;
         }
         switch(status_pack_marker){
           case 0x21:
             wfctaDecode->Getthresh(buf,offset+packSize,(short *)single_thresh, (short *)record_thresh);
             if(jdebug>4) printf("StatusDB::Fill:packetpos=%5ld, Fill threshold\n",offset+packSize);
             nfilled++;
             break;
           case 0x22:
             wfctaDecode->Deal22Pack(buf,offset+packSize,(long *)single_count);
             if(jdebug>4) printf("StatusDB::Fill:packetpos=%5ld, Fill single count\n",offset+packSize);
             nfilled++;
             break;
           case 0x23:
             wfctaDecode->Deal23Pack(buf,offset+packSize,(long *)single_count,(long *)single_time);
             if(jdebug>4) printf("StatusDB::Fill:packetpos=%5ld, Fill single count and time\n",offset+packSize);
             nfilled++;
             break;
           case 0x81:
             wfctaDecode->GetHV(buf,offset+packSize,(float *)HV);
             if(jdebug>4) printf("StatusDB::Fill:packetpos=%5ld, Fill HV\n",offset+packSize);
             nfilled++;
             break;
           case 0x82:
             wfctaDecode->GetPreTemp(buf,offset+packSize,(float *)PreTemp);
             if(jdebug>4) printf("StatusDB::Fill:packetpos=%5ld, Fill PreTemp\n",offset+packSize);
             nfilled++;
             break;
           case 0x83:
             //wfctaDecode->Deal83Package((float *)BigResistence);
             if(jdebug>5) printf("StatusDB::Fill:packetpos=%5ld, Fill 0x83 packet\n",offset+packSize);
             nfilled++;
             break;
           case 0x84:
             //wfctaDecode->Deal84Package((float *)SmallResistence);
             if(jdebug>5) printf("StatusDB::Fill:packetpos=%5ld, Fill 0x84 packet\n",offset+packSize);
             nfilled++;
             break;
           case 0x85:
             wfctaDecode->GetClbTemp(buf,offset+packSize,(float *)ClbTemp);
             if(jdebug>5) printf("StatusDB::Fill:packetpos=%5ld, Fill ClbTemp\n",offset+packSize);
             nfilled++;
             break;
           case 0x9:
             clb_initial_Time = wfctaDecode->GetclbInitialTime(buf,offset+packSize);
             clb_initial_time = wfctaDecode->GetclbInitialtime(buf,offset+packSize);
             fired_tube = wfctaDecode->GetFiredTube(buf,offset+packSize);
             status_readback_Time = wfctaDecode->GetStatusReadbackTime(buf,offset+packSize);
             status_readback_time = wfctaDecode->GetStatusReadbacktime(buf,offset+packSize);
             if(jdebug>4) printf("StatusDB::Fill:packetpos=%5ld, Fill F9 packet\n",offset+packSize);
             nfilled++;
             break;
         }
         offset+=packSize;
      }
      else break;
   }
   if(jdebug>3) printf("StatusDB::Fill: Fill %d packets for time=%d at %ld of file %s\n",nfilled,readtime,ftell(fp),FILENAME);
}
/*void StatusDB::Fill(){
   if(buflength<64) return;
   //const int nindex=17;
   //int index[nindex]={1,2,3,4,5,6,7,8,9,0x21,0x22,0x23,0x81,0x82,0x83,0x84,0x85};
   const int nindex=9;
   int index[nindex]={1,2,3,4,5,6,7,8,9};
   const int nindex2=8;
   int index2[nindex2]={0x21,0x22,0x23,0x81,0x82,0x83,0x84,0x85};
   for(int ii=0;ii<nindex;ii++){
      int exist=wfctaDecode->StatusPackCheck(buf,buflength,index[ii]);
      if(exist<0) continue;
      int64_t packSize0 = wfctaDecode->PackSize();
      if(jdebug>4) printf("StatusDB::Fill: F%d packsize0=%ld buflength=%ld\n",index[ii],packSize0,buflength);
      for(int i2=0;i2<(index[ii]<9?8:1);i2++){
         for(int i3=0;i3<((index[ii]==9)?1:(nindex2+1));i3++){
            int code=index[ii];
            int packSize=packSize0;
            if(i2>0){
               if(wfctaDecode->StatusPackCheck(buf+packSize0,buflength-(packSize0),index2[i2-1])<0) continue;
               code=index2[i2-1];
               packSize = wfctaDecode->PackSize()+packSize0;
               if(jdebug>4) printf("StatusDB::Fill: F%d 0x%x packsize=%ld\n",index[ii],index2[i2-1],packSize);
            }
            int imax;
            float max;
            switch(index[ii]){
               case 0x21:
                 wfctaDecode->Getthresh(buf,packSize,(short *)single_thresh, (short *)record_thresh);
                 if(jdebug>3) printf("StatusDB::Fill: Fill threshold\n");
                 continue;
               case 0x22:
                 wfctaDecode->Deal22Pack(buf,packSize,(long *)single_count);
                 if(jdebug>3) printf("StatusDB::Fill: Fill single count\n");
                 continue;
               case 0x23:
                 wfctaDecode->Deal23Pack(buf,packSize,(long *)single_count,(long *)single_time);
                 if(jdebug>3) printf("StatusDB::Fill: Fill single count and time\n");
                 continue;
               case 0x81:
                 wfctaDecode->GetHV(buf,packSize,(float *)HV);
                 imax=-1;
                 max=-1000;
                 for(int jj=0;jj<1024;jj++) if(HV[jj]>max) {imax=jj; max=HV[jj];}
                 if(jdebug>3) printf("StatusDB::Fill: Fill HV max=%f(%d)\n",max,imax);
                 continue;
               case 0x82:
                 wfctaDecode->GetPreTemp(buf,packSize,(float *)PreTemp);
                 imax=-1;
                 max=-1000;
                 for(int jj=0;jj<1024;jj++) if(PreTemp[jj]>max) {imax=jj; max=PreTemp[jj];}
                 if(jdebug>3) printf("StatusDB::Fill: Fill PreTemp max=%f(%d)\n",max,imax);
                 continue;
               //case 0x83:
               //  wfctaDecode->Deal83Package((float *)BigResistence);
               //  break;
               //case 0x84:
               //  wfctaDecode->Deal84Package((float *)SmallResistence);
               //  break;
               case 0x85:
                 wfctaDecode->GetClbTemp(buf,packSize,(float *)ClbTemp);
                 imax=-1;
                 max=-1000;
                 for(int jj=0;jj<1024;jj++) if(ClbTemp[jj]>max) {imax=jj; max=ClbTemp[jj];}
                 if(jdebug>3) printf("StatusDB::Fill: Fill ClbTemp max=%f(%d)\n",max,imax);
                 continue;
               case 0x9:
                 if(jdebug>3) printf("StatusDB::Fill: Fill F9 packet\n");
                 clb_initial_Time = wfctaDecode->GetclbInitialTime(buf,packSize);
                 clb_initial_time = wfctaDecode->GetclbInitialtime(buf,packSize);
                 fired_tube = wfctaDecode->GetFiredTube(buf,packSize);
                 status_readback_Time = wfctaDecode->GetStatusReadbackTime(buf,packSize);
                 status_readback_time = wfctaDecode->GetStatusReadbacktime(buf,packSize);
                 continue;
            }
         }
      }
   }
}*/

bool StatusDB::Locate(int iTel,int Time,double time){
   if(!LocateFile(iTel,Time,time)) return false;
   return LocateBlk(Time,time);
}

int StatusDB::GetVersion(int iTel,int Time,int i){
   if(!Locate(iTel,Time)) return -1000;
   return fpgaVersion[i];
}
long StatusDB::GetClbInitTime(int iTel,int Time){
   if(!Locate(iTel,Time)) return 0;
   return clb_initial_Time;
}
double StatusDB::GetClbInittime(int iTel,int Time){
   if(!Locate(iTel,Time)) return 0;
   return clb_initial_time;
}
int StatusDB::GetFiredTube(int iTel,int Time){
   if(!Locate(iTel,Time)) return 0;
   return fired_tube;
}
long StatusDB::GetReadbackTime(int iTel,int Time){
   if(!Locate(iTel,Time)) return 0;
   return status_readback_Time;
}
double StatusDB::GetReadbacktime(int iTel,int Time){
   if(!Locate(iTel,Time)) return 0;
   return status_readback_time;
}
short StatusDB::GetSingleThrd(int iTel,int Time,int i){
   if(!Locate(iTel,Time)) return -1000;
   return single_thresh[i];
}
short StatusDB::GetRecordThrd(int iTel,int Time,int i){
   if(!Locate(iTel,Time)) return -1000;
   return record_thresh[i];
}
long StatusDB::GetSingleCount(int iTel,int Time,int i){
   if(!Locate(iTel,Time)) return -1000;
   return single_count[i];
}
float StatusDB::GetDbTemp(int iTel,int Time,int i){
   if(!Locate(iTel,Time)) return -1000;
   return single_count[i];
}
long StatusDB::GetSingleTime(int iTel,int Time,int i){
   if(!Locate(iTel,Time)) return -1000;
   return single_time[i];
}
float StatusDB::GetHV(int iTel,int Time,int i){
   if(!Locate(iTel,Time)) return -1000;
   return HV[i];
}
float StatusDB::GetPreTemp(int iTel,int Time,int i){
   if(!Locate(iTel,Time)) return -1000;
   return PreTemp[i];
}
float StatusDB::GetBigR(int iTel,int Time,int i){
   if(!Locate(iTel,Time)) return -1000;
   return BigResistence[i];
}
float StatusDB::GetSmallR(int iTel,int Time,int i){
   if(!Locate(iTel,Time)) return -1000;
   return SmallResistence[i];
}
long StatusDB::GetClbTime(int iTel,int Time,int i){
   if(!Locate(iTel,Time)) return -1000;
   return ClbTime[i];
}
float StatusDB::GetClbTemp(int iTel,int Time,int i){
   if(!Locate(iTel,Time)) return -1000;
   return ClbTemp[i];
}
