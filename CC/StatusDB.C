#include "StatusDB.h"
#include "TString.h"
#include "common.h"

StatusDB* StatusDB::_Head=0;
char StatusDB::DBPath[200]="/eos/lhaaso/raw/wfcta";

uint8_t* StatusDB::buf=0;
int StatusDB::buflength=0;
char StatusDB::FILENAME[300]="";
FILE* StatusDB::fp=0;
int StatusDB::readtime=0;
vector<int> StatusDB::failtime;

WFCTADecode* StatusDB::wfctaDecode=0;

int StatusDB::timemargin=15;
void StatusDB::Init(){
   if(!buf){
      buf=new uint8_t[BUF_LEN];
      Reset();
   }
   if(!wfctaDecode) wfctaDecode=new WFCTADecode();
}
void StatusDB::Reset(){
   readtime=0;
   failtime.clear();
   for(int ii=0;ii<BUF_LEN;ii++) buf[ii]=0;
   buflength=0;
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
FILE* StatusDB::LocateFile(int Time,double time){
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
      printf("find file: %d %s time=%d\n",itar,namebuff.at(itar).data(),cfiletime);
      string cfilename(FILENAME);
      if(cfilename==namebuff.at(itar)){ //in the current file
         return fp;
      }
      else{ //in a new file
         if(fp) fclose(fp);
         strcpy(FILENAME,namebuff.at(itar).data());
         fp=fopen(FILENAME,"rb");
         fseek(fp,0,0);
         Reset();
         return fp;
      }
   }
   else return 0; //there is no file contains this Time
}
bool StatusDB::LocateBlk(int Time,double time){
   if(!fp) return false;
   for(int ii=0;ii<failtime.size();ii++){
      if(Time==failtime.at(ii)) return false;
   }

   //backup the read pointer and buffer
   long int clocation=ftell(fp);
   uint8_t buffer[BUF_LEN];
   for(int ii=0;ii<buflength;ii++) buffer[ii]=buf[ii];

   //
   int64_t packSize;
   size_t size_of_read;
   if(abs(Time-readtime)<=timemargin){ //just read the data from buf
      return true;
   }
   else if(Time>readtime){ //search from right
      fseek(fp,buflength,1);
      while((size_of_read=fread((uint8_t *)buf,1,STATUS_BUF_LEN,fp))!=0){
         bool find=wfctaDecode->StatusPackCheck(buf,size_of_read,9);
         uint64_t Time0=0;
         double time0=0;
         if(find){
            packSize = wfctaDecode->PackSize();
            Time0 = wfctaDecode->GetStatusReadbackTime(buf,packSize);
            time0 = wfctaDecode->GetStatusReadbacktime(buf,packSize);
         }
         if(abs(Time-Time0)<=timemargin){ //find it
            readtime=Time0;
            buflength=size_of_read;
            fseek(fp,-buflength,1);
            return true;
         }
         else if(Time<Time0){ //there is no coresponding time
            failtime.push_back(Time);
            fseek(fp,clocation,0);
            for(int ii=0;ii<buflength;ii++) buf[ii]=buffer[ii];
            return false;
         }
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
         bool find=wfctaDecode->StatusPackCheck(buf,-offset,9);
         uint64_t Time0=0;
         double time0=0;
         if(find){
            packSize = wfctaDecode->PackSize();
            Time0 = wfctaDecode->GetStatusReadbackTime(buf,packSize);
            time0 = wfctaDecode->GetStatusReadbacktime(buf,packSize);
         }
         if(abs(Time-Time0)<=timemargin){ //find it
            readtime=Time0;
            buflength=size_of_read;
            fseek(fp,-buflength,1);
            return true;
         }
         else if(Time>Time0){
            failtime.push_back(Time);
            fseek(fp,clocation,0);
            for(int ii=0;ii<BUF_LEN;ii++) buf[ii]=buffer[ii];
            return false;
         }
      }
   }
}
/*bool StatusDB::Fill(int Time,double time){
   if(!LocateFile(Time,time)) return false;
   bool getblk=LocateBlk(Time,time);
   if(!getblk) return false;

   const int nindex=17;
   int index[nindex]={1,2,3,4,5,6,7,8,9,0x21,0x22,0x23,0x81,0x82,0x83,0x84,0x85};
   for(int ii=0;ii<nindex;ii++){
      bool exist=wfctaDecode->StatusPackCheck(buf,buflength,index[ii]);
      if(!exist) continue;
      switch(index[ii]){
         case 0x21:
           wfctaDecode->Getthresh(buf,packSize,(short *)single_thresh, (short *)record_thresh);
           break;
         case 0x22:
           wfctaDecode->Deal22Pack(buf,packSize,(long *)single_count);
           break;
         case 0x23:
           wfctaDecode->Deal23Pack(buf,packSize,(long *)single_count,(long *)single_time);
           break;
         case 0x81:
           wfctaDecode->GetHV(buf,packSize,(float *)HV);
           break;
         case 0x82:
           wfctaDecode->GetPreTemp(buf,packSize,(float *)PreTemp);
           break;
         //case 0x83:
         //  wfctaDecode->Deal83Package((float *)BigResistence);
         //  break;
         //case 0x84:
         //  wfctaDecode->Deal84Package((float *)SmallResistence);
         //  break;
         case 0x85:
           wfctaDecode->GetClbTemp(buf,packSize,(float *)ClbTemp);
           break;
         case 0x9:
           clb_initial_Time = wfctaDecode->GetclbInitialTime(buf,packSize);
           clb_initial_time = wfctaDecode->GetclbInitialtime(buf,packSize);
           fired_tube = wfctaDecode->GetFiredTube(buf,packSize);
           status_readback_Time = wfctaDecode->GetStatusReadbackTime(buf,packSize);
           status_readback_time = wfctaDecode->GetStatusReadbacktime(buf,packSize);
      }
   }
}*/
