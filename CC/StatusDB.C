#include "StatusDB.h"
#include "TString.h"
#include "common.h"
#include <iostream>

StatusDB* StatusDB::_Head=0;
char StatusDB::DBPath0[2][200]={"/eos/lhaaso/raw/wfcta","/eos/lhaaso/decode/wfcta"};
char StatusDB::DBPath[200];

//char StatusDB::FILENAME[300]="";
//FILE* StatusDB::fp=0;
//uint8_t* StatusDB::buf=0;
//unsigned long StatusDB::buflength=0;
//
//TFile* StatusDB::fstatus=0;
//TTree* StatusDB::tree=0;
//
//int StatusDB::currentday=0;
//int StatusDB::nfiles[NCTMax];
//int StatusDB::telindex[NCTMax];
//int StatusDB::filetime[NCTMax][1000];
//int StatusDB::fileindex[NCTMax][1000];
//char StatusDB::filename[NCTMax][1000][200];
//
//int StatusDB::currentfile=-1;
//int StatusDB::currenttime=0;
//long int StatusDB::entryno=-1;
//vector<int> StatusDB::failday;
//vector<int> StatusDB::failtime;
//
//WFCTADecode* StatusDB::wfctaDecode=0;

bool StatusDB::UseDecodeData=true;
int StatusDB::jdebug=0;
int StatusDB::timemargin=15;
int StatusDB::MaxFileTime=3600;

////status variables
//Short_t StatusDB::iTel;
//int StatusDB::fpgaVersion[10];
//int StatusDB::f9mode;
//int StatusDB::f9pattern;
//int StatusDB::DbVersion[2][89];
//int StatusDB::ClbVersion[2][89];
//Long64_t StatusDB::clb_initial_Time;
//double StatusDB::clb_initial_time;
//int StatusDB::fired_tube;
//Long64_t StatusDB::status_readback_Time;
//double StatusDB::status_readback_time;
//int StatusDB::sipm[MAXPMT];
//int StatusDB::mask[MAXPMT];
//short StatusDB::single_thresh[MAXPMT];
//short StatusDB::record_thresh[MAXPMT];
//Long64_t StatusDB::single_count[MAXPMT];
//Long64_t StatusDB::single_time[MAXPMT];
//float StatusDB::DbTemp[MAXPMT];
//float StatusDB::HV[MAXPMT];
//float StatusDB::PreTemp[MAXPMT];
//float StatusDB::BigResistence[MAXPMT];
//float StatusDB::SmallResistence[MAXPMT];
//Long64_t StatusDB::ClbTime[MAXPMT];
//float StatusDB::ClbTemp[MAXPMT];

void StatusDB::Init(){
   if(UseDecodeData){
      fstatus=0;
      tree=0;
   }
   else{
      fp=0;
      buf=new uint8_t[BUF_LEN];
      wfctaDecode=new WFCTADecode();
   }
   Reset();
}
void StatusDB::Reset(){
   if(UseDecodeData){
      entryno=-1;
   }
   else{
      strcpy(FILENAME,"");
      if(fp) fseek(fp,0,0);
      if(buf) {for(int ii=0;ii<BUF_LEN;ii++) buf[ii]=0;}
      buflength=0;
   }

   if(tree) {delete tree; tree=0;}
   failday.clear();
   failtime.clear();

   for(int ii=0;ii<10;ii++) fpgaVersion[ii]=0;
   clb_initial_Time=0;
   clb_initial_time=0;
   fired_tube=0;
   status_readback_Time=0;
   status_readback_time=0;
   for(int ii=0;ii<MAXPMT;ii++){
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
   if(tree) delete tree;
   if(fstatus){fstatus->Close();}
   failday.clear();
   failtime.clear();
   _Head=0;
}
void StatusDB::SetDirectory(const char* dirname){
   if(dirname) strcpy(DBPath,dirname);
   else{
      strcpy(DBPath,DBPath0[UseDecodeData?1:0]);
   }
}
StatusDB* StatusDB::GetHead(){
   if(!_Head) _Head=new StatusDB();
   return _Head;
}
int StatusDB::LocateTel(int iTel0){
   int rectel=-1;
   for(int ii=0;ii<NCTMax;ii++){
      if(telindex[ii]==100||telindex[ii]==iTel0){
         rectel=ii;
         break;
      }
   }
   return rectel;
}
int StatusDB::LoadFile(int whichday,char* filename_in){
   int type=0;
   if(whichday>=20190101&&whichday<=20500000) type=1;
   if(whichday>=1300000000&&whichday<=2524579200) type=2;
   if(jdebug>2) printf("StatusDB::LoadFile: day=%d type=%d filename=%s\n",whichday,type,filename_in);
   if(type<1||type>2) return 0;
   if(whichday==currentday){
      if(type==1){
         int nsum=0;
         for(int ii=0;ii<NCTMax;ii) nsum+=nfiles[ii];
         return nsum;
      }
      if(type==2) return currentfile>=-1?1:0;
   }
   if(type==2){
      if(!filename_in) return 0;
      if(UseDecodeData){
         TFile* fbuff=TFile::Open(filename_in);
         if(!fbuff) return 0;
         TTree* tree0=(TTree*)fbuff->Get("Status");
         if(!tree0) {fbuff->Close(); return 0;}
         if(tree0->GetEntries()<=0) {fbuff->Close(); return 0;}
         fbuff->Close();
         currentday=whichday;
         return 1;
      }
      else{
         char filetype[30];
         CommonTools::GetFileType(filetype,filename_in);
         if(!strstr(filetype,"dat")) return 0;
         FILE* fp0=fopen(filename_in,"rb");
         if(fp0==NULL) return 0;
         fclose(fp0);
         currentday=whichday;
         return 1;
      }
   }
   int year=whichday/10000;
   int month=(whichday%10000)/100;
   int day=(whichday%100);
   char dirname[300]="";
   strcat(dirname,DBPath);
   strcat(dirname,Form("/%d",year));
   strcat(dirname,Form("/%02d%02d",month,day));
   string dirpath(dirname);

   int telindex0[NCTMax];
   for(int ii=0;ii<NCTMax;ii++){
      telindex0[ii]=telindex[ii];
      telindex[ii]=0;
   }
   int nfiles0[NCTMax];
   for(int ii=0;ii<NCTMax;ii++) nfiles0[ii]=0;
   //list all the status files
   vector<string> namebuff;
   CommonTools::getFiles(dirpath,namebuff);
   if(jdebug>2) printf("StatusDB::LoadFile: whichday=%d dirname=%s nfile=%d\n",whichday,dirname,namebuff.size());
   int nloaded=0;
   for(int ii=0;ii<namebuff.size();ii++){
      int cfiletime=CommonTools::GetTimeFromFileName(namebuff.at(ii).data());
      if(cfiletime<1300000000||cfiletime>2000000000) continue;
      char filetype[30];
      CommonTools::GetFileType(filetype,namebuff.at(ii).data());
      int itel=CommonTools::GetTelIndex(namebuff.at(ii).data());
      if(itel<=0) continue;
      else if(itel<=100){ //status data in status file
         if(UseDecodeData){
            if(!(strstr(namebuff.at(ii).data(),".status.root"))) continue;
            TFile* fbuff=TFile::Open(Form("root://eos01.ihep.ac.cn/%s",namebuff.at(ii).data()));
            if(!fbuff) continue;
            TTree* tree0=(TTree*)fbuff->Get("Status");
            if(!tree0) {fbuff->Close(); continue;}
            if(tree0->GetEntries()<=0) {fbuff->Close(); continue;}
            fbuff->Close();
         }
         else{
            if(!strstr(filetype,"dat")) continue;
            FILE* fp0=fopen(filename_in,"rb");
            if(fp0==NULL) continue;
            fclose(fp0);
         }
      }
      else{ //status data in one single tree inside the same file with science data
         itel=(itel%100);
         if(UseDecodeData){
            TFile* fbuff=TFile::Open(Form("root://eos01.ihep.ac.cn/%s",namebuff.at(ii).data()));
            if(!fbuff) continue;
            TTree* tree0=(TTree*)fbuff->Get("Status");
            if(!tree0) {fbuff->Close(); continue;}
            if(tree0->GetEntries()<=0) {fbuff->Close(); continue;}
            fbuff->Close();
         }
         else{
            if(!strstr(filetype,"dat")) continue;
            FILE* fp0=fopen(filename_in,"rb");
            if(fp0==NULL) continue;
            fclose(fp0);
         }
      }
      int rectel=(itel==100)?0:(itel-1);
      if(rectel<0||rectel>NCTMax) continue;
      //printf("index=%d %s %s time=%d\n",ii,dirpath.data(),namebuff.at(ii).data(),cfiletime);
      if(jdebug>2) printf("StatusDB::LoadFile: whichday=%d ii=%d(nfiles=%d) filename=%s\n",whichday,ii,nfiles0[rectel],namebuff.at(ii).data());

      //record all the status files
      if(nfiles0[rectel]>=1000){
         std::cerr<<"StatusDB::LoadFile: the max number of files is too small(n="<<namebuff.size()<<",max nfiles=1000)..."<<endl;
         break;
      }
      if(telindex[rectel]<=0) telindex[rectel]=itel;
      strcpy(filename[rectel][nfiles0[rectel]],namebuff.at(ii).data());
      //reorder according the file time
      int timeindex=-1;
      for(int ii=0;ii<nfiles0[rectel];ii++){
         int itime=filetime[rectel][ii];
         int itime2=filetime[rectel][ii+1];
         if(cfiletime>=itime&&cfiletime<itime2){
            timeindex=ii;
            break;
         }
      }
      if(jdebug>4){
         printf("StatusDB::LoadFile: cfiletime=%d itel=%d index=%d\n",cfiletime,rectel,timeindex);
         for(int ii=0;ii<nfiles0[rectel];ii++) printf("ii=%d time=%d\n",ii,filetime[rectel][ii]);
      }
      if(timeindex<0){
         if(jdebug>5) for(int ii=0;ii<nfiles0[rectel];ii++) printf("i1=%d time=%d\n",ii,filetime[rectel][ii]);
         if(nfiles0[rectel]-1>=0&&cfiletime>filetime[rectel][nfiles0[rectel]-1]){
            filetime[rectel][nfiles0[rectel]]=cfiletime;
            fileindex[rectel][nfiles0[rectel]]=nfiles0[rectel];
         }
         else{
            for(int ii=nfiles0[rectel];ii>=1;ii--){
               filetime[rectel][ii]=filetime[rectel][ii-1];
               fileindex[rectel][ii]=fileindex[rectel][ii-1];
            }
            filetime[rectel][0]=cfiletime;
            fileindex[rectel][0]=nfiles0[rectel];
         }
         if(jdebug>5) for(int ii=0;ii<nfiles0[rectel];ii++) printf("i2=%d time=%d\n",ii,filetime[rectel][ii]);
      }
      else{
         for(int ii=nfiles0[rectel];ii>=timeindex+2;ii--){
            filetime[rectel][ii]=filetime[rectel][ii-1];
            fileindex[rectel][ii]=fileindex[rectel][ii-1];
         }
         filetime[rectel][timeindex+1]=cfiletime;
         fileindex[rectel][timeindex+1]=nfiles0[rectel];
      }
      if(jdebug>4) for(int ii=0;ii<nfiles0[rectel];ii++) printf("i3=%d time=%d\n",ii,filetime[rectel][ii]);
      nfiles0[rectel]++;
      nloaded++;
   }
   if(nloaded>0){
      for(int ii=0;ii<NCTMax;ii++) nfiles[ii]=nfiles0[ii];
      currentday=year*10000+month*100+day;
      if(jdebug>3){
         for(int itel=0;itel<NCTMax;itel++){
            for(int ii=0;ii<nfiles[itel];ii++){
               printf("itel=%d ifile=%d(nfile=%d) whichday=%d time=%d filename=%s\n",itel,ii,nfiles[itel],whichday,filetime[itel][ii],filename[itel][ii]);
            }
         }
      }
   }
   else{
      for(int ii=0;ii<NCTMax;ii++) telindex[ii]=telindex0[ii];
   }
   return nloaded;
}
int StatusDB::LocateFile(int iTel0,int Time,char* filename_in,double time){
   if(Time<1300000000||Time>2000000000) return -1;
   int year=(CommonTools::TimeFlag(Time,1)%100)+2000;
   int month=CommonTools::TimeFlag(Time,2);
   int day=CommonTools::TimeFlag(Time,3);
   int hour=CommonTools::TimeFlag(Time,4);
   int minute=CommonTools::TimeFlag(Time,5);
   int second=CommonTools::TimeFlag(Time,6);
   int whichday=year*10000+month*100+day;
   if(filename_in){
      whichday=CommonTools::GetTimeFromFileName(filename_in);
   }

   //load the files if not loaded
   if(whichday!=currentday){
      if(jdebug>0) printf("StatusDB::LocateFile: Load new files for day=%d\n",whichday);
      for(int ii=0;ii<failday.size();ii++){
         if(whichday==failday.at(ii)) return -1;
      }
      int nloaded=LoadFile(whichday,filename_in);
      if(jdebug>1) printf("StatusDB::LocateFile: %d files Loaded for day=%d\n",nloaded,whichday);
      if(nloaded<=0){
         failday.push_back(whichday);
         return -1;
      }
      //reset some variables
      if(UseDecodeData){
         if(tree) {delete tree; tree=0;}
         if(fstatus) fstatus->Close();
      }
      else{
         if(fp) fclose(fp);
      }
      currentfile=-2;
      currenttime=0;
      entryno=-1;
      failtime.clear();
      if(jdebug>2) printf("StatusDB::LocateFile: variables reseted currentday=%d currentfile=%d currenttime=%d entryno=%ld\n",currentday,currentfile,currenttime,entryno);
   }
   
   if(Time==currenttime) return currentfile;
   for(int ii=0;ii<failtime.size();ii++){
      if(Time==failtime.at(ii)) return -1;
   }
   //search the time in the current files
   int located=-2;
   int rectel=0;
   if(!filename_in){
      rectel=LocateTel(iTel0);
      if(jdebug>2) printf("StatusDB::LocateFile: begin search file according to the time=%d\n",Time);
      if(rectel>=0){
         for(int ii=0;ii<nfiles[rectel];ii++){
            int time1=filetime[rectel][ii];
            int time2=(ii+1>=nfiles[rectel])?(filetime[rectel][nfiles[rectel]-1]+MaxFileTime):filetime[rectel][ii+1];
            if((Time>=time1&&Time<time2)&&(Time-time1)<MaxFileTime){
               located=ii; break;
            }
         }
      }
   }
   else located=-1;
   if(jdebug>2) printf("StatusDB::LocateFile: search file finished time=%d ifile=%d currentfile=%d filename_in=%s\n",Time,located,currentfile,filename_in);
   if(located<-1) failtime.push_back(Time);
   else if(located!=currentfile){
      int ifile=(!filename_in)?fileindex[rectel][located]:0;
      if(UseDecodeData){
         if(currentfile>=-1){
            if(tree) {delete tree; tree=0;}
            if(fstatus) fstatus->Close();
         }
         //TFile* fstatus;
         if(!filename_in) fstatus=TFile::Open(Form("root://eos01.ihep.ac.cn/%s",filename[rectel][ifile]),"READ");
         else fstatus=TFile::Open(filename_in,"READ");
         fstatus->cd();
         tree=(TTree*)fstatus->Get("Status");
         if(jdebug>3) printf("StatusDB::LocateFile: Open file and SetBranchAddress ifile=%d tree=%p filename=%s\n",ifile,tree,filename_in?filename_in:filename[rectel][ifile]);
         tree->SetBranchAddress("iTel", &iTel);
         tree->SetBranchAddress("fpgaVersion", &fpgaVersion);
         tree->SetBranchAddress("f9mode", &f9mode);
         tree->SetBranchAddress("f9pattern", &f9pattern);
         tree->SetBranchAddress("DbVersion", &DbVersion);
         tree->SetBranchAddress("ClbVersion", &ClbVersion);
         tree->SetBranchAddress("clb_initial_Time", (Long64_t*)&clb_initial_Time);
         tree->SetBranchAddress("clb_initial_time", &clb_initial_time);
         tree->SetBranchAddress("fired_tube", &fired_tube);
         tree->SetBranchAddress("status_readback_Time", (Long64_t*)&status_readback_Time);
         tree->SetBranchAddress("status_readback_time", &status_readback_time);
         tree->SetBranchAddress("sipm", &sipm);
         tree->SetBranchAddress("mask", &mask);
         tree->SetBranchAddress("single_thresh", &single_thresh);
         tree->SetBranchAddress("record_thresh", &record_thresh);
         tree->SetBranchAddress("single_count", &single_count);
         tree->SetBranchAddress("single_time", &single_time);
         tree->SetBranchAddress("DbTemp", &DbTemp);
         tree->SetBranchAddress("HV", &HV);
         tree->SetBranchAddress("PreTemp", &PreTemp);
         tree->SetBranchAddress("BigResistence", &BigResistence);
         tree->SetBranchAddress("SmallResistence", &SmallResistence);
         tree->SetBranchAddress("ClbTime", &ClbTime);
         tree->SetBranchAddress("ClbTemp", &ClbTemp);
         if(jdebug>3) printf("StatusDB::LocateFile: SetBranchAddress finished. ifile=%d tree=%p\n",ifile,tree);
         //tree->GetEntry(0);
         //if(jdebug>4) printf("StatusDB::LocateFile: check entry 0 of %d, time0=%d\n",tree->GetEntries(),status_readback_Time);
      }
      else{
         if(!filename_in) strcpy(FILENAME,filename[rectel][ifile]);
         else strcpy(FILENAME,filename_in);
         fp=fopen(FILENAME,"rb");
         if(jdebug>3) printf("StatusDB::LocateFile: dat file opened. ifile=%d pfile=%p\n",ifile,fp);
      }
      currentfile=located;
   }
   return located;
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
   if(fabs(Time-currenttime)<=timemargin){ //just read the data from buf
      if(jdebug>1) printf("StatusDB::LocateBlk: Time=%d in the current buffer(currenttime=%d)\n",Time,currenttime);
      return true;
   }
   else if(Time>currenttime){ //search from right
      fseek(fp,buflength,1);
      while((size_of_read=fread((uint8_t *)buf,1,STATUS_BUF_LEN,fp))!=0){
         bool find=wfctaDecode->statusPackCheck(buf,size_of_read,9)>=0;//change "StatusPackCheck" to "statusPackCheck" by youzhiyong --- 2019 08 22
         uint64_t Time0=0;
         double time0=0;
         if(find){
            packSize = wfctaDecode->PackSize();
            Time0 = wfctaDecode->GetStatusReadbackTime(buf,packSize);
            time0 = wfctaDecode->GetStatusReadbacktime(buf,packSize);
         }
         if(fabs(Time-(int)Time0)<=timemargin){ //find it
            if(jdebug>1) printf("StatusDB::LocateBlk: find Time=%d in the right(currenttime=%d currenttime_new=%d)\n",Time,currenttime,Time0);
            currenttime=Time0;
            buflength=(unsigned long)size_of_read;
            fseek(fp,-buflength,1);
            ResetBuffer();
            Fill();
            return true;
         }
         else if(Time<Time0){ //there is no coresponding time
            if(jdebug>1) printf("StatusDB::LocateBlk: couldn't find Time=%d in the right(currenttime=%d currenttime_new=%d)\n",Time,currenttime,Time0);
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
         bool find=wfctaDecode->statusPackCheck(buf,-offset,9)>=0;//change "StatusPackCheck" to "statusPackCheck" by youzhiyong --- 2019 08 22
         uint64_t Time0=0;
         double time0=0;
         if(find){
            packSize = wfctaDecode->PackSize();
            Time0 = wfctaDecode->GetStatusReadbackTime(buf,packSize);
            time0 = wfctaDecode->GetStatusReadbacktime(buf,packSize);
         }
         if(fabs(Time-(int)Time0)<=timemargin){ //find it
            if(jdebug>1) printf("StatusDB::LocateBlk: find Time=%d in the left(currenttime=%d currenttime_new=%d)\n",Time,currenttime,Time0);
            currenttime=Time0;
            buflength=(long int)size_of_read;
            fseek(fp,-buflength,1);
            ResetBuffer();
            Fill();
            return true;
         }
         else if(Time>Time0&&Time0>0){  //there is no coresponding time
            if(jdebug>1) printf("StatusDB::LocateBlk: couldn't find Time=%d in the left(currenttime=%d currenttime_new=%d)\n",Time,currenttime,Time0);
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
   int find9=wfctaDecode->statusPackCheck(buf,buflength,9)+63;//change "StatusPackCheck" to "statusPackCheck" by youzhiyong --- 2019 08 22
   int find1=wfctaDecode->statusPackCheck(buf,buflength,1);//change "StatusPackCheck" to "statusPackCheck" by youzhiyong --- 2019 08 22
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
      find9=wfctaDecode->statusPackCheck(buf,-offset,9)+63; //F9 packet should exist //change "StatusPackCheck" to "statusPackCheck" by youzhiyong --- 2019 08 22
      if(find9+1!=-offset) printf("StatusDB::ResetBuffer: Error, there is no F9 packet in the buffer\n");
      find1=wfctaDecode->statusPackCheck(buf,-offset,1);//change "StatusPackCheck" to "statusPackCheck" by youzhiyong --- 2019 08 22
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
      status_pack_marker=wfctaDecode->StatusPackCheck(buf+offset,buflength-offset,0);//add one more argument "0"---by youzhiyong 2019 08 22
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
             wfctaDecode->GetBigRes(buf,offset+packSize,(float *)BigResistence);
             if(jdebug>5) printf("StatusDB::Fill:packetpos=%5ld, Fill 0x83 packet\n",offset+packSize);
             nfilled++;
             break;
           case 0x84:
             wfctaDecode->GetSmallRes(buf,offset+packSize,(float *)SmallResistence);
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
   if(jdebug>3) printf("StatusDB::Fill: Fill %d packets for time=%d at %ld of file %s\n",nfilled,currenttime,ftell(fp),FILENAME);
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
                 for(int jj=0;jj<MAXPMT;jj++) if(HV[jj]>max) {imax=jj; max=HV[jj];}
                 if(jdebug>3) printf("StatusDB::Fill: Fill HV max=%f(%d)\n",max,imax);
                 continue;
               case 0x82:
                 wfctaDecode->GetPreTemp(buf,packSize,(float *)PreTemp);
                 imax=-1;
                 max=-1000;
                 for(int jj=0;jj<MAXPMT;jj++) if(PreTemp[jj]>max) {imax=jj; max=PreTemp[jj];}
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
                 for(int jj=0;jj<MAXPMT;jj++) if(ClbTemp[jj]>max) {imax=jj; max=ClbTemp[jj];}
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

long int StatusDB::LocateEntry(int Time){
   if(!tree) return -2;
   if(jdebug>4) printf("StatusDB::LocateEntry: Begin search entry for Time=%d\n",Time);
   long int maxentry=tree->GetEntries();
   if(maxentry<1) return -2;
   long int entry0=entryno;
   if(entry0>=0&&entry0<maxentry){
      entryno=entry0;
      tree->GetEntry(entryno);
      int timei=status_readback_Time;
      if(jdebug>4) printf("StatusDB::LocateEntry: First try Time=%d entry=%ld time=%d\n",Time,entryno,timei);
      if(fabs(Time-timei)<=timemargin) return entryno;
   }
   if(entry0+1>=0&&entry0+1<maxentry){
      entryno=entry0+1;
      tree->GetEntry(entryno);
      int timei=status_readback_Time;
      if(jdebug>4) printf("StatusDB::LocateEntry: Second try Time=%d entry=%ld time=%d\n",Time,entryno,timei);
      if(fabs(Time-timei)<=timemargin) return entryno;
   }
   if(entry0-1>=0&&entry0-1<maxentry){
      entryno=entry0-1;
      tree->GetEntry(entryno);
      int timei=status_readback_Time;
      if(jdebug>4) printf("StatusDB::LocateEntry: Third try Time=%d entry=%ld time=%d\n",Time,entryno,timei);
      if(fabs(Time-timei)<=timemargin) return entryno;
   }
   entryno=0;
   tree->GetEntry(entryno);
   int time1=status_readback_Time;
   entryno=maxentry-1;
   tree->GetEntry(entryno);
   int time2=status_readback_Time;
   if(Time<time1-timemargin){
      int retval;
      if(Time<time1-3*60){
         //entryno=-1;
         entryno=entry0;
         tree->GetEntry(entryno);
         retval=-1;
      }
      else{
         entryno=0;
         tree->GetEntry(entryno);
         retval=entryno;
      }
      if(jdebug>4) printf("StatusDB::LocateEntry: Time too small, Time=%d entry=%ld time=%d return=%d\n",Time,entryno,time1,retval);
      return retval;
   }
   else if(Time>time2+timemargin){
      int retval;
      if(Time>time2+3*60){
         //entryno=-1;
         entryno=entry0;
         tree->GetEntry(entryno);
         retval=-1;
      }
      else{
         entryno=maxentry-1;
         tree->GetEntry(entryno);
         retval=entryno;
      }
      if(jdebug>4) printf("StatusDB::LocateEntry: Time too large, Time=%d entry=%ld time=%d return=%d\n",Time,entryno,time1,retval);
      return retval;
   }
   else if(fabs(Time-time1)<=timemargin){
      entryno=0;
      tree->GetEntry(entryno);
      if(jdebug>5) printf("StatusDB::LocateEntry: try minimum entry, Time=%d entry=%ld time=%d\n",Time,entryno,time1);
      return entryno;
   }
   else if(fabs(Time-time2)<=timemargin){
      entryno=maxentry-1;
      tree->GetEntry(entryno);
      if(jdebug>5) printf("StatusDB::LocateEntry: try maximum entry, Time=%d entry=%ld time=%d\n",Time,entryno,time2);
      return entryno;
   }
   else{
      int entry1=0;
      int entry2=maxentry-1;
      int entryi=(entry1+entry2)/2;
      entryno=entryi;
      tree->GetEntry(entryno);
      int timei=status_readback_Time;
      while(fabs(Time-timei)>timemargin){
         if(Time<timei) {entry2=entryi; time2=timei;}
         else {entry1=entryi; time1=timei;}
         entryi=(entry1+entry2)/2;
         if(entryi==entry1||entryi==entry2){
            if(fabs(Time-time1)<fabs(Time-time2)) entryi=entry1;
            else entryi=entry2;
         }
         entryno=entryi;
         tree->GetEntry(entryno);
         timei=status_readback_Time;
         if(jdebug>5) printf("StatusDB::LocateEntry: ith try Time=%d entry=%ld(entry1=%ld entry2=%ld) time=%d\n",Time,entryno,entry1,entry2,timei);
         if(fabs(entry1-entry2)<1.5) break;
      }
      if(fabs(Time-timei)>timemargin){
         if(fabs(time1-time2)<3*60){
            if(jdebug>4) printf("StatusDB::LocateEntry: Missing a short time, Time=%d entry=%ld time=%d (entry1=%ld time1=%d,entry2=%ld time2=%d)\n",Time,entryno,timei,entry1,time1,entry2,time2);
            return entryi;
         }
         else{
            //entryno=-1;
            entryno=entry0;
            tree->GetEntry(entryno);
            if(jdebug>4) printf("StatusDB::LocateEntry: Missing some time, Time=%d entry=%ld time=%d (entry1=%ld time1=%d,entry2=%ld time2=%d)\n",Time,entryi,timei,entry1,time1,entry2,time2);
            return -1;
         }
      }
      else return entryno;
   }
}

long int StatusDB::Locate(int iTel0,int Time,char* filename_in,double time){
   int iindex=LocateFile(iTel0,Time,filename_in,time);
   if(iindex<-1) return -2;
   if(UseDecodeData){
      if(jdebug>5) printf("StatusDB::Locate: currentime=%d entryno=%ld Time=%d\n",currenttime,entryno,Time);
      if(fabs(Time-currenttime)<1){
         return entryno;
      }
      long int entry=LocateEntry(Time);
      if(entry>=0){
         //currentfile=iindex;
         currenttime=Time;
         entryno=entry;
      }
      else failtime.push_back(Time);
      if(jdebug>4) printf("StatusDB::Locate: entry=%ld located for iTel=%d Time=%d currentime=%d entryno=%ld\n",entry,iTel0,Time,currenttime,entryno);
      return entry;
   }
   else{
      return LocateBlk(Time,time)?1:-1;
   }
}

int StatusDB::GetVersion(int iTel0,int Time,int i,char* filename_in){
   if(Locate(iTel0,Time,filename_in)<0) return -1000;
   return fpgaVersion[i];
}
long StatusDB::GetClbInitTime(int iTel0,int Time,char* filename_in){
   if(Locate(iTel0,Time,filename_in)<0) return 0;
   return clb_initial_Time;
}
double StatusDB::GetClbInittime(int iTel0,int Time,char* filename_in){
   if(Locate(iTel0,Time,filename_in)<0) return 0;
   return clb_initial_time;
}
int StatusDB::GetFiredTube(int iTel0,int Time,char* filename_in){
   if(Locate(iTel0,Time,filename_in)<0) return 0;
   return fired_tube;
}
long StatusDB::GetReadbackTime(int iTel0,int Time,char* filename_in){
   if(Locate(iTel0,Time,filename_in)<0) return 0;
   return status_readback_Time;
}
double StatusDB::GetReadbacktime(int iTel0,int Time,char* filename_in){
   if(Locate(iTel0,Time,filename_in)<0) return 0;
   return status_readback_time;
}
short StatusDB::GetSingleThrd(int iTel0,int Time,int i,char* filename_in){
   if(Locate(iTel0,Time,filename_in)<0) return -1000;
   return single_thresh[i];
}
short StatusDB::GetRecordThrd(int iTel0,int Time,int i,char* filename_in){
   if(Locate(iTel0,Time,filename_in)<0) return -1000;
   return record_thresh[i];
}
long StatusDB::GetSingleCount(int iTel0,int Time,int i,char* filename_in){
   if(Locate(iTel0,Time,filename_in)<0) return -1000;
   return single_count[i];
}
float StatusDB::GetDbTemp(int iTel0,int Time,int i,char* filename_in){
   if(Locate(iTel0,Time,filename_in)<0) return -1000;
   return single_count[i];
}
long StatusDB::GetSingleTime(int iTel0,int Time,int i,char* filename_in){
   if(Locate(iTel0,Time,filename_in)<0) return -1000;
   return single_time[i];
}
float StatusDB::GetHV(int iTel0,int Time,int i,char* filename_in){
   if(Locate(iTel0,Time,filename_in)<0) return -1000;
   return HV[i];
}
float StatusDB::GetPreTemp(int iTel0,int Time,int i,char* filename_in){
   if(Locate(iTel0,Time,filename_in)<0) return -1000;
   return PreTemp[i];
}
float StatusDB::GetBigR(int iTel0,int Time,int i,char* filename_in){
   if(Locate(iTel0,Time,filename_in)<0) return -1000;
   return BigResistence[i];
}
float StatusDB::GetSmallR(int iTel0,int Time,int i,char* filename_in){
   if(Locate(iTel0,Time,filename_in)<0) return -1000;
   return SmallResistence[i];
}
long StatusDB::GetClbTime(int iTel0,int Time,int i,char* filename_in){
   if(Locate(iTel0,Time,filename_in)<0) return -1000;
   return ClbTime[i];
}
float StatusDB::GetClbTemp(int iTel0,int Time,int i,char* filename_in){
   if(Locate(iTel0,Time,filename_in)<0) return -1000;
   return ClbTemp[i];
}
