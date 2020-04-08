#include "WFCTASCDB.h"
#include "WFCTAEvent.h"
#include "TString.h"
int WFCTASCDB::TimeDelay=0;
int WFCTASCDB::jdebug=0;
WFCTASCDB* WFCTASCDB::_Head=0;
char WFCTASCDB::DirName[200]="/eos/user/w/wfcta/wfcta/wfcta/slowctrl/rootdata";
void WFCTASCDB::Init(){
   curTel=0;
   fin=0;
   tree=0;
   entries=-1;
   curentry=-1;
   Reset();
}
void WFCTASCDB::Reset(){
   T1MJD=0;
   T1HV=0;
   T1HVCurrent=0;
   T1NV=0;
   T1NVCurrent=0;
   T1P1V=0;
   T1P1VCurrent=0;
   T1P2V=0;
   T1P2VCurrent=0;
   T1P3V=0;
   T1P3VCurrent=0;
   T1Dledtem=-1000;
   T1DledDStem=-1000;
   T1DledDrtem=-1000;
   T1DledDrDStem=-1000;
   T1Rledtem=-1000;
   T1RledDStem=-1000;
   T1RledDrtem=-1000;
   T1RledDrDStem=-1000;
   T1CaImtem=-1000;
   T1CaExtem=-1000;
   T1Door1Deg=-1000;
   T1Door2Deg=-1000;
   T1HM=-1000;
   T1TemF=-1000;
   T1TemM=-1000;
   T1TemB=-1000;
   T1IncXDeg=-1000;
   T1IncYDeg=-1000;
   T1EleDeg=-1000;
   T1DledPW=0;
   T1DledPF=0;
   T1RledPW=0;
   T1RledPF=0;
}
void WFCTASCDB::Clear(){
   if(tree) {delete tree; tree=0;}
   if(fin) {fin->Close(); fin=0;}
}
void WFCTASCDB::SetDirName(char* dirname){
   strcpy(DirName,dirname);
}
WFCTASCDB* WFCTASCDB::GetHead(char* dirname){
   if(!_Head){
      if(dirname) WFCTASCDB::SetDirName(dirname);
      _Head=new WFCTASCDB();
   }
   return _Head;
}
void WFCTASCDB::SetBranchAddress(){
   Reset();
   if(!tree) return;
   if(tree->GetBranchStatus("T1MJD")) tree->SetBranchAddress("T1MJD",&T1MJD);
   if(tree->GetBranchStatus("T1HV")) tree->SetBranchAddress("T1HV",&T1HV);
   if(tree->GetBranchStatus("T1HVCurrent")) tree->SetBranchAddress("T1HVCurrent",&T1HVCurrent);
   if(tree->GetBranchStatus("T1NV")) tree->SetBranchAddress("T1NV",&T1NV);
   if(tree->GetBranchStatus("T1NVCurrent")) tree->SetBranchAddress("T1NVCurrent",&T1NVCurrent);
   if(tree->GetBranchStatus("T1P1V")) tree->SetBranchAddress("T1P1V",&T1P1V);
   if(tree->GetBranchStatus("T1P1VCurrent")) tree->SetBranchAddress("T1P1VCurrent",&T1P1VCurrent);
   if(tree->GetBranchStatus("T1P2V")) tree->SetBranchAddress("T1P2V",&T1P2V);
   if(tree->GetBranchStatus("T1P2VCurrent")) tree->SetBranchAddress("T1P2VCurrent",&T1P2VCurrent);
   if(tree->GetBranchStatus("T1P3V")) tree->SetBranchAddress("T1P3V",&T1P3V);
   if(tree->GetBranchStatus("T1P3VCurrent")) tree->SetBranchAddress("T1P3VCurrent",&T1P3VCurrent);
   if(tree->GetBranchStatus("T1Dledtem")) tree->SetBranchAddress("T1Dledtem",&T1Dledtem);
   if(tree->GetBranchStatus("T1DledDStem")) tree->SetBranchAddress("T1DledDStem",&T1DledDStem);
   if(tree->GetBranchStatus("T1DledDrtem")) tree->SetBranchAddress("T1DledDrtem",&T1DledDrtem);
   if(tree->GetBranchStatus("T1DledDrDStem")) tree->SetBranchAddress("T1DledDrDStem",&T1DledDrDStem);
   if(tree->GetBranchStatus("T1Rledtem")) tree->SetBranchAddress("T1Rledtem",&T1Rledtem);
   if(tree->GetBranchStatus("T1RledDStem")) tree->SetBranchAddress("T1RledDStem",&T1RledDStem);
   if(tree->GetBranchStatus("T1RledDrtem")) tree->SetBranchAddress("T1RledDrtem",&T1RledDrtem);
   if(tree->GetBranchStatus("T1RledDrDStem")) tree->SetBranchAddress("T1RledDrDStem",&T1RledDrDStem);
   if(tree->GetBranchStatus("T1CaImtem")) tree->SetBranchAddress("T1CaImtem",&T1CaImtem);
   if(tree->GetBranchStatus("T1CaExtem")) tree->SetBranchAddress("T1CaExtem",&T1CaExtem);
   if(tree->GetBranchStatus("T1Door1Deg")) tree->SetBranchAddress("T1Door1Deg",&T1Door1Deg);
   if(tree->GetBranchStatus("T1Door2Deg")) tree->SetBranchAddress("T1Door2Deg",&T1Door2Deg);
   if(tree->GetBranchStatus("T1HM")) tree->SetBranchAddress("T1HM",&T1HM);
   if(tree->GetBranchStatus("T1TemF")) tree->SetBranchAddress("T1TemF",&T1TemF);
   if(tree->GetBranchStatus("T1TemM")) tree->SetBranchAddress("T1TemM",&T1TemM);
   if(tree->GetBranchStatus("T1TemB")) tree->SetBranchAddress("T1TemB",&T1TemB);
   if(tree->GetBranchStatus("T1IncXDeg")) tree->SetBranchAddress("T1IncXDeg",&T1IncXDeg);
   if(tree->GetBranchStatus("T1IncYDeg")) tree->SetBranchAddress("T1IncYDeg",&T1IncYDeg);
   if(tree->GetBranchStatus("T1EleDeg")) tree->SetBranchAddress("T1EleDeg",&T1EleDeg);
   if(tree->GetBranchStatus("T1DledPW")) tree->SetBranchAddress("T1DledPW",&T1DledPW);
   if(tree->GetBranchStatus("T1DledPF")) tree->SetBranchAddress("T1DledPF",&T1DledPF);
   if(tree->GetBranchStatus("T1RledPW")) tree->SetBranchAddress("T1RledPW",&T1RledPW);
   if(tree->GetBranchStatus("T1RledPF")) tree->SetBranchAddress("T1RledPF",&T1RledPF);
}
long int WFCTASCDB::LoadDay(int Day,int iTel){
   //if(tree){
   //   if((!fin)||entries<=0){
   //      delete tree; tree=0; entries=-1; curentry=-1; curTel=0;
   //      if(fin) fin->Close();
   //      fin=0;
   //   }
   //   else{
   //      if(curentry<0) {curentry=0; tree->GetEntry(curentry);}
   //      int Time=CommonTools::ConvertMJD2Time(T1MJD);
   //      int curday=CommonTools::GetTelDay(Time);
   //      if(Day==curday&&curTel==iTel) return entries;
   //      else{ //load another day
   //         delete tree; tree=0; entries=-1; curentry=-1; curTel=0;
   //         fin->Close();
   //         fin=0;
   //      }
   //   }
   //}

   if(fin){
      char filename0[200]="";
      strcpy(filename0,fin->GetName());
      int size=strlen(filename0);
      int loc0=0;

      char telinfo[10];
      int ncount=0;
      int nchar=0;
      for(int ii=0;ii<size;ii++){
         if(filename0[ii]=='/') loc0=ii+1;
         if(filename0[ii]=='_'){
            ncount=1; nchar=0;
            continue;
         }
         if(ncount>0&&filename0[ii]=='.') ncount=0;
         if(ncount>0) telinfo[nchar++]=filename0[ii];
      }
      telinfo[nchar++]='\0';

      char year[10];
      char month[10];
      char day[10];
      ncount=0;
      int nchar1=0,nchar2=0,nchar3=0;
      for(int ii=loc0;ii<size;ii++){
         if(filename0[ii]=='-') {ncount++; continue;}
         if(filename0[ii]=='_') break;
         if(ncount==0){
            year[nchar1++]=filename0[ii];
         }
         else if(ncount==1){
            month[nchar2++]=filename0[ii];
         }
         else if(ncount==2){
            day[nchar3++]=filename0[ii];
         }
      }
      year[nchar1++]='\0';
      month[nchar2++]='\0';
      day[nchar3++]='\0';
      if(jdebug>0) printf("WFCTASCDB::LoadDay: iTel=%d(%d) Day=%04d-%02d-%02d(%d)\n",atoi(telinfo),iTel,atoi(year),atoi(month),atoi(day),Day);

      int curday=atoi(year)*10000+atoi(month)*100+atoi(day);
      if(Day==curday&&atoi(telinfo)==iTel) return entries;
      else{ //load another day
         if(tree){delete tree; tree=0;} entries=-1; curentry=-1; curTel=0;
         fin->Close();
         fin=0;
      }
   }

   char filename[200]="";
   TString ss(DirName); if(ss.BeginsWith("/eos/")) strcpy(filename,"root:://eos01.ihep.ac.cn/");
   strcat(filename,Form("%s/%04d-%d-%d_%d.root",DirName,Day/10000,(Day%10000)/100,(Day%100),iTel));
   TDirectory* gdir=gDirectory;
   fin=TFile::Open(filename,"READ");
   if(fin) fin->cd();
   tree=fin?(TTree*)fin->Get("Tel01_info"):0;
   SetBranchAddress();
   entries=tree?(tree->GetEntries()-1):-1;
   curentry=-1; curTel=iTel;
   if(entries>0){curentry=0; tree->GetEntry(curentry);}
   gdir->cd();
   return entries;
}
int WFCTASCDB::LoadTime(int time){
   if(entries<=0) return -10;
   if(curentry<0){curentry=0; tree->GetEntry(curentry);}
   time-=TimeDelay;
   long int entryi=curentry;
   int timei=(int)(CommonTools::ConvertMJD(T1MJD)+0.5);
   if(jdebug>2) printf("WFCTASCDB::LoadTime: test0,entry=%ld(%ld) Time=%d(%d) MJD=%.4lf\n",curentry,entries,timei,time,T1MJD);
   int timei2=0;
   long entryi2=-1;
   if(time==timei){
      if(jdebug>0) printf("WFCTASCDB::LoadTime: test0,return currrent entry. entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei,time);
      return 1;
   }
   else if(time>timei){
      if(curentry>=entries-1){
         if(jdebug>0) printf("WFCTASCDB::LoadTime: test0,over last entry. currrent entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei,time);
         return -3;
      }
      entryi2=entryi+1;
      curentry=entryi2;
      tree->GetEntry(curentry);
      timei2=(int)(CommonTools::ConvertMJD(T1MJD)+0.5);
      if(jdebug>2) printf("WFCTASCDB::LoadTime: test1,entry=%ld(%ld) Time=%d(%d) MJD=%.4lf\n",curentry,entries,timei2,time,T1MJD);
      if(time==timei2){
         if(jdebug>0) printf("WFCTASCDB::LoadTime: test1,return currrent entry. entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei2,time);
         return 1;
      }
      else if(time>timei2){
         entryi2=entries-1;
         curentry=entryi2;
         tree->GetEntry(curentry);
         timei2=(int)(CommonTools::ConvertMJD(T1MJD)+0.5);
         if(jdebug>2) printf("WFCTASCDB::LoadTime: test2,entry=%ld(%ld) Time=%d(%d) MJD=%.4lf\n",curentry,entries,timei2,time,T1MJD);
         if(time==timei2){
            if(jdebug>0) printf("WFCTASCDB::LoadTime: test2,return currrent entry. entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei2,time);
            return 1;
         }
         else if(time>timei2){
            if(jdebug>0) printf("WFCTASCDB::LoadTime: test2,over last entry. entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei2,time);
            return -3;
         }
      }
   }
   else{
      if(curentry<=0){
         if(jdebug>0) printf("WFCTASCDB::LoadTime: test0,before first entry. currrent entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei,time);
         return -2;
      }
      entryi2=entryi-1;
      curentry=entryi2;
      tree->GetEntry(curentry);
      timei2=(int)(CommonTools::ConvertMJD(T1MJD)+0.5);
      if(jdebug>2) printf("WFCTASCDB::LoadTime: test1,entry=%ld(%ld) Time=%d(%d) MJD=%.4lf\n",curentry,entries,timei2,time,T1MJD);
      if(time==timei2){
         if(jdebug>0) printf("WFCTASCDB::LoadTime: test1,return currrent entry. entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei2,time);
         return 1;
      }
      else if(time<timei2){
         entryi2=0;
         curentry=entryi2;
         tree->GetEntry(curentry);
         timei2=(int)(CommonTools::ConvertMJD(T1MJD)+0.5);
         if(jdebug>2) printf("WFCTASCDB::LoadTime: test2,entry=%ld(%ld) Time=%d(%d) MJD=%.4lf\n",curentry,entries,timei2,time,T1MJD);
         if(time==timei2){
            if(jdebug>0) printf("WFCTASCDB::LoadTime: test2,return currrent entry. entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei2,time);
            return 1;
         }
         else if(time<timei2){
            if(jdebug>0) printf("WFCTASCDB::LoadTime: test2,before first entry. entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei2,time);
            return -2;
         }
      }
   }

   int time1=entryi<entryi2?timei:timei2;
   long int entry1=entryi<entryi2?entryi:entryi2;
   int time2=entryi<entryi2?timei2:timei;
   long int entry2=entryi<entryi2?entryi2:entryi;

   int nloop=0;
   while(entry2>entry1+1){
      if(jdebug>1) printf("WFCTASCDB::LoadTime: loop%d time=%d entries=%ld ref1={%d,%ld} ref2={%d,%ld}\n",nloop,time,entries,time1,entry1,time2,entry2);
      long int entry_test=(entry1+entry2)/2;
      curentry=entry_test;
      tree->GetEntry(curentry);
      int time_test=(int)(CommonTools::ConvertMJD(T1MJD)+0.5);
      if(jdebug>3) printf("WFCTASCDB::LoadTime: loop%d,entry=%ld(%ld) Time=%d(%d) MJD=%.4lf\n",nloop,curentry,entries,time_test,time,T1MJD);
      if(time==time_test){
         if(jdebug>1) printf("WFCTASCDB::LoadTime: loop%d,return currrent entry. entry=%ld(%ld) Time=%d(%d)\n",nloop,curentry,entries,time_test,time);
         return 1;
      }
      else if(time>time_test){
         time1=time_test;
         entry1=entry_test;
      }
      else{
         time2=time_test;
         entry2=entry_test;
      }
      nloop++;
   }
   if(time==time1||time==time2) return 1;
   else if(time>time1&&time<time2) return -1;
   else if(time<time1) return -2;
   else return -3;
}
int WFCTASCDB::LoadFromrabbitTime(double time,int iTel){
   int daytime=(int)(time-TAI2UTC);
   int year=CommonTools::TimeFlag(daytime,1);
   int month=CommonTools::TimeFlag(daytime,2);
   int day=CommonTools::TimeFlag(daytime,3);
   int curday=year*10000+month*100+day;
   long int totentry=LoadDay(curday,iTel);
   return LoadTime((int)(time+0.5));
}
int WFCTASCDB::LoadEntry(long int entry){
   if(entry<0||entry>=entries) return 0;
   if(!tree) return -1;
   curentry=entry;
   tree->GetEntry(curentry);
   int Time=(int)(CommonTools::ConvertMJD(T1MJD)+0.5);
   return Time;
}

void WFCTASCDB::Dump(){
   int Time=(int)(CommonTools::ConvertMJD(T1MJD)+0.5)-TAI2UTC+TimeDelay;
   int year=CommonTools::TimeFlag(Time,1);
   int month=CommonTools::TimeFlag(Time,2);
   int day=CommonTools::TimeFlag(Time,3);
   int hour=CommonTools::TimeFlag(Time,4);
   int min=CommonTools::TimeFlag(Time,5);
   int sec=CommonTools::TimeFlag(Time,6);
   printf("FileName=%s tree=%p entries=%ld cur_entry=%ld iTel=%d cur_time=%d(%04d-%02d-%02d %02d:%02d:%02d)\n",fin?fin->GetName():"",tree,entries,curentry,curTel,Time,year,month,day,hour,min,sec);
   printf("HV=%.2lf HVCurrent=%.4lf\n",T1HV,T1HVCurrent);
   printf("NV=%.2lf NVCurrent=%.4lf\n",T1NV,T1NVCurrent);
   printf("P3V=%.2lf P3VCurrent=%.4lf\n",T1P3V,T1P3VCurrent);
   printf("DLed temp: %d %.2lf %.2lf %.2lf %.2lf\n",IsDLed(),T1Dledtem,T1DledDStem,T1DledDrtem,T1DledDrDStem);
   printf("RLed temp: %d %.2lf %.2lf %.2lf %.2lf\n",IsRLed(),T1Rledtem,T1RledDStem,T1RledDrtem,T1DledDrDStem);
   printf("DoorOpened: %.1lf %.1lf\n",T1Door1Deg,T1Door2Deg);
   printf("DoorTemp: %.2lf %.2lf %.2lf\n",T1TemF,T1TemM,T1TemB);
   printf("Inc Angle: %.2lf %.2lf %.2lf\n",T1IncXDeg,T1IncYDeg,T1EleDeg);
   printf("DLed: width=%d freq=%d\n",T1DledPW,T1DledPF);
   printf("RLed: width=%d freq=%d\n\n",T1RledPW,T1RledPF);
}

bool WFCTASCDB::IsDLed(WFCTAEvent* pev){
   if(pev) LoadFromrabbitTime(pev->rabbitTime+pev->rabbittime*2.0e-8,pev->iTel);
   if(!tree) return false;
   if(!tree->GetBranchStatus("T1DledDStem")) return false;
   if(!tree->GetBranchStatus("T1DledDrDStem")) return false;
   return (T1Dledtem<100);
}
bool WFCTASCDB::IsRLed(WFCTAEvent* pev){
   if(pev) LoadFromrabbitTime(pev->rabbitTime+pev->rabbittime*2.0e-8,pev->iTel);
   if(!tree) return false;
   if(!tree->GetBranchStatus("T1RledDStem")) return false;
   if(!tree->GetBranchStatus("T1RledDrDStem")) return false;
   return (T1Rledtem<100);
}
bool WFCTASCDB::DoorOpened(WFCTAEvent* pev){
   if(pev) LoadFromrabbitTime(pev->rabbitTime+pev->rabbittime*2.0e-8,pev->iTel);
   if(!tree) return false;
   if(!tree->GetBranchStatus("T1Door1Deg")) return false;
   if(!tree->GetBranchStatus("T1Door2Deg")) return false;
   return (T1Door1Deg>260.)&&(T1Door2Deg>260.);
}
