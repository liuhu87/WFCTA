#include "CalibLED.h"
#include "WFCTAEvent.h"
bool CalibLED::ForceCorr=false;
int CalibLED::TimeDelay=0;
int CalibLED::jdebug=0;
CalibLED* CalibLED::_Head=0;
char CalibLED::DirName[200]="/eos/lhaaso/rec/wfcta/LED_Driver_Correct_Factor_new";
void CalibLED::Init(){
   curTel=0;
   fin=0;
   tree=0;
   entries=-1;
   curentry=-1;
   Reset();
}
void CalibLED::Reset(){
   T1MJD=0;
   LED_DRF=1;
   T1HV=0;
   T1HVCurrent=0;
   T1DledDrDStem=-1000;
}
void CalibLED::Clear(){
   if(tree) {delete tree; tree=0;}
   if(fin) {fin->Close(); fin=0;}
}
void CalibLED::SetDirName(char* dirname){
   if(dirname) strcpy(DirName,dirname);
}
CalibLED* CalibLED::GetHead(char* dirname){
   if(!_Head){
      _Head=new CalibLED();
      if(dirname) SetDirName(dirname);
   }
   return _Head;
}
void CalibLED::SetBranchAddress(){
   Reset();
   if(!tree) return;
   if(tree->GetBranchStatus("T1MJD")) tree->SetBranchAddress("T1MJD",&T1MJD);
   if(tree->GetBranchStatus("LED_DRF")) tree->SetBranchAddress("LED_DRF",&LED_DRF);
   if(tree->GetBranchStatus("T1HVCurrent")) tree->SetBranchAddress("T1HVCurrent",&T1HVCurrent);
   if(tree->GetBranchStatus("T1HV")) tree->SetBranchAddress("T1HV",&T1HV);
   if(tree->GetBranchStatus("T1DledDrDStem")) tree->SetBranchAddress("T1DledDrDStem",&T1DledDrDStem);
}
long int CalibLED::LoadDay(int Day,int iTel){
   if(fin){
      char filename0[200]="";
      strcpy(filename0,fin->GetName());
      int size=strlen(filename0);
      int loc0=0;

      for(int ii=0;ii<size;ii++){
         if(filename0[ii]=='/') loc0=ii+1;
      }

      char telinfo[10];
      int nchar=0;
      char year[10];
      char month[10];
      char day[10];
      int ncount=0;
      int nchar1=0,nchar2=0,nchar3=0;
      for(int ii=loc0;ii<size;ii++){
         if(filename0[ii]=='_') {ncount++; continue;}
         if(filename0[ii]=='.') break;
         if(ncount==0){
            year[nchar1++]=filename0[ii];
         }
         else if(ncount==1){
            month[nchar2++]=filename0[ii];
         }
         else if(ncount==2){
            day[nchar3++]=filename0[ii];
         }
         else if(ncount==3){
            telinfo[nchar++]=filename0[ii];
         }
      }
      telinfo[nchar++]='\0';
      year[nchar1++]='\0';
      month[nchar2++]='\0';
      day[nchar3++]='\0';
      if(jdebug>0) printf("CalibLED::LoadDay: iTel=%d(%d) Day=%04d-%02d-%02d(%d)\n",atoi(telinfo),iTel,atoi(year),atoi(month),atoi(day),Day);

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
   strcat(filename,Form("%s/%04d_%02d_%02d_%02d.root",DirName,Day/10000,(Day%10000)/100,(Day%100),iTel));
   TDirectory* gdir=gDirectory;
   fin=TFile::Open(filename,"READ");
   if(fin) fin->cd();
   tree=fin?(TTree*)fin->Get("slow_control"):0;
   SetBranchAddress();
   entries=tree?(tree->GetEntries()-1):-1;
   curentry=-1; curTel=iTel;
   if(entries>0){curentry=0; tree->GetEntry(curentry);}
   gdir->cd();
   return entries;
}
int CalibLED::LoadTime(int time){
   if(entries<=0) return -10;
   if(curentry<0){curentry=0; tree->GetEntry(curentry);}
   time-=TimeDelay;
   long int entryi=curentry;
   int timei=(int)(CommonTools::ConvertMJD(T1MJD)+0.5);
   if(jdebug>2) printf("CalibLED::LoadTime: test0,entry=%ld(%ld) Time=%d(%d) MJD=%.4lf\n",curentry,entries,timei,time,T1MJD);
   int timei2=0;
   long entryi2=-1;
   if(time==timei){
      if(jdebug>0) printf("CalibLED::LoadTime: test0,return currrent entry. entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei,time);
      return 1;
   }
   else if(time>timei){
      if(curentry>=entries-1){
         if(jdebug>0) printf("CalibLED::LoadTime: test0,over last entry. currrent entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei,time);
         return -3;
      }
      entryi2=entryi+1;
      curentry=entryi2;
      tree->GetEntry(curentry);
      timei2=(int)(CommonTools::ConvertMJD(T1MJD)+0.5);
      if(jdebug>2) printf("CalibLED::LoadTime: test1,entry=%ld(%ld) Time=%d(%d) MJD=%.4lf\n",curentry,entries,timei2,time,T1MJD);
      if(time==timei2){
         if(jdebug>0) printf("CalibLED::LoadTime: test1,return currrent entry. entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei2,time);
         return 1;
      }
      else if(time>timei2){
         entryi2=entries-1;
         curentry=entryi2;
         tree->GetEntry(curentry);
         timei2=(int)(CommonTools::ConvertMJD(T1MJD)+0.5);
         if(jdebug>2) printf("CalibLED::LoadTime: test2,entry=%ld(%ld) Time=%d(%d) MJD=%.4lf\n",curentry,entries,timei2,time,T1MJD);
         if(time==timei2){
            if(jdebug>0) printf("CalibLED::LoadTime: test2,return currrent entry. entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei2,time);
            return 1;
         }
         else if(time>timei2){
            if(jdebug>0) printf("CalibLED::LoadTime: test2,over last entry. entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei2,time);
            return -3;
         }
      }
   }
   else{
      if(curentry<=0){
         if(jdebug>0) printf("CalibLED::LoadTime: test0,before first entry. currrent entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei,time);
         return -2;
      }
      entryi2=entryi-1;
      curentry=entryi2;
      tree->GetEntry(curentry);
      timei2=(int)(CommonTools::ConvertMJD(T1MJD)+0.5);
      if(jdebug>2) printf("CalibLED::LoadTime: test1,entry=%ld(%ld) Time=%d(%d) MJD=%.4lf\n",curentry,entries,timei2,time,T1MJD);
      if(time==timei2){
         if(jdebug>0) printf("CalibLED::LoadTime: test1,return currrent entry. entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei2,time);
         return 1;
      }
      else if(time<timei2){
         entryi2=0;
         curentry=entryi2;
         tree->GetEntry(curentry);
         timei2=(int)(CommonTools::ConvertMJD(T1MJD)+0.5);
         if(jdebug>2) printf("CalibLED::LoadTime: test2,entry=%ld(%ld) Time=%d(%d) MJD=%.4lf\n",curentry,entries,timei2,time,T1MJD);
         if(time==timei2){
            if(jdebug>0) printf("CalibLED::LoadTime: test2,return currrent entry. entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei2,time);
            return 1;
         }
         else if(time<timei2){
            if(jdebug>0) printf("CalibLED::LoadTime: test2,before first entry. entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei2,time);
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
      if(jdebug>1) printf("CalibLED::LoadTime: loop%d time=%d entries=%ld ref1={%d,%ld} ref2={%d,%ld}\n",nloop,time,entries,time1,entry1,time2,entry2);
      long int entry_test=(entry1+entry2)/2;
      curentry=entry_test;
      tree->GetEntry(curentry);
      int time_test=(int)(CommonTools::ConvertMJD(T1MJD)+0.5);
      if(jdebug>3) printf("CalibLED::LoadTime: loop%d,entry=%ld(%ld) Time=%d(%d) MJD=%.4lf\n",nloop,curentry,entries,time_test,time,T1MJD);
      if(time==time_test){
         if(jdebug>1) printf("CalibLED::LoadTime: loop%d,return currrent entry. entry=%ld(%ld) Time=%d(%d)\n",nloop,curentry,entries,time_test,time);
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
int CalibLED::LoadFromrabbitTime(double time,int iTel){
   int daytime=(int)(time-TAI2UTC);
   int curday=CommonTools::GetTelDay(daytime);
   long int totentry=LoadDay(curday,iTel);
   return LoadTime((int)(time+0.5));
}
int CalibLED::LoadEntry(long int entry){
   if(entry<0||entry>=entries) return 0;
   if(!tree) return -1;
   curentry=entry;
   tree->GetEntry(curentry);
   int Time=(int)(CommonTools::ConvertMJD(T1MJD)+0.5);
   return Time;
}

void CalibLED::Dump(){
   int Time=(int)(CommonTools::ConvertMJD(T1MJD)+0.5)-TAI2UTC+TimeDelay;
   int year=CommonTools::TimeFlag(Time,1);
   int month=CommonTools::TimeFlag(Time,2);
   int day=CommonTools::TimeFlag(Time,3);
   int hour=CommonTools::TimeFlag(Time,4);
   int min=CommonTools::TimeFlag(Time,5);
   int sec=CommonTools::TimeFlag(Time,6);
   printf("CalibLED::Dump: FileName=%s tree=%p entries=%ld cur_entry=%ld iTel=%d cur_time=%d(%04d-%02d-%02d %02d:%02d:%02d)\n",fin?fin->GetName():"",tree,entries,curentry,curTel,Time,year,month,day,hour,min,sec);
   printf("HV=%.2lf HVCurrent=%.4lf\n",T1HV,T1HVCurrent);
   printf("DLed Corr: %.4lf %.2lf\n",LED_DRF,T1DledDrDStem);
}

double CalibLED::DoLedDriveTempCorr(double input,int isipm,double time,int iTel){
   if(time>1300000000) LoadFromrabbitTime(time,iTel);
   double ImageX,ImageY;
   WFCTAEvent::GetImageXYCoo(isipm,ImageX,ImageY,-1,false);
   double angle=sqrt(pow(ImageX,2)+pow(ImageY,2));
   double corr=pow(cos(angle),4);
   double res=input/corr;
   if(!tree) return ForceCorr?-1.:res;
   if(!tree->GetBranchStatus("LED_DRF")) return ForceCorr?-1.:res;
   return (res*LED_DRF);
}
double CalibLED::DoLedDriveTempCorr(double input,int isipm,WFCTAEvent* pev){
   if(pev) return DoLedDriveTempCorr(input,isipm,pev->rabbitTime+pev->rabbittime*2.0e-8,pev->iTel);
   else return DoLedDriveTempCorr(input,isipm,0,0);
}
