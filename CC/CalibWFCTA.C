#include "CalibWFCTA.h"
//#include "camera.h"
extern void SC_Channel2SiPM(short F_DB, short mChannel, short *mSiPM);
#include "TFile.h"
#include "TString.h"
#include "CalibLED.h"
#include <iostream>
using namespace std;

CalibWFCTA* CalibWFCTA::_Head=0;
int CalibWFCTA::UseSiPMCalibVer=1;
bool CalibWFCTA::ForceCorr=true;
int CalibWFCTA::jdebug=0;
char CalibWFCTA::DirName[3][200]={"/eos/user/y/yangmj/wfcta/decode/led_calibrate/LED_Calibrate_End2End_Factor","/afs/ihep.ac.cn/users/h/hliu/public/WFDataDir","/workfs/ybj/yinlq/LHAASO/WFCTA/WFCTA_Info_Parameter"};

void CalibWFCTA::Init(){
   for(int ii=0;ii<NCTMax;ii++){
      sipmcalib_norm[ii]=0;
      sipmcalib_slope[ii]=0;
      for(int jj=0;jj<MAXPMT;jj++){
         unif_factor[ii][jj]=0;
         deltag_20[ii][jj]=0;
      }
   }

   fin=0;
   tree=0;
   tree2=0;
   entries=-1;
   curentry=-1;
   Reset();
}
void CalibWFCTA::Reset(){
   curTel=0;
   rabbitTime=0;
   for(int ii=0;ii<MAXPMT;ii++){
      Ratio_HL[ii]=22;
      Dead_Channel_Flag[ii]=0;
      mBaseH_RMS[ii]=0;
      mBaseL_RMS[ii]=0;
      mBaseH[ii]=0;
      mBaseL[ii]=0;
      mAdcH_RMS[ii]=0;
      mAdcL_RMS[ii]=0;
      mAdcH[ii]=0;
      mAdcL[ii]=0;
      Low_Flag[ii]=0;
      H_Gain_Factor[ii]=0;
      L_Gain_Factor[ii]=0;
      Gain_Factor[ii]=0;
   }
}
void CalibWFCTA::Release(){
   for(int ii=0;ii<NCTMax;ii++){
      if(sipmcalib_norm[ii]) {delete sipmcalib_norm[ii]; sipmcalib_norm[ii]=0;}
      if(sipmcalib_slope[ii]) {delete sipmcalib_slope[ii]; sipmcalib_slope[ii]=0;}
   }

   if(tree) {delete tree; tree=0;}
   if(tree2) {delete tree2; tree2=0;}
   if(fin) {fin->Close(); fin=0;}
}
int CalibWFCTA::LoadCalibSiPM(int version,char* dirname){
   int ndata=0;
   if(version==2||version<=0){
      char filename[200]="";
      TString ss(dirname?dirname:DirName[1]); if(ss.BeginsWith("/eos/")) strcpy(filename,"root:://eos01.ihep.ac.cn/");
      strcat(filename,Form("%s/temp_corr_v2.root",dirname?dirname:DirName[1]));
      TFile* fin0=TFile::Open(filename,"READ");
      if(!fin0) return 0;
      for(int ii=0;ii<NCTMax;ii++){
         sipmcalib_norm[ii]=(TGraphErrors*)fin0->Get(Form("Norm_Tel%d",ii+1));
         sipmcalib_slope[ii]=(TGraphErrors*)fin0->Get(Form("Slope_Tel%d",ii+1));
         if(sipmcalib_norm[ii]) ndata++;
      }
      fin0->Close();
   }
   if(version==3||version<=0){
      for(int itel=1;itel<=NCTMax;itel++){
         char filename[200]="";
         TString ss(dirname?dirname:DirName[2]); if(ss.BeginsWith("/eos/")) strcpy(filename,"root:://eos01.ihep.ac.cn/");
         strcat(filename,Form("%s/Uniform_Factor/WFCTA_%02d.txt",dirname?dirname:DirName[2],itel));
         char filename2[200]="";
         TString ss2(dirname?dirname:DirName[2]); if(ss2.BeginsWith("/eos/")) strcpy(filename2,"root:://eos01.ihep.ac.cn/");
         strcat(filename2,Form("%s/Tempature_Factor/W%02d_Tempature_Factor.txt",dirname?dirname:DirName[2],itel));
         FILE *fp_gt, *fp_unif;
         fp_unif=fopen(filename,"r");
         fp_gt=fopen(filename2,"r");
         if(fp_unif==NULL) { cerr<<"Can not open WFCTA_"<<itel<<".txt"<<endl;}
         if(fp_gt==NULL) { cerr<<"Can not open W"<<itel<<"_Tempature_Factor.txt"<<endl;}
         if(fp_unif==NULL||fp_gt==NULL) continue;

         //float deltag_20[MAXPMT];
         float temp = 0;
         short F_DB, mChannel, mSiPM;
         for(int i=1;i<9;i++){
          for(int j=1;j<9;j++){
            for(int k=1;k<17;k++){
              F_DB = j*10+i;
              mChannel = k;
              fscanf(fp_gt,"%f\n",&temp);
              SC_Channel2SiPM(F_DB, mChannel, &mSiPM);
              deltag_20[itel-1][mSiPM] = temp/100.;
            }
          }
         }
         fclose(fp_gt);
         ndata++;

         int ich;
         while(!feof(fp_unif))
         {
           fscanf(fp_unif,"%d %f %*lf %*d %*lf\n",&ich,&temp);
           unif_factor[itel-1][ich] = temp;
         }
         fclose(fp_unif);
      }
   }
   return ndata;
}

void CalibWFCTA::SetDirName(int version,char* dirname){
   if(version<1||version>3) return;
   else strcpy(DirName[version-1],dirname);
}
CalibWFCTA* CalibWFCTA::GetHead(){
   if(!_Head) _Head=new CalibWFCTA();
   return _Head;
}
void CalibWFCTA::SetBranchAddress(){
   Reset();
   if(!tree) return;
   if(tree->GetBranchStatus("iTel")) tree->SetBranchAddress("iTel",&curTel);
   if(tree->GetBranchStatus("rabbitTime")) tree->SetBranchAddress("rabbitTime",&rabbitTime);
   if(tree->GetBranchStatus("mBaseH_RMS")) tree->SetBranchAddress("mBaseH_RMS",&mBaseH_RMS);
   if(tree->GetBranchStatus("mBaseL_RMS")) tree->SetBranchAddress("mBaseL_RMS",&mBaseL_RMS);
   if(tree->GetBranchStatus("mBaseH")) tree->SetBranchAddress("mBaseH",&mBaseH);
   if(tree->GetBranchStatus("mBaseL")) tree->SetBranchAddress("mBaseL",&mBaseL);
   if(tree->GetBranchStatus("mAdcH")) tree->SetBranchAddress("mAdcH",&mAdcH);
   if(tree->GetBranchStatus("mAdcL")) tree->SetBranchAddress("mAdcL",&mAdcL);
   if(tree->GetBranchStatus("mAdcL_RMS")) tree->SetBranchAddress("mAdcL_RMS",&mAdcL_RMS);
   if(tree->GetBranchStatus("mAdcH_RMS")) tree->SetBranchAddress("mAdcH_RMS",&mAdcH_RMS);
   if(tree->GetBranchStatus("Low_Flag")) tree->SetBranchAddress("Low_Flag",&Low_Flag);
   if(tree->GetBranchStatus("H_Gain_Factor")) tree->SetBranchAddress("H_Gain_Factor",&H_Gain_Factor);
   if(tree->GetBranchStatus("L_Gain_Factor")) tree->SetBranchAddress("L_Gain_Factor",&L_Gain_Factor);
   if(tree->GetBranchStatus("Gain_Factor")) tree->SetBranchAddress("Gain_Factor",&Gain_Factor);
   if(tree2){
      if(tree2->GetBranchStatus("Ratio_HL")) tree2->SetBranchAddress("Ratio_HL",&Ratio_HL);
      if(tree2->GetBranchStatus("Dead_Channel_Flag")) tree2->SetBranchAddress("Dead_Channel_Flag",&Dead_Channel_Flag);
   }
}
long int CalibWFCTA::LoadDay(int Day,int iTel){
   if(fin){
      char filename0[200]="";
      strcpy(filename0,fin->GetName());
      int size=strlen(filename0);
      int loc0=0;

      for(int ii=0;ii<size;ii++){
         if(filename0[ii]=='/') loc0=ii+1;
      }
      int ncount=0;
      for(int ii=size-1;ii>=loc0;ii--){
         if(filename0[ii]=='_') ncount++;
         if(ncount==2){
            loc0=ii+1;
            break;
         }
      }

      ncount=0;
      int nchar=0;
      char telinfo[10];
      int nchar1=0;
      char dayinfo[10];
      int nchar2=0;
      char year[10];
      char month[10];
      char day[10];
      for(int ii=loc0;ii<size;ii++){
         if(filename0[ii]=='_') {ncount++; continue;}
         if(filename0[ii]=='.') break;
         if(ncount==0){
            dayinfo[nchar2++]=filename0[ii];
         }
         else if(ncount==1){
            telinfo[nchar1++]=filename0[ii];
         }
      }
      telinfo[nchar1++]='\0';
      nchar=0;
      for(int ii=0;ii<4;ii++) year[nchar++]=dayinfo[ii];
      year[nchar++]='\0';
      nchar=0;
      for(int ii=4;ii<6;ii++) month[nchar++]=dayinfo[ii];
      month[nchar++]='\0';
      nchar=0;
      for(int ii=6;ii<nchar2;ii++) day[nchar++]=dayinfo[ii];
      day[nchar++]='\0';
      if(jdebug>0) printf("CalibWFCTA::LoadDay: filename=%s iTel=%d(%d) Day=%04d-%02d-%02d(%d)\n",fin->GetName(),atoi(telinfo),iTel,atoi(year),atoi(month),atoi(day),Day);

      int curday=atoi(year)*10000+atoi(month)*100+atoi(day);
      if(Day==curday&&atoi(telinfo)==iTel) return entries;
      else{ //load another day
         if(tree){delete tree; tree=0;} entries=-1; curentry=-1;
         if(tree2){delete tree2; tree2=0;}
         fin->Close();
         fin=0;
      }
   }

   char filename[200]="";
   TString ss(DirName[0]); if(ss.BeginsWith("/eos/")) strcpy(filename,"root:://eos01.ihep.ac.cn/");
   strcat(filename,Form("%s/LED_Calibrate_Factor_%04d%02d%02d_%02d.root",DirName[0],Day/10000,(Day%10000)/100,(Day%100),iTel));
   TDirectory* gdir=gDirectory;
   fin=TFile::Open(filename,"READ");
   if(fin) fin->cd();
   tree=fin?(TTree*)fin->Get("LED_Signal"):0;
   tree2=fin?(TTree*)fin->Get("ChannelInfo"):0;
   SetBranchAddress();
   entries=tree?(tree->GetEntries()):-1;
   curentry=-1;
   if(entries>0){curentry=0; tree->GetEntry(curentry);}
   if(tree2&&tree2->GetEntries()>0) tree2->GetEntry(0);
   gdir->cd();
   return entries;
}
int CalibWFCTA::IsTimeEqual(int time1,int time2){
   int time_low=time1-15;
   int time_hig=time1+15;
   if(time2>=time_low&&time2<time_hig) return 0;
   else if(time2<time_low) return -1;
   else return 1;
}
int CalibWFCTA::LoadTime(int time){
   if(entries<=0) return -10;
   if(curentry<0){curentry=0; tree->GetEntry(curentry);}
   long int entryi=curentry;
   int timei=rabbitTime;
   if(jdebug>2) printf("CalibWFCTA::LoadTime: test0,entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei,time);
   int timei2=0;
   long entryi2=-1;
   if(IsTimeEqual(timei,time)==0){
      if(jdebug>0) printf("CalibWFCTA::LoadTime: test0,return currrent entry. entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei,time);
      return 1;
   }
   else if(IsTimeEqual(timei,time)>0){
      if(curentry>=entries-1){
         if(jdebug>0) printf("CalibWFCTA::LoadTime: test0,over last entry. currrent entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei,time);
         return -3;
      }
      entryi2=entryi+1;
      curentry=entryi2;
      tree->GetEntry(curentry);
      timei2=rabbitTime;
      if(jdebug>2) printf("CalibWFCTA::LoadTime: test1,entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei2,time);
      if(IsTimeEqual(timei2,time)==0){
         if(jdebug>0) printf("CalibWFCTA::LoadTime: test1,return currrent entry. entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei2,time);
         return 1;
      }
      else if(IsTimeEqual(timei2,time)>0){
         entryi2=entries-1;
         curentry=entryi2;
         tree->GetEntry(curentry);
         timei2=rabbitTime;
         if(jdebug>2) printf("CalibWFCTA::LoadTime: test2,entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei2,time);
         if(IsTimeEqual(timei2,time)==0){
            if(jdebug>0) printf("CalibWFCTA::LoadTime: test2,return currrent entry. entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei2,time);
            return 1;
         }
         else if(IsTimeEqual(timei2,time)>0){
            if(jdebug>0) printf("CalibWFCTA::LoadTime: test2,over last entry. entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei2,time);
            return -3;
         }
      }
   }
   else{
      if(curentry<=0){
         if(jdebug>0) printf("CalibWFCTA::LoadTime: test0,before first entry. currrent entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei,time);
         return -2;
      }
      entryi2=entryi-1;
      curentry=entryi2;
      tree->GetEntry(curentry);
      timei2=rabbitTime;
      if(jdebug>2) printf("CalibWFCTA::LoadTime: test1,entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei2,time);
      if(IsTimeEqual(timei2,time)==0){
         if(jdebug>0) printf("CalibWFCTA::LoadTime: test1,return currrent entry. entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei2,time);
         return 1;
      }
      else if(IsTimeEqual(timei2,time)<0){
         entryi2=0;
         curentry=entryi2;
         tree->GetEntry(curentry);
         timei2=rabbitTime;
         if(jdebug>2) printf("CalibWFCTA::LoadTime: test2,entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei2,time);
         if(IsTimeEqual(timei2,time)==0){
            if(jdebug>0) printf("CalibWFCTA::LoadTime: test2,return currrent entry. entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei2,time);
            return 1;
         }
         else if(IsTimeEqual(timei2,time)<0){
            if(jdebug>0) printf("CalibWFCTA::LoadTime: test2,before first entry. entry=%ld(%ld) Time=%d(%d)\n",curentry,entries,timei2,time);
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
      if(jdebug>1) printf("CalibWFCTA::LoadTime: loop%d time=%d entries=%ld ref1={%d,%ld} ref2={%d,%ld}\n",nloop,time,entries,time1,entry1,time2,entry2);
      long int entry_test=(entry1+entry2)/2;
      curentry=entry_test;
      tree->GetEntry(curentry);
      int time_test=rabbitTime;
      if(jdebug>3) printf("CalibWFCTA::LoadTime: loop%d,entry=%ld(%ld) Time=%d(%d)\n",nloop,curentry,entries,time_test,time);
      if(IsTimeEqual(time_test,time)==0){
         if(jdebug>1) printf("CalibWFCTA::LoadTime: loop%d,return currrent entry. entry=%ld(%ld) Time=%d(%d)\n",nloop,curentry,entries,time_test,time);
         return 1;
      }
      else if(IsTimeEqual(time_test,time)>0){
         time1=time_test;
         entry1=entry_test;
      }
      else{
         time2=time_test;
         entry2=entry_test;
      }
      nloop++;
   }
   if(IsTimeEqual(time1,time)==0||IsTimeEqual(time2,time)==0) return 1;
   else if(IsTimeEqual(time1,time)>0||IsTimeEqual(time2,time)<0) return -1;
   else if(IsTimeEqual(time1,time)<0) return -2;
   else return -3;
}
int CalibWFCTA::LoadFromrabbitTime(double time,int iTel){
   int daytime=(int)(time+0.5);
   int curday=CommonTools::GetTelDay(daytime);
   long int totentry=LoadDay(curday,iTel);
   return LoadTime(daytime);
}
int CalibWFCTA::LoadEntry(long int entry){
   if(entry<0||entry>=entries) return 0;
   if(!tree) return -1;
   curentry=entry;
   tree->GetEntry(curentry);
   int Time=rabbitTime;
   return Time;
}
void CalibWFCTA::Dump(int isipm){
   int Time=rabbitTime;
   int year=CommonTools::TimeFlag(Time,1);
   int month=CommonTools::TimeFlag(Time,2);
   int day=CommonTools::TimeFlag(Time,3);
   int hour=CommonTools::TimeFlag(Time,4);
   int min=CommonTools::TimeFlag(Time,5);
   int sec=CommonTools::TimeFlag(Time,6);
   printf("CalibWFCTA::Dump: FileName=%s tree={%p,%p} entries=%ld cur_entry=%ld iTel=%d cur_time=%d(%04d-%02d-%02d %02d:%02d:%02d)\n",fin?fin->GetName():"",tree,tree2,entries,curentry,curTel,Time,year,month,day,hour,min,sec);
   printf("isipm=%d\n",isipm);
   printf("BaseH={%.1lf %.3lf} BaseL={%.1lf %.3lf}\n",mBaseH[isipm],mBaseH_RMS[isipm],mBaseL[isipm],mBaseL_RMS[isipm]);
   printf("AdcH={%.1lf %.3lf} AdcL={%.1lf %.3lf} Flag=%d\n",mAdcH[isipm],mAdcH_RMS[isipm],mAdcL[isipm],mAdcL_RMS[isipm],Low_Flag[isipm]);
   printf("Gain Factor={%.1lf %.1lf %.1lf}\n",H_Gain_Factor[isipm],L_Gain_Factor[isipm],Gain_Factor[isipm]);
   printf("HL_Ratio=%.3lf Dead_Flag=%d\n",Ratio_HL[isipm],(int)Dead_Channel_Flag[isipm]);
}

double CalibWFCTA::DoCalibSiPM(int iTel,int isipm,double input,double temperature,double Time,int calibtype,int type){
   double res=-1;
   if(iTel<1||iTel>NCTMax) return res;
   if(isipm<0||isipm>1023) return res;
   double temp_ref=20.;
   if(UseSiPMCalibVer==1){
      int CalibType=(calibtype&0x3);
      if(CalibType==0x0) res=input;
      else{
         bool ishig=(bool)(type>10);
         int status=LoadFromrabbitTime((double)Time,iTel);
         if(status>0){
            double corr=ishig?H_Gain_Factor[isipm]:L_Gain_Factor[isipm];
            if(corr<=0) res=ForceCorr?res:input;
            else res=input/corr;
            ///additional changes due to sipm status
            if(Dead_Channel_Flag[isipm]==1) res=0;
            else if(Dead_Channel_Flag[isipm]==2){if(ishig) res=0;}
            else if(Dead_Channel_Flag[isipm]==3){if(!ishig) res=0;}
         }
         else res=ForceCorr?res:input;
      }
   }
   if(UseSiPMCalibVer==2){ //use sipm=529 as reference
      int CalibType=(calibtype&0x3);
      if(!sipmcalib_norm[iTel-1]) res=ForceCorr?res:input;
      else{
      int isipm_ref=529;
      double norm_ref=sipmcalib_norm[iTel-1]->Eval(isipm_ref);
      double slope_ref=sipmcalib_slope[iTel-1]->Eval(isipm_ref);
      double Qref_Tref=norm_ref+slope_ref*temp_ref;
      double Qref_Ti=norm_ref+slope_ref*temperature;
      double norm=sipmcalib_norm[iTel-1]->Eval(isipm);
      double slope=sipmcalib_slope[iTel-1]->Eval(isipm);
      double Qi_Ti=norm+slope*temperature;
      double Qi_Tref=norm+slope*temp_ref;
      //printf("CalibWFCTA::DoCalibSiPM1: iTel=%d isipm=%d temp=%lf, norm=%lf slope=%lf\n",iTel,isipm,temperature,norm,slope);
      if(CalibType==0x0) res=input;
      if(CalibType==0x1) res=input/Qi_Ti*Qref_Ti; //uniform correction
      if(CalibType==0x2) res=input/Qi_Ti*Qi_Tref; //temperature correction
      if(CalibType==0x3) res=input/Qi_Ti*Qref_Tref; //both correction
      if(Time<1573876800) res*=(8.e5/1.1e6);
      }
   }
   else if(UseSiPMCalibVer==3){
      if(unif_factor[iTel-1][isipm]==0||deltag_20[iTel-1][isipm]==0) res=ForceCorr?res:input;
      else{
      int CalibType=(calibtype&0x7);
      int Ncell = 360000;
      if(CalibType==0x0) res=input;
      if(CalibType==0x1) res=input/unif_factor[iTel-1][isipm]; //uniform correction
      if(CalibType==0x2) res=input*(1 + deltag_20[iTel-1][isipm]*(temperature-temp_ref)); //temperature correction
      if(CalibType==0x3) res=input/unif_factor[iTel-1][isipm]*(1 + deltag_20[iTel-1][isipm]*(temperature-temp_ref)); //both correction
      if(CalibType==0x7){
         res=input/unif_factor[iTel-1][isipm]*(1 + deltag_20[iTel-1][isipm]*(temperature-temp_ref));
         res=-Ncell*log(1-res/Ncell);
      }
      //printf("CalibWFCTA::DoCalibSiPM2: iTel=%d isipm=%d temp=%lf, unif=%lf delt=%lf\n",iTel,isipm,temperature,unif_factor[iTel-1][isipm],deltag_20[iTel-1][isipm]);
      if(Time<1573876800) res*=(8.e5/1.1e6);
      }
   }

   //additional specific correction for LED
   if((type%10)==2&&UseSiPMCalibVer>=1&&UseSiPMCalibVer<=3){
      int CalibType=(calibtype&0x8);
      if(CalibType==0x8) res=CalibLED::GetHead()->DoLedDriveTempCorr(res,isipm,(double)Time,iTel);
   }

   return res;
}
