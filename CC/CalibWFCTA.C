#include "CalibWFCTA.h"
//#include "camera.h"
extern void SC_Channel2SiPM(short F_DB, short mChannel, short *mSiPM);
#include "TFile.h"
#include <iostream>
using namespace std;

CalibWFCTA* CalibWFCTA::_Head=0;
int CalibWFCTA::UseSiPMCalibVer=1;

void CalibWFCTA::Init(){
   for(int ii=0;ii<CalibMaxTel;ii++){
      sipmcalib_norm[ii]=0;
      sipmcalib_slope[ii]=0;
      for(int jj=0;jj<1024;jj++){
         unif_factor[ii][jj]=0;
         deltag_20[ii][jj]=0;
      }
   }
}
void CalibWFCTA::Release(){
   for(int ii=0;ii<CalibMaxTel;ii++){
      if(sipmcalib_norm[ii]) {delete sipmcalib_norm[ii]; sipmcalib_norm[ii]=0;}
      if(sipmcalib_slope[ii]) {delete sipmcalib_slope[ii]; sipmcalib_slope[ii]=0;}
   }
   _Head=0;
}
int CalibWFCTA::LoadCalibSiPM(int version,char* dirname){
   int ndata=0;
   if(version==1||version<=0){
      TFile* fin=0;
      if(!dirname){
         fin=TFile::Open("/afs/ihep.ac.cn/users/h/hliu/public/WFDataDir/temp_corr_v2.root","READ");
      }
      else{
         fin=TFile::Open(Form("%s/temp_corr.root",dirname));
      }
      if(!fin) return 0;
      for(int ii=0;ii<CalibMaxTel;ii++){
         sipmcalib_norm[ii]=(TGraphErrors*)fin->Get(Form("Norm_Tel%d",ii+1));
         sipmcalib_slope[ii]=(TGraphErrors*)fin->Get(Form("Slope_Tel%d",ii+1));
         if(sipmcalib_norm[ii]) ndata++;
      }
      fin->Close();
   }
   if(version==2||version<=0){
      for(int itel=1;itel<=CalibMaxTel;itel++){
         FILE *fp_gt, *fp_unif;
         if(!dirname){
            fp_unif=fopen(Form("/workfs/ybj/yinlq/LHAASO/WFCTA/WFCTA_Info_Parameter/Uniform_Factor/WFCTA_%02d.txt",itel),"r");
            fp_gt=fopen(Form("/workfs/ybj/yinlq/LHAASO/WFCTA/WFCTA_Info_Parameter/Tempature_Factor/W%02d_Tempature_Factor.txt",itel),"r");
            if(fp_unif==NULL) { cerr<<"Can not open WFCTA_"<<itel<<".txt"<<endl;}
            if(fp_gt==NULL) { cerr<<"Can not open W"<<itel<<"_Tempature_Factor.txt"<<endl;}
            if(fp_unif==NULL||fp_gt==NULL) continue;
         }
         else{
            fp_unif=fopen(Form("%s/Uniform_Factor/WFCTA_%02d.txt",dirname,itel),"r");
            fp_gt=fopen(Form("%s/Tempature_Factor/W%02d_Tempature_Factor.txt",dirname,itel),"r");
            if(fp_unif==NULL) { cerr<<"Can not open WFCTA_"<<itel<<".txt"<<endl;}
            if(fp_gt==NULL) { cerr<<"Can not open W"<<itel<<"_Tempature_Factor.txt"<<endl;}
            if(fp_unif==NULL||fp_gt==NULL) continue;
         }

         //float deltag_20[1024];
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

CalibWFCTA* CalibWFCTA::GetHead(char* dirname){
   if(_Head) return _Head;
   else{
      _Head=new CalibWFCTA(dirname);
   }
   return _Head;
}
double CalibWFCTA::DoCalibSiPM(int iTel,int isipm,double input,double temperature,int calibtype){
   if(iTel<1||iTel>CalibMaxTel) return -1.;
   if(isipm<0||isipm>1023) return -1;
   double temp_ref=20.;
   if(UseSiPMCalibVer==1){ //use sipm=529 as reference
      int CalibType=(calibtype&0x3);
      if(!sipmcalib_norm[iTel-1]) return input;
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
      double res=-1;
      if(CalibType==0x0) res=input;
      if(CalibType==0x1) res=input/Qi_Ti*Qref_Ti; //uniform correction
      if(CalibType==0x2) res=input/Qi_Ti*Qi_Tref; //temperature correction
      if(CalibType==0x3) res=input/Qi_Ti*Qref_Tref; //both correction
      return res;
   }
   else if(UseSiPMCalibVer==2){
      int CalibType=(calibtype&0x7);
      int Ncell = 360000;
      double res=-1;
      if(CalibType==0x0) res=input;
      if(CalibType==0x1) res=input/unif_factor[iTel-1][isipm]; //uniform correction
      if(CalibType==0x2) res=input*(1 + deltag_20[iTel-1][isipm]*(temperature-temp_ref)); //temperature correction
      if(CalibType==0x3) res=input/unif_factor[iTel-1][isipm]*(1 + deltag_20[iTel-1][isipm]*(temperature-temp_ref)); //both correction
      if(CalibType==0x7){
         res=input/unif_factor[iTel-1][isipm]*(1 + deltag_20[iTel-1][isipm]*(temperature-temp_ref));
         res=-Ncell*log(1-res/Ncell);
      }
      //printf("CalibWFCTA::DoCalibSiPM2: iTel=%d isipm=%d temp=%lf, unif=%lf delt=%lf\n",iTel,isipm,temperature,unif_factor[iTel-1][isipm],deltag_20[iTel-1][isipm]);
      return res;
   }
   return -1;
}
