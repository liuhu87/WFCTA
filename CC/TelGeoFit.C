#include "TelGeoFit.h"
#include "Readparam.h"
#include "TMath.h"
#include "WFCTAEvent.h"
#include "RotateDB.h"
int TelGeoFit::jdebug=0;
void TelGeoFit::Init(){
   nevent=0;
   for(int ii=0;ii<NCTMax;ii++){
      telpos[ii][0]=telpos[ii][1]=0;
      telpos[ii][2]=-1.0e-8;
      teldir0[ii][0]=teldir0[ii][1]=0;
   }
   for(int ii=0;ii<MAXEvent;ii++){
      telindex[ii]=0;
      rabbitTime[ii]=0;
      iEvent[ii]=0;
      rotdir0[ii][0]=rotdir0[ii][1]=0;
      for(int jj=0;jj<MAXPMT;jj++) Npe_sipm[ii][jj]=0;
      Npe_sum[ii]=0;
   }
   minimizer=0;
}
void TelGeoFit::Clear(){
   if(minimizer) {delete minimizer; minimizer=0;}
}
void TelGeoFit::SetInfoInit(char* filename){
   char name_buff[200]="";
   if(!filename) strcpy(name_buff,"/afs/ihep.ac.cn/users/h/hliu/public/WFDataDir/default.inp");
   else strcpy(name_buff,filename);
   WReadConfig read;
   WReadConfig* readconfig=&read;
   readconfig->readparam(name_buff);
   double lhaaso_coo[3]={readconfig->GetLHAASOCoo(0),readconfig->GetLHAASOCoo(1),readconfig->GetLHAASOCoo(2)};
   int nict=0;
   for(int itel=1;itel<=NCTMax;itel++){
      if(nict>=readconfig-> GetCTNumber()) break;
      double flag=readconfig->GetCTPosition(itel,3);
      //if(jdebug>0) printf("TelGeoFit::SetInfoInit: iTel=%d flag=%.1lf\n",itel,flag);
      if(flag<-0.5) continue;
      telpos[itel-1][0]=readconfig->GetCTPosition(itel,0);
      telpos[itel-1][1]=readconfig->GetCTPosition(itel,1);
      telpos[itel-1][2]=readconfig->GetCTPosition(itel,2);
      teldir0[itel-1][0]=readconfig->GetCTPosition(itel,3) * TMath::DegToRad();
      teldir0[itel-1][1]=readconfig->GetCTPosition(itel,4) * TMath::DegToRad();
      nict++;
   }
   rotpos[0]=readconfig->GetLaserCoo(0)-lhaaso_coo[0];
   rotpos[1]=readconfig->GetLaserCoo(1)-lhaaso_coo[1];
   rotpos[2]=readconfig->GetLaserCoo(2)-lhaaso_coo[2];
}
void TelGeoFit::Dump(){
   printf("Tel and Rot Info:\n");
   for(int itel=0;itel<NCTMax;itel++){
      if(telpos[itel][2]<-9.e-9) continue;
      printf("iTel=%d: pos={%.1lf,%.1lf,%.1lf} ele=%.2lf azi=%.2lf\n",itel+1,telpos[itel][0],telpos[itel][1],telpos[itel][2],teldir0[itel][0]/PI*180,teldir0[itel][1]/PI*180);
   }
   printf("Rot: pos={%.1lf,%.1lf,%.1lf}\n",rotpos[0],rotpos[1],rotpos[2]);
   printf("Event Info: nevent=%d\n",nevent);
   for(int ievt=0;ievt<nevent;ievt++){
      printf("ievt=%d iTel=%d iEvent=%d rabbitTime=%d+%.8lf npe_sum=%.1lf rot_ele=%.2lf rot_azi=%.2lf\n",ievt,telindex[ievt],iEvent[ievt],(int)(rabbitTime[ievt]),rabbitTime[ievt]-(int)(rabbitTime[ievt]),Npe_sum[ievt],rotdir0[ievt][0]/PI*180,rotdir0[ievt][1]/PI*180);
   }
   printf("\n");
}
int TelGeoFit::AddEvent(WFCTAEvent* pev,int type,bool doclean){
   if(!pev) return nevent;
   //check wheather it is laser event
   //int index=(RotateDB::GetHead()->LaserIsFine(pev));
   int index=(RotateDB::GetHead()->GetEleAzi(pev));
   //if(index<=0) { printf("Not Laser Event. index=%d\n",index); return nevent;}
   telindex[nevent]=pev->iTel;
   rabbitTime[nevent]=pev->rabbitTime+pev->rabbittime*2.0e-8;
   iEvent[nevent]=pev->iEvent;
   rotdir0[nevent][0]=RotateDB::GetHead()->GetElevation()/180*PI;
   rotdir0[nevent][1]=RotateDB::GetHead()->GetAzimuth()/180*PI;
   int size=pev->iSiPM.size();
   Npe_sum[nevent]=0;
   for(int ii=0;ii<size;ii++){
      if(doclean&&(!pev->CleanImage(ii,pev->iTel,true,true))) continue;
      int sipm=pev->iSiPM.at(ii);
      Npe_sipm[nevent][sipm]=pev->GetContent(ii,pev->iTel,type,true,false);
      Npe_sum[nevent]+=Npe_sipm[nevent][sipm];
   }
   nevent++;
   return nevent;
}
bool TelGeoFit::GetImageCoo(double zenith,double azimuth,double dir_in[3],double &xx,double &yy,bool IsLocal){
   double norm=sqrt(dir_in[0]*dir_in[0]+dir_in[1]*dir_in[1]+dir_in[2]*dir_in[2]);
   double theta=IsLocal?(PI/2):(PI/2-zenith);
   double phi=IsLocal?0:azimuth;
   double zz0=(dir_in[0]*cos(phi)+dir_in[1]*sin(phi))/norm*cos(theta)+dir_in[2]/norm*sin(theta);
   if(zz0>=0) return false;
   double xx0=(dir_in[0]*cos(phi)+dir_in[1]*sin(phi))/norm*sin(theta)-dir_in[2]/norm*cos(theta);
   double yy0=(-dir_in[0]*sin(phi)+dir_in[1]*cos(phi))/norm;
   xx=-yy0/zz0/PI*180;
   yy=-xx0/zz0/PI*180;
   return true;
}
void TelGeoFit::GetOutDir(double zenith,double azimuth,double imagexy[2],double dir_out[3],bool IsLocal){
   double theta=IsLocal?(PI/2):(PI/2-zenith);
   double phi=IsLocal?0:azimuth;
   double xx0=imagexy[1]/180.*PI;
   double yy0=imagexy[0]/180.*PI;
   dir_out[0]=(-xx0*sin(theta)+cos(theta))*cos(phi)+yy0*sin(phi);
   dir_out[1]=(-xx0*sin(theta)+cos(theta))*sin(phi)-yy0*cos(phi);
   dir_out[2]=xx0*cos(theta)+sin(theta);
   double norm=sqrt(dir_out[0]*dir_out[0]+dir_out[1]*dir_out[1]+dir_out[2]*dir_out[2]);
   for(int ii=0;ii<3;ii++) dir_out[ii]/=norm;
}
double TelGeoFit::Interface(const double* par){
   double teldir_fit[NCTMax][2];
   for(int itel=0;itel<NCTMax;itel++){
      teldir_fit[itel][0]=teldir_fit[itel][1]=-1000;
   }
   double rotele=par[0];
   double rotazi=par[1];
   double rotrefele_bias=par[2];
   double rotrefazi_bias=par[3];
   int ntel=(int)(par[4]+0.5);
   int telstart=5;
   for(int itel=0;itel<ntel;itel++){
      int iTel=(int)(par[itel*3+0+telstart]+0.5);
      teldir_fit[iTel-1][0]=par[itel*3+1+telstart];
      teldir_fit[iTel-1][1]=par[itel*3+2+telstart];
      if(jdebug>1) printf("TelGeoFit::Interface: pars: iTel=%d dir={%.2lf,%.2lf}\n",iTel,teldir_fit[iTel-1][0]/PI*180,teldir_fit[iTel-1][1]/PI*180);
   }
   double Chi2=0;
   int Ndof=0;
   for(int ievt=0;ievt<nevent;ievt++){
      double chi2=0;
      int ndof=0;
      int itel=telindex[ievt]-1;
      if(itel<0||itel>=NCTMax) continue;
      if(teldir_fit[itel][0]<-500) continue;
      double rotdir_fit[2]={rotdir0[ievt][0]+rotrefele_bias,rotdir0[ievt][1]+rotrefazi_bias};
      WFCTAEvent::CalDir_out(rotele,rotazi,rotrefele_bias,rotrefazi_bias,rotdir0[ievt][0],rotdir0[ievt][1],rotdir_fit[0],rotdir_fit[1]);
      if(jdebug>1) printf("TelGeoFit::Interface: ievt=%d ele_in=%.2lf azi_in=%.2lf ele_out=%.2lf azi_out=%.2lf\n",ievt,rotdir0[ievt][0]/PI*180,rotdir0[ievt][1]/PI*180,rotdir_fit[0]/PI*180,rotdir_fit[1]/PI*180);
      double CC,phi;
      WFCTAEvent::GetCCphi(teldir_fit[itel][0],teldir_fit[itel][1],rotpos,rotdir_fit,CC,phi);
      if(jdebug>1) printf("TelGeoFit::Interface: ievt=%d(%d) itel=%d ele=%.2lf azi=%.2lf cc=%.2lf phi=%.2lf\n",ievt,nevent,itel+1,teldir_fit[itel][0]/PI*180,teldir_fit[itel][1]/PI*180,CC/PI*180,phi/PI*180);
      for(int isipm=0;isipm<MAXPMT;isipm++){
         if(Npe_sipm[ievt][isipm]<=0) continue;
         double ImageX,ImageY;
         double dcell=WFCTAEvent::GetImageXYCoo(isipm,ImageX,ImageY,-1,false);
         if(dcell<0) continue;
         double distance=(ImageX*sin(phi)-ImageY*cos(phi)+CC);
         double probi=Npe_sipm[ievt][isipm]/Npe_sum[ievt];
         chi2+=probi*pow(distance/(dcell/2.),2);
         ndof++;
         if(jdebug>2) printf("TelGeoFit::Interface: ievt=%d sipm=%d chi2=%.1lf ndof=%d\n",ievt,isipm,chi2,ndof);
      }
      if(jdebug>1) printf("TelGeoFit::Interface: ievt=%d chi2=%.1lf ndof=%d\n",ievt,chi2,ndof);
      Chi2+=chi2;
      Ndof+=ndof;
   }
   if(jdebug>0) printf("TelGeoFit::Interface: Chi2=%.1lf NDof=%d\n",Chi2,Ndof);
   return Chi2;
}
bool TelGeoFit::DoFit(int ntel,int* tellist,bool force){
   if(!tellist) return false;
   for(int itel=0;itel<ntel;itel++){
      if(tellist[itel]<1||tellist[itel]>NCTMax) return false;
   }
   if(!minimizer) minimizer=ROOT::Math::Factory::CreateMinimizer("Minuit","Migrad");
   minimizer->Clear();
   minimizer->SetMaxFunctionCalls(1000000);
   minimizer->SetMaxIterations(100000);
   minimizer->SetTolerance(0.001);
   minimizer->SetPrintLevel(0);
   #if defined(__CINT__)
   ROOT::Math::Functor f(this,"TelGeoFit","Interface");
   #else
   ROOT::Math::Functor f(this,&TelGeoFit::Interface,ntel*3+1+4);
   #endif
   minimizer->SetFunction(f);
   //minimizer->SetFixedVariable(0,"RotEle",PI/2);
   //minimizer->SetFixedVariable(1,"RotAzi",0);
   //minimizer->SetFixedVariable(2,"RotRefEle_bias",0);
   //minimizer->SetFixedVariable(3,"RotRefAzi_bias",0);
   double rotele_margin=5./180*PI;
   double rotazi_margin=180./180*PI;
   double rotrefele_margin=5./180*PI;
   double rotrefazi_margin=5./180*PI;
   minimizer->SetLimitedVariable(0,"RotEle",PI/2,fabs(0.1/180.*PI),TMath::Max(0.,PI/2-rotele_margin),TMath::Min(PI/2.,PI/2+rotele_margin));
   minimizer->SetLimitedVariable(1,"RotAzi",0,fabs(5./180.*PI),0.-rotazi_margin,0.+rotazi_margin);
   minimizer->SetLimitedVariable(2,"RotRefEle_bias",0,fabs(0.1/180.*PI),-rotrefele_margin,rotrefele_margin);
   minimizer->SetLimitedVariable(3,"RotRefAzi_bias",0,fabs(0.1/180.*PI),-rotrefazi_margin,rotrefazi_margin);
   minimizer->SetFixedVariable(4,"NTel",ntel);
   int telstart=5;
   double telele_margin=5./180*PI;
   double telazi_margin=5./180*PI;
   for(int itel=0;itel<ntel;itel++){
      int iTel=tellist[itel];
      minimizer->SetFixedVariable(itel*3+0+telstart,Form("Tel%d",iTel),iTel*1.);
      if(iTel==4||iTel==5){
      minimizer->SetLimitedVariable(itel*3+1+telstart,Form("Ele%d",iTel),teldir0[iTel-1][0],fabs(0.1/180.*PI),TMath::Max(0.,teldir0[iTel-1][0]-telele_margin),TMath::Min(PI/2,teldir0[iTel-1][0]+telele_margin));
      minimizer->SetLimitedVariable(itel*3+2+telstart,Form("Azi%d",iTel),teldir0[iTel-1][1],fabs(0.1/180.*PI),teldir0[iTel-1][1]-telazi_margin,teldir0[iTel-1][1]+telazi_margin);
      }
      else{
      minimizer->SetFixedVariable(itel*3+1+telstart,Form("Ele%d",iTel),teldir0[iTel-1][0]);
      minimizer->SetFixedVariable(itel*3+2+telstart,Form("Azi%d",iTel),teldir0[iTel-1][1]);
      }
      if(jdebug>1) printf("TelGeoFit::DoFit: iTel=%d Ele(%d)=%.2lf Azi(%d)=%.2lf\n",iTel,itel*3+1+telstart,teldir0[iTel-1][0]/PI*180,itel*3+2+telstart,teldir0[iTel-1][1]/PI*180);
   }
   minimizer->Minimize();
   minimizer->Hesse();
   return true;
}
//void TelGeoFit::PrintFit(char* filename){
//   for(int ievt=0;ievt<nevent;ievt++){
//      
//   }
//}
