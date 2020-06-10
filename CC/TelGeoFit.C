#include "TelGeoFit.h"
#include "Readparam.h"
#include "TMath.h"
#include "WFCTAEvent.h"
#include "RotateDB.h"
#include "Laser.h"
#include "slalib.h"
#include "WFCone.h"
int TelGeoFit::jdebug=0;
int TelGeoFit::FixRotDir[4]={0,0,0,0};
int TelGeoFit::FixRotPos[3]={0xFF,0xFF,0xFF};
int TelGeoFit::FixRotTime=0xF;
int TelGeoFit::FixTelDir[2]={0,0xFFFFF};
int TelGeoFit::FixTime=-1;
int TelGeoFit::FixLi=-1;
int TelGeoFit::FitLineWidth=1.5;
double TelGeoFit::InfNumber=-1000;
bool TelGeoFit::GoDown=false;
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
      image_pars[ii][0][0]=image_pars[ii][0][1]=0;
      image_pars[ii][1][0]=image_pars[ii][1][1]=0;
      for(int jj=0;jj<MAXPMT;jj++){
         Npe_sipm[ii][jj]=0;
         Time_sipm[ii][jj]=0;
      }
      Npe_sum[ii]=0;
      replaced[ii]=false;
   }
   minimizer=0;
   for(int ii=0;ii<MAXEvent;ii++){
      XYchi2_evt[ii]=0;
      XYndof_evt[ii]=0;
      Tchi2_evt[ii]=0;
      Tndof_evt[ii]=0;
   }
   Chi2_All=0;
   Ndof_All=0;
   for(int ii=0;ii<NCTMax;ii++){
      TelZen_fit[ii]=-1;
      TelAzi_fit[ii]=0;
   }
   for(int ii=0;ii<NCTMax;ii++){
      TelZen_int[ii]=-1;
      TelAzi_int[ii]=0;
   }
   for(int ii=0;ii<NRotMax;ii++){
      for(int jj=0;jj<3;jj++) rotpos[ii][jj]=0;
      rotpos[ii][2]=-1.0e15;
      RotEle_fit[ii]=PI/2;
      RotAzi_fit[ii]=0;
      RefEle_fit[ii]=0;
      RefAzi_fit[ii]=0;
      for(int jj=0;jj<3;jj++) RotPos_fit[ii][jj]=0;
      RotTime_fit[ii]=0;
      RotEle_int[ii]=PI/2;
      RotAzi_int[ii]=0;
      RefEle_int[ii]=0;
      RefAzi_int[ii]=0;
      for(int jj=0;jj<3;jj++) RotPos_int[ii][jj]=0;
      RotTime_int[ii]=0;
   }
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
   int nrot=0;
   for(int iRot=1;iRot<=NRotMax;iRot++){
      double flag1=readconfig->GetLaserCoo(iRot,2);
      double flag2=readconfig->GetLaserDir(iRot,0);
      //if(jdebug>0) printf("TelGeoFit::SetInfoInit: iRot=%d flag={%.1le,%.1lf}\n",iRot,flag1,flag2);
      if(flag1<-1.e10&&flag2<0) continue;
      rotpos[iRot-1][0]=readconfig->GetLaserCoo(iRot,0)-lhaaso_coo[0];
      rotpos[iRot-1][1]=readconfig->GetLaserCoo(iRot,1)-lhaaso_coo[1];
      rotpos[iRot-1][2]=readconfig->GetLaserCoo(iRot,2)-lhaaso_coo[2];
      nrot++;
   }
}
void TelGeoFit::Dump(){
   printf("Tel and Rot Info:\n");
   for(int itel=0;itel<NCTMax;itel++){
      if(telpos[itel][2]<-9.e-9) continue;
      printf("iTel=%d: pos={%.1lf,%.1lf,%.1lf} ele=%.2lf azi=%.2lf\n",itel+1,telpos[itel][0],telpos[itel][1],telpos[itel][2],teldir0[itel][0]/PI*180,teldir0[itel][1]/PI*180);
   }
   for(int irot=0;irot<NRotMax;irot++){
      if(rotpos[irot][2]<-1.e10) continue;
      printf("iRot=%d: pos={%.1lf,%.1lf,%.1lf}\n",irot+1,rotpos[irot][0],rotpos[irot][1],rotpos[irot][2]);
   }
   printf("\n");
   printf("Event Info: nevent=%d\n",nevent);
   for(int ievt=0;ievt<nevent;ievt++){
      double npe_max=0;
      for(int ii=0;ii<MAXPMT;ii++){
         if(Npe_sipm[ievt][ii]>npe_max) npe_max=Npe_sipm[ievt][ii];
      }
      printf("ievt=%d iTel=%d iEvent=%d rabbitTime=%d+%.8lf npe_sum=%.1lf npe_max=%.1lf rot_ele=%.2lf rot_azi=%.2lf replaced=%d\n",ievt,telindex[ievt],iEvent[ievt],(int)(rabbitTime[ievt]),rabbitTime[ievt]-(int)(rabbitTime[ievt]),Npe_sum[ievt],npe_max,90-rotdir0[ievt][0]/PI*180,rotdir0[ievt][1]/PI*180,replaced[ievt]);
      int isipm=30*16;
      printf("isipm=%d: npe_sipm=%.1lf time_sipm=%.2lf\n",isipm,Npe_sipm[ievt][isipm],Time_sipm[ievt][isipm]);
      printf("Image_Pars: CC=%.2lf+-%.3lf phi=%.2lf+-%.3lf\n",image_pars[ievt][0][0]/PI*180,image_pars[ievt][0][1]/PI*180,image_pars[ievt][1][0]/PI*180,image_pars[ievt][1][1]/PI*180);
   }
   printf("\n\n");
}
int TelGeoFit::AddEvent(WFCTAEvent* pev,int type,bool doclean){
   if(!pev) return nevent;
   if(nevent>=MAXEvent) return nevent;
   //check wheather it is laser event
   //int index=(RotateDB::GetHead()->LaserIsFine(pev));
   //int index=(RotateDB::GetHead()->GetEleAzi(pev));
   //if(index<=0) { printf("Not Laser Event. index=%d\n",index); return nevent;}

   //ignore some events
   int irot=RotateDB::GetHead()->GetLi(pev->rabbittime);
   if(irot<0) return nevent;
   if(irot>=0&&RotateDB::rotindex[irot]==2){
      if(pev->rabbitTime==1579981841&&pev->iEvent==52750) return nevent;
      if(pev->rabbitTime==1579966342&&pev->iEvent==50075) return nevent;
   }

   telindex[nevent]=pev->iTel;
   rabbitTime[nevent]=pev->rabbitTime+pev->rabbittime*2.0e-8;
   iEvent[nevent]=pev->iEvent;
   rotdir0[nevent][0]=PI/2-RotateDB::GetHead()->GetElevation()/180*PI;
   rotdir0[nevent][1]=RotateDB::GetHead()->GetAzimuth()/180*PI;
   int size=pev->iSiPM.size();
   Npe_sum[nevent]=0;
   for(int ii=0;ii<size;ii++){
      if(doclean&&(!pev->CleanImage(ii,pev->iTel,true,true))) continue;
      int sipm=pev->iSiPM.at(ii);
      Npe_sipm[nevent][sipm]=pev->GetContent(ii,pev->iTel,type,true,false);
      if(Npe_sipm[nevent][sipm]>30000) Npe_sipm[nevent][sipm]=0;
      Npe_sum[nevent]+=Npe_sipm[nevent][sipm];
      Time_sipm[nevent][sipm]=pev->GetContent(ii,pev->iTel,15,true,false);
   }
   //fitting parameters
   pev->DoFit(0,type,true);
   if(pev->minimizer){
      image_pars[nevent][0][0]=pev->minimizer->X()[2]; //CC
      image_pars[nevent][0][1]=pev->minimizer->Errors()[2]; //eCC
      image_pars[nevent][1][0]=pev->minimizer->X()[3]; //phi
      image_pars[nevent][1][1]=pev->minimizer->Errors()[3]; //ephi
      if(jdebug>2) printf("TelGeoFit::AddEvent: ievt=%d CC={%6.2lf+-%6.3lf} phi={%6.2lf+-%6.3lf}",nevent,image_pars[nevent][0][0]/PI*180,image_pars[nevent][0][1]/PI*180,image_pars[nevent][1][0]/PI*180,image_pars[nevent][1][1]/PI*180);
   }
   nevent++;
   return nevent;
}

void TelGeoFit::GetCrossCoor(double x,double y,double CC,double phi,double &x0,double &y0){
   x0=cos(phi)*(x*cos(phi)+y*sin(phi))-CC*sin(phi);
   y0=sin(phi)*(x*cos(phi)+y*sin(phi))+CC*cos(phi);
}
double TelGeoFit::GetImageXYCoo(int isipm,double &ImageX,double &ImageY,double focus,bool Isdegree){
   if(isipm<0||isipm>=NSIPM) return -1;
   if(WCamera::SiPMMAP[0][0]==0) WCamera::SetSiPMMAP();
   ImageX=-WCamera::GetSiPMX(isipm);
   ImageY=WCamera::GetSiPMY(isipm);
   double numcon=Isdegree?(180./PI):1.;
   double result=0;
   if(focus>0){
      ImageX*=numcon/focus;
      ImageY*=numcon/focus;
      result=SquareCone::D_ConeOut/focus*numcon;
   }
   else{
      ImageX*=numcon/WFTelescope::FOCUS;
      ImageY*=numcon/WFTelescope::FOCUS;
      result=SquareCone::D_ConeOut/WFTelescope::FOCUS*numcon;
   }
   //if(iTel==5) {ImageX*=-1; ImageY*=-1;}
   return result;
}
int TelGeoFit::IsInside(double xx,double yy){
   double ImageX0,ImageY0;
   double dcell=GetImageXYCoo(PIX/2,ImageX0,ImageY0,-1,false);
   int IndexY=-1,IndexX=-1;
   double xmin=100,xmax=-100;
   double ymin=100,ymax=-100;
   double mindistx=1000,mindisty=1000;
   int indexminx=-1,indexminy=-1;
   for(int Yi=0;Yi<PIX;Yi++){
      double ImageX,ImageY;
      GetImageXYCoo(PIX/2+Yi*PIX,ImageX,ImageY,-1,false);
      double y1=ImageY-dcell/2;
      double y2=ImageY+dcell/2;
      if(yy>=y1&&yy<=y2){
         IndexY=Yi;
         break;
      }
      if(y1<ymin) ymin=y1;
      if(y2>ymax) ymax=y2;
      if(fabs(yy-ImageY)<mindisty){
         mindisty=fabs(yy-ImageY);
         indexminy=Yi;
      }
   }
   if(IndexY<0){
      //printf("ymin=%.2lf ymax=%.2lf yy=%.2lf mindisty=%d(%.2lf)\n",ymin,ymax,yy,indexminy,mindisty);
      if((yy>ymax||yy<ymin)||indexminy<0) return -10;
      else{
         for(int Xi=0;Xi<PIX;Xi++){
            double ImageX,ImageY;
            GetImageXYCoo(Xi+indexminy*PIX,ImageX,ImageY,-1,false);
            double x1=ImageX-dcell/2;
            double x2=ImageX+dcell/2;
            if(xx>=x1&&xx<=x2){
               IndexX=Xi;
               break;
            }
            if(x1<xmin) xmin=x1;
            if(x2>xmax) xmax=x2;
            if(fabs(xx-ImageX)<mindistx){
               mindistx=fabs(xx-ImageX);
               indexminx=Xi;
            }
         }
         //printf("xmin=%.2lf xmax=%.2lf xx=%.2lf mindistx=%d(%.2lf)\n",xmin,xmax,xx,indexminx,mindistx);
         if((xx>xmax||xx<xmin)||indexminx<0) return -10;
         else return -1;
      }
   }
   for(int Xi=0;Xi<PIX;Xi++){
      double ImageX,ImageY;
      GetImageXYCoo(Xi+IndexY*PIX,ImageX,ImageY,-1,false);
      double x1=ImageX-dcell/2;
      double x2=ImageX+dcell/2;
      if(xx>=x1&&xx<=x2){
         IndexX=Xi;
         break;
      }
      if(x1<xmin) xmin=x1;
      if(x2>xmax) xmax=x2;
      if(fabs(xx-ImageX)<mindistx){
         mindistx=fabs(xx-ImageX);
         indexminx=Xi;
      }
   }
   if(IndexX<0){
      //printf("xmin=%.2lf xmax=%.2lf xx=%.2lf mindistx=%d(%.2lf)\n",xmin,xmax,xx,indexminx,mindistx);
      if((xx>xmax||xx<xmin)||indexminx<0) return -10;
      else return -1;
   }
   else return IndexY*PIX+IndexX;
}
bool TelGeoFit::GetRange(double CC,double phi,double XY1[4],double XY2[4]){
   if(jdebug>0) printf("TelGeoFit::GetRange:input,CC=%lf phi=%lf\n",CC/PI*180,phi/PI*180);
   double ImageX0,ImageY0;
   GetImageXYCoo(0+0*PIX,ImageX0,ImageY0,-1,false);
   double ImageXMax,ImageYMax;
   GetImageXYCoo(PIX-1+0*PIX,ImageXMax,ImageYMax,-1,false);
   bool xsign0=(ImageX0>ImageXMax);
   double margin=1.0e-5;
   if(fabs(phi)<(0.1/180.*PI)){
      double yy=CC/cos(phi);
      int IndexY=-1;
      for(int Yi=0;Yi<PIX;Yi++){
         double ImageX,ImageY;
         double dcell=GetImageXYCoo(PIX/2+Yi*PIX,ImageX,ImageY,-1,false);
         double y1=ImageY-dcell/2;
         double y2=ImageY+dcell/2;
         //printf("Yi=%d y12={%lf,%lf}\n",Yi,y1,y2);
         if(yy>=y1&&yy<=y2){
            IndexY=Yi;
            break;
         }
      }
      //printf("yy=%lf IndexY=%d\n",yy,IndexY);
      if(IndexY<0) return false;
      XY1[1]=yy;
      XY2[1]=yy;
      double dcell=GetImageXYCoo(0+IndexY*PIX,XY1[0],yy,-1,false);
      XY1[0]+=(xsign0?1:-1)*dcell/2.;
      XY1[2]=XY1[0]*cos(phi)+XY1[1]*sin(phi);
      XY1[3]=1;
      dcell=GetImageXYCoo((PIX-1)+IndexY*PIX,XY2[0],yy,-1,false);
      XY2[0]+=(xsign0?-1:+1)*dcell/2.;
      XY2[2]=XY2[0]*cos(phi)+XY2[1]*sin(phi);
      XY2[3]=1;
      //printf("xsign=%d XY1=%.5lf XY2=%.5lf\n",xsign0,XY1[0],XY2[0]);
   }
   else{
      int exist=0;
      int index=0;
      while(index>=0&&index<PIX){
         double ImageX,ImageY;
         double dcell=GetImageXYCoo(PIX/2+index*PIX,ImageX,ImageY,-1,false);
         double y1=ImageY+dcell/2;
         double y2=ImageY-dcell/2;
         double x1=(y1*cos(phi)-CC)/sin(phi);
         double x2=(y2*cos(phi)-CC)/sin(phi);
         bool inside1=IsInside(x1,y1-margin)>=-1;
         bool inside2=IsInside(x2,y2+margin)>=-1;
         if(jdebug>1) printf("TelGeoFit::GetRange: index=%d xy1={%lf,%lf} xy2={%lf,%lf} inside={%d,%d}\n",index,x1,y1,x2,y2,inside1,inside2);
         if(inside1){
            XY1[1]=y1;
            XY1[0]=x1;
            XY1[2]=XY1[0]*cos(phi)+XY1[1]*sin(phi);
            XY1[3]=1;
            exist=1;
            break;
         }
         else if(inside2){
            int xindex=((x1*ImageX0)>=0)?0:(PIX-1);
            dcell=GetImageXYCoo(xindex+index*PIX,ImageX,ImageY,-1,false);
            bool xsign=(x1>ImageX)?true:false;
            XY1[0]=ImageX+(xsign?1:-1)*dcell/2;
            XY1[1]=(XY1[0]*sin(phi)+CC)/cos(phi);
            XY1[2]=XY1[0]*cos(phi)+XY1[1]*sin(phi);
            XY1[3]=1;
            exist=1;
            break;
         }
         else if(x1*x2<0){
            dcell=GetImageXYCoo(0+index*PIX,ImageX,ImageY,-1,false);
            XY1[0]=ImageX+(xsign0?1:-1)*dcell/2;
            XY1[1]=(XY1[0]*sin(phi)+CC)/cos(phi);
            XY1[2]=XY1[0]*cos(phi)+XY1[1]*sin(phi);
            XY1[3]=1;
            dcell=GetImageXYCoo(PIX-1+index*PIX,ImageX,ImageY,-1,false);
            XY2[0]=ImageX+(xsign0?-1:1)*dcell/2;
            XY2[1]=(XY2[0]*sin(phi)+CC)/cos(phi);
            XY2[2]=XY2[0]*cos(phi)+XY2[1]*sin(phi);
            XY2[3]=1;
            if(XY1[2]>XY2[2]){
               double nbuff[4]={XY2[0],XY2[1],XY2[2],XY2[3]};
               for(int ii=0;ii<4;ii++) XY2[ii]=XY1[ii];
               for(int ii=0;ii<4;ii++) XY1[ii]=nbuff[ii];
            }
            exist=2;
            break;
         }
         index++;
      }
      if(exist<=0) return false;
      else if(exist>=2) return true;
      index=PIX-1;
      while(index>=0&&index<PIX){
         double ImageX,ImageY;
         double dcell=GetImageXYCoo(PIX/2+index*PIX,ImageX,ImageY,-1,false);
         double y1=ImageY-dcell/2;
         double y2=ImageY+dcell/2;
         double x1=(y1*cos(phi)-CC)/sin(phi);
         double x2=(y2*cos(phi)-CC)/sin(phi);
         bool inside1=IsInside(x1,y1+margin)>=-1;
         bool inside2=IsInside(x2,y2-margin)>=-1;
         if(inside1){
            XY2[1]=y1;
            XY2[0]=x1;
            XY2[2]=XY2[0]*cos(phi)+XY2[1]*sin(phi);
            XY2[3]=1;
            exist=2;
            break;
         }
         else if(inside2){
            int xindex=(x1*ImageX0>=0)?0:(PIX-1);
            dcell=GetImageXYCoo(xindex+index*PIX,ImageX,ImageY,-1,false);
            bool xsign=(x1>=ImageX)?true:false;
            XY2[0]=ImageX+(xsign?1:-1)*dcell/2;
            XY2[1]=(XY2[0]*sin(phi)+CC)/cos(phi);
            XY2[2]=XY2[0]*cos(phi)+XY2[1]*sin(phi);
            XY2[3]=1;
            exist=2;
            break;
         }
         else if(x1*x2<0){
            dcell=GetImageXYCoo(0+index*PIX,ImageX,ImageY);
            XY1[0]=ImageX+(xsign0?1:-1)*dcell/2;
            XY1[1]=(XY1[0]*sin(phi)+CC)/cos(phi);
            XY1[2]=XY1[0]*cos(phi)+XY1[1]*sin(phi);
            XY1[3]=1;
            dcell=GetImageXYCoo(PIX-1+index*PIX,ImageX,ImageY);
            XY2[0]=ImageX+(xsign0?-1:1)*dcell/2;
            XY2[1]=(XY2[0]*sin(phi)+CC)/cos(phi);
            XY2[2]=XY2[0]*cos(phi)+XY2[1]*sin(phi);
            XY2[3]=1;
            if(XY1[2]>XY2[2]){
               double nbuff[4]={XY2[0],XY2[1],XY2[2],XY2[3]};
               for(int ii=0;ii<4;ii++) XY2[ii]=XY1[ii];
               for(int ii=0;ii<4;ii++) XY1[ii]=nbuff[ii];
            }
            exist=2;
            break;
         }
         index--;
      }
      if(exist<2) return false;
   }
   if(jdebug>0) printf("TelGeoFit::GetRange:output,CC=%lf phi=%lf XY1={%lf,%lf} XY2={%lf,%lf}\n",CC/PI*180,phi/PI*180,XY1[0],XY1[1],XY2[0],XY2[1]);
   if(XY1[2]>XY2[2]){
      double nbuff[4]={XY2[0],XY2[1],XY2[2],XY2[3]};
      for(int ii=0;ii<4;ii++) XY2[ii]=XY1[ii];
      for(int ii=0;ii<4;ii++) XY1[ii]=nbuff[ii];
   }
   return true;
}
double TelGeoFit::GetWidth(double npe[MAXPMT],double CC,double phi){
   double cont_sum=0;
   TH1F hbuff("width_name",";Short Axis;Npe [pe]",26*6,-13./180*PI,13./180*PI);
   for(int ii=0;ii<MAXPMT;ii++){
      int isipm0=ii;
      double content=npe[isipm0];
      if(content<=0) continue;
      cont_sum+=content;

      double ImageX=0,ImageY=0;
      double dcell=TelGeoFit::GetImageXYCoo(isipm0,ImageX,ImageY,-1,false);
      double ImageX2=ImageX*cos(phi)+ImageY*sin(phi); //long axis
      double ImageY2=ImageY*cos(phi)-ImageX*sin(phi)-CC; //short axis
      hbuff.Fill(ImageY2,content);
   }
   if(cont_sum<=0) return -1;

   double sigm=-1;
   int ibin0=hbuff.GetXaxis()->FindBin(0.);
   double integral=0;
   for(int ibin=0;ibin<(hbuff.GetNbinsX()/2)+2;ibin++){
      integral+=hbuff.GetBinContent(ibin0+ibin);
      if(ibin>0) integral+=hbuff.GetBinContent(ibin0-ibin);
      if(fabs(integral/cont_sum)>0.95){
         sigm=hbuff.GetXaxis()->GetBinCenter(ibin0+ibin)-hbuff.GetXaxis()->GetBinCenter(ibin0);
         break;
      }
   }

   return sigm;
}
bool TelGeoFit::GetLongRange(double CC,double phi,double width,double LRange[2]){
   if(width<=0) return false;
   double XY1[4],XY2[4];
   bool inside=TelGeoFit::GetRange(CC,phi,XY1,XY2);
   if(!inside) return false;
   for(int ii=0;ii<2;ii++){
      double xx=(ii==0)?XY1[0]:XY2[0];
      double yy=(ii==0)?XY1[1]:XY2[1];
      double ImageX,ImageY1,ImageY2;
      double width1=TelGeoFit::GetImageXYCoo(PIX/2+0*PIX,ImageX,ImageY1,-1,false);
      double width2=TelGeoFit::GetImageXYCoo(PIX/2+(PIX-1)*PIX,ImageX,ImageY2,-1,false);
      if(ImageY1>ImageY2){
         ImageY1+=width1/2;
         ImageY2-=width2/2;
      }
      else{
         ImageY1-=width1/2;
         ImageY2+=width2/2;
      }
      double margin=0.1/180*PI;
      bool yside=(fabs(yy-ImageY1)<margin||fabs(yy-ImageY2)<margin)&&(fabs(sin(phi))>0.5/16.);
      double deltaL=0;
      if(yside) deltaL=fabs(width/tan(phi));
      else deltaL=fabs(width*tan(phi));
      if(ii==0) LRange[ii]=XY1[2]+deltaL;
      else LRange[ii]=XY2[2]-deltaL;
   }
   return (LRange[1]>LRange[0]);
}
TH1F* TelGeoFit::GetScatterAngle(int ievt,bool CleanEdge,bool UseFit){
   if(ievt<0||ievt>=nevent) return 0;
   double CC,phi;
   double zenith=0,azimuth=0;
   int itel=telindex[ievt]-1;
   if(itel<0){
      cerr<<"TelGeoFit::GetScatterAngle: No Telescope Information Loaded."<<endl;
      return 0;
   }
   bool isfit=TelZen_fit[itel]>=0;
   zenith=isfit?TelZen_fit[itel]:TelZen_int[itel];
   azimuth=isfit?TelAzi_fit[itel]:TelAzi_int[itel];
   int irot=RotateDB::GetHead()->GetLi((rabbitTime[ievt]-(int)(rabbitTime[ievt]))/2.0e-8);
   if(irot<0){
      cerr<<"TelGeoFit::GetScatterAngle: Not a laser event."<<endl;
      return 0;
   }
   int iRot=RotateDB::rotindex[irot];
   double rotcoo[3];
   for(int ii=0;ii<3;ii++) rotcoo[ii]=isfit?RotPos_fit[iRot-1][ii]:RotPos_int[iRot-1][ii];
   double rotdir[2];
   double ele_in=PI/2-rotdir0[ievt][0]; //from zenith to elevation
   double azi_in=rotdir0[ievt][1];
   double rotpar[4]={isfit?RotEle_fit[iRot-1]:RotEle_int[iRot-1],isfit?RotAzi_fit[iRot-1]:RotAzi_int[iRot-1],isfit?RefEle_fit[iRot-1]:RefEle_int[iRot-1],isfit?RefAzi_fit[iRot-1]:RefAzi_int[iRot-1]};
   CalDir_out(rotpar[0],rotpar[1],rotpar[2],rotpar[3],ele_in,azi_in,rotdir[0],rotdir[1]);
   rotdir[0]=PI/2-rotdir[0]; //from elevation to zenith

   if(UseFit){
      CC=image_pars[ievt][0][0];
      phi=image_pars[ievt][1][0];
   }
   else{
      TelGeoFit::GetCCphi(zenith,azimuth,telpos[itel],rotcoo,rotdir,CC,phi);
   }

   double LRange[2]={-30,30};
   if(CleanEdge){
      double width=GetWidth(Npe_sipm[ievt],CC,phi);
      bool longrange=TelGeoFit::GetLongRange(CC,phi,width,LRange);
      if(!longrange) return 0;
      LRange[0]*=180/PI;
      LRange[1]*=180/PI;
   }

   static int ihist=0;
   TH1F* hh=new TH1F(Form("scatter_longcoo_%d",ihist++),Form(";cos(scatter angle);%s","Npe [pe]"),100,-1,1);

   for(int ii=0;ii<MAXPMT;ii++){
      int isipm=ii;
      double content=Npe_sipm[ievt][ii];
      double econtent=10.;
      if(content<=0) continue;

      double ImageX,ImageY;
      double width=TelGeoFit::GetImageXYCoo(isipm,ImageX,ImageY,-1,false);
      double xrange[2]={ImageX-width/2,ImageX+width/2};
      double yrange[2]={ImageY-width/2,ImageY+width/2};

      int ndiv=5;
      for(int ibin1=0;ibin1<ndiv;ibin1++){
         double xx=xrange[0]+(xrange[1]-xrange[0])/ndiv*(ibin1+0.5);
         for(int ibin2=0;ibin2<ndiv;ibin2++){
            double yy=yrange[0]+(yrange[1]-yrange[0])/ndiv*(ibin2+0.5);
            double longaxis=xx*cos(phi)+yy*sin(phi); //long axis
            double ycoo=yy*cos(phi)-xx*sin(phi)-CC; //short axis
            if(CleanEdge) {if(longaxis<=LRange[0]||longaxis>=LRange[1]) continue;}
            double scatter_angle=TelGeoFit::GetScatAngleFromLongcoo(zenith,azimuth,telpos[itel],rotcoo,rotdir,longaxis);
            int ibin=hh->GetXaxis()->FindBin(cos(scatter_angle));
            hh->SetBinContent(ibin,hh->GetBinContent(ibin)+content/(ndiv*ndiv));
            hh->SetBinError(ibin,sqrt(pow(hh->GetBinError(ibin),2)+pow(econtent,2)/(ndiv*ndiv)));
         }
      }
   }
   if(hh->Integral()<=0){
      delete hh;
      hh=0;
   }
   return hh;
}

bool TelGeoFit::GetImageCoo(double zenith,double azimuth,double dir_in[3],double &xx,double &yy,bool IsLocal){
   double norm=sqrt(dir_in[0]*dir_in[0]+dir_in[1]*dir_in[1]+dir_in[2]*dir_in[2]);
   double theta=IsLocal?(PI/2):(PI/2-zenith);
   double phi=IsLocal?0:azimuth;
   double zz0=(dir_in[0]*cos(phi)+dir_in[1]*sin(phi))/norm*cos(theta)+dir_in[2]/norm*sin(theta);
   if(zz0>=0) return false;
   double xx0=(dir_in[0]*cos(phi)+dir_in[1]*sin(phi))/norm*sin(theta)-dir_in[2]/norm*cos(theta);
   double yy0=(-dir_in[0]*sin(phi)+dir_in[1]*cos(phi))/norm;
   xx=-(yy0/(-zz0));
   yy=(xx0/(-zz0));

   //double indir[2];
   //indir[0]=asin(-dir_in[2]/norm);
   //indir[1]=acos(-dir_in[0]/sqrt(dir_in[0]*dir_in[0]+dir_in[1]*dir_in[1]));
   //if(dir_in[1]>0) indir[1]=2*PI-indir[1];
   //double outdir[2];
   //slaDtp2s(xx/180*PI,yy/180*PI,azimuth,PI/2-zenith,&(outdir[1]),&(outdir[0]));
   //printf("Indir={%.2lf,%.2lf} ImageXY={%.2lf,%.2lf} Outdir={%.2lf,%.2lf}\n",indir[0]/PI*180,indir[1]/PI*180,xx,yy,outdir[0]/PI*180,outdir[1]/PI*180);
   return true;
}
void TelGeoFit::GetOutDir(double zenith,double azimuth,double imagexy[2],double dir_out[3],bool IsLocal){
   double theta=IsLocal?(PI/2):(PI/2-zenith);
   double phi=IsLocal?0:azimuth;
   double xx0=imagexy[1];
   double yy0=-imagexy[0];
   dir_out[0]=(-xx0*sin(theta)+cos(theta))*cos(phi)+yy0*sin(phi);
   dir_out[1]=(-xx0*sin(theta)+cos(theta))*sin(phi)-yy0*cos(phi);
   dir_out[2]=xx0*cos(theta)+sin(theta);
   double norm=sqrt(dir_out[0]*dir_out[0]+dir_out[1]*dir_out[1]+dir_out[2]*dir_out[2]);
   for(int ii=0;ii<3;ii++) dir_out[ii]/=norm;
}
void TelGeoFit::Getnz(double nz,double &nz1,double &nz2){
   nz1=nz/sqrt(1+nz);
   nz2=1./sqrt(1+nz);
   if(isinf(nz)){
      if(signbit(nz)) nz1=-1;
      else nz1=1;
      nz2=0;
   }
}
double TelGeoFit::GetImageCoo(double CC,double phi,double longcoo,int icoo){
   if(icoo==0) return (-CC*sin(phi)+longcoo*cos(phi));
   else return (CC*cos(phi)+longcoo*sin(phi));
}
void TelGeoFit::GetOutDir(double zenith,double azimuth,double CC,double phi,double longcoo,double dir_out[3]){
   double imagexy[2]={GetImageCoo(CC,phi,longcoo,0),GetImageCoo(CC,phi,longcoo,1)};
   GetOutDir(zenith,azimuth,imagexy,dir_out);
}

bool TelGeoFit::Getnz(double telcoo[3],double incoo[3],double indir[2],double &planephi,double &nz,double telzcoo[3]){ //indir[0] is the zenith angle
   double intheta=indir[0];
   double inphi=indir[1];
   if((incoo[2]!=telcoo[2])&&cos(intheta)!=0){
      double dz=(telcoo[2]-incoo[2])/cos(intheta);
      telzcoo[2]=telcoo[2];
      telzcoo[0]=incoo[0]+(sin(intheta)*cos(inphi))*dz;
      telzcoo[1]=incoo[1]+(sin(intheta)*sin(inphi))*dz;
   }
   else{
      for(int ii=0;ii<3;ii++) telzcoo[ii]=incoo[ii];
   }
   double xdir[3]={telzcoo[0]-telcoo[0],telzcoo[1]-telcoo[1],telzcoo[2]-telcoo[2]};
   double norm=sqrt(xdir[0]*xdir[0]+xdir[1]*xdir[1]+xdir[2]*xdir[2]);
   if(norm<=0) return false;
   double dist=sqrt(xdir[0]*xdir[0]+xdir[1]*xdir[1]);
   if(dist<=0) planephi=0;
   else{
      planephi=acos(xdir[0]/dist);
      if(xdir[1]<0) planephi=2*PI-planephi;
   }
   nz=sin(indir[0])*sin(indir[1]-planephi)/cos(indir[0]);
   return true;
}
double TelGeoFit::GetApar(double &Anz,double &phi,double zenith,double azimuth,double planephi,double nz){ //A or A*nz is a very important parameter
   double phi0=azimuth;
   double theta=PI/2-zenith;

   double AA=0;
   bool finitenz=isfinite(nz);

   double p1,p2,p3;
   if(finitenz){ //CC=AA*p1,cos(phi)=AA*p2,sin(phi)=AA*p3
      p1=sin(theta)*nz-cos(theta)*sin(phi0-planephi);
      p2=-cos(theta)*nz-sin(theta)*sin(phi0-planephi);
      p3=-cos(phi0-planephi);
   }
   else{ //CC=Anz*p1,cos(phi)=Anz*p2,sin(phi)=Anz*p3;
      p1=sin(theta);
      p2=-cos(theta);
      p3=0;
   }

   double norm=sqrt(p2*p2+p3*p3);
   if(norm==0){ //CC is infinite
      Anz=0;
      phi=0;
      return AA;
   }
   else{
      double phi_ref=acos(p3/norm);
      if(p2<0) phi_ref=2*PI-phi_ref;
      double phi1=CommonTools::ProcessAngle(-phi_ref+PI/2);
      double phi2=CommonTools::ProcessAngle(-phi_ref-PI/2);
      if(fabs(phi1-PI/2)==fabs(phi2-PI/2)) phi=(phi1>=0&&phi1<PI)?phi1:phi2;
      else if(fabs(phi1-PI/2)<fabs(phi2-PI/2)) phi=phi1;
      else phi=phi2;
      AA=fabs(p2)<fabs(p3)?(sin(phi)/p3):(cos(phi)/p2);
   }

   if(finitenz){
      Anz=AA*nz;
      return AA;
   }
   else{
      Anz=AA;
      return 0;
   }
}
void TelGeoFit::GetCCphi(double zenith,double azimuth,double planephi,double nz,double &CC,double &phi){
   double theta=PI/2-zenith;
   double phi0=azimuth;
   double Anz;
   double AA=GetApar(Anz,phi,zenith,azimuth,planephi,nz);
   CC=Anz*sin(theta)-AA*cos(theta)*sin(phi0-planephi);
}
void TelGeoFit::GetCCphi(double zenith,double azimuth,double telcoo[3],double incoo[3],double indir[2],double &CC,double &phi){ //indir[0] is zenith angle
   double planephi,nz;
   double telzcoo[3];
   Getnz(telcoo,incoo,indir,planephi,nz,telzcoo);
   return GetCCphi(zenith,azimuth,planephi,nz,CC,phi);
}

void TelGeoFit::CalPlane(double zenith,double azimuth,double CC,double phi,double &planephi,double &nz,int &signnz,bool sign){
   double phi0=azimuth;
   double theta=PI/2-zenith;

   double outdir[2][3];
   GetOutDir(zenith,azimuth,CC,phi,0./180*PI,outdir[0]);
   GetOutDir(zenith,azimuth,CC,phi,1./180*PI,outdir[1]);

   ///three equations/////
   //CC*cos(theta)+cos(phi)*sin(theta)=p1=-AA*sin(phi0-planephi)
   //sin(phi)=p2=-AA*cos(phi0-planephi)
   //CC*sin(theta)-cos(phi)*cos(theta)=p3=AA*nz
   //////////////////////
   double p1=CC*cos(theta)+cos(phi)*sin(theta);
   double p2=sin(phi);
   double p3=CC*sin(theta)-cos(phi)*cos(theta);
   double absAA=sqrt(p1*p1+p2*p2);

   //test the sign of AA
   bool signAA=true;
   double AA=signAA?absAA:(-absAA);

   double Acos_planephi=-p1*sin(phi0)-p2*cos(phi0);
   bool sign_cos_planephi=(!signbit(Acos_planephi))==signAA;
   double Asin_planephi=p1*cos(phi0)-p2*sin(phi0);
   bool sign_sin_planephi=(!signbit(Asin_planephi))==signAA;
   planephi=acos(absAA==0?(sign_cos_planephi?1:-1):(Acos_planephi/AA));
   if(!sign_sin_planephi) planephi=2*PI-planephi;
   //define the sign of AA and planephi
   double angle0=acos(cos(planephi)*outdir[0][0]+sin(planephi)*outdir[0][1]);
   double angle1=acos(cos(planephi)*outdir[1][0]+sin(planephi)*outdir[1][1]);
   bool time_increase=(GoDown==(angle1<angle0));
   if(time_increase!=sign){
      signAA=(!signAA);
      AA*=-1;
      planephi=CommonTools::ProcessAngle(planephi+PI,false);
   }

   nz=signAA?(p3/absAA):((-p3)/absAA);
   signnz=signbit(nz)?-1:1;
}
void TelGeoFit::CalPlane(double zenith,double azimuth,double CC,double phi,double &planephi,double &nz,int &signnz,double refplanephi){
   CalPlane(zenith,azimuth,CC,phi,planephi,nz,signnz,true);
   double outdir[3];
   GetOutDir(zenith,azimuth,CC,phi,0./180*PI,outdir);
   double spec_angle=acos(cos(planephi)*outdir[0]+sin(planephi)*outdir[1]);
   if(fabs(spec_angle-refplanephi)>PI/2){
      planephi=2*PI-planephi;
      signnz*=-1;
      nz*=-1;
   }
}
void TelGeoFit::CalPlane(double planephi,double nz,double xyzdir[3][3]){
   xyzdir[0][0]=cos(planephi);
   xyzdir[0][1]=sin(planephi);
   xyzdir[0][2]=0;
   xyzdir[2][0]=sin(planephi);
   xyzdir[2][1]=-cos(planephi);
   xyzdir[2][2]=nz;
   if(isinf(nz)){
      xyzdir[2][0]=0;
      xyzdir[2][1]=0;
      xyzdir[2][2]=(signbit(nz))?-1:1;
   }
   double norm=sqrt(pow(xyzdir[2][0],2)+pow(xyzdir[2][1],2)+pow(xyzdir[2][2],2));
   for(int icoo=0;icoo<3;icoo++) xyzdir[2][icoo]/=norm;
   Laser::cross(xyzdir[2],xyzdir[0],xyzdir[1]);
}

/*int TelGeoFit::CalTelDir(double CC,double phi,double planephi,double nz,int &nsol,double ele_out[4],double azi_out[4],int signA[4],double ele_in,double azi_in){
   nsol=0;
   bool finitenz=isfinite(nz);

   double margin=1.0e-6;
   //calculate elevation
   //when phi=PI/2,CC=0(or theta=0,planephi=azimuth), then we can't measure the elevation of telescope in this situation, but we can measure the planephi quite precisely
   if(!finitenz){ //nz is infinite
      if(fabs(sin(phi))>margin){
         printf("WFCTAEvent::CalTelDir: error of the input CC=%lf and phi=%lf when nz=%lf\n",CC/PI*180,phi/PI*180,nz);
         return -2;
      }
      else{
         printf("WFCTAEvent::CalTelDir: any azimuth is correct. Set azimuth to 8*PI degree\n");
         double p1=CC;
         double p2=cos(phi);
         double norm=sqrt(p1*p1+p2*p2);
         double angle_ref=acos(p1/norm);
         if(p2<0) angle_ref=2*PI-angle_ref;
         int result=-1;
         nsol=2;
         azi_out[0]=8*PI;
         ele_out[0]=CommonTools::ProcessAngle(angle_ref+PI/2);
         azi_out[1]=8*PI;
         ele_out[1]=CommonTools::ProcessAngle(angle_ref-PI/2);
         double maxcos=-1;
         for(int ii=0;ii<nsol;ii++){
            if(cos(ele_out[ii]-ele_in)>maxcos){
               result=ii;
               maxcos=cos(ele_out[ii]-ele_in);
            }
            double var=CC*sin(ele_out[ii])-cos(ele_out[ii])*cos(phi);
            int signnz=(signbit(nz))?-1:1;
            signA[ii]=(var/signnz>=0)?1:-1;
         }
         return result;
      }
   }
   else{
      double p2=CC;
      double p1=-cos(phi);
      double norm=sqrt(p1*p1+p2*p2);
      double absAnz=sqrt((1+CC*CC)/(1+nz*nz))*fabs(nz);
      if(norm<absAnz){
         printf("WFCTAEvent::CalTelDir: no theta is correct (abs(Anz)=%lf norm=%lf). Exiting...\n",absAnz,norm);
         return -3;
      }
      else if(norm<margin){
         printf("WFCTAEvent::CalTelDir: any theta is correct. Set theta to PI/2+8*PI degree\n");
         int result=-1;
         nsol=2;
         ele_out[0]=PI/2+8*PI;
         azi_out[0]=CommonTools::ProcessAngle(planephi+0);
         signA[0]=1;
         ele_out[1]=PI/2+8*PI;
         azi_out[1]=CommonTools::ProcessAngle(planephi+PI);
         signA[1]=-1;
         double maxcos=-1;
         for(int ii=0;ii<nsol;ii++){
            if(cos(azi_out[ii]-azi_in)>maxcos){
               result=ii;
               maxcos=cos(azi_out[ii]-azi_in);
            }
         }
         return result;
      }
      else{
         double angle_ref=acos(p1/norm);
         if(p2<0) angle_ref=2*PI-angle_ref;
         int signnz=(nz>=0)?1:-1;
         int result=-1;
         nsol=0;
         double maxcos=-1;
         for(int ii=0;ii<4;ii++){
            signA[nsol]=(ii/2<1)?1:-1;
            ele_out[nsol]=CommonTools::ProcessAngle(angle_ref+(((ii%2)==0)?1:-1)*acos(signA[nsol]*signnz*absAnz/norm));
            //calculate azimuth
            //when m=0 and n=0, then we can't measure the azimuth angle of telescope
            double p11=-(CC*cos(ele_out[nsol])+sin(ele_out[nsol])*cos(phi));
            double p22=sin(phi);
            double AA=signA[nsol]*sqrt(p11*p11+p22*p22);
            if(fabs(AA)<margin){
               printf("WFCTAEvent::CalTelDir: any azimuth is correct. Set azimuth to 0 degree\n");
               azi_out[nsol]=0;
            }
            else{
               double azi_ref=acos(p22/AA);
               if(p11/AA<0) azi_ref=2*PI-azi_ref;
               azi_out[ii]=CommonTools::ProcessAngle(planephi+azi_ref);
            }
            if(cos(ele_out[nsol]-ele_in)>maxcos){
               result=nsol;
               maxcos=cos(ele_out[nsol]-ele_in);
            }
            //printf("test: nsol=%d iele=%lf iazi=%lf isignA=%d maxcos=%lf result=%d\n",nsol,ele_out[ii]/PI*180,azi_out[ii]/PI*180,signA[ii],maxcos,result);
            nsol++;
         }
         return result;
      }
   }
}
int WFCTAEvent::CalTelDir(double CC,double phi,double* lasercoo,double* laserdir,int &nsol,double ele_out[4],double azi_out[4],int signA[4],double ele_in,double azi_in){
   double incoo[3];
   CooCorr(lasercoo,laserdir,incoo);
   double dist=sqrt(pow(incoo[0],2)+pow(incoo[1],2));
   if(dist<=0) return -1;
   double planephi=acos(incoo[0]/dist);
   if(incoo[1]<0) planephi=2*PI-planephi;
   double nz=(sin(laserdir[0])*sin(laserdir[1]-planephi))/cos(laserdir[0]);
   return CalTelDir(CC,phi,planephi,nz,nsol,ele_out,azi_out,signA,ele_in,azi_in);
}
int WFCTAEvent::CalTelDir(double CC,double phi,double planephi,double nz,double &elevation,double &azimuth,double ele_in,double azi_in){
   int nsol=0;
   double ele_out[4],azi_out[4];
   int signA[4];
   int result=CalTelDir(CC,phi,planephi,nz,nsol,ele_out,azi_out,signA,ele_in,azi_in);
   if(result<0) return result;
   int signA0=GetSign(CheckLaser(),CheckMC())>=0?1:-1;
   int signAnz=(signA0*nz>=0)?1:-1;
   int result2=-10;
   double maxcos=-1;
   double maxchi=-1;
   bool Isele=true;
   for(int ii=0;ii<nsol;ii++){
      if(fabs(azi_out[ii])>7*PI) {Isele=false; break;}
   }
   for(int ii=0;ii<nsol;ii++){
      int signA1=(signA[ii]>=0)?1:-1;
      if(signA0!=signA1) continue;
      //printf("Sign: ii=%d signA=%d(%d) signA0=%d result={%d,%d}\n",ii,signA[ii],signA1,signA0,result,result2);
      if(result==ii) {result2=ii; break;}
      else{
         //double icos=(Isele)?cos(ele_out[ii]-ele_in):cos(azi_out[ii]-azi_in);
         //if(icos>maxcos){
         //   result2=ii;
         //   maxcos=icos;
         //}
         //
      }
   }
   //printf("result2=%d nsol=%d result=%d\n",result2,nsol,result);
   if(result2<0) return result2;
   else{
      elevation=ele_out[result2]/PI*180;
      azimuth=azi_out[result2]/PI*180;
      return 1;
   }
}
int WFCTAEvent::CalTelDir(double CC,double phi,double* lasercoo,double* laserdir,double &elevation,double &azimuth,double ele_in,double azi_in){
   double incoo[3];
   CooCorr(lasercoo,laserdir,incoo);
   double dist=sqrt(pow(incoo[0],2)+pow(incoo[1],2));
   if(dist<=0) return -1;
   double planephi=acos(incoo[0]/dist);
   if(incoo[1]<0) planephi=2*PI-planephi;
   double nz=(sin(laserdir[0])*sin(laserdir[1]-planephi))/cos(laserdir[0]);
   return CalTelDir(CC,phi,planephi,nz,elevation,azimuth,ele_in,azi_in);
}
int WFCTAEvent::GetTelDir(double &elevation,double &errel,double &azimuth,double &erraz,double ele_in,double azi_in){
   if(!minimizer) return -1;
   double lasercoo[3]={laserevent.LaserCoo[0],laserevent.LaserCoo[1],laserevent.LaserCoo[2]};
   double lasertheta=laserevent.LaserDir[0]/180.*PI;
   double laserphi=laserevent.LaserDir[1]/180.*PI;
   double laserdir[2]={lasertheta,laserphi};
   double xdir[3];
   CooCorr(lasercoo,laserdir,xdir);

   double phi=minimizer->X()[3];
   double CC=minimizer->X()[2];
   int res=CalTelDir(CC,phi,xdir,laserdir,elevation,azimuth,ele_in,azi_in);
   if(res<=0) return res;

   double Cov[3]={minimizer->CovMatrix(2,2),minimizer->CovMatrix(2,3),minimizer->CovMatrix(3,3)};
   errel=0; erraz=0;
   double coeff[2][5];
   for(int iv=0;iv<5;iv++){
      coeff[0][iv]=0;
      coeff[1][iv]=0;
      double CC1=CC;
      double phi1=phi;
      double lasercoo1[3]={xdir[0],xdir[1],0};
      double laserdir1[2]={laserdir[0],laserdir[1]};
      double delta=0;
      if(iv==0){
         CC1=CC+sqrt(Cov[0]);
         delta=CC1-CC;
      }
      else if(iv==1){
         phi1=phi+sqrt(Cov[2]);
         if(phi1>0.999*PI) phi1=phi-sqrt(Cov[2]);
         delta=phi1-phi;
      }
      else if(iv==2){
         lasercoo1[0]=xdir[0]+Laser::LaserCooErr;
         delta=lasercoo1[0]-lasercoo[0];
      }
      else if(iv==3){
         lasercoo1[1]=xdir[1]+Laser::LaserCooErr;
         delta=lasercoo1[1]-lasercoo[1];
      }
      else if(iv==4){
         laserdir1[0]=laserdir[0]+Laser::LaserZenErr/180.*PI;
         delta=laserdir1[0]-laserdir[0];
      }
      else if(iv==5){
         laserdir1[1]=laserdir[1]+Laser::LaserAziErr/180.*PI;
         delta=laserdir1[1]-laserdir[1];
      }
      double elevation1,azimuth1;
      int res1=CalTelDir(CC1,phi1,lasercoo1,laserdir1,elevation1,azimuth1,ele_in,azi_in);
      if(fabs(azimuth1-azimuth)>300){
         azimuth1+=(azimuth1>azimuth)?(-360):(360);
      }
      if(res1<=0) continue;
      errel+=pow(elevation1-elevation,2);
      erraz+=pow(azimuth1-azimuth,2);
      if(delta!=0){
         coeff[0][iv]=(elevation1-elevation)/delta;
         coeff[1][iv]=(azimuth1-azimuth)/delta;
      }
      //printf("WFCTAEvent::GetTelDir: iv=%d inpar=%lf inpar2=%lf ele={%lf,%lf} azi={%lf,%lf}\n",iv,CC,CC1,elevation,elevation1,azimuth,azimuth1);
   }
   errel+=2*coeff[0][0]*coeff[0][1]*Cov[1];
   erraz+=2*coeff[1][0]*coeff[1][1]*Cov[1];
   errel=sqrt(errel);
   erraz=sqrt(erraz);
   //printf("WFCTAEvent::GetTelDir: CC={%lf,%lf} phi={%lf,%lf}\n",CC,sqrt(Cov[0]),phi,sqrt(Cov[2]));
   return 1;
}
bool WFCTAEvent::CalPHIRange0(double zenith,double azimuth,double planephi,double nz,double PHI_in,double* PHIRange){
   //the PHI range to have image ((x*cos(phi0)+y*sin(phi0))*cos(theta)+z*sin(theta)>=0)
   double theta=PI/2-zenith;
   double phi0=azimuth;
   double nz1,nz2;
   Getnz(nz,nz1,nz2);
   double p21=cos(theta)*cos(phi0-planephi);
   double p22=nz1*sin(phi0-planephi)*cos(theta)+nz2*sin(theta);
   double norm2=sqrt(p21*p21+p22*p22);
   double margin=1.0e-5;
   if(norm2<margin){ //any PHI is correct
      PHIRange[0]=0;
      PHIRange[1]=2*PI;
      return true;
   }
   else{
      double PHI_ref=acos(p21/norm2);
      if(p22<0) PHI_ref=2*PI-PHI_ref;
      PHIRange[0]=PHI_ref-PI/2;
      PHIRange[1]=PHI_ref+PI/2;
      return (cos(PHI_in-PHI_ref)>=margin);
   }
}
bool WFCTAEvent::GetImageXYCoo(double zenith,double azimuth,double planephi,double nz,double PHI_in,double &xx,double &yy){
   double PHIRange[2];
   if(!CalPHIRange0(zenith,azimuth,planephi,nz,PHI_in,PHIRange)) return false;
   double theta=PI/2-zenith;
   double phi0=azimuth;
   double nz1,nz2;
   Getnz(nz,nz1,nz2);
   double denum=(cos(PHI_in)*cos(phi0-planephi)+nz1*sin(PHI_in)*sin(phi0-planephi))*cos(theta)+nz2*sin(PHI_in)*sin(theta);
   double x2=-( (cos(PHI_in)*cos(phi0-planephi)+nz1*sin(PHI_in)*sin(phi0-planephi))*sin(theta)-nz2*sin(PHI_in)*cos(theta) )/denum;
   double y2=-( -cos(PHI_in)*sin(phi0-planephi)+nz1*sin(PHI_in)*cos(phi0-planephi) )/denum;
   xx=y2;
   yy=x2;
   return denum>=0;
}
bool WFCTAEvent::GetImageXYCoo(double zenith,double azimuth,double* incoo,double* indir,double PHI_in,double &xx,double &yy){
   double theta=PI/2-zenith;
   double phi0=azimuth;
   double planephi,nz;
   Getnz(incoo,indir,planephi,nz);                                                                                          {
   bool res=GetImageXYCoo(zenith,azimuth,planephi,nz,PHI_in,xx,yy);
   //printf("zenith=%lf azi=%lf PHI=%lf indir={%lf,%lf} planephi=%lf nz=%lf xy={%lf,%lf}\n",zenith/PI*180,azimuth/PI*180,PHI_in/PI*180,indir[0]/PI*180,indir[1]/PI*180,planephi/PI*180,nz,xx/PI*180,yy/PI*180);
   return res;

   //double theta=PI/2-zenith;
   //double phi0=azimuth;
   //double dist=sqrt(pow(incoo[0],2)+pow(incoo[1],2));
   //double zdir[3];
   //zdir[0]=incoo[1]/dist*cos(indir[0]);
   //zdir[1]=-incoo[0]/dist*cos(indir[0]);
   //zdir[2]=sin(indir[0])*(incoo[0]/dist*sin(indir[1])-incoo[1]/dist*cos(indir[1]));
   //double planephi=acos(incoo[0]/dist);
   //if(incoo[1]<0) planephi=2*PI-planephi;
   //double m2n2=sqrt(zdir[0]*zdir[0]+zdir[1]*zdir[1]);




   //double x1=cos(phi0-planephi)*cos(PHI_in)*m2n2+sin(phi0-planephi)*zdir[2]*sin(PHI_in);
   //double y1=-sin(phi0-planephi)*cos(PHI_in)*m2n2+cos(phi0-planephi)*zdir[2]*sin(PHI_in);

   //double x2=-(x1*sin(theta)-cos(theta)*m2n2*sin(PHI_in))/(x1*cos(theta)+sin(theta)*m2n2*sin(PHI_in));
   //double y2=-y1/(x1*cos(theta)+sin(theta)*m2n2*sin(PHI_in));
   //xx=y2;
   //yy=x2;
}
void WFCTAEvent::GetPHI(double zenith,double azimuth,double planephi,double nz,double* ImageCoo,double &PHI_in){
   double theta=PI/2-zenith;
   double phi0=azimuth;
   double nz1,nz2;
   Getnz(nz,nz1,nz2);
   double CC,phi;
   GetCCphi(zenith,azimuth,planephi,nz,CC,phi);

   //new method
   double x2=ImageCoo[1];
   double y2=ImageCoo[0];
   double dirx[3]={sin(theta)*cos(planephi-phi0),sin(planephi-phi0),cos(theta)*cos(planephi-phi0)};
   double diry[3]={-nz1*sin(theta)*sin(planephi-phi0)-nz2*cos(theta),nz1*cos(planephi-phi0),-nz1*cos(theta)*sin(planephi-phi0)+nz2*sin(theta)};
   double cosPHI=(-dirx[0]*x2-dirx[1]*y2+dirx[2])/sqrt(1+x2*x2+y2*y2);
   PHI_in=acos(cosPHI);
   if(-diry[0]*x2-diry[1]*y2+diry[2]<0) PHI_in=2*PI-PHI_in;
   //double x2=ImageCoo[1];
   //double y2=ImageCoo[0];
   //double px1=(x2*cos(theta)+sin(theta))*cos(phi0-planephi);
   //double px2=((x2*cos(theta)+sin(theta))*sin(phi0-planephi)*nz1+nz2*(x2*sin(theta)-cos(theta)));
   //double py1=y2*cos(theta)*cos(phi0-planephi)-sin(phi0-planephi);
   //double py2=(y2*cos(theta)*sin(phi0-planephi)+cos(phi0-planephi))*nz1+nz2*y2*sin(theta);

   //double norm=0;
   //      PHI_in=0;              {
   //      return;
   //   }          {
   //   else{
   //      double PhI_ref=acos(px1/norm);
   //      if(px2<0) PhI_ref=2*PI-PhI_ref;
   //   }  double PHI1=CommonTools::ProcessAngle(PhI_ref+PI/2);
   //      d{uble PHI2=CommonTools::ProcessAngle(PhI_ref-PI/2);
   //      if(fabs(PHI1-PI/2)<=fabs(PHI2-PI/2)) PHI_in=PHI1;
   //      else PHI_in=PHI2;
   //      return;




   //   }
   //}
   //else{ //use y to eval PHI
   //   norm=sqrt(py1*py1+py2*py2);
   //   if(norm<=0){
   //      printf("WFCTAEvent::GetPHI: no solution with y=%lf\n",y2);
   //      PHI_in=0;
   //      return;
   //   }
   //   else{
   //      double PhI_ref=acos(py1/norm);
   //      if(py2<0) PhI_ref=2*PI-PhI_ref;
   //      double PHI1=CommonTools::ProcessAngle(PhI_ref+PI/2);
   //      double PHI2=CommonTools::ProcessAngle(PhI_ref-PI/2);
   //      if(fabs(PHI1-PI/2)<=fabs(PHI2-PI/2)) PHI_in=PHI1;
   //      else PHI_in=PHI2;
   //      return;
   //   }
   //}
}
void WFCTAEvent::GetPHI(double zenith,double azimuth,double incoo[3],double indir[2],double* ImageCoo,double &PHI_in){
   double planephi,nz;
   Getnz(incoo,indir,planephi,nz);
   return GetPHI(zenith,azimuth,planephi,nz,ImageCoo,PHI_in);
}
void WFCTAEvent::GetPHI2(double zenith,double azimuth,double CC,double phi,double* ImageCoo,double &PHI_in){
   double planephi,nz;
   int signnz;
   CalPlane(CC,phi,zenith,azimuth,planephi,nz,signnz);
   return  GetPHI(zenith,azimuth,planephi,nz,ImageCoo,PHI_in);
   //double theta=PI/2-zenith;
   //double phi0=azimuth;
   //double AA=GetApar(CC,phi,zenith,azimuth);
   //double x2=ImageCoo[1];
   //double y2=ImageCoo[0];
   //double margin=1.0e-5;
   //if(fabs(x2*cos(theta)+sin(theta))<margin){
   //   double ff=(AA*sin(phi)*cos(theta));
   //   double x1=-(x2*sin(theta)-cos(theta));
   //   double y1=-y2;
   //   double f2=(x2*cos(theta)+sin(theta));
   //   double cosPHI1=(x1+AA*AA*(CC*cos(theta)+sin(theta)*cos(phi))*(CC*sin(theta)-cos(theta)*cos(phi))*f2)/(AA*sin(phi))*ff;
   //   double cosPHI2=(y1-AA*AA*sin(phi)*(CC*sin(theta)-cos(theta)*cos(phi))*f2)/(AA*(CC*cos(theta)+sin(theta)*cos(phi)))*ff;
   //   double PHI1=acos(cosPHI1);
   //   double PHI2=acos(cosPHI2);
   //   PHI_in=(PHI1+PHI2)/2./PI*180.;
   //}
   //else{
   //   double x1=-(x2*sin(theta)-cos(theta))/(x2*cos(theta)+sin(theta));
   //   double y1=-y2/(x2*cos(theta)+sin(theta));

   //   double ctanPHI1=(x1+AA*AA*(CC*cos(theta)+sin(theta)*cos(phi))*(CC*sin(theta)-cos(theta)*cos(phi)))/(AA*sin(phi));
   //   double ctanPHI2=(y1-AA*AA*sin(phi)*(CC*sin(theta)-cos(theta)*cos(phi)))/(AA*(CC*cos(theta)+sin(theta)*cos(phi)));
   //   double PHI1=acos(ctanPHI1/sqrt(1+pow(ctanPHI1,2)));
   //   double PHI2=acos(ctanPHI2/sqrt(1+pow(ctanPHI2,2)));
   //   PHI_in=(PHI1+PHI2)/2./PI*180.;
   //}
}
int WFCTAEvent::GetRange(double zenith,double azimuth,double planephi,double dirin[3],double phiin,double CCin,double* PHI,double* XX,double* YY){
   double theta=PI/2-zenith;
   double phi0=azimuth;
   bool IsLaser;
   double margin=1.0e-5;
   double nz,nz1,nz2;
   double CC,phi;
   double normdir=sqrt(dirin[0]*dirin[0]+dirin[1]*dirin[1]+dirin[2]*dirin[2]);
   double cosPHIL;
   bool UsePhi=(normdir<=0)&&(planephi==0);
   //printf("UsePhi=%d\n",UsePhi);
   if(UsePhi){
      bool useclosephi=true; //to determine the planephi, otherwise it need to be determined by the time information
      //use phiin and CCin as input parameter
      if(phiin<0){
         IsLaser=false;
         phi=-phiin;
         CC=CCin;
      }
      else{
         IsLaser=true;
         phi=phiin;
         CC=CCin;
      }
      double p1=CC*cos(theta)+sin(theta)*cos(phi);
      double p2=sin(phi);
      double p3=CC*sin(theta)-cos(theta)*cos(phi);
      double norm1=sqrt(p1*p1+p2*p2);
      if(norm1==0){ //any planephi is fine, set it to phi0
         planephi=phi0;
         nz2=0;
         nz1=p3>=0?1:-1;
         nz=nz1/nz2;
      }
      else{
         double phi_ref=acos(p2/norm1);
         if(p1<0) phi_ref=2*PI-phi_ref;
         double planephi1=CommonTools::ProcessAngle(phi0+phi_ref);
         phi_ref=acos(-p2/norm1);
         if(-p1<0) phi_ref=2*PI-phi_ref;
         double planephi2=CommonTools::ProcessAngle(phi0+phi_ref);
         double phi00=CommonTools::ProcessAngle(phi0);
         bool usephi1=useclosephi?(fabs(planephi1-phi00)<fabs(planephi2-phi00)):(fabs(planephi1-phi00)>fabs(planephi2-phi00));
         planephi=usephi1?planephi1:planephi2;
         double Apre=fabs(p1)<fabs(p2)?(p2/cos(planephi-phi0)):(p1/sin(planephi-phi0));
         nz=p3/Apre;
         Getnz(nz,nz1,nz2);
      }
   }
   else{
      IsLaser=dirin[2]>=0;
      if(!IsLaser) for(int icoo=0;icoo<3;icoo++) dirin[icoo]*=(-1);
      double dir[2];
      dir[0]=acos(dirin[2]/normdir);
      double normdir1=sqrt(dirin[0]*dirin[0]+dirin[1]*dirin[1]);
      if(normdir1<=0) dir[1]=0.;
      else{
         dir[1]=acos(dirin[0]/normdir1);
         if(dirin[1]<0) dir[1]=2*PI-dir[1];
      }

      cosPHIL=(sin(dir[0])*cos(dir[1]-planephi));
      if((1-fabs(cosPHIL))<margin){ //no image will be detected
         if(jdebug>0) printf("WFCTAEvent::GetXYRange: return=-2, error between the dirin and planephi, because angle=%lf\n",acos(cosPHIL)/PI*180);
         return -2;
      }

      double incoo[3]={cos(planephi)*1000,sin(planephi)*1000,0};
      GetCCphi(zenith,azimuth,incoo,dir,CC,phi);
      double planephi0;
      Getnz(incoo,dir,planephi0,nz);
      Getnz(nz,nz1,nz2);
      //printf("CC=%lf phi=%lf incoo={%lf,%lf,%lf} indir={%lf,%lf} zenith=%lf azimuth=%lf\n",CC/PI*180,phi/PI*180,incoo[0],incoo[1],incoo[2],dir[0]/PI*180,dir[1]/PI*180,zenith/PI*180,azimuth/PI*180);
   }

   double PHI1,PHI2;
   double xyboundary[2][4];
   bool inside=GetRange(CC,phi,xyboundary[0],xyboundary[1]);
   if(!inside){
      if(jdebug>0) printf("WFCTAEvent::GetXYRange: return=-4, Image is not in the field of view of telescope CC=%lf phi=%lf\n",CC/PI*180,phi/PI*180);
      return -4;
   }
   double xyboun[2][2]={{xyboundary[0][1]/180*PI,xyboundary[0][0]/180*PI},{xyboundary[1][1]/180*PI,xyboundary[1][0]/180*PI}};

   //from xy coor to PHI
   for(int ip=0;ip<2;ip++){
      double ImageCoo[2]={xyboun[ip][1],xyboun[ip][0]};
      double PHI_in;
      GetPHI(zenith,azimuth,planephi,nz,ImageCoo,PHI_in);
      if(ip==0) PHI1=PHI_in;
      else PHI2=PHI_in;
   }
   if(PHI1>PHI2){
      double buff=PHI1;
      PHI1=PHI2;
      PHI2=buff;
   }
   //printf("xyboun={%lf,%lf} {%lf,%lf}\n",xyboun[0][1],xyboun[0][0],xyboun[1][1],xyboun[1][0]);
   if(jdebug>0) printf("WFCTAEvent::GetRange: PHI1=%lf PHI2=%lf PHIL=%lf\n",PHI1/PI*180,PHI2/PI*180,acos(cosPHIL)/PI*180);
   PHI1=TMath::Max(0.,PHI1);
   PHI2=TMath::Min(UsePhi?PI:acos(cosPHIL),PHI2);
   if(PHI2<=PHI1){
      if(jdebug>0) printf("WFCTAEvent::GetXYRange: return=-6, No Image from 0 to acos(PHIL)=%lf\n",acos(cosPHIL)/PI*180);
      return -6;
   }

   //the PHI range to have image ((x*cos(phi0)+y*sin(phi0))*cos(theta)+z*sin(theta)>=0)
   double PHIRange[2];
   bool inside1=CalPHIRange0(zenith,azimuth,planephi,nz,PHI1,PHIRange);
   bool inside2=CalPHIRange0(zenith,azimuth,planephi,nz,PHI2,PHIRange);
   if((!inside1)&&(!inside2)){
      if(jdebug>0) printf("WFCTAEvent::GetXYRange: return=-7, error between the input parameters, because of PHI12={%lf,%lf} PHIRange={%lf,%lf}\n",PHI1/PI*180,PHI2/PI*180,PHIRange[0]/PI*180,PHIRange[1]/PI*180);
      return -7;
   }
   else if(inside1&&(!inside2)){
      PHI2=CommonTools::ProcessAngle(PHIRange[1]);
   }
   else if((!inside1)&&inside2){
      PHI1=CommonTools::ProcessAngle(PHIRange[0]);
   }
   double PHI_low=IsLaser?PHI1:PHI2;
   double PHI_hig=IsLaser?PHI2:PHI1;
   PHI[0]=PHI_low;
   PHI[1]=PHI_hig;

   //from PHI to coor
   bool res1=GetImageXYCoo(zenith,azimuth,planephi,nz,PHI[0],XX[0],YY[0]);
   bool res2=GetImageXYCoo(zenith,azimuth,planephi,nz,PHI[1],XX[1],YY[1]);

   return 1;
}
bool WFCTAEvent::Calnz_phiL(double planephi,double indir[2],double &nz,double &PHIL){
   double incoo[3]={cos(planephi),sin(planephi),0};
   Getnz(incoo,indir,planephi,nz);
   PHIL=acos(sin(indir[0])*cos(indir[1]-planephi));
   return true;
}
double WFCTAEvent::CalTime(double ele_tel,double azi_tel,double incoo[3],double indir[2],double ImageCooXY[2],double time0){
   double planephi,nz;
   Getnz(incoo,indir,planephi,nz);
   double PHIL;
   Calnz_phiL(planephi,indir,nz,PHIL);
   double nz1,nz2;
   Getnz(nz,nz1,nz2);
   double xdir[3]={sin(ele_tel)*cos(planephi-azi_tel),sin(planephi-azi_tel),cos(ele_tel)*cos(planephi-azi_tel)};
   double ydir[3]={-nz1*sin(ele_tel)*sin(planephi-azi_tel)-nz2*cos(ele_tel),nz1*cos(planephi-azi_tel),-nz1*cos(ele_tel)*sin(planephi-azi_tel)+nz2**
sin(ele_tel)};
   double norm=sqrt(1+ImageCooXY[0]*ImageCooXY[0]+ImageCooXY[1]*ImageCooXY[1]);
   double dir0[3]={-ImageCooXY[0]/norm,-ImageCooXY[1]/norm,1./norm};
   double cosPHI=dir0[0]*xdir[0]+dir0[1]*xdir[1]+dir0[2]*xdir[2];
   double sinPHI=dir0[0]*ydir[0]+dir0[1]*ydir[1]+dir0[2]*ydir[2];
   double PHI=acos(cosPHI);
   if(sinPHI<0) PHI=2*PI-PHI;
   if(PHI>=PHIL) return -1;
   double length=sqrt(incoo[0]*incoo[0]+incoo[1]*incoo[1]+incoo[2]*incoo[2]);
   double result=(sinPHI+sin(PHIL))/(sin(PHIL)*cosPHI-sinPHI*cos(PHIL));
   return time0+(length/vlight*1.0e9)*result; //in ns
}
double WFCTAEvent::CalTime(double CC,double phi,double incoo[3],double nz,double PHIL,double ImageCooXY[2],double time0,double ele_in,double azi_ii
n){
   double lengthxy=sqrt(incoo[0]*incoo[0]+incoo[1]*incoo[1]);
   double planephi=acos(incoo[0]/lengthxy);
   if(incoo[1]<0) planephi=2*PI-planephi;
   double xdir[3]={cos(planephi),sin(planephi),0};
   double nz1,nz2;
   Getnz(nz,nz1,nz2);
   double ydir[3]={-sin(planephi)*nz1,cos(planephi)*nz1,nz2};
   double dir0[3];
   for(int ii=0;ii<3;ii++) dir0[ii]=cos(PHIL)*xdir[ii]+sin(PHIL)*ydir[ii];
   double laserdir[2];
   laserdir[0]=acos(dir0[2]/sqrt(dir0[0]*dir0[0]+dir0[1]*dir0[1]+dir0[2]*dir0[2]));
   laserdir[1]=acos(dir0[0]/sqrt(dir0[0]*dir0[0]+dir0[1]*dir0[1]));
   if(dir0[1]<0) laserdir[1]=2*PI-laserdir[1];

   double incoo2[3];
   CooCorr(incoo,laserdir,incoo2);

   double ele_tel,azi_tel;
   int nsol;
   double ele_out[4],azi_out[4];
   int signA[4];
   int isol=CalTelDir(CC,phi,planephi,nz,nsol,ele_out,azi_out,signA,ele_in,azi_in);
   if(isol<0) return -2;
   else{
       ele_tel=ele_out[isol];
       azi_tel=azi_out[isol];
   }
   double time=CalTime(ele_tel,azi_tel,incoo,laserdir,ImageCooXY,time0);
   return time;
}*/


double TelGeoFit::GetOutAngle(double zenith,double azimuth,double CC,double phi,double longcoo,double refplanephi){
   double planephi,nz;
   int signnz;
   CalPlane(zenith,azimuth,CC,phi,planephi,nz,signnz,refplanephi);
   double outdir[3];
   GetOutDir(zenith,azimuth,CC,phi,longcoo,outdir);
   return acos(cos(planephi)*outdir[0]+sin(planephi)*outdir[1]);
}
double TelGeoFit::GetOutAngle(double zenith,double azimuth,double CC,double phi,double longcoo,bool sign){
   double planephi,nz;
   int signnz;
   CalPlane(zenith,azimuth,CC,phi,planephi,nz,signnz,sign);
   double outdir[3];
   GetOutDir(zenith,azimuth,CC,phi,longcoo,outdir);
   return acos(cos(planephi)*outdir[0]+sin(planephi)*outdir[1]);
}
double TelGeoFit::GetOutAngle(double zenith,double azimuth,double planephi,double nz,double longcoo){
   double CC,phi;
   GetCCphi(zenith,azimuth,planephi,nz,CC,phi);
   double outdir[3];
   GetOutDir(zenith,azimuth,CC,phi,longcoo,outdir);
   return acos(cos(planephi)*outdir[0]+sin(planephi)*outdir[1]);
}
double TelGeoFit::GetLongCoo(double zenith,double azimuth,double planephi,double nz,double OutAngle){
   double xyzdir[3][3];
   CalPlane(planephi,nz,xyzdir);
   double outdir[3];
   for(int ii=0;ii<3;ii++){
      outdir[ii]=-cos(OutAngle)*xyzdir[0][ii]-sin(OutAngle)*xyzdir[1][ii];
   }
   double xx,yy;
   bool inside=GetImageCoo(zenith,azimuth,outdir,xx,yy);
   if(!inside) return InfNumber;
   else{
      double CC,phi;
      GetCCphi(zenith,azimuth,planephi,nz,CC,phi);
      return xx*cos(phi)+yy*sin(phi);
   }
}
double TelGeoFit::GetInjAngle(double telcoo[3],double incoo[3],double indir[2],double telzcoo[3]){
   double planephi,nz;
   bool exist=Getnz(telcoo,incoo,indir,planephi,nz,telzcoo);
   if(!exist) return InfNumber;
   double injdir[3]={telzcoo[0]-telcoo[0],telzcoo[1]-telcoo[1],telzcoo[2]-telcoo[2]};
   double norm=sqrt(injdir[0]*injdir[0]+injdir[1]*injdir[1]+injdir[2]*injdir[2]);
   if(norm==0) return InfNumber;
   double outdir[3]={sin(indir[0])*cos(indir[1]),sin(indir[0])*sin(indir[1]),cos(indir[0])};
   if(outdir[2]>0){
      outdir[0]*=-1;
      outdir[1]*=-1;
      outdir[2]*=-1;
   }
   return acos((injdir[0]*outdir[0]+injdir[1]*outdir[1]+injdir[2]*outdir[2])/norm);
}
double TelGeoFit::GetOutAngleFromLength(double telcoo[3],double incoo[3],double indir[2],double length){
   double telzcoo[3];
   double planephi,nz;
   bool exist=Getnz(telcoo,incoo,indir,planephi,nz,telzcoo);
   if(!exist) return InfNumber;
   double outdir[3]={sin(indir[0])*cos(indir[1]),sin(indir[0])*sin(indir[1]),cos(indir[0])};
   double outcoo[3];
   for(int ii=0;ii<3;ii++){
      outcoo[ii]=incoo[ii]+outdir[ii]*length;
   }
   double dir0[3]={outcoo[0]-telcoo[0],outcoo[1]-telcoo[1],outcoo[2]-telcoo[2]};
   double norm0=0,norm1=0;
   double res=0;
   for(int ii=0;ii<3;ii++){
      norm0+=dir0[ii]*dir0[ii];
      norm1+=telzcoo[ii]*telzcoo[ii];
      res+=dir0[ii]*telzcoo[ii];
   }
   if(norm0==0||norm1==0) return InfNumber;
   res=acos(res/sqrt(norm0*norm1));
   return dir0[2]>0?res:(-res);
}
double TelGeoFit::GetOutAngleFromHeight(double telcoo[3],double incoo[3],double indir[2],double height){
   double telzcoo[3];
   double planephi,nz;
   bool exist=Getnz(telcoo,incoo,indir,planephi,nz,telzcoo);
   if(!exist) return InfNumber;
   double outdir[3]={sin(indir[0])*cos(indir[1]),sin(indir[0])*sin(indir[1]),cos(indir[0])};
   if(outdir[2]==0) return InfNumber;
   double dz=(telcoo[2]+height-incoo[2])/outdir[2];
   double outcoo[3];
   for(int ii=0;ii<3;ii++){
      outcoo[ii]=incoo[ii]+outdir[ii]*dz;
   }
   double dir0[3]={outcoo[0]-telcoo[0],outcoo[1]-telcoo[1],outcoo[2]-telcoo[2]};
   double norm0=0,norm1=0;
   double res=0;
   for(int ii=0;ii<3;ii++){
      norm0+=dir0[ii]*dir0[ii];
      norm1+=telzcoo[ii]*telzcoo[ii];
      res+=dir0[ii]*telzcoo[ii];
   }
   if(norm0==0||norm1==0) return InfNumber;
   res=acos(res/sqrt(norm0*norm1));
   return dir0[2]>0?res:(-res);
}
double TelGeoFit::GetLongCooFromLength(double zenith,double azimuth,double telcoo[3],double incoo[3],double indir[2],double length){
   double outangle=GetOutAngleFromLength(telcoo,incoo,indir,length);
   double planephi,nz;
   double telzcoo[3];
   bool exist=Getnz(telcoo,incoo,indir,planephi,nz,telzcoo);
   if(!exist) return InfNumber;
   return GetLongCoo(zenith,azimuth,planephi,nz,outangle);
}
double TelGeoFit::GetLongCooFromHeight(double zenith,double azimuth,double telcoo[3],double incoo[3],double indir[2],double height){
   double outangle=GetOutAngleFromHeight(telcoo,incoo,indir,height);
   double planephi,nz;
   double telzcoo[3];
   bool exist=Getnz(telcoo,incoo,indir,planephi,nz,telzcoo);
   if(!exist) return InfNumber;
   return GetLongCoo(zenith,azimuth,planephi,nz,outangle);
}
double TelGeoFit::GetScatAngleFromLength(double telcoo[3],double incoo[3],double indir[2],double length){
   double telzcoo[3];
   double injangle=GetInjAngle(telcoo,incoo,indir,telzcoo);
   double outangle=GetOutAngleFromLength(telcoo,incoo,indir,length);
   if(injangle<=InfNumber||outangle<=InfNumber) return InfNumber;
   double scatangle=fabs(injangle)+fabs(outangle);
   if(scatangle>=PI) return InfNumber;
   else return scatangle;
}
double TelGeoFit::GetScatAngleFromHeight(double telcoo[3],double incoo[3],double indir[2],double height){
   double telzcoo[3];
   double injangle=GetInjAngle(telcoo,incoo,indir,telzcoo);
   double outangle=GetOutAngleFromHeight(telcoo,incoo,indir,height);
   if(injangle<=InfNumber||outangle<=InfNumber) return InfNumber;
   double scatangle=fabs(injangle)+fabs(outangle);
   if(scatangle>=PI) return InfNumber;
   else return scatangle;
}
double TelGeoFit::GetScatAngleFromLongcoo(double zenith,double azimuth,double telcoo[3],double incoo[3],double indir[2],double longcoo){
   double telzcoo[3];
   double injangle=GetInjAngle(telcoo,incoo,indir,telzcoo);
   double planephi,nz;
   bool exist=Getnz(telcoo,incoo,indir,planephi,nz,telzcoo);
   if(!exist) return InfNumber;
   double outangle=GetOutAngle(zenith,azimuth,planephi,nz,longcoo);
   if(injangle<=InfNumber||outangle<=InfNumber) return InfNumber;
   double scatangle=fabs(injangle)+fabs(outangle);
   if(scatangle>=PI) return InfNumber;
   else return scatangle;
}
double TelGeoFit::GetGroundLength(double telcoo[3],double incoo[3],double indir[2]){
   double telzcoo[3];
   double injangle=GetInjAngle(telcoo,incoo,indir,telzcoo);
   if(injangle<=InfNumber) return sqrt(pow(incoo[0]-telcoo[0],2)+pow(incoo[1]-telcoo[1],2)+pow(incoo[2]-telcoo[2],2));
   else return sqrt(pow(telzcoo[0]-telcoo[0],2)+pow(telzcoo[1]-telcoo[1],2)+pow(telzcoo[2]-telcoo[2],2));
}
double TelGeoFit::GetRp(double telcoo[3],double incoo[3],double indir[2],double &length){
   double telzcoo[3];
   double injangle=GetInjAngle(telcoo,incoo,indir,telzcoo);
   if(injangle<=InfNumber) return InfNumber;
   double length0=GetGroundLength(telcoo,incoo,indir);
   double Rp=length0*sin(injangle);
   double len0=fabs(length0*cos(injangle));
   double len1=sqrt(pow(telzcoo[0]-incoo[0],2)+pow(telzcoo[1]-incoo[1],2)+pow(telzcoo[2]-incoo[2],2));
   if((incoo[2]-telzcoo[0])<0) length=len0+len1;
   else length=len0-len1;
   return Rp;
}
double TelGeoFit::GetTime(double outangle,double injangle,double groundlength,double delta_length){
   if(fabs(outangle)+fabs(injangle)>PI) return InfNumber;
   double Rp=groundlength*sin(injangle);
   double scat_angle=(PI-fabs(outangle)-fabs(injangle));
   //double length1=groundlength*cos(injangle)+Rp/tan(scat_angle/2); //go up
   double length1=groundlength*cos(injangle)+Rp/sin(scat_angle)*(1+cos(scat_angle)); //go up
   double length2=0+Rp*tan(scat_angle/2); //go down
   if(jdebug>9) printf("TelGeoFit::GetTime: outangle=%.2lf injangle=%.2lf scatangle=%.2lf groundlength=%.0lf delta_length=%.0lf Rp=%.0lf length={%.0lf,%.0lf}\n",outangle/PI*180,injangle/PI*180,scat_angle/PI*180,groundlength,delta_length,Rp,length1-delta_length,length2);
   return GoDown?(length2):(length1-delta_length);
}
double TelGeoFit::GetTime(double zenith,double azimuth,double CC,double phi,double longcoo,double groundlength,double delta_length,double injangle,bool sign){
   double outangle=GetOutAngle(zenith,azimuth,CC,phi,longcoo,sign);
   return GetTime(outangle,injangle,groundlength,delta_length);
}
double TelGeoFit::GetTimeFromLongcoo(double zenith,double azimuth,double telcoo[3],double incoo[3],double indir[2],double longcoo){
   double planephi,nz;
   double telzcoo[3];
   Getnz(telcoo,incoo,indir,planephi,nz,telzcoo);
   double outangle=GetOutAngle(zenith,azimuth,planephi,nz,longcoo);
   double injangle=GetInjAngle(telcoo,incoo,indir,telzcoo);
   double groundlength=GetGroundLength(telcoo,incoo,indir);
   double delta_length=sqrt(pow(incoo[0]-telzcoo[0],2)+pow(incoo[1]-telzcoo[1],2)+pow(incoo[2]-telzcoo[2],2));
   if(incoo[2]<telzcoo[2]) delta_length*=-1;
   return GetTime(outangle,injangle,groundlength,delta_length);
}
double TelGeoFit::GetTimeFromLength(double zenith,double azimuth,double telcoo[3],double incoo[3],double indir[2],double length,double &longcoo){
   longcoo=GetLongCooFromLength(zenith,azimuth,telcoo,incoo,indir,length);
   double outangle=GetOutAngleFromLength(telcoo,incoo,indir,length);
   double telzcoo[3];
   double injangle=GetInjAngle(telcoo,incoo,indir,telzcoo);
   double groundlength=GetGroundLength(telcoo,incoo,indir);
   double delta_length=sqrt(pow(incoo[0]-telzcoo[0],2)+pow(incoo[1]-telzcoo[1],2)+pow(incoo[2]-telzcoo[2],2));
   if(incoo[2]<telzcoo[2]) delta_length*=-1;
   return GetTime(outangle,injangle,groundlength,delta_length);
}
double TelGeoFit::GetTimeFromHeight(double zenith,double azimuth,double telcoo[3],double incoo[3],double indir[2],double height,double &longcoo){
   longcoo=GetLongCooFromHeight(zenith,azimuth,telcoo,incoo,indir,height);
   double outangle=GetOutAngleFromHeight(telcoo,incoo,indir,height);
   double telzcoo[3];
   double injangle=GetInjAngle(telcoo,incoo,indir,telzcoo);
   double groundlength=GetGroundLength(telcoo,incoo,indir);
   double delta_length=sqrt(pow(incoo[0]-telzcoo[0],2)+pow(incoo[1]-telzcoo[1],2)+pow(incoo[2]-telzcoo[2],2));
   if(incoo[2]<telzcoo[2]) delta_length*=-1;
   return GetTime(outangle,injangle,groundlength,delta_length);
}

double TelGeoFit::GetRotPos(int icoo,int Li){
   if(icoo<0||icoo>2) return InfNumber;
   if(Li==2){
      double ret_pos[3]={-58606,-12352,2715};
      return ret_pos[icoo];
   }
   else if(Li==3){
      double ret_pos[3]={-45006,99008,1107};
      return ret_pos[icoo];
   }
   else return InfNumber;
}
double TelGeoFit::GetRotDir(int ipar,int Li){
   if(ipar<0||ipar>3) return InfNumber;
   if(Li==2){
      //double ret_dir[4]={88.72/180*PI,14.81/180*PI,-4.10/180*PI,-1.22/180*PI};
      double ret_dir[4]={87.19/180*PI,17.39/180*PI,-3.50/180*PI,-1.23/180*PI};
      return ret_dir[ipar];
   }
   else if(Li==3){
      //double ret_dir[4]={89.10/180*PI,117.31/180*PI,0.62/180*PI,1.45/180*PI};
      double ret_dir[4]={88.72/180*PI,162.59/180*PI,1.96/180*PI,1.56/180*PI};
      return ret_dir[ipar];
   }
   else{
      double ret_dir[4]={90./180*PI,0./180*PI,0./180*PI,0./180*PI};
      return ret_dir[ipar];
   }
}
double TelGeoFit::GetRotTime(int Li){
   if(Li==2) return 990015724;
   else if(Li==3) return 990820158;
   else return InfNumber;
}
void TelGeoFit::SetRotTelPars(bool IsInit,bool IsFit){
   for(int irot=0;irot<NRotMax;irot++){
      if(IsInit){
         for(int ii=0;ii<3;ii++) RotPos_int[irot][ii]=rotpos[irot][ii];
         //RotPos_int[irot][0]=GetRotPos(0,irot+1);
         //RotPos_int[irot][1]=GetRotPos(1,irot+1);
         //RotPos_int[irot][2]=GetRotPos(2,irot+1);
         RotEle_int[irot]=GetRotDir(0,irot+1);
         RotAzi_int[irot]=GetRotDir(1,irot+1);
         RefEle_int[irot]=GetRotDir(2,irot+1);
         RefAzi_int[irot]=GetRotDir(3,irot+1);
         RotTime_int[irot]=GetRotTime(irot+1);
      }
      if(IsFit){
         for(int ii=0;ii<3;ii++) RotPos_fit[irot][ii]=rotpos[irot][ii];
         //RotPos_fit[irot][0]=GetRotPos(0,irot+1);
         //RotPos_fit[irot][1]=GetRotPos(1,irot+1);
         //RotPos_fit[irot][2]=GetRotPos(2,irot+1);
         RotEle_fit[irot]=GetRotDir(0,irot+1);
         RotAzi_fit[irot]=GetRotDir(1,irot+1);
         RefEle_fit[irot]=GetRotDir(2,irot+1);
         RefAzi_fit[irot]=GetRotDir(3,irot+1);
         RotTime_fit[irot]=GetRotTime(irot+1);
      }
   }
   if(IsInit){
      for(int itel=0;itel<NCTMax;itel++){
         TelZen_int[itel]=teldir0[itel][0];
         TelAzi_int[itel]=teldir0[itel][1];
      }
   }
   if(IsFit){
      for(int itel=0;itel<NCTMax;itel++){
         TelZen_fit[itel]=teldir0[itel][0];
         TelAzi_fit[itel]=teldir0[itel][1];
      }
   }
}
void TelGeoFit::SetRotPars(int Li,const double* par,bool IsInit,bool IsFit){
   int irot=Li-1;
   if(irot<0||irot>=NRotMax) return;
   if(IsInit){
      for(int ii=0;ii<3;ii++) RotPos_int[irot][ii]=par[ii];
      RotEle_int[irot]=par[0+3];
      RotAzi_int[irot]=par[1+3];
      RefEle_int[irot]=par[2+3];
      RefAzi_int[irot]=par[3+3];
      RotTime_int[irot]=par[4+3];
   }
   if(IsFit){
      for(int ii=0;ii<3;ii++) RotPos_fit[irot][ii]=par[ii];
      RotEle_fit[irot]=par[0+3];
      RotAzi_fit[irot]=par[1+3];
      RefEle_fit[irot]=par[2+3];
      RefAzi_fit[irot]=par[3+3];
      RotTime_fit[irot]=par[4+3];
   }
}
void TelGeoFit::SetTelPars(int iTel,const double* par,bool IsInit,bool IsFit){
   int itel=iTel-1;
   if(itel<0||itel>=NCTMax) return;
   if(IsInit){
      TelZen_int[itel]=par[0];
      TelAzi_int[itel]=par[1];
   }
   if(IsFit){
      TelZen_fit[itel]=par[0];
      TelAzi_fit[itel]=par[1];
   }
}
bool TelGeoFit::CalRotateZeroPos(double ele_rotate,double azi_rotate,double ele_ref,double azi_ref,double &ele0,double &azi0){
   double sintheta0=cos(ele_ref)*cos(azi_ref-azi_rotate)*cos(ele_rotate)+sin(ele_ref)*sin(ele_rotate);
   ele0=asin(sintheta0);

   double p1=cos(ele_ref)*cos(azi_ref-azi_rotate)*sin(ele_rotate)-sin(ele_ref)*cos(ele_rotate);
   double p2=cos(ele_ref)*sin(azi_ref-azi_rotate);
   double norm=sqrt(p1*p1+p2*p2);
   double margin=1.0e-5;
   if(fabs(norm)<margin){
      azi0=0.;
   }
   else{
      azi0=acos(p1/norm);
      if(p2<0) azi0=2*PI-azi0;
   }
   return true;
}
bool TelGeoFit::CalDir_out(double ele_rotate,double azi_rotate,double ele_ref,double azi_ref,double ele_in,double azi_in,double &ele_out,double &azi_out){
   double ele0,azi0;
   if(!CalRotateZeroPos(ele_rotate,azi_rotate,ele_ref,azi_ref,ele0,azi0)) return false;
   double p1=cos(ele0+ele_in)*cos(azi0+azi_in)*sin(ele_rotate)+sin(ele0+ele_in)*cos(ele_rotate);
   double p2=cos(ele0+ele_in)*sin(azi0+azi_in);
   double xx=p1*cos(azi_rotate)-p2*sin(azi_rotate);
   double yy=p1*sin(azi_rotate)+p2*cos(azi_rotate);
   double zz=-cos(ele0+ele_in)*cos(azi0+azi_in)*cos(ele_rotate)+sin(ele0+ele_in)*sin(ele_rotate);
   double norm=sqrt(xx*xx+yy*yy+zz*zz);
   ele_out=asin(zz/norm);
   azi_out=acos(xx/sqrt(xx*xx+yy*yy));
   if(yy<0) azi_out=2*PI-azi_out;
   return true;
}
bool TelGeoFit::CalDir_out(double ele_in,double azi_in,int Li,double &ele_out,double &azi_out){
   double ele_rotate=GetRotDir(0,Li);
   double azi_rotate=GetRotDir(1,Li);
   double ele_ref=GetRotDir(2,Li);
   double azi_ref=GetRotDir(3,Li);
   if(ele_rotate==InfNumber) return false;
   if(azi_rotate==InfNumber) return false;
   if(ele_ref==InfNumber) return false;
   if(azi_ref==InfNumber) return false;
   return CalDir_out(ele_rotate,azi_rotate,ele_ref,azi_ref,ele_in,azi_in,ele_out,azi_out);
}

int TelGeoFit::GetRotIndex(int ievt){
   if(ievt<0||ievt>=nevent) return -1;
   double rttime=rabbitTime[ievt]-(int)(rabbitTime[ievt]);
   return RotateDB::GetLi(rttime/2.0e-8);
}
double TelGeoFit::GetChi2XY(int &ndof,int ievt,double CC,double phi){
   double chi2=0;
   ndof=0;
   if(ievt<0||ievt>=nevent) return chi2;
   if(false){
      ///from CC and phi
      ndof=2;
      chi2=pow(fabs(CC-image_pars[ievt][0][0])/image_pars[ievt][0][1],2)+pow(fabs(phi-image_pars[ievt][1][0])/image_pars[ievt][1][1],2);
      //chi2=pow(fabs(CC-image_pars[ievt][0][0])/image_pars[ievt][0][1],2);
      //if(jdebug>2) printf("TelGeoFit::GetChi2XY: ievt=%d chi2=%.1lf ndof=%d\n",ievt,chi2,ndof);
   }
   else{
      ///from image
      for(int isipm=0;isipm<MAXPMT;isipm++){
         if(Npe_sipm[ievt][isipm]<=0) continue;
         double ImageX,ImageY;
         double dcell=TelGeoFit::GetImageXYCoo(isipm,ImageX,ImageY,-1,false);
         if(dcell<0) continue;
         double distance=(ImageX*sin(phi)-ImageY*cos(phi)+CC);
         double probi=Npe_sipm[ievt][isipm]/Npe_sum[ievt];
         chi2+=probi*pow(distance/(dcell/2.),2);
         ndof++;
         if(jdebug>2) printf("TelGeoFit::GetChi2XY: ievt=%d sipm=%d chi2=%.1lf ndof=%d\n",ievt,isipm,chi2,ndof);
      }
   }
   if(jdebug>1) printf("TelGeoFit::GetChi2XY: ievt=%d chi2=%.1lf ndof=%d\n",ievt,chi2,ndof);
   XYchi2_evt[ievt]=chi2;
   XYndof_evt[ievt]=ndof;
   return chi2;
}
double TelGeoFit::GetChi2Time(int &ndof,int ievt,double phi,double time0,double zenith,double azimuth,double telcoo[3],double incoo[3],double indir[2]){
   double chi2=0;
   ndof=0;
   if(ievt<0||ievt>=nevent) return chi2;
   double rabbittime=(rabbitTime[ievt]-((int)rabbitTime[ievt]))/2.0e-8;
   //int irot=RotateDB::GetLi(rabbittime);
   //if(irot<0) return chi2;
   //if(RotateDB::rotindex[irot]!=FixLi) return chi2;

   double time_ref=0;

   for(int isipm=0;isipm<MAXPMT;isipm++){
      if(Npe_sipm[ievt][isipm]<=0) continue;
      double ImageX,ImageY;
      double dcell=TelGeoFit::GetImageXYCoo(isipm,ImageX,ImageY,-1,false);
      if(dcell<0) continue;
      double datatime=(rabbittime*20-time0)+(Time_sipm[ievt][isipm]-time_ref)*80;
      double longcoo=ImageX*cos(phi)+ImageY*sin(phi);
      double modeltime=GetTimeFromLongcoo(zenith,azimuth,telcoo,incoo,indir,longcoo)/29.98;
      double probi=Npe_sipm[ievt][isipm]/Npe_sum[ievt];
      chi2+=probi*pow((datatime-modeltime)/80.,2);
      ndof++;
      if(jdebug>2) printf("TelGeoFit::GetChi2Time: ievt=%d sipm=%d time0=%9.0lf chi2=%.1lf ndof=%d\n",ievt,isipm,time0,chi2,ndof);
   }
   if(jdebug>1) printf("TelGeoFit::GetChi2Time: ievt=%d time0=%9.0lf chi2=%.1lf ndof=%d\n",ievt,time0,chi2,ndof);
   Tchi2_evt[ievt]=chi2;
   Tndof_evt[ievt]=ndof;
   return chi2;
}
double TelGeoFit::Interface(const double* par){
   int nrot=par[0];
   int rotstart=1;
   const int npars_rot=9;
   bool rotindex_fit[NRotMax];
   for(int ii=0;ii<NRotMax;ii++) rotindex_fit[ii]=false;
   ///8 pararmeters for rotate. 0-2: coordinate; 3:elevation of rotate axis; 4:azimuth of rotate axis; 5:bias of elevation zero; 6;bias of azimuth zero; 7: start time
   double rotpar_fit[NRotMax][npars_rot-1];
   for(int irot=0;irot<nrot;irot++){
      int iRot=(int)(par[irot*npars_rot+0+rotstart]+0.5);
      if(iRot>=1&&iRot<=NRotMax) rotindex_fit[iRot-1]=true;
      for(int ii=0;ii<npars_rot-1;ii++) rotpar_fit[iRot-1][ii]=par[irot*npars_rot+ii+1+rotstart];
      if(jdebug>1){
         bool nfixcoox=((FixRotPos[0]&(1<<(iRot-1)))==0);
         bool nfixcooy=((FixRotPos[1]&(1<<(iRot-1)))==0);
         bool nfixcooz=((FixRotPos[2]&(1<<(iRot-1)))==0);
         if( nfixcoox||nfixcooy||nfixcooz ) printf("TelGeoFit::Interface: pars: iRot=%d coo={%.2lf,%.2lf,%.2lf}\n",iRot,rotpar_fit[iRot-1][0],rotpar_fit[iRot-1][1],rotpar_fit[iRot-1][2]);

         bool nfixaxiszen=((FixRotDir[0]&(1<<(iRot-1)))==0);
         bool nfixaxisazi=((FixRotDir[1]&(1<<(iRot-1)))==0);
         bool nfixzen_bias=((FixRotDir[2]&(1<<(iRot-1)))==0);
         bool nfixazi_bias=((FixRotDir[3]&(1<<(iRot-1)))==0);
         if( (nfixaxiszen||nfixaxisazi) || (nfixzen_bias||nfixazi_bias) ) printf("TelGeoFit::Interface: pars: iRot=%d dir={%.2lf,%.2lf,%.2lf,%.2lf}\n",iRot,rotpar_fit[iRot-1][3]/PI*180,rotpar_fit[iRot-1][4]/PI*180,rotpar_fit[iRot-1][5]/PI*180,rotpar_fit[iRot-1][6]/PI*180);

         bool nfixtime=((FixRotTime&(1<<(iRot-1)))==0);
         if( nfixtime )  printf("TelGeoFit::Interface: pars: iRot=%d time=%9.0lf\n",iRot,rotpar_fit[iRot-1][7]);
      }
   }

   int ntel=(int)(par[nrot*npars_rot+rotstart]+0.5);
   int telstart=nrot*npars_rot+1+rotstart;
   const int npars_tel=3;
   bool telindex_fit[NCTMax];
   for(int ii=0;ii<NCTMax;ii++) telindex_fit[ii]=false;
   ///2 pararmeters for telescope. 0: zenith; 1: azimuth;
   double teldir_fit[NCTMax][2];
   for(int itel=0;itel<ntel;itel++){
      int iTel=(int)(par[itel*3+0+telstart]+0.5);
      if(iTel>=1&&iTel<=NCTMax) telindex_fit[iTel-1]=true;
      teldir_fit[iTel-1][0]=par[itel*3+1+telstart];
      teldir_fit[iTel-1][1]=par[itel*3+2+telstart];
      if(jdebug>1){
         bool nfixzen=((FixTelDir[0]&(1<<(iTel-1)))==0);
         bool nfixazi=((FixTelDir[1]&(1<<(iTel-1)))==0);
         if(nfixzen||nfixazi) printf("TelGeoFit::Interface: pars: iTel=%d dir={%.2lf,%.2lf}\n",iTel,teldir_fit[iTel-1][0]/PI*180,teldir_fit[iTel-1][1]/PI*180);
      }
   }

   double Chi2=0;
   int Ndof=0;
   for(int ievt=0;ievt<nevent;ievt++){
      int iTel=telindex[ievt];
      int itel=iTel-1;
      if(iTel<0||itel>=NCTMax) continue;
      if(!telindex_fit[itel]) continue;
      int irot=GetRotIndex(ievt);
      if(irot<0||irot>=NRotMax) continue;
      int iRot=RotateDB::rotindex[irot];
      if(!rotindex_fit[iRot-1]) continue;

      if(FixTime>1300000000&&(FixTime!=((int)rabbitTime[ievt]))) continue;
      if(replaced[ievt]) continue;
      if(FixLi>0&&FixLi!=iRot) continue;

      double rotdir_fit[2];
      double ele_in=PI/2-rotdir0[ievt][0]; //from zenith to elevation
      double azi_in=rotdir0[ievt][1];
      double rotele=rotpar_fit[iRot-1][3];
      double rotazi=rotpar_fit[iRot-1][4];
      double rotrefele_bias=rotpar_fit[iRot-1][5];
      double rotrefazi_bias=rotpar_fit[iRot-1][6];
      CalDir_out(rotele,rotazi,rotrefele_bias,rotrefazi_bias,ele_in,azi_in,rotdir_fit[0],rotdir_fit[1]);
      rotdir_fit[0]=PI/2-rotdir_fit[0]; //from elevation to zenith

      double rotpos0[3]={rotpar_fit[iRot-1][0],rotpar_fit[iRot-1][1],rotpar_fit[iRot-1][2]};
      double CC,phi;
      GetCCphi(teldir_fit[itel][0],teldir_fit[itel][1],telpos[itel],rotpos0,rotdir_fit,CC,phi);

      int xyndof=0;
      double xychi2=GetChi2XY(xyndof,ievt,CC,phi);
      Chi2+=xychi2;
      Ndof+=xyndof;
      
      int tndof=0,tchi2=0;
      bool nfixrottime=((FixRotTime&(1<<(iRot-1)))==0);
      if(nfixrottime){
         double time0=rotpar_fit[iRot-1][7];
         tchi2=GetChi2Time(tndof,ievt,phi,time0,teldir_fit[itel][0],teldir_fit[itel][1],telpos[itel],rotpos0,rotdir_fit);
      }
      Chi2+=tchi2;
      Ndof+=tndof;

      if(jdebug>1){
         bool nfixaxiszen=((FixRotDir[0]&(1<<(iRot-1)))==0);
         bool nfixaxisazi=((FixRotDir[1]&(1<<(iRot-1)))==0);
         bool nfixzen_bias=((FixRotDir[2]&(1<<(iRot-1)))==0);
         bool nfixazi_bias=((FixRotDir[3]&(1<<(iRot-1)))==0);
         bool nfixdir=( (nfixaxiszen||nfixaxisazi) || (nfixzen_bias||nfixazi_bias) );

         if(nfixdir) printf("TelGeoFit::Interface: ievt=%d rotzen_in=%.2lf rotazi_in=%.2lf rotzen_out=%.2lf rotazi_out=%.2lf pars={%.2lf,%.2lf,%.2lf,%.2lf}\n",ievt,rotdir0[ievt][0]/PI*180,rotdir0[ievt][1]/PI*180,rotdir_fit[0]/PI*180,rotdir_fit[1]/PI*180,rotele/PI*180,rotazi/PI*180,rotrefele_bias/PI*180,rotrefazi_bias/PI*180);
         printf("TelGeoFit::Interface: ievt=%d(%d) itel=%d telzen=%.2lf telazi=%.2lf cc=%.2lf phi=%.2lf image_pars={%.2lf,%.2lf}\n",ievt,nevent,itel+1,teldir_fit[itel][0]/PI*180,teldir_fit[itel][1]/PI*180,CC/PI*180,phi/PI*180,image_pars[ievt][0][0]/PI*180,image_pars[ievt][1][0]/PI*180);
      }
   }
   if(jdebug>0) printf("TelGeoFit::Interface: Chi2=%.1lf NDof=%d\n",Chi2,Ndof);
   Chi2_All=Chi2;
   Ndof_All=Ndof;
   return Chi2;
}

bool TelGeoFit::DoFit(int ntel,int* tellist,int nrot,int* rotlist,bool force){
   if(!tellist) return false;
   if(!rotlist) return false;
   for(int itel=0;itel<ntel;itel++){
      if(tellist[itel]<1||tellist[itel]>NCTMax) return false;
   }
   for(int irot=0;irot<nrot;irot++){
      if(rotlist[irot]<1||rotlist[irot]>NRotMax) return false;
   }
   if(!minimizer) minimizer=ROOT::Math::Factory::CreateMinimizer("Minuit","Migrad");
   minimizer->Clear();
   minimizer->SetMaxFunctionCalls(1000000);
   minimizer->SetMaxIterations(100000);
   minimizer->SetTolerance(0.001);
   minimizer->SetPrintLevel(0);
   int rotstart=1;
   const int npars_rot=9;
   int telstart=nrot*npars_rot+1+rotstart;
   const int npars_tel=3;
   int ntotal=telstart+ntel*npars_tel;
   #if defined(__CINT__)
   ROOT::Math::Functor f(this,"TelGeoFit","Interface");
   #else
   ROOT::Math::Functor f(this,&TelGeoFit::Interface,ntotal);
   #endif
   minimizer->SetFunction(f);

   double rotpos_margin=5.*1.0e2;
   double rotele_margin=5./180*PI;
   double rotazi_margin=180./180*PI;
   double rotrefele_margin=5./180*PI;
   double rotrefazi_margin=5./180*PI;
   double time_margin=10000.;
   minimizer->SetFixedVariable(0,"NRot",nrot);
   for(int irot=0;irot<nrot;irot++){
      int iRot=rotlist[irot]-1;
      minimizer->SetFixedVariable(rotstart+irot*npars_rot+0,Form("Rot%d",iRot+1),iRot+1.);

      bool nfixcoox=((FixRotPos[0]&(1<<(iRot)))==0);
      bool nfixcooy=((FixRotPos[1]&(1<<(iRot)))==0);
      bool nfixcooz=((FixRotPos[2]&(1<<(iRot)))==0);

      if(!nfixcoox) minimizer->SetFixedVariable(rotstart+irot*npars_rot+1,Form("RotPosx_%d",iRot+1),RotPos_int[iRot][0]);
      else minimizer->SetLimitedVariable(rotstart+irot*npars_rot+1,Form("RotPosx_%d",iRot+1),RotPos_int[iRot][0],50.,RotPos_int[iRot][0]-rotpos_margin,RotPos_int[iRot][0]+rotpos_margin);
      if(!nfixcooy) minimizer->SetFixedVariable(rotstart+irot*npars_rot+2,Form("RotPosy_%d",iRot+1),RotPos_int[iRot][1]);
      else minimizer->SetLimitedVariable(rotstart+irot*npars_rot+2,Form("RotPosx_%d",iRot+1),RotPos_int[iRot][1],50.,RotPos_int[iRot][1]-rotpos_margin,RotPos_int[iRot][1]+rotpos_margin);
      if(!nfixcooz) minimizer->SetFixedVariable(rotstart+irot*npars_rot+3,Form("RotPosz_%d",iRot+1),RotPos_int[iRot][2]);
      else minimizer->SetLimitedVariable(rotstart+irot*npars_rot+3,Form("RotPosx_%d",iRot+1),RotPos_int[iRot][2],50.,RotPos_int[iRot][2]-rotpos_margin,RotPos_int[iRot][2]+rotpos_margin);

      bool nfixaxiszen=((FixRotDir[0]&(1<<(iRot)))==0);
      bool nfixaxisazi=((FixRotDir[1]&(1<<(iRot)))==0);
      bool nfixzen_bias=((FixRotDir[2]&(1<<(iRot)))==0);
      bool nfixazi_bias=((FixRotDir[3]&(1<<(iRot)))==0);

      if(!nfixaxiszen) minimizer->SetFixedVariable(rotstart+irot*npars_rot+4,Form("RotEle_%d",iRot+1),RotEle_int[iRot]);
      else minimizer->SetLimitedVariable(rotstart+irot*npars_rot+4,Form("RotEle_%d",iRot+1),RotEle_int[iRot],fabs(0.5/180.*PI),TMath::Max(0.,RotEle_int[iRot]-rotele_margin),TMath::Min(PI/2.,RotEle_int[iRot]+rotele_margin));
      if(!nfixaxisazi) minimizer->SetFixedVariable(rotstart+irot*npars_rot+5,Form("RotAzi_%d",iRot+1),RotAzi_int[iRot]);
      else minimizer->SetLimitedVariable(rotstart+irot*npars_rot+5,Form("RotAzi_%d",iRot+1),RotAzi_int[iRot],fabs(1./180.*PI),RotAzi_int[iRot]-rotazi_margin,RotAzi_int[iRot]+rotazi_margin);
      if(!nfixzen_bias) minimizer->SetFixedVariable(rotstart+irot*npars_rot+6,Form("RotRefEle_bias_%d",iRot+1),RefEle_int[iRot]);
      else minimizer->SetLimitedVariable(rotstart+irot*npars_rot+6,Form("RotRefEle_bias_%d",iRot+1),RefEle_int[iRot],fabs(0.5/180.*PI),RefEle_int[iRot]-rotrefele_margin,RefEle_int[iRot]+rotrefele_margin);
      if(!nfixazi_bias) minimizer->SetFixedVariable(rotstart+irot*npars_rot+7,Form("RotRefAzi_bias_%d",iRot+1),RefAzi_int[iRot]);
      else minimizer->SetLimitedVariable(rotstart+irot*npars_rot+7,Form("RotRefAzi_bias_%d",iRot+1),RefAzi_int[iRot],fabs(0.5/180.*PI),RefAzi_int[iRot]-rotrefazi_margin,RefAzi_int[iRot]+rotrefazi_margin);

      bool nfixtime=((FixRotTime&(1<<(iRot)))==0);

      if(!nfixtime) minimizer->SetFixedVariable(rotstart+irot*npars_rot+8,Form("Time0_%d",iRot+1),RotTime_int[iRot]);
      else minimizer->SetLimitedVariable(rotstart+irot*npars_rot+8,Form("Time0_%d",iRot+1),RotTime_int[iRot],100,RotTime_int[iRot]-time_margin,RotTime_int[iRot]+time_margin);
   }

   minimizer->SetFixedVariable(telstart-1,"NTel",ntel);

   double telele_margin=2./180*PI;
   double telazi_margin=2./180*PI;
   for(int itel=0;itel<ntel;itel++){
      int iTel=tellist[itel]-1;
      minimizer->SetFixedVariable(itel*npars_tel+0+telstart,Form("Tel%d",iTel),iTel+1.);

      bool nfixzen=((FixTelDir[0]&(1<<(iTel)))==0);
      bool nfixazi=((FixTelDir[1]&(1<<(iTel)))==0);

      if(!nfixzen) minimizer->SetFixedVariable(itel*npars_tel+1+telstart,Form("Zen%d",iTel),TelZen_int[iTel]);
      else minimizer->SetLimitedVariable(itel*npars_tel+1+telstart,Form("Zen%d",iTel),TelZen_int[iTel],fabs(0.1/180.*PI),TMath::Max(0.,TelZen_int[iTel]-telele_margin),TMath::Min(PI/2,TelZen_int[iTel]+telele_margin));
      if(!nfixazi) minimizer->SetFixedVariable(itel*npars_tel+2+telstart,Form("Azi%d",iTel),TelAzi_int[iTel]);
      else minimizer->SetLimitedVariable(itel*npars_tel+2+telstart,Form("Azi%d",iTel),TelAzi_int[iTel],fabs(0.1/180.*PI),TelAzi_int[iTel]-telazi_margin,TelAzi_int[iTel]+telazi_margin);

      if(jdebug>1) printf("TelGeoFit::DoFit: iTel=%d Zen(%d)=%.2lf Azi(%d)=%.2lf\n",iTel,itel*npars_tel+1+telstart,TelZen_int[iTel]/PI*180,itel*npars_tel+2+telstart,TelAzi_int[iTel]/PI*180);
   }

   minimizer->Minimize();
   minimizer->Hesse();

   for(int irot=0;irot<nrot;irot++){
      int iRot=(int)(minimizer->X()[irot*npars_rot+rotstart+0]+0.5);
      for(int ii=0;ii<3;ii++) RotPos_fit[iRot-1][ii]=minimizer->X()[rotstart+irot*npars_rot+ii+1];
      RotEle_fit[iRot-1]=minimizer->X()[rotstart+irot*npars_rot+4];
      RotAzi_fit[iRot-1]=minimizer->X()[rotstart+irot*npars_rot+5];
      RefEle_fit[iRot-1]=minimizer->X()[rotstart+irot*npars_rot+6];
      RefAzi_fit[iRot-1]=minimizer->X()[rotstart+irot*npars_rot+7];
      RotTime_fit[iRot-1]=minimizer->X()[rotstart+irot*npars_rot+8];
   }
   for(int itel=0;itel<ntel;itel++){
      int iTel=(int)(minimizer->X()[itel*3+0+telstart]+0.5);
      TelZen_fit[iTel-1]=minimizer->X()[itel*3+1+telstart];
      TelAzi_fit[iTel-1]=minimizer->X()[itel*3+2+telstart];
   }
   return true;
}
bool TelGeoFit::FitProcedure1(int ntel,int* tellist,int nrot,int* rotlist){
return false;
   //bool fitted=false;
   //int tel_buff[10];

   //SetRotTelPars(true,false);
   //for(int itel=0;itel<ntel;itel++){
   //   tel_buff[0]=tellist[itel];
   //   if(itel==0){  ///fit the first tel
   //      FixRotPos=0;
   //      FixRotDir=0;
   //      fitted=DoFit(1,tel_buff);
   //   }
   //   else{  //fix Rot Pars, and fit other tels
   //      FixRotPos=0x7;
   //      FixRotDir=0xF;
   //      if(fitted&&minimizer){
   //      double pars_rot[7]={RotPos_fit[0],RotPos_fit[1],RotPos_fit[2],RotEle_fit,RotAzi_fit,RefEle_fit,RefAzi_fit};
   //      SetRotPars(pars_rot,true,false);
   //      fitted=DoFit(1,tel_buff);
   //      }
   //   }
   //}
   ////loose all the pars, and fit all tels
   //if(fitted){
   //FixRotPos=0;
   //FixRotDir=0;
   //for(int itel=0;itel<NCTMax;itel++){
   //   double pars_tel[2]={TelZen_fit[itel],TelAzi_fit[itel]};
   //   SetTelPars(itel+1,pars_tel,true,false);
   //}
   //fitted=DoFit(ntel,tellist);
   //}

   //return fitted;
}
bool TelGeoFit::FitProcedure2(int ntel,int* tellist,int nrot,int* rotlist){
return false;
   //bool fitted=false;
   //int tel_buff[10];

   //SetRotTelPars(true,false);
   //FixRotPos=0;
   //FixRotDir=0;
   //FixTelDir=0;
   //fitted=DoFit(ntel,tellist);

   //if(fitted&&minimizer){
   //   double pars_rot[7]={RotPos_fit[0],RotPos_fit[1],RotPos_fit[2],RotEle_fit,RotAzi_fit,RefEle_fit,RefAzi_fit};
   //   SetRotPars(pars_rot,true,false);
   //   FixRotPos=0x7;
   //   FixRotDir=0xF;
   //   for(int itel=0;itel<ntel;itel++){
   //      tel_buff[0]=tellist[itel];
   //      fitted=DoFit(1,tel_buff);
   //      if(fitted&&minimizer){
   //         double pars_tel[2]={TelZen_fit[tellist[itel]-1],TelAzi_fit[tellist[itel]-1]};
   //         SetTelPars(tellist[itel],pars_tel,true,false);
   //      }
   //   }
   //}

   //return fitted;
}
void TelGeoFit::DumpFit(){
   if(!minimizer){
      printf("TelGeoFit::DumpFit: Fit Failed\n");
   }
   else{
      printf("TelGeoFit::DumpFit: Fit Finished. Fitting pars:\n");
      int nrot=(int)(minimizer->X()[0]+0.5);
      int rotstart=1;
      const int npars_rot=9;
      int telstart=nrot*npars_rot+1+rotstart;
      const int npars_tel=3;
      printf("RotPars: nrot=%d\n",nrot);
      for(int irot=0;irot<nrot;irot++){
         int Li=(int)(minimizer->X()[rotstart+irot*npars_rot+0]+0.5);
         printf("Li=%d RotPosX=%7.2lf+-%6.3lf RotPosY=%7.2lf+-%6.3lf RotPosZ=%.2lf+-%.3lf\n",Li,minimizer->X()[rotstart+irot*npars_rot+1]/100,minimizer->Errors()[rotstart+irot*npars_rot+1]/100,minimizer->X()[rotstart+irot*npars_rot+2]/100,minimizer->Errors()[rotstart+irot*npars_rot+2]/100,minimizer->X()[rotstart+irot*npars_rot+3]/100,minimizer->Errors()[rotstart+irot*npars_rot+3]/100);
         printf("Li=%d RotEle=%6.2lf+-%6.3lf RotAzi=%6.2lf+-%6.3lf\n",Li,minimizer->X()[rotstart+irot*npars_rot+4]/PI*180,minimizer->Errors()[rotstart+irot*npars_rot+4]/PI*180,minimizer->X()[rotstart+irot*npars_rot+5]/PI*180,minimizer->Errors()[rotstart+irot*npars_rot+5]/PI*180);
         printf("Li=%d Ele_bias=%5.2lf+-%5.3lf Azi_bias=%5.2lf+-%5.3lf\n",Li,minimizer->X()[rotstart+irot*npars_rot+6]/PI*180,minimizer->Errors()[rotstart+irot*npars_rot+6]/PI*180,minimizer->X()[rotstart+irot*npars_rot+7]/PI*180,minimizer->Errors()[rotstart+irot*npars_rot+7]/PI*180);
         printf("Li=%d RotTime=%9.0lf+-%5.0lf\n",Li,minimizer->X()[rotstart+irot*npars_rot+8],minimizer->Errors()[rotstart+irot*npars_rot+8]);
      }

      int ntel=(int)(minimizer->X()[telstart-1]+0.5);
      printf("\n");
      printf("TelPars:\n");
      for(int itel=0;itel<ntel;itel++){
         int itel_fit=(int)(minimizer->X()[itel*npars_tel+0+telstart]+0.5);
         double izen_fit=minimizer->X()[itel*npars_tel+1+telstart];
         double err_zen_fit=minimizer->Errors()[itel*npars_tel+1+telstart];
         double iazi_fit=minimizer->X()[itel*npars_tel+2+telstart];
         double err_azi_fit=minimizer->Errors()[itel*npars_tel+2+telstart];
         printf("iTel=%2d: zen={%6.2lf+-%6.3lf} azi={%.2lf+-%.3lf}\n",itel_fit,izen_fit/PI*180,err_zen_fit/PI*180,iazi_fit/PI*180,err_azi_fit/PI*180);
      }
      printf("\n");
      printf("EventPars: Chi2/Ndof=%.1lf/%d\n",Chi2_All,Ndof_All);
      for(int ievt=0;ievt<nevent;ievt++){
         if(replaced[ievt]) continue;
         int iTel=telindex[ievt];
         int itel=iTel-1;
         if(iTel<0||itel>=NCTMax) continue;
         int irot=GetRotIndex(ievt);
         if(irot<0||irot>=NRotMax) continue;
         int iRot=RotateDB::rotindex[irot];

         double CC,phi;
         bool isfit_tel=TelZen_fit[itel]>=0;
         double teldir_fit[2]={isfit_tel?TelZen_fit[itel]:TelZen_int[itel],isfit_tel?TelAzi_fit[itel]:TelAzi_int[itel]};
         double rotdir_fit[2];
         double ele_in=PI/2-rotdir0[ievt][0]; //from zenith to elevation
         double azi_in=rotdir0[ievt][1];
         bool isfit_rot=TelZen_fit[itel]>=0;
         double rotpar_fit[4]={isfit_rot?RotEle_fit[iRot-1]:RotEle_int[iRot-1],isfit_rot?RotAzi_fit[iRot-1]:RotAzi_int[iRot-1],isfit_rot?RefEle_fit[iRot-1]:RefEle_int[iRot-1],isfit_rot?RefAzi_fit[iRot-1]:RefAzi_int[iRot-1]};
         CalDir_out(rotpar_fit[0],rotpar_fit[1],rotpar_fit[2],rotpar_fit[3],ele_in,azi_in,rotdir_fit[0],rotdir_fit[1]);
         rotdir_fit[0]=PI/2-rotdir_fit[0]; //from elevation to zenith
         GetCCphi(teldir_fit[0],teldir_fit[1],telpos[itel],RotPos_fit[iRot-1],rotdir_fit,CC,phi);
         printf("ievt=%5d iTel=%2d time=%.7lf iEvent=%5d Ele=%.0lf CC=%+6.2lf+-%6.2lf phi=%6.2lf+-%6.3lf CC_fit=%+6.2lf phi_fit=%6.2lf chi2/ndof=%.1lf/%d\n",ievt,telindex[ievt],rabbitTime[ievt],iEvent[ievt],90-rotdir0[ievt][0]/PI*180,image_pars[ievt][0][0]/PI*180,image_pars[ievt][0][1]/PI*180,image_pars[ievt][1][0]/PI*180,image_pars[ievt][1][1]/PI*180,CC/PI*180,phi/PI*180,XYchi2_evt[ievt],XYndof_evt[ievt]);
      }
      printf("\n");
   }

}
TCanvas* TelGeoFit::DrawFitXY(int iTel,int Li,int iEle,int Time,double plotrange){
   TCanvas* cc=new TCanvas();
   cc->SetTitle(Form("Chi2/Ndof=%.1lf/%d\n",Chi2_All,Ndof_All));
   double elelist[7]={10,20,30,40,50,55,60};
   int nplot=0;
   for(int ievt=0;ievt<nevent;ievt++){
      if(replaced[ievt]) continue;
      if(iTel>=1&&telindex[ievt]!=iTel) continue;
      if((iEle>=0&&iEle<7)&&fabs((90-rotdir0[ievt][0]/PI*180)-elelist[iEle])>0.5) continue;
      if(Time>=1300000000&&Time!=((int)rabbitTime[ievt])) continue;
      int irot=RotateDB::GetHead()->GetLi((rabbitTime[ievt]-(int)(rabbitTime[ievt]))/2.0e-8);
      if(irot<0) continue;
      if(Li>0&&RotateDB::rotindex[irot]!=Li) continue;
      nplot++;
   }
   int ndiv=sqrt(nplot);
   if(sqrt(nevent)-ndiv>0) ndiv++;
   cc->Divide(ndiv,ndiv);
   nplot=0;
   for(int ievt=0;ievt<nevent;ievt++){
      if(replaced[ievt]) continue;
      if(iTel>=1&&telindex[ievt]!=iTel) continue;
      if((iEle>=0&&iEle<7)&&fabs((90-rotdir0[ievt][0]/PI*180)-elelist[iEle])>0.5) continue;
      if(Time>=1300000000&&Time!=((int)rabbitTime[ievt])) continue;
      int irot=RotateDB::GetHead()->GetLi((rabbitTime[ievt]-(int)(rabbitTime[ievt]))/2.0e-8);
      if(irot<0) continue;
      int iRot=RotateDB::rotindex[irot];
      if(Li>0&&iRot!=Li) continue;

      double CC,phi;
      int itel=telindex[ievt]-1;
      bool isfit=TelZen_fit[itel]>=0;
      double teldir_fit[2]={isfit?TelZen_fit[itel]:TelZen_int[itel],isfit?TelAzi_fit[itel]:TelAzi_int[itel]};
      double rotdir_fit[2];
      double ele_in=PI/2-rotdir0[ievt][0]; //from zenith to elevation
      double azi_in=rotdir0[ievt][1];
      bool isfit_rot=TelZen_fit[itel]>=0;
      double rotpar_fit[4]={isfit_rot?RotEle_fit[iRot-1]:RotEle_int[iRot-1],isfit_rot?RotAzi_fit[iRot-1]:RotAzi_int[iRot-1],isfit_rot?RefEle_fit[iRot-1]:RefEle_int[iRot-1],isfit_rot?RefAzi_fit[iRot-1]:RefAzi_int[iRot-1]};
      CalDir_out(rotpar_fit[0],rotpar_fit[1],rotpar_fit[2],rotpar_fit[3],ele_in,azi_in,rotdir_fit[0],rotdir_fit[1]);
      rotdir_fit[0]=PI/2-rotdir_fit[0]; //from elevation to zenith
      GetCCphi(teldir_fit[0],teldir_fit[1],telpos[itel],RotPos_fit[iRot-1],rotdir_fit,CC,phi);
      if(jdebug>1) printf("TelGeoFit:DrawFitXY ievt=%d iTel=%d telzen=%.2lf telazi=%.2lf rotzen=%.2lf rotazi=%.2lf CC=%.2lf phi=%.2lf\n",ievt,itel+1,teldir_fit[0]/PI*180,teldir_fit[1]/PI*180,rotdir_fit[0]/PI*180,rotdir_fit[1]/PI*180,CC/PI*180,phi/PI*180);
      if(fabs(CC)>7./180*PI) continue;

      cc->cd(nplot+1);
      if(gPad){
         gPad->DrawFrame(-plotrange,-plotrange,plotrange,plotrange,Form("iTel=%d iEvent=%d Li=%d Ele=%.0lf time=%ld+%lf(%d-%02d-%02d %02d:%02d:%02d) Npe=%.0lf;X [degree];Y [degree]",telindex[ievt],iEvent[ievt],irot<0?0:RotateDB::rotindex[irot],(90-rotdir0[ievt][0]/PI*180),(long int)(rabbitTime[ievt]),rabbitTime[ievt]-(int)(rabbitTime[ievt]),CommonTools::TimeFlag((int)rabbitTime[ievt],1),CommonTools::TimeFlag((int)rabbitTime[ievt],2),CommonTools::TimeFlag((int)rabbitTime[ievt],3),CommonTools::TimeFlag((int)rabbitTime[ievt],4),CommonTools::TimeFlag((int)rabbitTime[ievt],5),CommonTools::TimeFlag((int)rabbitTime[ievt],6),Npe_sum[ievt]));
      }
      TH2Poly* hpoly=new TH2Poly();
      for(int ii=0;ii<NSIPM;ii++){
         double ImageX,ImageY;
         TelGeoFit::GetImageXYCoo(ii,ImageX,ImageY,-1,true);
         int ibin=hpoly->AddBin(ImageX-0.25,ImageY-0.25,ImageX+0.25,ImageY+0.25);
         if(Npe_sipm[ievt][ii]>0) hpoly->SetBinContent(ibin,Npe_sipm[ievt][ii]);
      }
      hpoly->Draw("colz same");
      hpoly->SetTitle(Form("iTel=%d iEvent=%d Li=%d Ele=%.0lf time=%ld+%lf(%d-%02d-%02d %02d:%02d:%02d) Npe=%.0lf;X [degree];Y [degree]",telindex[ievt],iEvent[ievt],irot<0?irot:RotateDB::rotindex[irot],(90-rotdir0[ievt][0]/PI*180),(long int)(rabbitTime[ievt]),rabbitTime[ievt]-(int)(rabbitTime[ievt]),CommonTools::TimeFlag((int)rabbitTime[ievt],1),CommonTools::TimeFlag((int)rabbitTime[ievt],2),CommonTools::TimeFlag((int)rabbitTime[ievt],3),CommonTools::TimeFlag((int)rabbitTime[ievt],4),CommonTools::TimeFlag((int)rabbitTime[ievt],5),CommonTools::TimeFlag((int)rabbitTime[ievt],6),Npe_sum[ievt]));

      TGraph* gr=new TGraph();
      for(int ii=0;ii<100;ii++){
         bool usex=fabs(phi-PI/2)<PI/4?false:true;
         double coo1=-10.+20/100.*(ii+0.5);
         double coo2=usex?(coo1*sin(phi)+CC/PI*180)/cos(phi):(coo1*cos(phi)-CC/PI*180)/sin(phi);
         gr->SetPoint(ii,usex?coo1:coo2,usex?coo2:coo1);
      }
      gr->SetLineColor(isfit?2:1);
      gr->SetLineWidth(FitLineWidth);
      gr->Draw("l");
      if(gPad) gPad->SetTitle(Form("chi2/ndof=%.1lf/%d",XYchi2_evt[ievt],XYndof_evt[ievt]));
      nplot++;
   }
   return cc;
}
TCanvas* TelGeoFit::DrawFitTime(int iTel,int Li,int iEle,int Time){
   TCanvas* cc=new TCanvas();
   cc->SetTitle(Form("Arrival Time vs Long Axis"));
   double elelist[7]={10,20,30,40,50,55,60};
   int nplot=0;
   for(int ievt=0;ievt<nevent;ievt++){
      if(replaced[ievt]) continue;
      if(iTel>=1&&telindex[ievt]!=iTel) continue;
      if((iEle>=0&&iEle<7)&&fabs((90-rotdir0[ievt][0]/PI*180)-elelist[iEle])>0.5) continue;
      if(Time>=1300000000&&Time!=((int)rabbitTime[ievt])) continue;
      int irot=RotateDB::GetHead()->GetLi((rabbitTime[ievt]-(int)(rabbitTime[ievt]))/2.0e-8);
      if(irot<0) continue;
      if(Li>0&&RotateDB::rotindex[irot]!=Li) continue;
      nplot++;
   }
   int ndiv=sqrt(nplot);
   if(sqrt(nevent)-ndiv>0) ndiv++;
   cc->Divide(ndiv,ndiv);
   static int ndraw=0;
   nplot=0;
   for(int ievt=0;ievt<nevent;ievt++){
      if(replaced[ievt]) continue;
      if(iTel>=1&&telindex[ievt]!=iTel) continue;
      int iele=-1;
      for(int ii=0;ii<7;ii++){
         if(fabs((90-rotdir0[ievt][0]/PI*180)-elelist[ii])<2) {iele=ii; break;}
      }
      if((iEle>=0&&iEle<7)&&iEle!=iele) continue;
      if(Time>=1300000000&&Time!=((int)rabbitTime[ievt])) continue;
      int irot=RotateDB::GetHead()->GetLi((rabbitTime[ievt]-(int)(rabbitTime[ievt]))/2.0e-8);
      if(irot<0) continue;
      int iRot=RotateDB::rotindex[irot];
      if(Li>0&&iRot!=Li) continue;
      cc->cd(nplot+1);
      int itel=telindex[ievt]-1;
      double rabbittime=(rabbitTime[ievt]-((int)rabbitTime[ievt]))/2.0e-8;
      if(itel<0||irot<0) continue;
      bool isfit=TelZen_fit[itel]>=0;
      double teldir_fit[2]={isfit?TelZen_fit[itel]:TelZen_int[itel],isfit?TelAzi_fit[itel]:TelAzi_int[itel]};
      double rotdir_fit[2];
      double ele_in=PI/2-rotdir0[ievt][0]; //from zenith to elevation
      double azi_in=rotdir0[ievt][1];
      bool isfit_rot=TelZen_fit[itel]>=0;
      double rotpar_fit[4]={isfit_rot?RotEle_fit[iRot-1]:RotEle_int[iRot-1],isfit_rot?RotAzi_fit[iRot-1]:RotAzi_int[iRot-1],isfit_rot?RefEle_fit[iRot-1]:RefEle_int[iRot-1],isfit_rot?RefAzi_fit[iRot-1]:RefAzi_int[iRot-1]};
      CalDir_out(rotpar_fit[0],rotpar_fit[1],rotpar_fit[2],rotpar_fit[3],ele_in,azi_in,rotdir_fit[0],rotdir_fit[1]);
      rotdir_fit[0]=PI/2-rotdir_fit[0]; //from elevation to zenith
      double CC,phi;
      GetCCphi(teldir_fit[0],teldir_fit[1],telpos[itel],RotPos_fit[iRot-1],rotdir_fit,CC,phi);

      TH1F* ncount=new TH1F(Form("count_%d",ievt),";Long Axis;Count",26,-13.,13.);
      TH1F* ncount_weight=new TH1F(Form("count_%d_weight",ievt),";Long Axis;Weighted Count",26,-13.,13.);
      TH1F* time_long_data=new TH1F(Form("Time_Long_Data_Evt%d_Tel%d_Ele%d_Draw%d",ievt,telindex[ievt],iele,ndraw),";Long Axis [degree];Peak Time [ns]",26,-13.,13.);
      double time_ref=0;
      //double maxpe=0;
      //for(int ii=0;ii<MAXPMT;ii++){
      //   if(Npe_sipm[ievt][ii]>maxpe){
      //      maxpe=Npe_sipm[ievt][ii];
      //      time_ref=Time_sipm[ievt][ii];
      //   }
      //}
      for(int ii=0;ii<MAXPMT;ii++){
         if(Npe_sipm[ievt][ii]<=0) continue;
         double ImageX,ImageY;
         TelGeoFit::GetImageXYCoo(ii,ImageX,ImageY,-1,true);
         double longaxis=ImageX*cos(phi)+ImageY*sin(phi);
         double probi=Npe_sipm[ievt][ii]/Npe_sum[ievt];
         ncount->Fill(longaxis);
         ncount_weight->Fill(longaxis,probi);
         double datatime=(rabbittime*20-RotTime_fit[iRot-1])+(Time_sipm[ievt][ii]-time_ref)*80;
         time_long_data->Fill(longaxis,datatime*probi);
      }
      ////remove first and last non-zero bins
      //for(int ii=1;ii<=ncount->GetNbinsX();ii++){
      //   if(ncount->GetBinContent(ii)>0) {time_long_data->SetBinContent(ii,0.); break;}
      //}
      //for(int ii=ncount->GetNbinsX();ii>=1;ii--){
      //   if(ncount->GetBinContent(ii)>0) {time_long_data->SetBinContent(ii,0.); break;}
      //}
      double minimumsipm=1;
      for(int ii=1;ii<=ncount->GetNbinsX();ii++){
         if(ncount->GetBinContent(ii)<minimumsipm||ncount_weight->GetBinContent(ii)==0) time_long_data->SetBinContent(ii,0);
      }

      time_long_data->Divide(ncount_weight);
      double tmaximum=-1.0e9,tminimum=1.0e9;
      for(int ii=1;ii<=time_long_data->GetNbinsX();ii++){
         double content=time_long_data->GetBinContent(ii);
         if(content!=0){
            time_long_data->SetBinError(ii,80.);
            if(content<tminimum) tminimum=content;
            if(content>tmaximum) tmaximum=content;
         }
      }
      if(jdebug>9){
         for(int ii=1;ii<=time_long_data->GetNbinsX();ii++){
            printf("ievt=%2d bin=%2d ncount=%.0lf weight=%.2le time_weight=%.2le\n",ievt,ii,ncount->GetBinContent(ii),ncount_weight->GetBinContent(ii),time_long_data->GetBinContent(ii));
         }
      }
      delete ncount;
      delete ncount_weight;
      //delete time_long_data;

      time_long_data->Draw("");
      time_long_data->SetTitle(Form("iTel=%d iEvent=%d Li=%d Ele=%.0lf time=%ld+%lf(%d-%02d-%02d %02d:%02d:%02d) Time=%.0lf",telindex[ievt],iEvent[ievt],irot<0?irot:RotateDB::rotindex[irot],(90-rotdir0[ievt][0]/PI*180),(long int)(rabbitTime[ievt]),rabbitTime[ievt]-(int)(rabbitTime[ievt]),CommonTools::TimeFlag((int)rabbitTime[ievt],1),CommonTools::TimeFlag((int)rabbitTime[ievt],2),CommonTools::TimeFlag((int)rabbitTime[ievt],3),CommonTools::TimeFlag((int)rabbitTime[ievt],4),CommonTools::TimeFlag((int)rabbitTime[ievt],5),CommonTools::TimeFlag((int)rabbitTime[ievt],6),time_long_data->GetMaximum()));
      time_long_data->GetYaxis()->SetRangeUser(tminimum-(iRot==2?800:1500),tmaximum+(iRot==2?800:1500));

      //model
      double zenith=TelZen_fit[itel];
      double azimuth=TelAzi_fit[itel];
      double telcoo[3]={telpos[itel][0],telpos[itel][1],telpos[itel][2]};
      double incoo[3]={RotPos_fit[iRot-1][0],RotPos_fit[iRot-1][1],RotPos_fit[iRot-1][2]};
      double indir[2]={rotdir_fit[0],rotdir_fit[1]};

      TGraph* gr=new TGraph();
      for(int ii=0;ii<100;ii++){
         double longcoo=-13.+26/100.*(ii+0.5);
         double dist=GetTimeFromLongcoo(zenith,azimuth,telcoo,incoo,indir,longcoo/180*PI);
         if(dist==TelGeoFit::InfNumber) continue;
         double sipmtime=dist/29.98;  //in nano-second
         gr->SetPoint(gr->GetN(),longcoo,sipmtime);
         if(jdebug>9){
            if(ii>=0) printf("ievt=%2d point=%d longcoo=%.2lf model_time=%.2le model_dist=%.2le\n",ievt,ii,longcoo,sipmtime,dist);
         }
      }
      gr->SetLineColor(FixRotTime==0?2:1);
      gr->SetLineWidth(FitLineWidth);
      gr->Draw("l");
      nplot++;
   }
   ndraw++;
   return cc;
}
TCanvas* TelGeoFit::DrawNpeScatA(int iTel,int Li,int iEle,int Time){
   TCanvas* cc=new TCanvas();
   cc->SetTitle(Form("Npe vs Scatter Angle"));
   double elelist[7]={10,20,30,40,50,55,60};
   double scatrange[2][7][2]={{{80,140},{95,150},{110,160},{120,170},{135,180},{140,180},{150,180}},
                              {{60,150},{70,160},{80,160},{100,170},{100,180},{110,180},{120,180}}
                             };
   int nplot=0;
   for(int ievt=0;ievt<nevent;ievt++){
      if(replaced[ievt]) continue;
      if(iTel>=1&&telindex[ievt]!=iTel) continue;
      if((iEle>=0&&iEle<7)&&fabs((90-rotdir0[ievt][0]/PI*180)-elelist[iEle])>0.5) continue;
      if(Time>=1300000000&&Time!=((int)rabbitTime[ievt])) continue;
      int irot=RotateDB::GetHead()->GetLi((rabbitTime[ievt]-(int)(rabbitTime[ievt]))/2.0e-8);
      if(irot<0) continue;
      if(Li>0&&RotateDB::rotindex[irot]!=Li) continue;
      nplot++;
   }
   int ndiv=sqrt(nplot);
   if(sqrt(nevent)-ndiv>0) ndiv++;
   cc->Divide(ndiv,ndiv);
   static int ndraw=0;
   nplot=0;
   for(int ievt=0;ievt<nevent;ievt++){
      if(replaced[ievt]) continue;
      if(iTel>=1&&telindex[ievt]!=iTel) continue;
      int iele=-1;
      for(int ii=0;ii<7;ii++){
         if(fabs((90-rotdir0[ievt][0]/PI*180)-elelist[ii])<2) {iele=ii; break;}
      }
      if((iEle>=0&&iEle<7)&&iEle!=iele) continue;
      if(Time>=1300000000&&Time!=((int)rabbitTime[ievt])) continue;
      int irot=RotateDB::GetHead()->GetLi((rabbitTime[ievt]-(int)(rabbitTime[ievt]))/2.0e-8);
      if(irot<0) continue;
      int iRot=RotateDB::rotindex[irot];
      if(Li>0&&iRot!=Li) continue;
      cc->cd(nplot+1);
      //if(gPad) gPad->SetLogy(1);
      int itel=telindex[ievt]-1;
      if(itel<0||irot<0) continue;
      //bool isfit=TelZen_fit[itel]>=0;
      //double teldir_fit[2]={isfit?TelZen_fit[itel]:TelZen_int[itel],isfit?TelAzi_fit[itel]:TelAzi_int[itel]};
      //double rotdir_fit[2];
      //double ele_in=PI/2-rotdir0[ievt][0]; //from zenith to elevation
      //double azi_in=rotdir0[ievt][1];
      //bool isfit_rot=TelZen_fit[itel]>=0;
      //double rotpar_fit[4]={isfit_rot?RotEle_fit[iRot-1]:RotEle_int[iRot-1],isfit_rot?RotAzi_fit[iRot-1]:RotAzi_int[iRot-1],isfit_rot?RefEle_fit[iRot-1]:RefEle_int[iRot-1],isfit_rot?RefAzi_fit[iRot-1]:RefAzi_int[iRot-1]};
      //CalDir_out(rotpar_fit[0],rotpar_fit[1],rotpar_fit[2],rotpar_fit[3],ele_in,azi_in,rotdir_fit[0],rotdir_fit[1]);
      //rotdir_fit[0]=PI/2-rotdir_fit[0]; //from elevation to zenith
      //double CC,phi;
      //GetCCphi(teldir_fit[0],teldir_fit[1],telpos[itel],RotPos_fit[iRot-1],rotdir_fit,CC,phi);
      //if(fabs(CC)>7./180*PI) continue;

      //TH1F* npe_scat_data=new TH1F(Form("Npe_ScatA_Data_Evt%d_Tel%d_Ele%d_Draw%d",ievt,telindex[ievt],iele,ndraw),";Scatter Angle [degree];Npe [pe]",100,0,180.);
      //for(int ii=0;ii<MAXPMT;ii++){
      //   if(Npe_sipm[ievt][ii]<=0) continue;
      //   double ImageX,ImageY;
      //   TelGeoFit::GetImageXYCoo(ii,ImageX,ImageY,-1,false);
      //   double longaxis=ImageX*cos(phi)+ImageY*sin(phi);
      //   double scat_angle=GetScatAngleFromLongcoo(teldir_fit[0],teldir_fit[1],telpos[itel],RotPos_fit[iRot-1],rotdir_fit,longaxis);
      //   if(scat_angle==InfNumber) continue;
      //   npe_scat_data->Fill(scat_angle/PI*180,Npe_sipm[ievt][ii]);
      //}

      TH1F* npe_scat_data=GetScatterAngle(ievt,true,false);

      if(jdebug>9){
         for(int ii=1;ii<=npe_scat_data->GetNbinsX();ii++){
            printf("ievt=%2d bin=%2d longcoo=%.2lf npe=%.2le\n",ievt,ii,npe_scat_data->GetXaxis()->GetBinCenter(ii),npe_scat_data->GetBinContent(ii));
         }
      }

      npe_scat_data->Draw("");
      npe_scat_data->SetTitle(Form("iTel=%d iEvent=%d Li=%d Ele=%.0lf time=%ld+%lf(%d-%02d-%02d %02d:%02d:%02d) Npe=%.0lf",telindex[ievt],iEvent[ievt],irot<0?irot:RotateDB::rotindex[irot],(90-rotdir0[ievt][0]/PI*180),(long int)(rabbitTime[ievt]),rabbitTime[ievt]-(int)(rabbitTime[ievt]),CommonTools::TimeFlag((int)rabbitTime[ievt],1),CommonTools::TimeFlag((int)rabbitTime[ievt],2),CommonTools::TimeFlag((int)rabbitTime[ievt],3),CommonTools::TimeFlag((int)rabbitTime[ievt],4),CommonTools::TimeFlag((int)rabbitTime[ievt],5),CommonTools::TimeFlag((int)rabbitTime[ievt],6),npe_scat_data->GetMaximum()));
      npe_scat_data->SetMarkerColor(4);
      npe_scat_data->SetMarkerSize(1.2);
      npe_scat_data->SetLineColor(4);
      npe_scat_data->SetLineWidth(2);
      npe_scat_data->GetXaxis()->SetRangeUser(scatrange[irot][iele][0],scatrange[irot][iele][1]);
      //npe_scat_data->GetXaxis()->SetRangeUser(60,180);
      //npe_scat_data->GetYaxis()->SetRangeUser(tminimum-2000,tmaximum+2000);

      nplot++;
   }
   ndraw++;
   return cc;
}

