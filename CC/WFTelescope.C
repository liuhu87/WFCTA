#include "WFTelescope.h"
#include "WFCTAEvent.h"
#include "WFMirror.h"
#include "WFCone.h"
#include "WFCamera.h"
#include "Readparam.h"
#include "CorsikaIO.h"
#include "TMath.h"
#include "TRandom.h"
#include "TFile.h"
#include "TGraph.h"
#include "TelGeoFit.h"
#include <vector>
using std::vector;

double WFTelescopeArray::lhaaso_coo[3]={0,0,0};
int WFTelescopeArray::jdebug=0;
bool WFTelescopeArray::DoSim=true;
int WFTelescopeArray::CTNumber=1;
WFTelescopeArray* WFTelescopeArray::_Head=0;
TRandom3* WFTelescopeArray::prandom=0;
double WFTelescopeArray::viewunit=0.5;
vector<int> WFTelescopeArray::TelID;
vector<TGraph*> WFTelescopeArray::gTelView;
bool WFTelescopeArray::LoadViewExt=true;
WFTelescopeArray* WFTelescopeArray::GetHead(char* filename,int seed){
   if(!DoSim) return 0;
   if(_Head) return _Head;
   else{
      if(filename) _Head=new WFTelescopeArray(filename,seed);
      return _Head;
   }
}
void WFTelescopeArray::Init(bool DoNew,int seed){
   if(!prandom) prandom=new TRandom3();
   prandom->SetSeed(seed);
   if(!DoNew) pct=0;
   else{
      pct=new WFTelescope*[CTNumber];
      for(int ii=0;ii<CTNumber;ii++){
         pct[ii]=new WFTelescope();
      }
   }
}
void WFTelescopeArray::Clear(){
   for(int ii=0;ii<CTNumber;ii++){
      if(pct[ii]) {delete pct[ii]; pct[ii]=0;}
   }
   delete []pct;
}
void WFTelescopeArray::ReadFromFile(char* filename){
  WReadConfig read;
  WReadConfig* readconfig=&read;
  readconfig->readparam(filename);

  int WLFlag, ThinFlag, FilterFlag,MirrorPointErrorFlag, MirrorGeometry;
  float MirrorSizeX, MirrorSizeY, MirrorSpot, MirrorPointError;
  char inputfilename[PATH_MAX_LENGTH],outputfilename[PATH_MAX_LENGTH];


  int Fadc_bins, Fadc_length; //The number of FADC bins and the bin length 
  float nsb;                  //night sky background 

  /*The telescope arrays settings*/
  int *CT_Index;
  float *CT_X, *CT_Y, *CT_Z, *CT_Azi, *CT_Zen;
  float *timefirst, *timelast,*timemean, *ncphoton;
  /*cos directions of the pointing of telescopes */
  float *CT_m, *CT_n, *CT_l;
  /*cos directions of the primary directions of showers */
  float prim_m, prim_n, prim_l;

  /*space angle between the arriving directions of showers and the pointing of the telescopes */
  float *Space_angle;

  float zenith, azimuth, Xmax, Nmax;
  double  mm, nn;

  WLFlag = readconfig->GetWLFlag();
  ThinFlag = readconfig->GetThinFlag();
  CorsikaIO::ThinFlag=ThinFlag;
  CorsikaIO::WLFlag=WLFlag;
  FilterFlag = readconfig->GetFilterFlag();

  MirrorGeometry = readconfig->GetMirrorGeometry();
  MirrorPointErrorFlag = readconfig->GetMirrorPointErrorFlag();
  MirrorPointError = readconfig->GetMirrorPointError();

  readconfig->GetMirrorSize(&MirrorSizeX , &MirrorSizeY);
  WFTelescope::MirrorSizeX=MirrorSizeX;
  WFTelescope::MirrorSizeY=MirrorSizeY;
  MirrorSpot =readconfig->GetMirrorSpot();

  //WFMirror::SetMirrorPointErrorFlag(MirrorPointErrorFlag);
  //WFMirror::SetMirrorPointError(MirrorPointError);
  WFMirror::SetMirrorSpot(MirrorSpot);
  WFMirror::SetMirrorGeometry(MirrorGeometry);

  lhaaso_coo[0]=readconfig->GetLHAASOCoo(0);
  lhaaso_coo[1]=readconfig->GetLHAASOCoo(1);
  lhaaso_coo[2]=readconfig->GetLHAASOCoo(2);

  CTNumber =readconfig->GetCTNumber();

  CT_Index = new int[CTNumber];
  CT_X = new float[CTNumber];
  CT_Y = new float[CTNumber];
  CT_Z = new float[CTNumber];
  CT_Azi = new float[CTNumber];
  CT_Zen = new float[CTNumber];
  CT_m = new float[CTNumber];
  CT_n = new float[CTNumber];
  CT_l = new float[CTNumber];
  Space_angle = new float[CTNumber];
  timefirst = new float[CTNumber];
  timelast = new float[CTNumber];
  timemean = new float[CTNumber];
  ncphoton = new float[CTNumber];

  int nict=0;
  for (int ict = 1; ict <=NCTMax; ict++){
     if(nict>=CTNumber) break;
     if(readconfig->GetCTPosition(ict,3)<0) continue;
     CT_Index[nict] = ict;
     CT_X[nict] =  readconfig->GetCTPosition(ict,0);
     CT_Y[nict] =  readconfig->GetCTPosition(ict,1);
     CT_Z[nict] =  readconfig->GetCTPosition(ict,2);
     CT_Zen[nict] =  readconfig->GetCTPosition(ict,3) * TMath::DegToRad();
     CT_Azi[nict] =  readconfig->GetCTPosition(ict,4) * TMath::DegToRad();
     CT_m[nict] = sin(CT_Zen[nict]) * cos(CT_Azi[nict]);
     CT_n[nict] = sin(CT_Zen[nict]) * sin(CT_Azi[nict]);
     CT_l[nict] = cos(CT_Zen[nict]);
     //printf("Read: itel=%d iTel=%d pos={%f,%f,%f} dir={%f,%f}\n",nict,ict,readconfig->GetCTPosition(ict,0),readconfig->GetCTPosition(ict,1),readconfig->GetCTPosition(ict,2),readconfig->GetCTPosition(ict,3),readconfig->GetCTPosition(ict,4));
     //printf("Read: itel=%d iTel=%d pos={%f,%f,%f} dir={%f,%f}\n",nict,CT_Index[nict],CT_X[nict],CT_Y[nict],CT_Z[nict],CT_Zen[nict]/PI*180,CT_Azi[nict]/PI*180);
     nict++;
  }
  //*The total intensty of nsb in the trigger window equlas Fadc_bins X Fabs_length X nsb *//
  Fadc_bins = readconfig->GetFadcBins();
  Fadc_length = readconfig->GetFadcLength();
  nsb = readconfig->GetNSB();
  WCamera::SetSiPMMAP();
  WCamera::SetNSB(nsb*Fadc_length*Fadc_bins);

  //* Init the telescope array *//
  if(!pct) pct=new WFTelescope*[CTNumber];
  for(int ict=0; ict<CTNumber; ict++){
     pct[ict]=new WFTelescope();
     pct[ict]->TelIndex_=CT_Index[ict];
     pct[ict]->SetXYZ(CT_X[ict],CT_Y[ict],CT_Z[ict]);
     pct[ict]->SetPointing(CT_Zen[ict],CT_Azi[ict]);
     pct[ict]->SetEulerMatrix(CT_Zen[ict],CT_Azi[ict]);
     //GetMirror(ict)->SetMirror();
     GetMirror(ict)->SetMirrorPointError(MirrorPointErrorFlag,MirrorPointError);
     //GetCamera(ict)->SetCTNumber(CTNumber);
     //GetCamera(ict)->Init();
     //printf("Adding new telescope: itel=%d iTel=%d pos={%lf,%lf,%lf} dir={%lf,%lf}\n",ict,CT_Index[ict],CT_X[ict],CT_Y[ict],CT_Z[ict],CT_Zen[ict]/PI*180,CT_Azi[ict]/PI*180);
  }

  delete []CT_Index;
  delete []CT_X;
  delete []CT_Y;
  delete []CT_Z;
  delete []CT_Azi;
  delete []CT_Zen;
  delete []CT_m;
  delete []CT_n;
  delete []CT_l;
  delete []Space_angle;
  delete []timefirst;
  delete []timelast;
  delete []timemean;
  delete []ncphoton;
}
int WFTelescopeArray::WhichTel(double x0, double y0,double z0,double m1,double n1,double l1){
   int whichct=-1;
   if(!CheckTelescope()) return whichct;
   if(m1*m1+n1*n1+l1*l1<=0) return whichct;
   double mindist=10000000000000;
   for(int ict=0;ict<CTNumber;ict++){
      double x1=(l1==0)?x0:(x0+(pct[ict]->Telz_+lhaaso_coo[2]-z0)/l1*m1);
      double y1=(l1==0)?y0:(y0+(pct[ict]->Telz_+lhaaso_coo[2]-z0)/l1*n1);
      double dist=sqrt(pow(x1-(pct[ict]->Telx_),2)+pow(y1-(pct[ict]->Tely_),2));
      if(dist<mindist){
         mindist=dist;
         whichct=ict;
      }
   }
   return whichct;
}
int WFTelescopeArray::RayTrace(double x0, double y0, double z0, double m1, double n1, double l1,double weight,double wavelength,int &itel,double &t,int &itube,int &icell){
   int whichct=WhichTel(x0,y0,z0,m1,n1,l1);
   if(whichct<0) return -1000;
   ///inside telescope
   if(jdebug>0) printf("WFTelescopeArray::RayTrace: Passing Telescope InCoo={%f,%f,%f} TelCoo={%f,%f,%f}\n",x0,y0,z0,pct[whichct]->Telx_,pct[whichct]->Tely_,pct[whichct]->Telz_+lhaaso_coo[2]);
   double xtel=(l1==0)?x0:(x0+(pct[whichct]->Telz_+lhaaso_coo[2]-z0)/l1*m1);
   double ytel=(l1==0)?y0:(y0+(pct[whichct]->Telz_+lhaaso_coo[2]-z0)/l1*n1);
   if(!pct[whichct]->IncidentTel(xtel,ytel)) return -1;
   if(jdebug>0) printf("WFTelescopeArray::RayTrace: Inside Telescope%d\n",whichct);

   ///apply all the efficiencies
   if(jdebug>0) printf("WFTelescopeArray::RayTrace: Before Efficiency\n");
   if(!pct[whichct]->Reflected(wavelength)) return -2;
   if(jdebug>0) printf("WFTelescopeArray::RayTrace: After Efficiency\n");

   ///convert coor.
   double x0new,y0new,z0new;
   double m1new,n1new,l1new;
   pct[whichct]->ConvertCoor(x0,y0,z0,m1,n1,l1,x0new,y0new,z0new,m1new,n1new,l1new);
   if(jdebug>0) printf("WFTelescopeArray::RayTrace(Convert Coor.): InCoor={%.3f,%.3f,%.3f} TelCoor={%.2f,%.2f} InDir={%.3f,%.3f,%.3f}\n",x0,y0,z0,pct[whichct]->Telx_,pct[whichct]->Tely_,m1,n1,l1);
   if(jdebug>0) printf("WFTelescopeArray::RayTrace(Convert Coor.): OtCoor={%.3f,%.3f,%.3f} TelCoor={%.2f,%.2f} OtDir={%.3f,%.3f,%.3f}\n",x0new,y0new,z0new,pct[whichct]->Telx_,pct[whichct]->Tely_,m1new,n1new,l1new);

   ///Ray Trace
   int result=pct[whichct]->RayTrace(x0new,y0new,z0new,m1new,n1new,l1new,wavelength,t,itube,icell);
   if(jdebug>0) printf("WFTelescopeArray::RayTrace: From Door to SiPM. result=%d\n",result);
   if(result>=0){
      GetCamera(whichct)->Fill(itube,t,weight);
      if(jdebug>0) printf("WFTelescopeArray::RayTrace: Filling the camera itube=%d weight=%le\n",itube,weight);
   }
   itel=whichct;
   return result>=0?result:(result-2);
}
bool WFTelescopeArray::CheckTelescope(){
   if(!pct) return false;
   bool exist=true;
   for(int ii=0;ii<WFTelescopeArray::CTNumber;ii++){
      if(!pct[ii]) {exist=false; break;}
   }
   return exist;
}
int WFTelescopeArray::GetTelescope(int iTel){
   int res=-1;
   for(int ii=0;ii<WFTelescopeArray::CTNumber;ii++){
      WFTelescope* pct0=pct[ii];
      if(!pct0) continue;
      if(iTel==pct0->TelIndex_) {res=ii; break;}
   }
   return res;
}
TGraph* WFTelescopeArray::TelViewCal(int iTel,bool IsLocal){
   WFTelescope* pct0=0;
   int itel=pct?GetTelescope(iTel):-1;
   if(itel>=0) pct0=pct[itel];
   bool exist=bool(pct0);
   if(!exist){
      if(IsLocal) pct0=new WFTelescope();
   }
   if(!pct0) return 0;
   int nbin=(int)(180/viewunit);
   vector<double> xlist;
   vector<double> ylist1;
   vector<double> ylist2;
   for(int ibinx=1;ibinx<=nbin;ibinx++){
      double xx=-90.+(90.-(-90.))/nbin*(ibinx-0.5);
      double ymin=100;
      double ymax=-100;
      bool inside_pre=false;
      for(int ibiny=1;ibiny<=nbin;ibiny++){
         double yy=-90.+(90.-(-90.))/nbin*(ibiny-0.5);
         double theta=sqrt(xx*xx+yy*yy)/180.*PI;
         if(theta>=PI/2) continue;
         bool inside=false;
         double x2=yy;
         double y2=-xx;
         double phi=0;
         if(x2!=0||y2!=0){
            phi=acos(x2/sqrt(x2*x2+y2*y2));
            if(y2<0) phi=2*PI-phi;
         }
         double dir_in[3]={-cos(phi)*sin(theta),-sin(phi)*sin(theta),-cos(theta)};
         double dir_tel[3]={0,0,-1};
         if(!IsLocal){
            dir_tel[0]=-cos(pct0->TelA_)*sin(pct0->TelZ_);
            dir_tel[1]=-sin(pct0->TelA_)*sin(pct0->TelZ_);
            dir_tel[2]=-cos(pct0->TelZ_);
         }
         if((dir_in[0]*dir_tel[0]+dir_in[1]*dir_tel[1]+dir_in[2]*dir_tel[2])<cos(13./180*PI)) continue;
         double imagexy[2];
         bool res=pct0->GetImageCoo(dir_in,imagexy[0],imagexy[1],IsLocal);
         if(res){
            int retval=TelGeoFit::IsInside(imagexy[0],imagexy[1]);
            if(retval>=-1) inside=true;
         }
         if(!inside_pre){
            if(inside&&ymin>99){
               ymin=yy;
            }
         }
         else{
            if(!inside) ymax=yy;
         }
         inside_pre=inside;
      }
      if(ymin<=ymax){
         xlist.push_back(xx);
         ylist1.push_back(ymin);
         ylist2.push_back(ymax);
      }
   }
   if((!exist)&&IsLocal) delete pct0;

   if(xlist.size()<=0) return 0;
   TGraph* gr=new TGraph();
   for(int ii=0;ii<xlist.size();ii++){
      gr->SetPoint(gr->GetN(),xlist.at(ii),ylist1.at(ii));
   }
   for(int ii=xlist.size()-1;ii>=0;ii--){
      gr->SetPoint(gr->GetN(),xlist.at(ii),ylist2.at(ii));
   }
   gr->SetPoint(gr->GetN(),xlist.at(0),ylist1.at(0));
   return gr;
}
TGraph* WFTelescopeArray::TelViewMC(int iTel,bool IsLocal){
   if(!pct) return 0;
   int itel=GetTelescope(iTel);
   if(itel<0) return 0;
   WFTelescope* pct0=pct[itel];
   if(!pct0) return 0;
   const int nbin=1000;
   double xbin[nbin+1];
   double ybin[nbin+1][2];
   for(int ibin=0;ibin<=nbin;ibin++){
      xbin[ibin]=-100.;
      double xx=-1.+2./nbin*ibin;
      bool inside_pre=false;
      for(int ibin2=0;ibin2<=nbin;ibin2++){
         double yy=-1.+2./nbin*ibin2;
         double rr=sqrt(xx*xx+yy*yy);
         if(rr>1.) continue;
         double zz=sqrt(1-rr*rr);
         double x0=0,y0=0,z0=0;
         double m1=-xx;
         double n1=-yy;
         double l1=-zz;

         //check wheather inside the field of view
         bool inside=false;
         bool quickpro=true;
         double x0new,y0new,z0new;
         double m1new,n1new,l1new;
         if(quickpro){
            double dir_in[3]={m1,n1,l1};
            double imagexy[2];
            bool res=pct0->GetImageCoo(dir_in,imagexy[0],imagexy[1],IsLocal);
            if(res){
               inside=TelGeoFit::IsInside(imagexy[0],imagexy[1])>=0;
            }
            if(jdebug>1) printf("WFTelescopeArray::TelView: res=%d inside=%d dir_in={%.2lf,%.2lf,%.2lf} imagexy={%.2lf,%.2lf}\n",res,inside,dir_in[0],dir_in[1],dir_in[2],imagexy[0],imagexy[1]);
         }
         else{
            if(IsLocal){
               x0new=x0*10;
               y0new=y0*10;
               z0new=z0*10;
               m1new=m1;
               n1new=n1;
               l1new=l1;
            }
            else{
               pct0->ConvertCoor(x0,y0,z0,m1,n1,l1,x0new,y0new,z0new,m1new,n1new,l1new);
            }
            z0new=WFTelescope::ZDOOR;
            double margin=100.;
            int ntest=60;
            for(int itest=0;itest<ntest*ntest;itest++){
               int ix=itest/ntest;
               int iy=itest%ntest;
               x0new=-(WFTelescope::D_DOOR-margin)/2.+(WFTelescope::D_DOOR-margin)/(ntest-1)*ix;
               y0new=-(WFTelescope::Hdoor-margin)/2.+(WFTelescope::Hdoor-margin)/(ntest-1)*iy;
               double t,xcluster,ycluster,m2,n2,l2;
               int result=pct0->RayTraceUpToCone(x0new,y0new,z0new,m1new,n1new,l1new,t,xcluster,ycluster,m2,n2,l2);
               if(result>=0){
                  double xc = ycluster;
                  double yc = xcluster;
                  int itube=(pct0->pcame)->GetCone(xc,yc);
                  if(itube>=0) inside=true;
                  if(inside) break;
               }
            }
         }
         if(jdebug>0&&inside) printf("WFTelescopeArray::TelView: xbin=%d(xx=%.2lf) ybin=%d(yy=%.2lf) zz=%lf inside=%d\n",ibin,xx,ibin2,yy,zz,inside);
         if(inside_pre==false){
            if(inside==true){
               xbin[ibin]=xx;
               ybin[ibin][0]=yy;
               inside_pre=inside;
            }
         }
         else{
            if(inside==false){
               xbin[ibin]=xx;
               ybin[ibin][1]=yy-2./nbin;
               break;
            }
         }
      }
   }

   TGraph* gr=new TGraph();
   int np=0;
   double x0,y0;
   for(int ibin=0;ibin<=nbin;ibin++){
      if(xbin[ibin]<-2) continue;
      if(jdebug>1) printf("WFTelescopeArray::TelView: np=%d xx=%.2lf yy=(%.2lf,%.2lf)\n",np,xbin[ibin],ybin[ibin][0],ybin[ibin][1]);
      if(np==0){x0=xbin[ibin]; y0=ybin[ibin][0];}
      gr->SetPoint(np,xbin[ibin],ybin[ibin][0]);
      np++;
   }
   for(int ibin=nbin;ibin>=0;ibin--){
      if(xbin[ibin]<-2) continue;
      gr->SetPoint(np,xbin[ibin],ybin[ibin][1]);
      np++;
   }
   gr->SetPoint(np,x0,y0);
   if(np<1) {delete gr; return 0;}
   else return gr;
}
int WFTelescopeArray::SaveTelView(char* filename){
   TDirectory* gdir=gDirectory;
   TFile* fout=TFile::Open(filename,"RECREATE");
   int ngr=0;
   for(int ii=0;ii<=CTNumber;ii++){
      TGraph* gr=ii>0?TelViewCal(ii,false):TelViewCal(1,true);
      if(!gr) continue;
      gr->Write(Form("Tel%d",ii));
      ngr++;
   }
   fout->Close();
   gdir->cd();
   return ngr;
}
int WFTelescopeArray::LoadTelView(char* filename){
   if(!filename) return 0;
   if(TelID.size()>0){
      TelID.clear();
      for(int ii=0;ii<gTelView.size();ii++){
         TGraph* gr=gTelView.at(ii);
         delete gr;
      }
      gTelView.clear();
   }
   TDirectory* gdir=gDirectory;
   TFile* fin=TFile::Open(filename,"READ");
   int ngr=0;
   for(int ii=0;ii<=CTNumber;ii++){
      int itel=GetTelescope(ii);
      if(ii>0&&itel<0) continue;
      TGraph* gr=(TGraph*)fin->Get(Form("Tel%d",ii));
      if(!gr) continue;
      TelID.push_back(ii);
      gTelView.push_back(gr);
      ngr++;
   }
   fin->Close();
   gdir->cd();
   return ngr;
}
TGraph* WFTelescopeArray::TelView(int iTel,bool IsLocal){
   if(LoadViewExt){
      for(int ii=0;ii<TelID.size();ii++){
         if(IsLocal){if(TelID.at(ii)==0) return gTelView.at(ii);}
         else{
            if(iTel>0&&iTel==TelID.at(ii)) return gTelView.at(ii);
         }
      }
   }
   //if not in the vector list, then just calculate
   TGraph* gr=TelViewCal(iTel,IsLocal);
   if(!gr) return 0;
   int iTel_rec=IsLocal?0:iTel;
   for(int ii=0;ii<TelID.size();ii++){
      if(iTel_rec==TelID.at(ii)) return gr;
   }
   //push this graph in the list
   TelID.push_back(iTel_rec);
   gTelView.push_back(gr);
   return gr;
}

WFMirrorArray* WFTelescopeArray::GetMirror(int iTel){
   if(!pct) return 0;
   if(iTel<0||iTel>=CTNumber) return 0;
   WFTelescope* p0=pct[iTel];
   if(!p0) return 0;
   return p0->pmirr;
}
SquareCone* WFTelescopeArray::GetCone(int iTel){
   if(!pct) return 0;
   if(iTel<0||iTel>=CTNumber) return 0;
   WFTelescope* p0=pct[iTel];
   if(!p0) return 0;
   return p0->pcone;
}
WCamera* WFTelescopeArray::GetCamera(int iTel){
   if(!pct) return 0;
   if(iTel<0||iTel>=CTNumber) return 0;
   WFTelescope* p0=pct[iTel];
   if(!p0) return 0;
   return p0->pcame;
}
double WFTelescopeArray::GetFocus(int iTel){
   if(iTel==1) return 2880;
   else if(iTel==2) return 2910;
   else if(iTel==3) return 2922;
   else if(iTel==4) return 2942;
   else if(iTel==5) return 2901;
   else if(iTel==6) return 2906;
   else return 2870;
}

double WFTelescope::MirrorSizeX=400.; //cm
double WFTelescope::MirrorSizeY=400.; //cm
double WFTelescope::ZDOOR=(2870+230+80); //mm
double WFTelescope::D_DOOR=2340;
double WFTelescope::Hdoor=2379;
double WFTelescope::TRANSPARENCY=0.87;
double WFTelescope::CLUSTER_X = 894.25;
double WFTelescope::CLUSTER_Y = 769.08;
double WFTelescope::FOCUS = 2870; //2870
double WFTelescope::ZCLUSTER0 = 2870;
double WFTelescope::ZCLUSTER1 =(2870+230);
TGraph2D *WFTelescope::Transmissivity=0;
double WFTelescope::filter_wl_max, WFTelescope::filter_wl_min;
double WFTelescope::filter_angle_max, WFTelescope::filter_angle_min;
double WFTelescope::Reflectivity[55];
TGraph* WFTelescope::quantumeff=0;
double WFTelescope::mirror_wl_max, WFTelescope::mirror_wl_min, WFTelescope::mirror_wl_step;
void WFTelescope::Init(){
   TelIndex_=0;
   Telx_=0;
   Tely_=0;
   Telz_=0;
   TelZ_=0;
   TelA_=0;
   for(int ii=0;ii<3;ii++){
      for(int jj=0;jj<3;jj++){
         matrix_[ii][jj]=0;
      }
   }
   pmirr=new WFMirrorArray();
   pcame=new WCamera();
   pcone=new SquareCone();
   SetTransmissivity();
   SetReflectivity();
   SetQuantumEff();
}
void WFTelescope::Clear(){
   if(pmirr) {delete pmirr; pmirr=0;}
   if(pcame) {delete pcame; pcame=0;}
   if(pcone) {delete pcone; pcone=0;}
   if(Transmissivity) {delete Transmissivity; Transmissivity=0;}
   if(quantumeff) {delete quantumeff; quantumeff=0;}
}
void WFTelescope::SetXYZ(double x, double y, double z)
{
    Telx_ = x;
    Tely_ = y;
    Telz_ = z;
}
void WFTelescope::SetPointing(double zenith,double azimuth)
{
    TelZ_ = zenith;
    TelA_ = azimuth;
    //printf("WFTelescope::SetPointing: TelZ %f TelA %f\n", TelZ_, TelA_);
}
void WFTelescope::SetEulerMatrix(double theta,double phi)
{
   double cosz, sinz, cosa, sina;
   cosz = cos(theta);
   sinz = sin(theta);
   cosa = cos(phi);
   sina = sin(phi);
   matrix_[0][0] = cosa*cosz;
   matrix_[0][1] = sina*cosz;
   matrix_[0][2] = -sinz;
   matrix_[1][0] = -sina;
   matrix_[1][1] = cosa;
   matrix_[1][2] = 0;
   matrix_[2][0] = cosa*sinz;
   matrix_[2][1] = sina*sinz;
   matrix_[2][2] = cosz;
}
void WFTelescope::Euler(double x0, double y0, double z0, double *x, double *y, double *z)
{
   *x = matrix_[0][0]*x0+matrix_[0][1]*y0+matrix_[0][2]*z0;
   *y = matrix_[1][0]*x0+matrix_[1][1]*y0+matrix_[1][2]*z0;
   *z = matrix_[2][0]*x0+matrix_[2][1]*y0+matrix_[2][2]*z0;
}
void WFTelescope::InverseEuler(double x0, double y0, double z0, double *x, double *y, double *z)
{
   *x = matrix_[0][0]*x0+matrix_[1][0]*y0+matrix_[2][0]*z0;
   *y = matrix_[0][1]*x0+matrix_[1][1]*y0+matrix_[2][1]*z0;
   *z = matrix_[0][2]*x0+matrix_[1][2]*y0+matrix_[2][2]*z0;
}
bool WFTelescope::IncidentTel(double x,double y)
{
    if(fabs(x-Telx_)<MirrorSizeX/2 && fabs(y-Tely_)<MirrorSizeY/2) return true;
    else return false;
}
void WFTelescope::SetTransmissivity()
{
  int wl;
  double p1,p2,p3,p4,p5,p6,p7,p8;
  double filter_wl_step = 5;  //nm
  double filter_angle_step = 5; //degree
  FILE *fp;
  filter_wl_min = 300;
  filter_wl_max = 600;
  filter_angle_max = 35;
  filter_angle_min = 0;
  int i=0;
  Transmissivity = new TGraph2D();
 // Transmissivity->SetName("t");
  char DatabaseDir[200]="";
  strcpy(DatabaseDir,getenv("WFCTADataDir"));
  fp = fopen(Form("%s/filter_transparency.txt",DatabaseDir),"r");
  while(!feof(fp)){
     fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %lf %lf\n",&wl,&p1,&p2,&p3,&p4,&p5,&p6,&p7,&p8);
     Transmissivity->SetPoint(i,wl,filter_angle_min+filter_angle_step*0,p1);
     i++;
     Transmissivity->SetPoint(i,wl,filter_angle_min+filter_angle_step*1,p2);
     i++;
     Transmissivity->SetPoint(i,wl,filter_angle_min+filter_angle_step*2,p3);
     i++;
     Transmissivity->SetPoint(i,wl,filter_angle_min+filter_angle_step*3,p4);
     i++;
     Transmissivity->SetPoint(i,wl,filter_angle_min+filter_angle_step*4,p5);
     i++;
     Transmissivity->SetPoint(i,wl,filter_angle_min+filter_angle_step*5,p6);
     i++;
     Transmissivity->SetPoint(i,wl,filter_angle_min+filter_angle_step*6,p7);
     i++;
     Transmissivity->SetPoint(i,wl,filter_angle_min+filter_angle_step*7,p8);
     i++;
  }
  if(Transmissivity->GetN()>0) printf("WFTelescopeArray::SetTransmissivity: load filter transparency from %s successfully\n",DatabaseDir);
  fclose(fp);

}
int WFTelescope::GetTransmissivity(double wavelength, double angle)
{
   // * updated by lingling Ma 2019.5.8
   int i;
   double p;
   double transmissivity;
   if(wavelength<filter_wl_min||wavelength>filter_wl_max||angle>filter_angle_max)  transmissivity = 0;
   else{
      transmissivity = Transmissivity-> Interpolate(wavelength, angle);
      //printf("WFTelescope::GetTransmissivity: wl=%lf angle=%lf transmissivity=%lf\n",wavelength,angle,transmissivity);
   }
   p = WFTelescopeArray::prandom->Rndm();
  if(p<transmissivity) return 1;
  else return 0;
}
bool WFTelescope::GetTransmissivity(double wavelength)
{

   /* The filter transmissivity is simulated for photons with wavelength
 *       from 280 nm to 415 nm with step 5 nm.
 *          */
   int i;
   double wl_max = 415; //nm
   double wl_min = 280; //nm
   double wl_step = 5;  //nm
   double wl0, wl1, t0, t1;

   double transmissivity;
   double Transmissivity[28] ={
  0.000600, 0.010656, 0.029956, 0.086233, 0.161489, 0.265211, 0.368200, 0.470700
, 0.552622, 0.629400, 0.679678, 0.726022, 0.756633, 0.781378, 0.797322, 0.806567
, 0.811344, 0.802433, 0.790100, 0.740322, 0.679389, 0.573833, 0.452278, 0.306722
, 0.188967, 0.089744, 0.044167, 0.012000
   };
  if(wavelength<wl_min||wavelength>wl_max) transmissivity = 0;

   else{
      i = int(wavelength-wl_min)/wl_step;
      wl0 = i*wl_step + wl_min;
      wl1 = (i+1)*wl_step + wl_min;
      t0 = Transmissivity[i];
      t1 = Transmissivity[i+1];
      transmissivity = t0 + (t1-t0)*(wavelength-wl0)/(wl1-wl0);
   }

   double p;
   p = WFTelescopeArray::prandom->Rndm();

  // 1 survived, 0 absorbed by the filter
  if(p<transmissivity) return true;
  else return false;
}
void WFTelescope::SetReflectivity()
{
  FILE *fp;
  double p;
  int wl;
  char DatabaseDir[200]="";
  strcpy(DatabaseDir,getenv("WFCTADataDir"));
  fp = fopen(Form("%s/reverseEfficiency.txt",DatabaseDir),"r");
  int i=0;
  mirror_wl_max = 800;
  mirror_wl_min = 250;
  mirror_wl_step = 10; //nm
  while(!feof(fp)){
     fscanf(fp,"%d %lf\n",&wl,&p);
     Reflectivity[i] = p;
     i++;
  }
  if(i>0) printf("WFTelescopeArray::SetReflectivity: load reflection efficiency from %s successfully\n",DatabaseDir);
  fclose(fp);
}
bool WFTelescope::Reflected(double wavelength)
{
  //Update by Lingling Ma 2019.05.08
  int i = int ((wavelength-mirror_wl_min)/mirror_wl_step);
  int j = i+1;
  double e1, e2, efficiency;
  double wl1, wl2;
  if((i>=0&&i<55)&&(j>=0&&j<55)){
  e1 = Reflectivity[i];
  e2 = Reflectivity[j];
  wl1 = mirror_wl_min+i*mirror_wl_step;
  wl2 = mirror_wl_min+j*mirror_wl_step;
  double k = (e1-e2)/(wl1-wl2);
  efficiency = (wavelength-wl2)*k+e2;
  }
  else if(wavelength<=0) efficiency = 1;
  else if(i<0) efficiency = Reflectivity[0];
  else if(j>54) efficiency = Reflectivity[54];

  //double EFFICIENCY=TRANSPARENCY*WFMirror::REFLECTIVITY*SquareCone::CONE_REFLECTIVE_EFFICIENCY;
  double EFFICIENCY=TRANSPARENCY*efficiency*SquareCone::CONE_REFLECTIVE_EFFICIENCY;
  double x;
  x = WFTelescopeArray::prandom->Rndm();
  if(x<EFFICIENCY)
     return true;
  else
     return false;
}
void WFTelescope::SetQuantumEff(){
  FILE *fp;
  double p;
  double wl;
  char DatabaseDir[200]="";
  strcpy(DatabaseDir,getenv("WFCTADataDir"));
  fp = fopen(Form("%s/quantumeff.txt",DatabaseDir),"r");
  quantumeff=new TGraph();
  int i=0;
  while(!feof(fp)){
     fscanf(fp,"%lf %lf\n",&wl,&p);
     quantumeff->SetPoint(i,wl,p);
     i++;
  }
  if(quantumeff->GetN()>0) printf("WFTelescopeArray::SetQuantumEff: load quantum efficiency from %s successfully\n",DatabaseDir);
  fclose(fp);
}
bool WFTelescope::GetQuantumEff(double wavelength){
  double eff=quantumeff->Eval(wavelength);
  double x = WFTelescopeArray::prandom->Rndm();
  if(x<eff)
     return true;
  else
     return false;
}
void WFTelescope::ConvertCoor(double x0,double y0,double z0,double m1,double n1,double l1,double &x0new,double &y0new,double &z0new,double &m1new,double &n1new,double &l1new){
   double x1=x0-Telx_;
   double y1=y0-Tely_;
   double z1=z0-(Telz_+WFTelescopeArray::lhaaso_coo[2]);
   Euler(x1,y1,z1,&x0new,&y0new,&z0new);
   //from cm to mm;
   x0new*=10;
   y0new*=10;
   z0new*=10;
   Euler(m1,n1,l1,&m1new,&n1new,&l1new);
}
void WFTelescope::Plane(double z,double x0,double y0,double z0,double m2,double n2,double l2,double *x,double *y)
{
  double k;
  k = (z-z0)/l2;
  *x = k*m2+x0;
  *y = k*n2+y0;
}
void WFTelescope::Sphere(double Z0,double R, double x0, double y0, double z0, double m, double n, double l, double *x, double *y, double *z)
{
  double A = 1.;
  double B = 2*m*x0+2*n*y0+2*l*z0-2*l*R;
  double C = x0*x0+y0*y0+z0*z0-2*z0*R;
  double delta = B*B-4*A*C;
  double k1, k2, z1, z2;

  if(delta<0){
    *x = -100000;
    *y = -100000;
    *z = -100000;
  }
  else{
    k1 = (-B + sqrt(delta))/(2*A);
    k2 = (-B - sqrt(delta))/(2*A);

    z1 = l*k1 + z0;
    z2 = l*k2 + z0;
    if(z1>Z0&&z2>Z0){
      *z = 1000000;
      *x = 1000000;
      *y = 1000000;
    }
    else if(z1<Z0&&z2<Z0){
      *z = 1000000;
      *x = 1000000;
      *y = 1000000;
    }
    else{
      if(z1<z2) {
        *z = z1;
        *x = m*k1+x0;
        *y = n*k1+y0;
      }
      if(z2<z1){
        *z = z2;
        *x = m*k2+x0;
        *y = n*k2+y0;
      }
    }
  }
}
bool WFTelescope::RayTraceDoor(double x0, double y0, double z0, double m1, double n1, double l1){
    double x,y;
    if(WFTelescopeArray::jdebug>3) printf("WFTelescope::RayTraceDoor: Above Door, InCoo={%.3f,%.3f,%.3f},InDir{%.3f,%.3f,%.3f}\n",x0,y0,z0,m1,n1,l1);
    Plane(ZDOOR,x0,y0,z0,m1,n1,l1,&x,&y);
    if(WFTelescopeArray::jdebug>3) printf("WFTelescope::RayTraceDoor: On Door, OtCoo={%.3f,%.3f} DoorSize={%.3f,%.3f}\n",x,y,D_DOOR/2,Hdoor/2);

    if(fabs(x)>D_DOOR/2.||fabs(y)>Hdoor/2.) return false; ///< shelded by the door; 13102015
    else return true;
}
bool WFTelescope::RayTraceCluster(double x0, double y0, double z0, double m1, double n1, double l1,bool IsCluster0,bool IsInside,double &x,double &y){
    if(WFTelescopeArray::jdebug>3) printf("WFTelescope::RayTraceCluster(%d,%d): Before Cluster, InCoo={%.3f,%.3f,%.3f},InDir{%.3f,%.3f,%.3f}\n",IsCluster0,IsInside,x0,y0,z0,m1,n1,l1);
    Plane(IsCluster0?ZCLUSTER0:ZCLUSTER1,x0,y0,z0,m1,n1,l1,&x,&y);
    if(WFTelescopeArray::jdebug>3) printf("WFTelescope::RayTraceCluster(%d,%d): Before Cluster, OtCoo={%.3f,%.3f},ClusterSize={%.3f,%.3f}\n",IsCluster0,IsInside,x,y,CLUSTER_X/2,CLUSTER_Y/2);

    if(!IsInside){
       if(fabs(x)<CLUSTER_X/2.&&fabs(y)<CLUSTER_Y/2.) return false; ///< shelded by the cluster; 13102015
       return true;
    }
    else{
       return true;
       if(fabs(x)>CLUSTER_X/2.||fabs(y)>CLUSTER_Y/2.) return false; /////out of the range of the cluster
       return true;
    }
}
bool WFTelescope::RayTraceMirror1(double x0, double y0, double z0, double m1, double n1, double l1,double &xmirror0,double &ymirror0, double &zmirror0){
    if(WFTelescopeArray::jdebug>3) printf("WFTelescope::RayTraceMirror1: Before Mirror, InCoo={%.3f,%.3f,%.3f},InDir{%.3f,%.3f,%.3f}\n",x0,y0,z0,m1,n1,l1);
    double ZMIRROR=WFMirror::CURVATURE - sqrt(WFMirror::CURVATURE*WFMirror::CURVATURE-D_DOOR*D_DOOR/4.);
    Sphere(ZMIRROR,WFMirror::CURVATURE,x0,y0,z0,m1,n1,l1,&xmirror0,&ymirror0,&zmirror0);
    if(WFTelescopeArray::jdebug>3) printf("WFTelescope::RayTraceMirror1: On Mirror, Zmirr,R={%.f,%.f} OtCoo={%.3f,%.3f,%.3f},MirrorXYSize={%.3f,%.3f}\n",ZMIRROR,WFMirror::CURVATURE,xmirror0,ymirror0,zmirror0,D_DOOR/2,Hdoor/2);

    if(fabs(xmirror0)>D_DOOR/2.||fabs(ymirror0)>Hdoor/2.) return false; //shelded by the container's wall
    else return true;
}
bool WFTelescope::RayTraceMirror2(double xmirror0,double ymirror0,double zmirror0,double m1, double n1, double l1,double &m2, double &n2, double &l2,int &ii,int &mm){
    double deltax = -10000;
    double deltay = -10000;
    double deltaz = -10000;
    if(WFTelescopeArray::jdebug>3){
       double delta0[3]={0,0,-1};
       double mindist=10000000;
       int minindex[2]={-1,-1};
       for(int i1=0;i1<NCircle;i1++){
          for(int jj=0;jj<pmirr->NMirror[i1];jj++){
             double dist=sqrt(pow(xmirror0-pmirr->pmirr[i1][jj]->mirrorx,2)+pow(ymirror0-pmirr->pmirr[i1][jj]->mirrory,2)+pow(zmirror0-pmirr->pmirr[i1][jj]->mirrorz,2));
             if(dist<mindist){
                mindist=dist;
                minindex[0]=i1;minindex[1]=jj;
                delta0[0]=xmirror0-pmirr->pmirr[i1][jj]->mirrorx;
                delta0[1]=ymirror0-pmirr->pmirr[i1][jj]->mirrory;
                delta0[2]=zmirror0-pmirr->pmirr[i1][jj]->mirrorz;
             }
          }
       }
       printf("WFTelescope::RayTraceMirror2: TelFrame: XYZ={%.3f,%.3f,%.3f} deltaxyz={%.3f,%.3f,%.3f}\n",xmirror0,ymirror0,zmirror0,delta0[0],delta0[1],delta0[2]);
    }
    pmirr->WhichMirror(xmirror0, ymirror0, zmirror0, &deltax, &deltay, &deltaz, &ii, &mm);
    if(WFTelescopeArray::jdebug>3) printf("WFTelescope::RayTraceMirror2: MirrorFrame: deltaxyz={%.3f,%.3f,%.3f} im={%d,%d}\n",deltax,deltay,deltaz,ii,mm);

    if(deltax==-10000) return false;  ///< out of mirror area

    double xmirror = WFTelescopeArray::prandom->Gaus(xmirror0,WFMirror::MirrorSpot);
    double ymirror = WFTelescopeArray::prandom->Gaus(ymirror0,WFMirror::MirrorSpot);
    double zmirror = zmirror0;
    /*//if(deltax==-10000){
       double x,y,z;
       x=0-xmirror;
       y=0-ymirror;
       z=WFMirror::CURVATURE-zmirror;
       double norm=sqrt(x*x+y*y+z*z);
       x/=norm;
       y/=norm;
       z/=norm;
       double costheta = m1*x+n1*y+l1*z;
       m2 = (m1-2.*costheta*x);
       n2 = (n1-2.*costheta*y);
       l2 = (l1-2.*costheta*z);
       return true;
    //}*/
    pmirr->pmirr[ii][mm]->GetReflected(m1,n1,l1,xmirror,ymirror,zmirror,&m2,&n2,&l2);
    if(WFTelescopeArray::jdebug>3) printf("WFTelescope::RayTraceMirror2: Reflection: RefCoo={%.3f,%.3f,%.3f} InDir={%.3f,%.3f,%.3f} OtDir={%.3f,%.3f,%.3f}\n",xmirror,ymirror,zmirror,m1,n1,l1,m2,n2,l2);

    return true;
}
int WFTelescope::RayTraceUpToCone(double x0, double y0, double z0, double m1, double n1, double l1,double &t,double &xcluster,double &ycluster,double &m2,double &n2,double &l2){
    double x,y;
    if(!RayTraceDoor(x0,y0,z0,m1,n1,l1)) return -1;
    if(WFTelescopeArray::jdebug>2) printf("WFTelescope::RayTraceUpToCone: Go Through Door\n");

    if(!RayTraceCluster(x0,y0,z0,m1,n1,l1,false,false,x,y)) return -2;
    if(!RayTraceCluster(x0,y0,z0,m1,n1,l1,true,false,x,y)) return -2;
    if(WFTelescopeArray::jdebug>2) printf("WFTelescope::RayTraceUpToCone: Go Through Cluster\n");

    double xmirror0,ymirror0,zmirror0; ///<impact point on the mirror
    if(!RayTraceMirror1(x0,y0,z0,m1,n1,l1,xmirror0,ymirror0,zmirror0)) return -3;
    if(WFTelescopeArray::jdebug>2) printf("WFTelescope::RayTraceUpToCone: Go Through Mirror1\n");

    int ii, mm;   ///< mirror index
    if(!RayTraceMirror2(xmirror0,ymirror0,zmirror0,m1,n1,l1,m2,n2,l2,ii,mm)) return -4;
    if(WFTelescopeArray::jdebug>2) printf("WFTelescope::RayTraceUpToCone: Go Through Mirror2\n");

    if(!RayTraceCluster(xmirror0,ymirror0,zmirror0,m2,n2,l2,true,true,xcluster,ycluster)) return -5;
    if(WFTelescopeArray::jdebug>2) printf("WFTelescope::RayTraceUpToCone: Go Through Cluster2\n");

    double dist = (xmirror0-x0)*(xmirror0-x0) + (ymirror0-y0)*(ymirror0-y0) + (zmirror0-z0)*(zmirror0-z0);
    dist = sqrt(dist);
    dist = sqrt((xmirror0-xcluster)*(xmirror0-xcluster)
         + (ymirror0-ycluster)*(ymirror0-ycluster)
         + (zmirror0-ZCLUSTER0)*(zmirror0-ZCLUSTER0))
         - dist;
    //printf("tcal: t1=%.9le t2=%.9le(dist=%lfmm dt=%.9le)\n\n",t,t+dist*0.1/vlight,dist,dist*0.1/vlight);
    t += dist*0.1/vlight;  //in second
    return 1;
}
int WFTelescope::RayTrace(double x0, double y0, double z0, double m1, double n1, double l1,double wavelength,double &t,int &itube,int &icell)
{
    double xcluster,ycluster;  ///< x,y on the pmt cluster plane
    double m2, n2, l2; ///< the direction of the reflected ray
    double deltax,deltay,deltaz;
    double Weight;

    int icone;
    int ix,iy;
    double x,y;
    double dircos[3];
    double hitpos[2];

    if(WFTelescopeArray::jdebug>1) printf("WFTelescope::RayTrace: Passing From Door to Cone InCoo={%.f,%.f,%.f} InDir={%.f,%.f,%.f}\n",x0,y0,z0,m1,n1,l1);
    double result=RayTraceUpToCone(x0,y0,z0,m1,n1,l1,t,xcluster,ycluster,m2,n2,l2);
    //l2=sqrt(1-m2*m2-n2*n2);
    if(!(l2>=0)) l2=0;
    if(result<0) return result;
    if(WFTelescopeArray::jdebug>1) printf("WFTelescope::RayTrace: Go Through From Door to Cone OtCoo={%.f,%.f} OtDir={%.f,%.f,%.f}\n",xcluster,ycluster,m2,n2,l2);

    ///pass filter
    if(wavelength>0){
       double angle=acos(fabs(l2))*TMath::RadToDeg();
       if(WFTelescopeArray::jdebug>0) printf("WFTelescope::RayTrace: Passing Filter wl=%lf angle=%lf\n",wavelength,angle);
       if(!GetTransmissivity(wavelength,angle)) return -6;
       if(WFTelescopeArray::jdebug>0) printf("WFTelescope::RayTrace: Go Through Filter\n");
    }

    //from Tel Coo to Cone Coo.
    double u,v,w,xc,yc;
    u = n2; //m2
    v = m2; //n2
    w = -l2;
    xc = ycluster; //-ycluster
    yc = xcluster;

    //*The photons that can enter the wenston cone*// 
    if(WFTelescopeArray::jdebug>1) printf("WFTelescope::RayTrace: Before Cone InCoo={%.f,%.f}\n",xc,yc);
    icone = pcame->GetCone(xc,yc);
    if(icone<0) return -7;
    //else{itube=icone; return 1;}
    if(WFTelescopeArray::jdebug>1) printf("WFTelescope::RayTrace: Inside Cone%d\n",icone);

    //*To Get the cone coordinates*//
    x = pcame->GetSiPMX(icone);
    y = pcame->GetSiPMY(icone);

    deltax =  xc - x;
    deltay =  yc - y;
    dircos[0] = u; //-u
    dircos[1] = v; //v
    dircos[2] = w; //-sqrt(1-u*u-v*v)

    if(WFTelescopeArray::jdebug>1) printf("WFTelescope::RayTrace: Passinfg Cone%d delta={%f.%f} dircos={%f,%f,%f}\n",icone,deltax,deltay,dircos[0],dircos[1],dircos[2]);
    //*The ray trace in the cone *//
    pcone -> SetInitDir(dircos);
    pcone -> SetInitPos(deltax,deltay);
    pcone -> SquareRaytracing();
    pcone -> GetHitPos(hitpos,Weight);
    if(!pcone -> GetStrike()) return -8;
    if(WFTelescopeArray::jdebug>1) printf("WFTelescope::RayTrace: Cone%d After GetStrike\n",icone);
    hitpos[0] += x;
    hitpos[1] += y;
    itube = pcame->GetTube(hitpos[0],hitpos[1]);
    if(itube<0) return -9;           //Throw away the bad point after ConeTracing.
    if(WFTelescopeArray::jdebug>1) printf("WFTelescope::RayTrace: Inside Tube%d hitpos={%f,%f}\n",itube,hitpos[0],hitpos[1]);
    if(itube!=icone) return -9;
    if(WFTelescopeArray::jdebug>1) printf("WFTelescope::RayTrace: Tube%d=Cone%d\n",itube,icone);
    double xx = WFTelescopeArray::prandom->Rndm();
    if(xx>Weight) return -10;

    hitpos[0] = hitpos[0] - x;
    hitpos[1] = hitpos[1] - y;
    ix = int((hitpos[0]+WCamera::D_SiPM/2)/WCamera::D_Cell);
    iy = int((hitpos[1]+WCamera::D_SiPM/2)/WCamera::D_Cell);
    icell = iy*WCamera::NCell+ix;
    if(WFTelescopeArray::jdebug>1) printf("WFTelescope::RayTrace: Passing All to SiPM Tube%d Cell%d\n",itube,icell);

    return TelIndex_>0?TelIndex_:1;
}
bool WFTelescope::GetImageCoo(double dir_in[3],double &xx,double &yy,bool IsLocal){
   double norm=sqrt(dir_in[0]*dir_in[0]+dir_in[1]*dir_in[1]+dir_in[2]*dir_in[2]);
   double theta=IsLocal?(PI/2):(PI/2-TelZ_);
   double phi=IsLocal?0:TelA_;
   double zz0=(dir_in[0]*cos(phi)+dir_in[1]*sin(phi))/norm*cos(theta)+dir_in[2]/norm*sin(theta);
   if(zz0>=0) return false;
   double xx0=(dir_in[0]*cos(phi)+dir_in[1]*sin(phi))/norm*sin(theta)-dir_in[2]/norm*cos(theta);
   double yy0=(-dir_in[0]*sin(phi)+dir_in[1]*cos(phi))/norm;
   xx=-yy0/zz0/PI*180;
   yy=-xx0/zz0/PI*180;
   return true;
}
void WFTelescope::GetOutDir(double imagexy[2],double dir_out[3],bool IsLocal){
   double theta=IsLocal?(PI/2):(PI/2-TelZ_);
   double phi=IsLocal?0:TelA_;
   double xx0=imagexy[1]/180.*PI;
   double yy0=imagexy[0]/180.*PI;
   dir_out[0]=(-xx0*sin(theta)+cos(theta))*cos(phi)+yy0*sin(phi);
   dir_out[1]=(-xx0*sin(theta)+cos(theta))*sin(phi)-yy0*cos(phi);
   dir_out[2]=xx0*cos(theta)+sin(theta);
   double norm=sqrt(dir_out[0]*dir_out[0]+dir_out[1]*dir_out[1]+dir_out[2]*dir_out[2]);
   for(int ii=0;ii<3;ii++) dir_out[ii]/=norm;
}
