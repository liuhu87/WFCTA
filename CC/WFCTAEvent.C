#include <stdlib.h>
#include <iostream>
#include "WFCTAEvent.h"
#include "WFCamera.h"
#include "TH2Poly.h"
#include "TMarker3DBox.h"
#include "TAxis3D.h"
#include <TCanvas.h>
#include <TView3D.h>
#include <TSystem.h>
#include "TF1.h"
#include "Laser.h"

using namespace std;

ClassImp(WFCTAEvent);

WFCTAEvent* WFCTAEvent::_Head=0;
TTree* WFCTAEvent::_Tree=0;
const char* WFCTAEvent::_Name="Event";
TBranch* WFCTAEvent::bAll=0;
TBranch* WFCTAEvent::bmcevent=0;
TBranch* WFCTAEvent::bledevent=0;
TBranch* WFCTAEvent::blaserevent=0;
WFCTAEvent::WFCTAEvent():TSelector()
{
   ievent.reserve(MAXPMT);
   ADC_Cut.reserve(MAXPMT);
   ImageBaseHigh.reserve(MAXPMT);
   ImageBaseLow.reserve(MAXPMT);
   ImageAdcHigh.reserve(MAXPMT);
   ImageAdcLow.reserve(MAXPMT);
   myImageBaseHigh.reserve(MAXPMT);
   myImageBaseLow.reserve(MAXPMT);
   myImageAdcHigh.reserve(MAXPMT);
   myImageAdcLow.reserve(MAXPMT);
   iSiPM.reserve(MAXPMT);
   Single_Threshold.reserve(MAXPMT);
   Record_Threshold.reserve(MAXPMT);
   peak.reserve(MAXPMT);
   mypeak.reserve(MAXPMT);
   peakamp.reserve(MAXPMT);
   gain_marker.reserve(MAXPMT);
   Over_Single_Marker.reserve(MAXPMT);
   Over_Record_Marker.reserve(MAXPMT);

   ADC_Cut.resize(MAXPMT);
   ImageBaseHigh.resize(MAXPMT);
   ImageBaseLow.resize(MAXPMT);
   ImageAdcHigh.resize(MAXPMT);
   ImageAdcLow.resize(MAXPMT);
   myImageBaseHigh.resize(MAXPMT);
   myImageBaseLow.resize(MAXPMT);
   myImageAdcHigh.resize(MAXPMT);
   myImageAdcLow.resize(MAXPMT);
   iSiPM.resize(MAXPMT);
   Single_Threshold.resize(MAXPMT);
   Record_Threshold.resize(MAXPMT);
   peak.resize(MAXPMT);
   mypeak.resize(MAXPMT);
   peakamp.resize(MAXPMT);
   gain_marker.resize(MAXPMT);
   Over_Single_Marker.resize(MAXPMT);
   Over_Record_Marker.resize(MAXPMT);

   Init();
}

WFCTAEvent::~WFCTAEvent()
{
   EventInitial();
   if(fModel) delete fModel;
   if(gDrawErr) delete gDrawErr;
   if(minimizer) delete minimizer;
}

void WFCTAEvent::Init()
{
   iTel=-1;
   iEvent=-1;
   rabbitTime=0;
   rabbittime=0;
   big_pack_lenth=-1;
   n_fired=-1;
   n_Channel=-1;
   iSiPM.clear();
   ievent.clear();
   gain_marker.clear();
   peak.clear();
   mypeak.clear();
   peakamp.clear();
   Single_Threshold.clear();
   Record_Threshold.clear();
   Over_Single_Marker.clear();
   Over_Record_Marker.clear();
   ADC_Cut.clear();
   ImageBaseHigh.clear();
   ImageBaseLow.clear();
   ImageAdcHigh.clear();
   ImageAdcLow.clear();
   myImageBaseHigh.clear();
   myImageBaseLow.clear();
   myImageAdcHigh.clear();
   myImageAdcLow.clear();

   for(int j=0;j<28;j++){
     Npoint[j]=j;
     for(int i=0;i<1024;i++){
       pulsehigh[i][j] = 0;
       pulselow[i][j] = 0;
     }
   }

   fModel=new TF1("ImageModel","[0]+[1]*x",-10,10);
   fModel->SetNpx(1000);
   fModel->SetLineColor(4);
   fModel->SetLineWidth(4);
   gDrawErr=0;
   minimizer=0;
  
   mcevent.Init();
   ledevent.Init();
   laserevent.Init();
}
void WFCTAEvent::EventInitial()
{
   iTel=-1;
   iEvent=-1;
   rabbitTime=0;
   rabbittime=0;
   big_pack_lenth=-1;
   n_fired=-1;
   n_Channel=-1;
   iSiPM.clear();
   ievent.clear();
   gain_marker.clear();
   peak.clear();
   mypeak.clear();
   peakamp.clear();
   Single_Threshold.clear();
   Record_Threshold.clear();
   Over_Single_Marker.clear();
   Over_Record_Marker.clear();
   ADC_Cut.clear();
   ImageBaseHigh.clear();
   ImageBaseLow.clear();
   ImageAdcHigh.clear();
   ImageAdcLow.clear();
   myImageBaseHigh.clear();
   myImageBaseLow.clear();
   myImageAdcHigh.clear();
   myImageAdcLow.clear();

   for(int j=0;j<28;j++){
     for(int i=0;i<1024;i++){
       pulsehigh[i][j] = 0;
       pulselow[i][j] = 0;
     }
   }

   fModel->SetParameters(0,0);
   if(minimizer) {delete minimizer; minimizer=0;}
   if(gDrawErr) {delete gDrawErr; gDrawErr=0;}
  
   mcevent.Reset();
   ledevent.Reset();
   laserevent.Reset();
}

void WFCTAEvent::InitTree(TTree *tree){
  //   Set branch addresses
  if (tree == 0) return;
  Tree() = tree;
  Head() = this;
  //Tree()->SetMakeClass(1);
  Tree()->SetBranchAddress(BranchName(),&Head());
}
void WFCTAEvent::CreateBranch(TTree *tree, int branchSplit){
   if(tree){
     Head()=this;
     tree->Branch(BranchName(),"WFCTAEvent",&Head(),32000,branchSplit);
     TBranch * branch=tree->GetBranch(BranchName());
     int clevel=branch->GetCompressionLevel();
     #ifdef __LZMA__
     if(clevel<6 && branch->GetCompressionAlgorithm()!=ROOT::kLZMA)clevel=6;
     #else
     if(clevel<6)clevel=6;
     #endif
     branch->SetCompressionLevel(clevel);

     cout <<" CompressionLevel "<<branch->GetCompressionLevel()<<" "<<branch->GetSplitLevel()<<endl;
     //tree->SetBranchStatus("TSelector",false);
     //tree->SetBranchStatus("mcevent",false);
     //tree->SetBranchStatus("ledevent",false);
     //tree->SetBranchStatus("laserevent",false);
   }
}
void WFCTAEvent::GetBranch(TTree *fChain){
   bAll=fChain->GetBranch(BranchName());
   bmcevent=fChain->GetBranch("mcevent");
   bledevent=fChain->GetBranch("ledevent");
   blaserevent=fChain->GetBranch("laserevent");
}
bool WFCTAEvent::GetAllContents(int _Entry){
   EventInitial();
   bool exist=true;
   if(!bAll) exist=false;
   //if(!bmcevent) exist=false;
   //if(!bledevent) exist=false;
   //if(!blaserevent) exist=false;
   if(!exist) return false;
   int ncount=bAll->GetEntry(_Entry);
   if(bmcevent) bmcevent->GetEntry(_Entry);
   if(bledevent) bledevent->GetEntry(_Entry);
   if(blaserevent) blaserevent->GetEntry(_Entry);
   return ncount>0;
}

void WFCTAEvent::CalculateADC(int itel){
   if(itel<0||itel>=NCTMax) return;
   for(int ii=0;ii<NSIPM;ii++){
      if(mcevent.TubeSignal[itel][ii]>0){
         iSiPM.push_back(ii);
         ImageAdcHigh.push_back(mcevent.TubeSignal[itel][ii]*WFCTAMCEvent::fAmpHig);
         ImageAdcLow.push_back(mcevent.TubeSignal[itel][ii]*WFCTAMCEvent::fAmpLow);
         myImageAdcHigh.push_back(mcevent.TubeSignal[itel][ii]*WFCTAMCEvent::fAmpHig);
         myImageAdcLow.push_back(mcevent.TubeSignal[itel][ii]*WFCTAMCEvent::fAmpLow);
      }
   }
}
int WFCTAEvent::GetMaxADCBin(int itel){
   int res=-1;
   double maxadc=-1;
   for(int ii=0;ii<iSiPM.size();ii++){
      double yy=ImageAdcHigh.at(ii);
      if(yy<=0) continue;
      if(yy>maxadc){
         maxadc=yy;
         res=ii;
      }
   }
   return res;
}
int WFCTAEvent::GetPeakADCBin(int isipm,int itel){
   if(itel<0||itel>=WFTelescopeArray::CTNumber) return -3;
   if(isipm<0||isipm>=NSIPM) return -2;
   double maxcontent=-1;
   int index=-1;
   for(int ii=0;ii<mcevent.NArrival[itel];ii++){
      double content=mcevent.ArrivalCount[itel][isipm][ii];
      if(content>maxcontent){
         maxcontent=content;
         index=ii;
      }
   }
   return index;
}
int WFCTAEvent::GetMaxTimeBin(int itel){
   if(itel<0||itel>=WFTelescopeArray::CTNumber) return -2;
   int res=-1;
   double maxtime=-1.e40;
   for(int ii=0;ii<mcevent.NArrival[itel];ii++){
      double yy=mcevent.ArrivalTime[itel][ii];
      if(yy>maxtime){
         maxtime=yy;
         res=ii;
      }
   }
   return res;
}
int WFCTAEvent::GetMinTimeBin(int itel){
   if(itel<0||itel>=WFTelescopeArray::CTNumber) return -2;
   int res=-1;
   double mintime=1.e40;
   for(int ii=0;ii<mcevent.NArrival[itel];ii++){
      double yy=mcevent.ArrivalTime[itel][ii];
      if(yy<mintime){
         mintime=yy;
         res=ii;
      }
   }
   return res;
}

bool WFCTAEvent::CleanImage(int isipm,int itel,int type){
   if(itel<0||itel>=WFTelescopeArray::CTNumber) return false;
   if(isipm<0||isipm>=NSIPM) return false;
   bool res=false;
   double trig0=(type==1||type==3)?(25*WFCTAMCEvent::fAmpHig):(25*WFCTAMCEvent::fAmpLow);
   double trig1=(type==1||type==3)?(45*WFCTAMCEvent::fAmpHig):(45*WFCTAMCEvent::fAmpLow);
   int size=iSiPM.size();
   for(int ii=0;ii<size;ii++){
      int isipm0=iSiPM.at(ii);
      if(isipm0!=isipm) continue;
      double content=(type==1||type==3)?myImageAdcHigh.at(ii):myImageAdcLow.at(ii);
      if(ADC_Cut.size()>ii) {if(ADC_Cut.at(ii)<trig0) continue;}
      else{if(content<trig1) continue;}
      double ImageXi=0,ImageYi=0;
      ImageXi=WCamera::GetSiPMX(ii)/WFTelescope::FOCUS/PI*180;
      ImageYi=WCamera::GetSiPMY(ii)/WFTelescope::FOCUS/PI*180;
      int nneigh=0;
      for(int jj=0;jj<size;jj++){
         if(jj==ii) continue;
         double ImageXj=0,ImageYj=0;
         ImageXj=WCamera::GetSiPMX(jj)/WFTelescope::FOCUS/PI*180;
         ImageYj=WCamera::GetSiPMY(jj)/WFTelescope::FOCUS/PI*180;
         double dist=sqrt(pow(ImageXi-ImageXj,2)+pow(ImageYi-ImageYj,2));
         if(dist>0.6) continue;
         double contentj=(type==1||type==3)?myImageAdcHigh.at(jj):myImageAdcLow.at(jj);
         if(ADC_Cut.size()>jj) {if(ADC_Cut.at(jj)<trig0) continue;}
         else{if(contentj<trig1) continue;}
         nneigh++;
      }
      if(nneigh<2) continue;
      res=true;
   }
   return res;
}
double WFCTAEvent::Interface(const double* par){
   fModel->SetParameter(0,par[2]);
   fModel->SetParameter(1,par[3]);
   int size=iSiPM.size();
   double chi2=0;
   int ndof=0;
   double sum=0;
   for(int ii=0;ii<size;ii++){
      int isipm=iSiPM.at(ii);
      if(!CleanImage(isipm,(int)(par[0]+0.5),(int)(par[1]+0.5))) continue;
      double ImageX=WCamera::GetSiPMX(isipm)/WFTelescope::FOCUS/PI*180;
      double ImageY=WCamera::GetSiPMY(isipm)/WFTelescope::FOCUS/PI*180;
      double content=(par[1]==1||par[1]==3)?myImageAdcHigh.at(ii):myImageAdcLow.at(ii);
      chi2+=pow((fModel->Eval(ImageX)-ImageY)/sqrt(1+par[3]*par[3])/0.25,2)*content;
      ndof++;
      sum+=content;
   }
   if(sum>0){
      chi2/=(sum/ndof);
   }
   return chi2;
}
bool WFCTAEvent::DoFit(int itel,int type,bool force){
   if(itel<0||itel>=WFTelescopeArray::CTNumber) return false;
   if(minimizer&&(!force)) return true;
   int size=iSiPM.size();
   int nbin=0;
   double mx,my,sx,sy,sxy,nn;
   for(int ii=0;ii<size;ii++){
      int isipm0=iSiPM.at(ii);
      if(!CleanImage(isipm0,itel,type)) continue;
      double content=(type==1||type==3)?myImageAdcHigh.at(ii):myImageAdcLow.at(ii);
      double ImageX,ImageY;
      ImageX=WCamera::GetSiPMX(isipm0)/WFTelescope::FOCUS/PI*180;
      ImageY=WCamera::GetSiPMY(isipm0)/WFTelescope::FOCUS/PI*180;

      mx += ImageX*content;
      my += ImageY*content;
      sx += ImageX*ImageX*content;
      sy += ImageY*ImageY*content;
      sxy += ImageX*ImageY*content;
      nn += content;
   }

   if(nn>0){
      mx/=nn;
      my/=nn;
      sx/=nn;
      sy/=nn;
      sxy/=nn;
   }
   else return false;
   double cx = sx - mx*mx;
   double cy = sy - my*my;
   double cxy = sxy - mx*my;
   double a = (cy-cx+sqrt((cy-cx)*(cy-cx)+4*cxy*cxy))/(2*cxy);
   double b = my-a*mx;
   //b = (my-DSourceY)-a*(mx-DSourceX);
   double ssx = (cx+2*a*cxy+a*a*cy)/(1+a*a);
   double ssy = (a*a*cx-2*a*cxy+cy)/(1+a*a);

   if(!minimizer) minimizer=ROOT::Math::Factory::CreateMinimizer("Minuit","Migrad");
   minimizer->Clear();
   minimizer->SetMaxFunctionCalls(1000000);
   minimizer->SetMaxIterations(100000);
   minimizer->SetTolerance(0.001);
   minimizer->SetPrintLevel(0);
   #if defined(__CINT__)
   ROOT::Math::Functor f(this,"WFCTAEvent","Interface");
   #else
   ROOT::Math::Functor f(this,&WFCTAEvent::Interface,4);
   #endif
   minimizer->SetFunction(f);
   fModel->SetParameters(b,a);
   //printf("LinearFit: kk=%lf bb=%lf\n",a,b);
   minimizer->SetFixedVariable(0,"iTel",0);
   minimizer->SetFixedVariable(1,"Type",3);
   minimizer->SetLimitedVariable(2,"bb",b,fabs(0.01*b),-1000,1000);
   minimizer->SetLimitedVariable(3,"kk",a,fabs(0.01*a),-1000,1000);
   minimizer->Minimize();
   minimizer->Hesse();
   return true;
}
bool WFCTAEvent::GetPlane(double xyzdir[3][3],double exyzdir[3][3],int itel,int type){
   WFTelescopeArray* pct=WFTelescopeArray::GetHead();
   if(!pct) return false;
   WFTelescope* pt=pct->pct[itel];
   if(!pt) return false;
   double phi0=pt->TelA_;
   double theta=PI/2-pt->TelZ_;
   if(!DoFit(itel,type)) return false;
   double fitpars[2]={minimizer->X()[2],minimizer->X()[3]};
   double AA=-1;
   double BB=fitpars[1];
   double CC=fitpars[0]/180.*PI;
   double zdir[3],xdir[3],ydir[3];
   double coeff[3][3];
   coeff[0][0]=sin(theta)*cos(phi0);
   coeff[0][1]=-sin(phi0);
   coeff[0][2]=-cos(theta)*cos(phi0)/180.*PI;
   coeff[1][0]=sin(theta)*sin(phi0);
   coeff[1][1]=cos(phi0);
   coeff[1][2]=-cos(theta)*sin(phi0)/180.*PI;
   coeff[2][0]=-cos(theta);
   coeff[2][1]=0;
   coeff[2][2]=-sin(theta)/180.*PI;
   //zdir[0]=AA*sin(theta)*cos(phi0)-BB*sin(phi0)-CC*cos(theta)*cos(phi0);
   //zdir[1]=AA*sin(theta)*sin(phi0)+BB*cos(phi0)-CC*cos(theta)*sin(phi0);
   //zdir[2]=-AA*cos(theta)-CC*sin(theta);
   zdir[0]=coeff[0][0]*AA+coeff[0][1]*BB+coeff[0][2]/PI*180.*CC;
   zdir[1]=coeff[1][0]*AA+coeff[1][1]*BB+coeff[1][2]/PI*180.*CC;
   zdir[2]=coeff[2][0]*AA+coeff[2][1]*BB+coeff[2][2]/PI*180.*CC;
   xdir[0]=-zdir[1];
   xdir[1]=zdir[0];
   xdir[2]=0;
   double cosangle=xdir[0]*cos(phi0)+xdir[1]*sin(phi0);
   if(cosangle<0){
      xdir[0]*=-1;
      xdir[1]*=-1;
   }
   if(zdir[0]==0&&zdir[1]==0){
      xdir[0]=cos(phi0)>=0?1:-1;
      xdir[1]=0;
   }
   Laser::cross(zdir,xdir,ydir);
   if(ydir[2]<0){
      for(int ii=0;ii<3;ii++){
         zdir[ii]*=-1;
         ydir[ii]*=-1;
      }
   }

   for(int ii=0;ii<3;ii++){
      xyzdir[0][ii]=xdir[ii];
      xyzdir[1][ii]=ydir[ii];
      xyzdir[2][ii]=zdir[ii];
   }
   for(int idir=0;idir<3;idir++){
      if(idir==1) continue;
      for(int icoo=0;icoo<3;icoo++){
         double coeff1;
         double coeff2;
         if(idir==2){
            coeff1=coeff[icoo][1];
            coeff2=coeff[icoo][2];
         }
         else if(idir==0){
            coeff1=(icoo==2)?0:(icoo==0?-coeff[1][1]:coeff[0][1]);
            coeff2=(icoo==2)?0:(icoo==0?-coeff[1][2]:coeff[0][2]);
         }
         exyzdir[idir][icoo]=sqrt(pow(coeff1,2)*minimizer->CovMatrix(3,3)+pow(coeff2,2)*minimizer->CovMatrix(2,2)+2*coeff1*coeff2*minimizer->CovMatrix(2,3));
      }
   }
   //calculate the error for ydir
   exyzdir[1][0]=sqrt(pow(xyzdir[2][0]*exyzdir[2][2],2)+pow(exyzdir[2][0]*xyzdir[2][2],2));
   exyzdir[1][1]=sqrt(pow(xyzdir[2][1]*exyzdir[2][2],2)+pow(exyzdir[2][1]*xyzdir[2][2],2));
   exyzdir[1][2]=sqrt(pow(2*xyzdir[2][0]*exyzdir[2][0],2)+pow(2*exyzdir[2][1]*xyzdir[2][1],2));
   //normalization
   for(int idir=0;idir<3;idir++){
      double norm=sqrt(xyzdir[idir][0]*xyzdir[idir][0]+xyzdir[idir][1]*xyzdir[idir][1]+xyzdir[idir][2]*xyzdir[idir][2]);
      for(int icoo=0;icoo<3;icoo++){
         xyzdir[idir][icoo]/=norm;
         exyzdir[idir][icoo]/=norm;
      }
   }
   return true;
}
TGraph* WFCTAEvent::DrawAxis(int iaxis,int itel,int type){
   if(iaxis<1||iaxis>3) return 0;
   double xyzdir[3][3];
   double exyzdir[3][3];
   if(!GetPlane(xyzdir,exyzdir,itel,type)) return 0;
   TGraph* gr=new TGraph();
   int np=100;
   for(int ii=0;ii<np;ii++){
      double xx=exp(log(0.5)+(log(1.0e7/0.5))/np*ii);
      double ymin=xx/xyzdir[iaxis-1][0]*xyzdir[iaxis-1][1];
      if(xyzdir[iaxis-1][0]==0){
         ymin=xx;
         xx=0;
      }
      double ymax=ymin;
      for(int jj=0;jj<4;jj++){
         double yy=xx/(xyzdir[iaxis-1][0]+exyzdir[iaxis-1][0]*(jj/2==0?1:-1))*(xyzdir[iaxis-1][1]+exyzdir[iaxis-1][1]*((jj%2)==0?1:-1));
         if(yy<ymin) ymin=yy;
         if(yy>ymax) ymax=yy;
      }
      gr->SetPoint(ii,xx,ymin);
      gr->SetPoint(2*np-1-ii,xx,ymax);
   }
   gr->SetFillStyle(3001);
   gr->SetFillColor(4);
   gr->Draw("F");
   return gr;
}

void WFCTAEvent::DrawFit(){
   if(!minimizer) return;
   double fitpars[2]={minimizer->X()[2],minimizer->X()[3]};
   double fitparse[2]={minimizer->Errors()[2],minimizer->Errors()[3]};
   if(gDrawErr) delete gDrawErr;
   gDrawErr=new TGraph();
   int np=100;
   for(int ii=0;ii<np;ii++){
      double xx=(-8.5+(8.5*2)/np*ii);
      fModel->SetParameters(fitpars[0],fitpars[1]);
      double ymin=fModel->Eval(xx);
      double ymax=ymin;
      for(int jj=0;jj<4;jj++){
         double pars[2];
         pars[0]=fitpars[0]+fitparse[0]*((jj/2==0)?1:-1);
         pars[1]=fitpars[1]+fitparse[1]*((jj%2)==0?1:-1);
         fModel->SetParameters(pars[0],pars[1]);
         double yy=fModel->Eval(xx);
         if(yy<ymin) ymin=yy;
         if(yy>ymax) ymax=yy;
      }
      gDrawErr->SetPoint(ii,xx,ymin);
      gDrawErr->SetPoint(2*np-1-ii,xx,ymax);
   }
   gDrawErr->SetFillStyle(3001);
   gDrawErr->SetFillColor(4);
   gDrawErr->Draw("F");
   fModel->SetParameters(fitpars[0],fitpars[1]);
   fModel->Draw("same");
}
TH2Poly* WFCTAEvent::Draw(int type,const char* opt,double threshold){
   TH2Poly* image=new TH2Poly();
   image->SetName("DrawPlot");
   image->SetTitle(";X;Y");
   for(int ii=0;ii<NSIPM;ii++){
      double ImageX,ImageY;
      ImageX=WCamera::GetSiPMX(ii)/WFTelescope::FOCUS/PI*180;
      ImageY=WCamera::GetSiPMY(ii)/WFTelescope::FOCUS/PI*180;
      image->AddBin(ImageX-0.25,ImageY-0.25,ImageX+0.25,ImageY+0.25);
   }
   for(int ii=0;ii<iSiPM.size();ii++){
      if(!CleanImage(iSiPM.at(ii),0,3)) continue;
      double content=0;
      if(type==0) content=ADC_Cut.at(ii)>threshold?1.:0.;
      if(type==1) content=ImageAdcHigh.at(ii)/WFCTAMCEvent::fAmpHig;
      if(type==2) content=ImageAdcLow.at(ii)/WFCTAMCEvent::fAmpLow;
      if(type==3) content=myImageAdcHigh.at(ii)/WFCTAMCEvent::fAmpHig;
      if(type==4) content=myImageAdcLow.at(ii)/WFCTAMCEvent::fAmpLow;
      image->SetBinContent(iSiPM.at(ii)+1,content>0?content:0);
   }
   image->Draw(opt);
   DrawFit();
   return image;
}
void WFCTAEvent::slaDtp2s(double xi, double eta, double raz, double decz, double &ra, double &dec ){
  double sdecz, cdecz, denom;

  sdecz = sin ( decz );
  cdecz = cos ( decz );
  denom = cdecz - eta * sdecz;
  ra = ( atan2 ( xi, denom ) + raz );
  dec = atan2 ( sdecz + eta * cdecz, sqrt ( xi * xi + denom * denom ) );
  //if(ra<0) ra=ra+2*PI;
  //printf("xi=%lf eta=%lf raz=%lf decz=%lf ra=%lf dec=%lf\n",xi/PI*180,eta/PI*180,raz/PI*180,decz/PI*180,ra/PI*180,dec/PI*180);
}
TH2Poly* WFCTAEvent::DrawGlobal(int type,const char* opt,double threshold){
   WFTelescopeArray* pct=WFTelescopeArray::GetHead();
   if(!pct) return 0;
   WFTelescope* pt=pct->pct[0];
   if(!pt) return 0;
   double raz=pt->TelA_;
   double decz=PI/2-pt->TelZ_;
   TH2Poly* image=new TH2Poly();
   image->SetName("DrawGlobalPlot");
   image->SetTitle(";Azimuth [degree];Elevation [degree]");
   double xrange[2]={1.0e5,-1.0e5};
   double yrange[2]={1.0e5,-1.0e5};
   for(int ii=0;ii<NSIPM;ii++){
      int PixI=ii/PIX;
      int PixJ=ii%PIX;
      double ImageX,ImageY;
      if(PixI%2==0) ImageX=PixJ+0.5-PIX/2.0;
      else ImageX=PixJ+1.0-PIX/2.0;
      ImageY=(PIX/2.0-PixI)-1/2.0;

      ImageX=ImageX*25.4/2870/PI*180;
      ImageY=ImageY*25.4/2870/PI*180;
      //ImageX-=0.31;
      //ImageY-=0.28;

      double theta0[5],phi0[5];
      for(int i2=0;i2<4;i2++){
         double xx,yy;
         if(i2==0){
            xx=ImageX-0.25;
            yy=ImageY-0.25;
         }
         else if(i2==1){
            xx=ImageX-0.25;
            yy=ImageY+0.25;
         }
         else if(i2==2){
            xx=ImageX+0.25;
            yy=ImageY+0.25;
         }
         else if(i2==3){
            xx=ImageX+0.25;
            yy=ImageY-0.25;
         }
         slaDtp2s(-xx/180*PI,yy/180*PI,raz,decz,phi0[i2],theta0[i2]);
         theta0[i2]=theta0[i2]/PI*180;
         phi0[i2]=phi0[i2]/PI*180;
         //printf("iSiPM=%d xx=%.1lf yy=%.1lf theta=%.1lf phi=%.1lf\n",ii,xx,yy,theta0[i2],phi0[i2]);
         if(phi0[i2]<xrange[0]) xrange[0]=phi0[i2];
         if(phi0[i2]>xrange[1]) xrange[1]=phi0[i2];
         if(theta0[i2]<yrange[0]) yrange[0]=theta0[i2];
         if(theta0[i2]>yrange[1]) yrange[1]=theta0[i2];
      }
      theta0[4]=theta0[0];
      phi0[4]=phi0[0];
      image->AddBin(5,phi0,theta0);
   }
   for(int ii=0;ii<iSiPM.size();ii++){
      double content=0;
      if(type==0) content=ADC_Cut.at(ii)>threshold?1.:0.;
      if(type==1) content=ImageAdcHigh.at(ii)/WFCTAMCEvent::fAmpHig;
      if(type==2) content=ImageAdcLow.at(ii)/WFCTAMCEvent::fAmpLow;
      if(type==3) content=myImageAdcHigh.at(ii)/WFCTAMCEvent::fAmpHig;
      if(type==4) content=myImageAdcLow.at(ii)/WFCTAMCEvent::fAmpLow;
      image->SetBinContent(iSiPM.at(ii)+1,content>0?content:0);
   }
   image->Draw(opt);
   image->GetXaxis()->SetRangeUser(xrange[0]-2.,xrange[1]+2.);
   image->GetYaxis()->SetRangeUser(yrange[0]-2.,yrange[1]+2.);
   return image;
}

TObjArray* WFCTAEvent::Draw3D(int type,const char* opt,double threshold,int ViewOpt){
   TObjArray* array=new TObjArray();
   double rmin[3]={1.0e10,1.0e10,1.0e100};
   double rmax[3]={-1.0e10,-1.0e10,-1.0e100};
   double tmin=mcevent.ArrivalTimeMin[0]*1.0e9; //in ns
   double tmax=mcevent.ArrivalTimeMax[0]*1.0e9; //in ns
   for(int ii=0;ii<NSIPM;ii++){
      double content=mcevent.TubeSignal[0][ii];
      if(content<=0) continue;
      if(tmax<=0) continue;
      int PixI=ii/PIX;
      int PixJ=ii%PIX;
      double ImageX,ImageY;
      if(PixI%2==0) ImageX=PixJ+0.5-PIX/2.0;
      else ImageX=PixJ+1.0-PIX/2.0;
      ImageY=(PIX/2.0-PixI)-1/2.0;

      ImageX=ImageX*16/32.0;
      ImageY=ImageY*16/32.0;
      ImageX-=0.31;
      ImageY-=0.28;

      if(ImageX<rmin[0]) rmin[0]=ImageX;
      if(ImageX>rmax[0]) rmax[0]=ImageX;
      if(ImageY<rmin[1]) rmin[1]=ImageY;
      if(ImageY>rmax[1]) rmax[1]=ImageY;
      if(tmin<rmin[2]) rmin[2]=tmin;
      if(tmax>rmax[2]) rmax[2]=tmax;
      //printf("Draw3D: ii=%d ImageX=%lf ImageY=%lf tmin=%le tmax=%le\n",ii,ImageX,ImageY,tmin,tmax);

      TMarker3DBox* box=new TMarker3DBox();
      box->SetDirection(0,0);
      box->SetPosition(ImageX,ImageY,(tmax+tmin)/2.);
      box->SetSize(0.5,0.5,(tmax-tmin)/2.);
      box->SetFillColor(2);
      array->Add((TObject*)box);
   }
   rmin[0]=-8.5; rmax[0]=8.5;
   rmin[1]=-8.5; rmax[1]=8.5;
   //printf("Draw3D range: xx={%lf,%lf},yy={%lf,%lf},zz={%le,%le}\n",rmin[0],rmax[0],rmin[1],rmax[1],rmin[2],rmax[2]);

   TCanvas* cc=new TCanvas();
   TView *view = TView::CreateView(1,rmin,rmax);
   view->ShowAxis();
   if(ViewOpt==1) view->Front();
   else if(ViewOpt==2) view->Side();
   else if(ViewOpt==3) view->Top();
   cc->SetView(view);

   array->Draw(opt);
   return array;
}

