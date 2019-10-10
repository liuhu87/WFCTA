#include <stdlib.h>
#include <cmath>
#include <iostream>
#include "WFCTAEvent.h"
#include "WFCamera.h"
#include "TH2Poly.h"
#include "TMarker3DBox.h"
#include "TAxis3D.h"
#include <TCanvas.h>
#include <TView3D.h>
#include <TSystem.h>

using namespace std;

ClassImp(WFCTAEvent);

WFCTAEvent* WFCTAEvent::_Head=0;
TTree* WFCTAEvent::_Tree=0;
const char* WFCTAEvent::_Name="WFCTAEvent";
TBranch* WFCTAEvent::bAll=0;
TBranch* WFCTAEvent::bmcevent=0;
TBranch* WFCTAEvent::bledevent=0;
TBranch* WFCTAEvent::blaserevent=0;
WFCTAEvent::WFCTAEvent():TSelector()
{
   eevent.reserve(MAXPMT);
   zipmod.reserve(MAXPMT);
   iSiPM.reserve(MAXPMT);
   winsum.reserve(MAXPMT);
   ADC_Cut.reserve(MAXPMT);
   eBaseH.reserve(MAXPMT);
   eBaseL.reserve(MAXPMT);
   eAdcH.reserve(MAXPMT);
   eAdcL.reserve(MAXPMT);
   eSatH.reserve(MAXPMT);
   eSatL.reserve(MAXPMT);
   BaseH.reserve(MAXPMT);
   BaseL.reserve(MAXPMT);
   AdcH.reserve(MAXPMT);
   AdcL.reserve(MAXPMT);
   SatH.reserve(MAXPMT);
   SatL.reserve(MAXPMT);
   Single_Threshold.reserve(MAXPMT);
   Record_Threshold.reserve(MAXPMT);
   peak.reserve(MAXPMT);
   PeakPosH.reserve(MAXPMT);
   PeakPosL.reserve(MAXPMT);
   PeakAmH.reserve(MAXPMT);
   PeakAmL.reserve(MAXPMT);
   gain_marker.reserve(MAXPMT);
   Over_Single_Marker.reserve(MAXPMT);
   Over_Record_Marker.reserve(MAXPMT);

   ADC_Cut.resize(MAXPMT);
   eevent.resize(MAXPMT);
   zipmod.resize(MAXPMT);
   iSiPM.resize(MAXPMT);
   winsum.resize(MAXPMT);
   eBaseH.resize(MAXPMT);
   eBaseL.resize(MAXPMT);
   eAdcH.resize(MAXPMT);
   eAdcL.resize(MAXPMT);
   eSatH.resize(MAXPMT);
   eSatL.resize(MAXPMT);
   BaseH.resize(MAXPMT);
   BaseL.resize(MAXPMT);
   AdcH.resize(MAXPMT);
   AdcL.resize(MAXPMT);
   SatH.resize(MAXPMT);
   SatL.resize(MAXPMT);
   Single_Threshold.resize(MAXPMT);
   Record_Threshold.resize(MAXPMT);
   peak.resize(MAXPMT);
   PeakPosH.resize(MAXPMT);
   PeakPosL.resize(MAXPMT);
   PeakAmH.resize(MAXPMT);
   PeakAmL.resize(MAXPMT);
   gain_marker.resize(MAXPMT);
   Over_Single_Marker.resize(MAXPMT);
   Over_Record_Marker.resize(MAXPMT);

   Init();
}

WFCTAEvent::~WFCTAEvent()
{
   EventInitial();
}

void WFCTAEvent::Init()
{
   iTel=-1;
   iEvent=-1;
   eEvent=-1;
   rabbitTime=0;
   rabbittime=0;
   big_pack_lenth=-1;
   n_fired=-1;
   n_Channel=-1;
   iSiPM.clear();
   eevent.clear();
   zipmod.clear();
   gain_marker.clear();
   peak.clear();
   PeakPosH.clear();
   PeakPosL.clear();
   PeakAmH.clear();
   PeakAmL.clear();
   Single_Threshold.clear();
   Record_Threshold.clear();
   Over_Single_Marker.clear();
   Over_Record_Marker.clear();
   winsum.clear();
   ADC_Cut.clear();
   eBaseH.clear();
   eBaseL.clear();
   eAdcH.clear();
   eAdcL.clear();
   eSatH.clear();
   eSatL.clear();
   BaseH.clear();
   BaseL.clear();
   AdcH.clear();
   AdcL.clear();
   SatH.clear();
   SatL.clear();

   for(int j=0;j<28;j++){
     Npoint[j]=j;
     for(int i=0;i<1024;i++){
       pulsehigh[i][j] = 0;
       pulselow[i][j] = 0;
     }
   }
  
   mcevent.Init();
   ledevent.Init();
   laserevent.Init();
}
void WFCTAEvent::EventInitial()
{
   iTel=-1;
   iEvent=-1;
   eEvent=-1;
   rabbitTime=0;
   rabbittime=0;
   big_pack_lenth=-1;
   n_fired=-1;
   n_Channel=-1;
   iSiPM.clear();
   zipmod.clear();
   eevent.clear();
   gain_marker.clear();
   peak.clear();
   PeakPosH.clear();
   PeakPosL.clear();
   PeakAmH.clear();
   PeakAmL.clear();
   Single_Threshold.clear();
   Record_Threshold.clear();
   Over_Single_Marker.clear();
   Over_Record_Marker.clear();
   winsum.clear();
   ADC_Cut.clear();
   eBaseH.clear();
   eBaseL.clear();
   eAdcH.clear();
   eAdcL.clear();
   eSatH.clear();
   eSatL.clear();
   BaseH.clear();
   BaseL.clear();
   AdcH.clear();
   AdcL.clear();
   SatH.clear();
   SatL.clear();

   for(int j=0;j<28;j++){
     for(int i=0;i<1024;i++){
       pulsehigh[i][j] = 0;
       pulselow[i][j] = 0;
     }
   }
  
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
         eAdcH.push_back(mcevent.TubeSignal[itel][ii]*WFCTAMCEvent::fAmpHig);
         eAdcL.push_back(mcevent.TubeSignal[itel][ii]*WFCTAMCEvent::fAmpLow);
         AdcH.push_back(mcevent.TubeSignal[itel][ii]*WFCTAMCEvent::fAmpHig);
         AdcL.push_back(mcevent.TubeSignal[itel][ii]*WFCTAMCEvent::fAmpLow);
      }
   }
}
int WFCTAEvent::GetMaxADCBin(){
   int res=-1;
   double maxadc=-1;
   for(int ii=0;ii<iSiPM.size();ii++){
      double yy=AdcH.at(ii);
      if(yy<=0) continue;
      if(yy>maxadc){
         maxadc=yy;
         res=ii;
      }
   }
   return res;
}
int WFCTAEvent::GetMaxTimeBin(){
   int res=-1;
   double maxtime=-1.e20;
   for(int ii=0;ii<1024;ii++){
      double yy=mcevent.ArrivalTimeMax[0][ii];
      if(yy<=0) continue;
      if(yy>maxtime){
         maxtime=yy;
         res=ii;
      }
   }
   return res;
}
int WFCTAEvent::GetMinTimeBin(){
   int res=-1;
   double mintime=1.e40;
   for(int ii=0;ii<1024;ii++){
      double yy=mcevent.ArrivalTimeMin[0][ii];
      if(yy<=0) continue;
      if(yy<mintime){
         mintime=yy;
         res=ii;
      }
   }
   return res;
}

TH2Poly* WFCTAEvent::Draw(int type,const char* opt,double threshold){
   TH2Poly* image=new TH2Poly();
   image->SetTitle(";X;Y");
   for(int ii=0;ii<NSIPM;ii++){
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

      image->AddBin(ImageX-0.25,ImageY-0.25,ImageX+0.25,ImageY+0.25);
   }
   for(int ii=0;ii<iSiPM.size();ii++){
      double content=0;
      if(type==0) content=ADC_Cut.at(ii)>threshold?1.:0.;
      if(type==1) content=eAdcH.at(ii)/WFCTAMCEvent::fAmpHig;
      if(type==2) content=eAdcL.at(ii)/WFCTAMCEvent::fAmpLow;
      if(type==3) content=AdcH.at(ii)/WFCTAMCEvent::fAmpHig;
      if(type==4) content=AdcL.at(ii)/WFCTAMCEvent::fAmpLow;
      image->SetBinContent(iSiPM.at(ii)+1,content>0?content:0);
   }
   image->Draw(opt);
   return image;
}

TObjArray* WFCTAEvent::Draw3D(int type,const char* opt,double threshold,int ViewOpt){
   TObjArray* array=new TObjArray();
   double rmin[3]={1.0e10,1.0e10,1.0e100};
   double rmax[3]={-1.0e10,-1.0e10,-1.0e100};
   for(int ii=0;ii<NSIPM;ii++){
      double content=mcevent.TubeSignal[0][ii];
      double tmin=mcevent.ArrivalTimeMin[0][ii]*1.0e9; //in ns
      double tmax=mcevent.ArrivalTimeMax[0][ii]*1.0e9; //in ns
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

/******************************************
 * read root file and do iamge processing *
 * ****************************************/
//set x and y, unit in degree
void WFCTAEvent::SetImage()
{
    int PixI,PixJ;
    for(int k=0;k<1024;k++)
    {
        PixI = (k) / PIX;
        PixJ = (k) % PIX;
        if(PixI%2==0)  {ImageX[k] = PixJ+0.5-PIX/2.0;}
        if(PixI%2==1)  {ImageX[k] = PixJ+1-PIX/2.0;}
        ImageY[k] = (PIX/2.0-PixI) - 1/2.0;

        ImageX[k] = ImageX[k]*16/32.;
        ImageY[k] = ImageY[k]*16/32.;

        ImageX[k] -= 0.31;
        ImageY[k] -= 0.28;
    }

}

//change rabbit time to local time
void WFCTAEvent::rabbittime2lt()
{
    double MJD19700101 = 40587;
    double TAI2UTC = 37;
    mjd = MJD19700101 + (rabbitTime + rabbittime*20/1000000000. - TAI2UTC)/86400;

    int j;
    double fd, d;
    long jd, n4, nd10;
    /* Check if date is acceptable */
    if ( ( mjd <= -2395520.0 ) || ( mjd >= 1e9 ) ) {
        j = -1;
    } else {
        j = 0;
    /* Separate day and fraction */
        fd = (mjd)>0.0?mjd-floor(mjd):mjd+floor(-mjd);
        if ( fd < 0.0 ) fd += 1.0;
        d = mjd - fd;
        d = d<0.0?ceil(d):floor(d);
    /* Express day in Gregorian calendar */
        jd = (long)d + 2400001;
        n4 = 4L*(jd+((6L*((4L*jd-17918L)/146097L))/4L+1L)/2L-37L);
        nd10 = 10L*(((n4-237L)%1461L)/4L)+5L;
        year = (int) (n4/1461L-4712L);
        month = (int) (((nd10/306L+2L)%12L)+1L);
        day = (int) ((nd10%306L)/10L+1L);
        j = 0;
        hour = int(fd * 24 + 8);
        minite = int((fd*24+8 - hour)*60);
        second = int(((fd*24+8-hour)*60-minite)*60);
		day += hour/24;
		hour = hour%24;
	//printf("time:%04d %02d%02d %02d:%02d:%02d\n",year,month,day,hour,minite,second);
    }

}

//initiate
void WFCTAEvent::InitImage()
{
	RawImagePe.clear();
	RawImageSiPM.clear();
	RawImageX.clear();
	RawImageY.clear();
	FullImagePe.clear();
	FullImageSiPM.clear();
	FullImageX.clear();
	FullImageY.clear();
	fNpixfriends.clear();
	CleanImagePe.clear();
	CleanImageSiPM.clear();
	CleanImageX.clear();
	CleanImageY.clear();
}

//change adc count to number of pe
void WFCTAEvent::AdcToPe()
{
   int isipm;
   double pe;
   double Ntotal = 360000;
   double theta;
   for(int ii=0;ii<iSiPM.size();ii++){
      isipm = (int)iSiPM.at(ii);
      if(SatH.at(ii)==0){  pe = AdcH.at(ii)/9.98;}
      else              {  pe = (AdcL.at(ii)*22)/9.98;}
      if(pe<0)          {  pe = 0;}
      else if(pe>Ntotal){  pe = Ntotal;}
      else              {  pe = -Ntotal*log(1-pe/Ntotal);}

      //pe = pe/factor[isipm];
      //pe = pe*(1 + deltag_20[isipm]*(correct_PreTemp[isipm]-T0));
      theta = pow(cos(sqrt(ImageX[isipm]*ImageX[isipm]+ImageY[isipm]*ImageY[isipm])/57.3),4);
	  //printf("%d theta:%lf\n",isipm,theta);
      pe = pe/theta;
      RawImagePe.push_back(pe);
      RawImageSiPM.push_back(isipm);
      RawImageX.push_back(ImageX[isipm]);
      RawImageY.push_back(ImageY[isipm]);
   }

}

//clean image preliminary
void WFCTAEvent::PrelimImageClean(double cut)
{
    for(int ii=0;ii<RawImagePe.size();ii++){
	if(RawImagePe.at(ii)<cut){continue;}
	FullImagePe.push_back(RawImagePe.at(ii));
        FullImageSiPM.push_back(RawImageSiPM.at(ii));
        FullImageX.push_back(RawImageX.at(ii));
        FullImageY.push_back(RawImageY.at(ii));
    }
}

//get neighbor trigger sipms of each sipm
void WFCTAEvent::GetNeighborPixs()
{
    int cnt;
    double distance;
    double MAXDIST=0.6;
    double x, y, x0, y0;
    for(int ii=0;ii<FullImagePe.size();ii++){
	cnt=0;
        x=FullImageX.at(ii);
        y=FullImageY.at(ii);
        for(int jj=0;jj<FullImagePe.size();jj++){
	    x0=FullImageX.at(jj);
            y0=FullImageY.at(jj);
            distance = sqrt((x0-x)*(x0-x)+(y0-y)*(y0-y));
            if(distance<=MAXDIST) {   //In degree
                cnt++;
            }
        }
        fNpixfriends.push_back(cnt);
    }
}

//calculate hillas parameters
int WFCTAEvent::CalcHillas()
{
    DNpix = 0;
    DSize = 0;
    DMeanX = 0;
    DMeanY = 0;
    Dslope = 0;
    Dintercept = 0;
    DLength = 0;
    DWidth = 0;
    double sx = 0;
    double sy = 0;
    double sxy = 0;
	for(int ii=0;ii<FullImagePe.size();ii++){
		if(fNpixfriends.at(ii)<=3){continue;}//use this to clean image
		DNpix++;
		DSize += FullImagePe.at(ii);
		DMeanX += FullImagePe.at(ii) * FullImageX.at(ii);
		DMeanY += FullImagePe.at(ii) * FullImageY.at(ii);
		sx += FullImagePe.at(ii) * FullImageX.at(ii) * FullImageX.at(ii);
		sy += FullImagePe.at(ii) * FullImageY.at(ii) * FullImageY.at(ii);
		sxy += FullImagePe.at(ii) * FullImageX.at(ii) * FullImageY.at(ii);
	}

    if(DSize==0.) return 2;
    DMeanX /= DSize;
    DMeanY /= DSize;
    sx /= DSize;
    sy /= DSize;
    sxy /= DSize;
    double cx = sx - DMeanX*DMeanX;
    double cy = sy - DMeanY*DMeanY;
    double cxy = sxy - DMeanX*DMeanY;
    double a = (cy-cx+sqrt((cy-cx)*(cy-cx)+4*cxy*cxy))/(2*cxy);
    double b = DMeanY-a*DMeanX;
    double ssx = (cx+2*a*cxy+a*a*cy)/(1+a*a);
    double ssy = (a*a*cx-2*a*cxy+cy)/(1+a*a);

    if(cx==0||cy==0) return 4;
    double delta = atan(a);
    Dslope = -a;
    Dintercept = -b;
    DLength = sqrt(ssx);
    DWidth = sqrt(ssy);

    return 0;
}

//clean image
void WFCTAEvent::GetCleanImage()
{
	for(int ii=0;ii<FullImagePe.size();ii++){
		if(fNpixfriends.at(ii)<=3){continue;}//use this to clean image
		CleanImagePe.push_back(FullImagePe.at(ii));
		CleanImageSiPM.push_back(FullImageSiPM.at(ii));
		CleanImageX.push_back(FullImageX.at(ii));
		CleanImageY.push_back(FullImageY.at(ii));
	}
}




