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
int WFCTAEvent::GetMaxADCBin(){
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
      if(type==1) content=ImageAdcHigh.at(ii)/WFCTAMCEvent::fAmpHig;
      if(type==2) content=ImageAdcLow.at(ii)/WFCTAMCEvent::fAmpLow;
      if(type==3) content=myImageAdcHigh.at(ii)/WFCTAMCEvent::fAmpHig;
      if(type==4) content=myImageAdcLow.at(ii)/WFCTAMCEvent::fAmpLow;
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
