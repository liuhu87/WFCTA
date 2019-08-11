#include <stdlib.h>
#include <iostream>
#include "WFCTAEvent.h"
#include "WFCamera.h"
#include "TH2Poly.h"

using namespace std;

ClassImp(WFCTAEvent);

WFCTAEvent* WFCTAEvent::_Head=0;
TTree* WFCTAEvent::_Tree=0;
const char* WFCTAEvent::_Name="Event";
TBranch* WFCTAEvent::bAll=0;
TBranch* WFCTAEvent::bEvent=0;
TBranch* WFCTAEvent::bTime=0;
TBranch* WFCTAEvent::btime=0;
TBranch* WFCTAEvent::bADC_Cut=0;
TBranch* WFCTAEvent::bBaseHigh=0;
TBranch* WFCTAEvent::bBaseLow=0;
TBranch* WFCTAEvent::bAdcHigh=0;
TBranch* WFCTAEvent::bAdcLow=0;
TBranch* WFCTAEvent::bmyBaseHigh=0;
TBranch* WFCTAEvent::bmyBaseLow=0;
TBranch* WFCTAEvent::bmyAdcHigh=0;
TBranch* WFCTAEvent::bmyAdcLow=0;
TBranch* WFCTAEvent::bSiPM=0;
TBranch* WFCTAEvent::bSThr=0;
TBranch* WFCTAEvent::bRThr=0;
TBranch* WFCTAEvent::bpeak=0;
TBranch* WFCTAEvent::bmypeak=0;
TBranch* WFCTAEvent::bgain_marker=0;
TBranch* WFCTAEvent::bOSMarker=0;
TBranch* WFCTAEvent::bORMarker=0;
TBranch* WFCTAEvent::bmcevent=0;
TBranch* WFCTAEvent::bledevent=0;
TBranch* WFCTAEvent::blaserevent=0;
WFCTAEvent::WFCTAEvent():TSelector()
{
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
   gain_marker.reserve(MAXPMT);
   Over_Single_Marker.reserve(MAXPMT);
   Over_Record_Marker.reserve(MAXPMT);

   Init();
}

WFCTAEvent::~WFCTAEvent()
{
   EventInitial();
}

void WFCTAEvent::Init()
{
   iEvent=-1;
   rabbitTime=0;
   rabbittime=0;
   iSiPM.clear();
   gain_marker.clear();
   peak.clear();
   mypeak.clear();
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
   
   mcevent.Init();
   ledevent.Init();
   laserevent.Init();
}
void WFCTAEvent::EventInitial()
{
   iEvent=-1;
   rabbitTime=0;
   rabbittime=0;
   iSiPM.clear();
   gain_marker.clear();
   peak.clear();
   mypeak.clear();
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
   bEvent=fChain->GetBranch("iEvent");
   bTime=fChain->GetBranch("rabbitTime");
   btime=fChain->GetBranch("rabbittime");
   bADC_Cut=fChain->GetBranch("ADC_Cut");
   bBaseHigh=fChain->GetBranch("ImageBaseHigh");
   bBaseLow=fChain->GetBranch("ImageBaseLow");
   bAdcHigh=fChain->GetBranch("ImageAdcHigh");
   bAdcLow=fChain->GetBranch("ImageAdcLow");
   bmyBaseHigh=fChain->GetBranch("myImageBaseHigh");
   bmyBaseLow=fChain->GetBranch("myImageBaseLow");
   bmyAdcHigh=fChain->GetBranch("myImageAdcHigh");
   bmyAdcLow=fChain->GetBranch("myImageAdcLow");
   bSiPM=fChain->GetBranch("iSiPM");
   bSThr=fChain->GetBranch("Single_Threshold");
   bRThr=fChain->GetBranch("Record_Threshold");
   bpeak=fChain->GetBranch("peak");
   bmypeak=fChain->GetBranch("mypeak");
   bgain_marker=fChain->GetBranch("gain_marker");
   bOSMarker=fChain->GetBranch("Over_Single_Marker");
   bORMarker=fChain->GetBranch("Over_Record_Marker");
   bmcevent=fChain->GetBranch("mcevent");
   bledevent=fChain->GetBranch("ledevent");
   blaserevent=fChain->GetBranch("laserevent");
}
bool WFCTAEvent::GetAllContents(int _Entry){
   EventInitial();
   bool exist=true;
   if(!bAll) exist=false;
   if(!bmcevent) exist=false;
   if(!bledevent) exist=false;
   if(!blaserevent) exist=false;
   if(!exist) return false;
   int ncount=bAll->GetEntry(_Entry);
   bEvent->GetEntry(_Entry);
   bTime->GetEntry(_Entry);
   btime->GetEntry(_Entry);
   bADC_Cut->GetEntry(_Entry);
   bBaseHigh->GetEntry(_Entry);
   bBaseLow->GetEntry(_Entry);
   bAdcHigh->GetEntry(_Entry);
   bAdcLow->GetEntry(_Entry);
   bmyBaseHigh->GetEntry(_Entry);
   bmyBaseLow->GetEntry(_Entry);
   bmyAdcHigh->GetEntry(_Entry);
   bmyAdcLow->GetEntry(_Entry);
   bSiPM->GetEntry(_Entry);
   bSThr->GetEntry(_Entry);
   bRThr->GetEntry(_Entry);
   bpeak->GetEntry(_Entry);
   bmypeak->GetEntry(_Entry);
   bgain_marker->GetEntry(_Entry);
   bOSMarker->GetEntry(_Entry);
   bORMarker->GetEntry(_Entry);
   bmcevent->GetEntry(_Entry);
   bledevent->GetEntry(_Entry);
   blaserevent->GetEntry(_Entry);
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

TH2Poly* WFCTAEvent::Draw(int type,const char* opt,double threshold){
   TH2Poly* image=new TH2Poly();
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
      if(type==0) {image->SetBinContent(iSiPM.at(ii)+1,ADC_Cut.at(ii)>threshold?1.:0.); /*printf("bin%d cont=%f\n",iSiPM.at(ii)+1,ADC_Cut.at(ii));*/}
      if(type==1) {image->SetBinContent(iSiPM.at(ii)+1,ImageAdcHigh.at(ii));}
      if(type==2) {image->SetBinContent(iSiPM.at(ii)+1,ImageAdcLow.at(ii));}
      if(type==3) {image->SetBinContent(iSiPM.at(ii)+1,myImageAdcHigh.at(ii));}
      if(type==4) {image->SetBinContent(iSiPM.at(ii)+1,myImageAdcLow.at(ii));}
   }
   image->Draw(opt);
   return image;
}
