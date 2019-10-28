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
#include "TF1.h"
#include "Laser.h"

using namespace std;

ClassImp(WFCTAEvent);

const int WFCTAEvent::cleanPix=3;
WFCTAEvent* WFCTAEvent::_Head=0;
TTree* WFCTAEvent::_Tree=0;
const char* WFCTAEvent::_Name="WFCTAEvent";
TBranch* WFCTAEvent::bAll=0;
TBranch* WFCTAEvent::bmcevent=0;
TBranch* WFCTAEvent::bledevent=0;
TBranch* WFCTAEvent::blaserevent=0;
int WFCTAEvent::jdebug=0;
bool WFCTAEvent::DoDraw=false;
double WFCTAEvent::adccuttrigger=350.;
double WFCTAEvent::npetrigger=45.;
int WFCTAEvent::nfiretrigger=3;
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
   WCamera::SetSiPMMAP();
}

WFCTAEvent::~WFCTAEvent()
{
   EventInitial();
   if(gDraw) delete gDraw;
   if(gDrawErr) delete gDrawErr;
   if(minimizer) delete minimizer;
   if(tlist){
      int size=tlist->Capacity();
      for(int ii=0;ii<size;ii++){
         delete tlist->At(ii);
      }
      for(int ii=0;ii<size;ii++){
         tlist->RemoveLast();
      }
      delete tlist;
   }
   graphlist.clear();
}

void WFCTAEvent::Init()
{
   tlist=new TList();

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

   gDraw=0;
   gDrawErr=0;
   minimizer=0;

   mcevent.Init();
   ledevent.Init();
   laserevent.Init();
}
void WFCTAEvent::EventInitial()
{
   int size=tlist->Capacity();
   for(int ii=0;ii<size;ii++){
      delete tlist->At(ii);
   }
   for(int ii=0;ii<size;ii++){
      tlist->RemoveLast();
   }
   for(int ii=0;ii<graphlist.size();ii++){
      delete (graphlist.at(ii));
   }
   graphlist.clear();

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

   if(gDraw) {delete gDraw; gDraw=0;}
   if(gDrawErr) {delete gDrawErr; gDrawErr=0;}
   if(minimizer) {delete minimizer; minimizer=0;}

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

bool WFCTAEvent::CheckLaser(){
   bool res=laserevent.Time>0;
   CommonTools::IsLaser=res;
   return res;
}
bool WFCTAEvent::CheckMC(){
   bool res=mcevent.Ngen>0;
   return res;
}

void WFCTAEvent::CalculateDataVar(int itel){
   if(itel<0||itel>=NCTMax) return;
   double tmin=mcevent.ArrivalTimeMin[itel];
   double tmax=mcevent.ArrivalTimeMax[itel];
   long int itmin=WCamera::FindBin(tmin);
   long int itmax=WCamera::FindBin(tmax);
   CommonTools::ResetHArrival();
   for(int ii=0;ii<NSIPM;ii++){
      if(mcevent.TubeSignal[itel][ii]>0){
         int peaktime=-1;
         double avetime=0;
         double nave=0;
         long int mintime=10000000000;
         long int maxtime=-1;
         double maxcount=0;
         for(int jj=0;jj<mcevent.NArrival[itel];jj++){
            if(mcevent.ArrivalCount[itel][ii][jj]>maxcount){
               maxcount=mcevent.ArrivalCount[itel][ii][jj];
               peaktime=jj;
            }
            if(mcevent.ArrivalCount[itel][ii][jj]>0){
               if(mcevent.ArrivalTime[itel][jj]<mintime) mintime=mcevent.ArrivalTime[itel][jj];
               if(mcevent.ArrivalTime[itel][jj]>maxtime) maxtime=mcevent.ArrivalTime[itel][jj];
               avetime+=mcevent.ArrivalTime[itel][jj]*mcevent.ArrivalCount[itel][ii][jj];
               nave+=mcevent.ArrivalCount[itel][ii][jj];
            }
         }
         if(nave>0) avetime/=nave;
         if(maxcount<1) continue;

         iSiPM.push_back(ii);
         eAdcH.push_back(mcevent.TubeSignal[itel][ii]*WFCTAMCEvent::fAmpHig);
         eAdcL.push_back(mcevent.TubeSignal[itel][ii]*WFCTAMCEvent::fAmpLow);
         AdcH.push_back(mcevent.TubeSignal[itel][ii]*WFCTAMCEvent::fAmpHig);
         AdcL.push_back(mcevent.TubeSignal[itel][ii]*WFCTAMCEvent::fAmpLow);

         //PeakPosH.push_back((mcevent.ArrivalTime[itel][peaktime]-itmin));
         PeakPosH.push_back(int(avetime)-itmin);
         PeakPosL.push_back(int(avetime)-itmin);
         PeakAmH.push_back(maxcount);
         PeakAmL.push_back(maxcount);

         //printf("WFCTAEvent::CalculateDataVar: iTel=%d PMT=%4d Signal=%5.0lf time={%5ld,%5ld,%5ld,%5ld} peakamp=%4.0lf tminmax={%ld,%ld} Overflow=%d\n",itel,ii,mcevent.TubeSignal[itel][ii],(mcevent.ArrivalTime[itel][peaktime]-itmin),(long int)(avetime)-itmin,(mintime-itmin),(maxtime-itmin),maxcount,itmin,itmax,mcevent.OverFlow[itel]);
         for(int jj=0;jj<MaxTimeBin;jj++){
            if(mcevent.ArrivalCount[itel][ii][jj]>0){
               //double timeref=mcevent.ArrivalTime[itel][jj]-mcevent.ArrivalTime[itel][peaktime];
               double timeref=mcevent.ArrivalTime[itel][jj]-int(avetime);
               CommonTools::HArrival[ii]->Fill(timeref,mcevent.ArrivalCount[itel][ii][jj]);
               //printf("itel=%d pmt=%d jj=%d time={%ld,%ld} maxcount=%lf\n",itel,ii,jj,mcevent.ArrivalTime[itel][jj],mcevent.ArrivalTime[itel][peaktime],maxcount);
            }
         }
      }
   }
}
int WFCTAEvent::GetMaxADCBin(int itel){
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

double WFCTAEvent::GetContent(int isipm,int itel,int type,bool IsIndex){
   itel=0;
   double content=0;
   int size=iSiPM.size();
   int start=IsIndex?isipm:0;
   int end=IsIndex?isipm:(size-1);
   if(IsIndex&&(isipm<0||isipm>=size)) return content;
   for(int ii=start;ii<=end;ii++){
      if((!IsIndex)&&(iSiPM.at(ii)!=isipm)) continue;
      if(type==0&&ii<ADC_Cut.size()) content=ADC_Cut.at(ii);
      if(type==1&&ii<eAdcH.size()) content=eAdcH.at(ii)/WFCTAMCEvent::fAmpHig;
      if(type==2&&ii<eAdcL.size()) content=eAdcL.at(ii)/WFCTAMCEvent::fAmpLow;
      if(type==3&&ii<AdcH.size()) content=AdcH.at(ii)/WFCTAMCEvent::fAmpHig;
      if(type==4&&ii<AdcL.size()) content=AdcL.at(ii)/WFCTAMCEvent::fAmpLow;
      if(type==5&&ii<PeakPosH.size()) content=PeakPosH.at(ii)+0.01;
      if(type==6&&ii<PeakPosL.size()) content=PeakPosL.at(ii)+0.01;
      if(type==7&&ii<PeakAmH.size()) content=PeakAmH.at(ii)+0.01;
      if(type==8&&ii<PeakAmL.size()) content=PeakAmL.at(ii)+0.01;
      if(type==9&&ii<SatH.size()) content=SatH.at(ii)?1.:-1.;
      if(type==10&&ii<SatL.size()) content=SatL.at(ii)?1.:-1.;
   }
   return content;
}
double WFCTAEvent::GetContentError(int isipm,int itel,int type,bool IsIndex){
   itel=0;
   bool ismc=CheckMC();
   double econtent=10.;
   int size=iSiPM.size();
   int start=IsIndex?isipm:0;
   int end=IsIndex?isipm:(size-1);
   if(IsIndex&&(isipm<0||isipm>=size)) return econtent;
   for(int ii=start;ii<=end;ii++){
      if((!IsIndex)&&(iSiPM.at(ii)!=isipm)) continue;
      if(ismc) econtent=mcevent.eTubeSignal[itel][iSiPM.at(ii)];
      if(type==5&&ii<PeakPosH.size()) econtent=1;
      if(type==6&&ii<PeakPosL.size()) econtent=1;
      if(type==9&&ii<SatH.size()) econtent=0;
      if(type==10&&ii<SatL.size()) econtent=0;
   }
   return econtent;
}
bool WFCTAEvent::CleanImage(int isipm,int itel,bool IsIndex){
   itel=0;
   if(itel<0||itel>=WFTelescopeArray::CTNumber) return false;
   int size=iSiPM.size();
   if(IsIndex&&(isipm<0||isipm>=size)) return false;
   if((!IsIndex)&&(isipm<0||isipm>=NSIPM)) return false;
   bool res=false;
   double trig0=adccuttrigger;
   double trig1=npetrigger;

   int start=IsIndex?isipm:0;
   int end=IsIndex?isipm:(size-1);
   for(int ii=start;ii<=end;ii++){
      int isipm0=iSiPM.at(ii);
      if((!IsIndex)&&(isipm0!=isipm)) continue;
      bool sat0=(GetContent(IsIndex?ii:isipm0,itel,9,IsIndex)>0.5);
      double content0=GetContent(IsIndex?ii:isipm0,itel,0,IsIndex);
      double content=sat0?GetContent(IsIndex?ii:isipm0,itel,4,IsIndex):GetContent(IsIndex?ii:isipm0,itel,3,IsIndex);
      //if(ADC_Cut.size()>ii) {if(content0<trig0) continue;}
      //else{if(content<trig1) continue;}
      if(content<trig1) continue;
      double ImageXi=0,ImageYi=0;
      ImageXi=WCamera::GetSiPMX(ii)/WFTelescope::FOCUS/PI*180;
      ImageYi=WCamera::GetSiPMY(ii)/WFTelescope::FOCUS/PI*180;
      if(iTel==5&&rabbitTime<1570680000) {ImageXi*=-1; ImageYi*=-1;}
      int nneigh=0;
      for(int jj=0;jj<size;jj++){
         if(jj==ii) continue;
         double ImageXj=0,ImageYj=0;
         ImageXj=WCamera::GetSiPMX(jj)/WFTelescope::FOCUS/PI*180;
         ImageYj=WCamera::GetSiPMY(jj)/WFTelescope::FOCUS/PI*180;
         if(iTel==5&&rabbitTime<1570680000) {ImageXj*=-1; ImageYj*=-1;}
         double dist=sqrt(pow(ImageXi-ImageXj,2)+pow(ImageYi-ImageYj,2));
         if(dist>0.6) continue;
         bool sat1=(GetContent(jj,itel,9,true)>0.5);
         double contentj0=GetContent(jj,itel,0,true);
         double contentj=sat1?GetContent(jj,itel,4,true):GetContent(jj,itel,3,true);
         //if(ADC_Cut.size()>jj) {if(contentj0<trig0) continue;}
         //else{if(contentj<trig1) continue;}
         if(contentj<trig1) continue;
         nneigh++;
      }
      if(nneigh<nfiretrigger) continue;
      res=true;
   }
   return res;
}
bool WFCTAEvent::PassClean(int itel,int nthreshold){
   int size=iSiPM.size();
   int ncontent=0;
   for(int ii=0;ii<size;ii++){
      if(CleanImage(ii,itel,true)) ncontent++;
   }
   if(ncontent<nthreshold) return false;
   else return true;
}

double WFCTAEvent::Interface(const double* par){
   int size=iSiPM.size();
   double chi2=0;
   int ndof=0;
   double sum=0;
   for(int ii=0;ii<size;ii++){
      int isipm=iSiPM.at(ii);
      if(!CleanImage(ii,(int)(par[0]+0.5),true)) continue;
      double ImageX=WCamera::GetSiPMX(isipm)/WFTelescope::FOCUS;///PI*180;
      double ImageY=WCamera::GetSiPMY(isipm)/WFTelescope::FOCUS;///PI*180;
      if(iTel==5&&rabbitTime<1570680000) {ImageX*=-1; ImageY*=-1;}
      double content=GetContent(ii,(int)(par[0]+0.5),(int)(par[1]+0.5),true);
      double err=0.25/180.*PI;
      chi2+=pow((ImageX*sin(par[3])-ImageY*cos(par[3])+par[2])/err,2)*content;
      ndof++;
      sum+=content;
   }
   ////for test
   //for(int ii=0;ii<9;ii++){
   //   chi2+=pow(((ii-4)*sin(par[3])-1.*cos(par[3])+par[2])/1.,2)*1;
   //   ndof++;
   //   sum+=1;
   //}
   if(sum>0){
      chi2/=(sum/ndof);
      //chi2/=sum;
   }
   return chi2;
}
bool WFCTAEvent::DoFit(int itel,int type,bool force){
   itel=0;
   if(itel<0||itel>=WFTelescopeArray::CTNumber) return false;
   if(minimizer&&(!force)) return true;
   int size=iSiPM.size();
   int nbin=0;
   double mx,my,sx,sy,sxy,nn;
   for(int ii=0;ii<size;ii++){
      int isipm0=iSiPM.at(ii);
      if(!CleanImage(ii,itel,true)) continue;
      double content=GetContent(ii,itel,type,true);
      double ImageX,ImageY;
      ImageX=WCamera::GetSiPMX(isipm0)/WFTelescope::FOCUS;///PI*180;
      ImageY=WCamera::GetSiPMY(isipm0)/WFTelescope::FOCUS;///PI*180;
      if(iTel==5&&rabbitTime<1570680000) {ImageX*=-1; ImageY*=-1;}

      mx += ImageX*content;
      my += ImageY*content;
      sx += ImageX*ImageX*content;
      sy += ImageY*ImageY*content;
      sxy += ImageX*ImageY*content;
      nn += content;

      nbin++;
   }

   if(nn>0&&nbin>5){
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
   minimizer->SetFixedVariable(0,"iTel",0);
   minimizer->SetFixedVariable(1,"Type",3);
   //minimizer->SetLimitedVariable(2,"bb",b,fabs(0.01*b),-1.e6,1.e6);
   //minimizer->SetLimitedVariable(3,"kk",a,fabs(0.01*a),-1.e6,1.e6);
   minimizer->SetLimitedVariable(2,"cc",b/sqrt(1+a*a),0.5/180.*PI,-25./180.*PI,25./180.*PI);
   minimizer->SetLimitedVariable(3,"phi",a>=0?atan(a):(PI+atan(a)),fabs(1./180.*PI),0,0.999*PI);
   //minimizer->SetLimitedVariable(2,"cc",1,0.01,-5,5);
   //minimizer->SetLimitedVariable(3,"phi",0,fabs(0.1/180.*PI),-PI/2,PI/2);
   //minimizer->SetFixedVariable(3,"phi",0.);
   minimizer->Minimize();
   minimizer->Hesse();
   //printf("LinearFit: kk=%lf bb=%lf cc={%lf,%lf} phi={%lf,%lf}\n",a,b,b/sqrt(1+a*a),minimizer->X()[2],(a>=0?atan(a):(PI+atan(a)))/PI*180,minimizer->X()[3]/PI*180.);
   //printf("LinearFit: phi={%lf,%lf} CC={%lf,%lf}\n",minimizer->X()[3],minimizer->Errors()[3],minimizer->X()[2],minimizer->Errors()[2]);
   return true;
}

bool WFCTAEvent::GetCrossCoor(double x,double y,double &x0,double &y0){
   if(!minimizer) return false;
   double CC=minimizer->X()[2]/PI*180.;
   double phi=minimizer->X()[3];
   x0=cos(phi)*(x*cos(phi)+y*sin(phi))-CC*sin(phi);
   y0=sin(phi)*(x*cos(phi)+y*sin(phi))+CC*cos(phi);
   return true;
}
TH1F* WFCTAEvent::GetLongDistribution(int itel,int type){
   itel=0;
   if(itel<0||itel>=WFTelescopeArray::CTNumber) return 0;
   if(!minimizer) return 0;
   double CC=minimizer->X()[2]/PI*180.;
   double phi=minimizer->X()[3];
   double margin=24.;
   const int nbin=24;
   TH1F* hh=new TH1F("longdis",";long dist [degree];Entries",nbin,-12.,12.);
   double nfill[nbin];
   for(int ii=0;ii<nbin;ii++) nfill[ii]=0;

   int size=iSiPM.size();
   for(int ii=0;ii<size;ii++){
      int isipm0=iSiPM.at(ii);
      if(!CleanImage(ii,itel,true)) continue;
      double ImageX,ImageY;
      ImageX=WCamera::GetSiPMX(isipm0)/WFTelescope::FOCUS/PI*180;
      ImageY=WCamera::GetSiPMY(isipm0)/WFTelescope::FOCUS/PI*180;
      if(iTel==5&&rabbitTime<1570680000) {ImageX*=-1; ImageY*=-1;}
      double dist=fabs(sin(phi)*ImageX-cos(phi)*ImageY+CC);
      if(dist>margin) continue;
      double x0,y0,xx,yy;
      bool res1=GetCrossCoor(0.,0.,x0,y0);
      bool res2=GetCrossCoor(ImageX,ImageY,xx,yy);
      if((!res1)||(!res2)) continue;
      double length=sqrt(pow(xx-x0,2)+pow(yy-y0,2));
      int sign;
      if(sin(phi)>sqrt(2.)/2) sign=(yy>y0)?1:(-1);
      else sign=(xx>x0)?1:(-1);
      int ibin=hh->GetXaxis()->FindBin(length*sign);
      hh->SetBinContent(ibin,hh->GetBinContent(ibin)+GetContent(ii,itel,type,true));
      hh->SetBinError(ibin,sqrt(pow(hh->GetBinError(ibin),2)+pow(GetContentError(ii,itel,type,true),2)));
      nfill[ibin-1]+=1;
      //printf("ii=%d ibin=%d length=%lf xy={%lf,%lf,%lf,%lf,%lf,%lf} sign=%d\n",ii,ibin,length,x0,y0,ImageX,ImageY,xx,yy,sign);
   }
   if(type==5){
      for(int ibin=1;ibin<=nbin;ibin++) {if(nfill[ibin-1]>0) hh->SetBinContent(ibin,hh->GetBinContent(ibin)/nfill[ibin-1]);}
   }
   if(hh->Integral()<=0){
      delete hh;
      return 0;
   }
   return hh;
}
TH1F* WFCTAEvent::GetShortDistribution(int itel,int type){
   itel=0;
   if(itel<0||itel>=WFTelescopeArray::CTNumber) return 0;
   if(!minimizer) return 0;
   double CC=minimizer->X()[2]/PI*180.;
   double phi=minimizer->X()[3];
   const int nbin=24;
   TH1F* hh=new TH1F("longdis",";long dist [degree];Entries",nbin,-6.,6.);
   double nfill[nbin];
   for(int ii=0;ii<nbin;ii++) nfill[ii]=0;

   int size=iSiPM.size();
   for(int ii=0;ii<size;ii++){
      int isipm0=iSiPM.at(ii);
      if(!CleanImage(ii,itel,true)) continue;
      double ImageX,ImageY;
      ImageX=WCamera::GetSiPMX(isipm0)/WFTelescope::FOCUS/PI*180;
      ImageY=WCamera::GetSiPMY(isipm0)/WFTelescope::FOCUS/PI*180;
      if(iTel==5&&rabbitTime<1570680000) {ImageX*=-1; ImageY*=-1;}
      double dist=fabs(sin(phi)*ImageX-cos(phi)*ImageY+CC);
      double xx,yy;
      bool res=GetCrossCoor(ImageX,ImageY,xx,yy);
      if(!res) continue;
      int sign;
      if(sin(phi)>sqrt(2.)/2) sign=(ImageX>xx)?1:(-1);
      else sign=(ImageY>yy)?1:(-1);
      int ibin=hh->GetXaxis()->FindBin(dist*sign);
      hh->SetBinContent(ibin,hh->GetBinContent(ibin)+GetContent(ii,itel,type,true));
      hh->SetBinError(ibin,sqrt(pow(hh->GetBinError(ibin),2)+pow(GetContentError(ii,itel,type,true),2)));
      nfill[ibin-1]+=1;
   }
   if(type==5){
      for(int ibin=1;ibin<=nbin;ibin++) {if(nfill[ibin-1]>0) hh->SetBinContent(ibin,hh->GetBinContent(ibin)/nfill[ibin-1]);}
   }
   if(hh->Integral()<=0){
      delete hh;
      return 0;
   }
   return hh;
}
int WFCTAEvent::GetSign(bool IsLaser,bool IsMC){
   int itel=0;
   //check wheather the time is increasing when x/y is increasing, and the sign of A parameter
   if(!minimizer) return 0;
   double CC=minimizer->X()[2]/PI*180.;
   double phi=minimizer->X()[3];

   const int nbin=24;
   double xminmax[2]={-12.,12.};
   double width=(xminmax[1]-xminmax[0])/nbin;
   int nfill[nbin];
   double avetime[nbin];
   double totadc[nbin];
   for(int ii=0;ii<nbin;ii++){
      nfill[ii]=0;
      avetime[ii]=0;
      totadc[ii]=0;
   }

   int size=iSiPM.size();
   int nn=0;
   for(int ii=0;ii<size;ii++){
      int isipm=iSiPM.at(ii);
      if(!CleanImage(ii,itel,true)) continue;
      double ImageX,ImageY;
      ImageX=WCamera::GetSiPMX(isipm)/WFTelescope::FOCUS/PI*180;
      ImageY=WCamera::GetSiPMY(isipm)/WFTelescope::FOCUS/PI*180;
      if(iTel==5&&rabbitTime<1570680000) {ImageX*=-1; ImageY*=-1;}
      double x0,y0,xx,yy;
      bool res1=GetCrossCoor(0.,0.,x0,y0);
      bool res2=GetCrossCoor(ImageX,ImageY,xx,yy);
      if((!res1)||(!res2)) continue;

      double length=sqrt(pow(xx-x0,2)+pow(yy-y0,2));
      int sign;
      if(sin(phi)>sqrt(2.)/2) sign=(yy>y0)?1:(-1);	//use the information along the y axis
      else sign=(xx>x0)?1:(-1);				//use the information along the x axis
      int ibin=(length*sign-xminmax[0])/width;
      if(ibin<0||ibin>=nbin) continue;

      double tt=GetContent(ii,itel,5,true);
      double content=GetContent(ii,itel,3,true);
      nfill[ibin-1]+=1;
      avetime[ibin-1]+=tt;
      totadc[ibin-1]+=content;
      nn++;
   }
   if(nn<=0) return 0;

   int peakbin=-1;
   double maxadc=-1;
   for(int ibin=0;ibin<nbin;ibin++){
      if(totadc[ibin]<=0) continue;
      if(nfill[ibin]>0) avetime[ibin]/=nfill[ibin];
      if(totadc[ibin]>maxadc){
         maxadc=totadc[ibin];
         peakbin=ibin;
      }
   }
   if(peakbin<0) return 0;
   int nfill1=0,nfill2=0;
   double pos1=0,time1=0,tot1=0;
   double pos2=0,time2=0,tot2=0;
   for(int ibin=0;ibin<nbin;ibin++){
      double pos=xminmax[0]+width*(ibin+0.5);
      if(ibin<peakbin){
         pos1+=pos*totadc[ibin];
         tot1+=totadc[ibin];
      }
      else if(ibin>peakbin){
         pos2+=pos*totadc[ibin];
         tot2+=totadc[ibin];
      }
   }
   if(tot1>0) pos1/=tot1;
   for(int ibin=0;ibin<(int)(pos1-xminmax[0]);ibin++){
      if(avetime[ibin]<=0) continue;
      nfill1++;
      time1+=avetime[ibin];
   }
   if(nfill1>0) time1/=nfill1;
   pos1=(xminmax[0]+width*(peakbin+0.5))-pos1;

   if(tot2>0) pos2/=tot2;
   for(int ibin=(int)(pos2-xminmax[0]);ibin<nbin;ibin++){
      if(avetime[ibin]<=0) continue;
      nfill2++;
      time2+=avetime[ibin];
   }
   if(nfill2>0) time2/=nfill2;
   pos2=pos2-(xminmax[0]+width*(peakbin+0.5));

   //printf("WFCTAEvent::GetSign: test: phi=%lf nfill={%d,%d} time={%lf,%lf} pos={%lf,%lf}\n",phi,nfill1,nfill2,time1,time2,pos1,pos2);
   if((nfill1<=0)||(nfill2<=0)) return 0;

   //the sign of A is determined by the time sequence,because:
   //dy_obs=1/A*sin(phi)*pow((y_obs*cos(theta)+sin(theta))/sin(PHI),2)*dPHI
   //dx_obs=1/A*cos(phi)*pow((y_obs*cos(theta)+sin(theta))/sin(PHI),2)*dPHI
   //and also it can be determined by the Image itself, because the image is not symmetric before and after Xmax
   int sign1,sign2;
   if(fabs(sin(phi))>sqrt(2.)/2){ //use the information along the y axis
      int dydt1=(time2>=time1)?1:(-1);
      int dydt2=(pos2>=pos1)?1:(-1);
      int dPHIdt=IsLaser?1:(-1);
      sign1=dydt1/dPHIdt/sin(phi);
      sign2=dydt2/dPHIdt/sin(phi);
   }
   else{ //use the information along the x axis
      int dxdt1=(time2>=time1)?1:(-1);
      int dxdt2=(pos2>=pos1)?1:(-1);
      int dPHIdt=IsLaser?1:(-1);
      sign1=dxdt1/dPHIdt/cos(phi);
      sign2=dxdt2/dPHIdt/cos(phi);
   }

   int res1,res2;
   double timeunit=IsMC?(CommonTools::timebinunit[IsLaser]):20.;
   if(fabs(time2-time1)>=2) res1=(sign1>0)?1:(-1);
   else res1=0;
   if(fabs(pos2/pos1-1)<0.1) res2=(sign2>0)?1:(-1);
   else res2=0;

   //if(res1*res2!=-4) printf("WFCTAEvent::GetSign: phi=%lf time={%lf,%lf} pos={%lf,%lf} res={%d,%d} LaserMC={%d,%d}\n",phi,time1,time2,pos1,pos2,res1,res2,IsLaser,IsMC);

   //use time information
   if(IsLaser){
      if(fabs(res1)<0.5) return res2;
      else return res1;
   }
   else{ //use image information
      if(fabs(res2)<0.5) return res1;
      else return res2;
   }
}
double WFCTAEvent::GetApar(double CC,double phi,double zenith,double azimuth){ //A is a very important parameter
   //A*A=(sin(phi)*sin(phi)+pow(CC*cos(theta)+sin(theta)*cos(phi),2))
   double phi0=azimuth;
   double theta=PI/2-zenith;
   double p1=CC*cos(theta)+sin(theta)*cos(phi);
   double p2=sin(phi);
   double Asq=(p1*p1+p2*p2);
   return (GetSign(CheckLaser(),CheckMC())>=0?1:-1)*sqrt(Asq);
}
double WFCTAEvent::GetApar(double &Anz,double zenith,double azimuth,double planephi,double nz){ //A or A*nz is a very important parameter
   double phi0=azimuth;
   double theta=PI/2-zenith;

   double AA=0;
   bool finitenz=isfinite(nz);

   double p1,p2,p3;
   if(finitenz){ //CC=AA*p1,cos(phi)=AA*p2,sin(phi)=AA*p3
      p1=sin(theta)*nz-cos(theta)*sin(phi0-planephi);
      p2=-cos(theta)*nz-sin(theta)*sin(phi0-planephi);
      p3=cos(phi0-planephi);
   }
   else{ //CC=Anz*p1,cos(phi)=Anz*p2,sin(phi)=Anz*p3;
      p1=sin(theta);
      p2=-cos(theta);
      p3=0;
   }

   double norm=sqrt(p2*p2+p3*p3);
   if(norm==0){ //CC is infinite
      Anz=0;
      return AA;
   }
   else{
      double phi_ref=acos(p3/norm);
      if(p2<0) phi_ref=2*PI-phi_ref;
      double phi1=CommonTools::ProcessAngle(-phi_ref+PI/2);
      double phi2=CommonTools::ProcessAngle(-phi_ref-PI/2);
      double phi;
      if(fabs(phi1-PI/2)==fabs(phi2-PI/2)) phi=phi1<phi2?phi1:phi2;
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
bool WFCTAEvent::CalPlane(double CC,double phi,double zenith,double azimuth,double &planephi,double &nz,int &signnz){
   double phi0=azimuth;
   double theta=PI/2-zenith;

   //when sin(phi)=0,CC=-sin(theta)/cos(theta)*cos(phi), then we can't measure the planephi in this situation, but we can get zdir[3]={0,0,1}
   double AA=GetApar(CC,phi,zenith,azimuth);
   int sign=GetSign(CheckLaser(),CheckMC())>=0?1:-1;
   if(AA==0){
      planephi=0; //any planephi is correct
      nz=sign>=0?1.0e8:-1.0e8;
      signnz=sign;
      return false;
   }
   else{
      double p1=(CC*cos(theta)+sin(theta)*cos(phi))/(AA);
      double p2=sin(phi)/AA;
      double phi_ref=acos(p1/fabs(AA));
      planephi=CommonTools::ProcessAngle(phi0+acos(p2));
      if(p1<0) planephi=CommonTools::ProcessAngle(phi0+2*PI-acos(p2));
      nz=(CC*sin(theta)-cos(theta)*cos(phi))/AA;
      int sign0=(CC*sin(theta)-cos(theta)*cos(phi))>=0?1:-1;
      signnz=sign*sign0;
      return true;
   }
}
bool WFCTAEvent::CalPlane(double CC,double phi,double zenith,double azimuth,double xyzdir[3][3]){
   double planephi,nz;
   int signnz;
   bool res=CalPlane(CC,phi,zenith,azimuth,planephi,nz,signnz);
   xyzdir[0][0]=cos(planephi);
   xyzdir[0][1]=sin(planephi);
   xyzdir[0][2]=0;
   xyzdir[2][0]=sin(planephi);
   xyzdir[2][1]=-cos(planephi);
   xyzdir[2][2]=nz;
   double norm=sqrt(pow(xyzdir[2][0],2)+pow(xyzdir[2][1],2)+pow(xyzdir[2][2],2));
   for(int icoo=0;icoo<3;icoo++) xyzdir[2][icoo]/=norm;
   Laser::cross(xyzdir[2],xyzdir[0],xyzdir[1]);
   return res;
}
bool WFCTAEvent::GetPlane(double &planephi,double &eplanephi,double &nz,double &enz,int itel,int type){
   WFTelescopeArray* pct=WFTelescopeArray::GetHead();
   if(!pct) return false;
   WFTelescope* pt=pct->pct[itel];
   if(!pt) return false;
   if(!DoFit(itel,type)) return false;
   double phi=minimizer->X()[3];
   double CC=minimizer->X()[2];

   int signnz;
   CalPlane(CC,phi,pt->TelZ_,pt->TelA_,planephi,nz,signnz);
   //error calculation
   double eele=1./180.*PI; //error of Tel zenith angle
   double eazi=1./180.*PI; //error of Tel azimuth angle
   double Cov[3]={minimizer->CovMatrix(2,2),minimizer->CovMatrix(2,3),minimizer->CovMatrix(3,3)};
   double coeff[2][4];
   eplanephi=0;
   enz=0;
   for(int iv=0;iv<4;iv++){
      coeff[0][iv]=0;
      coeff[1][iv]=0;
   }
   for(int iv=0;iv<4;iv++){
      double CC1=CC;
      double phi1=phi;
      double telz1=pt->TelZ_;
      double tela1=pt->TelA_;
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
         telz1=pt->TelZ_+eele;
         if(telz1>PI/2) telz1=pt->TelZ_-eele;
         delta=telz1-pt->TelZ_;
      }
      else if(iv==3){
         tela1=pt->TelA_+eazi;
         delta=tela1-pt->TelA_;
      }
      double planephi1,nz1;
      CalPlane(CC1,phi1,telz1,tela1,planephi1,nz1,signnz);
      eplanephi+=pow(planephi1-planephi,2);
      enz+=pow(nz1-nz,2);
      if(delta!=0){
         coeff[0][iv]=(planephi1-planephi)/delta;
         coeff[1][iv]=(nz1-nz)/delta;
      }
   }
   eplanephi+=2*coeff[0][0]*coeff[0][1]*Cov[1];
   eplanephi=sqrt(eplanephi);
   enz+=2*coeff[1][0]*coeff[1][1]*Cov[1];
   enz=sqrt(enz);
   return true;
}
bool WFCTAEvent::GetPlane(double xyzdir[3][3],double exyzdir[3][3],int itel,int type){
   WFTelescopeArray* pct=WFTelescopeArray::GetHead();
   if(!pct) return false;
   WFTelescope* pt=pct->pct[itel];
   if(!pt) return false;
   if(!DoFit(itel,type)) return false;
   double phi=minimizer->X()[3];
   double CC=minimizer->X()[2];

   CalPlane(CC,phi,pt->TelZ_,pt->TelA_,xyzdir);
   //error calculation
   double eele=1./180.*PI; //error of Tel zenith angle
   double eazi=1./180.*PI; //error of Tel azimuth angle
   double Cov[3]={minimizer->CovMatrix(2,2),minimizer->CovMatrix(2,3),minimizer->CovMatrix(3,3)};
   double coeff[3][3][4];
   for(int idir=0;idir<3;idir++){
      for(int icoo=0;icoo<3;icoo++){
         exyzdir[idir][icoo]=0;
         for(int iv=0;iv<4;iv++) coeff[idir][icoo][iv]=0;
      }
   }
   for(int iv=0;iv<4;iv++){
      double CC1=CC;
      double phi1=phi;
      double telz1=pt->TelZ_;
      double tela1=pt->TelA_;
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
         telz1=pt->TelZ_+eele;
         if(telz1>PI/2) telz1=pt->TelZ_-eele;
         delta=telz1-pt->TelZ_;
      }
      else if(iv==3){
         tela1=pt->TelA_+eazi;
         delta=tela1-pt->TelA_;
      }
      double xyzdir1[3][3];
      CalPlane(CC1,phi1,telz1,tela1,xyzdir1);
      for(int idir=0;idir<3;idir++){
         for(int icoo=0;icoo<3;icoo++){
            exyzdir[idir][icoo]+=pow(xyzdir1[idir][icoo]-xyzdir[idir][icoo],2);
            if(delta!=0) coeff[idir][icoo][iv]=(xyzdir1[idir][icoo]-xyzdir[idir][icoo])/delta;
         }
      }
   }
   for(int idir=0;idir<3;idir++){
      for(int icoo=0;icoo<3;icoo++){
         exyzdir[idir][icoo]+=2*coeff[idir][icoo][0]*coeff[idir][icoo][1]*Cov[1];
         exyzdir[idir][icoo]=sqrt(exyzdir[idir][icoo]);
      }
   }
   return true;
}

int WFCTAEvent::CalTelDir(double CC,double phi,double* lasercoo,double* laserdir,double &elevation,double &azimuth){
   //zdir should be (m,n,l)
   double dist=sqrt(pow(lasercoo[0],2)+pow(lasercoo[1],2));
   if(dist<=0) return -1;
   double planephi=acos(lasercoo[0]/dist);
   if(lasercoo[1]<0) planephi=2*PI-planephi;
   double zdir[3];
   zdir[0]=sin(planephi)*cos(laserdir[0]);
   zdir[1]=-cos(planephi)*cos(laserdir[0]);
   zdir[2]=sin(laserdir[0])*sin(laserdir[1]-planephi);

   double margin=1.0e-6;
   int signA=GetSign(CheckLaser(),CheckMC())>=0?1:-1;
   //calculate elevation
   //when phi=PI/2,CC=0(or theta=0,planephi=azimuth), then we can't measure the elevation of telescope in this situation, but we can measure the planephi quite precisely
   double nz=zdir[2]/cos(laserdir[0]);
   double Anz;
   if(cos(laserdir[0])==0){ //nz is infinite
      if(fabs(CC)<fabs(cos(phi))) Anz=(-cos(phi)>=0)?sqrt(1+CC*CC):(-sqrt(1+CC*CC));
      else Anz=(CC>=0)?sqrt(1+CC*CC):(-sqrt(1+CC*CC));
   }
   else{
      Anz=fabs(nz)>1/margin?sqrt(1./(1+1./nz/nz)):sqrt(nz*nz/(1+nz*nz));
      Anz*=sqrt(1+CC*CC);
      int signnz=nz>=0?1:-1;
      Anz*=(signA*signnz);
   }

   double p1=-cos(phi);
   double p2=CC;
   double norm=sqrt(p1*p1+p2*p2);
   if(norm==0){
      if(nz==0){
         printf("WFCTAEvent::CalTelDir: any theta is correct. Set theta to 90 degree\n");
         elevation=90;
         azimuth=CommonTools::ProcessAngle(planephi+signA>=0?0:PI)/PI*180;
         return 1;
      }
      else{
         printf("WFCTAEvent::CalTelDir: error of the input CC=%lf and phi=%lf when nz=%lf\n",CC/PI*180,phi/PI*180,nz);
         return -2;
      }
   }
   else if(fabs(norm)<fabs(Anz)){
      printf("WFCTAEvent::CalTelDir: no theta is correct (Anz=%lf norm=%lf). Exiting...\n",Anz,norm);
      return -3;
   }
   else{
      double theta_ref=acos(p1/norm);
      if(p2<0) theta_ref=2*PI-theta_ref;
      double sol1=CommonTools::ProcessAngle(theta_ref+acos(Anz/norm));
      double sol2=CommonTools::ProcessAngle(theta_ref+2*PI-acos(Anz/norm));
      if(fabs(sol1-PI/4)<fabs(sol2-PI/4)) elevation=sol1/PI*180;
      else elevation=sol2/PI*180;
   }

   //calculate azimuth
   //when m=0 and n=0, then we can't measure the azimuth angle of telescope
   double p11=-(CC*cos(elevation)+sin(elevation)*cos(phi));
   double p22=sin(phi);
   double AA=signA*sqrt(p11*p11+p22*p22);
   if(AA==0){
      printf("WFCTAEvent::CalTelDir: any azimuth is correct. Set azimuth to 0 degree\n");
      azimuth=0;
      return 1;
   }
   else{
      double azi_ref=acos(p22/AA);
      if(p11/AA<0) azi_ref=2*PI-azi_ref;
      azimuth=CommonTools::ProcessAngle(planephi-azi_ref)/PI*180;
   }

   return 1;
}
int WFCTAEvent::GetTelDir(double &elevation,double &errel,double &azimuth,double &erraz){
   if(!minimizer) return -1;
   double lasercoo[3]={laserevent.LaserCoo[0],laserevent.LaserCoo[1],laserevent.LaserCoo[2]};
   double xdir[3]={lasercoo[0]-0,lasercoo[1]-0,lasercoo[2]-0};
   double lasertheta=laserevent.LaserDir[0]/180.*PI;
   double laserphi=laserevent.LaserDir[1]/180.*PI;
   if(xdir[2]!=0){
      double dz=-xdir[2]/cos(lasertheta);
      xdir[0]+=(sin(lasertheta)*cos(laserphi))*dz;
      xdir[1]+=(sin(lasertheta)*sin(laserphi))*dz;
      xdir[2]=0;
   }
   double laserdir[2]={lasertheta,laserphi};
   double phi=minimizer->X()[3];
   double CC=minimizer->X()[2];
   int res=CalTelDir(CC,phi,xdir,laserdir,elevation,azimuth);
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
      int res1=CalTelDir(CC1,phi1,lasercoo1,laserdir1,elevation1,azimuth1);
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
}

void WFCTAEvent::Getnz(double incoo[3],double indir[2],double &planephi,double &nz,double &nz1,double &nz2){
   double xdir[3]={incoo[0]-0,incoo[1]-0,incoo[2]-0};
   double intheta=indir[0];
   double inphi=indir[1];
   if(incoo[2]!=0){
      double dz=-incoo[2]/cos(intheta);
      xdir[0]+=(sin(intheta)*cos(inphi))*dz;
      xdir[1]+=(sin(intheta)*sin(inphi))*dz;
      xdir[2]=0;
   }
   double dist=sqrt(xdir[0]*xdir[0]+xdir[1]*xdir[1]);
   if(dist<=0) planephi=0;
   else{
      planephi=acos(xdir[0]/dist);
      if(incoo[1]<0) planephi=2*PI-planephi;
   }
   double nz00=sin(indir[0])*sin(indir[1]-planephi);
   nz=nz00/cos(indir[0]);
   nz1=nz00/sqrt(cos(indir[0])*cos(indir[0])+nz00*nz00);
   nz2=cos(indir[0])/sqrt(cos(indir[0])*cos(indir[0])+nz00*nz00);
}
void WFCTAEvent::GetCCphi(double zenith,double azimuth,double planephi,double nz,double &CC,double &phi){
   double theta=PI/2-zenith;
   double phi0=azimuth;
   bool finitenz=isfinite(nz);
   double Anz;
   double AA=GetApar(Anz,zenith,azimuth,planephi,nz);
   CC=Anz*sin(theta)-AA*cos(theta)*sin(phi0-planephi);
   phi=acos(-Anz*cos(theta)-AA*sin(theta)*sin(phi0-planephi));
   if(AA*cos(phi0-planephi)<0) phi=2*PI-phi;
   //if(jdebug>0) printf("zen=%lf azi=%lf planephi=%lf nz=%lf AA=%lf Anz=%lf CC=%lf,phi=%lf\n",zenith/PI*180,azimuth/PI*180,planephi/PI*180,nz,AA,Anz,CC/PI*180,phi/PI*180);
}
void WFCTAEvent::GetCCphi(double zenith,double azimuth,double incoo[3],double indir[2],double &CC,double &phi){
   double planephi,nz,nz1,nz2;
   Getnz(incoo,indir,planephi,nz,nz1,nz2);
   return GetCCphi(zenith,azimuth,planephi,nz,CC,phi);
}

bool WFCTAEvent::GetImageXYCoo(double zenith,double azimuth,double planephi,double nz1,double nz2,double PHI_in,double &xx,double &yy){
   double theta=PI/2-zenith;
   double phi0=azimuth;
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
   double planephi,nz,nz1,nz2;
   Getnz(incoo,indir,planephi,nz,nz1,nz2);
   //printf("zenith=%lf azi=%lf PHI=%lf indir={%lf,%lf} planephi=%lf nz=%lf\n",zenith,azimuth,PHI_in,indir[0]/PI*180,indir[1]/PI*180,planephi/PI*180,nz);
   return GetImageXYCoo(zenith,azimuth,planephi,nz1,nz2,PHI_in,xx,yy);

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
void WFCTAEvent::GetPHI(double zenith,double azimuth,double planephi,double nz1,double nz2,double* ImageCoo,double &PHI_in){
   double theta=PI/2-zenith;
   double phi0=azimuth;

   double x2=ImageCoo[1];
   double y2=ImageCoo[0];
   double px1=(x2*cos(theta)+sin(theta))*cos(phi0-planephi);
   double px2=((x2*cos(theta)+sin(theta))*sin(phi0-planephi)*nz1+nz2*(x2*sin(theta)-cos(theta)));
   double py1=y2*cos(theta)*cos(phi0-planephi)-sin(phi0-planephi);
   double py2=(y2*cos(theta)*sin(phi0-planephi)+cos(phi0-planephi))*nz1+nz2*y2*sin(theta);

   double nz=nz1/nz2;
   double CC,phi;
   GetCCphi(zenith,azimuth,planephi,nz,CC,phi);
   double norm=0;
   if(fabs(sin(phi))<sqrt(2.)/2){ //use x to eval PHI
      norm=sqrt(px1*px1+px2*px2);
      if(norm<=0){
         printf("WFCTAEvent::GetPHI: no solution with x=%lf\n",x2);
         PHI_in=0;
         return;
      }
      else{
         double PhI_ref=acos(px1/norm);
         if(px2<0) PhI_ref=2*PI-PhI_ref;
         double PHI1=CommonTools::ProcessAngle(PhI_ref+PI/2);
         double PHI2=CommonTools::ProcessAngle(PhI_ref-PI/2);
         if(fabs(PHI1-PI/2)<=fabs(PHI2-PI/2)) PHI_in=PHI1;
         else PHI_in=PHI2;
         return;
      }
   }
   else{ //use y to eval PHI
      norm=sqrt(py1*py1+py2*py2);
      if(norm<=0){
         printf("WFCTAEvent::GetPHI: no solution with y=%lf\n",y2);
         PHI_in=0;
         return;
      }
      else{
         double PhI_ref=acos(py1/norm);
         if(py2<0) PhI_ref=2*PI-PhI_ref;
         double PHI1=CommonTools::ProcessAngle(PhI_ref+PI/2);
         double PHI2=CommonTools::ProcessAngle(PhI_ref-PI/2);
         if(fabs(PHI1-PI/2)<=fabs(PHI2-PI/2)) PHI_in=PHI1;
         else PHI_in=PHI2;
         return;
      }
   }
}
void WFCTAEvent::GetPHI(double zenith,double azimuth,double incoo[3],double indir[2],double* ImageCoo,double &PHI_in){
   double planephi,nz,nz1,nz2;
   Getnz(incoo,indir,planephi,nz,nz1,nz2);
   return GetPHI(zenith,azimuth,planephi,nz1,nz2,ImageCoo,PHI_in);
}
void WFCTAEvent::GetPHI(double zenith,double azimuth,double CC,double phi,double* ImageCoo,double &PHI_in){
   double planephi,nz;
   int signnz;
   CalPlane(CC,phi,zenith,azimuth,planephi,nz,signnz);
   double nz1,nz2;
   bool finitenz=isfinite(nz);
   if(finitenz){
      nz1=nz/sqrt(1+nz*nz);
      nz2=1./sqrt(1+nz*nz);
   }
   else{
      nz1=signnz>=0?1:-1;
      nz2=0;
   }
   return  GetPHI(zenith,azimuth,planephi,nz1,nz2,ImageCoo,PHI_in);

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
   double nz1,nz2;
   double CC,phi;
   double normdir=sqrt(dirin[0]*dirin[0]+dirin[1]*dirin[1]+dirin[2]*dirin[2]);
   double cosPHIL;
   bool UsePhi=(normdir<=0)&&(planephi==0);
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
         double nz=p3/Apre;
         nz1=nz/sqrt(1+nz*nz);
         nz2=1./sqrt(1+nz*nz);
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
      double nz,planephi0;
      Getnz(incoo,dir,planephi0,nz,nz1,nz2);
      /*//calculate phi and CC
      double p1=sin(dir[0])*cos(theta)*sin(dir[1]-planephi)+sin(theta)*cos(dir[0])*sin(phi0-planephi);
      double p2=cos(dir[0])*cos(phi0-planephi);
      double norm1=sqrt(p1*p1+p2*p2);
      if(norm1<margin){
         if(jdebug>0) printf("WFCTAEvent::GetXYRange: return=-3, error between the input parameters, because norm1={%lf,%lf,%lf}\n",p1,p2,norm1);
         return -3;
      }

      //calculate phi
      double phi_ref=acos(p2/norm1);
      if(p1<0) phi_ref=2*PI-phi_ref;
      phi=phi_ref+PI/2;
      if(phi>=(2*PI)) phi-=(2*PI);
      if(phi>=PI) phi-=PI;

      //calculate CC
      double absApar1=1/norm1;
      double Apar1;
      if(fabs(p2/norm1)<margin) Apar1=cos(phi)/(-p1);
      else Apar1=sin(phi)/p2;
      CC=Apar1*(sin(theta)*sin(dir[0])*sin(dir[1]-planephi)-cos(theta)*cos(dir[0])*sin(phi0-planephi));

      //calculate nz
      double nz00=sin(dir[0])*sin(dir[1]-planephi);
      double nz=nz00/cos(dir[0]);
      if(isfinite(nz)){
         nz1=nz/sqrt(1+nz*nz);
         nz2=1/sqrt(1+nz*nz);
      }
      else{
         nz1=nz00>=0?1:(-1);
         nz2=0;
      }*/
   }

   double PHI1,PHI2;

   //make sure the image is inside the field view of telescope
   double boundary=8.5/180.*PI;
   double xysol[4][2];
   int index1=-1,index2=-1;
   for(int iline=0;iline<4;iline++){
      xysol[iline][0]=xysol[iline][1]=1.0e5;
      double dist;
      if((iline/2)==0){
         xysol[iline][0]=((iline%2)==0?1:(-1))*boundary;
         if(fabs(cos(phi))<margin){
            if(cos(phi)*(xysol[iline][0]*sin(phi)+CC)>=0) xysol[iline][1]=1.0e5;
            else xysol[iline][1]=-1.0e5;
         }
         else{
            xysol[iline][1]=(xysol[iline][0]*sin(phi)+CC)/cos(phi);
         }
      }
      else{
         xysol[iline][1]=((iline%2)==0?1:(-1))*boundary;
         if(fabs(sin(phi))<margin){
            if(sin(phi)*(xysol[iline][1]*cos(phi)-CC)>=0) xysol[iline][0]=1.0e5;
            else xysol[iline][0]=-1.0e5;
         }
         else{
            xysol[iline][0]=(xysol[iline][1]*cos(phi)-CC)/sin(phi);
         }
      }
      if(fabs(xysol[iline][0])<=boundary+0.1/180.*PI && fabs(xysol[iline][1])<=boundary+0.1/180*PI){
         if(index1<0) index1=iline;
         else if(index2<0) index2=iline;
         else{
            if(jdebug>0) printf("WFCTAEvent::GetXYRange: return=-1, error in calculationg boundaries. line=%d calY=%d {%lf,%lf}\n",iline,(iline/2)==0,xysol[iline][0],xysol[iline][1]);
            //return -1;
         }
      }
   }
   if(index1<=0||index2<=0){
      if(jdebug>0) printf("WFCTAEvent::GetXYRange: return=-4, Image is not in the field of view of telescope CC=%lf phi=%lf\n",CC/PI*180,phi/PI*180);
      return -4;
   }
   double xyboun[2][2]={{xysol[index1][1],xysol[index1][0]},{xysol[index2][1],xysol[index2][0]}};

   //from xy coor to PHI
   for(int ip=0;ip<2;ip++){
      double ImageCoo[2]={xyboun[ip][1],xyboun[ip][0]};
      double PHI_in;
      GetPHI(zenith,azimuth,planephi,nz1,nz2,ImageCoo,PHI_in);
      if(ip==0) PHI1=PHI_in;
      else PHI2=PHI_in;
      /*double x0=xyboun[ip][0];
      double y0=xyboun[ip][1];
      double p1,p2;
      if(fabs(sin(phi))<sqrt(2.)/2){ //use x coordinate
         p1=(x0*cos(theta)+sin(theta))*cos(phi0-planephi);
         p2=(x0*cos(theta)+sin(theta))*sin(phi0-planephi)*nz1+(x0*sin(theta)-cos(theta))*nz2;
      }
      else{ //use y coordinate
         p1=y0*cos(theta)*cos(phi0-planephi)-sin(phi0-planephi);
         p2=(y0*cos(theta)*sin(phi0-planephi)+cos(phi0-planephi))*nz1+y0*sin(theta)*nz2;
      }
      double norm1=sqrt(p1*p1+p2*p2);
      if(fabs(norm1)<margin){
         if(jdebug>0) printf("WFCTAEvent::GetXYRange: return=-5, error between the input parameters, because of norm={%lf,%lf,%lf}\n",p1,p2,norm1);
         return -5;
      }
      else{
         double PHI_ref=acos(p1/norm1);
         if(p2<0) PHI_ref=2*PI-PHI_ref;
         double PHI00=PHI_ref-PI/2.;
         if(!(PHI00>=0&&PHI00<PI)) PHI00=PHI_ref+PI/2.;
         if(PHI00>=2*PI) PHI00-=2*PI;
         if(ip==0) PHI1=PHI00;
         else PHI2=PHI00;
      }*/
   }
   if(PHI1>PHI2){
      double buff=PHI1;
      PHI1=PHI2;
      PHI2=buff;
   }
   PHI1=TMath::Max(0.,PHI1);
   PHI2=TMath::Min(UsePhi?PI:acos(cosPHIL),PHI2);
   if(PHI2<=PHI1){
      if(jdebug>0) printf("WFCTAEvent::GetXYRange: return=-6, No Image from 0 to acos(PHIL)=%lf\n",acos(cosPHIL)/PI*180);
      return -6;
   }

   //the PHI range to have image ((x*cos(phi0)+y*sin(phi0))*cos(theta)+z*sin(theta)>=0)
   double p21=cos(theta)*cos(phi0-planephi);
   double p22=nz1*sin(phi0-planephi)*cos(theta)+nz2*sin(theta);
   double norm2=sqrt(p21*p21+p22*p22);
   if(norm2<margin){ //any PHI is correct
      ;
   }
   else{
      double PHI_ref=acos(p21/norm2);
      if(p22<0) PHI_ref=2*PI-PHI_ref;
      double angle1=CommonTools::ProcessAngle(PHI1-PHI_ref);
      double angle2=CommonTools::ProcessAngle(PHI2-PHI_ref);
      bool in1=( (angle1>=0&&angle1<=PI/2) || (angle1>=1.5*PI&&angle1<=2*PI) );
      bool in2=( (angle2>=0&&angle2<=PI/2) || (angle2>=1.5*PI&&angle2<=2*PI) );
      if(in1&&in2){
         ;
      }
      else if(in1&&(!in2)){
         PHI2=CommonTools::ProcessAngle(PI/2+PHI_ref);
      }
      else if((!in1)&&in2){
         PHI1=CommonTools::ProcessAngle(-PI/2+PHI_ref);
      }
      else{
         if(jdebug>0) printf("WFCTAEvent::GetXYRange: return=-7, error between the input parameters, because of PHI+PHI0={%lf,%lf}\n",angle1/PI*180,angle2/PI*180);
         return -7;
      }
   }

   double PHI_low=IsLaser?PHI1:PHI2;
   double PHI_hig=IsLaser?PHI2:PHI1;
   PHI[0]=PHI_low;
   PHI[1]=PHI_hig;

   //from PHI to coor
   bool res1=GetImageXYCoo(zenith,azimuth,planephi,nz1,nz2,PHI[0],XX[0],YY[0]);
   bool res2=GetImageXYCoo(zenith,azimuth,planephi,nz1,nz2,PHI[1],XX[1],YY[1]);

   /*double coox[2]={cos(PHI_low)*cos(planephi)-nz1*sin(PHI_low)*sin(planephi),cos(PHI_hig)*cos(planephi)-nz1*sin(PHI_hig)*sin(planephi)};
   double cooy[2]={cos(PHI_low)*sin(planephi)+nz1*sin(PHI_low)*cos(planephi),cos(PHI_hig)*sin(planephi)+nz1*sin(PHI_hig)*cos(planephi)};
   double cooz[2]={nz2*sin(PHI_low),nz2*sin(PHI_hig)};
   //interchange between XX and YY
   int index=0;
   YY[index]=-( (coox[index]*cos(phi0)+cooy[index]*sin(phi0))*sin(theta)-cooz[index]*cos(theta) )/( (coox[index]*cos(phi0)+cooy[index]*sin(phi0))*cos(theta)+cooz[index]*sin(theta) );
   index=1;
   YY[index]=-( (coox[index]*cos(phi0)+cooy[index]*sin(phi0))*sin(theta)-cooz[index]*cos(theta) )/( (coox[index]*cos(phi0)+cooy[index]*sin(phi0))*cos(theta)+cooz[index]*sin(theta) );
   index=0;
   XX[index]=-( -coox[index]*sin(phi0)+cooy[index]*cos(phi0) )/( (coox[index]*cos(phi0)+cooy[index]*sin(phi0))*cos(theta)+cooz[index]*sin(theta) );
   index=1;
   XX[index]=-( -coox[index]*sin(phi0)+cooy[index]*cos(phi0) )/( (coox[index]*cos(phi0)+cooy[index]*sin(phi0))*cos(theta)+cooz[index]*sin(theta) );*/
   return 1;
}

TGraph* WFCTAEvent::DrawCorePos(double* corepos,int itel,int type){
   double planephi,eplanephi,enz,nz;
   if(!GetPlane(planephi,eplanephi,nz,enz,itel,type)) return 0;
   TGraph* gr=new TGraph();
   gr->SetTitle(";X [cm];Y [cm]");
   double length=sqrt(corepos[0]*corepos[0]+corepos[1]*corepos[1])*1.5;
   int np=100;
   for(int ii=0;ii<np;ii++){
      double ll=exp(log(0.5)+(log(length/0.5))/np*ii);
      double phi=planephi;
      gr->SetPoint(ii,ll*cos(phi),ll*sin(phi));
   }
   gr->SetLineColor(4);
   gr->SetLineWidth(3);
   gr->SetTitle(";X [cm];Y [cm]");
   graphlist.push_back(gr);
   if(DoDraw) gr->Draw("l");
   return gr;
}
TGraph* WFCTAEvent::DrawCoreReg(double* corepos,int itel,int type){
   double planephi,eplanephi,enz,nz;
   if(!GetPlane(planephi,eplanephi,nz,enz,itel,type)) return 0;
   TGraph* gr=new TGraph();
   gr->SetTitle(";X [cm];Y [cm]");
   double length=sqrt(corepos[0]*corepos[0]+corepos[1]*corepos[1])*1.5;
   int np=100;
   for(int ii=0;ii<np;ii++){
      double ll=exp(log(0.5)+(log(length/0.5))/np*ii);
      double phi=planephi;
      double ephi=eplanephi;
      gr->SetPoint(ii,ll*cos(phi-ephi),ll*sin(phi-ephi));
      gr->SetPoint(2*np-1-ii,ll*cos(phi+ephi),ll*sin(phi+ephi));
   }
   gr->SetFillStyle(3001);
   gr->SetFillColor(4);
   gr->SetTitle(";X [cm];Y [cm]");
   graphlist.push_back(gr);
   if(DoDraw) gr->Draw("F");
   return gr;
}
TGraph* WFCTAEvent::DrawImageLine(double zenith,double azimuth,double incoo[3],double indir[2]){
   double xdir[3]={incoo[0]-0,incoo[1]-0,incoo[2]-0};
   double intheta=indir[0];
   double inphi=indir[1];
   if(incoo[2]!=0){
      double dz=-incoo[2]/cos(intheta);
      xdir[0]+=(sin(intheta)*cos(inphi))*dz;
      xdir[1]+=(sin(intheta)*sin(inphi))*dz;
      xdir[2]=0;
   }
   double PHI_max=acos( (xdir[0]*sin(intheta)*cos(inphi)+xdir[1]*sin(intheta)*sin(inphi))/sqrt(pow(xdir[0],2)+pow(xdir[1],2)+pow(xdir[2],2)) );

   double CC,phi;
   GetCCphi(zenith,azimuth,incoo,indir,CC,phi);
   printf("CC=%lf phi=%lf\n",CC/PI*180,phi/PI*180);

   int np=100;
   TGraph* gr=new TGraph();
   for(int ii=0;ii<np;ii++){
      double PHI=0.+(PHI_max-0.)/np*ii;
      double xx,yy;
      bool res=GetImageXYCoo(zenith,azimuth,incoo,indir,PHI,xx,yy);
      if(!res) continue;
      gr->SetPoint(gr->GetN(),xx,yy);
   }
   gr->SetLineColor(1);
   gr->SetLineWidth(4);
   gr->Draw("l");
   //graphlist.push_back(gr);
   return gr;
}
TGraph* WFCTAEvent::DrawImageLine(int itel){
   itel=0;
   WFTelescopeArray* pct=WFTelescopeArray::GetHead();
   if(!pct) return 0;
   WFTelescope* pt=pct->pct[itel];
   if(!pt) return 0;
   double zenith=pt->TelZ_;
   double azimuth=pt->TelA_;

   double incoo[3]={laserevent.LaserCoo[0],laserevent.LaserCoo[1],laserevent.LaserCoo[2]};
   double intheta=laserevent.LaserDir[0]/180.*PI;
   double inphi=laserevent.LaserDir[1]/180.*PI;
   double indir[2]={intheta,inphi};
   if(jdebug>-1) printf("WFCTAEvent::DrawImageLine: zenith=%lf azimuth=%lf incoo={%lf,%lf,%lf} indir={%lf,%lf}\n",zenith/PI*180,azimuth/PI*180,incoo[0],incoo[1],incoo[2],indir[0]/PI*180,indir[1]/PI*180);
   return DrawImageLine(zenith,azimuth,incoo,indir);
}

void WFCTAEvent::DrawFit(){
   if(!minimizer) return;
   double phi=minimizer->X()[3];
   double CC=minimizer->X()[2];
   double ephi=minimizer->Errors()[3];
   double eCC=minimizer->Errors()[2];
   double Cov[3]={minimizer->CovMatrix(2,2),minimizer->CovMatrix(2,3),minimizer->CovMatrix(3,3)};
   //printf("WFCTAEvent::DrawFit: Fitting Error Information, fitpars={%lf,%lf} fitparse={%lf,%lf} fitcov={%lf,%lf,%lf}\n",CC,phi,eCC,ephi,sqrt(Cov[0]),Cov[1],sqrt(Cov[2]));
   if(gDraw) delete gDraw;
   gDraw=new TGraph();
   if(gDrawErr) delete gDrawErr;
   gDrawErr=new TGraph();
   int np=100;
   for(int ii=0;ii<np;ii++){
      double xx,yy,err;
      if(fabs(cos(phi))<sqrt(2.)/2){
         yy=(-8.5+(8.5*2)/np*ii);
         xx=(yy*cos(phi)-CC/PI*180.)/sin(phi);
         double ceff1=-180./PI/sin(phi);
         double ceff2=(CC/PI*180.*cos(phi)-yy)/pow(sin(phi),2);
         err=sqrt(pow(ceff1,2)*Cov[0]+2*ceff1*ceff2*Cov[1]+pow(ceff2,2)*Cov[2]);
         gDrawErr->SetPoint(ii,xx-err,yy);
         gDrawErr->SetPoint(2*np-1-ii,xx+err,yy);
      }
      else{
         xx=(-8.5+(8.5*2)/np*ii);
         yy=(xx*sin(phi)+CC/PI*180.)/cos(phi);
         double ceff1=180./PI/cos(phi);
         double ceff2=(CC/PI*180.*sin(phi)+xx)/pow(cos(phi),2);
         err=sqrt(pow(ceff1,2)*Cov[0]+2*ceff1*ceff2*Cov[1]+pow(ceff2,2)*Cov[2]);
         gDrawErr->SetPoint(ii,xx,yy-err);
         gDrawErr->SetPoint(2*np-1-ii,xx,yy+err);
      }
      gDraw->SetPoint(ii,xx,yy);
   }
   gDrawErr->SetFillStyle(3001);
   gDrawErr->SetFillColor(1);
   if(DoDraw) gDrawErr->Draw("F");
   gDraw->SetLineColor(1);
   gDraw->SetLineWidth(4);
   if(DoDraw) gDraw->Draw("l");
   printf("WFCTAEvent::DrawFit: phi={%lf,%lf} cc={%lf,%lf}\n",phi/PI*180,sqrt(Cov[2])/PI*180,CC/PI*180,eCC/PI*180);
}
TH2Poly* WFCTAEvent::Draw(int type,const char* opt,bool DoClean,double threshold){
   TH2Poly* image=new TH2Poly();
   image->SetName("DrawPlot");
   image->SetTitle(Form("iTel=%d iEvent=%d time=%ld+%lf(%d-%02d-%02d %02d:%02d:%02d);X [degree];Y [degree]",iTel,iEvent,rabbitTime,rabbittime*20*1.0e-9,CommonTools::TimeFlag((int)rabbitTime,1),CommonTools::TimeFlag((int)rabbitTime,2),CommonTools::TimeFlag((int)rabbitTime,3),CommonTools::TimeFlag((int)rabbitTime,4),CommonTools::TimeFlag((int)rabbitTime,5),CommonTools::TimeFlag((int)rabbitTime,6)));
   for(int ii=0;ii<NSIPM;ii++){
      double ImageX,ImageY;
      ImageX=WCamera::GetSiPMX(ii)/WFTelescope::FOCUS/PI*180;
      ImageY=WCamera::GetSiPMY(ii)/WFTelescope::FOCUS/PI*180;
      //if(iTel==5&&rabbitTime<1570680000) {ImageX*=-1; ImageY*=-1;}
      image->AddBin(ImageX-0.25,ImageY-0.25,ImageX+0.25,ImageY+0.25);
      //printf("WFCTAEvent::Draw: SiPM=%d ImageX=%lf ImageY=%lf\n",ii,ImageX,ImageY);
   }
   int ncontent=0;
   for(int ii=0;ii<iSiPM.size();ii++){
      if(DoClean) {if(!CleanImage(ii,iTel,true)) continue;}
      double content=GetContent(ii,iTel,type,true);
      image->SetBinContent(iSiPM.at(ii)+1,content>0?content:0);
      //printf("WFCTAEvent::Draw: SiPM=%d content=%lf\n",iSiPM.at(ii),AdcH.at(ii));
      ncontent++;
   }
   if(type==-1){
      for(int ii=0;ii<NSIPM;ii++){
         double content=((ii%2)==1)?ii:0;
         image->SetBinContent(ii+1,content>0?content:0);
      }
   }
   else{
      if(ncontent<5){
         delete image;
         return 0;
      }
   }
   image->Draw(opt);
   DrawFit();
   //TGraph* gr=DrawImageLine(0);
   //graphlist.push_back(gr);
   tlist->Add(image);
   return image;
}
TGraph* WFCTAEvent::DrawPulse(int isipm,const char* opt,bool IsHigh,bool DoClean){
   if(isipm<0||isipm>=NSIPM) return 0;
   if(DoClean&&(!PassClean(iTel))) return 0;
   TGraph* gr=new TGraph();
   for(int ii=0;ii<28;ii++){
      double xx=Npoint[ii];
      double yy=IsHigh?pulsehigh[isipm][ii]:pulselow[isipm][ii];
      if(yy>0) gr->SetPoint(gr->GetN(),xx,yy);
   }
   if(gr->GetN()<3){
      delete gr;
      return 0;
   }
   gr->SetTitle(Form("iTel=%d iEvent=%d iSiPM=%d time=%ld+%lf;Time Window Index;Amplitude [ADC]",iTel,iEvent,isipm,rabbitTime,rabbittime*20*1.0e-9));
   gr->SetMarkerStyle(20);
   gr->SetMarkerColor(4);
   gr->SetMarkerSize(1.3);
   gr->SetLineColor(4);
   gr->SetLineWidth(2);
   graphlist.push_back(gr);
   if(DoDraw) gr->Draw(opt);
   return gr;
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
TH2Poly* WFCTAEvent::DrawGlobal(int type,const char* opt,bool DoClean,double threshold){
   WFTelescopeArray* pct=WFTelescopeArray::GetHead();
   if(!pct) return 0;
   int itel=pct->GetTelescope(iTel);
   if(itel<0) return 0;
   WFTelescope* pt=pct->pct[itel];
   if(!pt) return 0;
   double raz=pt->TelA_;
   double decz=PI/2-pt->TelZ_;
   TH2Poly* image=new TH2Poly();
   image->SetName("DrawGlobalPlot");
   image->SetTitle(";Azimuth [degree];Elevation [degree]");
   double xrange[2]={1.0e5,-1.0e5};
   double yrange[2]={1.0e5,-1.0e5};
   for(int ii=0;ii<NSIPM;ii++){
      double ImageX,ImageY;
      ImageX=WCamera::GetSiPMX(ii)/WFTelescope::FOCUS/PI*180;
      ImageY=WCamera::GetSiPMY(ii)/WFTelescope::FOCUS/PI*180;
      if(iTel==5&&rabbitTime<1570680000) {ImageX*=-1; ImageY*=-1;}

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
      if(DoClean) {if(!CleanImage(ii,iTel,true)) continue;}
      double content=GetContent(ii,iTel,type,true);
      if(type==0&&ii<ADC_Cut.size()) content=content>adccuttrigger?1.:0.;
      image->SetBinContent(iSiPM.at(ii)+1,content>0?content:0);
   }
   if(DoDraw) image->Draw(opt);
   image->GetXaxis()->SetRangeUser(xrange[0]-2.,xrange[1]+2.);
   image->GetYaxis()->SetRangeUser(yrange[0]-2.,yrange[1]+2.);
   tlist->Add(image);
   return image;
}

TH2Poly* WFCTAEvent::DrawCloudFormat(int type,const char* opt,bool DoClean,double threshold){
   WFTelescopeArray* pct=WFTelescopeArray::GetHead();
   if(!pct) return 0;
   int itel=pct->GetTelescope(iTel);
   if(itel<0) return 0;
   WFTelescope* pt=pct->pct[itel];
   if(!pt) return 0;
   double raz=pt->TelA_;
   double decz=PI/2-pt->TelZ_;
   TH2Poly* image=new TH2Poly();
   image->SetName("DrawCloudFormatPlot");
   image->SetTitle(Form("iTel=%d iEvent=%d time=%ld+%lf(%d-%02d-%02d %02d:%02d:%02d);West<-->East [degree];South<-->North [degree]",iTel,iEvent,rabbitTime,rabbittime*20*1.0e-9,CommonTools::TimeFlag((int)rabbitTime,1),CommonTools::TimeFlag((int)rabbitTime,2),CommonTools::TimeFlag((int)rabbitTime,3),CommonTools::TimeFlag((int)rabbitTime,4),CommonTools::TimeFlag((int)rabbitTime,5),CommonTools::TimeFlag((int)rabbitTime,6)));
   double xrange[2]={1.0e5,-1.0e5};
   double yrange[2]={1.0e5,-1.0e5};
   for(int ii=0;ii<NSIPM;ii++){
      double ImageX,ImageY;
      ImageX=WCamera::GetSiPMX(ii)/WFTelescope::FOCUS/PI*180;
      ImageY=WCamera::GetSiPMY(ii)/WFTelescope::FOCUS/PI*180;
      if(iTel==5&&rabbitTime<1570680000) {ImageX*=-1; ImageY*=-1;}

      double x0[5],y0[5];
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
         double theta0,phi0;
         slaDtp2s(-xx/180*PI,yy/180*PI,raz,decz,phi0,theta0);
         theta0=PI/2-theta0;
         x0[i2]=-theta0/PI*180*sin(phi0);
         y0[i2]=theta0/PI*180*cos(phi0);
         //printf("iSiPM=%d xx=%.1lf yy=%.1lf theta=%.1lf phi=%.1lf xy={%.1lf,%.1lf}\n",ii,xx,yy,theta0,phi0,x0[i2],y0[i2]);
         if(x0[i2]<xrange[0]) xrange[0]=x0[i2];
         if(x0[i2]>xrange[1]) xrange[1]=x0[i2];
         if(y0[i2]<yrange[0]) yrange[0]=y0[i2];
         if(y0[i2]>yrange[1]) yrange[1]=y0[i2];
      }
      x0[4]=x0[0];
      y0[4]=y0[0];
      image->AddBin(5,x0,y0);
   }
   for(int ii=0;ii<iSiPM.size();ii++){
      if(DoClean) {if(!CleanImage(ii,iTel,true)) continue;}
      double content=GetContent(ii,iTel,type,true);
      if(type==0&&ii<ADC_Cut.size()) content=content>adccuttrigger?1.:0.;
      image->SetBinContent(iSiPM.at(ii)+1,content>0?content:0);
   }
   if(DoDraw) image->Draw(opt);
   image->GetXaxis()->SetRangeUser(xrange[0]-2.,xrange[1]+2.);
   image->GetYaxis()->SetRangeUser(yrange[0]-2.,yrange[1]+2.);
   //image->GetXaxis()->SetRangeUser(-90.,90.);
   //image->GetYaxis()->SetRangeUser(-90.,90.);
   tlist->Add(image);
   return image;
}

TGraph2D* WFCTAEvent::Draw3D(int type,const char* opt,double threshold,int ViewOpt){
   TGraph2D* array=new TGraph2D();
   double rmin[3]={1.0e10,1.0e10,1.0e100};
   double rmax[3]={-1.0e10,-1.0e10,-1.0e100};

   int size=iSiPM.size();
   for(int ii=0;ii<size;ii++){
      int isipm=iSiPM.at(ii);
      if(!CleanImage(ii,iTel,true)) continue;
      double ImageX,ImageY;
      ImageX=WCamera::GetSiPMX(isipm)/WFTelescope::FOCUS/PI*180;
      ImageY=WCamera::GetSiPMY(isipm)/WFTelescope::FOCUS/PI*180;
      if(iTel==5&&rabbitTime<1570680000) {ImageX*=-1; ImageY*=-1;}
      double content=GetContent(ii,iTel,type,true);
      if(type==0&&ii<ADC_Cut.size()) content=content>adccuttrigger?1.:0.;
      array->SetPoint(array->GetN(),ImageX,ImageY,content);

      if(ImageX<rmin[0]) rmin[0]=ImageX;
      if(ImageX>rmax[0]) rmax[0]=ImageX;
      if(ImageY<rmin[1]) rmin[1]=ImageY;
      if(ImageY>rmax[1]) rmax[1]=ImageY;
      if(content<rmin[2]) rmin[2]=content;
      if(content>rmax[2]) rmax[2]=content;
      //printf("Draw3D: ii=%d ImageX=%lf ImageY=%lf tmin=%le tmax=%le\n",ii,ImageX,ImageY,tmin,tmax);
   }
   rmin[0]=-8.5; rmax[0]=8.5;
   rmin[1]=-8.5; rmax[1]=8.5;

   //TObjArray* array=new TObjArray();
   //double rmin[3]={1.0e10,1.0e10,1.0e100};
   //double rmax[3]={-1.0e10,-1.0e10,-1.0e100};
   //double tmin=1.0e20,tmax=-1;
   //int size=iSiPM.size();
   //for(int ii=0;ii<size;ii++){
   //   if(mypeak.at(ii)<tmin) tmin=mypeak.at(ii);
   //   if(mypeak.at(ii)>tmax) tmax=mypeak.at(ii);
   //}

   //for(int ii=0;ii<size;ii++){
   //   int isipm=iSiPM.at(ii);
   //   if(!CleanImage(isipm,0,3)) continue;
   //   double ImageX,ImageY;
   //   ImageX=WCamera::GetSiPMX(isipm)/WFTelescope::FOCUS/PI*180;
   //   ImageY=WCamera::GetSiPMY(isipm)/WFTelescope::FOCUS/PI*180;

   //   if(ImageX<rmin[0]) rmin[0]=ImageX;
   //   if(ImageX>rmax[0]) rmax[0]=ImageX;
   //   if(ImageY<rmin[1]) rmin[1]=ImageY;
   //   if(ImageY>rmax[1]) rmax[1]=ImageY;
   //   if(tmin<rmin[2]) rmin[2]=tmin;
   //   if(tmax>rmax[2]) rmax[2]=tmax;
   //   //printf("Draw3D: ii=%d ImageX=%lf ImageY=%lf tmin=%le tmax=%le\n",ii,ImageX,ImageY,tmin,tmax);

   //   TMarker3DBox* box=new TMarker3DBox();
   //   box->SetDirection(0,0);
   //   box->SetPosition(ImageX,ImageY,mypeak.at(ii));
   //   box->SetSize(0.25,0.25,0.1);
   //   box->SetFillColor(2);
   //   array->Add((TObject*)box);
   //}
   //rmin[0]=-8.5; rmax[0]=8.5;
   //rmin[1]=-8.5; rmax[1]=8.5;
   ////printf("Draw3D range: xx={%lf,%lf},yy={%lf,%lf},zz={%le,%le}\n",rmin[0],rmax[0],rmin[1],rmax[1],rmin[2],rmax[2]);

   TCanvas* cc=new TCanvas();
   TView *view = TView::CreateView(1,rmin,rmax);
   view->ShowAxis();
   if(ViewOpt==1) view->Front();
   else if(ViewOpt==2) view->Side();
   else if(ViewOpt==3) view->Top();
   cc->SetView(view);

   tlist->Add(cc);
   graphlist.push_back((TGraph*)array);
   if(DoDraw) array->Draw(opt);
   return array;
}

bool WFCTAEvent::IsLed(int nfire_threshold){
   int nfire=0;
   int size=iSiPM.size();
   for(int ii=0;ii<size;ii++){
      if(!CleanImage(ii,iTel,true)) continue;
      nfire++;
   }
   //printf("event=%d nfire=%d\n",iEvent,nfire);
   return nfire>=nfire_threshold;
}
bool WFCTAEvent::IsLaser(){
   return false;
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
		if(FullImageSiPM.at(ii)==0||FullImageSiPM.at(ii)==1023)
		{
			if(fNpixfriends.at(ii)<=2){continue;}
		}
		else
		{
			if(fNpixfriends.at(ii)<=cleanPix){continue;}//use this to clean image
		}
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
		if(FullImageSiPM.at(ii)==0||FullImageSiPM.at(ii)==1023)
		{
			if(fNpixfriends.at(ii)<=2){continue;}
		}
		else
		{
			if(fNpixfriends.at(ii)<=cleanPix){continue;}//use this to clean image
		}
		CleanImagePe.push_back(FullImagePe.at(ii));
		CleanImageSiPM.push_back(FullImageSiPM.at(ii));
		CleanImageX.push_back(FullImageX.at(ii));
		CleanImageY.push_back(FullImageY.at(ii));
	}
}
