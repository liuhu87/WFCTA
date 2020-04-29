#include <cmath>
#include <iostream>
#include "WFCTAEvent.h"
#include "WFCamera.h"
#include "WFCone.h"
#include "TH2Poly.h"
#include "TMarker3DBox.h"
#include "TAxis3D.h"
#include <TCanvas.h>
#include <TView3D.h>
#include <TSystem.h>
#include "TF1.h"
#include "Laser.h"
#include "slalib.h"
#include "RotateDB.h"
#include "StatusDB.h"
#include "CalibWFCTA.h"
#include "LHChain.h"
#include "common.h"

using namespace std;

ClassImp(WFCTAEvent);

const int WFCTAEvent::cleanPix=3;
WFCTAEvent* WFCTAEvent::_Head=0;
TTree* WFCTAEvent::_Tree=0;
LHChain* WFCTAEvent::_Chain=0;
const char* WFCTAEvent::_Name="WFCTAEvent";
TBranch* WFCTAEvent::bAll=0;
TBranch* WFCTAEvent::bmcevent=0;
TBranch* WFCTAEvent::bledevent=0;
TBranch* WFCTAEvent::blaserevent=0;
int WFCTAEvent::_Entry=-1;
int WFCTAEvent::jdebug=0;
bool WFCTAEvent::DoDraw=false;
double WFCTAEvent::adccuttrigger=350.;
double WFCTAEvent::npetrigger=45.;
int WFCTAEvent::nfiretrigger=3;
int WFCTAEvent::CalibType=0;
WFCTAEvent::WFCTAEvent():TSelector()
{
	packCheck.reserve(MAXPMT);
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
	LaserBaseH.reserve(MAXPMT);
	LaserBaseL.reserve(MAXPMT);
	AdcH.reserve(MAXPMT);
	AdcL.reserve(MAXPMT);
	LaserAdcH.reserve(MAXPMT);
	LaserAdcL.reserve(MAXPMT);
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
	LaserBaseH.resize(MAXPMT);
	LaserBaseL.resize(MAXPMT);
	AdcH.resize(MAXPMT);
	AdcL.resize(MAXPMT);
	LaserAdcH.resize(MAXPMT);
	LaserAdcL.resize(MAXPMT);
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
	   tlist=0;
	}
	graphlist.clear();
}

void WFCTAEvent::Init()
{
   tlist=0;

   iTel=-1;
   iEvent=-1;
   eEvent=-1;
   rabbitTime=0;
   rabbittime=0;
   big_pack_lenth=-1;
   n_fired=-1;
   n_Channel=-1;
   iSiPM.clear();
   packCheck.clear();
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
   htlong=0;

   Type=0;

   mcevent.Init();
   ledevent.Init();
   laserevent.Init();
}
void WFCTAEvent::EventInitial()
{
	if(tlist){
	   int size=tlist->Capacity();
	   for(int ii=0;ii<size;ii++){
	      delete tlist->At(ii);
	   }
	   for(int ii=0;ii<size;ii++){
	      tlist->RemoveLast();
	   }
	}
	for(int ii=0;ii<graphlist.size();ii++){
	   delete (graphlist.at(ii));
	}
	if(graphlist.size()>0) graphlist.clear();

	iTel=-1;
	merge_size=-1;
	iEvent=-1;
	eEvent=-1;
	rabbitTime=0;
	rabbittime=0;
	big_pack_lenth=-1;
	n_fired=-1;
	n_Channel=-1;
	iSiPM.clear();
	packCheck.clear();
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
	LaserBaseH.clear();
	LaserBaseL.clear();
	AdcH.clear();
	AdcL.clear();
	LaserAdcH.clear();
	LaserAdcL.clear();
	SatH.clear();
	SatL.clear();

	for(int j=0;j<28;j++){
		Npoint[j]=j;
		for(int i=0;i<1024;i++){
			pulsehigh[i][j] = 0;
			pulselow[i][j] = 0;
		}
	}

	if(gDraw) {delete gDraw; gDraw=0;}
	if(gDrawErr) {delete gDrawErr; gDrawErr=0;}
	if(minimizer) {delete minimizer; minimizer=0;}
        if(htlong) {delete htlong; htlong=0;}

	Type=0;

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
void WFCTAEvent::SetLHChain(LHChain* chain){
	_Chain=chain;
}
void WFCTAEvent::CreateBranch(TTree *tree, int branchSplit){
	if(tree){
		Head()=this;
		tree->Branch(BranchName(),"WFCTAEvent",&Head(),320000,branchSplit);
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
bool WFCTAEvent::GetAllContents(int entry){
	if(entry<0) entry=_Entry;
	else _Entry=entry;
	if(_Entry<0) return false;
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
        Type=0;
	return ncount>0;
}
const char* WFCTAEvent::GetFileName(){
   return _Chain?_Chain->GetFileName():0;
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
         LaserAdcH.push_back(mcevent.TubeSignal[itel][ii]*WFCTAMCEvent::fAmpHig);
         LaserAdcL.push_back(mcevent.TubeSignal[itel][ii]*WFCTAMCEvent::fAmpLow);
         eSatH.push_back(false);
         eSatL.push_back(false);
         SatH.push_back(false);
         SatL.push_back(false);

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
double WFCTAEvent::GetImageXYCoo(int isipm,double &ImageX,double &ImageY,double focus,bool Isdegree){
   if(isipm<0||isipm>=NSIPM) return -1;
   if(WCamera::SiPMMAP[0][0]==0) WCamera::SetSiPMMAP();
   ImageX=WCamera::GetSiPMX(isipm);
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
int WFCTAEvent::GetMaxADCBin(int itel){
	int res=-1;
	double maxadc=-1;
	for(int ii=0;ii<iSiPM.size();ii++){
		double yy=LaserAdcH.at(ii);
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

bool WFCTAEvent::IsSaturated(int isipm,int itel,bool IsHigh,bool IsIndex){
   int size=iSiPM.size();
   int start=IsIndex?isipm:0;
   int end=IsIndex?isipm:(size-1);
   double content=2;
   for(int ii=start;ii<=end;ii++){
      if((!IsIndex)&&(iSiPM.at(ii)!=isipm)) continue;
      if(IsHigh){
         double baseh=LaserBaseH.at(ii);
         if(ii<LaserAdcH.size()&&ii<LaserBaseH.size()&&ii<winsum.size()){
            content=(LaserAdcH.at(ii)>6000||(winsum.at(ii)+LaserBaseH.at(ii)*4)>9000)?1.:-1.;
         }
         //if(ii<LaserAdcH.size()&&ii<eSatH.size()) content=(LaserAdcH.at(ii)>6000||eSatH.at(ii))?1.:-1.;
      }
      else{
         if(ii<eSatL.size()) content=(eSatL.at(ii))?1.:-1.;
      }
   }
   return (content>0);
}
double WFCTAEvent::GetContent(int isipm,int itel,int type,bool IsIndex,bool IsFit){
   itel=0;
   double content=0;
   int size=iSiPM.size();
   int start=IsIndex?isipm:0;
   int end=IsIndex?isipm:(size-1);
   if(IsIndex&&(isipm<0||isipm>=size)) return content;
   for(int ii=start;ii<=end;ii++){
      if((!IsIndex)&&(iSiPM.at(ii)!=isipm)) continue;
      bool ishig=false;
      if(type==0){
         if(ii<ADC_Cut.size()) content=ADC_Cut.at(ii);
         else content=-1;
      }
      if(type==1){
         if(ii<eAdcH.size()) content=eAdcH.at(ii)/WFCTAMCEvent::fAmpHig;
         else content=-1;
         ishig=true;
      }
      if(type==2){
         if(ii<eAdcL.size()) content=eAdcL.at(ii)/WFCTAMCEvent::fAmpLow;
         else content=-1;
      }
      if(type==3){
         if(ii<LaserAdcH.size()) content=LaserAdcH.at(ii)/WFCTAMCEvent::fAmpHig;
         else content=-1;
         ishig=true;
      }
      if(type==4){
         if(ii<LaserAdcL.size()) content=LaserAdcL.at(ii)/WFCTAMCEvent::fAmpLow;
         else content=-1;
      }
      if(type==5){
         if(ii<PeakPosH.size()) content=PeakPosH.at(ii);
         else content=-1;
         ishig=true;
      }
      if(type==6){
         if(ii<PeakPosL.size()) content=PeakPosL.at(ii);
         else content=-1;
      }
      if(type==7){
         if(ii<PeakAmH.size()) content=PeakAmH.at(ii);
         else content=-1;
         ishig=true;
      }
      if(type==8){
         if(ii<PeakAmL.size()) content=PeakAmL.at(ii);
         else content=-1;
      }
      if(type==9){
         content=IsSaturated(ii,itel,true,true)?1:-1;
         ishig=true;
      }
      if(type==10){
         content=IsSaturated(ii,itel,false,true)?1:-1;
      }
      if(type==11){
         ishig=(!IsSaturated(ii,itel,true,true));
         if(ishig) content=(ii<LaserAdcH.size())?(LaserAdcH.at(ii)/WFCTAMCEvent::fAmpHig):-1;
         else      content=(ii<LaserAdcL.size())?(LaserAdcL.at(ii)/WFCTAMCEvent::fAmpLow):-1;
      }
      if(type==12){
         ishig=(!IsSaturated(ii,itel,true,true));
         bool satl=IsSaturated(ii,itel,false,true);
         if((!ishig)&&satl) content=-1;
         else if(!ishig) content=(ii<LaserAdcL.size())?(LaserAdcL.at(ii)/WFCTAMCEvent::fAmpLow):-1;
         else content=(ii<LaserAdcH.size())?(LaserAdcH.at(ii)/WFCTAMCEvent::fAmpHig):-1;
      }
      if(type==13){
         content=(ii<AdcH.size())?(AdcH.at(ii)/WFCTAMCEvent::fAmpHig):-1;
         ishig=true;
      }
      if(type==14){
         content=(ii<AdcL.size())?(AdcL.at(ii)/WFCTAMCEvent::fAmpLow):-1;
      }
      if(type==15){
         ishig=(!IsSaturated(ii,itel,true,true));
         bool satl=IsSaturated(ii,itel,false,true);
         if((!ishig)&&satl) content=-1;
         else if(!ishig) content=(ii<PeakPosL.size())?PeakPosL.at(ii):-1;
         else content=(ii<PeakPosH.size())?PeakPosH.at(ii):-1;
      }
      if(((CalibType&0x3)!=0&&content>0)&&(!IsFit)){
         if((type>=1&&type<=4)||(type>=7&&type<=8)||(type>=11&&type<=14)){
            TDirectory* gdir=gDirectory;
            if(CalibWFCTA::UseSiPMCalibVer==1){
               if(jdebug>5) printf("WFCTAEvent::GetContent: iTel=%d Time=%ld SiPM=%d(%d of %d) cont_before=%lf\n",iTel,rabbitTime,iSiPM.at(ii),ii,iSiPM.size(),content);
               content*=(ishig?WFCTAMCEvent::fAmpHig:WFCTAMCEvent::fAmpLow);
               content=CalibWFCTA::GetHead()->DoCalibSiPM(iTel,iSiPM.at(ii),content,0,rabbitTime+rabbittime*2.0e-8,CalibType,ishig?(Type+10):Type);
               if(jdebug>5) printf("WFCTAEvent::GetContent: iTel=%d Time=%ld SiPM=%d(%d of %d) cont_after=%lf\n",iTel,rabbitTime,iSiPM.at(ii),ii,iSiPM.size(),content);
            }
            else if(CalibWFCTA::UseSiPMCalibVer>1){
               char filename[200]="";
               bool exist=CommonTools::GetStatusFile(filename,(char*)GetFileName());
               double temperature=StatusDB::GetHead()->GetPreTemp(iTel,rabbitTime,iSiPM.at(ii),exist?filename:0);
               if(CalibWFCTA::UseSiPMCalibVer==2){
                  double sum=0;
                  int nn=0;
                  for(int i0=0;i0<1024;i0++){
                     double itemp=StatusDB::GetHead()->PreTemp[i0];
                     if(itemp<-100) continue;
                     sum+=itemp;
                     nn++;
                  }
                  if(nn>0) temperature=sum/nn;
               }
               if(jdebug>5) printf("WFCTAEvent::GetContent: iTel=%d Time=%ld SiPM=%d(%d of %d) temp=%lf cont_before=%lf\n",iTel,rabbitTime,iSiPM.at(ii),ii,iSiPM.size(),temperature,content);
               if(temperature>-100) content=CalibWFCTA::GetHead()->DoCalibSiPM(iTel,iSiPM.at(ii),content,temperature,rabbitTime,CalibType,Type);
               if(jdebug>5) printf("WFCTAEvent::GetContent: iTel=%d Time=%ld SiPM=%d(%d of %d) temp=%lf cont_after=%lf\n",iTel,rabbitTime,iSiPM.at(ii),ii,iSiPM.size(),temperature,content);
            }
            if(gdir) gdir->cd();
         }
      }
   }
   return content;
}
double WFCTAEvent::GetContentError(int isipm,int itel,int type,bool IsIndex,bool IsFit){
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
      if(type==5){
         if(ii<PeakPosH.size()) econtent=1;
         else econtent=0;
      }
      if(type==6){
         if(ii<PeakPosL.size()) econtent=1;
         else econtent=0;
      }
      if(type==9){
         if(ii<SatH.size()&&ii<eSatH.size()) econtent=0;
         else if(ii<SatH.size()) econtent=0;
         else econtent=0;
      }
      if(type==10){
         if(ii<SatL.size()&&ii<eSatL.size()) econtent=0;
         else if(ii<SatL.size()) econtent=0;
         else econtent=0;
      }
   }
   return econtent;
}
double WFCTAEvent::GetTotalPe(double &error,int &ncontent,int itel,int type,bool DoClean,bool IsFit){
   ncontent=0;
   double sum=0;
   double esum=0;
   for(int ii=0;ii<iSiPM.size();ii++){
      if(DoClean) {if(!CleanImage(ii,iTel,true)) continue;}
      double content=GetContent(ii,itel,type,true,IsFit);
      double econtent=GetContentError(ii,itel,type,true,IsFit);
      if(content<=0) continue;
      sum+=content;
      esum=sqrt(pow(esum,2)+pow(econtent,2));
      ncontent++;
   }
   error=esum;
   return sum;
}

bool WFCTAEvent::CleanImage(int isipm,int itel,bool IsIndex,bool IsFit){
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
      double content0=GetContent(IsIndex?ii:isipm0,itel,0,IsIndex,IsFit);
      //bool sat0=(GetContent(IsIndex?ii:isipm0,itel,9,IsIndex,IsFit)>0.5);
      double content=GetContent(IsIndex?ii:isipm0,itel,11,IsIndex,IsFit);
      //if(ADC_Cut.size()>ii) {if(content0>0&&content0<trig0) continue;}
      //else{if(content<trig1) continue;}
      if(content<trig1) continue;
      double ImageXi=0,ImageYi=0;
      GetImageXYCoo(ii,ImageXi,ImageYi);
      int nneigh=0;
      for(int jj=0;jj<size;jj++){
         if(jj==ii) continue;
         double ImageXj=0,ImageYj=0;
         GetImageXYCoo(jj,ImageXj,ImageYj);
         double dist=sqrt(pow(ImageXi-ImageXj,2)+pow(ImageYi-ImageYj,2));
         if(dist>0.6) continue;
         bool sat1=(GetContent(jj,itel,9,true,IsFit)>0.5);
         double contentj0=GetContent(jj,itel,0,true,IsFit);
         double contentj=GetContent(jj,itel,11,true,IsFit);
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
   double content[1024];
   double econtent[1024];
   double allcontent=0;
   double eallcontent=0;
   int itel=(int)(par[0]+0.5);
   int type=(int)(par[1]+0.5);
   bool IsFit=((int)(par[4]+0.5));
   for(int ii=0;ii<size;ii++){
      int isipm=iSiPM.at(ii);
      if(CleanImage(ii,itel,true,IsFit)){
         content[ii]=GetContent(ii,itel,type,true,IsFit);
         econtent[ii]=GetContentError(ii,itel,type,true,IsFit);
      }
      else{
         content[ii]=0;
         econtent[ii]=0;
      }
      if(content[ii]>=0){
         allcontent+=content[ii];
         eallcontent=sqrt(pow(eallcontent,2)+pow(econtent[ii],2));
      }
   }
   double chi2=0;
   int ndof=0;
   //printf("\n");
   for(int ii=0;ii<size&&(allcontent>0);ii++){
      int isipm=iSiPM.at(ii);
      double ImageX,ImageY;
      double dcell=GetImageXYCoo(isipm,ImageX,ImageY,WFTelescope::FOCUS,false);
      ImageX*=-1;
      if(dcell<0) continue;
      double distance=(ImageX*sin(par[3])-ImageY*cos(par[3])+par[2]);
      double probi=content[ii]/allcontent;
      //printf("isipm=%d dcell=%.2lf XY={%.2lf,%.2lf} dist=%.2lf prob=%.2e\n",isipm,dcell/PI*180,ImageX/PI*180,ImageY/PI*180,distance/PI*180,probi);
      chi2+=probi*pow(distance/(dcell/2.),2);
      ndof++;
   }
   //printf("chi2=%.2e(par=%.2lf,%.2lf)\n\n",chi2,par[3]/PI*180,par[2]/PI*180);
   ////for test
   //for(int ii=0;ii<9;ii++){
   //   chi2+=pow(((ii-4)*sin(par[3])-1.*cos(par[3])+par[2])/1.,2)*1;
   //   ndof++;
   //   sum+=1;
   //}

   //chi2/=(chi2/ndof);
   return chi2;
}
/*double WFCTAEvent::InterfaceT(const double* par){
   if(!htlong){
      htlong=GetDistribution(true,(int)(par[0]+0.5),5,false);
      htlong=CorrTimeDist(htlong,CheckLaser(),CheckMC());
   }
   int size=iSiPM.size();
   double content[1024];
   double econtent[1024];
   double allcontent=0;
   double eallcontent=0;
   for(int ii=0;ii<size;ii++){
      int isipm=iSiPM.at(ii);
      if(CleanImage(ii,(int)(par[0]+0.5),true)){
         content[ii]=GetContent(ii,(int)(par[0]+0.5),(int)(par[1]+0.5),true);
         econtent[ii]=GetContentError(ii,(int)(par[0]+0.5),(int)(par[1]+0.5),true);
      }
      else{
         content[ii]=0;
         econtent[ii]=0;
      }
      allcontent+=content[ii];
      eallcontent=sqrt(pow(eallcontent,2)+pow(econtent[ii],2));
   }
   double chi2=0;
   int ndof=0;
   for(int ii=0;ii<size&&(allcontent>0);ii++){
      int isipm=iSiPM.at(ii);
      double ImageX,ImageY;
      double dcell=GetImageXYCoo(isipm,ImageX,ImageY,WFTelescope::FOCUS,false);
      if(dcell<0) continue;
      double distance=(ImageX*sin(par[3])-ImageY*cos(par[3])+par[2]);
      double probi=content[ii]/allcontent;
      chi2+=probi*pow(distance/(dcell/2.),2);
      ndof++;
   }

   return chi2;
}*/
bool WFCTAEvent::DoFit(int itel,int type,bool force){
   itel=0;
   if(itel<0||itel>=WFTelescopeArray::CTNumber) return false;
   if(minimizer&&(!force)) return true;
   int size=iSiPM.size();
   int nbin=0;
   double mx,my,sx,sy,sxy,nn;
   bool IsFit=true;
   for(int ii=0;ii<size;ii++){
      int isipm0=iSiPM.at(ii);
      if(!CleanImage(ii,itel,true,IsFit)) continue;
      double content=GetContent(ii,itel,type,true,IsFit);
      if(content<=0) continue;
      double ImageX,ImageY;
      GetImageXYCoo(isipm0,ImageX,ImageY,-1,false);
      ImageX*=-1;

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
   ROOT::Math::Functor f(this,&WFCTAEvent::Interface,5);
   #endif
   minimizer->SetFunction(f);
   minimizer->SetFixedVariable(0,"iTel",0);
   minimizer->SetFixedVariable(1,"Type",4);
   //minimizer->SetLimitedVariable(2,"bb",b,fabs(0.01*b),-1.e6,1.e6);
   //minimizer->SetLimitedVariable(3,"kk",a,fabs(0.01*a),-1.e6,1.e6);
   double cc0=a>=0?b/sqrt(1+a*a):-b/sqrt(1+a*a);
   double phi0=a>=0?atan(a):(PI+atan(a));
   minimizer->SetLimitedVariable(2,"cc",cc0,0.5/180.*PI,-12./180.*PI,12./180.*PI);
   minimizer->SetLimitedVariable(3,"phi",phi0,fabs(1./180.*PI),0,0.999*PI);
   //minimizer->SetLimitedVariable(2,"cc",1,0.01,-5,5);
   //minimizer->SetLimitedVariable(3,"phi",0,fabs(0.1/180.*PI),-PI/2,PI/2);
   //minimizer->SetFixedVariable(3,"phi",0.);
   minimizer->SetFixedVariable(4,"IsFit",IsFit?1.:0.);
   minimizer->Minimize();
   minimizer->Hesse();
   //printf("LinearFit: kk=%lf bb=%lf cc={%lf,%lf} phi={%lf,%lf}\n",a,b/PI*180,cc0/PI*180,minimizer->X()[2]/PI*180,phi0/PI*180,minimizer->X()[3]/PI*180.);
   //printf("LinearFit: phi={%lf,%lf} CC={%lf,%lf}\n",minimizer->X()[3],minimizer->Errors()[3],minimizer->X()[2],minimizer->Errors()[2]);
   return true;
}

bool WFCTAEvent::GetCrossCoor(double x,double y,double CC,double phi,double &x0,double &y0){
   CC=CC/PI*180;
   x0=cos(phi)*(x*cos(phi)+y*sin(phi))-CC*sin(phi);
   y0=sin(phi)*(x*cos(phi)+y*sin(phi))+CC*cos(phi);
   return true;
}
bool WFCTAEvent::GetCrossCoor(double x,double y,double &x0,double &y0){
   if(!minimizer) return false;
   double CC=minimizer->X()[2];
   double phi=minimizer->X()[3];
   return GetCrossCoor(x,y,CC,phi,x0,y0);
}

int WFCTAEvent::IsInside(double xx,double yy){
   double ImageX0,ImageY0;
   double dcell=GetImageXYCoo(PIX/2,ImageX0,ImageY0);
   int IndexY=-1,IndexX=-1;
   double xmin=100,xmax=-100;
   double ymin=100,ymax=-100;
   double mindistx=1000,mindisty=1000;
   int indexminx=-1,indexminy=-1;
   for(int Yi=0;Yi<PIX;Yi++){
      double ImageX,ImageY;
      GetImageXYCoo(PIX/2+Yi*PIX,ImageX,ImageY);
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
      if((yy>ymax||yy<ymin)||indexminy<0) return -1;
      else{
         for(int Xi=0;Xi<PIX;Xi++){
            double ImageX,ImageY;
            GetImageXYCoo(Xi+indexminy*PIX,ImageX,ImageY);
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
         if((xx>xmax||xx<xmin)||indexminx<0) return -1;
         else return -2;
      }
   }
   for(int Xi=0;Xi<PIX;Xi++){
      double ImageX,ImageY;
      GetImageXYCoo(Xi+IndexY*PIX,ImageX,ImageY);
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
      if((xx>xmax||xx<xmin)||indexminx<0) return -1;
      else return -2;
   }
   else return IndexY*PIX+IndexX;
}
bool WFCTAEvent::GetRange(double CC,double phi,double XY1[4],double XY2[4]){
   if(jdebug>0) printf("WFCTAEvent::GetRange:input,CC=%lf phi=%lf\n",CC/PI*180,phi/PI*180);
   CC*=180./PI;
   double margin=1.0e-5;
   if(fabs(phi)<(0.1/180.*PI)){
      double yy=CC/cos(phi);
      int IndexY=-1;
      for(int Yi=0;Yi<PIX;Yi++){
         double ImageX,ImageY;
         double dcell=GetImageXYCoo(PIX/2+Yi*PIX,ImageX,ImageY);
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
      double dcell=GetImageXYCoo(0+IndexY*PIX,XY1[0],yy);
      XY1[0]-=dcell/2.;
      XY1[2]=XY1[0]*cos(phi)+XY1[1]*sin(phi);
      XY1[3]=1;
      dcell=GetImageXYCoo((PIX-1)+IndexY*PIX,XY2[0],yy);
      XY2[0]+=dcell/2.;
      XY2[2]=XY2[0]*cos(phi)+XY2[1]*sin(phi);
      XY2[3]=1;
   }
   else{
      int exist=0;
      int index=0;
      while(index>=0&&index<PIX){
         double ImageX,ImageY;
         double dcell=GetImageXYCoo(PIX/2+index*PIX,ImageX,ImageY);
         double y1=ImageY+dcell/2;
         double y2=ImageY-dcell/2;
         double x1=(y1*cos(phi)-CC)/sin(phi);
         double x2=(y2*cos(phi)-CC)/sin(phi);
         bool inside1=IsInside(x1,y1-margin)>=0;
         bool inside2=IsInside(x2,y2+margin)>=0;
         if(jdebug>1) printf("WFCTAEvent::GetRange: index=%d xy1={%lf,%lf} xy2={%lf,%lf} inside={%d,%d}\n",index,x1,y1,x2,y2,inside1,inside2);
         if(inside1){
            XY1[1]=y1;
            XY1[0]=x1;
            XY1[2]=XY1[0]*cos(phi)+XY1[1]*sin(phi);
            XY1[3]=1;
            exist=1;
            /*if(!inside2) exist=1;
            else{
               XY2[1]=y2;
               XY2[0]=x2;
               XY2[2]=XY2[0]*cos(phi)+XY2[1]*sin(phi);
               XY2[3]=1;
               if(XY1[2]>XY2[2]){
                  double nbuff[4]={XY2[0],XY2[1],XY2[2],XY2[3]};
                  for(int ii=0;ii<4;ii++) XY2[ii]=XY1[ii];
                  for(int ii=0;ii<4;ii++) XY1[ii]=nbuff[ii];
               }
               exist=2;
            }*/
            break;
         }
         else if(inside2){
            dcell=GetImageXYCoo(((x1>0)?(PIX-1):0)+index*PIX,ImageX,ImageY);
            XY1[0]=ImageX+((x1>0)?1:-1)*dcell/2;
            XY1[1]=(XY1[0]*sin(phi)+CC)/cos(phi);
            XY1[2]=XY1[0]*cos(phi)+XY1[1]*sin(phi);
            XY1[3]=1;
            exist=1;
            break;
         }
         else if(x1*x2<0){
            dcell=GetImageXYCoo(0+index*PIX,ImageX,ImageY);
            XY1[0]=ImageX-dcell/2;
            XY1[1]=(XY1[0]*sin(phi)+CC)/cos(phi);
            XY1[2]=XY1[0]*cos(phi)+XY1[1]*sin(phi);
            XY1[3]=1;
            dcell=GetImageXYCoo(PIX-1+index*PIX,ImageX,ImageY);
            XY2[0]=ImageX+dcell/2;
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
         double dcell=GetImageXYCoo(PIX/2+index*PIX,ImageX,ImageY);
         double y1=ImageY-dcell/2;
         double y2=ImageY+dcell/2;
         double x1=(y1*cos(phi)-CC)/sin(phi);
         double x2=(y2*cos(phi)-CC)/sin(phi);
         bool inside1=IsInside(x1,y1+margin)>=0;
         bool inside2=IsInside(x2,y2-margin)>=0;
         if(inside1){
            XY2[1]=y1;
            XY2[0]=x1;
            XY2[2]=XY2[0]*cos(phi)+XY2[1]*sin(phi);
            XY2[3]=1;
            exist=2;
            break;
         }
         else if(inside2){
            dcell=GetImageXYCoo(((x1>0)?(PIX-1):0)+index*PIX,ImageX,ImageY);
            XY2[0]=ImageX+((x1>0)?1:-1)*dcell/2;
            XY2[1]=(XY2[0]*sin(phi)+CC)/cos(phi);
            XY2[2]=XY2[0]*cos(phi)+XY2[1]*sin(phi);
            XY2[3]=1;
            exist=2;
            break;
         }
         else if(x1*x2<0){
            dcell=GetImageXYCoo(0+index*PIX,ImageX,ImageY);
            XY1[0]=ImageX-dcell/2;
            XY1[1]=(XY1[0]*sin(phi)+CC)/cos(phi);
            XY1[2]=XY1[0]*cos(phi)+XY1[1]*sin(phi);
            XY1[3]=1;
            dcell=GetImageXYCoo(PIX-1+index*PIX,ImageX,ImageY);
            XY2[0]=ImageX+dcell/2;
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

      /*double ImageX,ImageY;
      double dcell=GetImageXYCoo(PIX/2+0*PIX,ImageX,ImageY);
      XY1[1]=ImageY+dcell/2;
      XY1[0]=(XY1[1]*cos(phi)-CC)/sin(phi);
      XY1[2]=XY1[0]*cos(phi)+XY1[1]*sin(phi);
      XY1[3]=-1;
      bool inside1=IsInside(XY1[0],XY1[1]-margin)>=0;
      dcell=GetImageXYCoo(PIX/2+(PIX-1)*PIX,ImageX,ImageY);
      XY2[1]=ImageY-dcell/2;
      XY2[0]=(XY2[1]*cos(phi)-CC)/sin(phi);
      XY2[2]=XY2[0]*cos(phi)+XY2[1]*sin(phi);
      XY2[3]=-1;
printf("XY1={%lf,%lf,%lf,%lf} XY2={%lf,%lf,%lf,%lf}\n",XY1[0],XY1[1],XY1[2],XY1[3],XY2[0],XY2[1],XY2[2],XY2[3]);
      bool inside2=IsInside(XY2[0],XY2[1]+margin)>=0;
      //printf("inside1=%d inside2=%d\n",inside1,inside2);
      if((!inside1)&&(!inside2)) return false;
      else if(inside1&&inside2){
         ;
      }
      else{
         bool useupper=(inside1&&(!inside2));
         int index=useupper?(PIX-1):0;
         while(index>=0&&index<PIX){
            double ImageX,ImageY;
            GetImageXYCoo(PIX/2+index*PIX,ImageX,ImageY);
            double y1=ImageY+(useupper?-1:1)*dcell/2;
            double x1=(y1*cos(phi)-CC)/sin(phi);
            if(IsInside(x1,y1-(useupper?-1:1)*margin)>=0){
               if(useupper){
                  XY2[0]=x1;
                  XY2[1]=y1;
                  XY2[2]=XY2[0]*cos(phi)+XY2[1]*sin(phi);
                  XY2[3]=1;
               }
               else{
                  XY1[0]=x1;
                  XY1[1]=y1;
                  XY1[2]=XY1[0]*cos(phi)+XY1[1]*sin(phi);
                  XY1[3]=1;
               }
               break;
            }
            double y2=ImageY-(useupper?-1:1)*dcell/2;
            double x2=(y2*cos(phi)-CC)/sin(phi);
            if(IsInside(x1,y1+(useupper?-1:1)*margin)>=0){
               if(useupper){
                  XY2[0]=x2;
                  XY2[1]=y2;
                  XY2[2]=XY2[0]*cos(phi)+XY2[1]*sin(phi);
                  XY2[3]=1;
               }
               else{
                  XY1[0]=x2;
                  XY1[1]=y2;
                  XY1[2]=XY1[0]*cos(phi)+XY1[1]*sin(phi);
                  XY1[3]=1;
               }
               break;
            }

            if(useupper) index--;
            else index++;
         }
      }*/
   }
   if(jdebug>0) printf("WFCTAEvent::GetRange:output,CC=%lf phi=%lf XY1={%lf,%lf} XY2={%lf,%lf}\n",CC,phi/PI*180,XY1[0],XY1[1],XY2[0],XY2[1]);
   if(XY1[2]>XY2[2]){
      double nbuff[4]={XY2[0],XY2[1],XY2[2],XY2[3]};
      for(int ii=0;ii<4;ii++) XY2[ii]=XY1[ii];
      for(int ii=0;ii<4;ii++) XY1[ii]=nbuff[ii];
   }
   return true;
}
TH1F* WFCTAEvent::GetDistribution(bool IsLong,int itel,int type,bool IsWidth,bool CleanEdge){
   itel=0;
   if(itel<0||itel>=WFTelescopeArray::CTNumber) return 0;
   if(!minimizer) return 0;
   double CC=minimizer->X()[2]/PI*180.;
   double phi=minimizer->X()[3];
   double XY1[4],XY2[4];
   bool inside=GetRange(CC/180.*PI,phi,XY1,XY2);
   if(!inside) return 0;
   double margin=24.;
   const int nbin=48;
   static int ihist=0;
   TH1F* hh=new TH1F(Form("%sdis_%s%s_h%d",IsLong?"long":"short",((type>=5&&type<=6)||type==15)?"time":"npe",IsWidth?"_width":"",ihist++),Form(";%s axis [degree];%s",IsLong?"Long":"Short",IsWidth?"Width [degree]":(((type>=5&&type<=6)||type==15)?"Time [ns]":"Npe [pe]")),IsWidth?24:nbin,IsLong?-12.:-6.,IsLong?12.:6.);
   double nfill[nbin];
   //double nmark[nbin];
   for(int ii=0;ii<nbin;ii++){
      nfill[ii]=0;
      //nmark[ii]=false;
   }

   int size=iSiPM.size();
   TH2F* h2d=new TH2F(Form("h2d"),";Long Axis [degree];Short Axis [degree]",24,-12.,12.,128,-8,8);
   for(int ii=0;ii<size;ii++){
      int isipm0=iSiPM.at(ii);
      if(!CleanImage(ii,itel,true)) continue;
      double ImageX=0,ImageY=0;
      double dcell=GetImageXYCoo(isipm0,ImageX,ImageY);
      double ImageX2=ImageX*cos(phi)+ImageY*sin(phi); //long axis
      double ImageY2=-ImageX*sin(phi)+ImageY*cos(phi)-CC; //short axis
      double content=GetContent(ii,itel,3,true);
      double econtent=GetContentError(ii,itel,3,true);
      if(content<=0) continue;
      int ibinx=h2d->GetXaxis()->FindBin(ImageX2);
      int ibiny=h2d->GetYaxis()->FindBin(ImageY2);
      h2d->SetBinContent(ibinx,ibiny,h2d->GetBinContent(ibinx,ibiny)+content);
      h2d->SetBinError(ibinx,ibiny,sqrt(pow(h2d->GetBinError(ibinx,ibiny),2)+pow(econtent,2)));
   }

   double sigm=-1;
   for(int ibin=0;ibin<=h2d->GetNbinsX();ibin++){
      TH1D* hbuff=0;
      if(ibin==0) hbuff=h2d->ProjectionY("hbuff",1,h2d->GetNbinsX());
      else hbuff=h2d->ProjectionY(Form("hbuff_bin%d",ibin),ibin,ibin);
      double sigmi=-1;
      if(hbuff->Integral()>50){
         double allcont=hbuff->Integral();
         int ibin0=hbuff->GetXaxis()->FindBin(0.);
         double acccont=0;
         for(int ibin=0;ibin<(hbuff->GetNbinsX()/2)+2;ibin++){
            acccont+=hbuff->GetBinContent(ibin0+ibin);
            if(ibin>0) acccont+=hbuff->GetBinContent(ibin0-ibin);
            if(fabs(acccont/allcont)>0.99){
               sigmi=hbuff->GetXaxis()->GetBinCenter(ibin0+ibin)-hbuff->GetXaxis()->GetBinCenter(ibin0);
               break;
            }
         }
         //hbuff->Fit("gaus","QS0");
         //TF1* f1=hbuff->GetFunction("gaus");
         //if(f1) sigmi=f1->GetParameter(2)*2;
      }
      delete hbuff;
      if(ibin==0) sigm=sigmi;
      if(IsWidth){
         if(ibin>0&&sigmi>0){
            hh->SetBinContent(ibin,sigmi);
            hh->SetBinError(ibin,0);
         }
      }
      else if(ibin==0){
         break;
      }
   }
   if(h2d) {delete h2d; h2d=0;}
   if(sigm<=0){
      delete hh;
      return 0;
   }
   if(jdebug>5) printf("WFCTAEvent::GetDistribution: sigm=%lf\n",sigm);
   double LRange[2]={-30,30};
   for(int ii=0;ii<3;ii++){
      double xx=XY1[0];
      double yy=XY1[1];
      if(XY1[3]<0) xx=XY1[0]+sigm*(ii-1);
      else yy=XY1[1]+sigm*(ii-1);
      double LL=xx*cos(phi)+yy*sin(phi);
      if(LL>LRange[0]) LRange[0]=LL;
   }
   for(int ii=0;ii<3;ii++){
      double xx=XY2[0];
      double yy=XY2[1];
      if(XY2[3]<0) xx=XY2[0]+sigm*(ii-1);
      else yy=XY2[1]+sigm*(ii-1);
      double LL=xx*cos(phi)+yy*sin(phi);
      if(LL<LRange[1]) LRange[1]=LL;
   }
   if(jdebug>0) printf("WFCTAEvent::GetDistribution: sigm=%lf LRange={%lf,%lf}\n",sigm,LRange[0],LRange[1]);
   if(IsWidth){ //clean
      for(int ibin=1;ibin<=hh->GetNbinsX();ibin++){
         if(hh->GetBinContent(ibin)>0){
            hh->SetBinContent(ibin,0);
            hh->SetBinError(ibin,0);
            break;
         }
      }
      for(int ibin=hh->GetNbinsX();ibin>=1;ibin--){
         if(hh->GetBinContent(ibin)>0){
            hh->SetBinContent(ibin,0);
            hh->SetBinError(ibin,0);
            break;
         }
      }
      for(int ibin=1;ibin<=hh->GetNbinsX();ibin++){
         double xcenter=hh->GetXaxis()->GetBinCenter(ibin);
         if(xcenter<=LRange[0]||xcenter>=LRange[1]){
            hh->SetBinContent(ibin,0);
            hh->SetBinError(ibin,0);
         }
      }
      return hh;
   }

   bool IsTime=((type>=5&&type<=6)||type==15);
   const int ndivide=30;
   for(int ii=0;ii<size;ii++){
      int isipm0=iSiPM.at(ii);
      if(!CleanImage(ii,itel,true)) continue;
      double content=GetContent(ii,itel,type,true);
      double econtent=GetContentError(ii,itel,type,true);
      if(content<=0) continue;
      double norm=GetContent(ii,itel,3,true);
      if(norm<=0) continue;
      double ImageX,ImageY;
      double dcell=GetImageXYCoo(isipm0,ImageX,ImageY);
      bool isedge=WCamera::IsEdge(isipm0);
      for(int idiv1=0;idiv1<ndivide;idiv1++){
         double ImageXi=ImageX-dcell/2.+dcell/ndivide*(idiv1+0.5);
         for(int idiv2=0;idiv2<ndivide;idiv2++){
            double ImageYi=ImageY-dcell/2.+dcell/ndivide*(idiv2+0.5);
            double ImageX2=ImageXi*cos(phi)+ImageYi*sin(phi); //long axis
            double ImageY2=-ImageXi*sin(phi)+ImageYi*cos(phi)-CC; //short axis
            if(IsLong&&fabs(ImageY2)>margin) continue;
            if(CleanEdge) {if(ImageX2<=LRange[0]||ImageX2>=LRange[1]) continue;}
            int ibin=hh->GetXaxis()->FindBin(IsLong?ImageX2:ImageY2);
            if(ibin<1||ibin>nbin) continue;
            if(IsTime){
            hh->SetBinContent(ibin,hh->GetBinContent(ibin)+content*norm/(ndivide*ndivide));
            hh->SetBinError(ibin,sqrt(pow(hh->GetBinError(ibin),2)+pow(econtent,2)*norm/(ndivide*ndivide)));
            nfill[ibin-1]+=norm/(ndivide*ndivide);
            //if(ImageX2>4) printf("isipm=%d idiv={%d,%d} cont={%lf,%lf,%lf}\n",ii,idiv1,idiv2,content,content*norm/(ndivide*ndivide),norm/(ndivide*ndivide));
            }
            else{
            hh->SetBinContent(ibin,hh->GetBinContent(ibin)+content/(ndivide*ndivide));
            hh->SetBinError(ibin,sqrt(pow(hh->GetBinError(ibin),2)+pow(econtent,2)/(ndivide*ndivide)));
            }
         }
      }
   }
   if(IsTime){
      for(int ibin=1;ibin<=nbin;ibin++){
         if(nfill[ibin-1]<=0) continue;
         hh->SetBinContent(ibin,hh->GetBinContent(ibin)/nfill[ibin-1]);
         hh->SetBinError(ibin,hh->GetBinError(ibin)/nfill[ibin-1]);
      }
   }
   //clean edge bins
   if(CleanEdge){
   for(int ibin=1;ibin<hh->GetNbinsX();ibin++){
      //if(nmark[ibin-1]){
      //   hh->SetBinContent(ibin,0);
      //   hh->SetBinError(ibin,0);
      //   printf("unused bin: %d %lf\n",ibin,hh->GetXaxis()->GetBinCenter(ibin));
      //}
      if(hh->GetBinContent(ibin)>0){
         hh->SetBinContent(ibin,0);
         hh->SetBinError(ibin,0);
         break;
      }
   }
   for(int ibin=hh->GetNbinsX();ibin>=1;ibin--){
      if(hh->GetBinContent(ibin)>0){
         hh->SetBinContent(ibin,0);
         hh->SetBinError(ibin,0);
         break;
      }
   }
   }
   if(hh->Integral()<=0){
      delete hh;
      hh=0;
   }
   return hh;
}
TH1F* WFCTAEvent::CorrTimeDist(TH1F* hist,int IsLaser,int IsMC){
   if(IsLaser==0) IsLaser=(CheckLaser())?1:-1;
   if(IsMC==0) IsMC=(CheckMC())?1:-1;
   if(!hist) return hist;
   double timeunit=IsMC?(CommonTools::timebinunit[IsLaser]):20.;
   for(int ibin=0;ibin<=hist->GetNbinsX()+1;ibin++){
      hist->SetBinContent(ibin,hist->GetBinContent(ibin)*timeunit);
      hist->SetBinError(ibin,hist->GetBinError(ibin)*timeunit);
   }
   return hist;
}
int WFCTAEvent::GetSign(bool IsLaser,bool IsMC){
   if(!minimizer) return 0;
   double CC=minimizer->X()[2]/PI*180.;
   double phi=minimizer->X()[3];
   TH1F* hwidth=GetDistribution(true,0,3,true);
   TH1F* htime=GetDistribution(true,0,5,false);
   if((!hwidth)||(!htime)){
      printf("WFCTAEvent::GetSign: hwidth=%p htime=%p\n",hwidth,htime);
      if(hwidth) delete hwidth;
      if(htime) delete htime;
      return 0;
   }
   //use the information from width
   int maxbin=-1;
   double maxcontent=0;
   for(int ibin=1;ibin<hwidth->GetNbinsX();ibin++){
      if(hwidth->GetBinContent(ibin)>maxcontent){
          maxbin=ibin;
          maxcontent=hwidth->GetBinContent(ibin);
      }
   }
   if(maxbin<0){
      if(!hwidth) delete hwidth;
      if(!htime) delete htime;
      return 0;
   }
   int nleft=0,nright=0;
   double aveleft=0,averight=0;
   for(int ibin=1;ibin<hwidth->GetNbinsX();ibin++){
      if(hwidth->GetBinContent(ibin)<=0) continue;
      if(ibin<maxbin){
         nleft++;
         aveleft+=fabs(hwidth->GetXaxis()->GetBinCenter(ibin)-hwidth->GetXaxis()->GetBinCenter(maxbin));
      }
      else if(ibin>maxbin){
         nright++;
         averight+=fabs(hwidth->GetXaxis()->GetBinCenter(ibin)-hwidth->GetXaxis()->GetBinCenter(maxbin));
      }
   }
   if(nleft>0) aveleft/=nleft;
   if(nright>0) averight/=nright;
   //use the information from time
   htime->Fit("pol1","QS0");
   TF1* f1=htime->GetFunction("pol1");
   double slope=f1?(f1->GetParameter(1)):0;
   if(true){
      if(!hwidth) delete hwidth;
      if(!htime) delete htime;
   }

   //the sign of A is determined by the time sequence,because:
   //dL=pow((y_obs*cos(theta)+sin(theta))/sin(PHI)/sin(PHI),2)/(A*nz2)*dPHI
   //and also it can be determined by the Image itself, because the image is not symmetric before and after Xmax
   double slope0=0.1,avedL0=1.;
   int sign_shape,sign_time;
   int dLdt1=(fabs(slope)<slope0)?0:(slope>0?1:-1);
   int dLdt2=(fabs(aveleft-averight)<avedL0)?0:(aveleft<averight?1:-1);
   int dPHIdt=IsLaser?1:(-1);
   sign_time=dLdt1/dPHIdt;
   sign_shape=dLdt1/dPHIdt;

   //printf("slope=%lf ave={%lf,%lf} sign={%lf,%lf}\n",slope,aveleft,averight,sign_time,sign_shape);

   //use time information
   if(IsLaser){
      if(fabs(sign_time)<0.5) return sign_shape;
      else return sign_time;
   }
   else{ //use image information
      if(fabs(sign_shape)<0.5) return sign_time;
      else return sign_shape;
   }
}
/*int WFCTAEvent::GetSign(bool IsLaser,bool IsMC){
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
      GetImageXYCoo(isipm,ImageX,ImageY);
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
}*/
int WFCTAEvent::GetSign(double zenith,double azimuth,double planephi,double nz){
   double Anz;
   double AA=GetApar(Anz,zenith,azimuth,planephi,nz);
   if(AA==0) return 0;
   else return (AA>0)?1:-1;
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

void WFCTAEvent::CooCorr(double incoo[3],double indir[2],double outcoo[3]){
   double margin=1.0e-5;
   double dirout[3]={sin(indir[0])*cos(indir[1]),sin(indir[0])*sin(indir[1]),cos(indir[0])};
   if(fabs(incoo[2])>margin&&fabs(dirout[2])>margin){
      double dz=(0-incoo[2])/dirout[2];
      outcoo[0]=incoo[0]+dirout[0]*dz;
      outcoo[1]=incoo[1]+dirout[1]*dz;
      outcoo[2]=0;
   }
   else{
      outcoo[0]=incoo[0];
      outcoo[1]=incoo[1];
      outcoo[2]=incoo[2];
   }
}
//void WFCTAEvent::CooCorr(double incoo[3],double indir[3],double outcoo[3]){
//   double margin=1.0e-5;
//   if(fabs(incoo[2])>margin&&fabs(indir[2])>margin){
//      double dz=(0-incoo[2])/indir[2];
//      outcoo[0]=incoo[0]+indir[0]*dz;
//      outcoo[1]=incoo[1]+indir[1]*dz;
//      outcoo[2]=0;
//   }
//   else{
//      outcoo[0]=incoo[0];
//      outcoo[1]=incoo[1];
//      outcoo[2]=incoo[2];
//   }
//}
int WFCTAEvent::CalTelDir(double CC,double phi,double planephi,double nz,int &nsol,double ele_out[4],double azi_out[4],int signA[4],double ele_in,double azi_in){
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
         /*double icos=(Isele)?cos(ele_out[ii]-ele_in):cos(azi_out[ii]-azi_in);
         if(icos>maxcos){
            result2=ii;
            maxcos=icos;
         }*/
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

void WFCTAEvent::Getnz(double nz,double &nz1,double &nz2){
   nz1=nz/sqrt(1+nz);
   nz2=1./sqrt(1+nz);
   if(isinf(nz)){
      if(signbit(nz)) nz1=-1;
      else nz1=1;
      nz2=0;
   }
}
void WFCTAEvent::Getnz(double incoo[3],double indir[2],double &planephi,double &nz){
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
   nz=sin(indir[0])*sin(indir[1]-planephi)/cos(indir[0]);
}
void WFCTAEvent::GetCCphi(double zenith,double azimuth,double planephi,double nz,double &CC,double &phi){
   double theta=PI/2-zenith;
   double phi0=azimuth;
   double Anz;
   double AA=GetApar(Anz,zenith,azimuth,planephi,nz);
   CC=Anz*sin(theta)-AA*cos(theta)*sin(phi0-planephi);
   phi=acos(-Anz*cos(theta)-AA*sin(theta)*sin(phi0-planephi));
   if(AA*cos(phi0-planephi)<0) phi=2*PI-phi;
   //if(jdebug>0) printf("zen=%lf azi=%lf planephi=%lf nz=%lf AA=%lf Anz=%lf CC=%lf,phi=%lf\n",zenith/PI*180,azimuth/PI*180,planephi/PI*180,nz,AA,Anz,CC/PI*180,phi/PI*180);
}
void WFCTAEvent::GetCCphi(double zenith,double azimuth,double incoo[3],double indir[2],double &CC,double &phi){
   double planephi,nz;
   Getnz(incoo,indir,planephi,nz);
   return GetCCphi(zenith,azimuth,planephi,nz,CC,phi);
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
   Getnz(incoo,indir,planephi,nz);
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

   /*double x2=ImageCoo[1];
   double y2=ImageCoo[0];
   double px1=(x2*cos(theta)+sin(theta))*cos(phi0-planephi);
   double px2=((x2*cos(theta)+sin(theta))*sin(phi0-planephi)*nz1+nz2*(x2*sin(theta)-cos(theta)));
   double py1=y2*cos(theta)*cos(phi0-planephi)-sin(phi0-planephi);
   double py2=(y2*cos(theta)*sin(phi0-planephi)+cos(phi0-planephi))*nz1+nz2*y2*sin(theta);

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
   }*/
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
   /*double boundary=8.5/180.*PI;
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
            if(jdebug>0) printf("WFCTAEvent::GetXYRange: return=-1, error in calculationg boundaries. line=%d calY=%d {%lf,%lf}\n",iline,(iline/2)==0,xysol[iline][0]/PI*180,xysol[iline][1]/PI*180);
            //return -1;
         }
      }
      if(jdebug>0) printf("WFCTAEvent::GetXYRange: iline=%d index1=%d index2=%d\n",iline,index1,index2);
   }
   if(index1<0||index2<0){
      if(jdebug>0) printf("WFCTAEvent::GetXYRange: return=-4, Image is not in the field of view of telescope CC=%lf phi=%lf\n",CC/PI*180,phi/PI*180);
      return -4;
   }
   double xyboun[2][2]={{xysol[index1][1],xysol[index1][0]},{xysol[index2][1],xysol[index2][0]}};*/

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
   /*double p21=cos(theta)*cos(phi0-planephi);
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
   }*/

   double PHI_low=IsLaser?PHI1:PHI2;
   double PHI_hig=IsLaser?PHI2:PHI1;
   PHI[0]=PHI_low;
   PHI[1]=PHI_hig;

   //from PHI to coor
   bool res1=GetImageXYCoo(zenith,azimuth,planephi,nz,PHI[0],XX[0],YY[0]);
   bool res2=GetImageXYCoo(zenith,azimuth,planephi,nz,PHI[1],XX[1],YY[1]);

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
bool WFCTAEvent::CalRotateZeroPos(double ele_rotate,double azi_rotate,double ele_ref,double azi_ref,double &ele0,double &azi0){
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
bool WFCTAEvent::CalDir_out(double ele_rotate,double azi_rotate,double ele_ref,double azi_ref,double ele_in,double azi_in,double &ele_out,double &azi_out){
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
   double ydir[3]={-nz1*sin(ele_tel)*sin(planephi-azi_tel)-nz2*cos(ele_tel),nz1*cos(planephi-azi_tel),-nz1*cos(ele_tel)*sin(planephi-azi_tel)+nz2*sin(ele_tel)};
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
//double WFCTAEvent::CalTime(double ele_tel,double azi_tel,double incoo[3],double indir[2],double Llong,double time0){
//   double CC,phi;
//   GetCCphi(PI/2-ele_tel,azi_tel,incoo,indir,CC,phi);
//   double ImageCooXY[2];
//   return CalTime();
//}
double WFCTAEvent::CalTime(double CC,double phi,double incoo[3],double nz,double PHIL,double ImageCooXY[2],double time0,double ele_in,double azi_in){
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
}
TGraph* WFCTAEvent::DrawTimeLine(double CC,double phi,double incoo[3],double nz,double PHIL,double time0,double ele_in,double azi_in){
   double x0,y0;
   GetCrossCoor(0.,0.,CC,phi,x0,y0);
   TGraph* gr=new TGraph();
   int np=100;
   for(int ii=0;ii<np;ii++){
      double ll=-10.+(10.+10.)/np*ii;
      double xx=x0+ll*sin(phi);
      double yy=y0+ll*cos(phi);
      double ImageCooXY[2]={xx,yy};
      double time=CalTime(CC,phi,incoo,nz,PHIL,ImageCooXY,time0,ele_in,azi_in);
      if(time<0) continue;
      gr->SetPoint(gr->GetN(),ll,time);
   }
   gr->SetTitle(";Coordinate [degree];Time [ns]");
   gr->SetLineColor(4);
   gr->SetLineWidth(3);
   if(DoDraw) gr->Draw("l");
   return gr;
}
TGraph* WFCTAEvent::DrawTimeLine(double ele_rotate,double azi_rotate,double ele_ref,double azi_ref,double ele_in,double azi_in,double ele_tel,double azi_tel,double incoo[3],double time0){
   double ele_out,azi_out;
   if(!CalDir_out(ele_rotate,azi_rotate,ele_ref,azi_ref,ele_in,azi_in,ele_out,azi_out)) return 0;
   double indir[2]={PI/2-ele_out,azi_out};
   double dir0[3]={sin(indir[0])*cos(indir[1]),sin(indir[0])*sin(indir[1]),cos(indir[0])};
   double incoo0[3]={incoo[0],incoo[1],incoo[2]};
   if(fabs(incoo[2])>1.0e-5&&fabs(dir0[2])>1.0e-5){
      double dl=(0-incoo[2])/dir0[2];
      incoo0[0]=incoo[0]+dir0[0]*dl;
      incoo0[1]=incoo[1]+dir0[1]*dl;
      incoo0[2]=0;
   }
   DoFit(0,4);
   double phi=minimizer->X()[3];
   double CC=minimizer->X()[2];
   TGraph* gr=new TGraph();
   int np=100;
   for(int ii=0;ii<np;ii++){
      double ll=-10.+(10.+10.)/np*ii;
      double xx=ll*sin(phi)+CC/PI*180*cos(phi);
      double yy=ll*cos(phi)-CC/PI*180*sin(phi);
      double ImageCooXY[2]={xx,yy};
      double time=CalTime(ele_tel,azi_tel,incoo0,indir,ImageCooXY,time0);
      if(time<0) continue;
      gr->SetPoint(gr->GetN(),ll,time);
   }
   gr->SetTitle(";Coordinate [degree];Time [ns]");
   gr->SetLineColor(4);
   gr->SetLineWidth(3);
   graphlist.push_back(gr);
   if(DoDraw) gr->Draw("l");
   return gr;
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
   //printf("CC=%.2lf phi=%.2lf\n",CC/PI*180,phi/PI*180);

   int np=100;
   TGraph* gr=new TGraph();
   for(int ii=0;ii<np;ii++){
      double PHI=0.+(PHI_max-0.)/np*ii;
      double xx,yy;
      bool res=GetImageXYCoo(zenith,azimuth,incoo,indir,PHI,xx,yy);
      if(!res) continue;
      printf("line: %d,xy={%lf,%lf}\n",gr->GetN(),xx/PI*180,yy/PI*180);
      gr->SetPoint(gr->GetN(),xx/PI*180,yy/PI*180);
   }
   gr->SetLineColor(1);
   gr->SetLineWidth(4);
   gr->Draw("l");
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
   if(jdebug>0) printf("WFCTAEvent::DrawImageLine: zenith=%lf azimuth=%lf incoo={%lf,%lf,%lf} indir={%lf,%lf}\n",zenith/PI*180,azimuth/PI*180,incoo[0],incoo[1],incoo[2],indir[0]/PI*180,indir[1]/PI*180);
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
   int Liindex=RotateDB::GetLi((double)rabbittime);
   for(int ii=0;ii<NSIPM;ii++){
      double ImageX,ImageY;
      GetImageXYCoo(ii,ImageX,ImageY);
      image->AddBin(-ImageX-0.25,ImageY-0.25,-ImageX+0.25,ImageY+0.25);
      //printf("WFCTAEvent::Draw: SiPM=%d ImageX=%lf ImageY=%lf\n",ii,ImageX,ImageY);
   }
   int ncontent=0;
   double sum,esum;
   sum=GetTotalPe(esum,ncontent,iTel,type,DoClean);
   for(int ii=0;ii<iSiPM.size();ii++){
      if(DoClean) {if(!CleanImage(ii,iTel,true)) continue;}
      double content=GetContent(ii,iTel,type,true);
      image->SetBinContent(iSiPM.at(ii)+1,content>0?content:0);
      if(jdebug>8) printf("WFCTAEvent::Draw: SiPM=%d content=%lf satl=%d\n",iSiPM.at(ii),content,(int)(eSatL.at(ii)));
   }
   image->SetTitle(Form("iTel=%d iEvent=%d Li=%d time=%ld+%lf(%d-%02d-%02d %02d:%02d:%02d) Npe=%.0lf;X [degree];Y [degree]",iTel,iEvent,Liindex<0?Liindex:RotateDB::rotindex[Liindex],rabbitTime,rabbittime*20*1.0e-9,CommonTools::TimeFlag((int)rabbitTime,1),CommonTools::TimeFlag((int)rabbitTime,2),CommonTools::TimeFlag((int)rabbitTime,3),CommonTools::TimeFlag((int)rabbitTime,4),CommonTools::TimeFlag((int)rabbitTime,5),CommonTools::TimeFlag((int)rabbitTime,6),sum));
   if(type==-1){
      for(int ii=0;ii<NSIPM;ii++){
         double content=((ii%2)==1)?ii:0;
         image->SetBinContent(ii+1,content>0?content:0);
      }
   }
   else{
      if((type>=1&&type<=4)||(type>=7&&type<=8)||(type>=11&&type<=14)){
         if(ncontent<5){
            delete image;
            return 0;
         }
      }
   }
   image->Draw(opt);
   DrawFit();
   TGraph* gr=DrawImageLine(0);
   if(!tlist) tlist=new TList();
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
/*void WFCTAEvent::slaDtp2s(double xi, double eta, double raz, double decz, double &ra, double &dec ){
  double sdecz, cdecz, denom;

  sdecz = sin ( decz );
  cdecz = cos ( decz );
  denom = cdecz - eta * sdecz;
  ra = ( atan2 ( xi, denom ) + raz );
  dec = atan2 ( sdecz + eta * cdecz, sqrt ( xi * xi + denom * denom ) );
  //if(ra<0) ra=ra+2*PI;
  //printf("xi=%lf eta=%lf raz=%lf decz=%lf ra=%lf dec=%lf\n",xi/PI*180,eta/PI*180,raz/PI*180,decz/PI*180,ra/PI*180,dec/PI*180);
}*/
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
      GetImageXYCoo(ii,ImageX,ImageY);

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
         slaDtp2s(-xx/180*PI,yy/180*PI,raz,decz,&(phi0[i2]),&(theta0[i2]));
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
      if(content<=0) continue;
      if(type==0&&ii<ADC_Cut.size()) content=content>adccuttrigger?1.:0.;
      image->SetBinContent(iSiPM.at(ii)+1,content>0?content:0);
   }
   if(DoDraw) image->Draw(opt);
   image->GetXaxis()->SetRangeUser(xrange[0]-2.,xrange[1]+2.);
   image->GetYaxis()->SetRangeUser(yrange[0]-2.,yrange[1]+2.);
   if(!tlist) tlist=new TList();
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
      GetImageXYCoo(ii,ImageX,ImageY);

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
         slaDtp2s(-xx/180*PI,yy/180*PI,raz,decz,&phi0,&theta0);
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
      if(content<=0) continue;
      if(type==0&&ii<ADC_Cut.size()) content=content>adccuttrigger?1.:0.;
      image->SetBinContent(iSiPM.at(ii)+1,content>0?content:0);
   }
   if(DoDraw) image->Draw(opt);
   image->GetXaxis()->SetRangeUser(xrange[0]-2.,xrange[1]+2.);
   image->GetYaxis()->SetRangeUser(yrange[0]-2.,yrange[1]+2.);
   //image->GetXaxis()->SetRangeUser(-90.,90.);
   //image->GetYaxis()->SetRangeUser(-90.,90.);
   if(!tlist) tlist=new TList();
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
                GetImageXYCoo(isipm,ImageX,ImageY);
		double content=GetContent(ii,iTel,type,true);
		if(content<=0) continue;
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

	if(!tlist) tlist=new TList();
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
   int irot=RotateDB::GetLi(rabbittime);
   return (irot>=0);
}
void WFCTAEvent::CalInfo(double result[100]){
   double MAX=0;
   double TNum=0;
   double TSum=0;
   double TSum2=0;
   double Tvariance=0;
   const int ncut=9;
   double value[ncut];
   double Num[ncut];
   double Sum[ncut];
   double SumX[ncut];
   double SumY[ncut];
   double SumX2[ncut];
   double SumY2[ncut];
   double Xaverage[ncut];
   double Yaverage[ncut];
   double SXvariance[ncut];
   double SYvariance[ncut];

   int size=iSiPM.size();
   for(int ii=0;ii<size;ii++){
      double npe=GetContent(ii,0,3,true);
      double ntime=GetContent(ii,0,5,true);
      if(npe>MAX) MAX=npe;
      if(npe>0){
         TNum++;
         TSum+=ntime;
         TSum2+=ntime*ntime;
      }
   }
   if(TNum>0) Tvariance=sqrt(TSum2/TNum-pow(TSum/TNum,2));
   int np=0;
   result[np++]=MAX;
   result[np++]=ncut;
   for(int icut=0;icut<ncut;icut++){
      value[icut]=0.1*(ncut-icut)*MAX;
      Num[icut]=0;
      Sum[icut]=0;
      SumX[icut]=0;
      SumY[icut]=0;
      SumX2[icut]=0;
      SumY2[icut]=0;
      Xaverage[icut]=0;
      Yaverage[icut]=0;
      SXvariance[icut]=0;
      SYvariance[icut]=0;
      for(int ii=0;ii<size;ii++){
         double npe=GetContent(ii,0,3,true);
         if(npe>value[icut]){
            Num[icut]++;
            Sum[icut]+=npe;
            double ImageXi,ImageYi;
            GetImageXYCoo(iSiPM.at(ii),ImageXi,ImageYi);
            SumX[icut]+=ImageXi*npe;
            SumY[icut]+=ImageYi*npe;
            SumX2[icut]+=ImageXi*ImageXi*npe;
            SumY2[icut]+=ImageYi*ImageYi*npe;
         }
      }
      if(Sum[icut]>0){
         Xaverage[icut]=SumX[icut]/Sum[icut];
         Yaverage[icut]=SumY[icut]/Sum[icut];
         SXvariance[icut]=SumX2[icut]/Sum[icut]-pow(Xaverage[icut],2);
         SYvariance[icut]=SumY2[icut]/Sum[icut]-pow(Yaverage[icut],2);
      }
      result[np++]=value[icut];
      result[np++]=Num[icut];
      result[np++]=Xaverage[icut];
      result[np++]=Yaverage[icut];
      result[np++]=SXvariance[icut];
      result[np++]=SYvariance[icut];
   }
   result[np++]=Tvariance;
}
bool WFCTAEvent::IsNoise(double* pars) {
   double result[100];
   CalInfo(result);
   if(result[0]>400) return false;
   int p0=(int)(pars[0]+0.5);
   int p2=(int)(pars[2]+0.5);
   double p1=pars[1];
   double p3=pars[3];
   double num=result[2+p0*6+1];
   double xvar=sqrt(result[2+p2*6+4]);
   double yvar=sqrt(result[2+p2*6+5]);
   return (num>p1)&&(xvar>=p3&&yvar>=p3);
}
bool WFCTAEvent::IsCR(double* pars) {
   double result[100];
   CalInfo(result);
   if(result[0]>400) return false;
   int p0=(int)(pars[0]+0.5);
   int p2=(int)(pars[2]+0.5);
   double p1=pars[1];
   double p3=pars[3];
   double num=result[2+p0*6+1];
   double xvar=sqrt(result[2+p2*6+4]);
   double yvar=sqrt(result[2+p2*6+5]);
   return (num<p1)&&(xvar<=p3&&yvar<=p3);
}

/******************************************
 * read root file and do iamge processing *
 * ****************************************/
//set x and y, unit in degree
void WFCTAEvent::SetImage()
{
	double D_ConeOut=25.8 ;// mm
	//double interval = 0.36;//1.0; //mm, gaps between subclusters 
	//int Interx,Intery;
	double centerx, centery;
	centerx = 414.3;
	centery = 414.3;
	for(int k=0;k<1024;k++){
		int i = k/32;
		int j = k%32;
		//Intery = i/4;
		//Interx = j/4;
		if(i%2==0)
			ImageX[k] = ((j+0.5)*D_ConeOut - centerx);// + interval*Interx-centerx);
		if(i%2==1)
			ImageX[k] = ((j+1)*D_ConeOut - centerx);// + interval*Interx-centerx);
		ImageY[k] = ((PIX-i)*D_ConeOut - centery);// + interval*Intery-centery);
		//ImageY[k] = ((PIX-i)*D_ConeOut + interval*Intery-centery);
	}

	/*
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
	   */
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
//void WFCTAEvent::AdcToPe()
void WFCTAEvent::AdcToPe(float *deltag_20, float *correct_PreTemp, int isledevent)
{
	int isipm;
	double pe;
	double Ntotal = 360000;
	double theta;
	double T0=20.;
	for(int ii=0;ii<iSiPM.size();ii++){
		isipm = (int)iSiPM.at(ii);
		if(SatH.at(ii)==0){  pe = AdcH.at(ii)/9.98;}
		else              {  pe = (AdcL.at(ii)*22)/9.98;}
		if(pe<0)          {  pe = 0;}
		else if(pe>Ntotal){  pe = Ntotal;}
		else              {  pe = -Ntotal*log(1-pe/Ntotal);}

		if(isledevent){
			theta = pow(cos(sqrt(ImageX[isipm]*ImageX[isipm]+ImageY[isipm]*ImageY[isipm])/2870),4);
			pe = pe/theta;
		}
		//if(iEvent%1000==0)
		//  printf("%d deltag_20:%f,coTemp:%f,coTemp-20:%f,theta:%lf\n",isipm,*(deltag_20+isipm),*(correct_PreTemp+isipm),*(correct_PreTemp+isipm)-T0,theta);
		pe = pe/(1 + (*(deltag_20+isipm))*(*(correct_PreTemp+isipm)-T0));

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
	double MAXDIST=30.4;//0.6;
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
	Dslope = a;
	Dintercept = b;
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
