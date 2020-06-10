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
#include "TelGeoFit.h"

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
int WFCTAEvent::ndiv=5;
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
               if(CommonTools::HArrival[ii]) CommonTools::HArrival[ii]->Fill(timeref,mcevent.ArrivalCount[itel][ii][jj]);
               //if(peaktime>=0) printf("itel=%d pmt=%d jj=%d time={%ld,%ld} maxcount=%lf\n",itel,ii,jj,mcevent.ArrivalTime[itel][jj],mcevent.ArrivalTime[itel][peaktime],maxcount);
            }
         }
      }
   }
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
         double baseh=ii<LaserBaseH.size()?LaserBaseH.at(ii):0;
         double wins=ii<winsum.size()?winsum.at(ii):0;
         if(ii<LaserAdcH.size()){
            content=(LaserAdcH.at(ii)>6000||(wins+baseh*4)>9000)?1.:-1.;
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
      bool ismc=CheckMC();
      if(((CalibType&0x3)!=0&&content>0)&&(!IsFit)&&(!ismc)){
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
   if(itel<1||itel>=NCTMax) itel=iTel;
   int telindex=WFTelescopeArray::GetHead()->GetTelescope(itel);
   if(telindex<0) return 0;
   bool ismc=CheckMC();
   double econtent=10.;
   int size=iSiPM.size();
   int start=IsIndex?isipm:0;
   int end=IsIndex?isipm:(size-1);
   if(IsIndex&&(isipm<0||isipm>=size)) return econtent;
   for(int ii=start;ii<=end;ii++){
      if((!IsIndex)&&(iSiPM.at(ii)!=isipm)) continue;
      if(ismc) econtent=mcevent.eTubeSignal[telindex][iSiPM.at(ii)];
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
      double ImageXi=0,ImageYi=0;
      TelGeoFit::GetImageXYCoo(iSiPM.at(ii),ImageXi,ImageYi,-1,true);
      if(jdebug>10) printf("WFCTAEvent::CleanImage: iTel=%d isipm=%d content=%.1lf trig1=%.1lf XYCoo={%.2lf,%.2lf}\n",iTel,isipm0,content,trig1,ImageXi,ImageYi);
      if(content<trig1) continue;
      int nneigh=0;
      for(int jj=0;jj<size;jj++){
         if(jj==ii) continue;
         double ImageXj=0,ImageYj=0;
         TelGeoFit::GetImageXYCoo(iSiPM.at(jj),ImageXj,ImageYj,-1,true);
         double dist=sqrt(pow(ImageXi-ImageXj,2)+pow(ImageYi-ImageYj,2));
         if(dist>0.6) continue;
         bool sat1=(GetContent(jj,itel,9,true,IsFit)>0.5);
         double contentj0=GetContent(jj,itel,0,true,IsFit);
         double contentj=GetContent(jj,itel,11,true,IsFit);
         //if(ADC_Cut.size()>jj) {if(contentj0<trig0) continue;}
         //else{if(contentj<trig1) continue;}
         if(contentj<trig1) continue;
         if(jdebug>10) printf("WFCTAEvent::CleanImage: iTel=%d isipm=%d isipm2=%d pass trigger (content=%.1lf XYCoo={%.2lf,%.2lf})\n",iTel,isipm0,iSiPM.at(jj),contentj,ImageXj,ImageYj);
         nneigh++;
      }
      if(jdebug>10) printf("WFCTAEvent::CleanImage: iTel=%d isipm=%d nneigh=%d nfiretrigger=%d\n",iTel,isipm0,nneigh,nfiretrigger);
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
      double dcell=TelGeoFit::GetImageXYCoo(isipm,ImageX,ImageY,WFTelescope::FOCUS,false);
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
   if(itel<1||itel>WFTelescopeArray::CTNumber) itel=iTel;
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
      TelGeoFit::GetImageXYCoo(isipm0,ImageX,ImageY,-1,false);

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

TH2F* WFCTAEvent::GetRotImage(int itel,int type,double CC,double phi){
   if(phi<0){
      if(!minimizer) return 0;
      CC=minimizer->X()[2];
      phi=minimizer->X()[3];
   }

   bool IsTime=((type>=5&&type<=6)||type==15);

   static int ihist=0;
   TH2F* h2d=new TH2F(Form("rot_image+hist%d",ihist++),";Long Axis [degree];Short Axis [degree]",13,-13.,13.,32,-8.,8.);
   TH2F* hweight=0;
   if(IsTime){
      hweight=new TH2F(Form("hweight"),";Long Axis [degree];Short Axis [degree]",13,-13.,13.,32,-8.,8.);
   }

   int size=iSiPM.size();
   for(int ii=0;ii<size;ii++){
      int isipm0=iSiPM.at(ii);
      if(!CleanImage(ii,itel,true)) continue;
      double content=GetContent(ii,itel,type,true);
      double econtent=GetContentError(ii,itel,type,true);
      double weight=IsTime?GetContent(ii,itel,12,true):1.;
      if(weight<=0) continue;

      double ImageX=0,ImageY=0;
      double width=TelGeoFit::GetImageXYCoo(isipm0,ImageX,ImageY,-1,true);
      double xrange[2]={ImageX-width/2,ImageX+width/2};
      double yrange[2]={ImageY-width/2,ImageY+width/2};

      for(int ibin1=0;ibin1<ndiv;ibin1++){
         double xx=xrange[0]+(xrange[1]-xrange[0])/ndiv*(ibin1+0.5);
         for(int ibin2=0;ibin2<ndiv;ibin2++){
            double yy=yrange[0]+(yrange[1]-yrange[0])/ndiv*(ibin2+0.5);
            double ImageX2=xx*cos(phi)+yy*sin(phi); //long axis
            double ImageY2=yy*cos(phi)-xx*sin(phi)-CC/PI*180; //short axis
            int ibinx=h2d->GetXaxis()->FindBin(ImageX2);
            int ibiny=h2d->GetYaxis()->FindBin(ImageY2);
            h2d->SetBinContent(ibinx,ibiny,h2d->GetBinContent(ibinx,ibiny)+content*weight/(ndiv*ndiv));
            h2d->SetBinError(ibinx,ibiny,sqrt(pow(h2d->GetBinError(ibinx,ibiny),2)+pow(econtent,2)*weight/(ndiv*ndiv)));
            if(IsTime){
               hweight->SetBinContent(ibinx,ibiny,hweight->GetBinContent(ibinx,ibiny)+weight/(ndiv*ndiv));
            }
         }
      }
   }
   if(IsTime){
      h2d->Divide(hweight);
      delete hweight;
   }
   return h2d;
}
TH1F* WFCTAEvent::GetWidth(int itel,int type,double CC,double phi){
   bool IsTime=((type>=5&&type<=6)||type==15);

   static int ihist=0;
   TH1F* hh=new TH1F(Form("longdis_%swidth_h%d",IsTime?"time_":"npe_",ihist++),Form(";long axis [degree];%s",IsTime?"Time [ns]":"Npe [pe]"),13,-13,13);

   TH2F* h2d=GetRotImage(itel,type,CC,phi);
   for(int ibin=1;ibin<=h2d->GetNbinsX();ibin++){
      TH1D* hbuff=h2d->ProjectionY(Form("hbuff_bin%d",ibin),ibin,ibin);
      if(IsTime){
         double time_mean=0;
         double time_rms=0;
         double etime_rms=0;
         int nn=0;
         for(int ii=1;ii<=hbuff->GetNbinsX();ii++){
            double timei=hbuff->GetBinContent(ii);
            double etimei=hbuff->GetBinError(ii);
            if(timei<=0) continue;
            time_mean+=timei;
            time_rms+=pow(timei,2);
            etime_rms+=pow(timei*etimei,2);
            nn++;
         }
         if(nn>0){
            time_mean/=nn;
            time_rms=sqrt(time_rms/nn-time_mean*time_mean);
            etime_rms=sqrt(etime_rms)/nn/time_rms;
            hh->SetBinContent(ibin,time_rms);
            hh->SetBinError(ibin,etime_rms);
         }
      }
      else{
         double sigmi=-1;
         if(hbuff->Integral()>50){
            double allcont=hbuff->Integral();
            int ibin0=hbuff->GetXaxis()->FindBin(0.);
            double acccont=0;
            for(int ibin=0;ibin<(hbuff->GetNbinsX()/2)+2;ibin++){
               acccont+=hbuff->GetBinContent(ibin0+ibin);
               if(ibin>0) acccont+=hbuff->GetBinContent(ibin0-ibin);
               if(fabs(acccont/allcont)>0.95){
                  sigmi=hbuff->GetXaxis()->GetBinCenter(ibin0+ibin)-hbuff->GetXaxis()->GetBinCenter(ibin0);
                  break;
               }
            }
         }
         if(sigmi>0){
            hh->SetBinContent(ibin,sigmi);
            hh->SetBinError(ibin,0);
         }
      }
      delete hbuff;
   }
   delete h2d;

   return hh;
}
TH1F* WFCTAEvent::GetDistribution(bool IsLong,int itel,int type,bool CleanEdge){
   itel=-1;
   if(itel<0||itel>=WFTelescopeArray::CTNumber) itel=iTel;
   //if(!minimizer) return 0;
   //double CC=minimizer->X()[2];
   //double phi=minimizer->X()[3];
   double CC,phi;
   if(!GetCCphi(CC,phi)) return 0;
   double XY1[4],XY2[4];
   bool inside=TelGeoFit::GetRange(CC,phi,XY1,XY2);
   if(!inside) return 0;
   double margin=24.;
   const int nbin=48;

   TH1F* hh=0;
   TH2F* h2d=GetRotImage(itel,type,CC,phi);
   if(!h2d) return hh;
   static int ihist=0;
   if(IsLong){
      hh=(TH1F*)h2d->ProjectionX(Form("%sdis_%s_h%d",IsLong?"long":"short",((type>=5&&type<=6)||type==15)?"time":"npe",ihist++),0,h2d->GetNbinsY()+1);
      if(CleanEdge){
         double LRange[2]={-30,30};
         double npe[MAXPMT];
         for(int ii=0;ii<MAXPMT;ii++) npe[ii]=0;
         for(int ii=0;ii<iSiPM.size();ii++){
            if(!CleanImage(ii,itel,true)) continue;
            npe[iSiPM.at(ii)]=GetContent(ii,itel,12,true);
         }
         double width=TelGeoFit::GetWidth(npe,CC,phi);
         bool longrange=TelGeoFit::GetLongRange(CC,phi,width,LRange);
         if(!longrange) {delete hh; hh=0;}
         LRange[0]*=180/PI;
         LRange[1]*=180/PI;
         //clean edge
         for(int ibin=1;ibin<=hh->GetNbinsX();ibin++){
            double xcenter=hh->GetXaxis()->GetBinCenter(ibin);
            if(xcenter<LRange[0]||xcenter>LRange[1]){
               hh->SetBinContent(ibin,0);
               hh->SetBinError(ibin,0);
            }
         }
      }
   }
   else{
      hh=(TH1F*)h2d->ProjectionY(Form("%sdis_%s_h%d",IsLong?"long":"short",((type>=5&&type<=6)||type==15)?"time":"npe",ihist++),0,h2d->GetNbinsX()+1);
   }
   if(hh) hh->SetTitle(Form(";%s axis [degree];%s",IsLong?"Long":"Short",((type>=5&&type<=6)||type==15)?"Time [ns]":"Npe [pe]"));
   delete h2d;
   return hh;
}
TH1F* WFCTAEvent::GetScatterAngle(int itel,int type,bool CleanEdge,bool UseFit){
   itel=0;
   if(itel<1||itel>WFTelescopeArray::CTNumber) itel=iTel;

   WFTelescopeArray* pta=WFTelescopeArray::GetHead();
   int telindex=pta->GetTelescope(itel);
   if(telindex<0){
      cerr<<"WFCTAEvent::GetScatterAngle: No Telescope Information Loaded."<<endl;
      return 0;
   }
   WFTelescope* pt=pta->pct[telindex];
   double zenith=pt->TelZ_;
   double azimuth=pt->TelA_;
   double telcoo[3];
   telcoo[0]=pt->Telx_;
   telcoo[1]=pt->Tely_;
   telcoo[2]=pt->Telz_+WFTelescopeArray::lhaaso_coo[2];

   int irot=RotateDB::GetHead()->GetLi(rabbittime);
   if(irot<0){
      cerr<<"WFCTAEvent::GetScatterAngle: Not a Laser Event."<<endl;
      return 0;
   }
   double rotcoo[3];
   for(int ii=0;ii<3;ii++) rotcoo[ii]=TelGeoFit::GetRotPos(ii,RotateDB::rotindex[irot])+(ii==2?WFTelescopeArray::lhaaso_coo[2]:0);

   int index=-1;
   bool ismc=CheckMC();
   double ele_in,azi_in;
   if(ismc){
      double laserdir[2]={90-laserevent.LaserDir[0],laserevent.LaserDir[1]};
      index=RotateDB::IsFineAngle(laserdir[0],laserdir[1],RotateDB::rotindex[irot],itel);
      if(index>0) index+=RotateDB::rotindex[irot]*100;
      ele_in=laserdir[0];
      azi_in=laserdir[1];
      index=1;
      //printf("MC: ele_in=%.2lf azi=%.2lf index=%d\n",laserdir[0],laserdir[1],index);
   }
   else{
      index=RotateDB::GetHead()->GetEleAzi(this);
      ele_in=RotateDB::GetHead()->GetElevation();
      azi_in=RotateDB::GetHead()->GetAzimuth();
   }
   if(index<=0){
      cerr<<"WFCTAEvent::GetScatterAngle: Not Rotate Log Found."<<endl;
      return 0;
   }
   double rotdir[2];
   TelGeoFit::CalDir_out(ele_in/180*PI,azi_in/180*PI,RotateDB::rotindex[irot],rotdir[0],rotdir[1]);
   rotdir[0]=PI/2-rotdir[0];
   double CC,phi;
   if(UseFit){
      if(!minimizer){
         cerr<<"WFCTAEvent::GetScatterAngle: No Fit Exist."<<endl;
         return 0;
      }
      CC=minimizer->X()[2];
      phi=minimizer->X()[3];
   }
   else{
      TelGeoFit::GetCCphi(zenith,azimuth,telcoo,rotcoo,rotdir,CC,phi);
   }
   printf("zen=%.2lf azi=%.2lf telcoo={%.1lf,%.1lf,%.1lf} rotcoo={%.1lf,%.1lf,%.1lf} rotdir={%.2lf,%.2lf} CC=%.2lf phi=%.2lf\n",zenith/PI*180,azimuth/PI*180,telcoo[0],telcoo[1],telcoo[2],rotcoo[0],rotcoo[1],rotcoo[2],rotdir[0]/PI*180,rotdir[1]/PI*180,CC/PI*180,phi/PI*180);

   double LRange[2]={-30,30};
   if(CleanEdge){
      double npe[MAXPMT];
      for(int ii=0;ii<MAXPMT;ii++) npe[ii]=0;
      for(int ii=0;ii<iSiPM.size();ii++){
         int isipm=iSiPM.at(ii);
         npe[isipm]=GetContent(ii,itel,type,true);
      }
      double width=TelGeoFit::GetWidth(npe,CC,phi);
      bool longrange=TelGeoFit::GetLongRange(CC,phi,width,LRange);
      LRange[0]*=180/PI;
      LRange[1]*=180/PI;
      printf("width=%.2lf LRange={%.2lf,%.2lf} return=%d\n",width/PI*180,LRange[0],LRange[1],longrange);
      if(!longrange) return 0;
   }

   static int ihist=0;
   TH1F* hh=new TH1F(Form("scatter_longcoo_%d",ihist++),Form(";cos(scatter angle);%s",((type>=5&&type<=6)||type==15)?"Time [ns]":"Npe [pe]"),100,-1,1);

   int size=iSiPM.size();
   for(int ii=0;ii<size;ii++){
      int isipm=iSiPM.at(ii);
      if(!CleanImage(ii,itel,true)) continue;
      double content=GetContent(ii,itel,type,true);
      double econtent=GetContentError(ii,itel,type,true);
      if(content<=0) continue;

      double ImageX,ImageY;
      double width=TelGeoFit::GetImageXYCoo(isipm,ImageX,ImageY,-1,false);
      double xrange[2]={ImageX-width/2,ImageX+width/2};
      double yrange[2]={ImageY-width/2,ImageY+width/2};

      for(int ibin1=0;ibin1<ndiv;ibin1++){
         double xx=xrange[0]+(xrange[1]-xrange[0])/ndiv*(ibin1+0.5);
         for(int ibin2=0;ibin2<ndiv;ibin2++){
            double yy=yrange[0]+(yrange[1]-yrange[0])/ndiv*(ibin2+0.5);
            double longaxis=xx*cos(phi)+yy*sin(phi); //long axis
            double ycoo=yy*cos(phi)-xx*sin(phi)-CC; //short axis
            if(CleanEdge) {if(longaxis<=LRange[0]||longaxis>=LRange[1]) continue;}
            double scatter_angle=TelGeoFit::GetScatAngleFromLongcoo(zenith,azimuth,telcoo,rotcoo,rotdir,longaxis);
            int ibin=hh->GetXaxis()->FindBin(cos(scatter_angle));
            hh->SetBinContent(ibin,hh->GetBinContent(ibin)+content/(ndiv*ndiv));
            hh->SetBinError(ibin,sqrt(pow(hh->GetBinError(ibin),2)+pow(econtent,2)/(ndiv*ndiv)));
         }
      }
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

TGraph* WFCTAEvent::DrawImageLine(double zenith,double azimuth,double telcoo[3],double incoo[3],double indir[2]){
   double planephi,nz;
   double telzcoo[3];
   TelGeoFit::Getnz(telcoo,incoo,indir,planephi,nz,telzcoo);
   double norm=1.;
   double xdir[3]={telzcoo[0]-telcoo[0],telzcoo[1]-telcoo[1],telzcoo[2]-telcoo[2]};
   norm=sqrt(pow(xdir[0],2)+pow(xdir[1],2)+pow(xdir[2],2));
   for(int ii=0;ii<3;ii++) xdir[ii]/=norm;
   double dir0[3]={incoo[0]-telcoo[0],incoo[1]-telcoo[1],incoo[2]-telcoo[2]};
   norm=sqrt(pow(dir0[0],2)+pow(dir0[1],2)+pow(dir0[2],2));
   for(int ii=0;ii<3;ii++) dir0[ii]/=norm;
   double dir1[3]={sin(zenith)*cos(azimuth),sin(zenith)*sin(azimuth),cos(zenith)};
   double PHI_min=acos(dir0[0]*xdir[0]+dir0[1]*xdir[1]+dir0[2]*xdir[2]);
   if(incoo[2]-telzcoo[2]<0) PHI_min*=-1;
   double PHI_max=acos(dir1[0]*xdir[0]+dir1[1]*xdir[1]+dir1[2]*xdir[2]);

   double CC,phi;
   TelGeoFit::GetCCphi(zenith,azimuth,telcoo,incoo,indir,CC,phi);
   //printf("telcoo={%.1lf,%.1lf,%.1lf} incoo={%.1lf,%.1lf,%.1lf}\n",telcoo[0],telcoo[1],telcoo[2],incoo[0],incoo[1],incoo[2]);
   //printf("CC=%.2lf phi=%.2lf PHI_min=%.2lf PHI_max=%.2lf\n",CC/PI*180,phi/PI*180,PHI_min/PI*180,PHI_max/PI*180);

   int np=100;
   TGraph* gr=new TGraph();
   for(int ii=0;ii<np;ii++){
      double PHI=PHI_min+(PHI_max-PHI_min)/np*(ii+0.5);
      double longcoo=TelGeoFit::GetLongCoo(zenith,azimuth,planephi,nz,PHI);
      double xx=TelGeoFit::GetImageCoo(CC,phi,longcoo,0);
      double yy=TelGeoFit::GetImageCoo(CC,phi,longcoo,1);
      //printf("line: %d,xy={%lf,%lf}\n",gr->GetN(),xx/PI*180,yy/PI*180);
      gr->SetPoint(gr->GetN(),xx/PI*180,yy/PI*180);
   }
   gr->SetLineColor(2);
   gr->SetLineWidth(2);
   gr->Draw("l");
   return gr;
}
bool WFCTAEvent::GetLaserPosDir(double coo[3],double &zen,double &azi,bool TrueDir){
   int irot=RotateDB::GetLi(rabbittime);
   if(irot<0) return false;
   bool ismc=CheckMC();
   if(ismc){
      for(int ii=0;ii<3;ii++) coo[ii]=laserevent.LaserCoo[ii];
      zen=laserevent.LaserDir[0]/180.*PI;
      azi=laserevent.LaserDir[1]/180.*PI;
   }
   else{
      for(int ii=0;ii<3;ii++) coo[ii]=TelGeoFit::GetRotPos(ii,RotateDB::rotindex[irot])+(ii==2?WFTelescopeArray::lhaaso_coo[2]:0);
      int index=RotateDB::GetHead()->GetEleAzi(this);
      zen=PI/2-RotateDB::GetHead()->GetElevation()/180*PI;
      azi=RotateDB::GetHead()->GetAzimuth()/180*PI;
      if(RotateDB::GetHead()->GetElevation()==0) return false;
   }
   if(!TrueDir) return true;
   else{
      double indir[2]={zen,azi};
      TelGeoFit::CalDir_out(PI/2-indir[0],indir[1],RotateDB::rotindex[irot],zen,azi);
      zen=PI/2-zen;
      return true;
   }
}
bool WFCTAEvent::GetCCphi(double &CC,double &phi){
   WFTelescopeArray* pct=WFTelescopeArray::GetHead();
   if(!pct) return false;
   int telindex=pct->GetTelescope(iTel);
   if(telindex<0) return false;
   WFTelescope* pt=pct->pct[telindex];
   if(!pt) return false;
   double zenith=pt->TelZ_;
   double azimuth=pt->TelA_;
   double telcoo[3]={pt->Telx_,pt->Tely_,pt->Telz_+WFTelescopeArray::lhaaso_coo[2]};
   double incoo[3];
   double indir[2];
   if(!GetLaserPosDir(incoo,indir[0],indir[1],true)) return false;
   TelGeoFit::GetCCphi(zenith,azimuth,telcoo,incoo,indir,CC,phi);
   return true;
}
TGraph* WFCTAEvent::DrawImageLine(int itel){
   if(itel<1||itel>NCTMax) itel=iTel;
   WFTelescopeArray* pct=WFTelescopeArray::GetHead();
   if(!pct) return 0;
   int telindex=pct->GetTelescope(itel);
   if(telindex<0) return 0;
   WFTelescope* pt=pct->pct[telindex];
   if(!pt) return 0;
   double zenith=pt->TelZ_;
   double azimuth=pt->TelA_;
   double telcoo[3]={pt->Telx_,pt->Tely_,pt->Telz_+WFTelescopeArray::lhaaso_coo[2]};

   int irot=RotateDB::GetLi(rabbittime);
   if(irot<0) return 0;

   double incoo[3];
   double indir[2];
   if(!GetLaserPosDir(incoo,indir[0],indir[1],true)) return 0;
   if(jdebug>0) printf("WFCTAEvent::DrawImageLine: zenith=%lf azimuth=%lf incoo={%lf,%lf,%lf} indir={%lf,%lf}\n",zenith/PI*180,azimuth/PI*180,incoo[0],incoo[1],incoo[2],indir[0]/PI*180,indir[1]/PI*180);
   return DrawImageLine(zenith,azimuth,telcoo,incoo,indir);
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
   gDrawErr->SetFillColor(14);
   if(DoDraw) gDrawErr->Draw("F");
   gDraw->SetLineColor(1);
   gDraw->SetLineWidth(2);
   if(DoDraw) gDraw->Draw("l");
   printf("WFCTAEvent::DrawFit: phi={%lf,%lf} cc={%lf,%lf}\n",phi/PI*180,sqrt(Cov[2])/PI*180,CC/PI*180,eCC/PI*180);
}
TH2Poly* WFCTAEvent::Draw(int type,const char* opt,bool DoClean,int index,double threshold){
   TH2Poly* image=new TH2Poly();
   image->SetName("DrawPlot");
   int Liindex=RotateDB::GetLi((double)rabbittime);
   for(int ii=0;ii<NSIPM;ii++){
      double ImageX,ImageY;
      TelGeoFit::GetImageXYCoo(ii,ImageX,ImageY,-1,true);
      image->AddBin(ImageX-0.25,ImageY-0.25,ImageX+0.25,ImageY+0.25);
      //printf("WFCTAEvent::Draw: SiPM=%d ImageX=%lf ImageY=%lf\n",ii,ImageX,ImageY);
   }
   int ncontent=0;
   double sum,esum;
   sum=GetTotalPe(esum,ncontent,iTel,type,DoClean);
   for(int ii=0;ii<iSiPM.size();ii++){
      //if(jdebug>8||true) printf("WFCTAEvent::Draw: iTel=%d SiPM=%d size=%d Adc=%.2lf\n",iTel,iSiPM.at(ii),iSiPM.size(),LaserAdcH.at(ii));
      if(DoClean) {if(!CleanImage(ii,iTel,true)) continue;}
      double content=GetContent(ii,iTel,type,true);
      image->SetBinContent(iSiPM.at(ii)+1,content>0?content:0);
      if(jdebug>8) printf("WFCTAEvent::Draw: iTel=%d SiPM=%d content=%lf satl=%d\n",iTel,iSiPM.at(ii),content,(int)(eSatL.at(ii)));
   }
   image->SetTitle(Form("iTel=%d iEvent=%d Li=%d time=%ld+%lf(%d-%02d-%02d %02d:%02d:%02d) Npe=%.0lf;X [degree];Y [degree]",iTel,iEvent,Liindex<0?Liindex:(index==-1?RotateDB::rotindex[Liindex]:index),rabbitTime,rabbittime*20*1.0e-9,CommonTools::TimeFlag((int)rabbitTime,1),CommonTools::TimeFlag((int)rabbitTime,2),CommonTools::TimeFlag((int)rabbitTime,3),CommonTools::TimeFlag((int)rabbitTime,4),CommonTools::TimeFlag((int)rabbitTime,5),CommonTools::TimeFlag((int)rabbitTime,6),sum));
   if(type==-1){
      for(int ii=0;ii<NSIPM;ii++){
         double content=((ii%2)==1)?ii:0;
         image->SetBinContent(ii+1,content>0?content:0);
      }
   }
   else{
      if((type>=1&&type<=4)||(type>=7&&type<=8)||(type>=11&&type<=14)){
         if(ncontent<1){
            delete image;
            return 0;
         }
      }
   }
   image->Draw(opt);
   DrawFit();
   if(DoDraw) DrawImageLine();
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
      TelGeoFit::GetImageXYCoo(ii,ImageX,ImageY,-1,true);

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
      TelGeoFit::GetImageXYCoo(ii,ImageX,ImageY,-1,true);

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
                TelGeoFit::GetImageXYCoo(isipm,ImageX,ImageY,-1,true);
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
            TelGeoFit::GetImageXYCoo(iSiPM.at(ii),ImageXi,ImageYi,-1,true);
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
