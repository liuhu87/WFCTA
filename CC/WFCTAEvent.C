#include <stdlib.h>
#include <cmath>
#include <iostream>
#include "WFCTAEvent.h"
#include "TMarker3DBox.h"
#include "TAxis3D.h"
#include <TCanvas.h>
#include <TView3D.h>
#include <TSystem.h>
#include "TF1.h"

using namespace std;

ClassImp(WFCTAEvent);

const int WFCTAEvent::cleanPix=3;
WFCTAEvent* WFCTAEvent::_Head=0;
TTree* WFCTAEvent::_Tree=0;
const char* WFCTAEvent::_Name="WFCTAEvent";
TBranch* WFCTAEvent::bAll=0;
TBranch* WFCTAEvent::bledevent=0;
int WFCTAEvent::jdebug=0;
double WFCTAEvent::adccuttrigger=350.;
double WFCTAEvent::npetrigger=45.;
int WFCTAEvent::nfiretrigger=3;
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
	packCheck.resize(MAXPMT);
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
   //packCheck=-1;
   packCheck.clear();
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


   ledevent.Init();
}
void WFCTAEvent::EventInitial()
{
	iTel=-1;
	merge_size=-1;
	iEvent=-1;
	eEvent=-1;
	rabbitTime=0;
	rabbittime=0;
	big_pack_lenth=-1;
	n_fired=-1;
	n_Channel=-1;
	//packCheck=-1;
	packCheck.clear();
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

	ledevent.Reset();
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
		//tree->SetBranchStatus("ledevent",false);
	}
}
void WFCTAEvent::GetBranch(TTree *fChain){
	bAll=fChain->GetBranch(BranchName());
	bledevent=fChain->GetBranch("ledevent");
}
bool WFCTAEvent::GetAllContents(int _Entry){
	EventInitial();
	bool exist=true;
	if(!bAll) exist=false;
	//if(!bledevent) exist=false;
	if(!exist) return false;
	int ncount=bAll->GetEntry(_Entry);
	if(bledevent) bledevent->GetEntry(_Entry);
	return ncount>0;
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

/******************************************
 * read root file and do iamge processing *
 * ****************************************/
//set x and y, unit in degree
void WFCTAEvent::SetImage()
{
	int PIX=32;
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
void WFCTAEvent::PrelimImageClean(double cut, int laserEvent, int ledEvent)
{
	if(ledEvent){
		for(int ii=0;ii<RawImagePe.size();ii++){
			FullImagePe.push_back(RawImagePe.at(ii));
			FullImageSiPM.push_back(RawImageSiPM.at(ii));
			FullImageX.push_back(RawImageX.at(ii));
			FullImageY.push_back(RawImageY.at(ii));
		}
	}
	else if(laserEvent){
		for(int ii=0;ii<RawImagePe.size();ii++){
			if(RawImagePe.at(ii)<cut){continue;}
			FullImagePe.push_back(RawImagePe.at(ii));
			FullImageSiPM.push_back(RawImageSiPM.at(ii));
			FullImageX.push_back(RawImageX.at(ii));
			FullImageY.push_back(RawImageY.at(ii));
		}
	}
	else{
		cut =20;
		int Maxpro_Peakpos_H=-1;
		int Maxpro_Peakpos_H_Times=-1;
		map<short,int> maxpro_peakposh;
		map<short,int>::iterator maxpro_peakposh_iter;
		maxpro_peakposh.clear();
		for(int ii=0;ii<RawImagePe.size();ii++){
			if(RawImagePe.at(ii)<cut){continue;}
			maxpro_peakposh[PeakPosH.at(ii)]++;
		}
		for(maxpro_peakposh_iter=maxpro_peakposh.begin(); maxpro_peakposh_iter!=maxpro_peakposh.end(); maxpro_peakposh_iter++){
			if(Maxpro_Peakpos_H_Times<maxpro_peakposh_iter->second){
				Maxpro_Peakpos_H_Times = maxpro_peakposh_iter->second;
				Maxpro_Peakpos_H = maxpro_peakposh_iter->first;
			}
		}
		for(int ii=0;ii<RawImagePe.size();ii++){
			if(PeakPosH.at(ii)<Maxpro_Peakpos_H-2 || PeakPosH.at(ii)>Maxpro_Peakpos_H+2){/*printf("%d\n",RawImagePe.at(ii));*/continue;}
			if(RawImagePe.at(ii)<cut){continue;}
			FullImagePe.push_back(RawImagePe.at(ii));
			FullImageSiPM.push_back(RawImageSiPM.at(ii));
			FullImageX.push_back(RawImageX.at(ii));
			FullImageY.push_back(RawImageY.at(ii));
		}
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
int WFCTAEvent::CalcHillas(int ledEvent)
{
	DNpix = 0;
	DSize = 0;
	DMeanX = 0;
	DMeanY = 0;
	Dslope = 0;
	Dintercept = 0;
	DLength = 0;
	DWidth = 0;
	if(ledEvent){
		for(int ii=0;ii<FullImagePe.size();ii++){
			DNpix++;
			DSize += FullImagePe.at(ii);
		}
	}
	else{
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
	}
	return 0;
}

//clean image
void WFCTAEvent::GetCleanImage(int LedEvent)
{
	if(LedEvent){
		for(int ii=0;ii<FullImagePe.size();ii++){
			CleanImagePe.push_back(FullImagePe.at(ii));
			CleanImageSiPM.push_back(FullImageSiPM.at(ii));
			CleanImageX.push_back(FullImageX.at(ii));
			CleanImageY.push_back(FullImageY.at(ii));
		}
	}
	else{
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
}
