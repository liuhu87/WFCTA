#ifndef WFCTAEVENT_H
#define WFCTAEVENT_H

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include "TObject.h"
#include "TObjArray.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH2Poly.h"
#include "common.h"

#include "WFCTAMCEvent.h"
#include "WFCTALedEvent.h"
#include "WFCTALaserEvent.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

using namespace std;
class LHChain;
class WFCTAEvent : public TSelector
{
	protected:
		static const int cleanPix;
		static WFCTAEvent* _Head;
		static TTree* _Tree;
		static LHChain* _Chain;
		static const char * _Name;
		static TBranch* bAll;
		static TBranch* bmcevent;
		static TBranch* bledevent;
		static TBranch* blaserevent;
		static int _Entry;

   public:
   static int jdebug;
   static bool DoDraw;
   static double adccuttrigger;
   static double npetrigger;
   static int nfiretrigger;
   static int CalibType;

   TList* tlist; //!
   vector<TGraph*> graphlist; //!

	public:
		short iTel;
		short merge_size;
		long iEvent;
		long eEvent;
		long rabbitTime;
		double rabbittime;
		int big_pack_lenth;
		short n_fired;
		short n_Channel;

		vector<short> packCheck;
		vector<long> eevent;
		vector<short> zipmod;
		vector<short> iSiPM;
		vector<float> winsum;
		vector<float> ADC_Cut;  //!
		vector<float> eBaseH;  //!
		vector<float> eBaseL;  //!
		vector<float> eAdcH;  //!
		vector<float> eAdcL;  //!
		vector<bool> eSatH;
		vector<bool> eSatL;
		vector<float> BaseH;
		vector<float> BaseL;
		vector<float> LaserBaseH;
		vector<float> LaserBaseL;
		vector<float> AdcH;
		vector<float> AdcL;
		vector<float> LaserAdcH;
		vector<float> LaserAdcL;
		vector<bool> SatH;
		vector<bool> SatL;
		vector<short> Single_Threshold;  //!
		vector<short> Record_Threshold;  //!
		vector<char> peak;  //!
		vector<short> PeakPosH;
		vector<short> PeakPosL;
		vector<int> PeakAmH;
		vector<int> PeakAmL;
		vector<bool> gain_marker;  //!
		vector<bool> Over_Single_Marker;
		vector<bool> Over_Record_Marker;

		int Npoint[28];  //!
		int pulsehigh[1024][28]; //!
		int pulselow[1024][28];  //!
		double ImageX[1024];  //!
		double ImageY[1024];  //!

		WFCTAMCEvent mcevent;
		WFCTALedEvent ledevent;
		WFCTALaserEvent laserevent;

                /// for the fitting of the image
                TGraph* gDraw; //!
                TGraph* gDrawErr; //!
                ROOT::Math::Minimizer* minimizer; //!
                /// for the time vs long axis
                TH1F* htlong; //!

                ///Event Type: 1:laser 2:led 3:CRs
                int Type; //!

	public:
		double mjd;  //!
		int year;  //!
		int month;  //!
		int day;  //!
		int hour;  //!
		int minite;  //!
		int second;  //!
		int DNpix;  //!
		double DSize;  //!
		double DMeanX;  //!
		double DMeanY;  //!
		double Dslope;  //!
		double Dintercept;  //!
		double DLength;  //!
		double DWidth;  //!
		vector<double> RawImagePe;  //!
		vector<int> RawImageSiPM;  //!
		vector<double> RawImageX;  //!
		vector<double> RawImageY;  //!
		vector<double> FullImagePe;  //!
		vector<int> FullImageSiPM;  //!
		vector<double> FullImageX;  //!
		vector<double> FullImageY;  //!
		vector<int> CleanImagePe;  //!
		vector<int> CleanImageSiPM;  //!
		vector<double> CleanImageX;  //!
		vector<double> CleanImageY;  //!

		vector<int> fNpixfriends;  //!

	public:
		WFCTAEvent();
		~WFCTAEvent();
		void Init();
		void EventInitial();
		void InitTree(TTree* tree);
		static const char *  BranchName() {return _Name;}
		static WFCTAEvent* & Head()  {return _Head;}
		static TTree* & Tree()  {return _Tree;}
                static void SetLHChain(LHChain* chain);
		void CreateBranch(TTree *tree, int branchSplit);
		void GetBranch(TTree *fChain);
		bool GetAllContents(int entry=-1);
                const char* GetFileName();
                bool CheckLaser();
                bool CheckMC();
                void CalculateDataVar(int itel=0);
                static double GetImageXYCoo(int isipm,double &ImageX,double &ImageY,double focus=-1,bool Isdegree=true);
                int GetMaxADCBin(int itel=0);
                int GetPeakADCBin(int isipm,int itel=0);
                int GetMaxTimeBin(int itel=0);
                int GetMinTimeBin(int itel=0);
                bool IsSaturated(int isipm,int itel=0,bool IsHigh=true,bool IsIndex=false);
                double GetContent(int isipm,int itel=0,int type=3,bool IsIndex=false,bool IsFit=false);
                double GetContentError(int isipm,int itel=0,int type=3,bool IsIndex=false,bool IsFit=false);
                double GetTotalPe(double &error,int &ncontent,int itel,int type,bool DoClean=true,bool IsFit=false);
                bool CleanImage(int isipm,int itel=0,bool IsIndex=false,bool IsFit=false);
                bool PassClean(int itel=0,int nthreshold=5);
                double Interface(const double* par);
                bool DoFit(int itel=0,int type=3,bool force=false);
                static bool GetCrossCoor(double x, double y, double CC, double phi, double &x0, double &y0);
                bool GetCrossCoor(double x, double y, double &x0, double &y0);
                static int IsInside(double xx,double yy);
                static bool GetRange(double CC,double phi,double XY1[4],double XY2[4]);
                TH1F* GetDistribution(bool IsLong,int itel=0,int type=3,bool IsWidth=false,bool CleanEdge=true);
                TH1F* CorrTimeDist(TH1F* hist,int IsLaser=0,int IsMC=0);

                int GetSign(bool IsLaser,bool IsMC=false);
                static int GetSign(double zenith,double azimuth,double planephi,double nz);
                double GetApar(double CC,double phi,double zenith,double azimuth);
                static double GetApar(double &Anz,double zenith,double azimuth,double planephi,double nz);
                bool CalPlane(double CC,double phi,double zenith,double azimuth,double &planephi,double &nz,int &signnz);
                bool CalPlane(double CC,double phi,double zenith,double azimuth,double xyzdir[3][3]);
                bool GetPlane(double &planephi,double &eplanephi,double &nz,double &enz,int itel=0,int type=3);
                bool GetPlane(double xyzdir[3][3],double exyzdir[3][3],int itel=0,int type=3);
                static void CooCorr(double incoo[3],double indir[2],double outcoo[3]);
                //static void CooCorr(double incoo[3],double indir[3],double outcoo[3]);
                static int CalTelDir(double CC,double phi,double planephi,double nz,int &nsol,double ele_out[4],double azi_out[4],int signA[4],double ele_in=PI/4,double azi_in=0);
                static int CalTelDir(double CC,double phi,double* lasercoo,double* laserdir,int &nsol,double ele_out[4],double azi_out[4],int signA[4],double ele_in=PI/4,double azi_in=0);
                int CalTelDir(double CC,double phi,double planephi,double nz,double &elevation,double &azimuth,double ele_in=PI/4,double azi_in=0);
                int CalTelDir(double CC,double phi,double* lasercoo,double* laserdir,double &elevation,double &azimuth,double ele_in=PI/4,double azi_in=0);
                int GetTelDir(double &elevation,double &errel,double &azimuth,double &erraz,double ele_in=PI/4,double azi_in=0);

                static void Getnz(double nz,double &nz1,double &nz2);
                static void Getnz(double incoo[3],double indir[2],double &planephi,double &nz);
                static void GetCCphi(double zenith,double azimuth,double planephi,double nz,double &CC,double &phi);
                static void GetCCphi(double zenith,double azimuth,double incoo[3],double indir[2],double &CC,double &phi);
                static bool CalPHIRange0(double zenith,double azimuth,double planephi,double nz,double PHI_in,double* PHIRange);
                static bool GetImageXYCoo(double zenith,double azimuth,double planephi,double nz,double PHI_in,double &xx,double &yy);
                static bool GetImageXYCoo(double zenith,double azimuth,double* incoo,double* indir,double PHI_in,double &xx,double &yy);
                static void GetPHI(double zenith,double azimuth,double planephi,double nz,double* ImageCoo,double &PHI_in);
                static void GetPHI(double zenith,double azimuth,double incoo[3],double indir[2],double* ImageCoo,double &PHI_in);
                void GetPHI2(double zenith,double azimuth,double CC,double phi,double* ImageCoo,double &PHI_in);
                static int GetRange(double zenith,double azimuth,double planephi,double dirin[3],double phiin,double CCin,double* PHI,double* XX,double* YY);
                static bool CalRotateZeroPos(double ele_rotate,double azi_rotate,double ele_ref,double azi_ref,double &ele0,double &azi0);
                static bool CalDir_out(double ele_rotate,double azi_rotate,double ele_ref,double azi_ref,double ele_in,double azi_in,double &ele_out,double &azi_out);
                static bool Calnz_phiL(double planephi,double indir[2],double &nz,double &PHIL);
                static double CalTime(double ele_tel,double azi_tel,double incoo[3],double indir[2],double ImageCooXY[2],double time0=0);
                static double CalTime(double CC,double phi,double incoo[3],double nz,double PHIL,double ImageCooXY[2],double time0=0,double ele_in=PI/4,double azi_in=0);
                static TGraph* DrawTimeLine(double CC,double phi,double incoo[3],double nz,double PHIL,double time0=0,double ele_in=PI/4,double azi_in=0);
                TGraph* DrawTimeLine(double ele_rotate,double azi_rotate,double ele_ref,double azi_ref,double ele_in,double azi_in,double ele_tel,double azi_tel,double incoo[3],double time0=0);

	        TH2Poly* Draw(int type=0,const char* opt="scat colz",bool DoClean=true,double threshold=500.);
                TGraph* DrawPulse(int iSiPM,const char* opt="apl",bool IsHigh=true,bool DoClean=true);
                void DrawFit();
                static TGraph* DrawImageLine(double zenith,double azimuth,double incoo[3],double indir[2]);
                TGraph* DrawImageLine(int itel=0);
                TGraph* DrawCorePos(double* corepos,int itel=0,int type=3);
                TGraph* DrawCoreReg(double* corepos,int itel=0,int type=3);
                //static void slaDtp2s ( double xi, double eta, double raz, double decz, double &ra, double &dec );
	        TH2Poly* DrawGlobal(int type=0,const char* opt="colz",bool DoClean=true,double threshold=500.);
                TH2Poly* DrawCloudFormat(int type=0,const char* opt="colz",bool DoClean=true,double threshold=500.);
                TGraph2D* Draw3D(int type,const char* opt,double threshold,int ViewOpt=0);

                bool IsLed(int nfire_threshold=1000);
                bool IsLaser();
                bool IsNoise(double* pars);
                bool IsCR(double* pars);
                void CalInfo(double result[100]);

		void rabbittime2lt();
		void InitImage();
		void SetImage();
		void AdcToPe(float *deltag_20, float *correct_PreTemp, int isledevent);
		//void AdcToPe();
		void PrelimImageClean(double cut);
		void GetNeighborPixs();
		int CalcHillas();
		void GetCleanImage();


		ClassDef(WFCTAEvent,3);
};

//#if !defined(CINT)
ClassImp(WFCTAEvent);
//#endif

#endif // WFCTAEVENT_H
