#ifndef WFCTAEVENT_H
#define WFCTAEVENT_H

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include "TObject.h"
#include "TObjArray.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH2Poly.h"

#include "WFCTAMCEvent.h"
#include "WFCTALedEvent.h"
#include "WFCTALaserEvent.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

using namespace std;
const int MAXPMT=1024;
class WFCTAEvent : public TSelector
{
protected:
   static WFCTAEvent* _Head;
   static TTree* _Tree;
   static const char * _Name;
   static TBranch* bAll;
   static TBranch* bmcevent;
   static TBranch* bledevent;
   static TBranch* blaserevent;

public:
        short iTel;
	long iEvent;
	long rabbitTime;
	double rabbittime;
	int big_pack_lenth;
	short n_fired;
	short n_Channel;

	vector<long> ievent;
	vector<float> ADC_Cut;  //!
	vector<float> ImageBaseHigh;  
	vector<float> ImageBaseLow;  //!
	vector<float> ImageAdcHigh;  //!
	vector<float> ImageAdcLow;  //!
	vector<float> myImageBaseHigh;
	vector<float> myImageBaseLow;
	vector<float> myImageAdcHigh;
	vector<float> myImageAdcLow;
        vector<short> iSiPM;
	vector<short> Single_Threshold;  //!
	vector<short> Record_Threshold;  //!
	vector<int> peak;  //!
	vector<int> mypeak;
	vector<int> peakamp;
	vector<bool> gain_marker;  //!
	vector<bool> Over_Single_Marker;
	vector<bool> Over_Record_Marker;

	int Npoint[28];  //!
	int pulsehigh[1024][28];  //!
	int pulselow[1024][28];  //!

        /// for the fitting of the image
        TGraph* gDraw; //!
        TGraph* gDrawErr; //!
        ROOT::Math::Minimizer* minimizer; //!

        WFCTAMCEvent mcevent;  
        WFCTALedEvent ledevent;  //!
        WFCTALaserEvent laserevent;

public:
        WFCTAEvent();
        ~WFCTAEvent();
        void Init();
        void EventInitial();
        void InitTree(TTree* tree);
        static const char *  BranchName() {return _Name;}
        static WFCTAEvent* & Head()  {return _Head;}
        static TTree* & Tree()  {return _Tree;}
        void CreateBranch(TTree *tree, int branchSplit);
        void GetBranch(TTree *fChain);
        bool GetAllContents(int _Entry);
        bool CheckLaser();
        bool CheckMC();
        void CalculateDataVar(int itel=0);
        int GetMaxADCBin(int itel=0);
        int GetPeakADCBin(int isipm,int itel=0);
        int GetMaxTimeBin(int itel=0);
        int GetMinTimeBin(int itel=0);
        double GetContent(int isipm,int itel=0,int type=3,bool IsIndex=false);
        bool CleanImage(int isipm,int itel=0,int type=3,bool IsIndex=false);
        double Interface(const double* par);
        bool DoFit(int itel=0,int type=3,bool force=false);
        bool GetCrossCoor(double x, double y, double &x0, double &y0);
        TH1F* GetLongDistribution(int itel=0,int type=3);
        TH1F* GetShortDistribution(int itel=0,int type=3);
        int GetSign(bool IsLaser,bool IsMC);
        double GetApar(double CC,double phi,double zenith,double azimuth);
        bool CalPlane(double CC,double phi,double zenith,double azimuth,double &planephi,double &nz);
        bool CalPlane(double CC,double phi,double zenith,double azimuth,double xyzdir[3][3]);
        bool GetPlane(double &planephi,double &eplanephi,double &nz,double &enz,int itel=0,int type=3);
        bool GetPlane(double xyzdir[3][3],double exyzdir[3][3],int itel=0,int type=3);
        int CalTelDir(double CC,double phi,double* lasercoo,double* laserdir,double &elevation,double &azimuth);
        int GetTelDir(double &elevation,double &errel,double &azimuth,double &erraz);
        static void GetImageXYCoo(double zenith,double azimuth,double* lasercoo,double* laserdir,double PHI_in,double &xx,double &yy);
        void GetPHI(double zenith,double azimuth,double CC,double phi,double* ImageCoo,double &PHI_in);
	TH2Poly* Draw(int type=0,const char* opt="scat colz",double threshold=500.);
        void DrawFit();
        TGraph* DrawImageLine(int itel=0);
        TGraph* DrawCorePos(double* corepos,int itel=0,int type=3);
        TGraph* DrawCoreReg(double* corepos,int itel=0,int type=3);
        static void slaDtp2s ( double xi, double eta, double raz, double decz, double &ra, double &dec );
	TH2Poly* DrawGlobal(int type=0,const char* opt="scat colz",double threshold=500.);
        TGraph2D* Draw3D(int type,const char* opt,double threshold,int ViewOpt=0);

   ClassDef(WFCTAEvent,3);
};

//#if !defined(CINT)
ClassImp(WFCTAEvent);
//#endif

#endif // WFCTAEVENT_H
