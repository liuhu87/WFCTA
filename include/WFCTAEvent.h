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
   static int jdebug;
   static bool DoDraw;

        short iTel;
	long iEvent;
        long eEvent;
	long rabbitTime;
	double rabbittime;
	int big_pack_lenth;
	short n_fired;
	short n_Channel;

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
	vector<float> AdcH;
	vector<float> AdcL;
        vector<bool> SatH;
        vector<bool> SatL;
	vector<short> Single_Threshold;  //!
	vector<short> Record_Threshold;  //!

	vector<int> peak;  //!
	vector<int> mypeak;
	vector<int> peakamp;

	//vector<char> peak;  //!
	vector<char> PeakPosH;
        vector<char> PeakPosL;
	vector<int> PeakAmH;
        vector<int> PeakAmL;

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

        double ImageX[1024];  //!
        double ImageY[1024];  //!

        WFCTAMCEvent mcevent;  
        WFCTALedEvent ledevent;  //!
        WFCTALaserEvent laserevent;

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
        static double GetApar(double &Anz,double zenith,double azimuth,double planephi,double nz);
        bool CalPlane(double CC,double phi,double zenith,double azimuth,double &planephi,double &nz,int &signnz);
        bool CalPlane(double CC,double phi,double zenith,double azimuth,double xyzdir[3][3]);
        bool GetPlane(double &planephi,double &eplanephi,double &nz,double &enz,int itel=0,int type=3);
        bool GetPlane(double xyzdir[3][3],double exyzdir[3][3],int itel=0,int type=3);
        int CalTelDir(double CC,double phi,double* lasercoo,double* laserdir,double &elevation,double &azimuth);
        int GetTelDir(double &elevation,double &errel,double &azimuth,double &erraz);

        static void Getnz(double incoo[3],double indir[2],double &planephi,double &nz,double &nz1,double &nz2);
        static void GetCCphi(double zenith,double azimuth,double planephi,double nz,double &CC,double &phi);
        static void GetCCphi(double zenith,double azimuth,double incoo[3],double indir[2],double &CC,double &phi);
        static bool GetImageXYCoo(double zenith,double azimuth,double planephi,double nz1,double nz2,double PHI_in,double &xx,double &yy);
        static bool GetImageXYCoo(double zenith,double azimuth,double* incoo,double* indir,double PHI_in,double &xx,double &yy);
        static void GetPHI(double zenith,double azimuth,double planephi,double nz1,double nz2,double* ImageCoo,double &PHI_in);
        static void GetPHI(double zenith,double azimuth,double incoo[3],double indir[2],double* ImageCoo,double &PHI_in);
        void GetPHI(double zenith,double azimuth,double CC,double phi,double* ImageCoo,double &PHI_in);
        static int GetRange(double zenith,double azimuth,double planephi,double dirin[3],double phiin,double CCin,double* PHI,double* XX,double* YY);

	TH2Poly* Draw(int type=0,const char* opt="scat colz",double threshold=500.);
        void DrawFit();
        static TGraph* DrawImageLine(double zenith,double azimuth,double incoo[3],double indir[2]);
        TGraph* DrawImageLine(int itel=0);
        TGraph* DrawCorePos(double* corepos,int itel=0,int type=3);
        TGraph* DrawCoreReg(double* corepos,int itel=0,int type=3);
        static void slaDtp2s ( double xi, double eta, double raz, double decz, double &ra, double &dec );
	TH2Poly* DrawGlobal(int type=0,const char* opt="scat colz",double threshold=500.);
        TGraph2D* Draw3D(int type,const char* opt,double threshold,int ViewOpt=0);

        void rabbittime2lt();
        void InitImage();
	void SetImage();
	void AdcToPe();
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
