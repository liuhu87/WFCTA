#ifndef WFCTAEVENT_H
#define WFCTAEVENT_H

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include "TSelector.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TTree.h"
#include "TBranch.h"

//#include "WFCTAMCEvent.h"
#include "WFCTALedEvent.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

using namespace std;
const int MAXPMT=1024;
class WFCTAEvent : public TSelector
{
	protected:
		static const int cleanPix;
		static WFCTAEvent* _Head;
		static TTree* _Tree;
		static const char * _Name;
		static TBranch* bAll;
		static TBranch* bledevent;
		static TBranch* blaserevent;

	public:
		static int jdebug;
		static double adccuttrigger;
		static double npetrigger;
		static int nfiretrigger;

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
		//short packCheck;
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

		WFCTALedEvent ledevent;  //!

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
		int GetMaxADCBin(int itel=0);

		static void slaDtp2s ( double xi, double eta, double raz, double decz, double &ra, double &dec );

		bool IsNoise(int p0=5,double p1=150,int p2=4,double p3=3.5);

		void rabbittime2lt();
		void InitImage();
		void SetImage();
		void AdcToPe(float *deltag_20, float *correct_PreTemp, int isledevent);
		//void AdcToPe();
		void PrelimImageClean(double cut, int laserEvent, int ledEvent);
		void GetNeighborPixs();
		int CalcHillas(int ledEvent);
		void GetCleanImage(int LedEvent);


		ClassDef(WFCTAEvent,3);
};

//#if !defined(CINT)
ClassImp(WFCTAEvent);
//#endif

#endif // WFCTAEVENT_H
