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
	vector<float> eBaseH;  
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
	vector<char> peak;  //!
	vector<char> PeakPosH;
        vector<char> PeakPosL;
	vector<int> PeakAmH;
        vector<int> PeakAmL;
	vector<bool> gain_marker;  //!
	vector<bool> Over_Single_Marker;
	vector<bool> Over_Record_Marker;

	int Npoint[28]; //! 
	int pulsehigh[1024][28]; //!
	int pulselow[1024][28];  //!
        double ImageX[1024];  //!
        double ImageY[1024];  //!

        WFCTAMCEvent mcevent;  //!
        WFCTALedEvent ledevent;  //!
        WFCTALaserEvent laserevent;  //!

public:
        int year;  //!
        int month;  //!
        int day;  //!
        int hour;  //!
        int minite;  //!
        int second;  //!
        vector<double> FullImagePe;  //!
        vector<double> FullImageX;  //!
	vector<double> FullImageY;  //!
	vector<double> CleanImagePe;  //!
	vector<double> CleanImageX;  //!
	vector<double> CleanImageY;  //!

	vector<double>::iterator it_pe;  //!
	vector<double>::iterator it_x;  //!
	vector<double>::iterator it_y;  //!
	vector<double>::iterator it_pe0;  //!
	vector<double>::iterator it_x0;  //!
	vector<double>::iterator it_y0;  //!

        vector<int> fNpixfriends;  //!
	vector<int>::iterator it_npix;  //!
	int Npix;  //!

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
        void CalculateADC(int itel=0);
        int GetMaxADCBin();
        int GetMaxTimeBin();
        int GetMinTimeBin();
	TH2Poly* Draw(int type=0,const char* opt="scat colz",double threshold=500.);
        TObjArray* Draw3D(int type,const char* opt,double threshold,int ViewOpt=0);

	void SetImage();
	void rabbittime2lt();
	void InitImage();
	void AdcToPe();
	void ImageClean(double cut);

   ClassDef(WFCTAEvent,3);
};

//#if !defined(CINT)
ClassImp(WFCTAEvent);
//#endif

#endif // WFCTAEVENT_H
