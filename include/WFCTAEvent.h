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
   static TBranch* bEvent;
   static TBranch* bTime;
   static TBranch* btime;
   static TBranch* bADC_Cut;
   static TBranch* bBaseHigh;
   static TBranch* bBaseLow;
   static TBranch* bAdcHigh;
   static TBranch* bAdcLow;
   static TBranch* bmyBaseHigh;
   static TBranch* bmyBaseLow;
   static TBranch* bmyAdcHigh;
   static TBranch* bmyAdcLow;
   static TBranch* bSiPM;
   static TBranch* bSThr;
   static TBranch* bRThr;
   static TBranch* bpeak;
   static TBranch* bmypeak;
   static TBranch* bgain_marker;
   static TBranch* bOSMarker;
   static TBranch* bORMarker;
   static TBranch* bmcevent;
   static TBranch* bledevent;
   static TBranch* blaserevent;

public:
	long iEvent;
	long rabbitTime;
	double rabbittime;
	vector<float> ADC_Cut;
	vector<float> ImageBaseHigh;
	vector<float> ImageBaseLow;
	vector<float> ImageAdcHigh;
	vector<float> ImageAdcLow;
	vector<float> myImageBaseHigh;
	vector<float> myImageBaseLow;
	vector<float> myImageAdcHigh;
	vector<float> myImageAdcLow;
        //vector<int> iSiPM;
	//vector<int> Single_Threshold;
	//vector<int> Record_Threshold;
	//vector<int> peak;
	//vector<int> mypeak;
	//vector<int> gain_marker;
	//vector<int> Over_Single_Marker;
	//vector<int> Over_Record_Marker;
        
        vector<short> iSiPM;
	vector<short> Single_Threshold;
	vector<short> Record_Threshold;
	vector<char> peak;
	vector<char> mypeak;
	vector<bool> gain_marker;
	vector<bool> Over_Single_Marker;
	vector<bool> Over_Record_Marker;

        WFCTAMCEvent mcevent;
        WFCTALedEvent ledevent;
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
        void CalculateADC(int itel=0);
	TH2Poly* Draw(int type=0,const char* opt="scat colz",double threshold=500.);
        TObjArray* Draw3D(int type,const char* opt,double threshold,int ViewOpt=0);

   ClassDef(WFCTAEvent,3);
};

//#if !defined(CINT)
ClassImp(WFCTAEvent);
//#endif

#endif // WFCTAEVENT_H
