#ifndef WFCTAEVENT_H
#define WFCTAEVENT_H

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include "TObject.h"

#include "WFCTAMCEvent.h"
#include "WFCTALedEvent.h"
#include "WFCTALaserEvent.h"

using namespace std;
class WFCTAEvent : public TSelector
{
public:
        long iEvent;
        long rabbitTime;
        double rabbittime;
        vector<short> iSiPM;
	vector<bool> gain_marker;
	vector<char> peak;
	vector<char> mypeak;
	vector<short> Single_Threshold;
	vector<short> Record_Threshold;
	vector<bool> Over_Single_Marker;
	vector<bool> Over_Record_Marker;
	vector<float> ADC_Cut;
	vector<float> ImageBaseHigh;
	vector<float> ImageBaseLow;
	vector<float> ImageAdcHigh;
	vector<float> ImageAdcLow;
	vector<float> myImageBaseHigh;
	vector<float> myImageBaseLow;
	vector<float> myImageAdcHigh;
	vector<float> myImageAdcLow;

        WFCTAMCEvent mcevent;
        WFCTALedEvent ledevent;
        WFCTALaserEvent laserevent;

public:
        WFCTAEvent();
        ~WFCTAEvent();
        void Init();
        void EventInitial();

   ClassDef(WFCTAEvent,2);
};

#endif // WFCTAEVENT_H
