#ifndef WFCTAEVENT_H
#define WFCTAEVENT_H

#include <vector>
#include <string>
#include <map>
#include "TObject.h"

using namespace std;
class WFCTAEvent //: public TObject 
{
public:
        WFCTAEvent();
        ~WFCTAEvent();
        void EventInitial();
	void StatusInitial();
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
	void SetEvent(long event);
	void SetrabbitTime(long rabTime);
	void Setrabbittime(double rabtime);
        void SetSiPM(short SiPM);
        void Setgain_marker(bool w_gain_marker);
        void Setpeak(char w_peak);
        void Setmypeak(char w_mypeak);
        void SetSingle_Threshold(short w_Single_Threshold);
        void SetRecord_Threshold(short w_Record_Threshold);
        void SetOver_Single_Marker(bool w_Over_Single_Marker);
        void SetOver_Record_Marker(bool w_Over_Record_Marker);
        void SetADC_Cut(float w_ADC_Cut);
        void SetImageBaseHigh(float w_ImageBaseHigh);
        void SetImageBaseLow(float w_ImageBaseLow);
        void SetImageAdcHigh(float w_ImageAdcHigh);
        void SetImageAdcLow(float w_ImageAdcLow);
        void SetmyImageBaseHigh(float w_myImageBaseHigh);
        void SetmyImageBaseLow(float w_myImageBaseLow);
        void SetmyImageAdcHigh(float w_myImageAdcHigh);
        void SetmyImageAdcLow(float w_myImageAdcLow);

   ClassDef(WFCTAEvent,2);
};

#endif // WFCTAEVENT_H
