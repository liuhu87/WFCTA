#include <stdlib.h>
#include <iostream>
#include "../include/WFCTAEvent.h"

using namespace std;

ClassImp(WFCTAEvent);

WFCTAEvent::WFCTAEvent()
{
}

WFCTAEvent::~WFCTAEvent()
{
}

void WFCTAEvent::EventInitial()
{
	iSiPM.clear();
	gain_marker.clear();
	peak.clear();
	mypeak.clear();
	Single_Threshold.clear();
	Record_Threshold.clear();
	Over_Single_Marker.clear();
	Over_Record_Marker.clear();
	ADC_Cut.clear();
	ImageAdcHigh.clear();
	ImageAdcLow.clear();
	ImageBaseHigh.clear();
	ImageBaseLow.clear();
	myImageAdcHigh.clear();
	myImageAdcLow.clear();
	myImageBaseHigh.clear();
	myImageBaseLow.clear();
}

void WFCTAEvent::SetEvent(long event)
{
	iEvent = event;
}

void WFCTAEvent::SetrabbitTime(long rabTime)
{
        rabbitTime = rabTime;
}
void WFCTAEvent::Setrabbittime(double rabtime)
{
        rabbittime = rabtime;
}

void WFCTAEvent::SetSiPM(short SiPM)
{
        iSiPM.push_back(SiPM);
}

void WFCTAEvent::Setgain_marker(bool w_gain_marker)
{
        gain_marker.push_back(w_gain_marker); 
}

void WFCTAEvent::Setpeak(char w_peak)
{
        peak.push_back(w_peak); 
}

void WFCTAEvent::Setmypeak(char w_mypeak)
{
        mypeak.push_back(w_mypeak); 
}

void WFCTAEvent::SetSingle_Threshold(short w_Single_Threshold)
{
        Single_Threshold.push_back(w_Single_Threshold); 
}

void WFCTAEvent::SetRecord_Threshold(short w_Record_Threshold)
{
        Record_Threshold.push_back(w_Record_Threshold); 
}

void WFCTAEvent::SetOver_Single_Marker(bool w_Over_Single_Marker)
{
        Over_Single_Marker.push_back(w_Over_Single_Marker); 
}

void WFCTAEvent::SetOver_Record_Marker(bool w_Over_Record_Marker)
{
        Over_Record_Marker.push_back(w_Over_Record_Marker); 
}

void WFCTAEvent::SetADC_Cut(float w_ADC_Cut)
{
        ADC_Cut.push_back(w_ADC_Cut); 
}

void WFCTAEvent::SetImageBaseHigh(float w_ImageBaseHigh)
{
        ImageBaseHigh.push_back(w_ImageBaseHigh); 
}

void WFCTAEvent::SetImageBaseLow(float w_ImageBaseLow)
{
        ImageBaseLow.push_back(w_ImageBaseLow); 
}

void WFCTAEvent::SetImageAdcHigh(float w_ImageAdcHigh)
{
        ImageAdcHigh.push_back(w_ImageAdcHigh); 
}

void WFCTAEvent::SetImageAdcLow(float w_ImageAdcLow)
{
        ImageAdcLow.push_back(w_ImageAdcLow); 
}

void WFCTAEvent::SetmyImageBaseHigh(float w_myImageBaseHigh)
{
        myImageBaseHigh.push_back(w_myImageBaseHigh); 
}

void WFCTAEvent::SetmyImageBaseLow(float w_myImageBaseLow)
{
        myImageBaseLow.push_back(w_myImageBaseLow);
}

void WFCTAEvent::SetmyImageAdcHigh(float w_myImageAdcHigh)
{
        myImageAdcHigh.push_back(w_myImageAdcHigh);
}

void WFCTAEvent::SetmyImageAdcLow(float w_myImageAdcLow)
{
        myImageAdcLow.push_back(w_myImageAdcLow);
}


void WFCTAEvent::StatusInitial()
{
//	DbTempCount = 0;
//	ClbTimeCount = 0;
//	ClbTempCount = 0;
}

