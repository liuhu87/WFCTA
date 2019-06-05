#include <stdlib.h>
#include <iostream>
#include "WFCTAEvent.h"

using namespace std;

//ClassImp(WFCTAEvent);

WFCTAEvent::WFCTAEvent()
{
   Init();
}

WFCTAEvent::~WFCTAEvent()
{
   Reset();
}

void WFCTAEvent::Init()
{
   iEvent=-1;
   rabbitTime=0;
   rabbittime=0;
   iSiPM.clear();
   gain_marker.clear();
   peak.clear();
   mypeak.clear();
   Single_Threshold.clear();
   Record_Threshold.clear();
   Over_Single_Marker.clear();
   Over_Record_Marker.clear();
   ADC_Cut.clear();
   ImageBaseHigh.clear();
   ImageBaseLow.clear();
   ImageAdcHigh.clear();
   ImageAdcLow.clear();
   myImageBaseHigh.clear();
   myImageBaseLow.clear();
   myImageAdcHigh.clear();
   myImageAdcLow.clear();
   
   mcevent.Init();
   ledevent.Init();
   laserevent.Init();
}
void WFCTAEvent::Reset()
{
   iEvent=-1;
   rabbitTime=0;
   rabbittime=0;
   iSiPM.clear();
   gain_marker.clear();
   peak.clear();
   mypeak.clear();
   Single_Threshold.clear();
   Record_Threshold.clear();
   Over_Single_Marker.clear();
   Over_Record_Marker.clear();
   ADC_Cut.clear();
   ImageBaseHigh.clear();
   ImageBaseLow.clear();
   ImageAdcHigh.clear();
   ImageAdcLow.clear();
   myImageBaseHigh.clear();
   myImageBaseLow.clear();
   myImageAdcHigh.clear();
   myImageAdcLow.clear();
   
   mcevent.Reset();
   ledevent.Reset();
   laserevent.Reset();
}


