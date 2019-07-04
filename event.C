#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include "include/dumpPack.h"
#include "include/WFCTAEvent.h"
#include "include/WFCTADecode.h"
#include <TFile.h>
#include <TTree.h>

#define BUF_LEN 200000
#define STATUS_BUF_LEN 200000

using namespace std;

int main(int argc, char**argv)
{
  FILE *fp;
  uint8_t *buf = new uint8_t[BUF_LEN];
  size_t size_of_read;
  fp = fopen(argv[1],"rb");

  int64_t packSize = 0;
  map<short, int>* sipm_position;
  map<short, int>::iterator sipm_position_iter;
  int pulsehigh[32];
  int pulselow[32];

  WFCTAEvent *wfctaEvent = new WFCTAEvent();
  TFile *rootfile = new TFile(argv[2],"recreate");
  /*********************************************************************/
  TTree *eventShow = new TTree("eventShow","info of evnets");
  wfctaEvent -> CreateBranch(eventShow,1);
  /*********************************************************************/

  WFCTADecode *wfctaDecode = new WFCTADecode();

  //Events Initial//
  wfctaEvent->EventInitial();
  for(int j=0;j<32;j++){
    pulsehigh[j] = 0;
    pulselow[j] = 0;
  }

  fp = fopen(argv[1],"rb");
  while(true)
  {
      size_of_read = fread((uint8_t *)buf,1,BUF_LEN,fp);
      fseek(fp,-size_of_read,1);
      if(size_of_read==0){break;}
      //dumpPacket(buf,size_of_read,16);

      if(wfctaDecode->bigPackCheck(buf,int(size_of_read)))
      {
	  //get info eventID and rabbit_time//
          packSize = wfctaDecode->PackSize();

          wfctaEvent->iEvent=wfctaDecode->eventId(buf);
          wfctaEvent->rabbitTime=wfctaDecode->RabbitTime(buf);
          wfctaEvent->rabbittime=wfctaDecode->Rabbittime(buf);

	  //find sipms and their position in this pack//
	  wfctaDecode->Find_SiPMs(buf,packSize);
	  sipm_position = &(wfctaDecode->GetSiPM_Position());

	  //get info of each sipm: q, base, peakposition...//
	  for(sipm_position_iter=sipm_position->begin(); sipm_position_iter!=sipm_position->end(); sipm_position_iter++){
	    wfctaEvent->iSiPM.push_back( sipm_position_iter->first );
	    wfctaEvent->peak.push_back( wfctaDecode->GetPeak(buf,sipm_position_iter->first) );
	    wfctaEvent->Single_Threshold.push_back( wfctaDecode->GetSingle_Thresh(buf,sipm_position_iter->first) );
	    wfctaEvent->Record_Threshold.push_back( wfctaDecode->GetRecord_Thresh(buf,sipm_position_iter->first) );
	    wfctaEvent->Over_Single_Marker.push_back( wfctaDecode->GetOver_Single_Mark(buf,sipm_position_iter->first) );
	    wfctaEvent->Over_Record_Marker.push_back( wfctaDecode->GetOver_Record_Mark(buf,sipm_position_iter->first) );
	    wfctaEvent->ImageBaseHigh.push_back( wfctaDecode->BaseHigh(buf,sipm_position_iter->first) );
            wfctaEvent->ImageBaseLow.push_back( wfctaDecode->BaseLow(buf,sipm_position_iter->first) );
            wfctaEvent->ImageAdcHigh.push_back( wfctaDecode->AdcHigh(buf,sipm_position_iter->first) );
            wfctaEvent->ImageAdcLow.push_back( wfctaDecode->AdcLow(buf,sipm_position_iter->first) );
	    //wfctaDecode->GetWaveForm(buf,sipm_position_iter->first,(int *)pulsehigh, (int *)pulselow);
	    wfctaEvent->gain_marker.push_back( wfctaDecode->Getgain_marker(buf,sipm_position_iter->first) );
	    wfctaEvent->mypeak.push_back( wfctaDecode->Getmypeak(buf,sipm_position_iter->first) );
            wfctaEvent->ADC_Cut.push_back( wfctaDecode->GetADC_Cut(buf,sipm_position_iter->first) );
	    wfctaEvent->myImageBaseHigh.push_back( wfctaDecode->GetmyImageBaseHigh(buf,sipm_position_iter->first) );
            wfctaEvent->myImageBaseLow.push_back( wfctaDecode->GetmyImageBaseLow(buf,sipm_position_iter->first) );
            wfctaEvent->myImageAdcHigh.push_back( wfctaDecode->GetmyImageAdcHigh(buf,sipm_position_iter->first) );
            wfctaEvent->myImageAdcLow.push_back( wfctaDecode->GetmyImageAdcLow(buf,sipm_position_iter->first) );
	  }

          ///just do some test to exam the reading program
          //(wfctaEvent->mcevent).iuse=wfctaEvent->iEvent;
          //(wfctaEvent->mcevent).RayTrace.push_back(wfctaEvent->iEvent*2);
          //(wfctaEvent->mcevent).TelTrigger[0]=wfctaEvent->iEvent*3;
          //(wfctaEvent->ledevent).Time=wfctaEvent->rabbitTime;
          //(wfctaEvent->ledevent).Frequency=wfctaEvent->iEvent;
          //(wfctaEvent->ledevent).DoorOpen=(wfctaEvent->iEvent%2);
          //(wfctaEvent->laserevent).Time=wfctaEvent->rabbitTime;
          //(wfctaEvent->laserevent).Frequency=wfctaEvent->iEvent;
          //(wfctaEvent->laserevent).flux=wfctaEvent->iEvent*1.0;

          eventShow->Fill();
	  wfctaEvent->EventInitial();
          for(int j=0;j<32;j++){
            pulsehigh[j] = 0;
            pulselow[j] = 0;
          }
          fseek(fp,packSize,1);
      }
      else
      {
          fseek(fp,1,1);
      }
  }
  fclose(fp);

/******************************************************************************/
  //eventShow->Write();
  rootfile->Write();
  rootfile->Close();

}
