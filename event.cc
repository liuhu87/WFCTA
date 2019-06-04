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
  eventShow -> Branch("Event","WFCTAEvent",wfctaEvent,32000,1);
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

      if(wfctaDecode->bigPackCheck(buf,BUF_LEN))
      {
	  //get info eventID and rabbit_time//
          packSize = wfctaDecode->PackSize();

          wfctaEvent->SetEvent(wfctaDecode->eventId(buf));
          wfctaEvent->SetrabbitTime(wfctaDecode->RabbitTime(buf));
          wfctaEvent->Setrabbittime(wfctaDecode->Rabbittime(buf));

	  //find sipms and their position in this pack//
	  wfctaDecode->Find_SiPMs(buf,packSize);
	  sipm_position = &(wfctaDecode->GetSiPM_Position());

	  //get info of each sipm: q, base, peakposition...//
	  for(sipm_position_iter=sipm_position->begin(); sipm_position_iter!=sipm_position->end(); sipm_position_iter++){
	    wfctaEvent->SetSiPM( sipm_position_iter->first );
	    wfctaEvent->Setpeak( wfctaDecode->GetPeak(buf,sipm_position_iter->first) );
	    wfctaEvent->SetSingle_Threshold( wfctaDecode->GetSingle_Thresh(buf,sipm_position_iter->first) );
	    wfctaEvent->SetRecord_Threshold( wfctaDecode->GetRecord_Thresh(buf,sipm_position_iter->first) );
	    wfctaEvent->SetOver_Single_Marker( wfctaDecode->GetOver_Single_Mark(buf,sipm_position_iter->first) );
	    wfctaEvent->SetOver_Record_Marker( wfctaDecode->GetOver_Record_Mark(buf,sipm_position_iter->first) );
	    wfctaEvent->SetImageBaseHigh( wfctaDecode->BaseHigh(buf,sipm_position_iter->first) );
            wfctaEvent->SetImageBaseLow( wfctaDecode->BaseLow(buf,sipm_position_iter->first) );
            wfctaEvent->SetImageAdcHigh( wfctaDecode->AdcHigh(buf,sipm_position_iter->first) );
            wfctaEvent->SetImageAdcLow( wfctaDecode->AdcLow(buf,sipm_position_iter->first) );
	    //wfctaDecode->GetWaveForm(buf,sipm_position_iter->first,(int *)pulsehigh, (int *)pulselow);
	    wfctaEvent->Setgain_marker( wfctaDecode->Getgain_marker(buf,sipm_position_iter->first) );
	    wfctaEvent->Setmypeak( wfctaDecode->Getmypeak(buf,sipm_position_iter->first) );
            wfctaEvent->SetADC_Cut( wfctaDecode->GetADC_Cut(buf,sipm_position_iter->first) );
	    wfctaEvent->SetmyImageBaseHigh( wfctaDecode->GetmyImageBaseHigh(buf,sipm_position_iter->first) );
            wfctaEvent->SetmyImageBaseLow( wfctaDecode->GetmyImageBaseLow(buf,sipm_position_iter->first) );
            wfctaEvent->SetmyImageAdcHigh( wfctaDecode->GetmyImageAdcHigh(buf,sipm_position_iter->first) );
            wfctaEvent->SetmyImageAdcLow( wfctaDecode->GetmyImageAdcLow(buf,sipm_position_iter->first) );
	  }

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
