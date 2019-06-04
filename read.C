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

  uint8_t status_pack_marker;

  int fpgaVersion[10];
  long clb_initial_Time;
  double clb_initial_time;
  int fired_tube;
  long status_readback_Time;
  double status_readback_time;
  //int sipm[1024];  for(int i=0;i<1024;i++) {sipm[i]=i;}
  short single_thresh[1024];
  short record_thresh[1024];
  long single_count[1024];
  float DbTemp[1024];
  long single_time[1024];
  float HV[1024];
  float PreTemp[1024];
  float BigResistence[1024];
  float SmallResistence[1024];
  long ClbTime[1024];
  float ClbTemp[1024];

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
/*
  eventShow -> Branch("rabbitTime",&wfctaEvent->rabbitTime);
  eventShow -> Branch("rabbittime",&wfctaEvent->rabbittime);
  eventShow -> Branch("iSiPM",&wfctaEvent->iSiPM);
  eventShow -> Branch("gain_marker",&wfctaEvent->gain_marker);
  eventShow -> Branch("peak",&wfctaEvent->peak);
  eventShow -> Branch("mypeak",&wfctaEvent->mypeak);
  eventShow -> Branch("Single_Threshold",&wfctaEvent->Single_Threshold);
  eventShow -> Branch("Record_Threshold",&wfctaEvent->Record_Threshold);
  eventShow -> Branch("Over_Single_Marker",&wfctaEvent->Over_Single_Marker);
  eventShow -> Branch("Over_Record_Marker",&wfctaEvent->Over_Record_Marker);
  eventShow -> Branch("ADC_Cut",&wfctaEvent->ADC_Cut);
  eventShow -> Branch("ImageBaseHigh",&wfctaEvent->ImageBaseHigh);
  eventShow -> Branch("ImageBaseLow",&wfctaEvent->ImageBaseLow);
  eventShow -> Branch("ImageAdcHigh",&wfctaEvent->ImageAdcHigh);
  eventShow -> Branch("ImageAdcLow",&wfctaEvent->ImageAdcLow);
  eventShow -> Branch("myImageBaseHigh",&wfctaEvent->myImageBaseHigh);
  eventShow -> Branch("myImageBaseLow",&wfctaEvent->myImageBaseLow);
  eventShow -> Branch("myImageAdcHigh",&wfctaEvent->myImageAdcHigh);
  eventShow -> Branch("myImageAdcLow",&wfctaEvent->myImageAdcLow);
*/
  /*********************************************************************/

  WFCTADecode *wfctaDecode = new WFCTADecode();
  for(int i=0;i<1024;i++){
    single_thresh[i] = -1000;
    record_thresh[i] = -1000;
    single_count[i] = -1000;
    single_time[i] = -1000;
    DbTemp[i] = -1000; 
    HV[i] = -1000;
    PreTemp[i] = -1000;
    BigResistence[i] = -1000;
    SmallResistence[i] = -1000;
    ClbTime[i] = -1000;
    ClbTemp[i] = -1000;
  }
  for(int i=0;i<10;i++){
    fpgaVersion[i] = -1000;
  }
  while(true)
  {
      size_of_read = fread((uint8_t *)buf,1,STATUS_BUF_LEN,fp);
      fseek(fp,-size_of_read,1);
      if(size_of_read==0){break;}
      if(wfctaDecode->StatusPackCheck(buf,STATUS_BUF_LEN))
      {
	  status_pack_marker = wfctaDecode->StatusPackCheck(buf,STATUS_BUF_LEN);
          packSize = wfctaDecode->PackSize();
	  //printf("packSize:%lld | status_pack_marker:%x\n",packSize,status_pack_marker);
          //dumpPacket(buf,packSize,16);

            switch(status_pack_marker){
              case 0x21:
		wfctaDecode->Getthresh(buf,packSize,(short *)single_thresh, (short *)record_thresh);
                break;
              case 0x22:
		wfctaDecode->Deal22Pack(buf,packSize,(long *)single_count);
                break;
              case 0x23:
                wfctaDecode->Deal23Pack(buf,packSize,(long *)single_count,(long *)single_time);
                break;
              case 0x81:
		wfctaDecode->GetHV(buf,packSize,(float *)HV);
                break;
              case 0x82:
		wfctaDecode->GetPreTemp(buf,packSize,(float *)PreTemp);
                break;
              //case 0x83:
              //  wfctaDecode->Deal83Package((float *)BigResistence);
              //  break;
              //case 0x84:
              //  wfctaDecode->Deal84Package((float *)SmallResistence);
              //  break;
              case 0x85:
		wfctaDecode->GetClbTemp(buf,packSize,(float *)ClbTemp);
                break;
	      case 0x9:
		clb_initial_Time = wfctaDecode->GetclbInitialTime(buf,packSize);
		clb_initial_time = wfctaDecode->GetclbInitialtime(buf,packSize);
		fired_tube = wfctaDecode->GetFiredTube(buf,packSize);
		status_readback_Time = wfctaDecode->GetStatusReadbackTime(buf,packSize);
		status_readback_time = wfctaDecode->GetStatusReadbacktime(buf,packSize);

		printf("a status:\n");
		printf("%d %d\n",single_thresh[0],record_thresh[0]);
		printf("%ld %ld\n",single_count[0],single_time[0]);
		printf("%f %f %f\n",HV[0],PreTemp[0],ClbTemp[0]);
		printf("%ld %lf | %ld %lf\n\n",clb_initial_Time,clb_initial_time,status_readback_Time,status_readback_time);

		//Status->Fill();
		for(int i=0;i<1024;i++){
		  single_thresh[i] = -1000;
                  record_thresh[i] = -1000;
                  single_count[i] = -1000;
                  single_time[i] = -1000;
                  DbTemp[i] = -1000;
                  HV[i] = -1000;
                  PreTemp[i] = -1000;
                  BigResistence[i] = -1000;
                  SmallResistence[i] = -1000;
                  ClbTime[i] = -1000;
		  ClbTemp[i] = -1000;
                }
                for(int i=0;i<10;i++){
                  fpgaVersion[i] = -1000;
                }
		break;
	    }
            //if(status_pack_marker>0&&status_pack_marker<10){
            //  wfctaDecode->DealFPGAPackage((int *)fpgaVersion);
            //}

          fseek(fp,packSize,1);
      }
      else
      {
          fseek(fp,size_of_read,1);
      }
  }
  fclose(fp);


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
  eventShow->Write();
  //rootfile->Write();
  rootfile->Close();

}
