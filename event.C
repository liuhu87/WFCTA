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

#define BUF_LEN 2000000
#define STATUS_BUF_LEN 2000000

using namespace std;

int main(int argc, char**argv)
{
  if(argc!=3)
  {
      printf("Use %s inputfile outfile\n",argv[0]);
      return 0;
  }

  FILE *fp;
  uint8_t *buf = NULL;// = new uint8_t[BUF_LEN];
  int32_t slicelength;
  size_t size_of_read;
  short ITEL;
  int FEEDataHead;

  int64_t packStart = 0;
  map<short, int>* sipm_position;
  map<short, int>::iterator sipm_position_iter;

  TFile *rootfile = new TFile(argv[2],"recreate");
  WFCTAEvent *wfctaEvent = new WFCTAEvent();
  /*********************************************************************/
  TTree *eventShow = new TTree("eventShow","info of evnets");
  wfctaEvent -> CreateBranch(eventShow,1);
  /*********************************************************************/

  WFCTADecode *wfctaDecode = new WFCTADecode();

  //Events Initial//
  wfctaEvent->EventInitial();

  fp = fopen(argv[1],"rb");
  int nevent=0;
  while(true)
  {
      if(buf!=NULL)
        {delete buf;}
      buf = new uint8_t[40];
      size_of_read = fread((uint8_t *)buf,1,40,fp);
      if(size_of_read==0){break;}
      //dumpPacket(buf,24,16);
      if(wfctaDecode->FEEDataFragment(buf))
      {
	FEEDataHead = wfctaDecode->feeDataHead();
        slicelength = wfctaDecode->sliceLength(buf,FEEDataHead); 
	ITEL = wfctaDecode->Telid(buf,FEEDataHead);
        fseek(fp,-size_of_read+FEEDataHead,1);

	delete buf;
	buf = new uint8_t[slicelength];
	size_of_read = fread((uint8_t *)buf,1,slicelength,fp);
	//printf("slicelength:%lld\n",slicelength);
	packStart = 0;
	while(1)
	{
	  //dumpPacket(buf,24,16);
          //printf("packStart:%lld | size_of_read:%d\n",packStart,size_of_read);
          if(wfctaDecode->bigPackCheck(buf,int(size_of_read),packStart))
          {
	    //dumpPacket(buf,24,16);
	    wfctaEvent->iTel = ITEL;
	    //printf("ITEL%d:\n",ITEL);
	    //get info eventID and rabbit_time//
	    //printf("packStart:%lld | size_of_read:%d\n",packStart,size_of_read);
	    wfctaEvent->big_pack_lenth = wfctaDecode->bigpackLen();

            wfctaEvent->iEvent=wfctaDecode->eventId(buf);
            wfctaEvent->rabbitTime=wfctaDecode->RabbitTime(buf);
            wfctaEvent->rabbittime=wfctaDecode->Rabbittime(buf);
            wfctaEvent->n_fired = wfctaDecode->nFired(buf);
            printf("iEvent:%d:\n\n",wfctaDecode->eventId(buf));

  	    //find sipms and their position in this pack//
  	    wfctaDecode->Find_SiPMs(buf);//,0);
	    sipm_position = &(wfctaDecode->GetSiPM_Position());

	    //get info of each sipm: q, base, peakposition...//
	    wfctaEvent->n_Channel = 0;
	    for(sipm_position_iter=sipm_position->begin(); sipm_position_iter!=sipm_position->end(); sipm_position_iter++){
	      wfctaEvent->n_Channel++;
	      wfctaEvent->iSiPM.push_back( sipm_position_iter->first );
  	      wfctaEvent->ievent.push_back( wfctaDecode->eventId_in_channel(buf,sipm_position_iter->first) );
	      wfctaEvent->Over_Single_Marker.push_back( wfctaDecode->GetOver_Single_Mark(buf,sipm_position_iter->first) );
	      wfctaEvent->Over_Record_Marker.push_back( wfctaDecode->GetOver_Record_Mark(buf,sipm_position_iter->first) );
	      wfctaEvent->ImageBaseHigh.push_back( wfctaDecode->BaseHigh(buf,sipm_position_iter->first) );
              //printf("%d %d:\n",wfctaDecode->eventId_in_channel(buf,sipm_position_iter->first),sipm_position_iter->first);
	      wfctaDecode->GetWaveForm(buf,sipm_position_iter->first,(int *)(wfctaEvent->pulsehigh), (int *)(wfctaEvent->pulselow));
	      wfctaEvent->mypeak.push_back( wfctaDecode->Getwavepeak(buf,sipm_position_iter->first) );
              wfctaEvent->peakamp.push_back( wfctaDecode->GetpeakAmp(buf,sipm_position_iter->first) );
	      wfctaEvent->myImageBaseHigh.push_back( wfctaDecode->GetwaveImageBaseHigh(buf,sipm_position_iter->first) );
              wfctaEvent->myImageBaseLow.push_back( wfctaDecode->GetwaveImageBaseLow(buf,sipm_position_iter->first) );
              wfctaEvent->myImageAdcHigh.push_back( wfctaDecode->GetwaveImageAdcHigh(buf,sipm_position_iter->first) );
              wfctaEvent->myImageAdcLow.push_back( wfctaDecode->GetwaveImageAdcLow(buf,sipm_position_iter->first) );
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

	    packStart = wfctaDecode->PackSize();
          }
          else
          {
	      break;
          }
        } 
      }
      else
      {
	  fseek(fp,-size_of_read+20,1);
      }
  }
  fclose(fp);

/******************************************************************************/
  eventShow->Write();
  //rootfile->Write();
  rootfile->Close();

}
