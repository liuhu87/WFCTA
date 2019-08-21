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
  int32_t big_pack_lenth;
  //int16_t n_fired; 
  int n_Channel;
  int Npoint[28] = {0};  for(int i=0;i<28;i++) {Npoint[i]=i;}
  vector<int>* peakamp = new vector<int>();
  int pulsehigh[1024][28];
  int pulselow[1024][28];

  WFCTAEvent *wfctaEvent = new WFCTAEvent();
  TFile *rootfile = new TFile(argv[2],"recreate");
  /*********************************************************************/
  TTree *eventShow = new TTree("eventShow","info of evnets");
  wfctaEvent -> CreateBranch(eventShow,1);
  eventShow -> Branch("big_pack_lenth",&big_pack_lenth,"big_pack_lenth/I");
  //eventShow -> Branch("n_fired",&n_fired,"n_fired/S");
  eventShow -> Branch("n_Channel",&n_Channel,"n_Channel/I");
  eventShow -> Branch("peakamp","vector<int>",&peakamp);
  eventShow -> Branch("Npoint",Npoint,"Npoint[28]/I");
  eventShow -> Branch("pulsehigh",pulsehigh,"pulsehigh[1024][28]/I");
  eventShow -> Branch("pulselow",pulselow,"pulselow[1024][28]/I");
  /*********************************************************************/

  WFCTADecode *wfctaDecode = new WFCTADecode();

  //Events Initial//
  wfctaEvent->EventInitial();
  peakamp->clear();
  for(int i=0;i<1024;i++){
    for(int j=0;j<28;j++){
      pulsehigh[i][j] = 0;
      pulselow[i][j] = 0;
    }
  }

  fp = fopen(argv[1],"rb");
  while(true)
  {
      size_of_read = fread((uint8_t *)buf,1,BUF_LEN,fp);
      fseek(fp,-size_of_read,1);
      if(size_of_read==0){break;}
      //dumpPacket(buf,size_of_read,16);
//printf("test\n");
      if(wfctaDecode->bigPackCheck(buf,int(size_of_read)))
      {
	  //get info eventID and rabbit_time//
          packSize = wfctaDecode->PackSize();
	  big_pack_lenth = wfctaDecode->bigpackLen();
          printf("%d:\n",wfctaDecode->eventId(buf));

          wfctaEvent->iEvent=wfctaDecode->eventId(buf);
          wfctaEvent->rabbitTime=wfctaDecode->RabbitTime(buf);
          wfctaEvent->rabbittime=wfctaDecode->Rabbittime(buf);
          wfctaEvent->n_fired = wfctaDecode->nFired(buf);

	  //find sipms and their position in this pack//
	  wfctaDecode->Find_SiPMs(buf,packSize);
	  sipm_position = &(wfctaDecode->GetSiPM_Position());

	  //get info of each sipm: q, base, peakposition...//
	  n_Channel = 0;
	  for(sipm_position_iter=sipm_position->begin(); sipm_position_iter!=sipm_position->end(); sipm_position_iter++){
	    n_Channel++;
	    wfctaEvent->iSiPM.push_back( sipm_position_iter->first );
	    wfctaEvent->ievent.push_back( wfctaDecode->eventId_in_channel(buf,sipm_position_iter->first) );
	    wfctaEvent->Over_Single_Marker.push_back( wfctaDecode->GetOver_Single_Mark(buf,sipm_position_iter->first) );
	    wfctaEvent->Over_Record_Marker.push_back( wfctaDecode->GetOver_Record_Mark(buf,sipm_position_iter->first) );
	    wfctaEvent->ImageBaseHigh.push_back( wfctaDecode->BaseHigh(buf,sipm_position_iter->first) );
            //printf("%d %d:\n",wfctaDecode->eventId_in_channel(buf,sipm_position_iter->first),sipm_position_iter->first);
	    wfctaDecode->GetWaveForm(buf,sipm_position_iter->first,(int *)pulsehigh, (int *)pulselow);
	    wfctaEvent->mypeak.push_back( wfctaDecode->Getwavepeak(buf,sipm_position_iter->first) );
            peakamp->push_back( wfctaDecode->GetpeakAmp(buf,sipm_position_iter->first) );
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
	  peakamp->clear();
	  for(int i=0;i<1024;i++){
            for(int j=0;j<28;j++){
              pulsehigh[i][j] = 0;
              pulselow[i][j] = 0;
            }
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
