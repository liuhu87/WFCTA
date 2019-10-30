#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <vector>
#include <algorithm>
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
	int32_t slicelength = 2000000;
	size_t size_of_read;
	short ITEL;
	int FEEDataHead;
	long FEEPos;

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
	int nevent[20]={0};
	vector<long>* rb_Time = new vector<long>();
	vector<double>* rb_time = new vector<double>();
	vector<long>* sort_time = new vector<long>();
	vector<long>* pack_pos = new vector<long>();
	vector<long>* pack_len = new vector<long>();
	vector<long>* iEVT = new vector<long>();
	vector<long>* eEVT = new vector<long>();
	rb_Time->clear();
	rb_time->clear();
	sort_time->clear();
	pack_pos->clear();
	pack_len->clear();
	iEVT->clear();
	eEVT->clear();
	while(true)
	{
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
			FEEPos = ftell(fp);

			delete[] buf;
			buf = new uint8_t[slicelength];
			size_of_read = fread((uint8_t *)buf,1,slicelength,fp);
			//printf("slicelength:%lld\n",slicelength);
			packStart = 0;
			while(1)
			{
				if(wfctaDecode->bigPackCheck(buf,int(size_of_read),packStart))
				{
					wfctaEvent->iTel = ITEL;
					//get info eventID and rabbit_time//
					wfctaEvent->big_pack_lenth = wfctaDecode->bigpackLen();
					nevent[ITEL]++;
					wfctaEvent->iEvent=nevent[ITEL];
					wfctaEvent->eEvent=wfctaDecode->eventId(buf);
					wfctaEvent->rabbitTime=wfctaDecode->RabbitTime(buf);
					wfctaEvent->rabbittime=wfctaDecode->Rabbittime(buf);
					//printf("iEvent:%d Time:%lld time:%lf\n\n",wfctaEvent->iEvent,wfctaEvent->rabbitTime,wfctaEvent->rabbittime);

					rb_Time->push_back( wfctaEvent->rabbitTime);
					rb_time->push_back( wfctaEvent->rabbittime);
					pack_pos->push_back( FEEPos + packStart );
					pack_len->push_back( wfctaEvent->big_pack_lenth );
					iEVT->push_back( nevent[ITEL] );
					eEVT->push_back( wfctaDecode->eventId(buf) );

					//eventShow->Fill();
					wfctaEvent->EventInitial();

					packStart = wfctaDecode->PackSize();
				}
				else
				{
					break;
				}
			} 
			delete[] buf;
		}
		else
		{
			delete[] buf;
			fseek(fp,-size_of_read+20,1);
		}
	}
	fclose(fp);

	long rb_TimeMin = *min_element(rb_Time->begin(),rb_Time->end());
	printf("TimeMin:%lld\n\n",rb_TimeMin);
	for(int ii=0;ii<rb_Time->size();ii++){
		sort_time->push_back( long((rb_Time->at(ii)-rb_TimeMin)*1000000000+rb_time->at(ii)) );
	}
	map<long, long> time_position;// = nullptr;
	map<long, long>::iterator time_position_iter;
	for(int ii=0;ii<sort_time->size();ii++){
		time_position.insert(pair<long,long>(sort_time->at(ii),pack_pos->at(ii)));
	}

	float adch;
	float adcl;
	fp = fopen(argv[1],"rb");
	for(time_position_iter=time_position.begin(); time_position_iter!=time_position.end(); time_position_iter++)
	{
		printf("Time:%lld pos:%lld\n\n",time_position_iter->first,time_position_iter->second);
		fseek(fp,time_position_iter->second,0);
		buf = new uint8_t[slicelength];
		fread((uint8_t *)buf,1,slicelength,fp);
		
		//dumpPacket(buf,48,16);
		//printf("packStart:%lld | size_of_read:%d\n",packStart,size_of_read);
		if(wfctaDecode->bigPackCheck(buf,int(slicelength),0))
		{
			//dumpPacket(buf,24,16);
			wfctaEvent->iTel = ITEL;
			//printf("ITEL%d:\n",ITEL);
			//get info eventID and rabbit_time//
			//printf("packStart:%lld | size_of_read:%d\n",packStart,size_of_read);
			wfctaEvent->big_pack_lenth = wfctaDecode->bigpackLen();
			nevent[ITEL]++;
			wfctaEvent->iEvent=nevent[ITEL];
			wfctaEvent->eEvent=wfctaDecode->eventId(buf);
			wfctaEvent->rabbitTime=wfctaDecode->RabbitTime(buf);
			wfctaEvent->rabbittime=wfctaDecode->Rabbittime(buf);
			wfctaEvent->n_fired = wfctaDecode->nFired(buf);
			printf("iEvent:%d:\n\n",wfctaEvent->iEvent);

			//find sipms and their position in this pack//
			wfctaDecode->Find_SiPMs(buf);//,0);
			sipm_position = &(wfctaDecode->GetSiPM_Position());

			//get info of each sipm: q, base, peakposition...//
			wfctaEvent->n_Channel = 0;
			for(sipm_position_iter=sipm_position->begin(); sipm_position_iter!=sipm_position->end(); sipm_position_iter++){
				wfctaEvent->n_Channel++;
				if(ITEL==5)
				{
					wfctaEvent->iSiPM.push_back( 1023 - (sipm_position_iter->first) );
				}
				else
				{
					wfctaEvent->iSiPM.push_back( sipm_position_iter->first );
				}
				wfctaEvent->eevent.push_back( wfctaDecode->eventId_in_channel(buf,sipm_position_iter->first) );
				wfctaEvent->zipmod.push_back( wfctaDecode->zipMode(buf,sipm_position_iter->first) );
				wfctaEvent->Over_Single_Marker.push_back( wfctaDecode->GetOver_Single_Mark(buf,sipm_position_iter->first) );
				wfctaEvent->Over_Record_Marker.push_back( wfctaDecode->GetOver_Record_Mark(buf,sipm_position_iter->first) );
				wfctaEvent->winsum.push_back( wfctaDecode->Getwinsum(buf,sipm_position_iter->first) );
				wfctaDecode->GetWaveForm(buf,sipm_position_iter->first,(int *)(wfctaEvent->pulsehigh), (int *)(wfctaEvent->pulselow));
				wfctaEvent->PeakPosH.push_back( wfctaDecode->GetPeakPosH(buf,sipm_position_iter->first) );
				wfctaEvent->PeakPosL.push_back( wfctaDecode->GetPeakPosL(buf,sipm_position_iter->first) );
				wfctaEvent->PeakAmH.push_back( wfctaDecode->GetPeakAmH(buf,sipm_position_iter->first) );
				wfctaEvent->PeakAmL.push_back( wfctaDecode->GetPeakAmL(buf,sipm_position_iter->first) );
				wfctaEvent->BaseH.push_back( wfctaDecode->GetwaveImageBaseHigh(buf,sipm_position_iter->first) );
				wfctaEvent->BaseL.push_back( wfctaDecode->GetwaveImageBaseLow(buf,sipm_position_iter->first) );
				adch = wfctaDecode->GetwaveImageAdcHigh(buf,sipm_position_iter->first);
				adcl = wfctaDecode->GetwaveImageAdcLow(buf,sipm_position_iter->first);
				wfctaEvent->AdcH.push_back( adch );
				wfctaEvent->AdcL.push_back( adcl );
				wfctaEvent->eSatH.push_back( wfctaDecode->eSaturationHigh(buf,sipm_position_iter->first) );
				wfctaEvent->eSatL.push_back( wfctaDecode->eSaturationLow(buf,sipm_position_iter->first) );
				if(adch>6000){	wfctaEvent->SatH.push_back(1);}
				else         {	wfctaEvent->SatH.push_back(0);}
				if(adcl>6000){    wfctaEvent->SatL.push_back(1);}
				else         {    wfctaEvent->SatL.push_back(0);}
			}

			eventShow->Fill();
			wfctaEvent->EventInitial();

			packStart = wfctaDecode->PackSize();
		}
		//else
		//{
		//	break;
		//}
		
		delete[] buf;
	} 
	fclose(fp);
	/******************************************************************************/
	rootfile->Write();
	rootfile->Close();

}
