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
#include "include/WFCTAMerge.h"
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
		printf("Use %s iptpath/inputfile outpath/outfile\n",argv[0]);
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

	char Name1[300]="root://eos01.ihep.ac.cn/";
	char Name2[300];
	strcpy(Name2,Name1);
	strcat(Name2,argv[2]);
	TFile *rootfile = TFile::Open(Name2,"recreate");
	WFCTAEvent *wfctaEvent = new WFCTAEvent();
	/*********************************************************************/
	TTree *eventShow = new TTree("eventShow","info of evnets");
	wfctaEvent -> CreateBranch(eventShow,1);
	/*********************************************************************/

	WFCTADecode *wfctaDecode = new WFCTADecode();
	//Events Initial//
	wfctaEvent->EventInitial();

	fp = fopen(argv[1],"rb");
	vector<long>* rb_Time = new vector<long>();
	vector<double>* rb_time = new vector<double>();
	vector<long>* pack_pos = new vector<long>();
	vector<long>* pack_len = new vector<long>();
	rb_Time->clear();
	rb_time->clear();
	pack_pos->clear();
	pack_len->clear();
	while(true)
	{
		buf = new uint8_t[40];
		size_of_read = fread((uint8_t *)buf,1,40,fp);
		if(size_of_read==0){break;}
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
			packStart = 0;
			while(1)
			{
				if(wfctaDecode->bigPackCheck(buf,int(size_of_read),packStart))
				{
					//get info rabbit_time and position in file//
					rb_Time->push_back( wfctaDecode->RabbitTime(buf) );
					rb_time->push_back( wfctaDecode->Rabbittime(buf) );
					pack_pos->push_back( FEEPos + packStart );
					pack_len->push_back( wfctaDecode->bigpackLen() );
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
	map<long, long> time_position;
	map<long, long>::iterator time_position_iter;
	map<long, long>::iterator next_time_position_iter;
	for(int ii=0;ii<rb_Time->size();ii++){
		time_position.insert( pair<long,long>( long((rb_Time->at(ii)-rb_TimeMin)*1000000000+rb_time->at(ii)*20), pack_pos->at(ii) ) );
	}

	float adch;
	float adcl;
	long deltaTime;
	int nevent[20]={0};
	short ISIPM;
	WFCTAMerge merge_ev;
	vector<WFCTAMerge> merge_evs;
	merge_evs.clear();

	slicelength = 2000000;
	fp = fopen(argv[1],"rb");
	next_time_position_iter=time_position.begin();
	next_time_position_iter++;
	for(time_position_iter=time_position.begin(); time_position_iter!=time_position.end(); time_position_iter++)
	{
		fseek(fp,time_position_iter->second,0);
		buf = new uint8_t[slicelength];
		fread((uint8_t *)buf,1,slicelength,fp);

		//dumpPacket(buf,36,16);
		merge_ev.EventInitial();
		wfctaDecode->bigPackCheck(buf,int(slicelength),0);
		merge_ev.eEvent=wfctaDecode->eventId(buf);
		merge_ev.rabbitTime=wfctaDecode->RabbitTime(buf);
		merge_ev.rabbittime=wfctaDecode->Rabbittime(buf);
		merge_ev.big_pack_lenth = wfctaDecode->bigpackLen();
		merge_ev.n_fired = wfctaDecode->nFired(buf);
		//find sipms and their position in this pack//
		wfctaDecode->Find_SiPMs(buf);//,0);
		sipm_position = &(wfctaDecode->GetSiPM_Position());
		//get info of each sipm: q, base, peakposition...//
		for(sipm_position_iter=sipm_position->begin(); sipm_position_iter!=sipm_position->end(); sipm_position_iter++){
			merge_ev.n_Channel++;
			ISIPM = sipm_position_iter->first;
			merge_ev.IsData[ISIPM] = 1;
			merge_ev.eevent[ISIPM] = wfctaDecode->eventId_in_channel(buf,sipm_position_iter->first);
			merge_ev.zipmod[ISIPM] = wfctaDecode->zipMode(buf,sipm_position_iter->first);
			merge_ev.Over_Single_Marker[ISIPM] = wfctaDecode->GetOver_Single_Mark(buf,sipm_position_iter->first);
			merge_ev.Over_Record_Marker[ISIPM] = wfctaDecode->GetOver_Record_Mark(buf,sipm_position_iter->first);
			merge_ev.winsum[ISIPM] = wfctaDecode->Getwinsum(buf,sipm_position_iter->first);
			//wfctaDecode->GetWaveForm(buf,sipm_position_iter->first,(int *)(merge_ev.pulsehigh), (int *)(merge_ev.pulselow));
			//wfctaDecode->GeteSaturation(buf,sipm_position_iter->first,(int *)(merge_ev.saturationH), (int *)(merge_ev.saturationL));
			wfctaDecode->GetWaveForm(buf,short(ISIPM),(int *)(merge_ev.pulsehigh), (int *)(merge_ev.pulselow));
			wfctaDecode->GeteSaturation(buf,short(ISIPM),(int *)(merge_ev.saturationH), (int *)(merge_ev.saturationL));
		}

		merge_evs.push_back(merge_ev);
		deltaTime = next_time_position_iter->first - time_position_iter->first;
		printf("TimeNext:%lld    TimeThis:%lld    deltaTime:%lld\n",next_time_position_iter->first,time_position_iter->first,deltaTime);

		if(deltaTime!=1600)
		{
			wfctaEvent->iTel = ITEL;
			wfctaEvent->merge_size = merge_evs.size();
			nevent[ITEL]++;
			wfctaEvent->iEvent=nevent[ITEL];
			wfctaEvent->eEvent=WFCTAMerge::GeteEvent(merge_evs);
			wfctaEvent->rabbitTime=WFCTAMerge::RabbitTime(merge_evs);
			wfctaEvent->rabbittime=WFCTAMerge::Rabbittime(merge_evs);
			wfctaEvent->big_pack_lenth=WFCTAMerge::GetBigPackLen(merge_evs);
			wfctaEvent->n_fired=WFCTAMerge::GetNFired(merge_evs);
			wfctaEvent->n_Channel=WFCTAMerge::GetNChannel(merge_evs);
			for(int isipm=0;isipm<1024;isipm++)
			{
				if(!WFCTAMerge::IsData_Merge(isipm,merge_evs)){continue;}
				if(ITEL==5)	{wfctaEvent->iSiPM.push_back( 1023 - isipm );}
				else		{wfctaEvent->iSiPM.push_back( isipm );}
				wfctaEvent->eevent.push_back( WFCTAMerge::eevent_Merge(isipm,merge_evs) );
				wfctaEvent->zipmod.push_back( WFCTAMerge::zipmod_Merge(isipm,merge_evs) );
				wfctaEvent->Over_Single_Marker.push_back( WFCTAMerge::OvSigMarker_Merge(isipm,merge_evs) );
				wfctaEvent->Over_Record_Marker.push_back( WFCTAMerge::OvRecMarker_Merge(isipm,merge_evs) );
				wfctaEvent->winsum.push_back( WFCTAMerge::WimSum_Merge(isipm,merge_evs) );
				wfctaEvent->eSatH.push_back( WFCTAMerge::eSatH_Merge(isipm,merge_evs) );
				wfctaEvent->eSatL.push_back( WFCTAMerge::eSatL_Merge(isipm,merge_evs) );
				wfctaEvent->PeakPosH.push_back( WFCTAMerge::GetPeakPosH(isipm,merge_evs) );
				wfctaEvent->PeakPosL.push_back( WFCTAMerge::GetPeakPosL(isipm,merge_evs) );
				wfctaEvent->PeakAmH.push_back( WFCTAMerge::GetPeakAmpH(isipm,merge_evs) );
				wfctaEvent->PeakAmL.push_back( WFCTAMerge::GetPeakAmpL(isipm,merge_evs) );
				WFCTAMerge::Calc_Q_Base(isipm,merge_evs,0);
				wfctaEvent->BaseH.push_back( WFCTAMerge::GetBaseH(isipm,merge_evs) );
				wfctaEvent->BaseL.push_back( WFCTAMerge::GetBaseL(isipm,merge_evs) );
				adch = WFCTAMerge::GetAdcH(isipm,merge_evs);
				adcl = WFCTAMerge::GetAdcL(isipm,merge_evs);
				wfctaEvent->AdcH.push_back( adch );
				wfctaEvent->AdcL.push_back( adcl );
				if(adch>6000){  wfctaEvent->SatH.push_back(1);}
				else         {  wfctaEvent->SatH.push_back(0);}
				if(adcl>6000){    wfctaEvent->SatL.push_back(1);}
				else         {    wfctaEvent->SatL.push_back(0);}
				WFCTAMerge::Calc_Q_Base(isipm,merge_evs,1);
				wfctaEvent->LaserBaseH.push_back( WFCTAMerge::GetLaserBaseH(isipm,merge_evs) );
				wfctaEvent->LaserBaseL.push_back( WFCTAMerge::GetLaserBaseL(isipm,merge_evs) );
				wfctaEvent->LaserAdcH.push_back( WFCTAMerge::GetLaserAdcH(isipm,merge_evs) );
				wfctaEvent->LaserAdcL.push_back( WFCTAMerge::GetLaserAdcL(isipm,merge_evs) );
				if(isipm==1){
					printf("winsum_merge:%f\n",WFCTAMerge::WimSum_Merge(isipm,merge_evs));
				}
			}
			printf("merge event:%d\n\n",merge_evs.size());
			eventShow->Fill();
			merge_evs.clear();
			wfctaEvent->EventInitial();
		}

		delete[] buf;
		if(next_time_position_iter!=time_position.end())
		{
			next_time_position_iter++;
		}
	} 
	fclose(fp);
	/******************************************************************************/
	rootfile->Write();
	rootfile->Close();

}
