#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include "WFCTAMerge.h"

using namespace std;

WFCTAMerge::WFCTAMerge()
{
	EventInitial();
}

WFCTAMerge::~WFCTAMerge()
{
	EventInitial();
}

void WFCTAMerge::EventInitial()
{
	big_pack_lenth=-1;
	eEvent=-1;
	rabbitTime=0;
	rabbittime=0;
	n_fired=-1;
	n_Channel=0;
	for(int i=0;i<1024;i++)
	{
		eevent[i] = -1;
		zipmod[i] = -1;
		winsum[i] = 0;
		Over_Single_Marker[i] = 0;
		Over_Record_Marker[i] = 0;
	}
	for(int j=0;j<28;j++){
		for(int i=0;i<1024;i++){
			pulsehigh[i][j] = 0;
			pulselow[i][j] = 0;
		}
	}

}

long WFCTAMerge::GeteEvent(vector<WFCTAMerge> &evs)
{
	long rb_Time=0;
	long e_Evt;
	vector<WFCTAMerge>::iterator evs_iter;
	for(evs_iter=evs.begin(); evs_iter!=evs.end(); evs_iter++)
	{
		if(rb_Time<(*evs_iter).rabbitTime)
			e_Evt = (*evs_iter).eEvent;
	}
	return e_Evt;
}

long WFCTAMerge::RabbitTime(vector<WFCTAMerge> &evs)
{
	long rb_Time=0;
	vector<WFCTAMerge>::iterator evs_iter;
	for(evs_iter=evs.begin(); evs_iter!=evs.end(); evs_iter++)
	{
		if(rb_Time<(*evs_iter).rabbitTime)
			rb_Time = (*evs_iter).rabbitTime;
	}
	return rb_Time;
}
long WFCTAMerge::Rabbittime(vector<WFCTAMerge> &evs)
{
	long rb_time=0;
	vector<WFCTAMerge>::iterator evs_iter;
	for(evs_iter=evs.begin(); evs_iter!=evs.end(); evs_iter++)
	{
		if(rb_time<(*evs_iter).rabbittime)
			rb_time = (*evs_iter).rabbittime;
	}
	return rb_time;
}


int WFCTAMerge::GetBigPackLen(vector<WFCTAMerge> &evs)
{
	int big_pac_len=0;
	vector<WFCTAMerge>::iterator evs_iter;
	for(evs_iter=evs.begin(); evs_iter!=evs.end(); evs_iter++)
	{
		big_pac_len += (*evs_iter).big_pack_lenth;
	}
	return big_pac_len;
}

short WFCTAMerge::GetNFired(vector<WFCTAMerge> &evs)
{
	short NFired=0;
	vector<WFCTAMerge>::iterator evs_iter;
	for(evs_iter=evs.begin(); evs_iter!=evs.end(); evs_iter++)
	{
		NFired += (*evs_iter).n_fired;
	}
	return NFired;
}


short WFCTAMerge::GetNChannel(vector<WFCTAMerge> &evs)
{
	short NChannel=0;
	vector<WFCTAMerge>::iterator evs_iter;
	for(evs_iter=evs.begin(); evs_iter!=evs.end(); evs_iter++)
	{
		NChannel += (*evs_iter).n_Channel;
	}
	return NChannel;
}

long WFCTAMerge::eevent_Merge(int isipm, vector<WFCTAMerge> &evs)
{
	long rb_Time=0;
	long e_evt;
	vector<WFCTAMerge>::iterator evs_iter;
	for(evs_iter=evs.begin(); evs_iter!=evs.end(); evs_iter++)
	{
		if(rb_Time<(*evs_iter).rabbitTime)
			e_evt = (*evs_iter).eevent[isipm];
	}
	return e_evt;
}


short WFCTAMerge::zipmod_Merge(int isipm, vector<WFCTAMerge> &evs)
{
	long rb_Time=0;
	short z_mod;
	vector<WFCTAMerge>::iterator evs_iter;
	for(evs_iter=evs.begin(); evs_iter!=evs.end(); evs_iter++)
	{
		if(rb_Time<(*evs_iter).rabbitTime)
			z_mod = (*evs_iter).zipmod[isipm];
	}
	return z_mod;
}

bool WFCTAMerge::OvSigMarker_Merge(int isipm, vector<WFCTAMerge> &evs)
{
	bool merged_OvSigMarker = 0;
	vector<WFCTAMerge>::iterator evs_iter;
	for(evs_iter=evs.begin(); evs_iter!=evs.end(); evs_iter++)
	{
		if((*evs_iter).Over_Single_Marker[isipm]==1)
			merged_OvSigMarker = 1;
	}
	return merged_OvSigMarker;
}

bool WFCTAMerge::OvRecMarker_Merge(int isipm, vector<WFCTAMerge> &evs)
{
	bool merged_OvRecMarker = 0;
	vector<WFCTAMerge>::iterator evs_iter;
	for(evs_iter=evs.begin(); evs_iter!=evs.end(); evs_iter++)
	{
		if((*evs_iter).Over_Record_Marker[isipm]==1)
			merged_OvRecMarker = 1;
	}
	return merged_OvRecMarker;
}


float WFCTAMerge::WimSum_Merge(int isipm, vector<WFCTAMerge> &evs)
{
	float merged_winsum = -1000;
	vector<WFCTAMerge>::iterator evs_iter;
	for(evs_iter=evs.begin(); evs_iter!=evs.end(); evs_iter++)
	{
		if(merged_winsum<(*evs_iter).winsum[isipm])
			merged_winsum = (*evs_iter).winsum[isipm];
	}
	return merged_winsum;
}

void WFCTAMerge::WaveForm_Merge(int isipm, vector<WFCTAMerge> &evs)
{
	vector<WFCTAMerge>::iterator evs_iter;
	for(evs_iter=evs.begin(); evs_iter!=evs.end(); evs_iter++)
	{
		/*
		printf("waveform: ");
		for(int ipoint=0;ipoint<28;ipoint++)
		{
			printf("%d ",(*evs_iter).pulsehigh[isipm][ipoint]);
		}
		printf("\n");
		*/
	}
}
