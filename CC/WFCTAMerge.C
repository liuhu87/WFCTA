#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include "WFCTAMerge.h"

using namespace std;

int32_t WFCTAMerge::peakAmpH=0;
int32_t WFCTAMerge::peakAmpL=0;
uint8_t WFCTAMerge::peakPosH=0;
uint8_t WFCTAMerge::peakPosL=0;
float WFCTAMerge::m_Basehigh=0;
float WFCTAMerge::m_Baselow=0;
float WFCTAMerge::m_Adchigh=0;
float WFCTAMerge::m_Adclow=0;
vector<int> WFCTAMerge::merged_pulsehigh(0,0);
vector<int> WFCTAMerge::merged_pulselow(0,0);

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
		IsData[i] = 0;
		eevent[i] = -1;
		zipmod[i] = -1;
		winsum[i] = -1;
		Over_Single_Marker[i] = 0;
		Over_Record_Marker[i] = 0;
	}
	for(int j=0;j<28;j++){
		for(int i=0;i<1024;i++){
			saturationH[i][j] = 0;
			saturationL[i][j] = 0;
			pulsehigh[i][j] = 0;
			pulselow[i][j] = 0;
		}
	}

}

long WFCTAMerge::GeteEvent(vector<WFCTAMerge> &evs)
{
	long e_Evt=-1;
	vector<WFCTAMerge>::iterator evs_iter;
	evs_iter=evs.begin();
	e_Evt = (*evs_iter).eEvent;
	return e_Evt;
}

long WFCTAMerge::RabbitTime(vector<WFCTAMerge> &evs)
{
	long rb_Time=0;
	vector<WFCTAMerge>::iterator evs_iter;
	evs_iter=evs.begin();
	rb_Time = (*evs_iter).rabbitTime;
	return rb_Time;
}
long WFCTAMerge::Rabbittime(vector<WFCTAMerge> &evs)
{
	long rb_time=0;
	vector<WFCTAMerge>::iterator evs_iter;
	evs_iter=evs.begin();
	rb_time = (*evs_iter).rabbittime;
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
	int nchannel[1024] = {0};
	for(int isipm=0;isipm<1024;isipm++)
	{
		for(evs_iter=evs.begin(); evs_iter!=evs.end(); evs_iter++)
		{
			nchannel[isipm] += (*evs_iter).IsData[isipm];
		}
		if(nchannel[isipm]!=0)
			NChannel++;
	}
	return NChannel;
}

bool WFCTAMerge::IsData_Merge(int isipm, vector<WFCTAMerge> &evs)
{
	bool isdata = 0;
	vector<WFCTAMerge>::iterator evs_iter;
	for(evs_iter=evs.begin(); evs_iter!=evs.end(); evs_iter++)
	{
		if((*evs_iter).IsData[isipm]==1)
			isdata = 1;
	}
	return isdata;
}

long WFCTAMerge::eevent_Merge(int isipm, vector<WFCTAMerge> &evs)
{
	long e_evt;
	vector<WFCTAMerge>::iterator evs_iter;
	evs_iter=evs.begin();
	e_evt = (*evs_iter).eevent[isipm];
	return e_evt;
}


short WFCTAMerge::zipmod_Merge(int isipm, vector<WFCTAMerge> &evs)
{
	short z_mod;
	vector<WFCTAMerge>::iterator evs_iter;
	evs_iter=evs.begin();
	z_mod = (*evs_iter).zipmod[isipm];
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

int WFCTAMerge::eSatH_Merge(int isipm, vector<WFCTAMerge> &evs)
{
	int merged_eSath = 0;
	vector<WFCTAMerge>::iterator evs_iter;
	for(evs_iter=evs.begin(); evs_iter!=evs.end(); evs_iter++)
	{
		for(int ipoint=0;ipoint<28;ipoint++)
		{
			if((*evs_iter).saturationH[isipm][ipoint]==1)
				merged_eSath = 1;
		}
	}
	return merged_eSath;
}
int WFCTAMerge::eSatL_Merge(int isipm, vector<WFCTAMerge> &evs)
{
	int merged_eSatl = 0;
	vector<WFCTAMerge>::iterator evs_iter;
	for(evs_iter=evs.begin(); evs_iter!=evs.end(); evs_iter++)
	{
		for(int ipoint=0;ipoint<28;ipoint++)
		{
			if((*evs_iter).saturationL[isipm][ipoint]==1)
				merged_eSatl = 1;
		}
	}
	return merged_eSatl;
}


char WFCTAMerge::GetPeakPosH(int isipm, vector<WFCTAMerge> &evs)
{
	WFCTAMerge::FindPeak(isipm,evs);
	if(isipm==1)
	{
		printf("merged_waveform_h: ");
		for(int ii=0;ii<merged_pulsehigh.size();ii++){
			printf("%d ",merged_pulsehigh.at(ii));
		}
		printf("\n");
		printf("merged_waveform_l: ");
		for(int ii=0;ii<merged_pulsehigh.size();ii++){
			printf("%d ",merged_pulselow.at(ii));
		}
		printf("\n");
	}
	return peakPosH;
}

char WFCTAMerge::GetPeakPosL(int isipm, vector<WFCTAMerge> &evs)
{
	WFCTAMerge::FindPeak(isipm,evs);
	return peakPosL;
}
int WFCTAMerge::GetPeakAmpH(int isipm, vector<WFCTAMerge> &evs)
{
	WFCTAMerge::FindPeak(isipm,evs);
	return peakAmpH;
}
int WFCTAMerge::GetPeakAmpL(int isipm, vector<WFCTAMerge> &evs)
{
	WFCTAMerge::FindPeak(isipm,evs);
	return peakAmpL;
}
float WFCTAMerge::GetBaseH(int isipm, vector<WFCTAMerge> &evs)
{
	WFCTAMerge::Calc_Q_Base(isipm,evs);
	return m_Basehigh;
}
float WFCTAMerge::GetBaseL(int isipm, vector<WFCTAMerge> &evs)
{
	WFCTAMerge::Calc_Q_Base(isipm,evs);
	return m_Baselow;
}
float WFCTAMerge::GetAdcH(int isipm, vector<WFCTAMerge> &evs)
{
	WFCTAMerge::Calc_Q_Base(isipm,evs);
	return m_Adchigh;
}
float WFCTAMerge::GetAdcL(int isipm, vector<WFCTAMerge> &evs)
{
	WFCTAMerge::Calc_Q_Base(isipm,evs);
	return m_Adclow;
}


void WFCTAMerge::Calc_Q_Base(int isipm, vector<WFCTAMerge> &evs)
{
	int wave_len = merged_pulsehigh.size();
	WFCTAMerge::FindPeak(isipm,evs);
	m_Basehigh = 0;
	m_Baselow = 0;
	m_Adchigh = 0;
	m_Adclow = 0;

	if(peakPosH<3)
	{
		for(int i=6;i<wave_len;i++)  { m_Basehigh += merged_pulsehigh.at(i);}
		for(int i=0;i<6;i++)		 { m_Adchigh += merged_pulsehigh.at(i); }
	}
	else if(peakPosH>wave_len-5)
	{
		for(int i=0;i<wave_len-6;i++)		 { m_Basehigh += merged_pulsehigh.at(i);}
		for(int i=wave_len-6;i<wave_len;i++) { m_Adchigh += merged_pulsehigh.at(i); }
	}
	else
	{
		for(int i=0;i<peakPosH-2;i++)           { m_Basehigh += merged_pulsehigh.at(i);}
		for(int i=peakPosH+4;i<wave_len;i++)	{ m_Basehigh += merged_pulsehigh.at(i);}
		for(int i=peakPosH-2;i<peakPosH+4;i++)  { m_Adchigh += merged_pulsehigh.at(i); }
	}

	if(peakPosL<3)
	{
		for(int i=6;i<wave_len;i++)  { m_Baselow += merged_pulselow.at(i);}
		for(int i=0;i<6;i++)		 { m_Adclow += merged_pulselow.at(i); }
	}
	else if(peakPosL>wave_len-5)
	{
		for(int i=0;i<wave_len-6;i++)		 { m_Baselow += merged_pulselow.at(i);}
		for(int i=wave_len-6;i<wave_len;i++) { m_Adclow += merged_pulselow.at(i); }
	}
	else
	{
		for(int i=0;i<peakPosL-2;i++)			{ m_Baselow += merged_pulselow.at(i);}
		for(int i=peakPosL+4;i<wave_len;i++)    { m_Baselow += merged_pulselow.at(i);}
		for(int i=peakPosL-2;i<peakPosL+4;i++)  { m_Adclow += merged_pulselow.at(i); }
	}

	m_Basehigh = m_Basehigh/(4*(wave_len-6));
	m_Baselow = m_Baselow/(4*(wave_len-6));
	m_Adchigh -= m_Basehigh*24;
	m_Adclow -= m_Baselow*24;
}

void WFCTAMerge::FindPeak(int isipm, vector<WFCTAMerge> &evs)
{
	WFCTAMerge::WaveForm_Merge(isipm,evs);
	double sumhighmax = -1000;
	double sumhigh;
	for(int ii=0;ii<merged_pulsehigh.size();ii++){
		sumhigh = merged_pulsehigh.at(ii);
		if(sumhighmax<sumhigh) {sumhighmax = sumhigh; peakPosH = ii; peakAmpH = sumhigh;}
	}
	double sumlowmax = -1000;
	double sumlow;
	for(int ii=0;ii<merged_pulsehigh.size();ii++){
		sumlow = merged_pulselow.at(ii);
		if(sumlowmax<sumlow) {sumlowmax = sumlow; peakPosL = ii; peakAmpL = sumlow;}
	}
}

void WFCTAMerge::WaveForm_Merge(int isipm, vector<WFCTAMerge> &evs)
{
	vector<WFCTAMerge>::iterator evs_iter;
	merged_pulsehigh.clear();
	merged_pulselow.clear();
	int first_event=0;
	for(evs_iter=evs.begin(); evs_iter!=evs.end(); evs_iter++)
	{
		if(!first_event)
		{
			for(int ipoint=0;ipoint<28;ipoint++)
			{
				if((*evs_iter).pulsehigh[isipm][ipoint]==0){continue;}
				merged_pulsehigh.push_back((*evs_iter).pulsehigh[isipm][ipoint]);
				merged_pulselow.push_back((*evs_iter).pulselow[isipm][ipoint]);
				first_event++;
			}
		}
		else
		{
			for(int ipoint=8;ipoint<28;ipoint++)
			{
				if((*evs_iter).pulsehigh[isipm][ipoint]==0){continue;}
				merged_pulsehigh.push_back((*evs_iter).pulsehigh[isipm][ipoint]);
				merged_pulselow.push_back((*evs_iter).pulselow[isipm][ipoint]);
			}
		}

	/*	
		if(isipm==1)
		{
			printf("%02d:waveform_h:     ",first_event);
			for(int ipoint=0;ipoint<28;ipoint++)
			{
				printf("%d ",(*evs_iter).pulsehigh[isipm][ipoint]);
			}
			printf("\n");
			printf("%02d:waveform_l:     ",first_event);
			for(int ipoint=0;ipoint<28;ipoint++)
			{
				printf("%d ",(*evs_iter).pulselow[isipm][ipoint]);
			}
			printf("\n");
		}
		*/
		
	}
}
