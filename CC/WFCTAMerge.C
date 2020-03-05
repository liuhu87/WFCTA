#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include "WFCTAMerge.h"

using namespace std;

short WFCTAMerge::emptyPulse=0;
int32_t WFCTAMerge::peakAmpH=0;
int32_t WFCTAMerge::peakAmpL=0;
int16_t WFCTAMerge::peakPosH=0;
int16_t WFCTAMerge::peakPosL=0;
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

/*
void WFCTAMerge::packCheck_Merge(vector<WFCTAMerge> &evs)
{
	vector<short> pack_check;
	pack_check.clear();
	vector<WFCTAMerge>::iterator evs_iter;
	for(evs_iter=evs.begin(); evs_iter!=evs.end(); evs_iter++)
	{
		pack_check.push_back( (*evs_iter).packCheck );
	}
	//return &pack_check;
}
*/
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


short WFCTAMerge::GetPeakPosH(int isipm, vector<WFCTAMerge> &evs)
{
	WFCTAMerge::FindPeak(isipm,evs);
	/*if(isipm==1)
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
	}*/
	short peak_pos_H = short(peakPosH+(emptyPulse/28)*20);
	//if(isipm==1){
	//	printf("peak_pos_H:%d emptyPulse:%d shift:%d\n",peak_pos_H,emptyPulse,(emptyPulse/28)*20);
	//}
	return peak_pos_H;
}

short WFCTAMerge::GetPeakPosL(int isipm, vector<WFCTAMerge> &evs)
{
	WFCTAMerge::FindPeak(isipm,evs);
	short peak_pos_L = short(peakPosL+(emptyPulse/28)*20);
	return peak_pos_L;
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
	//WFCTAMerge::Calc_Q_Base(isipm,evs,laserCalc);
	return m_Basehigh;
}
float WFCTAMerge::GetBaseL(int isipm, vector<WFCTAMerge> &evs)
{
	//WFCTAMerge::Calc_Q_Base(isipm,evs,laserCalc);
	return m_Baselow;
}
float WFCTAMerge::GetAdcH(int isipm, vector<WFCTAMerge> &evs)
{
	//WFCTAMerge::Calc_Q_Base(isipm,evs,laserCalc);
	return m_Adchigh;
}
float WFCTAMerge::GetAdcL(int isipm, vector<WFCTAMerge> &evs)
{
	//WFCTAMerge::Calc_Q_Base(isipm,evs,laserCalc);
	return m_Adclow;
}

float WFCTAMerge::GetLaserBaseH(int isipm, vector<WFCTAMerge> &evs)
{
	//WFCTAMerge::Calc_Q_Base(isipm,evs,laserCalc);
	return m_Basehigh;
}
float WFCTAMerge::GetLaserBaseL(int isipm, vector<WFCTAMerge> &evs)
{
	//WFCTAMerge::Calc_Q_Base(isipm,evs,laserCalc);
	return m_Baselow;
}
float WFCTAMerge::GetLaserAdcH(int isipm, vector<WFCTAMerge> &evs)
{
	//WFCTAMerge::Calc_Q_Base(isipm,evs,laserCalc);
	return m_Adchigh;
}
float WFCTAMerge::GetLaserAdcL(int isipm, vector<WFCTAMerge> &evs)
{
	//WFCTAMerge::Calc_Q_Base(isipm,evs,laserCalc);
	return m_Adclow;
}

void WFCTAMerge::Calc_Q_Base(int isipm, vector<WFCTAMerge> &evs, int laserCalc)
{
	int waveStart,waveEnd;
	int wave_calc_len;
	int base_calc_len;

	int wave_len = merged_pulsehigh.size();
	WFCTAMerge::FindPeak(isipm,evs);
	//if(isipm==1){
	//	printf("peakPosH:%d peakAmH:%d\n",peakPosH,peakAmpH);
	//}
	m_Basehigh = 0;
	m_Baselow = 0;
	m_Adchigh = 0;
	m_Adclow = 0;

	if(laserCalc==0)
	{
		//calc high gain cosmic ray adc and base
		base_calc_len=0;
		wave_calc_len=0;
		if(peakPosH>0)	{	waveStart = peakPosH-1;}
		else			{	waveStart = 0;}
		if(peakPosH<wave_len-2)	{	waveEnd = peakPosH+2;}
		else					{	waveEnd = wave_len-1;}
		for(int i=waveStart;i<=waveEnd;i++){
			m_Adchigh += merged_pulsehigh.at(i);
			wave_calc_len++;
		}
		for(int i=0;i<waveStart-2;i++){
			m_Basehigh += merged_pulsehigh.at(i);
			base_calc_len++;
		}
		for(int i=waveEnd+3;i<wave_len;i++){
			m_Basehigh += merged_pulsehigh.at(i);
			base_calc_len++;
		}
		m_Basehigh = m_Basehigh/(4*base_calc_len);
		m_Adchigh -= m_Basehigh*4*wave_calc_len;
		//if(isipm==1){
		//	printf("BaseH:%f AdcH:%f waveStart:%d waveEnd:%d wave_calc_len:%d base_calc_len:%d\n",m_Basehigh,m_Adchigh,waveStart,waveEnd,wave_calc_len,base_calc_len);
		//}

		//calc low gain cosmic ray adc and base
		base_calc_len=0;
		wave_calc_len=0;
		if(peakPosL>0)	{	waveStart = peakPosL-1;}
		else			{	waveStart = 0;}
		if(peakPosL<wave_len-2)	{	waveEnd = peakPosL+2;}
		else					{	waveEnd = wave_len-1;}
		for(int i=waveStart;i<=waveEnd;i++){
			m_Adclow += merged_pulselow.at(i);
			wave_calc_len++;
		}
		for(int i=0;i<waveStart-2;i++){
			m_Baselow += merged_pulselow.at(i);
			base_calc_len++;
		}
		for(int i=waveEnd+3;i<wave_len;i++){
			m_Baselow += merged_pulselow.at(i);
			base_calc_len++;
		}
		m_Baselow = m_Baselow/(4*base_calc_len);
		m_Adclow -= m_Baselow*4*wave_calc_len;
	}
	else
	{
		double preBaseH=0;
		double preBaseL=0;
		if(peakPosH<6)			{	for(int i=14;i<wave_len;i++)	{preBaseH += merged_pulsehigh.at(i);}}
		if(peakPosH>wave_len-9)	{	for(int i=0;i<wave_len-14;i++)	{preBaseH += merged_pulsehigh.at(i);}}
		else					{
									for(int i=0;i<peakPosH-5;i++)         { preBaseH += merged_pulsehigh.at(i);}
									for(int i=peakPosH+9;i<wave_len;i++)	{ preBaseH += merged_pulsehigh.at(i);}
								}
		if(peakPosL<6)			{	for(int i=14;i<wave_len;i++)	{preBaseL += merged_pulselow.at(i);}}
		if(peakPosL>wave_len-9)	{	for(int i=0;i<wave_len-14;i++)	{preBaseL += merged_pulselow.at(i);}}
		else					{
									for(int i=0;i<peakPosL-5;i++)         { preBaseL += merged_pulselow.at(i);}
									for(int i=peakPosL+9;i<wave_len;i++)	{ preBaseL += merged_pulselow.at(i);}
								}
		preBaseH /= (wave_len-14);
		preBaseL /= (wave_len-14);

		int ipulseH;
		//calc high gain laser adc and base
		base_calc_len=0;
		wave_calc_len=0;
		ipulseH = peakPosH;
		m_Adchigh += merged_pulsehigh.at(ipulseH);
		wave_calc_len++;
		ipulseH -= 1;
		while(1){
			if(ipulseH<0){break;}
			m_Adchigh += merged_pulsehigh.at(ipulseH);
			wave_calc_len++;
			ipulseH -= 1;
			if(ipulseH<=0 || merged_pulsehigh.at(ipulseH)<preBaseH+200){
				if(ipulseH<0){break;}
				m_Adchigh += merged_pulsehigh.at(ipulseH);
				wave_calc_len++;
				break;
			}
		}
		if(ipulseH>=0)	{	waveStart = ipulseH;}
		else			{	waveStart = 0;}
		for(int i=0;i<waveStart-2;i++){
			m_Basehigh += merged_pulsehigh.at(i);
			base_calc_len++;
		}
		ipulseH = peakPosH+1;
		if(ipulseH<wave_len){
			m_Adchigh += merged_pulsehigh.at(ipulseH);
			wave_calc_len++;
			ipulseH += 1;
		}
		while(1){
			if(ipulseH>wave_len-1){break;}
			m_Adchigh += merged_pulsehigh.at(ipulseH);
			wave_calc_len++;
			ipulseH += 1;
			if(ipulseH>=wave_len-1 || merged_pulsehigh.at(ipulseH)<preBaseH+200){
				if(ipulseH>wave_len-1){break;}
				m_Adchigh += merged_pulsehigh.at(ipulseH);
				wave_calc_len++;
				break;
			}
		}
		if(ipulseH<=wave_len-1) {   waveEnd = ipulseH;}
		else					{   waveEnd = wave_len-1;}
		for(int i=waveEnd+3;i<wave_len;i++){
			m_Basehigh += merged_pulsehigh.at(i);
			base_calc_len++;
		}
		m_Basehigh = m_Basehigh/(4*base_calc_len);
		m_Adchigh -= m_Basehigh*4*wave_calc_len;
		//if(isipm==1){
		//	printf("BaseH:%f AdcH:%f waveStart:%d waveEnd:%d wave_calc_len:%d base_calc_len:%d preBaseH:%lf\n",m_Basehigh,m_Adchigh,waveStart,waveEnd,wave_calc_len,base_calc_len,preBaseH);
		//}

		int ipulseL;
		//calc low gain laser adc and base
		base_calc_len=0;
		wave_calc_len=0;
		ipulseL = peakPosL;
		m_Adclow += merged_pulselow.at(ipulseL);
		wave_calc_len++;
		ipulseL -= 1;
		while(1){
			if(ipulseL<0){break;}
			m_Adclow += merged_pulselow.at(ipulseL);
			wave_calc_len++;
			ipulseL -= 1;
			if(ipulseL<=0 || merged_pulselow.at(ipulseL)<preBaseL+9){
				if(ipulseL<0){break;}
				m_Adclow += merged_pulselow.at(ipulseL);
				wave_calc_len++;
				break;
			}
		}
		if(ipulseL>=0)	{	waveStart = ipulseL;}
		else			{	waveStart = 0;}
		for(int i=0;i<waveStart-2;i++){
			m_Baselow += merged_pulselow.at(i);
			base_calc_len++;
		}
		ipulseL = peakPosL+1;
		if(ipulseL<wave_len){
			m_Adclow += merged_pulselow.at(ipulseL);
			wave_calc_len++;
			ipulseL += 1;
		}
		while(1){
			if(ipulseL>wave_len-1){break;}
			m_Adclow += merged_pulselow.at(ipulseL);
			wave_calc_len++;
			ipulseL += 1;
			if(ipulseL>=wave_len-1 || merged_pulselow.at(ipulseL)<preBaseL+9){
				if(ipulseL>wave_len-1){break;}
				m_Adclow += merged_pulselow.at(ipulseL);
				wave_calc_len++;
				break;
			}
		}
		if(ipulseL<=wave_len-1) {   waveEnd = ipulseL;}
		else					{   waveEnd = wave_len-1;}
		for(int i=waveEnd+3;i<wave_len;i++){
			m_Baselow += merged_pulselow.at(i);
			base_calc_len++;
		}
		m_Baselow = m_Baselow/(4*base_calc_len);
		m_Adclow -= m_Baselow*4*wave_calc_len;
	}
}

void WFCTAMerge::FindPeak(int isipm, vector<WFCTAMerge> &evs)
{
	WFCTAMerge::WaveForm_Merge(isipm,evs);
	double sumhighmax = -1000;
	double sumhigh;
	double sumlowmax = -1000;
	double sumlow;
	for(int ii=0;ii<merged_pulsehigh.size()-1;ii++){
		sumhigh = merged_pulsehigh.at(ii)+merged_pulsehigh.at(ii+1);
		sumlow = merged_pulselow.at(ii)+merged_pulselow.at(ii+1);
		if(sumhighmax<sumhigh){
			sumhighmax = sumhigh;
			if(merged_pulsehigh.at(ii)>merged_pulsehigh.at(ii+1)){
				peakPosH = ii; peakAmpH = merged_pulsehigh.at(ii);
			}
			else{
				peakPosH = ii+1; peakAmpH = merged_pulsehigh.at(ii+1);
			}
		}
		if(sumlowmax<sumlow){
			sumlowmax = sumlow;
			if(merged_pulselow.at(ii)>merged_pulselow.at(ii+1)){
				peakPosL = ii; peakAmpL = merged_pulselow.at(ii);
			}
			else{
				peakPosL = ii+1; peakAmpL = merged_pulselow.at(ii+1);
			}
		}
	}
}

void WFCTAMerge::WaveForm_Merge(int isipm, vector<WFCTAMerge> &evs)
{
	vector<WFCTAMerge>::iterator evs_iter;
	merged_pulsehigh.clear();
	merged_pulselow.clear();
	int first_event=0;
	emptyPulse=0;
	for(evs_iter=evs.begin(); evs_iter!=evs.end(); evs_iter++)
	{
		if(!first_event)
		{
			for(int ipoint=0;ipoint<28;ipoint++)
			{
				if((*evs_iter).pulsehigh[isipm][ipoint]==0){emptyPulse++;continue;}
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
