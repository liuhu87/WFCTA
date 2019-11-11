#ifndef WFCTAMERGE_H
#define WFCTAMERGE_H

#include <stdint.h>
#include <iostream>
#include <vector>
#include <string>
#include <map>

using namespace std;
class WFCTAMerge
{
	protected:
		static int32_t peakAmpH;
		static int32_t peakAmpL;
		static uint8_t peakPosH;
		static uint8_t peakPosL;
		static float m_Basehigh;
		static float m_Baselow;
		static float m_Adchigh;
		static float m_Adclow;
		static vector<int> merged_pulsehigh;
		static vector<int> merged_pulselow;

		static void FindPeak(int isipm, vector<WFCTAMerge> &evs);
		static void WaveForm_Merge(int isipm, vector<WFCTAMerge> &evs);
	public:
		int big_pack_lenth;
		long eEvent;
		long rabbitTime;
		double rabbittime;
		short n_fired;
		short n_Channel;
		long eevent[1024];
		short zipmod[1024];
		float winsum[1024];
		bool Over_Single_Marker[1024];
		bool Over_Record_Marker[1024];

		int Npoint[28];
		int IsData[1024];
		int saturationH[1024][28];
		int saturationL[1024][28];
		int pulsehigh[1024][28];
		int pulselow[1024][28];

	public:
		WFCTAMerge();
		~WFCTAMerge();
		void EventInitial();
		static long GeteEvent(vector<WFCTAMerge> &evs);
		static long RabbitTime(vector<WFCTAMerge> &evs);
		static long Rabbittime(vector<WFCTAMerge> &evs);
		static int GetBigPackLen(vector<WFCTAMerge> &evs);
		static short GetNFired(vector<WFCTAMerge> &evs);
		static short GetNChannel(vector<WFCTAMerge> &evs);
		static long eevent_Merge(int isipm, vector<WFCTAMerge> &evs);
		static short zipmod_Merge(int isipm, vector<WFCTAMerge> &evs);
		static bool OvSigMarker_Merge(int isipm, vector<WFCTAMerge> &evs);
		static bool OvRecMarker_Merge(int isipm, vector<WFCTAMerge> &evs);
		static float WimSum_Merge(int isipm, vector<WFCTAMerge> &evs);
		static char GetPeakPosH(int isipm, vector<WFCTAMerge> &evs);
		static char GetPeakPosL(int isipm, vector<WFCTAMerge> &evs);
		static int GetPeakAmpH(int isipm, vector<WFCTAMerge> &evs);
		static int GetPeakAmpL(int isipm, vector<WFCTAMerge> &evs);
		static void Calc_Q_Base(int isipm, vector<WFCTAMerge> &evs, int laserCalc);
		static float GetBaseH(int isipm, vector<WFCTAMerge> &evs);
		static float GetBaseL(int isipm, vector<WFCTAMerge> &evs);
		static float GetAdcH(int isipm, vector<WFCTAMerge> &evs);
		static float GetAdcL(int isipm, vector<WFCTAMerge> &evs);
		static float GetLaserBaseH(int isipm, vector<WFCTAMerge> &evs);
		static float GetLaserBaseL(int isipm, vector<WFCTAMerge> &evs);
		static float GetLaserAdcH(int isipm, vector<WFCTAMerge> &evs);
		static float GetLaserAdcL(int isipm, vector<WFCTAMerge> &evs);
		static int eSatH_Merge(int isipm, vector<WFCTAMerge> &evs);
		static int eSatL_Merge(int isipm, vector<WFCTAMerge> &evs);

		static bool IsData_Merge(int isipm, vector<WFCTAMerge> &evs);

};

#endif // WFCTAMERGE_H
