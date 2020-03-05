#ifndef WFCTAMERGE_H
#define WFCTAMERGE_H

#include <stdint.h>
#include <iostream>
#include <vector>
#include <string>

class WFCTAMerge
{
	protected:
		static short emptyPulse;
		static int32_t peakAmpH;
		static int32_t peakAmpL;
		static int16_t peakPosH;
		static int16_t peakPosL;
		static float m_Basehigh;
		static float m_Baselow;
		static float m_Adchigh;
		static float m_Adclow;
		static std::vector<int> merged_pulsehigh;
		static std::vector<int> merged_pulselow;

		static void FindPeak(int isipm, std::vector<WFCTAMerge> &evs);
		static void WaveForm_Merge(int isipm, std::vector<WFCTAMerge> &evs);
	public:
		int big_pack_lenth;
		long eEvent;
		long rabbitTime;
		double rabbittime;
		short n_fired;
		short n_Channel;
		short packCheck;
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
		static long GeteEvent(std::vector<WFCTAMerge> &evs);
		static long RabbitTime(std::vector<WFCTAMerge> &evs);
		static long Rabbittime(std::vector<WFCTAMerge> &evs);
		//static void packCheck_Merge(std::vector<WFCTAMerge> &evs);
		static int GetBigPackLen(std::vector<WFCTAMerge> &evs);
		static short GetNFired(std::vector<WFCTAMerge> &evs);
		static short GetNChannel(std::vector<WFCTAMerge> &evs);
		static long eevent_Merge(int isipm, std::vector<WFCTAMerge> &evs);
		static short zipmod_Merge(int isipm, std::vector<WFCTAMerge> &evs);
		static bool OvSigMarker_Merge(int isipm, std::vector<WFCTAMerge> &evs);
		static bool OvRecMarker_Merge(int isipm, std::vector<WFCTAMerge> &evs);
		static float WimSum_Merge(int isipm, std::vector<WFCTAMerge> &evs);
		static short GetPeakPosH(int isipm, std::vector<WFCTAMerge> &evs);
		static short GetPeakPosL(int isipm, std::vector<WFCTAMerge> &evs);
		static int GetPeakAmpH(int isipm, std::vector<WFCTAMerge> &evs);
		static int GetPeakAmpL(int isipm, std::vector<WFCTAMerge> &evs);
		static void Calc_Q_Base(int isipm, std::vector<WFCTAMerge> &evs, int laserCalc);
		static float GetBaseH(int isipm, std::vector<WFCTAMerge> &evs);
		static float GetBaseL(int isipm, std::vector<WFCTAMerge> &evs);
		static float GetAdcH(int isipm, std::vector<WFCTAMerge> &evs);
		static float GetAdcL(int isipm, std::vector<WFCTAMerge> &evs);
		static float GetLaserBaseH(int isipm, std::vector<WFCTAMerge> &evs);
		static float GetLaserBaseL(int isipm, std::vector<WFCTAMerge> &evs);
		static float GetLaserAdcH(int isipm, std::vector<WFCTAMerge> &evs);
		static float GetLaserAdcL(int isipm, std::vector<WFCTAMerge> &evs);
		static int eSatH_Merge(int isipm, std::vector<WFCTAMerge> &evs);
		static int eSatL_Merge(int isipm, std::vector<WFCTAMerge> &evs);

		static bool IsData_Merge(int isipm, std::vector<WFCTAMerge> &evs);

};

#endif // WFCTAMERGE_H
