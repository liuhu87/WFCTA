#ifndef WFCTAMERGE_H
#define WFCTAMERGE_H

#include <iostream>
#include <vector>
#include <string>
#include <map>

using namespace std;
class WFCTAMerge
{
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
		static void WaveForm_Merge(int isipm, vector<WFCTAMerge> &evs);

};

#endif // WFCTAMERGE_H
