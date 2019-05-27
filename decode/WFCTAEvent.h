#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include "TObject.h"

using namespace std;
class WFCTAEvent //: public TObject 
{
public:
	WFCTAEvent();
	~WFCTAEvent();
	bool OpenFile();
	bool OpenFile(const char* s);
	inline void CloseFile(){fInputStream.close();};
	int FirstBigpackage();
	bool Big_Package_Head();
        bool Little_Package();
        bool Big_Package_Tail();

	bool Status_Package(int *fpga_marker, int *package_marker);

        bool EndFile();

        void BranchInitial();
	void StatusInitial();
	void DealBigPackageHead();
        void ReadLittlePackage();
	void LittlePackageCalc(int *pulseh, int *pulsel);
        void DealBigPackageTail();

	void Deal21Package(short *single_thresh, short *record_thresh);
	void Deal22Package(long *single_count,float *DbTemp);
	void Deal23Package(long *single_count,long *single_time);
	void Deal81Package(float *HV);
	void Deal82Package(float *PreTemp);
	void Deal83Package(float *BigResistence);
	void Deal84Package(float *SmallResistence);
	void Deal85Package(long *ClbTime,float *ClbTemp);
	void DealFPGAPackage(int *fpgaVer);
	void DealF9Package(int *fpgaVer);

public:
	void SC_Channel2SiPM(short F_DB, short mChannel, short *mSiPM);

        inline vector<short>&  GetSiPM() {return WSiPM;};
	long GetEvent();
        inline vector<short>&  GetSC() {return WSC;};
        inline vector<short>&  GetChannel() {return WChannel;};
        long GetrabbitTime();
        double Getrabbittime();
        short Getpackagenum();
        inline vector<bool>&  Getgain_marker() {return Wgain_marker;};
        inline vector<char>&  Getpeak() {return Wpeak;};
        inline vector<char>&  Getmypeak() {return Wmypeak;};
        inline vector<short>&  GetSingle_Threshold() {return WSingle_Threshold;};
        inline vector<short>&  GetRecord_Threshold() {return WRecord_Threshold;};
        inline vector<bool>&  GetOver_Single_Marker() {return WOver_Single_Marker;};
        inline vector<bool>&  GetOver_Record_Marker() {return WOver_Record_Marker;};
        inline vector<float>&  GetADC_Cut() {return WADC_Cut;};
//	inline vector<double>&  GetImageX() {return WImageX;};
//	inline vector<double>&  GetImageY() {return WImageY;};
        inline vector<float>&  GetImageBaseHigh() {return WImageBaseHigh;};
        inline vector<float>&  GetImageBaseLow() {return WImageBaseLow;};
        inline vector<float>&  GetImageAdcHigh() {return WImageAdcHigh;};
        inline vector<float>&  GetImageAdcLow() {return WImageAdcLow;};
        inline vector<float>&  GetmyImageBaseHigh() {return WmyImageBaseHigh;};
        inline vector<float>&  GetmyImageBaseLow() {return WmyImageBaseLow;};
        inline vector<float>&  GetmyImageAdcHigh() {return WmyImageAdcHigh;};
        inline vector<float>&  GetmyImageAdcLow() {return WmyImageAdcLow;};

	long GetclbInitialTime() {return Wclb_initial_Time;};
        double GetclbInitialtime() {return Wclb_initial_time;};
	int GetFiredTube() {return Wfired_tube;};
	long GetStatusReadbackTime() {return Wstatus_readback_Time;};
	double GetStatusReadbacktime() {return Wstatus_readback_time;};

private:
	string fFileName;
	ifstream fInputStream;
	short temp;
	vector<long> bigpackagehead;
	vector<int> littlepackage;//[126];
        vector<long> bigpackagetail;//[24];
	vector<long> fpgapackage;//[64]
	vector<long> clb85package;//[74]
	vector<long> clb_db_package;//[72]
	bool jude;
private:
        vector<short> WSiPM;
        long WEvent;
        vector<short> WSC;
        vector<short> WChannel;
	long WrabbitTime;
	double Wrabbittime;
	int Wpackagenum;
	vector<bool> Wgain_marker;
	vector<char> Wpeak;
	vector<char> Wmypeak;
	vector<short> WSingle_Threshold;
	vector<short> WRecord_Threshold;
	vector<bool> WOver_Single_Marker;
	vector<bool> WOver_Record_Marker;
	vector<float> WADC_Cut;
//        vector<double> WImageX;
//        vector<double> WImageY;
	vector<float> WImageBaseHigh;
	vector<float> WImageBaseLow;
	vector<float> WImageAdcHigh;
	vector<float> WImageAdcLow;
	vector<float> WmyImageBaseHigh;
	vector<float> WmyImageBaseLow;
	vector<float> WmyImageAdcHigh;
	vector<float> WmyImageAdcLow;

	int s_SC;
	long Wclb_initial_Time;
	double Wclb_initial_time;
	int Wfired_tube;
	long Wstatus_readback_Time;
	double Wstatus_readback_time;
	float WDbTemp;
	long Wsingle_time;
	long WClbTime;
	float WClbTemp;

   ClassDef(WFCTAEvent,1);
};
