#ifndef WFCTADECODE_H
#define WFCTADECODE_H

#include <stdint.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include "TObject.h"

using namespace std;

class WFCTADecode//: public TObject
{
public:
    WFCTADecode();
    ~WFCTADecode();

    uint8_t StatusPackCheck(uint8_t *begin, int bufsize);
    bool StatusPackCheck(uint8_t *begin, int bufsize,int type);

    bool bigPackCheck(uint8_t *begin, int bufsize);
    void Find_SiPMs(uint8_t *begin, int packsize);

    uint64_t PackSize() {return packSize;};
    map<short, int>& GetSiPM_Position() {return m_sipm_position;};

    uint64_t eventId(uint8_t *begin);
    uint64_t RabbitTime(uint8_t *begin);
    //uint64_t Rabbittime(uint8_t *begin);
    double Rabbittime(uint8_t *begin);

    uint8_t GetPeak(uint8_t *begin, short isipm);
    bool Getgain_marker(uint8_t *begin, short isipm);
    uint8_t Getmypeak(uint8_t *begin, short isipm);
    uint16_t GetSingle_Thresh(uint8_t *begin, short isipm);
    uint16_t GetRecord_Thresh(uint8_t *begin, short isipm);
    bool GetOver_Single_Mark(uint8_t *begin, short isipm);
    bool GetOver_Record_Mark(uint8_t *begin, short isipm);
    float GetADC_Cut(uint8_t *begin, short isipm);
    float AdcHigh(uint8_t *begin, short isipm);
    float AdcLow(uint8_t *begin, short isipm);
    float BaseHigh(uint8_t *begin, short isipm);
    float BaseLow(uint8_t *begin, short isipm);
    float GetmyImageBaseHigh(uint8_t *begin, short isipm);
    float GetmyImageBaseLow(uint8_t *begin, short isipm);
    float GetmyImageAdcHigh(uint8_t *begin, short isipm);
    float GetmyImageAdcLow(uint8_t *begin, short isipm);
    void GetWaveForm(uint8_t *begin, short isipm, int *pulseh, int *pulsel);
    

    void Getthresh(uint8_t *begin, int packsize, short *single_thresh, short *record_thresh);//Deal21Package
    void Deal22Pack(uint8_t *begin, int packsize, long *single_count);//Deal22Package
    void Deal23Pack(uint8_t *begin, int packsize, long *single_count, long *single_time);//Deal23Package
    void GetHV(uint8_t *begin, int packsize, float *HV);//Deal81Package
    void GetPreTemp(uint8_t *begin, int packsize, float *PreTemp);//Deal82Package
    void GetClbTemp(uint8_t *begin, int packsize, float *ClbTemp);//Deal85pack

    uint64_t GetclbInitialTime(uint8_t *begin, int packsize);
    double GetclbInitialtime(uint8_t *begin, int packsize);
    int GetFiredTube(uint8_t *begin, int packsize);
    uint64_t GetStatusReadbackTime(uint8_t *begin, int packsize);
    //uint64_t GetStatusReadbacktime(uint8_t *begin, int packsize);
    double GetStatusReadbacktime(uint8_t *begin, int packsize);

private:

    void waveform(uint8_t *begin, short isipm);
    void Calc_Q_Base(uint8_t *begin, short isipm);

    //uint8_t *m_buffer;
    //size_t m_size;
    int64_t head,tail;
    int64_t readPos;

    uint8_t status_pack_mark;
    int64_t packSize;
    map<short, int> m_sipm_position;
    map<short, int>::iterator m_sipm_position_iter;

    //char m_peak;
    float m_adc_high;
    float m_adc_low;
    float m_base_high;
    float m_base_low;

    float m_Basehigh;
    float m_Baselow;
    float m_Adchigh;
    float m_Adclow;

    uint8_t m_mypeak;

    int pulsehigh[32];
    int pulselow[32];

  ClassDef(WFCTADecode,1);
};

#endif // WFCTADECODE_H
