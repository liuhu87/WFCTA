#include <stdlib.h>
#include <iostream>
#include "WFCTADecode.h"
#include "dumpPack.h"
#include "camera.h"

using namespace std;

ClassImp(WFCTADecode);

WFCTADecode::WFCTADecode()
{
}

WFCTADecode::~WFCTADecode()
{
}


/*******************
 * *****************
 * ***status data***
 * *****************
 * *****************/

/**********************
 * **find status pack**
 * ********************/
uint8_t WFCTADecode::StatusPackCheck(uint8_t *begin, int bufsize)
{
    readPos = 0;
    while(readPos<bufsize)
    {
        if( *(begin+readPos+0)==0x12 && *(begin+readPos+1)==0x34 && *(begin+readPos+62)==0xab && *(begin+readPos+63)==0xcd ){
            packSize = readPos+64;
            status_pack_mark = *(begin+readPos+2);//FPGA 1-9 PACK
	    if(status_pack_mark==0){status_pack_mark = 9;}
            return status_pack_mark;
        }

        if( *(begin+readPos+0)==0x12 && *(begin+readPos+1)==0x34 && *(begin+readPos+70)==0xab && *(begin+readPos+71)==0xcd){
            packSize = readPos+72;
	    status_pack_mark = *(begin+readPos+3);
            return status_pack_mark;
        }

        if( *(begin+readPos+0)==0x12 && *(begin+readPos+1)==0x34 && *(begin+readPos+72)==0xab && *(begin+readPos+73)==0xcd){
            packSize = readPos+74;
            status_pack_mark = *(begin+readPos+3);
            return status_pack_mark;
        }
        readPos++;
    }
    return 0;
}

/*****************************************
 * **get single_thresh and record_thresh**
 * ***************************************/
void WFCTADecode::Getthresh(uint8_t *begin, int packsize, short *single_thresh, short *record_thresh)//Deal21Package
{
    head = packsize - 72;

    short sipm;
    short fpga = (short)begin[head+2]&0x0f;
    short db = ((short)begin[head+2])>>4;
    short sc = db*10+fpga;

    for(int i=0;i<15;i++)
    {
        SC_Channel2SiPM(sc,i+1,&sipm);
        *(single_thresh+sipm) = ((int16_t)begin[head+4+4*i]<<8) | (int16_t)begin[head+5+4*i];
    }
        SC_Channel2SiPM(sc,16,&sipm);
        *(single_thresh+sipm) = ((int16_t)begin[head+66]<<8) | (int16_t)begin[head+67];

    for(int i=0;i<14;i++)
    {
        SC_Channel2SiPM(sc,i+1,&sipm);
        *(record_thresh+sipm) = ((int16_t)begin[head+6+4*i]<<8) | (int16_t)begin[head+7+4*i];
    }
        SC_Channel2SiPM(sc,15,&sipm);
        *(record_thresh+sipm) = ((int16_t)begin[head+64]<<8) | (int16_t)begin[head+65];
        SC_Channel2SiPM(sc,16,&sipm);
        *(record_thresh+sipm) = ((int16_t)begin[head+68]<<8) | (int16_t)begin[head+69];
}

/*******************************
 * **deal status package of 22**
 * *****************************/
void WFCTADecode::Deal22Pack(uint8_t *begin, int packsize, long *single_count)//Deal22Package
{
    head = packsize - 72;
    short sipm;
    short fpga = (short)begin[head+2]&0x0f;
    short db = ((short)begin[head+2])>>4;
    short sc = db*10+fpga;

    for(int i=0;i<8;i++)
    {       
        SC_Channel2SiPM(sc,i+1,&sipm);
        *(single_count+sipm) =  ((int64_t)begin[head+9+i*5]<<32) | 
				((int64_t)begin[head+10+i*5]<<24) | 
				((int64_t)begin[head+11+i*5]<<16) | 
				((int64_t)begin[head+12+i*5]<<8) | 
				((int64_t)begin[head+13+i*5]);
    }        

        //WDbTemp = (clb_db_package[4]<<8)+clb_db_package[5];
        //*(DbTemp+DbTempCount) = s_SC;
        //*(DbTemp+64+DbTempCount) = WDbTemp;
        //DbTempCount++;
}

/*******************************
 * **deal status package of 23**
 * *****************************/
void WFCTADecode::Deal23Pack(uint8_t *begin, int packsize, long *single_count, long *single_time)//Deal23Package
{
    head = packsize - 72;
    short sipm;
    short fpga = (short)begin[head+2]&0x0f;
    short db = ((short)begin[head+2])>>4;
    short sc = db*10+fpga;

    long m_single_time = ((int64_t)begin[head+44]<<32) |
			 ((int64_t)begin[head+45]<<24) | 
			 ((int64_t)begin[head+48]<<16) |
			 ((int64_t)begin[head+49]<<8) |
			 ((int64_t)begin[head+50]);
    for(int i=0;i<16;i++)
    {
        SC_Channel2SiPM(sc,i+1,&sipm);
        *(single_time+sipm) = m_single_time;
    }
    for(int i=8;i<16;i++)
    {
        SC_Channel2SiPM(sc,i+1,&sipm);
        *(single_count+sipm) =  ((int64_t)begin[head+4+(i-8)*5]<<32) |
				((int64_t)begin[head+5+(i-8)*5]<<24) |
				((int64_t)begin[head+6+(i-8)*5]<<16) |
				((int64_t)begin[head+7+(i-8)*5]<<8) |
				((int64_t)begin[head+8+(i-8)*5]);
    }
}

/*******************************
 * **deal status package of 81**
 * *****************************/
void WFCTADecode::GetHV(uint8_t *begin, int packsize, float *HV)//Deal81Package
{
    head = packsize - 72;
    //dumpPacket(begin+head,72,16);

    short sipm;
    short fpga = (short)begin[head+2]&0x0f;
    short db = ((short)begin[head+2])>>4;
    short sc = db*10+fpga;

    //printf("fpga:%d, db:%d, sipm:%d\n",fpga,db,sipm);
 
    SC_Channel2SiPM(sc,1,&sipm);
    //SC_Channel2eSiPM(fpga, db, 1, &sipm); sipm -= 1;
    *(HV+sipm) = (float)( ((uint32_t)begin[head+13]<<19) | ((uint32_t)begin[head+14]<<11) | ((uint32_t)begin[head+15]<<3) | ((uint32_t)begin[head+16]>>5) );

    SC_Channel2SiPM(sc,2,&sipm);
    //SC_Channel2eSiPM(fpga, db, 2, &sipm); sipm -= 1;
    *(HV+sipm) = (float)( (((uint32_t)begin[head+16]&0x1f)<<22) | ((uint32_t)begin[head+17]<<14) | ((uint32_t)begin[head+18]<<6) | ((uint32_t)begin[head+19]>>2) );

    SC_Channel2SiPM(sc,3,&sipm);
    //SC_Channel2eSiPM(fpga, db, 3, &sipm); sipm -= 1;
    *(HV+sipm) = (float)( (((uint32_t)begin[head+19]&0x02)<<25) | ((uint32_t)begin[head+20]<<17) | ((uint32_t)begin[head+21]<<9) | ((uint32_t)begin[head+22]<<1) | ((uint32_t)begin[head+23]>>7) );

    SC_Channel2SiPM(sc,4,&sipm);
    //SC_Channel2eSiPM(fpga, db, 4, &sipm); sipm -= 1;
    *(HV+sipm) = (float)( (((uint32_t)begin[head+23]&0x7f)<<20) | ((uint32_t)begin[head+24]<<12) | ((uint32_t)begin[head+25]<<4) | ((uint32_t)begin[head+26]>>4) );

    SC_Channel2SiPM(sc,5,&sipm);
    //SC_Channel2eSiPM(fpga, db, 5, &sipm); sipm -= 1;
    *(HV+sipm) = (float)( (((uint32_t)begin[head+26]&0x0f)<<23) | ((uint32_t)begin[head+27]<<15) | ((uint32_t)begin[head+28]<<7) | ((uint32_t)begin[head+29]>>1) );

    SC_Channel2SiPM(sc,6,&sipm);
    //SC_Channel2eSiPM(fpga, db, 6, &sipm); sipm -= 1;
    *(HV+sipm) = (float)( (((uint32_t)begin[head+29]&0x01)<<26) | ((uint32_t)begin[head+30]<<18) | ((uint32_t)begin[head+31]<<10) | ((uint32_t)begin[head+32]<<2) | ((uint32_t)begin[head+33]>>6) );

    SC_Channel2SiPM(sc,7,&sipm);
    //SC_Channel2eSiPM(fpga, db, 7, &sipm); sipm -= 1;
    *(HV+sipm) = (float)( (((uint32_t)begin[head+33]&0x3f)<<21) | ((uint32_t)begin[head+34]<<13) | ((uint32_t)begin[head+35]<<5) | ((uint32_t)begin[head+36]>>3) );

    SC_Channel2SiPM(sc,8,&sipm);
    //SC_Channel2eSiPM(fpga, db, 8, &sipm); sipm -= 1;
    *(HV+sipm) = (float)( (((uint32_t)begin[head+36]&0x07)<<24) | ((uint32_t)begin[head+37]<<16) | ((uint32_t)begin[head+40]<<8) | ((uint32_t)begin[head+41]) );

    SC_Channel2SiPM(sc,9,&sipm);
    //SC_Channel2eSiPM(fpga, db, 9, &sipm); sipm -= 1;
    *(HV+sipm) = (float)( ((uint32_t)begin[head+42]<<19) | ((uint32_t)begin[head+43]<<11) | ((uint32_t)begin[head+44]<<3) | ((uint32_t)begin[head+45]>>5) );

    SC_Channel2SiPM(sc,10,&sipm);
    //SC_Channel2eSiPM(fpga, db, 10, &sipm); sipm -= 1;
    *(HV+sipm) = (float)( (((uint32_t)begin[head+45]&0x1f)<<22) | ((uint32_t)begin[head+46]<<14) | ((uint32_t)begin[head+47]<<6) | ((uint32_t)begin[head+48]>>2) );

    SC_Channel2SiPM(sc,11,&sipm);
    //SC_Channel2eSiPM(fpga, db, 11, &sipm); sipm -= 1;
    *(HV+sipm) = (float)( (((uint32_t)begin[head+48]&0x02)<<25) | ((uint32_t)begin[head+49]<<17) | ((uint32_t)begin[head+50]<<9) | ((uint32_t)begin[head+51]<<1) | ((uint32_t)begin[head+52]>>7) );

    SC_Channel2SiPM(sc,12,&sipm);
    //SC_Channel2eSiPM(fpga, db, 12, &sipm); sipm -= 1;
    *(HV+sipm) = (float)( (((uint32_t)begin[head+52]&0x7f)<<20) | ((uint32_t)begin[head+53]<<12) | ((uint32_t)begin[head+54]<<4) | ((uint32_t)begin[head+55]>>4) );

    SC_Channel2SiPM(sc,13,&sipm);
    //SC_Channel2eSiPM(fpga, db, 13, &sipm); sipm -= 1;
    *(HV+sipm) = (float)( (((uint32_t)begin[head+55]&0x0f)<<23) | ((uint32_t)begin[head+56]<<15) | ((uint32_t)begin[head+57]<<7) | ((uint32_t)begin[head+58]>>1) );

    SC_Channel2SiPM(sc,14,&sipm);
    //SC_Channel2eSiPM(fpga, db, 14, &sipm); sipm -= 1;
    *(HV+sipm) = (float)( (((uint32_t)begin[head+58]&0x01)<<26) | ((uint32_t)begin[head+59]<<18) | ((uint32_t)begin[head+60]<<10) | ((uint32_t)begin[head+61]<<2) | ((uint32_t)begin[head+62]>>6) );

    SC_Channel2SiPM(sc,15,&sipm);
    //SC_Channel2eSiPM(fpga, db, 15, &sipm); sipm -= 1;
    *(HV+sipm) = (float)( (((uint32_t)begin[head+62]&0x3f)<<21) | ((uint32_t)begin[head+63]<<13) | ((uint32_t)begin[head+64]<<5) | ((uint32_t)begin[head+65]>>3) );

    SC_Channel2SiPM(sc,16,&sipm);
    //SC_Channel2eSiPM(fpga, db, 16, &sipm); sipm -= 1;
    *(HV+sipm) = (float)( (((uint32_t)begin[head+65]&0x07)<<24) | ((uint32_t)begin[head+66]<<16) | ((uint32_t)begin[head+67]<<8) | ((uint32_t)begin[head+68]) );

    for(int i=0;i<16;i++)
    {
	SC_Channel2SiPM(sc,i+1,&sipm);
        //SC_Channel2eSiPM(fpga, db, i+1, &sipm); sipm -= 1;
        *(HV+sipm) /= (512*427.4087);
    }
}

/*******************************
 * **deal status package of 82**
 * *****************************/
void WFCTADecode::GetPreTemp(uint8_t *begin, int packsize, float *PreTemp)//Deal82Package
{
    head = packsize - 72;
    //dumpPacket(begin+head,72,16);

    double A = 0.00433;
    double B = 13.582;
    double C;

    short sipm;
    short fpga = (short)begin[head+2]&0x0f;
    short db = ((short)begin[head+2])>>4;
    short sc = db*10+fpga;

    for(int i=0;i<8;i++)
    {
	SC_Channel2SiPM(sc,i+1,&sipm);
        //SC_Channel2eSiPM(fpga, db, i+1, &sipm); sipm -= 1;
        *(PreTemp+sipm) = (float)( ((int16_t)begin[head+13+i*2]<<8) | ((int16_t)begin[head+14+i*2]) );
    }
	SC_Channel2SiPM(sc,9,&sipm);
        //SC_Channel2eSiPM(fpga, db, 9, &sipm); sipm -= 1;
        *(PreTemp+sipm) = (float)( ((int16_t)begin[head+13+8*2]<<8) | ((int16_t)begin[head+14+9*2]) );
    for(int i=9;i<16;i++)
    {
	SC_Channel2SiPM(sc,i+1,&sipm);
        //SC_Channel2eSiPM(fpga, db, i+1, &sipm); sipm -= 1;
        *(PreTemp+sipm) = (float)( ((int16_t)begin[head+13+(i+1)*2]<<8) | ((int16_t)begin[head+14+(i+1)*2]) );
    }
    for(int i=0;i<16;i++)
    {
        SC_Channel2SiPM(sc,i+1,&sipm);
        //SC_Channel2eSiPM(fpga, db, i+1, &sipm); sipm -= 1;
        C = *(PreTemp+sipm)*10000./32768-2230.8;
        *(PreTemp+sipm) = (-1*B+sqrt(B*B-4*A*C))/(2*A);
        *(PreTemp+sipm) += 30;
    }

}

/********************************
 * ***get clb board temperature**
 * ******************************/
void WFCTADecode::GetClbTemp(uint8_t *begin, int packsize, float *ClbTemp)//Deal85pack
{
    head = packsize - 74;

    short sipm;
    short fpga = (short)begin[head+2]&0x0f;
    short db = ((short)begin[head+2])>>4;
    short sc = db*10+fpga;

    float m_ClbTemp = (float)( ((int16_t)begin[head+15]<<8) | (int16_t)begin[head+16] );
    if(m_ClbTemp<40960) {m_ClbTemp /= 256.;}
    else                {m_ClbTemp = (65536-m_ClbTemp)/256.;}

    for(int i=0;i<16;i++)
    {
	SC_Channel2SiPM(sc,i+1,&sipm);
        //SC_Channel2eSiPM(fpga, db, i+1, &sipm); sipm -= 1;
        *(ClbTemp+sipm) = m_ClbTemp;
    }
}

/**********************
 * **get initial time**
 * ********************/
uint64_t WFCTADecode::GetclbInitialTime(uint8_t *begin, int packsize)
{
    uint64_t m_clb_initial_Time = ((uint64_t)begin[packsize-26]<<30)|
                                  ((uint64_t)begin[packsize-25]<<22)|
                                  ((uint64_t)begin[packsize-24]<<14)|
                                  ((uint64_t)begin[packsize-23]<<6)|
                                  ((uint64_t)begin[packsize-22]>>2&0x3f);
    return m_clb_initial_Time;
}
double WFCTADecode::GetclbInitialtime(uint8_t *begin, int packsize)
{
    double m_clb_initial_time = (double)( (((uint64_t)begin[packsize-22]&0x03)<<24)|
                                          ((uint64_t)begin[packsize-21]<<16)|
                                          ((uint64_t)begin[packsize-20]<<8)|
                                          ((uint64_t)begin[packsize-19]<<0) );
    return m_clb_initial_time;
}

/**********************************
 * **get setted fired tube number**
 * ********************************/
int WFCTADecode::GetFiredTube(uint8_t *begin, int packsize)
{
    int m_fired_tube = (int)( ((int32_t)begin[packsize-18]<<8) | (int32_t)begin[packsize-17]);
    return m_fired_tube;
}

/******************************
 * **get status readback time**
 * ****************************/
uint64_t  WFCTADecode::GetStatusReadbackTime(uint8_t *begin, int packsize)
{
    uint64_t m_status_readback_Time = ((uint64_t)begin[packsize-11]<<30)|
                        	      ((uint64_t)begin[packsize-10]<<22)|
                        	      ((uint64_t)begin[packsize-9]<<14)|
                        	      ((uint64_t)begin[packsize-8]<<6)|
	                              ((uint64_t)begin[packsize-7]>>2&0x3f);
    return m_status_readback_Time;
}
double  WFCTADecode::GetStatusReadbacktime(uint8_t *begin, int packsize)
{
    double m_status_readback_time = (double)( (((uint64_t)begin[packsize-7]&0x03)<<24)|
    		                              ((uint64_t)begin[packsize-6]<<16)|
                  		              ((uint64_t)begin[packsize-5]<<8)|
                        	              ((uint64_t)begin[packsize-4]<<0) );
    return m_status_readback_time;
}



/*******************
 * *****************
 * ***events data***
 * *****************
 * *****************/

/**********************
 * **find big package**
 * ********************/
bool WFCTADecode::bigPackCheck(uint8_t *begin, int bufsize)
{
    int64_t big_pack_len;
    head = 0;
    tail = 0;

    readPos = 0;
    while(readPos<bufsize)
    {
	if(   *(begin+readPos+0)==0xcc && *(begin+readPos+1)==0xcc 
	   && *(begin+readPos+2)==0xdd && *(begin+readPos+3)==0xdd
	   && *(begin+readPos+4)==0xee && *(begin+readPos+5)==0xee
	   && *(begin+readPos+6)==0xff && *(begin+readPos+7)==0xff)
	{
	    head = readPos;
	    break;
	}
        readPos++;
    }

    readPos = 0;
    while(readPos<bufsize)
    {
	if(   *(begin+readPos+0)==0x11 && *(begin+readPos+1)==0x11 
	   && *(begin+readPos+2)==0x22 && *(begin+readPos+3)==0x22  
	   && *(begin+readPos+4)==0x33 && *(begin+readPos+5)==0x33 
	   && *(begin+readPos+6)==0x44 && *(begin+readPos+7)==0x44)
	{
	    tail = readPos+7;
	    break;
	}
	readPos++;
    }

    big_pack_len = tail - head;
    packSize = tail+1;
    if(big_pack_len>0)
    {
	    //dumpPacket(begin+head,24);printf("%lld ** %lld ** packsize:%lld | \n",head,tail,packSize);
	    return true;
    }
    else
    {
	return false;
    }
}

/**********************
 * **   find sipms   **
 * ********************/
void WFCTADecode::Find_SiPMs(uint8_t *begin, int packsize)
{
    short fpga,db;
    short sc,channel,sipm;
    m_sipm_position.clear();
    readPos = 0;
    while(readPos<packsize)
    {
        if(   *(begin+readPos+0)==0xaa && *(begin+readPos+1)==0xaa
           && *(begin+readPos+124)==0xbb && *(begin+readPos+125)==0xbb)
        {
	    fpga = *(begin+readPos+5)&0x0f;
	    db = *(begin+readPos+5)>>4;
	    sc = db*10+fpga;
	    channel = *(begin+readPos+4);
	    SC_Channel2SiPM(sc,channel,&sipm);
	    m_sipm_position.insert(pair<short,int>(sipm,(int)readPos));
        }
        readPos++;
    }
    //dumpPacket(begin,packsize,16);
}

/**********************
 * **  get event ID  **
 * ********************/
uint64_t WFCTADecode::eventId(uint8_t *begin)
{
    uint64_t evtid = ((uint64_t)begin[head+8]<<24)|
                      ((uint64_t)begin[head+9]<<16)|
                      ((uint64_t)begin[head+10]<<8)|
                      ((uint64_t)begin[head+11]);
    //dumpPacket(begin,24,16);
    //printf("event:%llu\n",evtid);
    return evtid;
}

/**********************
 * ** get rabbitTime **
 * ********************/
uint64_t WFCTADecode::RabbitTime(uint8_t *begin)
{
    uint64_t rab_Time = ((uint64_t)begin[tail-15]<<30)|
                        ((uint64_t)begin[tail-14]<<22)|
                        ((uint64_t)begin[tail-13]<<14)|
                        ((uint64_t)begin[tail-12]<<6)|
                        ((uint64_t)begin[tail-11]>>2&0x3f);
    return rab_Time;
}
/**********************
 * ** get rabbittime **
 * ********************/
double WFCTADecode::Rabbittime(uint8_t *begin)
{
    double rab_time = (double)( (((uint64_t)begin[tail-11]&0x03)<<24)|
                                ((uint64_t)begin[tail-10]<<16)|
                                ((uint64_t)begin[tail-9]<<8)|
                                ((uint64_t)begin[tail-8]<<0) );
    return rab_time;
}

/****************************
 * **get wave peak position**
 * **************************/
uint8_t WFCTADecode::GetPeak(uint8_t *begin, short isipm)
{
    m_sipm_position_iter = m_sipm_position.find(isipm);
    uint8_t m_peak = ((uint8_t)begin[m_sipm_position_iter->second+2]*64)|
	     ((uint8_t)begin[m_sipm_position_iter->second+3]>>2);
    //dumpPacket(begin+m_sipm_position_iter->second,128,16);
    return m_peak;
}

/**********************
 * **single threshold**
 * ********************/
uint16_t WFCTADecode::GetSingle_Thresh(uint8_t *begin, short isipm)
{
    m_sipm_position_iter = m_sipm_position.find(isipm);
    uint16_t m_Single_Threshold = ((uint16_t)begin[m_sipm_position_iter->second+120]<<8)|
			       ((uint16_t)begin[m_sipm_position_iter->second+121]);
    return m_Single_Threshold;
}

/**********************
 * **record threshold**
 * ********************/
uint16_t WFCTADecode::GetRecord_Thresh(uint8_t *begin, short isipm)
{
    m_sipm_position_iter = m_sipm_position.find(isipm);
    uint16_t m_Record_Threshold = ((uint16_t)begin[m_sipm_position_iter->second+122]<<8)|
                               ((uint16_t)begin[m_sipm_position_iter->second+123]);
    return m_Record_Threshold;
}

/************************
 * **over single marker**
 * **********************/
bool WFCTADecode::GetOver_Single_Mark(uint8_t *begin, short isipm)
{
    m_sipm_position_iter = m_sipm_position.find(isipm);
    bool m_Over_Single_Mark = begin[m_sipm_position_iter->second+3]&0x2;
    return m_Over_Single_Mark;
}

/************************
 * **over record marker**
 * **********************/
bool WFCTADecode::GetOver_Record_Mark(uint8_t *begin, short isipm)
{
    m_sipm_position_iter = m_sipm_position.find(isipm);
    bool m_Over_Record_Mark = begin[m_sipm_position_iter->second+3]&0x1;
    return m_Over_Record_Mark;
}

/***********************
 * **   high gain Q   **
 * *********************/
float WFCTADecode::AdcHigh(uint8_t *begin, short isipm)
{
    m_sipm_position_iter = m_sipm_position.find(isipm);
    m_adc_high = (float)( ((uint64_t)begin[m_sipm_position_iter->second+104]<<24)|
		          ((uint64_t)begin[m_sipm_position_iter->second+105]<<16)|
		          ((uint64_t)begin[m_sipm_position_iter->second+106]<<8)|
		          ((uint64_t)begin[m_sipm_position_iter->second+107]) ) - m_base_high*64.;
    return m_adc_high;
}
/**********************
 * **   low gian Q   **
 * ********************/
float WFCTADecode::AdcLow(uint8_t *begin, short isipm)
{
    m_sipm_position_iter = m_sipm_position.find(isipm);
    m_adc_low = (float)( ((uint64_t)begin[m_sipm_position_iter->second+108]<<24)|
                         ((uint64_t)begin[m_sipm_position_iter->second+109]<<16)|
                         ((uint64_t)begin[m_sipm_position_iter->second+110]<<8)|
                         ((uint64_t)begin[m_sipm_position_iter->second+111]) ) - m_base_low*64.;
    return m_adc_low;
}
/**********************
 * ** high gain base **
 * ********************/
float WFCTADecode::BaseHigh(uint8_t *begin, short isipm)
{
    m_sipm_position_iter = m_sipm_position.find(isipm);
    m_base_high = (float)( ((uint64_t)begin[m_sipm_position_iter->second+112]<<24)|
                          ((uint64_t)begin[m_sipm_position_iter->second+113]<<16)|
                          ((uint64_t)begin[m_sipm_position_iter->second+114]<<8)|
                          ((uint64_t)begin[m_sipm_position_iter->second+115]) ) / 256;
    return m_base_high;
}
/***********************
 * **  low gain base  **
 * *********************/
float WFCTADecode::BaseLow(uint8_t *begin, short isipm)
{
    m_sipm_position_iter = m_sipm_position.find(isipm);
    m_base_low = (float)( ((uint64_t)begin[m_sipm_position_iter->second+116]<<24)|
                         ((uint64_t)begin[m_sipm_position_iter->second+117]<<16)|
                         ((uint64_t)begin[m_sipm_position_iter->second+118]<<8)|
                         ((uint64_t)begin[m_sipm_position_iter->second+119]) ) / 256;
    return m_base_low;
}

/***********************
 * ** H|L gain marker **
 * *********************/
bool WFCTADecode::Getgain_marker(uint8_t *begin, short isipm)
{
    WFCTADecode::waveform(begin,isipm);
    bool m_gain_marker = 0;
    for(int i=0;i<29;i++){
	if(pulsehigh[i]>4000){m_gain_marker = 1;}
    }
    return m_gain_marker;
}


/****************************
 * ** get peak in waveform **
 * **************************/
uint8_t WFCTADecode::Getmypeak(uint8_t *begin, short isipm)
{
    WFCTADecode::waveform(begin,isipm);
    double sumhighmax = -1000;
    double sumhigh;
    for(int i=0;i<29;i++){
	sumhigh = pulsehigh[i]+pulsehigh[i+1]+pulsehigh[i+2]+pulsehigh[i+3];
        sumhigh /=4.;
	if(sumhighmax<sumhigh) {sumhighmax = sumhigh; m_mypeak = i+1;}
    }
    return m_mypeak;
}

/*********************************************
 * ** get maximum calculation of four point **
 * *******************************************/
float WFCTADecode::GetADC_Cut(uint8_t *begin, short isipm)
{
    WFCTADecode::waveform(begin,isipm);
    double Four_Point_Q = -1000;
    double sum_4;
    for(int i=0;i<29;i++){
        sum_4 = pulsehigh[i]+pulsehigh[i+1]+pulsehigh[i+2]+pulsehigh[i+3];
        if(Four_Point_Q<sum_4) {Four_Point_Q = sum_4;}
    }
    Four_Point_Q -= (4. * m_base_high);
    return Four_Point_Q;

}

/*****************************************************************
 * ** get qhigh/qlow/basehigh/baselow, which calc from waveform **
 * ***************************************************************/
float WFCTADecode::GetmyImageBaseHigh(uint8_t *begin, short isipm)
{
    WFCTADecode::Calc_Q_Base(begin,isipm);
    return m_Basehigh;
}

float WFCTADecode::GetmyImageBaseLow(uint8_t *begin, short isipm)
{   
    WFCTADecode::Calc_Q_Base(begin,isipm);
    return m_Baselow;
}

float WFCTADecode::GetmyImageAdcHigh(uint8_t *begin, short isipm)
{
    WFCTADecode::Calc_Q_Base(begin,isipm);
    return m_Adchigh;
}

float WFCTADecode::GetmyImageAdcLow(uint8_t *begin, short isipm)
{
    WFCTADecode::Calc_Q_Base(begin,isipm);
    return m_Adclow;
}

/******************************
 * ** get wave form [public] **
 * ****************************/
void WFCTADecode::GetWaveForm(uint8_t *begin, short isipm, int *pulseh, int *pulsel)
{
//dumpPacket(begin+m_sipm_position_iter->second+6,54,3);
//dumpPacket(begin+m_sipm_position_iter->second+60,5);
//dumpPacket(begin+m_sipm_position_iter->second+65,39,3);
}




/*************************************
 * ** calc q and base from waveform **
 * ***********************************/
void WFCTADecode::Calc_Q_Base(uint8_t *begin, short isipm)
{
    WFCTADecode::Getmypeak(begin,isipm);
    m_Basehigh = 0;
    m_Baselow = 0;
    m_Adchigh = 0;
    m_Adclow = 0;
    if(m_mypeak<8)
    {
        for(int i=31;i>21;i--)  { m_Basehigh += pulsehigh[i]; m_Baselow += pulselow[i];}
        m_Basehigh = m_Basehigh/10.;
        m_Baselow = m_Baselow/10.;
	for(int i=0;i<m_mypeak+9;i++) { m_Adchigh += pulsehigh[i]-m_Basehigh; m_Adclow += pulselow[i]-m_Baselow;}
    }
    else if(m_mypeak>22)
    {
        for(int i=0;i<10;i++)  { m_Basehigh += pulsehigh[i]; m_Baselow += pulselow[i];}
        m_Basehigh = m_Basehigh/10.;
        m_Baselow = m_Baselow/10.;
        for(int i=m_mypeak-6;i<32;i++) { m_Adchigh += pulsehigh[i]-m_Basehigh; m_Adclow += pulselow[i]-m_Baselow;}
    }
    else
    {
        for(int i=0;i<5;i++)  { m_Basehigh += pulsehigh[i]; m_Baselow += pulselow[i];}
        for(int i=31;i>26;i--)  { m_Basehigh += pulsehigh[i]; m_Baselow += pulselow[i];}
        m_Basehigh = m_Basehigh/10.;
        m_Baselow = m_Baselow/10.;
        for(int i=m_mypeak-6;i<m_mypeak+9;i++) { m_Adchigh += pulsehigh[i]-m_Basehigh; m_Adclow += pulselow[i]-m_Baselow;}
    }
}

/******************************
 * ** get wave form [privat] **
 * ****************************/
void WFCTADecode::waveform(uint8_t *begin, short isipm)
{   
    m_sipm_position_iter = m_sipm_position.find(isipm);
    int waveStart1 = 6;
    int waveEnd1 = 58;
    int waveStart2 = 65;
    int waveEnd2 = 104;
    int middleWave = 60;
    int ipoint = 0;
    
    for(int i=waveStart1; i<waveEnd1; i = i+3)
    {   
        pulsehigh[ipoint] = ((int)begin[m_sipm_position_iter->second+i]<<4)|((int)begin[m_sipm_position_iter->second+i+1]>>4);
        pulselow[ipoint] = (((int)begin[m_sipm_position_iter->second+i+1]&0x0f)<<8)|((int)begin[m_sipm_position_iter->second+i+2]);
        ipoint++;
    }   
        pulsehigh[ipoint] = ((int)begin[m_sipm_position_iter->second+middleWave]<<4)|((int)begin[m_sipm_position_iter->second+middleWave+1]>>4);
        pulselow[ipoint] = (((int)begin[m_sipm_position_iter->second+middleWave+1]&0x0f)<<8)|((int)begin[m_sipm_position_iter->second+middleWave+4]);
        ipoint++;
    for(int i=waveStart2; i<waveEnd2; i = i+3)
    {   
        pulsehigh[ipoint] = ((int)begin[m_sipm_position_iter->second+i]<<4)|((int)begin[m_sipm_position_iter->second+i+1]>>4);
        pulselow[ipoint] = (((int)begin[m_sipm_position_iter->second+i+1]&0x0f)<<8)|((int)begin[m_sipm_position_iter->second+i+2]);
        ipoint++;
    }
}







