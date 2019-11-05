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
bool WFCTADecode::StatusPack(uint8_t *begin, int bufsize, int64_t packStart)
{
	//head = 0;
	//tail = 0;

	readPos = 0+packStart;
	while(readPos<bufsize)
	{   
		if(   *(begin+readPos+0)==0xbb && *(begin+readPos+1)==0xbb
				&& *(begin+readPos+2)==0xbb && *(begin+readPos+3)==0xbb
				&& *(begin+readPos+4)==0xee && *(begin+readPos+5)==0xee
				&& *(begin+readPos+6)==0xee && *(begin+readPos+7)==0xee
				&& *(begin+readPos+41544)==0xaa && *(begin+readPos+41545)==0xaa
				&& *(begin+readPos+41546)==0xaa && *(begin+readPos+41547)==0xaa
				&& *(begin+readPos+41548)==0xdd && *(begin+readPos+41549)==0xdd
				&& *(begin+readPos+41550)==0xdd && *(begin+readPos+41551)==0xdd)
		{   
			//printf("a status pack\n");
			return true;
		}
		readPos++;
	}
	return false;
}

/***************************
 * **find each status pack**
 * *************************/
uint8_t WFCTADecode::StatusPackCheck(uint8_t *begin, int bufsize, int64_t packStart)
{
	readPos = 0+packStart;
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
	packSize = readPos;
	return 100;
}

int WFCTADecode::statusPackCheck(uint8_t *begin, int bufsize,int type)
{
	int find=-1;
	readPos = 0;
	if(type>=1&&type<=9){
		while(readPos+63<bufsize){
			if( *(begin+readPos+0)==0x12 && *(begin+readPos+1)==0x34 && *(begin+readPos+62)==0xab && *(begin+readPos+63)==0xcd && *(begin+readPos+2)==(type==9?0:type) ){
				packSize = readPos+64;
				find=readPos;
				break;
			}
			readPos++;
		}
	}
	else if((type>=0x21&&type<=0x23)||(type>=0x81&&type<=0x84)){
		while(readPos+71<bufsize){
			if( *(begin+readPos+0)==0x12 && *(begin+readPos+1)==0x34 && *(begin+readPos+70)==0xab && *(begin+readPos+71)==0xcd && *(begin+readPos+3)==type ){
				packSize = readPos+72;
				find=readPos;
				break;
			}
			readPos++;
		}
	}
	else if(type==0x85){
		while(readPos+73<bufsize){
			if( *(begin+readPos+0)==0x12 && *(begin+readPos+1)==0x34 && *(begin+readPos+72)==0xab && *(begin+readPos+73)==0xcd && *(begin+readPos+3)==type ){
				packSize = readPos+74;
				find=readPos;
				break;
			}
			readPos++;
		}
	}
	return find;
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
void WFCTADecode::Deal22Pack(uint8_t *begin, int packsize, long *single_count)//Deal22Package get single count [channel 1-8]
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
void WFCTADecode::Deal23Pack(uint8_t *begin, int packsize, long *single_count, long *single_time)//Deal23Package get single [channel 9-16] count and single time
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

void WFCTADecode::GetBigRes(uint8_t *begin, int packsize, float *BigResistence)//Deal83Package
{
	head = packsize - 72;

	short sipm;
	short fpga = (short)begin[head+2]&0x0f;
	short db = ((short)begin[head+2])>>4;
	short sc = db*10+fpga;

	for(int i=0;i<4;i++)
	{
		SC_Channel2SiPM(sc,i+1,&sipm);
		*(BigResistence+sipm) = (float)( ((int16_t)begin[head+13+i*2]<<8) | (int16_t)begin[head+14+i*2] ) *50/256.;
	}
	SC_Channel2SiPM(sc,5,&sipm);
	*(BigResistence+sipm) = (float)( ((int16_t)begin[head+13+4*2]<<8) | (int16_t)begin[head+14+5*2] ) *50/256.;
	for(int i=5;i<16;i++)
	{
		SC_Channel2SiPM(sc,i+1,&sipm);
		*(BigResistence+sipm) = (float)( ((int16_t)begin[head+13+(i+1)*2]<<8) | (int16_t)begin[head+14+(i+1)*2] ) *50/256.;
	}

}

void WFCTADecode::GetSmallRes(uint8_t *begin, int packsize, float *SmallResistence)//Deal84Package
{
	head = packsize - 72;

	short sipm;
	short fpga = (short)begin[head+2]&0x0f;
	short db = ((short)begin[head+2])>>4;
	short sc = db*10+fpga;

	SC_Channel2SiPM(sc,1,&sipm);
	*(SmallResistence+sipm) = (float)( ((int16_t)begin[head+13]<<8) | (int16_t)begin[head+14+2] )  *10/1024.;
	for(int i=1;i<16;i++)
	{
		SC_Channel2SiPM(sc,i+1,&sipm);
		*(SmallResistence+sipm) = (float)( ((int16_t)begin[head+13+(i+1)*2]<<8) | (int16_t)begin[head+14+(i+1)*2] ) *10/1024.;
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

void WFCTADecode::GetClbTime(uint8_t *begin, int packsize, long *ClbTime)//Deal85pack
{
	head = packsize - 74;

	short sipm;
	short fpga = (short)begin[head+2]&0x0f;
	short db = ((short)begin[head+2])>>4;
	short sc = db*10+fpga;

	long m_ClbTime = (long)( ((uint64_t)begin[head+10]<<32)|
			((uint64_t)begin[head+11]<<24)|
			((uint64_t)begin[head+12]<<16)|
			((uint64_t)begin[head+13]<<8)|
			(uint64_t)begin[head+14] );
	for(int i=0;i<16;i++)
	{
		SC_Channel2SiPM(sc,i+1,&sipm);
		*(ClbTime+sipm) = m_ClbTime;
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

/*********************************
 * **get mode info in F1-F8 pack**
 * *******************************/
int WFCTADecode::GetF18Version(uint8_t *begin, int packsize)
{
	int m_f18version = (int)( (int32_t)begin[packsize-3]);
	printf("m_f18version:%0x\n",m_f18version);
	return m_f18version;
}

int WFCTADecode::GetDBVersion(uint8_t *begin, int packsize)
{
	int m_dbversion = (int)( (int32_t)begin[packsize-3]);
	printf("m_dbversion:%0x\n",m_dbversion);
	return m_dbversion;
}
int WFCTADecode::GetDBNumber(uint8_t *begin, int packsize)
{
	int db = (int)( (int8_t)(begin[packsize-70]>>4)&0xf);
	int fpga = (int)( (int8_t)(begin[packsize-70])&0xf);
	int m_dbnumber = db*10 + fpga;
	printf("m_dbnumber:%d\n",m_dbnumber);
	return m_dbnumber;
}

int WFCTADecode::GetClbVersion(uint8_t *begin, int packsize)
{
	int m_clbversion = (int)( ((int32_t)begin[packsize-6]<<8) | (int32_t)begin[packsize-5]);
	printf("m_clbversion:%0x\n",m_clbversion);
	return m_clbversion;
}
int WFCTADecode::GetClbNumber(uint8_t *begin, int packsize)
{
	int db = (int)( (int8_t)(begin[packsize-72]>>4)&0xf);
	int fpga = (int)( (int8_t)(begin[packsize-72])&0xf);
	int m_clbnumber = db*10 + fpga;
	printf("m_clbnumber:%d\n",m_clbnumber);
	return m_clbnumber;
}

/******************************
 * **get mode info in F9 pack**
 * ****************************/
int WFCTADecode::GetF9Mode(uint8_t *begin, int packsize)
{
	int m_f9_mode = (int)( (int32_t)begin[packsize-60]);
	printf("m_f9_mode:%0x\n",m_f9_mode);
	return m_f9_mode;
}
int WFCTADecode::GetPattern(uint8_t *begin, int packsize)
{
	int m_pattern = (int)( (int32_t)begin[packsize-59]);
	printf("m_pattern:%0x\n",m_pattern);
	return m_pattern;
}
int WFCTADecode::GetF9Version(uint8_t *begin, int packsize)
{
	int m_f9version = (int)( (int32_t)begin[packsize-3]);
	printf("m_f9version:%0x\n",m_f9version);
	return m_f9version;
}
int WFCTADecode::GetF9PlusVersion(uint8_t *begin, int packsize)
{
	int m_f9plusversion = (int)( (int32_t)begin[packsize-16]);
	printf("m_f9pversion:%0x\n",m_f9plusversion);
	return m_f9plusversion;
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


/****************************
 * **find FEE Data Fragment**
 * **************************/
bool WFCTADecode::FEEDataFragment(uint8_t *begin)
{
	//dumpPacket(begin,4);
	FEEDataHead = 0;
	while(FEEDataHead<20)
	{
		if( *(begin+FEEDataHead)==0xbb && *(begin+FEEDataHead+1)==0x34 && *(begin+FEEDataHead+2)==0x12 && *(begin+FEEDataHead+3)==0xbb )
			return true;
		else
			FEEDataHead++;
	}
	return false;
}

/**********************
 * **find big package**
 * ********************/
bool WFCTADecode::bigPackCheck(uint8_t *begin, int bufsize, int64_t packStart)
{
	//int64_t big_pack_len;
	head = 0;
	tail = 0;

	readPos = 0+packStart;
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

	readPos = 0+packStart;
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

	big_pack_len = tail + 1 - head;
	packSize = tail+1;
	if(big_pack_len>1)
	{
		//dumpPacket(begin+head,24);printf("%lld ** %lld ** packsize:%lld | \n",head,tail,packSize);
		return true;
	}
	else
	{
		return false;
	}
}

/***********************
 * **  get slice size  **
 * *********************/
int32_t WFCTADecode::sliceLength(uint8_t *begin, int feedatahead)
{
	int32_t sliceLength = ((int32_t)begin[feedatahead+7]<<24)|
		((int32_t)begin[feedatahead+6]<<16)|
		((int32_t)begin[feedatahead+5]<<8)|
		((int32_t)begin[feedatahead+4]);
	//printf("totle size:%d\n",sliceLength);
	return sliceLength;
}

/**********************
 * **   get tel id   **
 * ********************/
short WFCTADecode::Telid(uint8_t *begin, int feedatahead)
{
	short telId = (short)begin[feedatahead+8];
	//printf("telid:%d\n",telId);
	return telId;
}

/**********************
 * **   find sipms   **
 * ********************/
void WFCTADecode::Find_SiPMs(uint8_t *begin)//, int packStart)
{
	short fpga,db;
	short sc,channel,sipm;
	m_sipm_position.clear();
	readPos = head;

	while(readPos<tail)
	{
		if(   *(begin+readPos+0)==0xaa && *(begin+readPos+1)==0xaa
				&& *(begin+readPos+124)==0xbb && *(begin+readPos+125)==0xbb)
		{
			fpga = *(begin+readPos+5)&0x0f;
			db = (*(begin+readPos+5)>>4)&0x0f;
			sc = db*10+fpga;
			channel = *(begin+readPos+4);
			SC_Channel2SiPM(sc,channel,&sipm);
			m_sipm_position.insert(pair<short,int>(sipm,(int)readPos));
			readPos += 126;
		}
		else
		{
			readPos++;
		}
	}
	//dumpPacket(begin,packSize,16);
}

/**********************
 * **  get event ID  **
 * ********************/
uint64_t WFCTADecode::eventId(uint8_t *begin)
{
	uint64_t evtid = ((uint64_t)begin[head+12]<<24)|
		((uint64_t)begin[head+13]<<16)|
		((uint64_t)begin[head+14]<<8)|
		((uint64_t)begin[head+15]);
	//dumpPacket(begin,24,16);
	//printf("event:%llu\n",evtid);
	return evtid;
}

/************************
 * ** get Npix in data **
 * **********************/
int16_t WFCTADecode::nFired(uint8_t *begin)
{
	int16_t nfired = ((int16_t)begin[head+22]<<8) | ((int16_t)begin[head+23]);
	return nfired;
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
	//printf("%02x %02x %02x %02x %02x %02x %02x %02x %02x %02x\n",begin[tail-15],begin[tail-14],begin[tail-13],begin[tail-12],begin[tail-11],begin[tail-10],begin[tail-9],begin[tail-8],begin[tail-7],begin[tail-6]);
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

/*-------------------------------------------------------------------------------*/
//little package message//
/*-------------------------------------------------------------------------------*/

/********************************************
 * ** get event id in each channel package **
 * ******************************************/
uint64_t WFCTADecode::eventId_in_channel(uint8_t *begin, short isipm)
{   
	m_sipm_position_iter = m_sipm_position.find(isipm);
	int packposition = m_sipm_position_iter->second;

	uint64_t evtid_in_channel = ((uint64_t)begin[packposition+2]<<2)|
		((uint64_t)begin[packposition+3]>>6);
	return evtid_in_channel;
}

/****************************************************************
 * ** get zip mode, 00 represent 4 point zipped into one point **
 * **************************************************************/
int16_t WFCTADecode::zipMode(uint8_t *begin, short isipm)
{
	m_sipm_position_iter = m_sipm_position.find(isipm);
	int packposition = m_sipm_position_iter->second;

	int16_t zip_mode = ((int16_t)((begin[packposition+3]>>4)&0x03));
	return zip_mode;
}

/**********************************************
 * ** get trigger marker, over_single_marker **
 * ********************************************/
bool WFCTADecode::GetOver_Single_Mark(uint8_t *begin, short isipm)
{
	m_sipm_position_iter = m_sipm_position.find(isipm);
	int packposition = m_sipm_position_iter->second;

	bool m_Over_Single_Mark = (begin[packposition+3]&0x1); 
	return m_Over_Single_Mark;
}

/*********************************************
 * ** get trigger marker, over record marker**
 * *******************************************/
bool WFCTADecode::GetOver_Record_Mark(uint8_t *begin, short isipm)
{
	m_sipm_position_iter = m_sipm_position.find(isipm);
	int packposition = m_sipm_position_iter->second;

	bool m_Over_Record_Mark = 0;
	int m_over_record_mark = begin[m_sipm_position_iter->second+3]&0x0e;
	//printf("m_over_record_mark:%d\n",m_over_record_mark);
	if(m_over_record_mark==0x0e)
	{
		m_Over_Record_Mark=1;
	}
	//bool m_Over_Record_Mark = begin[m_sipm_position_iter->second+3]&0xe;
	return m_Over_Record_Mark;
}

/**********************
 * ** high gain base **
 * ********************/
float WFCTADecode::Getwinsum(uint8_t *begin, short isipm)
{
	m_sipm_position_iter = m_sipm_position.find(isipm);
	int packposition = m_sipm_position_iter->second;

	float m_winsum = (float)( ((uint64_t)begin[packposition+64]<<24)|
			((uint64_t)begin[packposition+65]<<16)|
			((uint64_t)begin[packposition+66]<<8)|
			((uint64_t)begin[packposition+67]) ) ;
	return m_winsum;
}


/***********************************************************************
 * ** get PeakPosH,PeakPosL,PeakAmH,PeakAmL, which calc from waveform **
 * *********************************************************************/
uint8_t WFCTADecode::GetPeakPosH(uint8_t *begin, short isipm)
{
	WFCTADecode::wavepeak(begin,isipm);
	return m_wavepeakH;
}

uint8_t WFCTADecode::GetPeakPosL(uint8_t *begin, short isipm)
{
	WFCTADecode::wavepeak(begin,isipm);
	return m_wavepeakL;
}

int32_t WFCTADecode::GetPeakAmH(uint8_t *begin, short isipm)
{
	WFCTADecode::wavepeak(begin,isipm);
	return peakAmpH;
}

int32_t WFCTADecode::GetPeakAmL(uint8_t *begin, short isipm)
{
	WFCTADecode::wavepeak(begin,isipm);
	return peakAmpL;
}

/*****************************************************************
 * ** get qhigh/qlow/basehigh/baselow, which calc from waveform **
 * ***************************************************************/
float WFCTADecode::GetwaveImageBaseHigh(uint8_t *begin, short isipm)
{
	WFCTADecode::Calc_Q_Base(begin,isipm);
	return m_Basehigh;
}

float WFCTADecode::GetwaveImageBaseLow(uint8_t *begin, short isipm)
{   
	WFCTADecode::Calc_Q_Base(begin,isipm);
	return m_Baselow;
}

float WFCTADecode::GetwaveImageAdcHigh(uint8_t *begin, short isipm)
{
	WFCTADecode::Calc_Q_Base(begin,isipm);
	return m_Adchigh;
}

float WFCTADecode::GetwaveImageAdcLow(uint8_t *begin, short isipm)
{
	WFCTADecode::Calc_Q_Base(begin,isipm);
	return m_Adclow;
}

bool WFCTADecode::eSaturationHigh(uint8_t *begin, short isipm)
{
	WFCTADecode::Stauration(begin,isipm);
	return eSatH;
}

bool WFCTADecode::eSaturationLow(uint8_t *begin, short isipm)
{
	WFCTADecode::Stauration(begin,isipm);
	return eSatL;
}

/******************************
 * ** get wave form [public] **
 * ****************************/
void WFCTADecode::GetWaveForm(uint8_t *begin, short isipm, int *pulseh, int *pulsel)
{
	WFCTADecode::waveform(begin,isipm);
	for(int i=0;i<28;i++)
	{
		*(pulseh+isipm*28+i) = pulsehigh[i];
		*(pulsel+isipm*28+i) = pulselow[i];
	}
	//printf("%5d ",pulsehigh[i]);
	//printf("\n");
	//for(int i=0;i<28;i++)
	//{
	//printf("%5d ",pulselow[i]);
	//}
	//printf("\n");
	//dumpPacket(begin+m_sipm_position_iter->second+6,112,16);
	//dumpPacket(begin+m_sipm_position_iter->second+60,5);
	//dumpPacket(begin+m_sipm_position_iter->second+65,39,3);
}


void WFCTADecode::GeteSaturation(uint8_t *begin, short isipm, int *esath, int *esatl)
{
	WFCTADecode::waveform(begin,isipm);
	for(int i=0;i<28;i++)
	{
		*(esath+isipm*28+i) = saturationH[i];
		*(esatl+isipm*28+i) = saturationL[i];
	}
}


/************************************************
 * ** get saturation marker of high & low gain **
 * **********************************************/
void WFCTADecode::Stauration(uint8_t *begin, short isipm)
{
	WFCTADecode::waveform(begin,isipm);
	eSatH = 0;
	eSatL = 0;
	for(int i=0;i<28;i++)
	{
		if(saturationH[i]==1) {eSatH = 1;}
		if(saturationL[i]==1) {eSatL = 1;}
	}
}

/*************************************
 * ** calc q and base from waveform **
 * ***********************************/
void WFCTADecode::Calc_Q_Base(uint8_t *begin, short isipm)
{
	WFCTADecode::wavepeak(begin,isipm);
	m_Basehigh = 0;
	m_Baselow = 0;
	m_Adchigh = 0;
	m_Adclow = 0;
	if(m_wavepeakH<3)
	{
		for(int i=6;i<28;i++)  { m_Basehigh += pulsehigh[i];}
		for(int i=0;i<6;i++)   { m_Adchigh += pulsehigh[i]; }
	}
	else if(m_wavepeakH>23)
	{
		for(int i=0;i<22;i++)  { m_Basehigh += pulsehigh[i];}
		for(int i=22;i<28;i++) { m_Adchigh += pulsehigh[i]; }
	}
	else
	{
		for(int i=0;i<m_wavepeakH-2;i++)             { m_Basehigh += pulsehigh[i];}
		for(int i=m_wavepeakH+4;i<28;i++)            { m_Basehigh += pulsehigh[i];}
		for(int i=m_wavepeakH-2;i<m_wavepeakH+4;i++) { m_Adchigh += pulsehigh[i]; }
	}

	if(m_wavepeakL<3)
	{
		for(int i=6;i<28;i++)  { m_Baselow += pulselow[i];}
		for(int i=0;i<6;i++)   { m_Adclow += pulselow[i]; }
	}
	else if(m_wavepeakL>23)
	{
		for(int i=0;i<22;i++)  { m_Baselow += pulselow[i];}
		for(int i=22;i<28;i++) { m_Adclow += pulselow[i]; }
	}
	else
	{
		for(int i=0;i<m_wavepeakL-2;i++)             { m_Baselow += pulselow[i];}
		for(int i=m_wavepeakL+4;i<28;i++)            { m_Baselow += pulselow[i];}
		for(int i=m_wavepeakL-2;i<m_wavepeakL+4;i++) { m_Adclow += pulselow[i]; }
	}

	m_Basehigh = m_Basehigh/88.;
	m_Baselow = m_Baselow/88.;
	m_Adchigh -= m_Basehigh*24;
	m_Adclow -= m_Baselow*24;

}

/******************************
 * ** get wave form [privat] **
 * ****************************/
void WFCTADecode::waveform(uint8_t *begin, short isipm)
{   
	m_sipm_position_iter = m_sipm_position.find(isipm);
	int packposition = m_sipm_position_iter->second;
	int waveStart1 = packposition+6;
	int waveStart2 = packposition+12;

	for(int i=0; i<14; i++)
	{   
		pulsehigh[i] = ((int)(begin[waveStart1+i*4]&0x7f)<<8)|((int)begin[waveStart1+i*4+1]);
		saturationH[i] = ((int)((begin[waveStart1+i*4]>>7)&0x01));
		pulselow[i]  = ((int)(begin[waveStart1+i*4+2]&0x7f)<<8)|((int)begin[waveStart1+i*4+3]);
		saturationL[i] = ((int)((begin[waveStart1+i*4+2]>>7)&0x01));
	}   

	for(int i=14; i<28; i++)
	{   
		pulsehigh[i] = ((int)(begin[waveStart2+i*4]&0x7f)<<8)|((int)begin[waveStart2+i*4+1]);
		saturationH[i] = ((int)((begin[waveStart1+i*4]>>7)&0x01));
		pulselow[i]  = ((int)(begin[waveStart2+i*4+2]&0x7f)<<8)|((int)begin[waveStart2+i*4+3]);
		saturationL[i] = ((int)((begin[waveStart1+i*4+2]>>7)&0x01));
	}
}


/***************************************
 * ** get peak and peakam in waveform **
 * *************************************/
void WFCTADecode::wavepeak(uint8_t *begin, short isipm)
{
	WFCTADecode::waveform(begin,isipm);
	double sumhighmax = -1000;
	double sumhigh;
	for(int i=0;i<28;i++){
		sumhigh = pulsehigh[i];
		if(sumhighmax<sumhigh) {sumhighmax = sumhigh; m_wavepeakH = i; peakAmpH = sumhigh;}
	}
	double sumlowmax = -1000;
	double sumlow;
	for(int i=0;i<28;i++){
		sumlow = pulselow[i];
		if(sumlowmax<sumlow) {sumlowmax = sumlow; m_wavepeakL = i; peakAmpL = sumlow;}
	}
}





