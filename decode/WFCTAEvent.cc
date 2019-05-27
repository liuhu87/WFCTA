#include <iostream>
#include "WFCTAEvent.h"
#include "WFCTAparameter.h"
using namespace std;

ClassImp(WFCTAEvent);

WFCTAEvent::WFCTAEvent()
{
	bigpackagehead.resize(16,0);
	littlepackage.resize(126,0);
	bigpackagetail.resize(24,0);
	fpgapackage.resize(64,0);
	clb85package.resize(74,0);
	clb_db_package.resize(72,0);
}

WFCTAEvent::~WFCTAEvent()
{

}

bool WFCTAEvent::OpenFile(const char* s)
{
	fFileName = s;
	return OpenFile();
}

bool WFCTAEvent::OpenFile()
{
	fInputStream.open(fFileName.c_str(),ios::in);
	if(!fInputStream) 
	{
		cerr<<"Error opening file "<<fFileName<<endl;
		return false;
	}
	else
	{
		return true;
	}
}

int WFCTAEvent::FirstBigpackage()
{
	while(true)
	{
		fInputStream.read((char*)&temp,1);
		if(temp == 204)
		{
			bigpackagehead[0] = 204;
			for(int i=1;i<16;i++)
			{
				fInputStream.read((char*)&temp,1);
				bigpackagehead[i] = temp;
			}
			if(bigpackagehead[0]==204&&bigpackagehead[1]==204&&bigpackagehead[2]==221&&bigpackagehead[3]==221&&bigpackagehead[4]==238&&bigpackagehead[5]==238&&bigpackagehead[6]==255&&bigpackagehead[7]==255)
			{
				fInputStream.seekg(-16,ios::cur);
				return fInputStream.tellg();
			}
			else
			{
				fInputStream.seekg(-15,ios::cur);
			}
		}
		if(!fInputStream)
		{
	                cerr<<"Error opening file "<<fFileName<<endl;
        	        return fInputStream.tellg();
		}
	}
}

void WFCTAEvent::BranchInitial()
{
	WSiPM.clear();
	WSC.clear();
	WChannel.clear();
	Wgain_marker.clear();
	Wpeak.clear();
	Wmypeak.clear();
	WSingle_Threshold.clear();
	WRecord_Threshold.clear();
	WOver_Single_Marker.clear();
	WOver_Record_Marker.clear();
	WADC_Cut.clear();
//	ImageX.clear();
//	ImageY.clear();
	WImageAdcHigh.clear();
	WImageAdcLow.clear();
	WImageBaseHigh.clear();
	WImageBaseLow.clear();
	WmyImageAdcHigh.clear();
	WmyImageAdcLow.clear();
	WmyImageBaseHigh.clear();
	WmyImageBaseLow.clear();
}
void WFCTAEvent::StatusInitial()
{
//map <int, double>::iterator clbTemp_Iter;
//for(clbTemp_Iter=clbTemp.begin(); clbTemp_Iter!=clbTemp.end(); clbTemp_Iter++){
//   printf("%d %lf\n",clbTemp_Iter->first,clbTemp_Iter->second);
//}
	DbTempCount = 0;
	ClbTimeCount = 0;
	ClbTempCount = 0;
}

bool WFCTAEvent::Big_Package_Head()
{
	if(temp == 204)
	{
		fInputStream.read((char*)&temp,1);
	if(temp == 204)
	{
		fInputStream.read((char*)&temp,1);
	if(temp == 221)
	{
		bigpackagehead[0] = 204;
		bigpackagehead[1] = 204;
		for(int i=2;i<16;i++)
		{
			bigpackagehead[i] = temp;
			fInputStream.read((char*)&temp,1);
		}
		if(bigpackagehead[0]==204&&bigpackagehead[1]==204&&bigpackagehead[2]==221&&bigpackagehead[3]==221&&bigpackagehead[4]==238&&bigpackagehead[5]==238&&bigpackagehead[6]==255&&bigpackagehead[7]==255)
		{
			return true;
		}
	}
	}
	}
	return false;
}

bool WFCTAEvent::Little_Package()
{
	jude = false;
	if(temp == 170)
	{
		fInputStream.read((char*)&temp,1);
	if(temp == 170)
	{
		littlepackage[0] = 170;
		for(int i=1;i<126;i++)
		{
			littlepackage[i] = temp;
			fInputStream.read((char*)&temp,1);
			if(fInputStream.eof()) {jude = true; break;}
		}
		if(jude)  {jude = false;}
		if(littlepackage[0]==170&&littlepackage[1]==170&&littlepackage[124]==187&&littlepackage[125]==187)
		{
			return true;
		}
		else
		{
//			Littlepackagenumhead++;
			fInputStream.seekg(-124,ios::cur);
		}
	}
	}
	else
	{
	        fInputStream.read((char*)&temp,1);
	}
	return false;
}

bool WFCTAEvent::Big_Package_Tail()
{
	if(temp == 17)
	{
		fInputStream.read((char*)&temp,1);
	if(temp == 17)
	{
		fInputStream.read((char*)&temp,1);
	if(temp == 34)
	{
		fInputStream.seekg(-19,ios::cur);
		for(int i=0;i<24;i++)
		{
			fInputStream.read((char*)&temp,1);
			bigpackagetail[i] = temp;
		}
		if(bigpackagetail[16]==17&&bigpackagetail[17]==17&&bigpackagetail[18]==34&&bigpackagetail[19]==34&&bigpackagetail[20]==51&&bigpackagetail[21]==51&&bigpackagetail[22]==68&&bigpackagetail[23]==68)
		{
			return true;
		}
	}
	}
	}
	return false;
}

bool WFCTAEvent::Status_Package(int *fpga_marker, int *package_marker)
{
	jude = false;
	if(temp == 18)
	{
		fInputStream.read((char*)&temp,1);
	if(temp == 52)
	{
		fInputStream.read((char*)&temp,1);
		*fpga_marker = temp;
		fInputStream.read((char*)&temp,1);
		*package_marker = temp;
//cout<<hex<<"fpga:"<<*fpga_marker<<" package:"<<*package_marker<<endl;
		if(*fpga_marker <10)
		{
			fpgapackage[0] = 18;
			fpgapackage[1] = 52;
			fpgapackage[2] = *fpga_marker;
			for(int i=3;i<64;i++)
			{
				fpgapackage[i] = temp;
				fInputStream.read((char*)&temp,1);
				if(fInputStream.eof()) {jude = true; break;}
			}
			if(jude)  {jude = false;}
			if(fpgapackage[0]==18&&fpgapackage[1]==52&&fpgapackage[62]==171&&fpgapackage[63]==205)
			{
				fInputStream.seekg(-1,ios::cur);
				return true;
			}
			else
			{
				fInputStream.seekg(-62,ios::cur);
			}
		}
		else if(*package_marker==0x85)
		{
			clb85package[0] = 18;
			clb85package[1] = 52;
			clb85package[2] = *fpga_marker;
			for(int i=3;i<74;i++)
			{
				clb85package[i] = temp;
				fInputStream.read((char*)&temp,1);
				if(fInputStream.eof()) {jude = true; break;}
			}
			if(jude)  {jude = false;}
			if(clb85package[0]==0x12&&clb85package[1]==0x34&&clb85package[72]==0xab&&clb85package[73]==0xcd)
			{
				fInputStream.seekg(-1,ios::cur);
				return true;
			}
			else
			{
				fInputStream.seekg(-72,ios::cur);
			}
		}
		else
		{
			clb_db_package[0] = 18;
			clb_db_package[1] = 52;
			clb_db_package[2] = *fpga_marker;
			for(int i=3;i<72;i++)
			{
				clb_db_package[i] = temp;
				fInputStream.read((char*)&temp,1);
				if(fInputStream.eof()) {jude = true; break;}
			}
			if(jude)  {jude = false;}
			if(clb_db_package[0]==0x12&&clb_db_package[1]==0x34&&clb_db_package[70]==0xab&&clb_db_package[71]==0xcd)
			{
				fInputStream.seekg(-1,ios::cur);
				return true;
			}
			else
			{
				fInputStream.seekg(-70,ios::cur);
			}
		}
	}
	}
	return false;
}
void WFCTAEvent::Deal21Package(short *single_thresh, short *record_thresh)
{
	s_SC = (clb_db_package[2]>>4)*10+clb_db_package[2]-(clb_db_package[2]>>4)*16;

	for(int i=0;i<15;i++)
	{
		SC_Channel2SiPM(s_SC,i+1,&SiPM);
		*(single_thresh+SiPM) = (clb_db_package[4+4*i]<<8)+clb_db_package[5+4*i];
	}
	SC_Channel2SiPM(s_SC,16,&SiPM);
	*(single_thresh+SiPM) = (clb_db_package[66]<<8)+clb_db_package[67];

	for(int i=0;i<14;i++)
	{
		SC_Channel2SiPM(s_SC,i+1,&SiPM);
		*(record_thresh+SiPM) = (clb_db_package[6+4*i]<<8)+clb_db_package[7+4*i];
	}
	SC_Channel2SiPM(s_SC,15,&SiPM);
	*(record_thresh+SiPM) = (clb_db_package[64]<<8)+clb_db_package[65];
	SC_Channel2SiPM(s_SC,16,&SiPM);
	*(record_thresh+SiPM) = (clb_db_package[68]<<8)+clb_db_package[69];
}
void WFCTAEvent::Deal22Package(long *single_count,float *DbTemp)
{
	s_SC = (clb_db_package[2]>>4)*10+clb_db_package[2]-(clb_db_package[2]>>4)*16;

	WDbTemp = (clb_db_package[4]<<8)+clb_db_package[5];
	*(DbTemp+DbTempCount) = s_SC;
	*(DbTemp+64+DbTempCount) = WDbTemp;
	DbTempCount++;

	for(int i=0;i<8;i++)
	{
		SC_Channel2SiPM(s_SC,i+1,&SiPM);
		*(single_count+SiPM) = (clb_db_package[9+i*5]<<32)+(clb_db_package[10+i*5]<<24)+(clb_db_package[11+i*5]<<16)+(clb_db_package[12+i*5]<<8)+clb_db_package[13+i*5];
	}
}
void WFCTAEvent::Deal23Package(long *single_count,long *single_time)
{
	s_SC = (clb_db_package[2]>>4)*10+clb_db_package[2]-(clb_db_package[2]>>4)*16;
	Wsingle_time = (clb_db_package[44]<<32)+(clb_db_package[45]<<24)+(clb_db_package[48]<<16)+(clb_db_package[49]<<8)+clb_db_package[50];
	for(int i=0;i<16;i++)
	{
		SC_Channel2SiPM(s_SC,i+1,&SiPM);
		*(single_time+SiPM) = Wsingle_time;
	}
	for(int i=8;i<16;i++)
	{
		SC_Channel2SiPM(s_SC,i+1,&SiPM);
		*(single_count+SiPM) = (clb_db_package[4+(i-8)*5]<<32)+(clb_db_package[5+(i-8)*5]<<24)+(clb_db_package[6+(i-8)*5]<<16)+(clb_db_package[7+(i-8)*5]<<8)+clb_db_package[8+(i-8)*5];
	}
}
void WFCTAEvent::Deal81Package(float *HV)
{
	s_SC = (clb_db_package[2]>>4)*10+clb_db_package[2]-(clb_db_package[2]>>4)*16;
	SC_Channel2SiPM(s_SC,1,&SiPM);
	*(HV+SiPM) = (clb_db_package[13]<<19)+(clb_db_package[14]<<11)+(clb_db_package[15]<<3)+(clb_db_package[16]>>5);
	SC_Channel2SiPM(s_SC,2,&SiPM);
	*(HV+SiPM) = (clb_db_package[16]-(clb_db_package[16]>>5)*32)*4194304+(clb_db_package[17]<<14)+(clb_db_package[18]<<6)+(clb_db_package[19]>>2);
	SC_Channel2SiPM(s_SC,3,&SiPM);
	*(HV+SiPM) = (clb_db_package[19]-(clb_db_package[19]>>2)*4)*33554432+(clb_db_package[20]<<17)+(clb_db_package[21]<<9)+(clb_db_package[22]<<1)+(clb_db_package[23]>>7);
	SC_Channel2SiPM(s_SC,4,&SiPM);
	*(HV+SiPM) = (clb_db_package[23]-(clb_db_package[23]>>7)*128)*1048576+(clb_db_package[24]<<12)+(clb_db_package[25]<<4)+(clb_db_package[26]>>4);
	SC_Channel2SiPM(s_SC,5,&SiPM);
	*(HV+SiPM) = (clb_db_package[26]-(clb_db_package[26]>>4)*16)*8388608+(clb_db_package[27]<<15)+(clb_db_package[28]<<7)+(clb_db_package[29]>>1);
	SC_Channel2SiPM(s_SC,6,&SiPM);
	*(HV+SiPM) = (clb_db_package[29]-(clb_db_package[29]>>1)*2)*67108864+(clb_db_package[30]<<18)+(clb_db_package[31]<<10)+(clb_db_package[32]<<2)+(clb_db_package[33]>>6);
	SC_Channel2SiPM(s_SC,7,&SiPM);
	*(HV+SiPM) = (clb_db_package[33]-(clb_db_package[33]>>6)*64)*2097152+(clb_db_package[34]<<13)+(clb_db_package[35]<<5)+(clb_db_package[36]>>3);
	SC_Channel2SiPM(s_SC,8,&SiPM);
 	*(HV+SiPM) = (clb_db_package[36]-(clb_db_package[36]>>3)*8)*16777216+(clb_db_package[37]<<16)+(clb_db_package[40]<<8)+(clb_db_package[41]);
	SC_Channel2SiPM(s_SC,9,&SiPM);
	*(HV+SiPM) = (clb_db_package[42]<<19)+(clb_db_package[43]<<11)+(clb_db_package[44]<<3)+(clb_db_package[45]>>5);
	SC_Channel2SiPM(s_SC,10,&SiPM);
	*(HV+SiPM) = (clb_db_package[45]-(clb_db_package[45]>>5)*32)*4194304+(clb_db_package[46]<<14)+(clb_db_package[47]<<6)+(clb_db_package[48]>>2);
	SC_Channel2SiPM(s_SC,11,&SiPM);
	*(HV+SiPM) = (clb_db_package[48]-(clb_db_package[48]>>2)*4)*33554432+(clb_db_package[49]<<17)+(clb_db_package[50]<<9)+(clb_db_package[51]<<1)+(clb_db_package[52]>>7);
	SC_Channel2SiPM(s_SC,12,&SiPM);
	*(HV+SiPM) = (clb_db_package[52]-(clb_db_package[52]>>7)*128)*1048576+(clb_db_package[53]<<12)+(clb_db_package[54]<<4)+(clb_db_package[55]>>4);
	SC_Channel2SiPM(s_SC,13,&SiPM);
	*(HV+SiPM) = (clb_db_package[55]-(clb_db_package[55]>>4)*16)*8388608+(clb_db_package[56]<<15)+(clb_db_package[57]<<7)+(clb_db_package[58]>>1);
	SC_Channel2SiPM(s_SC,14,&SiPM);
	*(HV+SiPM) = (clb_db_package[58]-(clb_db_package[58]>>1)*2)*67108864+(clb_db_package[59]<<18)+(clb_db_package[60]<<10)+(clb_db_package[61]<<2)+(clb_db_package[62]>>6);
	SC_Channel2SiPM(s_SC,15,&SiPM);
	*(HV+SiPM) = (clb_db_package[62]-(clb_db_package[62]>>6)*64)*2097152+(clb_db_package[63]<<13)+(clb_db_package[64]<<5)+(clb_db_package[65]>>3);
	SC_Channel2SiPM(s_SC,16,&SiPM);
	*(HV+SiPM) = (clb_db_package[65]-(clb_db_package[65]>>3)*8)*16777216+(clb_db_package[66]<<16)+(clb_db_package[67]<<8)+(clb_db_package[68]);

	for(int i=0;i<16;i++)
	{
		SC_Channel2SiPM(s_SC,i+1,&SiPM);
		*(HV+SiPM) /= (512*427.4087);
	}
}
void WFCTAEvent::Deal82Package(float *PreTemp)
{
	double A = 0.00433;
	double B = 13.582;
	double C;

	s_SC = (clb_db_package[2]>>4)*10+clb_db_package[2]-(clb_db_package[2]>>4)*16;
	for(int i=0;i<8;i++)
	{
		SC_Channel2SiPM(s_SC,i+1,&SiPM);
		*(PreTemp+SiPM) = (clb_db_package[13+i*2]<<8)+clb_db_package[14+i*2];
	}
	SC_Channel2SiPM(s_SC,9,&SiPM);
	*(PreTemp+SiPM) = (clb_db_package[13+8*2]<<8)+clb_db_package[14+9*2];
	for(int i=9;i<16;i++)
	{
		SC_Channel2SiPM(s_SC,i+1,&SiPM);
		*(PreTemp+SiPM) = (clb_db_package[13+(i+1)*2]<<8)+clb_db_package[14+(i+1)*2];
	}
	for(int i=0;i<16;i++)
	{
		SC_Channel2SiPM(s_SC,i+1,&SiPM);
		C = *(PreTemp+SiPM)*10000/32768-2230.8;
		*(PreTemp+SiPM) = (-1*B+sqrt(B*B-4*A*C))/(2*A);
		*(PreTemp+SiPM) += 30;
	}
}
void WFCTAEvent::Deal83Package(float *BigResistence)
{
	s_SC = (clb_db_package[2]>>4)*10+clb_db_package[2]-(clb_db_package[2]>>4)*16;
	for(int i=0;i<4;i++)
	{
		SC_Channel2SiPM(s_SC,i+1,&SiPM);
		*(BigResistence+SiPM) = ((clb_db_package[13+i*2]<<8)+clb_db_package[14+i*2])*50/256.;
	}
	SC_Channel2SiPM(s_SC,5,&SiPM);
	*(BigResistence+SiPM) = ((clb_db_package[13+4*2]<<8)+clb_db_package[14+5*2])*50/256.;
	for(int i=5;i<16;i++)
	{
		SC_Channel2SiPM(s_SC,i+1,&SiPM);
		*(BigResistence+SiPM) = ((clb_db_package[13+(i+1)*2]<<8)+clb_db_package[14+(i+1)*2])*50/256.;
	}
}
void WFCTAEvent::Deal84Package(float *SmallResistence)
{
	s_SC = (clb_db_package[2]>>4)*10+clb_db_package[2]-(clb_db_package[2]>>4)*16;
	SC_Channel2SiPM(s_SC,1,&SiPM);
	*(SmallResistence+SiPM) = ((clb_db_package[13]<<8)+clb_db_package[14+2])*10/1024.;
	for(int i=1;i<16;i++)
	{
		SC_Channel2SiPM(s_SC,i+1,&SiPM);
		*(SmallResistence+SiPM) = ((clb_db_package[13+(i+1)*2]<<8)+clb_db_package[14+(i+1)*2])*10/1024.;
	}

}
void WFCTAEvent::Deal85Package(long *ClbTime,float *ClbTemp)
{
	s_SC = (clb85package[2]>>4)*10+clb85package[2]-(clb85package[2]>>4)*16;
	WClbTime = 0;
	for(int i=10;i<15;i++)
	{
		WClbTime += (clb85package[i]<<(14-i)*8);
	}
	*(ClbTime+ClbTimeCount) = s_SC;
	*(ClbTime+64+ClbTimeCount) = WClbTime;
	ClbTimeCount++;


	WClbTemp = (clb85package[15]<<8)+clb85package[16];
	if(WClbTemp<40960)  {WClbTemp /= 256.;}
	else 			  {WClbTemp = (65536-WClbTemp)/256.;}
	*(ClbTemp+ClbTempCount) = s_SC;
	*(ClbTemp+64+ClbTempCount) = WClbTemp;
	ClbTempCount++;
/*
  if(s_SC==78){
  for(int i=0;i<74;i++){
  cout<<hex<<clb85package[i]<<" ";
  }
  cout<<endl;
  }
*/
}

void WFCTAEvent::DealFPGAPackage(int *fpgaVer)
{
	int FpgaNumber = fpgapackage[2];
        *(fpgaVer+FpgaNumber-1) = fpgapackage[61];
}
void WFCTAEvent::DealF9Package(int *fpgaVer)
{
	*(fpgaVer+8) = fpgapackage[61];
	*(fpgaVer+9) = fpgapackage[48];

	Wclb_initial_time = 0;
        for(int i=38;i<42;i++)
        {
                Wclb_initial_Time += (fpgapackage[i]<<(41-i)*8+6);
        }
        Wclb_initial_Time += (fpgapackage[42]>>2);
        Wclb_initial_time = ((fpgapackage[42]-(fpgapackage[42]>>2)*4)<<24)+(fpgapackage[43]<<16)+(fpgapackage[44]<<8)+(fpgapackage[45]<<0);

	//for(int i=38;i<46;i++)
	//{
	//	Wclb_initial_time += (fpgapackage[i]<<((45-i)*8));
	//}

	Wfired_tube = (fpgapackage[46]<<8)+fpgapackage[47];

        Wstatus_readback_Time = 0;
        for(int i=53;i<57;i++)
        {
                Wstatus_readback_Time += (fpgapackage[i]<<(56-i)*8+6);
        }
        Wstatus_readback_Time += (fpgapackage[57]>>2);
        Wstatus_readback_time = ((fpgapackage[57]-(fpgapackage[57]>>2)*4)<<24)+(fpgapackage[58]<<16)+(fpgapackage[59]<<8)+(fpgapackage[60]<<0);	

}

bool WFCTAEvent::EndFile()
{
	if(fInputStream.eof())
	{
		return true;
	}
	return false;
}

void WFCTAEvent::DealBigPackageHead()
{
	WEvent = 0;
	for(int i=8;i<12;i++)
	{
		WEvent += (bigpackagehead[i]<<((11-i)*8));
	}
}

void WFCTAEvent::ReadLittlePackage()
{
	pulseend = 6;
	portnumber = (littlepackage[pulseend-1]>>4);
	fpga = littlepackage[pulseend-1]-(littlepackage[pulseend-1]>>4)*16;
	SC = subcluster[portnumber-1][fpga-1];
	Channel = littlepackage[pulseend-2];
	basehigh = ((littlepackage[112]<<24)+(littlepackage[113]<<16)+(littlepackage[114]<<8)+(littlepackage[115]))/256.;
	baselow = ((littlepackage[116]<<24)+(littlepackage[117]<<16)+(littlepackage[118]<<8)+(littlepackage[119]))/256.;
	qhigh = ((littlepackage[104]<<24)+(littlepackage[105]<<16)+(littlepackage[106]<<8)+(littlepackage[107]))-basehigh*64;
	qlow = ((littlepackage[108]<<24)+(littlepackage[109]<<16)+(littlepackage[110]<<8)+(littlepackage[111]))-baselow*64;
	threshold_single = (littlepackage[120]<<8)+littlepackage[121];
	threshold_record = (littlepackage[122]<<8)+littlepackage[123];
	over_record = littlepackage[3]%2;
	over_single = (littlepackage[3]%4-littlepackage[3]%2)/2;
	peakpoint = littlepackage[2]*64+(littlepackage[3]>>2);

	WSC.push_back(SC);
	WChannel.push_back(Channel);
	Wpeak.push_back(peakpoint);
	WSingle_Threshold.push_back(threshold_single);
	WRecord_Threshold.push_back(threshold_record);
	WOver_Single_Marker.push_back(over_single);
	WOver_Record_Marker.push_back(over_record);
        WImageBaseHigh.push_back(basehigh);
        WImageBaseLow.push_back(baselow);
	WImageAdcHigh.push_back(qhigh);
	WImageAdcLow.push_back(qlow);
}

//void WFCTAEvent::LittlePackageCalc(double *pulseh, double *pulsel, int Telscope)
void WFCTAEvent::LittlePackageCalc(int *pulseh, int *pulsel)
{
	int pulsepointsize = 32;
        int iPulsehigh = 0;
        int iPulselow = 0;
        double sumhigh;
        double sumhighmax = 0;
	int Wpulsehigh[1024][32];//temp
	int Wpulselow[1024][32];//temp

        GAIN_MARKER = 0;
        Four_Point_ADC = 0;
        pulseend = 6;
        start = pulseend;
        maxpoint = 0;
        Qhigh = 0;
        Qlow = 0;
        Basehigh = 0;
        Baselow = 0;
        baselineBasehigh = 0;
        baselineBaselow = 0;

	portnumber = (littlepackage[pulseend-1]>>4);
	fpga = littlepackage[pulseend-1]-(littlepackage[pulseend-1]>>4)*16;
	SC = subcluster[portnumber-1][fpga-1];
	Channel = littlepackage[pulseend-2];
	SC_Channel2SiPM(SC,Channel,&SiPM);

/////////////////////////reade raw_data pulse/////////////////////////
	for(int i=start;i<start+52;i=i+3)
	{
		Wpulsehigh[SiPM][iPulsehigh] = (littlepackage[i]<<4)+(littlepackage[i+1]>>4);
		Wpulselow[SiPM][iPulselow] = (littlepackage[i+1]-(littlepackage[i+1]>>4)*16)*256+littlepackage[i+2];
		iPulsehigh++;iPulselow++;
		end = i+3;
	}
	Wpulsehigh[SiPM][iPulsehigh] = (littlepackage[end]<<4)+(littlepackage[end+1]>>4);
	Wpulselow[SiPM][iPulselow] = (littlepackage[end+1]-(littlepackage[end+1]>>4)*16)*256+littlepackage[end+4];
	iPulsehigh++;iPulselow++;
	for(int i=end+5;i<end+42;i=i+3)
	{
		Wpulsehigh[SiPM][iPulsehigh] = (littlepackage[i]<<4)+(littlepackage[i+1]>>4);
		Wpulselow[SiPM][iPulselow] = (littlepackage[i+1]-(littlepackage[i+1]>>4)*16)*256+littlepackage[i+2];
		iPulsehigh++;iPulselow++;
		pulseend = i+33;
	}
	end = pulseend;
//////////////////////////////////////////////put pulse into pulsehigh and pulselow//////////////////////////////////////////////
	for(int i=0;i<pulsepointsize;i++)
	{
		*(pulseh+SiPM*pulsepointsize+i) = Wpulsehigh[SiPM][i];
		*(pulsel+SiPM*pulsepointsize+i) = Wpulselow[SiPM][i];
	}

//////////////////////////////////////////////use High_Gain or Low_Gain//////////////////////////////////////////////
	for(int i=0;i<31;i++)
	{
		if(Wpulsehigh[SiPM][i]>4000) {GAIN_MARKER = 1;}
	}
//////////////////////////////////////////////trigger by myself//////////////////////////////////////////////
	for(int i=0;i<29;i++)
	{
		sumhigh=Wpulsehigh[SiPM][i]+Wpulsehigh[SiPM][i+1]+Wpulsehigh[SiPM][i+2]+Wpulsehigh[SiPM][i+3];
		sumhigh /= 4.;
		if(sumhighmax<sumhigh) {sumhighmax=sumhigh;maxpoint=i;}
	}
//	if(sumhighmax*4-4*basehigh>100){TUBETRIGGER = 1;WNpix++;}
	Four_Point_ADC = sumhighmax*4-4*basehigh;
//////////////////////////////////////////////find signal pulse Q_calculation//////////////////////////////////////////////
	if(maxpoint<8)
	{
		for(int i=31;i>21;i--)  {Basehigh += Wpulsehigh[SiPM][i];Baselow += Wpulselow[SiPM][i];}
		Basehigh = Basehigh/10.;
		Baselow = Baselow/10.;
//		for(int i=0;i<maxpoint+9;i++) {Qhigh += Wpulsehigh[SiPM][i]-basehigh;Qlow += Wpulselow[SiPM][i]-baselow;}
		for(int i=0;i<maxpoint+9;i++) {Qhigh += Wpulsehigh[SiPM][i]-Basehigh;Qlow += Wpulselow[SiPM][i]-Baselow;}
	}
	else if(maxpoint>22)
	{
		for(int i=0;i<10;i++)  {Basehigh += Wpulsehigh[SiPM][i];Baselow += Wpulselow[SiPM][i];}
		Basehigh = Basehigh/10.;
		Baselow = Baselow/10.;
//		for(int i=maxpoint-6;i<31;i++) {Qhigh += Wpulsehigh[SiPM][i]-basehigh;Qlow += Wpulselow[SiPM][i]-baselow;}
		for(int i=maxpoint-6;i<31;i++) {Qhigh += Wpulsehigh[SiPM][i]-Basehigh;Qlow += Wpulselow[SiPM][i]-Baselow;}
	}
	else
	{
		for(int i=0;i<5;i++)  {Basehigh += Wpulsehigh[SiPM][i];Baselow += Wpulselow[SiPM][i];}
		for(int i=31;i>26;i--)  {Basehigh += Wpulsehigh[SiPM][i];Baselow += Wpulselow[SiPM][i];}
		Basehigh = Basehigh/10.;
		Baselow = Baselow/10.;
//		for(int i=maxpoint-6;i<maxpoint+9;i++)  {Qhigh += Wpulsehigh[SiPM][i]-basehigh;Qlow += Wpulselow[SiPM][i]-baselow;}
		for(int i=maxpoint-6;i<maxpoint+9;i++)  {Qhigh += Wpulsehigh[SiPM][i]-Basehigh;Qlow += Wpulselow[SiPM][i]-Baselow;}
	}

	WSiPM.push_back(SiPM);
	Wgain_marker.push_back(GAIN_MARKER);
	Wmypeak.push_back(maxpoint);
	WADC_Cut.push_back(Four_Point_ADC);
//        WImageX.push_back(imagex);
//        WImageY.push_back(imagey);
	WmyImageBaseHigh.push_back(Basehigh);
	WmyImageBaseLow.push_back(Baselow);
	WmyImageAdcHigh.push_back(Qhigh);
	WmyImageAdcLow.push_back(Qlow);
}

void WFCTAEvent::DealBigPackageTail()
{
      WrabbitTime = 0;
      for(int i=8;i<12;i++)
      {
        WrabbitTime += (bigpackagetail[i]<<(11-i)*8+6);
      }
      WrabbitTime += (bigpackagetail[12]>>2);
      Wrabbittime = ((bigpackagetail[12]-(bigpackagetail[12]>>2)*4)<<24)+(bigpackagetail[13]<<16)+(bigpackagetail[14]<<8)+(bigpackagetail[15]<<0);
      Wpackagenum = bigpackagetail[6]*256+bigpackagetail[7];
}

void WFCTAEvent::SC_Channel2SiPM(short F_DB, short mChannel, short *mSiPM)
{   
    int FPGA;
    int DB;
    double SC_X;
    double SC_Y;
    double Channel_X;
    double Channel_Y;
    FPGA = F_DB%10;
    DB = F_DB/10;
    SC_X = ADRESS[DB-1][FPGA-1]%10;
    SC_Y = ADRESS[DB-1][FPGA-1]/10;

    if(mChannel<9)
    { 
      if((mChannel)%2==1){
        Channel_Y = (SC_X-1)*4.+1;
        Channel_X = 33-((SC_Y-1)*4.+(mChannel+1)/2.);
      }
      if((mChannel)%2==0){
        Channel_Y = (SC_X-1)*4.+2;
        Channel_X = 32.5-((SC_Y-1)*4.+(mChannel)/2.);
      }
    }
    if(mChannel>=9)
    { 
      if((mChannel)%2==1){
        Channel_Y = (SC_X-1)*4.+3;
        Channel_X = 33-((SC_Y-1)*4.+(mChannel-7)/2.);
      }
      if((mChannel)%2==0){
        Channel_Y = (SC_X-1)*4.+4;
        Channel_X = 32.5-((SC_Y-1)*4.+(mChannel-8)/2.);
      }
    }
    if(int(Channel_X*2)%2==1){Channel_X = Channel_X+0.5;}
    else {Channel_X = Channel_X;}
    Channel_Y = 32-Channel_Y;
    *mSiPM = int(1023-(Channel_Y*32+Channel_X-1));
}

long WFCTAEvent::GetEvent()
{
	return WEvent;
}

long WFCTAEvent::GetrabbitTime()
{
        return WrabbitTime;
}

double WFCTAEvent::Getrabbittime()
{
        return Wrabbittime;
}

short WFCTAEvent::Getpackagenum()
{
        return Wpackagenum;
}
