#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <TFile.h>
#include "WFCTAEvent.h"
#include <TTree.h>
#include <cmath>
#include <map>

using namespace std;
int main(int argc,char** argv)
{
  if(argc<3) {
        printf("Usage: %s infile outfile  \n", argv[0]);
        return 1;
    }
  vector<short>* iSiPM = new vector<short>();
  long iEvent;
//  vector<int>* iSC = new vector<int>();
//  vector<int>* iChannel = new vector<int>();
  long rabbitTime;
  double rabbittime;
  vector<bool>* gain_marker = new vector<bool>();
  vector<char>* peak = new vector<char>();
  vector<char>* mypeak = new vector<char>();
  vector<short>* Single_Threshold = new vector<short>();
  vector<short>* Record_Threshold = new vector<short>();
  vector<bool>* Over_Single_Marker = new vector<bool>();
  vector<bool>* Over_Record_Marker = new vector<bool>();
  vector<float>* ADC_Cut = new vector<float>();//not electronics program, this trigger is done in this program
  vector<float>* ImageBaseHigh = new vector<float>();
  vector<float>* ImageBaseLow = new vector<float>();
  vector<float>* ImageAdcHigh = new vector<float>();
  vector<float>* ImageAdcLow = new vector<float>();
  vector<float>* myImageBaseHigh = new vector<float>();
  vector<float>* myImageBaseLow = new vector<float>();
  vector<float>* myImageAdcHigh = new vector<float>();
  vector<float>* myImageAdcLow = new vector<float>();
  int Npoint[32] = {0};  for(int i=0;i<32;i++) {Npoint[i]=i;}
  int pulsehigh[1024][32];
  int pulselow[1024][32];


  TFile *rootfile = new TFile(Form("%s.root",argv[argc-1]),"recreate");
/*********************************************************************/
                                            ////////////////////creat data trees////////////////////
  TTree *eventShow = new TTree("eventShow","Show Event");

  eventShow -> Branch("iEvent",&iEvent,"iEvent/L");
  eventShow -> Branch("rabbitTime",&rabbitTime,"rabbitTime/L");
  eventShow -> Branch("rabbittime",&rabbittime,"rabbittime/D");
  eventShow -> Branch("iSiPM","vector<short>",&iSiPM);
//  eventShow -> Branch("iSC","vector<int>",&iSC);
//  eventShow -> Branch("iChannel","vector<int>",&iChannel);
  eventShow -> Branch("gain_marker","vector<bool>",&gain_marker);
  eventShow -> Branch("peak","vector<char>",&peak);
  eventShow -> Branch("mypeak","vector<char>",&mypeak);
  eventShow -> Branch("Single_Threshold","vector<short>",&Single_Threshold);
  eventShow -> Branch("Record_Threshold","vector<short>",&Record_Threshold);
  eventShow -> Branch("Over_Single_Marker","vector<bool>",&Over_Single_Marker);
  eventShow -> Branch("Over_Record_Marker","vector<bool>",&Over_Record_Marker);
  eventShow -> Branch("ADC_Cut","vector<float>",&ADC_Cut);
  eventShow -> Branch("ImageBaseHigh","vector<float>",&ImageBaseHigh);
  eventShow -> Branch("ImageBaseLow","vector<float>",&ImageBaseLow);
  eventShow -> Branch("ImageAdcHigh","vector<float>",&ImageAdcHigh);
  eventShow -> Branch("ImageAdcLow","vector<float>",&ImageAdcLow);
  eventShow -> Branch("myImageBaseHigh","vector<float>",&myImageBaseHigh);
  eventShow -> Branch("myImageBaseLow","vector<float>",&myImageBaseLow);
  eventShow -> Branch("myImageAdcHigh","vector<float>",&myImageAdcHigh);
  eventShow -> Branch("myImageAdcLow","vector<float>",&myImageAdcLow);
//  eventShow -> Branch("Npoint",Npoint,"Npoint[32]/I");
//  eventShow -> Branch("pulsehigh",pulsehigh,"pulsehigh[1024][32]/I");
//  eventShow -> Branch("pulselow",pulselow,"pulselow[1024][32]/I");


  int fpgaVersion[10];
  long clb_initial_Time;
  double clb_initial_time;
  int fired_tube;
  long status_readback_Time;
  double status_readback_time;
  int sipm[1024];  for(int i=0;i<1024;i++) {sipm[i]=i;}
  short single_thresh[1024];
  short record_thresh[1024];
  long single_count[1024];
  float DbTemp[2][64];
  long single_time[1024];
  float HV[1024];
  float PreTemp[1024];
  float BigResistence[1024];
  float SmallResistence[1024];
  long ClbTime[2][64];
  float ClbTemp[2][64];

  TTree *Status = new TTree("Status","Status Tree");
  Status -> Branch("fpgaVersion",fpgaVersion,"fpgaVersion[10]/I");
  Status -> Branch("clb_initial_Time",&clb_initial_Time,"clb_initial_Time/L");
  Status -> Branch("clb_initial_time",&clb_initial_time,"clb_initial_time/D");
  Status -> Branch("fired_tube",&fired_tube,"fired_tube/I");
  Status -> Branch("status_readback_Time",&status_readback_Time,"status_readback_Time/L");
  Status -> Branch("status_readback_time",&status_readback_time,"status_readback_time/D");
  Status -> Branch("sipm",sipm,"sipm[1024]/I");
  Status -> Branch("single_thresh",single_thresh,"single_thresh[1024]/I");
  Status -> Branch("record_thresh",record_thresh,"record_thresh[1024]/I");
  Status -> Branch("single_count",single_count,"single_count[1024]/L");
  Status -> Branch("single_time",single_time,"single_time[1024]/L");
  Status -> Branch("DbTemp",DbTemp,"DbTemp[2][64]/F");
  Status -> Branch("HV",HV,"HV[1024]/F");
  Status -> Branch("PreTemp",PreTemp,"PreTemp[1024]/F");
  Status -> Branch("BigResistence",BigResistence,"BigResistence[1024]/F");
  Status -> Branch("SmallResistence",SmallResistence,"SmallResistence[1024]/F");
  Status -> Branch("ClbTime",ClbTime,"ClbTime[2][64]/L");
  Status -> Branch("ClbTemp",ClbTemp,"ClbTemp[2][64]/F");

  int FirstBigpackageAdress;
  int FPGA_marker;
  int Package_marker;
  WFCTAEvent* event = new WFCTAEvent();
/********************************************************************************/
//////////////////////////////////read file//////////////////////////////////
/********************************************************************************/
//open file and deal this file//
  event->OpenFile(argv[1]);
  FirstBigpackageAdress = event->FirstBigpackage();
  cout<<"FirstBigpackageAdress: "<<FirstBigpackageAdress<<endl;

  event->BranchInitial();
  event->StatusInitial();
  for(int i=0;i<1024;i++){
    single_thresh[i] = -1000;
    record_thresh[i] = -1000;
    single_count[i] = -1000;
    single_time[i] = -1000;
    HV[i] = -1000;
    PreTemp[i] = -1000;
    BigResistence[i] = -1000;
    SmallResistence[i] = -1000;
    for(int j=0;j<32;j++){
      pulsehigh[i][j] = 0;
      pulselow[i][j] = 0;
    }
  }
  for(int i=0;i<10;i++){
    fpgaVersion[i] = -1000;
  }

  while(true)
  {
    if(event->EndFile()){break;}

//deal big package head//
    if(event->Big_Package_Head())
    {
      event->DealBigPackageHead();
      iEvent = event->GetEvent();
		//      Littlepackagenum = 0;      Littlepackagenumhead = 0;      Littlepackagenumtail = 0;
    }
//deal little package//
    if(event->Little_Package())
    {
      event->ReadLittlePackage();
//      iSC = &(event->GetSC());
//      iChannel = &(event->GetChannel());
      peak = &(event->Getpeak());
      Single_Threshold = &(event->GetSingle_Threshold());
      Record_Threshold = &(event->GetRecord_Threshold());
      Over_Single_Marker = &(event->GetOver_Single_Marker());
      Over_Record_Marker = &(event->GetOver_Record_Marker());
      ImageBaseHigh = &(event->GetImageBaseHigh());
      ImageBaseLow = &(event->GetImageBaseLow());
      ImageAdcHigh = &(event->GetImageAdcHigh());
      ImageAdcLow = &(event->GetImageAdcLow());

      event->LittlePackageCalc((int *)pulsehigh, (int *)pulselow);
      iSiPM = &(event->GetSiPM());
      gain_marker = &(event->Getgain_marker());
      mypeak = &(event->Getmypeak());
      ADC_Cut = &(event->GetADC_Cut());
      myImageBaseHigh = &(event->GetmyImageBaseHigh());
      myImageBaseLow = &(event->GetmyImageBaseLow());
      myImageAdcHigh = &(event->GetmyImageAdcHigh());
      myImageAdcLow = &(event->GetmyImageAdcLow());
    }
//deal big package tail//
    if(event->Big_Package_Tail())
    {
      event->DealBigPackageTail();
      rabbitTime = event->GetrabbitTime();
      rabbittime = event->Getrabbittime();
      eventShow->Fill();
      event->BranchInitial();
      for(int i=0;i<1024;i++){
        for(int j=0;j<32;j++){
          pulsehigh[i][j] = 0;
          pulselow[i][j] = 0;
        }
      }
    }

//status package//
    if(event->Status_Package(&FPGA_marker,&Package_marker))
    {
//cout<<hex<<"fpga:"<<FPGA_marker<<" package:"<<Package_marker<<endl;
      if(FPGA_marker>10){
        switch(Package_marker){
	  case 0x21:
	    event->Deal21Package((short *)single_thresh, (short *)record_thresh);
	    break;
	  case 0x22:
            event->Deal22Package((long *)single_count,(float *)DbTemp);
	    break;
	  case 0x23:
            event->Deal23Package((long *)single_count,(long *)single_time);
	    break;
	  case 0x81:
            event->Deal81Package((float *)HV);
	    break;
          case 0x82:
            event->Deal82Package((float *)PreTemp);
            break;
          case 0x83:
            event->Deal83Package((float *)BigResistence);
            break;
          case 0x84:
            event->Deal84Package((float *)SmallResistence);
            break;
          case 0x85:
            event->Deal85Package((long *)ClbTime,(float *)ClbTemp);
	    break;
	}
      }
      if(FPGA_marker>0&&FPGA_marker<10){
        event->DealFPGAPackage((int *)fpgaVersion);
      }
      if(FPGA_marker==0){
        event->DealF9Package((int *)fpgaVersion);
        clb_initial_Time = event->GetclbInitialTime();
        clb_initial_time = event->GetclbInitialtime();
        fired_tube = event->GetFiredTube();
	status_readback_Time = event->GetStatusReadbackTime();
        status_readback_time = event->GetStatusReadbacktime();
        Status->Fill();
	event->StatusInitial();
        for(int i=0;i<1024;i++){
          single_thresh[i] = -1000;
          record_thresh[i] = -1000;
	  single_count[i] = -1000;
	  single_time[i] = -1000;
	  HV[i] = -1000;
	  PreTemp[i] = -1000;
	  BigResistence[i] = -1000;
	  SmallResistence[i] = -1000;
        }
        for(int i=0;i<10;i++){
          fpgaVersion[i] = -1000;
        }
      }
    }
  }
  event->CloseFile();
/******************************************************************************/
//////////////////////////////////write file//////////////////////////////////
/******************************************************************************/
  rootfile->Write();
  rootfile->Close();
}
