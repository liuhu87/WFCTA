#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include "include/dumpPack.h"
#include "include/WFCTAEvent.h"
#include "include/WFCTADecode.h"
#include <TFile.h>
#include <TTree.h>

#define BUF_LEN 200000
#define STATUS_BUF_LEN 200000

using namespace std;

int main(int argc, char**argv)
{
  if(argc!=3)
  {
      printf("Use %s inputfile outfile\n",argv[0]);
      return 0;
  }

  FILE *fp;
  uint8_t *buf;// = new uint8_t[BUF_LEN];
  size_t size_of_read;
  fp = fopen(argv[1],"rb");
  bool statuspackloop;
  int32_t slicelength;
  short ITEL;

  int64_t packStart = 0;
  int64_t packSize = 0;
  int FEEDataHead;
  uint8_t status_pack_marker;

  short iTel;
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
  long single_time[1024];
  float DbTemp[1024];
  float HV[1024];
  float PreTemp[1024];
  float BigResistence[1024];
  float SmallResistence[1024];
  long ClbTime[1024];
  float ClbTemp[1024];


  //WFCTAEvent *wfctaEvent = new WFCTAEvent();
  TFile *rootfile = new TFile(argv[2],"recreate");
  /*********************************************************************/
  TTree *Status = new TTree("Status","Status Tree");
  Status -> Branch("iTel",&iTel,"iTel/S");
  Status -> Branch("fpgaVersion",fpgaVersion,"fpgaVersion[10]/I");
  Status -> Branch("clb_initial_Time",&clb_initial_Time,"clb_initial_Time/L");
  Status -> Branch("clb_initial_time",&clb_initial_time,"clb_initial_time/D");
  Status -> Branch("fired_tube",&fired_tube,"fired_tube/I");
  Status -> Branch("status_readback_Time",&status_readback_Time,"status_readback_Time/L");
  Status -> Branch("status_readback_time",&status_readback_time,"status_readback_time/D");
  Status -> Branch("sipm",sipm,"sipm[1024]/I");
  Status -> Branch("single_thresh",single_thresh,"single_thresh[1024]/S");
  Status -> Branch("record_thresh",record_thresh,"record_thresh[1024]/S");
  Status -> Branch("single_count",single_count,"single_count[1024]/L");
  Status -> Branch("single_time",single_time,"single_time[1024]/L");
  Status -> Branch("DbTemp",DbTemp,"DbTemp[1024]/F");
  Status -> Branch("HV",HV,"HV[1024]/F");
  Status -> Branch("PreTemp",PreTemp,"PreTemp[1024]/F");
  Status -> Branch("BigResistence",BigResistence,"BigResistence[1024]/F");
  Status -> Branch("SmallResistence",SmallResistence,"SmallResistence[1024]/F");
  Status -> Branch("ClbTime",ClbTime,"ClbTime[1024]/L");
  Status -> Branch("ClbTemp",ClbTemp,"ClbTemp[1024]/F");
  /*********************************************************************/

  WFCTADecode *wfctaDecode = new WFCTADecode();
  clb_initial_Time = -1000;
  clb_initial_time = -1000;
  fired_tube = -1000;
  status_readback_Time = -1000;
  status_readback_time = -1000;
  for(int i=0;i<1024;i++){
    single_thresh[i] = -1000;
    record_thresh[i] = -1000;
    single_count[i] = -1000;
    single_time[i] = -1000;
    DbTemp[i] = -1000; 
    HV[i] = -1000;
    PreTemp[i] = -1000;
    BigResistence[i] = -1000;
    SmallResistence[i] = -1000;
    ClbTime[i] = -1000;
    ClbTemp[i] = -1000;
  }
  for(int i=0;i<10;i++){
    fpgaVersion[i] = -1000;
  }
  while(true)
  {
      delete buf;
      buf = new uint8_t[40];
      size_of_read = fread((uint8_t *)buf,1,40,fp);
      if(size_of_read==0){break;}
      if(wfctaDecode->FEEDataFragment(buf))
      {
	FEEDataHead = wfctaDecode->feeDataHead();
        slicelength = wfctaDecode->sliceLength(buf,FEEDataHead);
        ITEL = wfctaDecode->Telid(buf,FEEDataHead);
        fseek(fp,-size_of_read+FEEDataHead,1);

        delete buf;
        buf = new uint8_t[slicelength];
        size_of_read = fread((uint8_t *)buf,1,slicelength,fp);
	//printf("slicelength:%lld\n",slicelength);
        packStart = 0;
	while(1)
	{
          //dumpPacket(buf,24,16);
          if(wfctaDecode->StatusPack(buf,int(size_of_read),packStart))
          {
	    statuspackloop = true;
	    while(statuspackloop)
	    {
  	      status_pack_marker = wfctaDecode->StatusPackCheck(buf,int(size_of_read),packStart);
              packSize = wfctaDecode->PackSize();
	      packStart = packSize;
	      printf("packSize:%lld | sizeofread:%lld | status_pack_marker:%x\n",packSize,size_of_read,status_pack_marker);
              //dumpPacket(buf,10,16);

              switch(status_pack_marker){
                case 0x21:
   		  wfctaDecode->Getthresh(buf,packSize,(short *)single_thresh, (short *)record_thresh);
                  break;
                case 0x22:
		  wfctaDecode->Deal22Pack(buf,packSize,(long *)single_count);
                  break;
                case 0x23:
                  wfctaDecode->Deal23Pack(buf,packSize,(long *)single_count,(long *)single_time);
                  break;
                case 0x81:
	  	  wfctaDecode->GetHV(buf,packSize,(float *)HV);
                  break;
                case 0x82:
		  wfctaDecode->GetPreTemp(buf,packSize,(float *)PreTemp);
                  break;
                case 0x83:
                  wfctaDecode->GetBigRes(buf,packSize,(float *)BigResistence);
                //  wfctaDecode->Deal83Package((float *)BigResistence);
                  break;
                case 0x84:
                  wfctaDecode->GetSmallRes(buf,packSize,(float *)SmallResistence);
                //  wfctaDecode->Deal84Package((float *)SmallResistence);
                  break;
                case 0x85:
	  	  wfctaDecode->GetClbTemp(buf,packSize,(float *)ClbTemp);
                  break;
	        case 0x9:
		  iTel = ITEL;
		  printf("itel:%d-----\n\n",iTel);
		  //dumpPacket(buf,24,16);
	  	  clb_initial_Time = wfctaDecode->GetclbInitialTime(buf,packSize);
		  clb_initial_time = wfctaDecode->GetclbInitialtime(buf,packSize);
		  fired_tube = wfctaDecode->GetFiredTube(buf,packSize);
		  status_readback_Time = wfctaDecode->GetStatusReadbackTime(buf,packSize);
		  status_readback_time = wfctaDecode->GetStatusReadbacktime(buf,packSize);

		  Status->Fill();
		  //iTel = -1;
		  clb_initial_Time = -1000;
		  clb_initial_time = -1000;
		  fired_tube = -1000;
		  status_readback_Time = -1000;
		  status_readback_time = -1000;
		  for(int i=0;i<1024;i++){
		    single_thresh[i] = -1000;
                    record_thresh[i] = -1000;
                    single_count[i] = -1000;
                    single_time[i] = -1000;
                    DbTemp[i] = -1000;
                    HV[i] = -1000;
                    PreTemp[i] = -1000;
                    BigResistence[i] = -1000;
                    SmallResistence[i] = -1000;
                    ClbTime[i] = -1000;
		    ClbTemp[i] = -1000;
                  }
                  for(int i=0;i<10;i++){
                    fpgaVersion[i] = -1000;
                  }
		  statuspackloop = false;
		  //packStart = packSize;
		  break;
		case 100:
		  iTel = ITEL;
                  Status->Fill();
                  clb_initial_Time = -1000;
                  clb_initial_time = -1000;
                  fired_tube = -1000;
                  status_readback_Time = -1000;
                  status_readback_time = -1000;
                  for(int i=0;i<1024;i++){
                    single_thresh[i] = -1000;
                    record_thresh[i] = -1000;
                    single_count[i] = -1000;
                    single_time[i] = -1000;
                    DbTemp[i] = -1000;
                    HV[i] = -1000;
                    PreTemp[i] = -1000;
                    BigResistence[i] = -1000;
                    SmallResistence[i] = -1000;
                    ClbTime[i] = -1000;
                    ClbTemp[i] = -1000;
                  }
                  for(int i=0;i<10;i++){
                    fpgaVersion[i] = -1000;
                  }
                  statuspackloop = false;
		  printf("statuspackloop:%d itel:%d---\n\n",statuspackloop,iTel);
		  break;
	      }
              //if(status_pack_marker>0&&status_pack_marker<10){
              //  wfctaDecode->DealFPGAPackage((int *)fpgaVersion);
              //}      }
	    }
            packStart = wfctaDecode->PackSize();
          }
          else
          {
		break;
              //fseek(fp,size_of_read,1);
          }
	}
      }
      else
      {
          fseek(fp,-size_of_read+20,1);
      }
  }
  fclose(fp);


/******************************************************************************/
  rootfile->Write();
  rootfile->Close();

}
