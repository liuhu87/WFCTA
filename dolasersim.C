#include "WFTelescope.h"
#include <stdlib.h>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "Laser.h"
#include "WFCTAEvent.h"
#include "TelGeoFit.h"

using namespace std;

int main(int argc, char**argv)
{
   // arguments
   if (argc<3) {
      printf("   Usage %s <nevent> <outname> <seed>\n",argv[0]);
      return 0;
      //exit(0);
   }
   // vars 
   int nevent = atoi(argv[1]);
   const char* outname = argv[2];
   int seed = argc>3?atoi(argv[3]):0;
   cout<< nevent << "  " << outname << "  " << seed <<endl;

   TH1::SetDefaultSumw2();
   CommonTools::IsLaser=true;
   CommonTools::InitHArrival();

   TFile *fout = new TFile(outname,"RECREATE");
   TH1D* hNevt = new TH1D("NEvent","",2,0,2);
   TH1D* hraytrace = (TH1D*) WFCTAMCEvent::hRayTrace->Clone("RayTraceRes");
   hraytrace->Reset();
   TH1D* hPMTs = new TH1D("PMTs","",1024,-0.5,1023.5);
   TH1D* hweight = (TH1D*) WFCTALaserEvent::hweight->Clone("WeightRes");
   hweight->Reset();
   TTree *tree = new TTree("eventShow","info of laser evnets");

   WFTelescopeArray::jdebug=0;
   WFTelescopeArray::DoSim=true;
   //WFTelescopeArray::GetHead(Form("%s/default.inp",getenv("WFCTADataDir")));
   WFTelescopeArray::GetHead(Form("default.inp"));
   WFTelescope* pt=0;
   if(WFTelescopeArray::GetHead()){
      pt=WFTelescopeArray::GetHead()->pct[0];
   }
   if(!pt) return 0;
   //Atmosphere::SetParameters(Form("%s/default.inp",getenv("WFCTADataDir")));
   Atmosphere::SetParameters(Form("default.inp"));
   Atmosphere::GetHead()->AddATMModel(Form("/afs/ihep.ac.cn/users/h/hliu/Documents/LHAASO/WFCTA/ATMModel.txt"));
   Atmosphere::ATMRayModel=0;
   Atmosphere::ATMMieModel=-1;
   Atmosphere::scale=500.;

   Laser::UseTestScat=true;
   Laser::WhichRot=2;
   Laser::WhichTel=-1;
   Laser::scale=1.0e-8;
   //Laser::Doigen=9632;
   Laser::DoPlot=false;
   Laser::jdebug=3;
   Laser::TelSimDist=200.;
   Laser::IniRange[0][0]=-1.e3;
   Laser::IniRange[0][1]=2.0e5;
   Laser::IniRange[1][0]=-1.e4;
   Laser::IniRange[1][1]=1.0e4;
   Laser::IniRange[2][0]=-10;
   Laser::IniRange[2][1]=5.e5; //1.3e5
   //Laser::IniRange[3][0]=-1;
   //Laser::IniRange[3][1]=-1;

   Laser* pl=new Laser(seed);
   if(!pl->pwfc) pl->pwfc=new WFCTAEvent();
   WFCTAEvent* pevt=(pl->pwfc);
   printf("WFCTAEvent: %p\n",pevt);
   //pl->SetParameters(Form("%s/default.inp",getenv("WFCTADataDir")));
   pl->SetParameters(Form("default.inp"));
   double lasercoo0[3]={pl->lasercoo[0],pl->lasercoo[1],pl->lasercoo[2]};
   double laserdir0[2]={pl->laserdir[0],pl->laserdir[1]};
   pevt->CreateBranch(tree,1);

   int nlaserdir=10;

   int Time=1;
   double time=TelGeoFit::GetRotTime(Laser::WhichRot)/20.;
   for(int ii=0;ii<nevent;ii++){
      printf("ievent=%d Time=%d time=%lf\n",ii,Time,time);
      //double angle=-90.+180./nevent*ii;
      //pl->lasercoo[0]=1.0e5*cos(angle/180.*PI);
      //pl->lasercoo[1]=1.0e5*sin(angle/180.*PI);
      //pl->laserdir[0]=laserdir0[0]-5.+10./nlaserdir*(ii/nevent);
      //pl->laserdir[1]=laserdir0[1]-5.+10./nlaserdir*(ii%nevent);
      long int ngentel=pl->EventGen(Time,time,true);
      for(int itel=0;itel<pl->tellist.size();itel++){
         WFTelescope* pt=WFTelescopeArray::GetHead()->pct[pl->tellist.at(itel)];
         bool dosim=pl->DoWFCTASim(pl->tellist.at(itel));
         
         for(int ipmt=0;ipmt<NSIPM;ipmt++){
            double content=hPMTs->GetBinContent(ipmt+1);
            double econtent=hPMTs->GetBinError(ipmt+1);
            hPMTs->SetBinContent(ipmt+1,content+pevt->mcevent.TubeSignal[0][ipmt]);
            hPMTs->SetBinError(ipmt+1,sqrt(pow(econtent,2)+pow(content+pevt->mcevent.eTubeSignal[0][ipmt],2)));
         }

         double smax=0;
         //for(int jj=0;jj<NSIPM;jj++){
         //   if(pevt->mcevent.TubeSignal[0][jj]>smax) smax=pevt->mcevent.TubeSignal[0][jj];
         //}
         for(int jj=0;jj<pevt->iSiPM.size();jj++){
            if(pevt->LaserAdcH.at(jj)>smax) smax=pevt->LaserAdcH.at(jj);
         }
         printf("Evt=%d Ngen=%lf(%lf) iTel=%d ngentel=%ld evtn=%d size=%d maxsignal=%lf\n\n",ii,pevt->mcevent.Ngen,pl->count_gen,pt->TelIndex_,ngentel,pevt->iEvent,pevt->mcevent.RayTrace.size(),smax);
         tree->Fill();
         pevt->EventInitial();
      }
      if(Laser::DoPlot) pl->Draw("al",0,"./");
      //fill the event
      hNevt->Fill(0.5,pl->count_gen);
      hNevt->Fill(1.5,ngentel/Laser::scale);
      hraytrace->Add(WFCTAMCEvent::hRayTrace);
      hweight->Add(WFCTALaserEvent::hweight);
   }

   fout->cd();
   tree->Write();
   hNevt->Write();
   hraytrace->Write();
   hPMTs->Write();
   hweight->Write();
   for(int ii=0;ii<NSIPM;ii++){
      if(CommonTools::HArrival[ii]->Integral()>1000) CommonTools::HArrival[ii]->Write();
   }

   if(pl->hlength) pl->hlength->Write();
   if(pl->htheta) pl->htheta->Write();
   if(pl->hphi) pl->hphi->Write();

   fout->Close();
   //delete pl;
   //delete WFTelescopeArray::GetHead();

   return 0;
}
