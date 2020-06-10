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
      printf("   Usage %s <nevent> <outname> <seed> <Flat_Scatter_Angle>\n",argv[0]);
      return 0;
      //exit(0);
   }
   // vars 
   int nevent = atoi(argv[1]);
   const char* outname = argv[2];
   int seed = argc>3?atoi(argv[3]):0;
   int flat_scatter = argc>4?atoi(argv[4]):-1;
   cout<< nevent << "  " << outname << "  " << seed <<endl;

   TH1::SetDefaultSumw2();
   CommonTools::IsLaser=true;
   CommonTools::InitHArrival();

   char cdir[200]="/afs/ihep.ac.cn/users/h/hliu/Documents/LHAASO/WFCTA";

   WFTelescopeArray::jdebug=0;
   WFTelescopeArray::DoSim=true;
   //WFTelescopeArray::GetHead(Form("%s/default.inp",getenv("WFCTADataDir")));
   WFTelescopeArray::GetHead(Form("%s/default.inp",cdir));
   WFTelescope* pt=0;
   if(WFTelescopeArray::GetHead()){
      pt=WFTelescopeArray::GetHead()->pct[0];
   }
   if(!pt) return 0;

   TFile *fout = new TFile(outname,"RECREATE");
   TH1D* hraytrace = (TH1D*) WFCTAMCEvent::hRayTrace->Clone("RayTraceRes");
   hraytrace->Reset();
   TH1D* hweightraytrace = (TH1D*) WFCTAMCEvent::hRayTrace->Clone("WeightRayTraceRes");
   hweightraytrace->Reset();
   TH1D* hPMTs[NCTMax];
   for(int itel=0;itel<NCTMax;itel++){
      int iTel=itel+1;
      int telindex=WFTelescopeArray::GetHead()->GetTelescope(iTel);
      if(telindex<0) hPMTs[itel]=0;
      else hPMTs[itel] = new TH1D(Form("PMTs_Tel%d",iTel),"",1024,-0.5,1023.5);
   }
   TH1D* hweight = (TH1D*) WFCTALaserEvent::hweight->Clone("WeightRes");
   hweight->Reset();
   TTree *tree = new TTree("eventShow","info of laser evnets");

   //Atmosphere::SetParameters(Form("%s/default.inp",getenv("WFCTADataDir")));
   Atmosphere::SetParameters(Form("%s/default.inp",cdir));
   Atmosphere::GetHead()->AddATMModel(Form("%s/ATMModel.txt",cdir));
   Atmosphere::ATMRayModel=0;
   Atmosphere::ATMMieModel=-1;
   Atmosphere::scale=500.;

   Laser::UseTestScat=(flat_scatter>0);
   Laser::WhichRot=2;
   Laser::WhichTel=-1;
   Laser::scale=1.0e-8;//-5.0e6;//4.0e-8;
   //Laser::Doigen=9632;
   Laser::DoPlot=false;
   Laser::jdebug=1;
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
   pl->SetParameters(Form("%s/default.inp",cdir));
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
         int telindex=pl->tellist.at(itel);
         if(telindex<0) continue;
         WFTelescope* pt=WFTelescopeArray::GetHead()->pct[telindex];
         int iTel=pt->TelIndex_;
         bool dosim=pl->DoWFCTASim(telindex);
         
         //fill the event
         hraytrace->Add(WFCTAMCEvent::hRayTrace);
         hweightraytrace->Add(WFCTAMCEvent::hWeightRayTrace);
         hweight->Add(WFCTALaserEvent::hweight);

         for(int ipmt=0;ipmt<NSIPM;ipmt++){
            double content=hPMTs[iTel-1]->GetBinContent(ipmt+1);
            double econtent=hPMTs[iTel-1]->GetBinError(ipmt+1);
            hPMTs[iTel-1]->SetBinContent(ipmt+1,content+pevt->mcevent.TubeSignal[telindex][ipmt]);
            hPMTs[iTel-1]->SetBinError(ipmt+1,sqrt(pow(econtent,2)+pow(pevt->mcevent.eTubeSignal[telindex][ipmt],2)));
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
      hraytrace->Fill(-20.);
      hweightraytrace->Fill(-20.);
      if(Laser::DoPlot) pl->Draw("al",0,"./");
   }

   fout->cd();
   tree->Write();
   hraytrace->Write();
   hweightraytrace->Write();
   for(int ii=0;ii<NCTMax;ii++){
      if(!hPMTs[ii]) continue;
      if(hPMTs[ii]->Integral()<=0) continue;
      hPMTs[ii]->Write();
   }
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
