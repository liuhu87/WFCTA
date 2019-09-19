#include "WFTelescope.h"
#include <stdlib.h>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "Laser.h"
#include "WFCTAEvent.h"

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
   WFTelescopeArray::GetHead(Form("%s/default.inp",getenv("WFCTADataDir")));
   WFTelescope* pt=0;
   if(WFTelescopeArray::GetHead()){
      pt=WFTelescopeArray::GetHead()->pct[0];
   }
   if(!pt) return 0;
   Atmosphere::SetParameters(Form("%s/default.inp",getenv("WFCTADataDir")));
   Atmosphere::scale=1.0e3;
   Laser::scale=1.0e-8;
   //Laser::Doigen=9632;
   Laser::DoPlot=false;
   Laser::jdebug=3;
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
   pl->SetParameters(Form("%s/default.inp",getenv("WFCTADataDir")));
   pevt->CreateBranch(tree,1);

   int nlaserdir=10;
   double lasertheta=pl->laserdir[0];
   double laserphi=pl->laserdir[1];

   int Time=1;
   double time=0;
   for(int ii=0;ii<nevent;ii++){
      printf("ievent=%d Time=%d time=%lf\n",ii,Time,time);
      //pl->laserdir[0]=lasertheta+(ii/10)*3.;
      //pl->laserdir[1]=laserphi+(ii%10)*1.;
      long int ngentel=pl->EventGen(Time,time,true);
      if(Laser::DoPlot) pl->Draw("al",0,"./");
      //fill the event
      hNevt->Fill(0.5,pl->count_gen);
      hNevt->Fill(1.5,ngentel/Laser::scale);
      hraytrace->Add(WFCTAMCEvent::hRayTrace);
      hweight->Add(WFCTALaserEvent::hweight);
      for(int ipmt=0;ipmt<1024;ipmt++){
         double content=hPMTs->GetBinContent(ipmt+1);
         double econtent=hPMTs->GetBinError(ipmt+1);
         hPMTs->SetBinContent(ipmt+1,content+pevt->mcevent.TubeSignal[0][ipmt]);
         hPMTs->SetBinError(ipmt+1,sqrt(pow(econtent,2)+pow(content+pevt->mcevent.eTubeSignal[0][ipmt],2)));
      }

      double smax=0;
      for(int jj=0;jj<NSIPM;jj++){
         if(pevt->mcevent.TubeSignal[0][jj]>smax) smax=pevt->mcevent.TubeSignal[0][jj];
      }
      printf("Evt=%d Ngen=%lf(%lf) ngentel=%ld evtn=%d size=%d maxsignal=%lf\n\n",ii,pevt->mcevent.Ngen,pl->count_gen,ngentel,pevt->iEvent,pevt->mcevent.RayTrace.size(),smax);
      tree->Fill();
      pevt->EventInitial();
   }

   fout->cd();
   tree->Write();
   hNevt->Write();
   hraytrace->Write();
   hPMTs->Write();
   hweight->Write();
   fout->Close();
   //delete pl;
   //delete WFTelescopeArray::GetHead();

   return 0;
}
