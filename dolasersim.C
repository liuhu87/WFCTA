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
   TTree *tree = new TTree("eventShow","info of laser evnets");

   WFTelescopeArray::jdebug=0;
   WFTelescopeArray::DoSim=true;
   WFTelescopeArray::GetHead(Form("/afs/ihep.ac.cn/users/h/hliu/Documents/LHAASO/WFCTA/default.inp"));
   Atmosphere::SetParameters();
   Atmosphere::scale=1.0e6;
   Laser::scale=1.0e-9;
   //Laser::spotrange = 0;//0.001;//0.001;//mm
   Laser::divergence = 1.;//0.0573;
   Laser::jdebug=3;
   Laser* pl=new Laser(seed);
   if(!pl->pwfc) pl->pwfc=new WFCTAEvent();
   WFCTAEvent* pevt=(pl->pwfc);
   printf("WFCTAEvent: %p\n",pevt);
   pl->SetParameters();
   pevt->CreateBranch(tree,1);

   int Time=10000;
   double time=0;
   for(int ii=0;ii<nevent;ii++){
      printf("ievent=%d Time=%d time=%lf\n",ii,Time,time);
      long int ngentel=pl->EventGen(Time,time,true);
      //fill the event
      hNevt->Fill(0.5,pl->count_gen);
      hNevt->Fill(1.5,ngentel/Laser::scale);
      hraytrace->Add(WFCTAMCEvent::hRayTrace);
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
   fout->Close();
   delete pl;
   //delete WFTelescopeArray::GetHead();

   return 0;
}
