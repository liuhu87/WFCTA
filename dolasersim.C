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

   TFile *fout = new TFile(outname,"RECREATE");
   TTree *tree = new TTree("eventShow","info of laser evnets");

   WFTelescopeArray::jdebug=0;
   WFTelescopeArray::DoSim=true;
   WFTelescopeArray::GetHead(Form("/afs/ihep.ac.cn/users/c/chenqh/sqn/WFCTA/default.inp"));
   Atmosphere::SetParameters();
   Laser::scale=1.0e-10;
   //Laser::spotrange = 0;//0.001;//0.001;//mm
   Laser::divergence = 1.;//0.0573;
   Laser::jdebug=2;
   Laser* pl=new Laser(seed);
   if(!pl->pwfc) pl->pwfc=new WFCTAEvent();
   WFCTAEvent* pevt=(pl->pwfc);
   printf("WFCTAEvent: %p\n",pevt);
   pl->SetParameters();
   pevt->CreateBranch(tree,1);

   int Time=10000;
   double time=0;
   for(int ii=0;ii<nevent;ii++){
      //printf("Event=%d\n",pl->ievent_gen);
      long int ngentel=pl->EventGen(Time,time);
      //fill the event
      double smax=0;
      for(int jj=0;jj<NSIPM;jj++){
         if(pevt->mcevent.TubeSignal[0][jj]>smax) smax=pevt->mcevent.TubeSignal[0][jj];
      }
      printf("Evt=%d Ngen=%ld(%ld) ngentel=%ld evtn=%d size=%d maxsignal=%lf\n\n",ii,pevt->mcevent.Ngen,pl->count_gen,ngentel,pevt->iEvent,pevt->mcevent.RayTrace.size(),smax);
      tree->Fill();
      pevt->EventInitial();
   }

   fout->cd();
   tree->Write();
   fout->Close();
   delete pl;
   //delete WFTelescopeArray::GetHead();

   return 0;
}
