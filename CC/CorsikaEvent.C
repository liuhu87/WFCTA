#include "CorsikaEvent.h"
#include "EventNtuple.h"
#include "ReadTrack.h"
void CorsikaEvent::Init(){
   run=-1;
   date=0;
   version=0;
   oheight=0;
   //nlevel=0;

   //EVTH
   event=-1;
   partp=0;
   ep=0;
   for(int ii=0;ii<3;ii++) pp[ii]=0;
   thetap=0;
   phip=0;
   stheight=0;
   for(int ii=0;ii<NuseMax;ii++) {corex[ii]=0; corey[ii]=0;}
   //EVTE
   nphoton=0;
   nelectron=0;
   nhadron=0;
   nmuon=0;
   nparticle=0;

   //RUNE
   nevent=0;

   int nsize=1000000;
   //CERENKOV Light
   nclight=0;
   //cx.resize(nsize);
   //cy.resize(nsize);
   //ct.resize(nsize);
   //cu.resize(nsize);
   //cv.resize(nsize);
   //height.resize(nsize);
   //wavelength.resize(nsize);
   //cweight.resize(nsize);

   cx.clear();
   cy.clear();
   ct.clear();
   cu.clear();
   cv.clear();
   height.clear();
   wavelength.clear();
   cweight.clear();
   //PARTICLE
   //ipart.resize(nsize);
   //igen.resize(nsize);
   //ilevel.resize(nsize);
   //px.resize(nsize);
   //py.resize(nsize);
   //pt.resize(nsize);
   //for(int ii=0;ii<3;ii++) ppp[ii].resize(nsize);
   //pweight.resize(nsize);

   ipart.clear();
   igen.clear();
   ilevel.clear();
   px.clear();
   py.clear();
   pt.clear();
   for(int ii=0;ii<3;ii++) ppp[ii].clear();
   pweight.clear();

   pwfc=0;
}
void CorsikaEvent::Reset(){
   run=-1;
   date=0;
   version=0;
   oheight=0;
   //nlevel=0;

   //EVTH
   event=-1;
   partp=0;
   ep=0;
   for(int ii=0;ii<3;ii++) pp[ii]=0;
   thetap=0;
   phip=0;
   stheight=0;
   for(int ii=0;ii<NuseMax;ii++) {corex[ii]=0; corey[ii]=0;}
   //EVTE
   nphoton=0;
   nelectron=0;
   nhadron=0;
   nmuon=0;
   nparticle=0;

   //RUNE
   nevent=0;
   //CERENKOV Light
   nclight=0;
   cx.clear();  
   cy.clear();
   ct.clear();
   cu.clear();
   cv.clear();
   height.clear();
   wavelength.clear();
   cweight.clear();

   //cx.swap(vector<float>(cx));
   //cy.swap(vector<float>(cy));
   //ct.swap(vector<float>(ct));
   //cu.swap(vector<float>(cu));
   //cv.swap(vector<float>(cv));
   //height.swap(vector<float>(height));
   //wavelength.swap(vector<float>(wavelength));
   //cweight.swap(vector<float>(cweight));

   //PARTICLE
   ipart.clear();
   igen.clear();
   ilevel.clear();
   px.clear();
   py.clear();
   pt.clear();
   for(int ii=0;ii<3;ii++) ppp[ii].clear();
   pweight.clear();

   //vector<int>(ipart).swap(ipart);
   //vector<int>(igen).swap(igen);
   //vector<int>(ilevel).swap(ilevel);
   //vector<int>(px).swap(px);
   //vector<int>(py).swap(py);
   //vector<int>(pt).swap(pt);
   //for(int ii=0;ii<3;ii++) vector<int>(ppp[ii]).swap(ppp[ii]);
   //vector<int>(pweight).swap(pweight);

   if(pwfc) pwfc->EventInitial();
}
void CorsikaEvent::Clear(){
   if(pwfc) delete pwfc;
}
void CorsikaEvent::Copy(CorsikaIO* pcorio){
   if(!pcorio) return;
   run=(pcorio->Evt).irun;
   date=(pcorio->Evt).idate;
   version=(pcorio->Evt).version;
   //nlevel=(pcorio->Evt).nlevel;
   oheight=(pcorio->Evt).height;

   event=(pcorio->Evt).ievent;
   partp=(pcorio->Evt).ipartp;
   ep=(pcorio->Evt).ep;
   pp[0]=(pcorio->Evt).pxp;
   pp[1]=(pcorio->Evt).pyp;
   pp[2]=(pcorio->Evt).pzp;
   thetap=(pcorio->Evt).thetap;
   phip=(pcorio->Evt).phip;
   stheight=(pcorio->Evt).stheight;
   for(int ii=0;ii<NuseMax;ii++) {corex[ii]=(pcorio->Evt).corex[ii]; corey[ii]=(pcorio->Evt).corey[ii];}

   nphoton=(pcorio->Evt).nphoton;
   nelectron=(pcorio->Evt).nelectron;
   nhadron=(pcorio->Evt).nhadron;
   nmuon=(pcorio->Evt).nmuon;
   nparticle=(pcorio->Evt).nparticle;

   nevent=(pcorio->Evt).nevent;

   if(pcorio->ftype==0||pcorio->ftype==2){
   nclight=(pcorio->Cer).nclight;
   cx.push_back((pcorio->Cer).x);
   cy.push_back((pcorio->Cer).y);
   ct.push_back((pcorio->Cer).t);
   cu.push_back((pcorio->Cer).u);
   cv.push_back((pcorio->Cer).v);
   height.push_back((pcorio->Cer).height);
   wavelength.push_back((pcorio->Cer).wavelength);
   cweight.push_back((pcorio->Cer).weight);
   }

   if(pcorio->ftype==1||pcorio->ftype==2){
   ipart.push_back((pcorio->Par).ipart);
   igen.push_back((pcorio->Par).igen);
   ilevel.push_back((pcorio->Par).ilevel);

   px.push_back((pcorio->Par).x);
   py.push_back((pcorio->Par).y);
   pt.push_back((pcorio->Par).t);
   ppp[0].push_back((pcorio->Par).px);
   ppp[1].push_back((pcorio->Par).py);
   ppp[2].push_back((pcorio->Par).pz);
   pweight.push_back((pcorio->Par).weight);
   }
}
int CorsikaEvent::WhichCore(double x0,double y0){
   int whichcore=-1;
   double mindist=1000000000;
   for(int ii=0;ii<Nuse;ii++){
      double dist=sqrt(pow(x0-corex[ii],2)+pow(y0-corey[ii],2));
      if(dist<mindist) { mindist=dist; whichcore=ii; }
   }
   return whichcore;
}
bool CorsikaEvent::DoWFCTASim(int iuse){
   if((iuse<0||iuse>=Nuse)&&(iuse!=100)){
      cerr<<"CorsikaEvent::DoWFCTASim: iuse ou out of range(iuse="<<iuse<<" Nuse="<<Nuse<<"), exiting..."<<endl;
      return false;
   }
   WFTelescopeArray* pct=WFTelescopeArray::GetHead();
   if(!pct) return false;
   if(!pct->CheckTelescope()) return false;
   else{//do the WFCTA simulation
      if(!pwfc) pwfc=new WFCTAEvent();
      for(int itel=0;itel<WFTelescopeArray::CTNumber;itel++) pct->GetCamera(itel)->ReSet();
      pwfc->EventInitial();
      int CERSize=cx.size();
      bool findcore=false;
      for(int icer=0;icer<CERSize;icer++){
         double x0,y0,z0;
         double m1,n1,l1;
         double weight,wave;
         double t=ct.at(icer)*1.0e-9; //from nano second to second
         int itube,icell;
         int whichtel;

         int whichcore=WhichCore(cx.at(icer),cy.at(icer));
         if(whichcore<0) continue;
         if((whichcore!=iuse)&&(iuse!=100)) continue;
         findcore=true;

         x0=cx.at(icer)-corex[whichcore];
         y0=cy.at(icer)-corey[whichcore];
         z0=0;
         ///the direction are verified with the shower image plot
         m1=cu.at(icer);
         n1=cv.at(icer);
         l1=(1-m1*m1-n1*n1>=0)?-sqrt(1-m1*m1-n1*n1):-1;
         weight=cweight.at(icer);
         wave=wavelength.at(icer);
         int res=pct->RayTrace(x0,y0,z0,m1,n1,l1,weight,wave,whichtel,t,itube,icell);
         if(WFCTAMCEvent::RecordRayTrace) (pwfc->mcevent).RayTrace.push_back(res);
         (pwfc->mcevent).hRayTrace->Fill(res,weight);
         (pwfc->mcevent).Ngen+=weight;

         //printf("CorsikaEvent::DoWFCTASim: run=%d event=%d ilight=%d cxy={%.0lf,%.0lf} corexy={%.0lf,%.0lf} coo={%lf,%lf} res=%d\n",run,event,icer,cx.at(icer),cy.at(icer),corex[whichcore],corey[whichcore],x0,y0,res);
         //printf("Corepos(Nuse=%d): ",Nuse);
         for(int ii=0;ii<Nuse;ii++){
            //printf("{%.0lf,%.0lf} ",corex[ii],corey[ii]);
         }
         //printf("\n");
      }
      if(!findcore) return false;
      if(pwfc){
         (pwfc->mcevent).Copy(pct);
         pwfc->CalculateDataVar();
         (pwfc->mcevent).GetTubeTrigger();
         (pwfc->mcevent).GetTelescopeTrigger(pct);
      }
      return true;
   }
}
void CorsikaEvent::Fill(){
   if(CorsikaIO::jdebug>0||true) printf("CorsikaEvent::Fill: Processing run=%d event=%d\n",run,event);
   if(EventNtuple::GetHead()){ //fill the data in EventNtuple
      EventNtuple::GetHead()->Fill(this,0);
      for(int iuse=0;iuse<Nuse;iuse++){
         if(DoWFCTASim(iuse)){
            pwfc->iEvent=event;
            pwfc->rabbitTime=CommonTools::Convert(date*1000000.+iuse);
            EventNtuple::GetHead()->Fill(0,pwfc,iuse);
         }
      }
   }
   if(ReadTrack::GetHead()&&ReadTrack::DoPlot){ //Copy Shower information to ReadTrack to do plot
      ReadTrack::GetHead()->Copy(this);
   }
   if(CorsikaIO::jdebug>=10) printf("CorsikaEvent::Fill: totally %d cer lights and %d particles\n",cx.size(),px.size());
   Reset();
}
int CorsikaEvent::CountParticle(){
   return TMath::Max(cx.size(),px.size());
}
