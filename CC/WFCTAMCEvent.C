#include "TMath.h"
#include "TRandom.h"
#include <set>
#include "WFCTAMCEvent.h"
#include "EventNtuple.h"
bool WFCTAMCEvent::RecordRayTrace=false;
double WFCTAMCEvent::fAmpLow=(5.6/22);
double WFCTAMCEvent::fAmpHig=5.6;
TH1D* WFCTAMCEvent::hRayTrace=new TH1D("MCRayTrace","Ray Trace Result",35,-25-0.5,10-0.5);
void WFCTAMCEvent::Init(int size){
   //hRayTrace=new TH1D("MCRayTrace","Ray Trace Result",35,-25-0.5,10-0.5);
   //for(int ii=0;ii<3;ii++){
   //   Coogen[ii].resize(size>0?size:1000000);
   //   Dirgen[ii].resize(size>0?size:1000000);
   //}
   //Wavegen.resize(size>0?size:1000000);
   RayTrace.resize(size>0?size:1000000);
   Reset();
   //Ngen=0;
}
void WFCTAMCEvent::Reset(){
   Ngen=0;
   Timegen=0;
   //memset(this,0,sizeof(*this));
   if(hRayTrace) {if(hRayTrace->Integral()>0) hRayTrace->Reset();}
   //for(int ii=0;ii<3;ii++){
   //   Coogen[ii].clear();
   //   Dirgen[ii].clear();
   //}
   //Wavegen.clear();
   /*RayTrace.clear();
   for(int ict=0;ict<NCTMax;ict++){
      for(int ii=0;ii<NSIPM;ii++){
         TubeSignal[ict][ii]=0;
         eTubeSignal[ict][ii]=0;
         TubeTrigger[ict][ii]=0;
      }
      TelTrigger[ict]=0;

      NArrival[ict]=0;
      ArrivalTimeMin[ict]=1.0e20;
      ArrivalTimeMax[ict]=-1.0e20;
      OverFlow[ict]=false;
      for(int jj=0;jj<MaxTimeBin;jj++){
         ArrivalTime[ict][jj]=0;
         for(int ii=0;ii<NSIPM;ii++){
            ArrivalCount[ict][ii][jj]=0;
            ArrivalCountE[ict][ii][jj]=0;
         }
      }
   }*/
}
void WFCTAMCEvent::Copy(WFTelescopeArray* pct){
   if(!pct) pct=WFTelescopeArray::GetHead();
   if(!pct) return;
   for(int ict=0;ict<WFTelescopeArray::CTNumber;ict++){
      WCamera* pcame=pct->GetCamera(ict);
      if(!pcame){
         for(int ii=0;ii<NSIPM;ii++){
            TubeSignal[ict][ii]=0;
            eTubeSignal[ict][ii]=0;
         }
         continue;
      }
      pcame->AddNSB();
      NArrival[ict]=pcame->NArrival;
      ArrivalTimeMin[ict]=pcame->ArrivalTimeMin;
      ArrivalTimeMax[ict]=pcame->ArrivalTimeMax;
      OverFlow[ict]=pcame->OverFlow;
      for(int jj=0;jj<MaxTimeBin;jj++){
         ArrivalTime[ict][jj]=pcame->ArrivalTime[jj];
      }
      for(int ii=0;ii<NSIPM;ii++){
         TubeSignal[ict][ii]=(pcame->TubeSignal).at(ii);
         eTubeSignal[ict][ii]=(pcame->eTubeSignal).at(ii);
         int peaktime=-1;
         long int itmin=10000000000;
         long int itmax=-1;
         double maxcount=-1;
         for(int jj=0;jj<MaxTimeBin;jj++){
            ArrivalCount[ict][ii][jj]=pcame->ArrivalCount[ii][jj];
            ArrivalCountE[ict][ii][jj]=pcame->ArrivalCountE[ii][jj];
            if(ArrivalCount[ict][ii][jj]>maxcount){
               maxcount=ArrivalCount[ict][ii][jj];
               peaktime=jj;
            }
            if(ArrivalCount[ict][ii][jj]>0){
               if(ArrivalTime[ict][jj]<itmin) itmin=ArrivalTime[ict][jj];
               if(ArrivalTime[ict][jj]>itmax) itmax=ArrivalTime[ict][jj];
            }
         }
         //if(maxcount>5) printf("WFCTAMCEvent::Copy: T%d PMT%d signal=%lf ipeak=%d time={%ld,%ld,%ld} peakcount=%lf\n",ict,ii,TubeSignal[ict][ii],peaktime,ArrivalTime[ict][peaktime],itmin,itmax,maxcount);
      }
   }
}
void WFCTAMCEvent::GetTubeTrigger()
{
 for(int ict=0; ict<WFTelescopeArray::CTNumber; ict++){
    for(int itube=0; itube<NSIPM; itube++){
       if((TubeSignal[ict][itube]-WCamera::NSB)/sqrt(WCamera::NSB)>4) TubeTrigger[ict][itube] = 1;
         else TubeTrigger[ict][itube] = 0;
    }
 }
}
struct Dir{
  int ipmt;
  double u;
  double v;
  double w;
  
  //this line overloads '<' because the structure need be sorted in std::set!
   bool operator<(const Dir &a) const {return ipmt<a.ipmt;}
 };//This structure is only used below.
void WFCTAMCEvent::GetTelescopeTrigger(WFTelescopeArray* pct)
{
  if(!pct) pct=WFTelescopeArray::GetHead();
  if(!pct) return;
  double u0,v0,w0;
  double u1,v1,w1;
  double x0,y0,z0;
  double COS = cos(0.6*TMath::DegToRad());
  int NTrigger;
  struct Dir direction;
  vector<struct Dir> dir_collection;
  set<struct Dir> Pattern;
  set<struct Dir>::iterator itor;
  int Trigger_single = 3;
  int Trigger_multiple = 3;
  int addpmt;

  //Get directions of all triggered tubes.   if it is triggered push (u,v,w) into the vector[itel]
  //Calculate trigger for each telescope 

  for(int iTel=0; iTel<WFTelescopeArray::CTNumber; iTel++)
  {
    WFTelescope* p0=pct->pct[iTel];
    if(!p0) continue;
    for(int ipmt = 0;ipmt<NSIPM;ipmt++)
    {
      // push all triggered pmts into a vector
      if(TubeTrigger[iTel][ipmt]==0) continue;
      x0=WCamera::SiPMMAP[ipmt][0];
      y0=WCamera::SiPMMAP[ipmt][1];
      z0 = WFTelescope::FOCUS;
      u0=-x0/sqrt(x0*x0+y0*y0+z0*z0);
      v0=-y0/sqrt(x0*x0+y0*y0+z0*z0);
      w0=z0/sqrt(x0*x0+y0*y0+z0*z0);
      p0->InverseEuler(u0,v0,w0,&u1,&v1,&w1);
      direction.ipmt = iTel*NSIPM+ipmt;
      direction.u=u1;
      direction.v=v1;
      direction.w=w1;
      dir_collection.push_back(direction);

      // trigger for iTel
      Pattern.clear();
      Pattern.insert(direction);
      while(1)
      {
        addpmt = 0;
        for(int jpmt = 0;jpmt<NSIPM;jpmt++)
        {
           if(jpmt==ipmt||TubeTrigger[iTel][jpmt]==0) continue;
          x0=WCamera::SiPMMAP[jpmt][0];
          y0=WCamera::SiPMMAP[jpmt][1];
          z0 = WFTelescope::FOCUS;
          u0=-x0/sqrt(x0*x0+y0*y0+z0*z0);
          v0=-y0/sqrt(x0*x0+y0*y0+z0*z0);
          w0=z0/sqrt(x0*x0+y0*y0+z0*z0);
          p0->InverseEuler(u0,v0,w0,&u1,&v1,&w1) ;
          direction.ipmt = iTel*NSIPM+jpmt;
          direction.u=u1;
          direction.v=v1;
          direction.w=w1;
          for(itor=Pattern.begin();itor!=Pattern.end();itor++)
          {
            if(itor->ipmt==jpmt) continue;
            if(itor->u*u1+itor->v*v1+itor->w*w1>COS)
            {
              if(Pattern.insert(direction).second)
              {
                addpmt++;
                break;
              }
            }
          } // for(itor=Pattern.begin();itor!=Pattern.end();itor++)
          if(Pattern.size()>=Trigger_single) break;
        }//  for(int jpmt = 0;jpmt<NPMT;jpmt++)
        if(addpmt==0||Pattern.size()>=Trigger_single) break;
      } //while(1)
     if(Pattern.size()>=Trigger_single) {TelTrigger[iTel]=1; break;}

    }//  for(int ipmt = 0;ipmt<NPMT;ipmt++)
  }//  for(int iTel=0; iTel<ntel; iTel++)
  for(int iTel=0;iTel<WFTelescopeArray::CTNumber;iTel++) if(TelTrigger[iTel]!=0) return;

  // multiple trigger mode is applied only all telescopes are not triggered

  for(int ipmt=0;ipmt<dir_collection.size();ipmt++)
  {
    Pattern.clear();
    Pattern.insert(dir_collection[ipmt]);
    while(1)
    {
      addpmt = 0;
      for(int jpmt=0;jpmt<dir_collection.size();jpmt++)
      {
        if(jpmt==ipmt) continue;
        for(itor=Pattern.begin();itor!=Pattern.end();itor++)
        {
          if(itor->u*dir_collection[ipmt].u+itor->v*dir_collection[ipmt].v+itor->w*dir_collection[ipmt].w>COS)
          {
            if(Pattern.insert(dir_collection[ipmt]).second)
            {
              addpmt++;
              break;
            }
          }
        } // for(itor=Pattern.begin();itor!=Pattern.end();itor++)
        if(Pattern.size()>=Trigger_multiple) break;
      }//  for(int jpmt = 0;jpmt<NPMT;jpmt++)

      if(addpmt==0||Pattern.size()>=Trigger_multiple) break;
    } //while(1)
    if(Pattern.size()>=Trigger_multiple) for(int iTel=0;iTel<WFTelescopeArray::CTNumber;iTel++) TelTrigger[iTel]=2;
  }//  for(int ipmt = 0;ipmt<NPMT;ipmt++)
}
