#include "wcamera.h"
#include "TMath.h"
#include "TRandom.h"
#include <set>
WCamera::WCamera()
{

}

WCamera::~WCamera()
{

}

void WCamera::SetSiPMMAP()
{
  int  k;
  for(int i=0; i<PIX; i++){
     for(int j=0; j<PIX; j++){
         k = i*PIX + j;
         if(i%2==0)
            SiPMMAP[k][0] = (j+0.5-PIX/2.0)*D_ConeOut;

         if(i%2==1)
            SiPMMAP[k][0] = (j+1-PIX/2.0)*D_ConeOut;

         SiPMMAP[k][1] = (PIX/2.0-i)*D_ConeOut - D_ConeOut/2.0;
     }
  }
}

double WCamera::GetSiPMX(int itube)
{
  return SiPMMAP[itube][0];
}

double WCamera::GetSiPMY(int itube)
{
  return SiPMMAP[itube][1];
}

int WCamera::GetCone(double clusterx, double clustery)
{
  double deltax, deltay;
  int  itube;

  itube = -100;
  for(int k=0; k<NSIPM; k++){
     deltax = clusterx-SiPMMAP[k][0];
     deltay = clustery-SiPMMAP[k][1];
     if(fabs(deltax)<D_ConeIn/2.0&&fabs(deltay)<D_ConeIn/2.0){  //the gap between PMT are considered
        itube = k;
        break;
     }
  }
  return itube;
}

int WCamera::GetTube(double clusterx, double clustery)
{
  double deltax, deltay;
  int itube;

  itube = -100; 
  for(int k=0; k<NSIPM; k++)
  {
    deltax = clusterx-SiPMMAP[k][0];
    deltay = clustery-SiPMMAP[k][1];
    if(fabs(deltax)<D_SiPM/2.0&&fabs(deltay)<D_SiPM/2.0)
    {  //the gap between PMT are considered
       itube = k;
       break;
    }
  }
  return itube;
}

void WCamera::SetNSB(float nsb)
{
   NSB = nsb;
}

void WCamera::SetTriggerSigma(float triggersigma)
{
  TriggerSigma = triggersigma;
}

void WCamera::SetFixTriggerThreshold(float fixtriggerthreshold)
{
   FixTriggerThreshold = fixtriggerthreshold;
}

void WCamera::SetCTNumber(int ctnumber)
{  
   CTNumber = ctnumber;
}


void WCamera::Init()
{
  TubeSignal.resize(CTNumber*NSIPM);
  TubeTrigger.resize(CTNumber*NSIPM); 
  TubeSignalIntoCone.resize(CTNumber*NSIPM);
  TubeSignalAfterConeTracing.resize(CTNumber*NSIPM);
  TelTrigger.resize(CTNumber);
  //** Added by linglingMa on 2018-1-22 **//
  //** To Get the arrive time on each SiPM **//
  PeakTime.resize(CTNumber*NSIPM);
  TubeSignalInTriggerWindow.resize(CTNumber*NSIPM); 
}


void WCamera::ReSet(int ict)
{
   for(int i=0; i<NSIPM; i++){
      TubeSignal[ict*NSIPM+i] = 0;
      TubeTrigger[ict*NSIPM+i] = 0;
      TubeSignalIntoCone[ict*NSIPM+i] = 0;
      TubeSignalAfterConeTracing[ict*NSIPM+i] =0;
      //** Added by linglingMa on 2018-1-22 **//
      //** To Get the arrive time on each SiPM **//
      PeakTime[ict*NSIPM+i] = 0;
      TubeSignalInTriggerWindow[ict*NSIPM+i] = 0;  
   }
   TelTrigger[ict] = 0;

   cell.clear();
   ArriveTime.clear();
}


void WCamera::PhotonCellToTube()
{
  // TubeSignal[ict*NSIPM+itube]+= outpe ;
  int ict, itube;
  for(ct_iter=cell.begin();ct_iter!=cell.end();ct_iter++){
     ict = ct_iter->first;
     for(tube_iter=ct_iter->second.begin();tube_iter!=ct_iter->second.end();tube_iter++){
        itube = tube_iter->first;
        for(cell_iter=tube_iter->second.begin();cell_iter!=tube_iter->second.end();cell_iter++){
           TubeSignal[ict*NSIPM+itube] += cell_iter->second;
        } 
     }
  } 
}

void WCamera::GetPeakTime()
{
  //** Added by linglingMa on 2018-1-22 **//
  //** To Get the arrive time on each SiPM **//  
  int ict, itube, maxpe;
  for(ct_time_iter=ArriveTime.begin(); ct_time_iter!=ArriveTime.end(); ct_time_iter++){
     ict = ct_time_iter->first;
     for(tube_time_iter=ct_time_iter->second.begin();tube_time_iter!=ct_time_iter->second.end();tube_time_iter++){
        itube = tube_time_iter->first;
        maxpe=0;
        for(time_time_iter=tube_time_iter->second.begin();time_time_iter!=tube_time_iter->second.end();time_time_iter++){
           if(time_time_iter->second > maxpe) {
              maxpe = time_time_iter->second;
              PeakTime[ict*NSIPM+itube] = time_time_iter->first;
           }
           
        }  
    }
  }
  

}

void WCamera::GetPhotonInTriggerWindow(double trigger_window)
{ 
  int ict, itube;
  float delta_time;
  int Nphoton_in_triggerwindow;

  for(ct_time_iter=ArriveTime.begin(); ct_time_iter!=ArriveTime.end(); ct_time_iter++){
     ict = ct_time_iter->first;
     for(tube_time_iter=ct_time_iter->second.begin();tube_time_iter!=ct_time_iter->second.end();tube_time_iter++){
        itube = tube_time_iter->first;
        Nphoton_in_triggerwindow = 0;
        if( PeakTime[ict*NSIPM+itube] ==0) TubeSignalInTriggerWindow[ict*NSIPM+itube] = 0;
        else{
           for(time_time_iter=tube_time_iter->second.begin();time_time_iter!=tube_time_iter->second.end();time_time_iter++){
              delta_time = time_time_iter->first - PeakTime[ict*NSIPM+itube];
              if(fabs(delta_time)<=trigger_window/2.){
                 Nphoton_in_triggerwindow+=time_time_iter->second;
              }
           }
           TubeSignalInTriggerWindow[ict*NSIPM+itube] = Nphoton_in_triggerwindow;
        }
     }
   } 

}

void WCamera::PhotonIntoCone(int ict,int itube, int outpe)
{  
   TubeSignalIntoCone[ict*NSIPM+itube]+= outpe ;
}

void WCamera::GetArriveTime(int ict, int itube,int icell, int itime, int iphoton)
{
  if(cell[ict][itube][icell]==0) 
     ArriveTime[ict][itube][itime]+=iphoton; 
  else ArriveTime[ict][itube][itime] += 0;
}
void WCamera::PhotonIntoCell(int ict, int itube, int icell, int outpe)
{
   cell[ict][itube][icell] += outpe;
}

void WCamera::PhotonAfterConeTracing(int ict,int itube, int outpe)
{
   TubeSignalAfterConeTracing[ict*NSIPM+itube] += outpe;
}

void WCamera::AddNSB(int fadcflag)
{
  double nsb;
  for(int ict=0; ict<CTNumber; ict++){
     for(int itube=0; itube<NSIPM; itube++){
       nsb = gRandom->Poisson(NSB);
       if(fadcflag)  TubeSignalInTriggerWindow[ict*NSIPM+itube] +=int(nsb);
       else  TubeSignal[ict*NSIPM+itube] +=int(nsb);
     }
  }
}


void WCamera::GetTubeTrigger(int nsbflag, int fadcflag)
{
   int tubepe;
   if(nsbflag){
     for(int ict=0; ict<CTNumber; ict++){ 
        for(int itube=0; itube<NSIPM; itube++){
           if(fadcflag){
             tubepe = TubeSignalInTriggerWindow[ict*NSIPM+itube];
           } 
           else tubepe = TubeSignal[ict*NSIPM+itube];
           //Update By Lingling Ma 2019-5-9
           tubepe += gRandom->Poisson(NSB);
           if((tubepe-NSB)/sqrt(NSB)>TriggerSigma) TubeTrigger[ict*NSIPM+itube] = 1;
           else TubeTrigger[ict*NSIPM+itube] = 0;
        }
     }
   }
   else{
      for(int ict=0; ict<CTNumber; ict++){
        for(int itube=0; itube<NSIPM; itube++){
           if(fadcflag){
             tubepe = TubeSignalInTriggerWindow[ict*NSIPM+itube];
           }
           else tubepe = TubeSignal[ict*NSIPM+itube];
           if(tubepe>FixTriggerThreshold) TubeTrigger[ict*NSIPM+itube] = 1;
           else TubeTrigger[ict*NSIPM+itube] = 0;
        }
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


void WCamera::GetEulerMatrix(float TelZ,float TelA)
{
   float cosz, sinz, cosa, sina;
   cosz = cos(TelZ);
   sinz = sin(TelZ);
   cosa = cos(TelA);
   sina = sin(TelA);
   matrix_[0][0] = cosa*cosz;
   matrix_[0][1] = sina*cosz;
   matrix_[0][2] = -sinz;
   matrix_[1][0] = -sina;
   matrix_[1][1] = cosa;
   matrix_[1][2] = 0;
   matrix_[2][0] = cosa*sinz;
   matrix_[2][1] = sina*sinz;
   matrix_[2][2] = cosz;
}

void WCamera::InverseEuler(double x0, double y0, double z0, double *x, double *y, double *z)
{
   *x = matrix_[0][0]*x0+matrix_[1][0]*y0+matrix_[2][0]*z0;
   *y = matrix_[0][1]*x0+matrix_[1][1]*y0+matrix_[2][1]*z0;
   *z = matrix_[0][2]*x0+matrix_[1][2]*y0+matrix_[2][2]*z0;
}

void WCamera::GetTelescopeTrigger(int CTNumber,float *CT_Zen, float *CT_Azi)
{

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

  for(int iTel=0; iTel<CTNumber; iTel++)
  {
    GetEulerMatrix(CT_Zen[iTel],CT_Azi[iTel]);
    for(int ipmt = 0;ipmt<NSIPM;ipmt++)
    {
      // push all triggered pmts into a vector
      if(TubeTrigger[iTel*NSIPM+ipmt]==0) continue;
      x0=SiPMMAP[ipmt][0];
      y0=SiPMMAP[ipmt][1];
      z0 = FOCUS;
      u0=-x0/sqrt(x0*x0+y0*y0+z0*z0);
      v0=-y0/sqrt(x0*x0+y0*y0+z0*z0);
      w0=z0/sqrt(x0*x0+y0*y0+z0*z0);
      InverseEuler(u0,v0,w0,&u1,&v1,&w1);
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
           if(jpmt==ipmt||TubeTrigger[iTel*NSIPM+jpmt]==0) continue;
          x0=SiPMMAP[jpmt][0];
          y0=SiPMMAP[jpmt][1];
          z0 = FOCUS;
          u0=-x0/sqrt(x0*x0+y0*y0+z0*z0);
          v0=-y0/sqrt(x0*x0+y0*y0+z0*z0);
          w0=z0/sqrt(x0*x0+y0*y0+z0*z0);
          InverseEuler(u0,v0,w0,&u1,&v1,&w1) ;
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
  for(int iTel=0;iTel<CTNumber;iTel++) if(TelTrigger[iTel]!=0) return;

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
    if(Pattern.size()>=Trigger_multiple) for(int iTel=0;iTel<CTNumber;iTel++) TelTrigger[iTel]=2;
  }//  for(int ipmt = 0;ipmt<NPMT;ipmt++)
}
