#include "WFCamera.h"
#include "WFCone.h"
#include "TMath.h"
#include "TRandom.h"
#include "WFTelescope.h"
#include <set>
double WCamera::D_SiPM=15.0;
double WCamera::D_Cell = 2.5e-2;
int WCamera::NCell=600;
double WCamera::SiPMMAP[NSIPM][2];
float WCamera::NSB=0;
WCamera::WCamera()
{
   Init();
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
            SiPMMAP[k][0] = (j+0.5-PIX/2.0)*SquareCone::D_ConeOut;

         if(i%2==1)
            SiPMMAP[k][0] = (j+1-PIX/2.0)*SquareCone::D_ConeOut;

         SiPMMAP[k][1] = (PIX/2.0-i)*SquareCone::D_ConeOut - SquareCone::D_ConeOut/2.0;
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
     if(fabs(deltax)<SquareCone::D_ConeIn/2.0&&fabs(deltay)<SquareCone::D_ConeIn/2.0){  //the gap between PMT are considered
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

//void WCamera::SetCTNumber(int ctnumber)
//{  
//   CTNumber = ctnumber;
//}


void WCamera::Init()
{
  TubeSignal.resize(NSIPM);
  eTubeSignal.resize(NSIPM);
  ArrivalTimeMin.resize(NSIPM);
  ArrivalTimeMax.resize(NSIPM);
  ArrivalAccTime.resize(NSIPM);
  NArrival.resize(NSIPM);
  TubeTrigger.resize(NSIPM); 
}

void WCamera::ReSet()
{
   for(int i=0; i<NSIPM; i++){
      TubeSignal[i] = 0;
      eTubeSignal[i] = 0;
      ArrivalTimeMin[i]=0;
      ArrivalTimeMax[i]=0;
      ArrivalAccTime[i]=0;
      NArrival[i]=0;
      TubeTrigger[i] = 0;
   }
   TelTrigger = 0;
}

void WCamera::AddNSB()
{
  double nsb;
  for(int itube=0; itube<NSIPM; itube++){
    nsb = WFTelescopeArray::prandom->Poisson(NSB);
    //TubeSignal[itube] +=int(nsb);
    TubeSignal[itube] +=nsb;
  }
}

void WCamera::Fill(int itube,double time,double weight){
   if(itube<0||itube>=NSIPM) return;
   else{
      TubeSignal[itube] += weight;
      eTubeSignal[itube] = sqrt(pow(eTubeSignal[itube],2)+pow(weight,2));
      if(time!=0){
         if(ArrivalTimeMin[itube]!=0){
            if(time<ArrivalTimeMin[itube]) ArrivalTimeMin[itube]=time;
         }
         else ArrivalTimeMin[itube]=time;
         if(ArrivalTimeMax[itube]!=0){
            if(time>ArrivalTimeMax[itube]) ArrivalTimeMax[itube]=time;
         }
         else ArrivalTimeMax[itube]=time;

         ArrivalAccTime[itube]+=time;
         NArrival[itube]+=1;
      }
   }
}

void WCamera::GetTubeTrigger()
{
    for(int itube=0; itube<NSIPM; itube++){
       if((TubeSignal[itube]-NSB)/sqrt(NSB)>4) TubeTrigger[itube] = 1;
         else TubeTrigger[itube] = 0;
    }
}

/*struct Dir{
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

void WCamera::GetTelescopeTrigger(float CT_Zen, float *CT_Azi)
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
}*/
