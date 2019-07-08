#include "Laser.h"
#include "stdio.h"
#include <iostream>
#include "math.h"
#include "common.h"
#include "WFTelescope.h"
#include "WFCTAEvent.h"
using namespace std;

TRandom3* Atmosphere::prandom = 0;
double Atmosphere::aod_air = 0;
double Atmosphere::aod_aerosol = 0;
double Atmosphere::scat_air = 0;
double Atmosphere::scat_aerosol = 0;

void Atmosphere::Init(int seed){
   if(!prandom){
      prandom = new TRandom3();
   }
   prandom -> SetSeed(seed);
}
void Atmosphere::Release(){
   if(prandom) delete prandom;
}

void Atmosphere::SetParameters(char* filename){
;
}

void Atmosphere::RayScatterAngle(double wavelength, double &theta, double &phi){
///
;
}

void Atmosphere::MieScatterAngle(double wavelength, double &theta, double &phi){
///
;
}

double Atmosphere::ZDependence(double z,int type){
   return 1;
}
double Atmosphere::DeltaZ(double z){
   return 10.;
}
double Atmosphere::FreeIntgLength(){
   if(!prandom) return 0;
   if((aod_air+aod_aerosol)<0) return 0;
   else if(aod_air+aod_aerosol==0) return 1.0e10;
   double xx=prandom->Uniform();
   return log(1/(1-xx));
}
double Atmosphere::FreePathLength(double z0,double dir0[3]){
   double intglength=FreeIntgLength();
   if(dir0[2]==0){
      return intglength/((aod_air+aod_aerosol)*ZDependence(z0));
   }
   double norm=sqrt(pow(dir0[0],2)+pow(dir0[1],2)+pow(dir0[2],2));
   double integ=0;
   double length=0;
   while(integ<intglength){
      double z1=z0+DeltaZ(z0);
      double dintg=(aod_air+aod_aerosol)*ZDependence((z0+z1)/2)*fabs((z1-z0)/dir0[2]*norm);
      if(integ+dintg>=intglength){
         length+=(intglength-integ)/dintg*fabs((z1-z0)/dir0[2]*norm);
         return length;
      }
      else{
         integ+=dintg;
         length+=fabs((z1-z0)/dir0[2]*norm);
         z0=z1;
      }
   }
   return -1;
}

int Atmosphere::IsScattering(double z0){
   if(!prandom) return 0;
   double sum=aod_air*ZDependence(z0,1)+aod_aerosol*ZDependence(z0,2);
   double sum_ext=aod_air*ZDependence(z0,1)-scat_air*ZDependence(z0,3)+aod_aerosol*ZDependence(z0,2)-scat_aerosol*ZDependence(z0,4);
   if(sum_ext<0) sum_ext=0;
   double sum_scat=scat_air*ZDependence(z0,3)+scat_aerosol*ZDependence(z0,4);
   sum=sum_ext+sum_scat;
   if(sum<0) return 0;
   else if(sum==0) return 100;
   double xx=prandom->Uniform();
   if(xx<sum_ext/sum) return 0;  //absorbed
   else if(xx<sum_ext/sum+scat_air*ZDependence(z0,3)/sum) return 1; //Rayleigh Scattering
   else return 2;  //Mie Scattering
}

//double Atmosphere::probability(){
//   prob(distance) = (aod_air + aod_aerosol) * exp((-(aod_air + aod_aerosol)) * distance );//distance is the variable
//}


/*\
The class for laser photon propagation
*/

TRandom3* Laser::prandom = 0;
double Laser::TelSimDist=500.; //in cm
double Laser::unittime=1600.; //in ns
double Laser::intensity = 2;//mj
double Laser::intensity_err = 0;
double Laser::wavelength0 = 355;//nm
double Laser::wavelength0_err = 0;
double Laser::frequency = 1;
double Laser::pulsetime = 0.01;
double Laser::spotrange = 0.001;//mm, the range of the initial laser spot
double Laser::divergence = 0.0573;

void Laser::Init(int seed){
   if(!prandom) prandom = new TRandom3();
   prandom->SetSeed(seed);
   pwfc=0;
   Reset();
}
void Laser::Release(){
   if(prandom) delete prandom;
   if(pwfc) delete pwfc;
}
void Laser::Reset(){
   for(int ii=0;ii<3;ii++) lasercoo[ii]=0;
   for(int ii=0;ii<3;ii++) laserdir[ii]=0;
   count_gen=0;
   Time_gen=10000;
   time_gen=0;
   ievent_gen=0;
   wavelength_gen=0;
   vgwav.clear();
   for(int ii=0;ii<3;ii++){
      coor_gen[ii]=0;
      dir_gen[ii]=0;
      vgcoo[ii].clear();
      vgdir[ii].clear();
   }
   Telindex=-2;
   votel.clear();
   for(int ii=0;ii<3;ii++){
      coor_out[ii]=0;
      dir_out[ii]=0;
      vocoo[ii].clear();
      vodir[ii].clear();
   }
   if(pwfc) pwfc->EventInitial();
}

void Laser::cross(double dir1[3],double dir2[3],double dir3[3]){
   dir3[0]=dir1[1]*dir2[2]-dir1[2]*dir2[1];
   dir3[1]=dir1[2]*dir2[0]-dir1[0]*dir2[2];
   dir3[2]=dir1[0]*dir2[1]-dir1[1]*dir2[0];
}
bool Laser::CartesianFrame(double zero[3],double coor_in[3],double dir_in[3],double xdir[3],double ydir[3],double zdir[3]){
   //z axis
   for(int ii=0;ii<3;ii++) zdir[ii]=dir_in[ii];
   double norm=sqrt(zdir[0]*zdir[0]+zdir[1]*zdir[1]+zdir[2]*zdir[2]);
   if(norm<=0) return false;
   for(int ii=0;ii<3;ii++) zdir[ii]/norm;

   double dir1[3]={coor_in[0]-zero[0],coor_in[1]-zero[1],coor_in[2]-zero[2]};
   norm=sqrt(dir1[0]*dir1[0]+dir1[1]*dir1[1]+dir1[2]*dir1[2]);
   if(norm<=0) return false;

   //temporary x and y axis
   double xdir0[3],ydir0[3];
   cross(dir1,zdir,xdir0);
   norm=sqrt(xdir0[0]*xdir0[0]+xdir0[1]*xdir0[1]+xdir0[2]*xdir0[2]);
   if(norm<=0){
      double mult=dir1[0]*zdir[0]+dir1[1]*zdir[1]+dir1[2]*zdir[2];
      if(mult>0) return false;
      double zaxis[3]={0,0,1};
      cross(dir1,zaxis,xdir0);
      norm=sqrt(xdir0[0]*xdir0[0]+xdir0[1]*xdir0[1]+xdir0[2]*xdir0[2]);
   }
   if(norm<=0) return false;
   for(int ii=0;ii<3;ii++) xdir0[ii]/=norm;
   cross(zdir,xdir0,ydir0);
   norm=sqrt(ydir0[0]*ydir0[0]+ydir0[1]*ydir0[1]+ydir0[2]*ydir0[2]);
   for(int ii=0;ii<3;ii++) ydir0[ii]/=norm;

   //the new x and y axis, for which x point to the zero[3], and ydir=zdir X xdir
   double coor_min[3],dir_min[3];
   bool decrease;
   mindist(zero,coor_in,dir_in,coor_min,decrease);
   double mult=(zero[0]-coor_min[0])*ydir0[0]+(zero[1]-coor_min[1])*ydir0[1]+(zero[2]-coor_min[2])*ydir0[2];
   for(int ii=0;ii<3;ii++) xdir[ii]=(mult>=0?1:-1)*ydir0[ii];
   cross(zdir,xdir,ydir);
   return true;
}
double Laser::mindist(double zero[3],double coor_in[3],double dir_in[3],double coor_min[3],bool &decrease){
   double length=pow(dir_in[0],2)+pow(dir_in[1],2)+pow(dir_in[2],2);
   double kk=((zero[0]-coor_in[0])*dir_in[0]+(zero[1]-coor_in[1])*dir_in[1]+(zero[2]-coor_in[2])*dir_in[2])/length;
   decrease=(kk>=0);
   for(int ii=0;ii<3;ii++) coor_min[ii]=coor_in[ii]+kk*dir_in[ii];
   return sqrt(pow(zero[0]-coor_min[0],2)+pow(zero[1]-coor_min[1],2)+pow(zero[2]-coor_min[2],2));
}
double Laser::mindist(double coor_in[3],double dir_in[3],int &whichtel,double coor_min[3],bool &decrease){
   int nct;
   double telcoor[NCTMax][2];
   WFTelescopeArray* pta=WFTelescopeArray::GetHead();
   if(!pta){
      nct=1; telcoor[0][0]=telcoor[0][1]=0;
   }
   else{
      nct=pta->CTNumber;
      for(int ii=0;ii<nct;ii++){
         WFTelescope* pt=pta->pct[ii];
         telcoor[ii][0]=pt->Telx_;
         telcoor[ii][1]=pt->Tely_;
      }
   }
   whichtel=-2;
   //loop the telescope
   double min=1.0e10;
   for(int ii=0;ii<nct;ii++){
      double coor_min[3];
      double zero[3]={telcoor[ii][0],telcoor[ii][1],0};
      double dis=mindist(zero,coor_in,dir_in,coor_min,decrease);
      if(dis<min){
         whichtel=ii;
         min=dis;
      }
   }
   if((!pta) && whichtel>=0) whichtel=-1;
   return min;
}
void Laser::PositionDis(double &xx,double &yy){
   double ran1=prandom->Uniform(0,spotrange);
   double rr=ran1;
   double phi=prandom->Uniform(0,2*PI);
   xx=rr*cos(phi);
   yy=rr*sin(phi);
}
void Laser::DirectionDis(double &theta,double &phi){
   double ran1=prandom->Uniform(1-divergence,1);
   theta=acos(sqrt(ran1));
   phi=prandom->Uniform(0,2*PI);
}
bool Laser::InitialGen(){
   double norm=1;
   //direction
   double theta,phi;
   DirectionDis(theta,phi);
   double xdir[3],ydir[3],zdir[3];
   double zero[3]={0,0,0};
   if(!CartesianFrame(zero,lasercoo,laserdir,xdir,ydir,zdir)) return false;

   double rr=tan(theta);
   for(int ii=0;ii<3;ii++){
      dir_gen[ii]=zdir[ii]+rr*(cos(phi)*xdir[ii]+sin(phi)*ydir[ii]);
   }
   norm=sqrt(dir_gen[0]*dir_gen[0]+dir_gen[1]*dir_gen[1]+dir_gen[2]*dir_gen[2]);
   for(int ii=0;ii<3;ii++) dir_gen[ii]/=norm;

   //position
   double xx,yy;
   PositionDis(xx,yy);
   for(int ii=0;ii<3;ii++){
      coor_gen[ii]=lasercoo[ii]+xx*xdir[ii]+yy*ydir[ii];
   }
   
   //wavelength
   wavelength_gen=WaveLengthGen();
   return true;
}
double Laser::WaveLengthGen(){
   return prandom->Gaus(wavelength0,wavelength0_err);
}

int Laser::EventGen(int &Time,double &time){
   double ngen0=intensity*1.0e-3/(hplank*vlight/(wavelength0*1.0e-7))/pulsetime;
   double acctime=0;
   double time1=Time+time*20*1.0e-9;
   double time2=Time+time*20*1.0e-9+unittime;
   int nn1=(int)(time1*frequency);
   int nn2=(int)(time2*frequency);
   if(nn1==nn2){
      double xx1=TMath::Max(nn1/frequency,time1);
      double xx2=TMath::Min(nn1/frequency+pulsetime,time2);
      if(xx2>xx1) acctime+=(xx2-xx1);
   }
   else{
      acctime+=nn1/frequency+pulsetime-TMath::Min(time1,nn1/frequency+pulsetime);
      acctime+=TMath::Min(nn2/frequency+pulsetime,time2)-nn2/frequency;
      if(nn2>nn1+1) acctime+=(nn2-nn1-1)*pulsetime;
   }
   if(acctime<=0) return 0;

   Reset();
   int ngen=(int)(acctime*ngen0);
   count_gen=ngen;
   Time_gen=Time;
   time_gen=time;
   for(int igen=0;igen<ngen;igen++){
      bool dogen=InitialGen();
      if(!dogen) continue;
      count_gen++;
      vgwav.push_back(wavelength_gen);
      for(int ii=0;ii<3;ii++){
         vgcoo[ii].push_back(coor_gen[ii]);
         vgdir[ii].push_back(dir_gen[ii]);
      }
      bool res=Propagate();
      if(res<0) Telindex=-10;
      votel.push_back(Telindex);
      for(int ii=0;ii<3;ii++){
         vocoo[ii].push_back(res>=0?coor_out[ii]:0);
         vodir[ii].push_back(res>=0?dir_out[ii]:0);
      }
      DoWFCTASim();
   }
   ievent_gen++;

   time+=unittime;
   if(time>=1.0){
      Time++;
      time-=1.0;
   }
   return ngen;
}

int Laser::Propagate(){
   double coor_min[3];
   int whichtel;
   bool decrease;
   double freelength=Atmosphere::FreePathLength(coor_gen[2],dir_gen);
   double znew=coor_gen[2]+dir_gen[2]/sqrt(pow(dir_gen[0],2)+pow(dir_gen[1],2)+pow(dir_gen[2],2))*freelength;
   int scatter=Atmosphere::IsScattering(znew);
   if(scatter<=0){ //absorbed
      double dist=mindist(coor_gen,dir_gen,whichtel,coor_min,decrease);
      if(decrease&&dist<TelSimDist){
         double dist2=sqrt(pow(coor_gen[0]-coor_min[0],2)+pow(coor_gen[1]-coor_min[1],2)+pow(coor_gen[2]-coor_min[2],2));
         if(dist2<freelength){
            Telindex=whichtel;
            for(int ii=0;ii<3;ii++){
               coor_out[ii]=coor_min[ii];
               dir_out[ii]=dir_gen[ii];
            }
            return 0;
         }
         else return -2;
      }
      else return -1;
   }
   else if(scatter>2){ //no absorbtion, no scatter
      double dist=mindist(coor_gen,dir_gen,whichtel,coor_min,decrease);
      if(decrease&&dist<TelSimDist){
         Telindex=whichtel;
         for(int ii=0;ii<3;ii++){
            coor_out[ii]=coor_min[ii];
            dir_out[ii]=dir_gen[ii];
         }
         return 0;
      }
      else return -1;
   }
   else{ //scattering
      double coor_scat[3];
      double dir_scat[3];

      for(int ii=0;ii<3;ii++) coor_scat[ii]=coor_gen[ii]+dir_gen[ii]*freelength;

      double theta,phi;
      if(scatter==1){	//Rayleigh scattering
         Atmosphere::RayScatterAngle(wavelength_gen,theta,phi);
      }
      else if(scatter==2){	//Mie scattering
         Atmosphere::MieScatterAngle(wavelength_gen,theta,phi);
      }
      double xdir[3],ydir[3],zdir[3];
      double zero[3]={0,0,0};
      CartesianFrame(zero,coor_gen,dir_gen,xdir,ydir,zdir);
      double rr=tan(theta);
      for(int ii=0;ii<3;ii++){
         dir_scat[ii]=zdir[ii]+rr*(cos(phi)*xdir[ii]+sin(phi)*ydir[ii]);
      }
      double norm=sqrt(dir_scat[0]*dir_scat[0]+dir_scat[1]*dir_scat[1]+dir_scat[2]*dir_scat[2]);
      for(int ii=0;ii<3;ii++) dir_scat[ii]/=norm;

      //the distance to the closest telescope
      double dist=mindist(coor_scat,dir_scat,whichtel,coor_min,decrease);
      if(decrease&&dist<TelSimDist){ //inside the field of view of one telescope after scattering
         double freelength2=Atmosphere::FreePathLength(coor_scat[2],dir_scat);
         if(dist<freelength2){
            Telindex=whichtel;
            for(int ii=0;ii<3;ii++){
               coor_out[ii]=coor_min[ii];
               dir_out[ii]=dir_scat[ii];
            }
            return scatter;
         }
         else return -2;
      }
      else return -1;
   }
}

bool Laser::DoWFCTASim(){
   WFTelescopeArray* pct=WFTelescopeArray::GetHead();
   if(!pct) return false;
   if(!pct->CheckTelescope()) return false;
   else{//do the WFCTA simulation
      if(!pwfc) pwfc=new WFCTAEvent();
      for(int itel=0;itel<WFTelescopeArray::CTNumber;itel++) pct->GetCamera(itel)->ReSet();
      pwfc->EventInitial();
      bool findtel=false;
      int LaserSize=(vocoo[0]).size();
      for(int il=0;il<LaserSize;il++){
         double x0,y0,z0;
         double m1,n1,l1;
         double weight,wave;
         double t;
         int itube,icell;
         int whichtel;

         whichtel=votel.at(il);
         if(whichtel<0||whichtel>=WFTelescopeArray::CTNumber) continue;
         WFTelescope* pt=pct->pct[whichtel];
         if(!pt) continue;
         findtel=true;

         x0=vocoo[0].at(il)-pt->Telx_;
         y0=vocoo[1].at(il)-pt->Tely_;
         z0=vocoo[2].at(il);
         ///not quite sure about this direction
         m1=vodir[0].at(il);
         n1=vodir[1].at(il);
         l1=vodir[2].at(il);
         weight=1.;
         wave=vgwav.at(il);
         int res=pct->RayTrace(x0,y0,z0,m1,n1,l1,weight,wave,whichtel,t,itube,icell);
         (pwfc->mcevent).RayTrace.push_back(res);
      }
      if(!findtel) return false;
      if(pwfc){
         (pwfc->mcevent).Copy(pct);
         (pwfc->mcevent).GetTubeTrigger();
         (pwfc->mcevent).GetTelescopeTrigger(pct);
         (pwfc->laserevent).Frequency=frequency;
      }
      return true;
   }

   return false;
}
