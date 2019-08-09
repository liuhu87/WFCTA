#include "Laser.h"
#include "stdio.h"
#include <iostream>
#include "math.h"
#include "common.h"
#include "WFTelescope.h"
#include "WFCTAEvent.h"
using namespace std;

double Atmosphere::aod_air = 0;
double Atmosphere::aod_aerosol = 0;
double Atmosphere::scat_air = 0;
double Atmosphere::scat_aerosol = 0;

void Atmosphere::Init(int seed){
;
}
void Atmosphere::Release(){
;
}

void Atmosphere::SetParameters(char* filename){
   aod_air = 0.04*1.0e-5;
   scat_air = 0.03*1.0e-5;
}

bool Atmosphere::RayScatterAngle(double wavelength, double &theta, double &phi){
   if(!Laser::prandom) return false;
   phi=Laser::prandom->Uniform(0,2*PI);
   double costheta=Laser::prandom->Uniform(-1,1);
   theta=acos(costheta);
   return true;
}

bool Atmosphere::MieScatterAngle(double wavelength, double &theta, double &phi){
   if(!Laser::prandom) return false;
   phi=Laser::prandom->Uniform(0,2*PI);
   double costheta=Laser::prandom->Uniform(-1,1);
   theta=acos(costheta);
   return true;
}

double Atmosphere::ZDependence(double z,int type){
   return 1;
}
double Atmosphere::DeltaZ(double z){
   return 1000000.;
}
double Atmosphere::FreeIntgLength(){
   if(!Laser::prandom) return 0;
   if((aod_air+aod_aerosol)<0) return 0;
   else if(aod_air+aod_aerosol==0) return 1.0e10;
   double xx=Laser::prandom->Uniform();
   return log(1/(1-xx));
}
double Atmosphere::FreePathLength(double z0,double dir0[3]){
   double intglength=FreeIntgLength();
   if(intglength>0.9e10) return 1.0e10;
   if(dir0[2]==0){
      return intglength/((aod_air+aod_aerosol)*ZDependence(z0));
   }
   double norm=sqrt(pow(dir0[0],2)+pow(dir0[1],2)+pow(dir0[2],2));
   double integ=0;
   double length=0;
   int nstep=0;
   while(integ<intglength){
      //printf("step=%d: intglength=%le integ=%le z0=%le dirz=%lf\n",nstep,intglength,integ,z0,dir0[2]);
      //if(nstep>10) break;
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
      nstep++;
   }
   return -1;
}

int Atmosphere::IsScattering(double z0){
   if(!Laser::prandom) return 0;
   double sum=aod_air*ZDependence(z0,1)+aod_aerosol*ZDependence(z0,2);
   double sum_ext=aod_air*ZDependence(z0,1)-scat_air*ZDependence(z0,3)+aod_aerosol*ZDependence(z0,2)-scat_aerosol*ZDependence(z0,4);
   if(sum_ext<0) sum_ext=0;
   double sum_scat=scat_air*ZDependence(z0,3)+scat_aerosol*ZDependence(z0,4);
   sum=sum_ext+sum_scat;
   if(sum<0) return 0;
   else if(sum==0) return 100;
   double xx=Laser::prandom->Uniform();
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

int Laser::jdebug=0;
TRandom3* Laser::prandom = 0;
double Laser::TelSimDist=400.; //in cm
double Laser::scale=1.0e-4;
double Laser::unittime=1600.; //in ns
double Laser::intensity = 2;//mj
double Laser::intensity_err = 0;
double Laser::wavelength0 = 355;//nm
double Laser::wavelength0_err = 0;
double Laser::frequency = 1;
double Laser::pulsetime = 7.7e-9;
double Laser::spotrange = 1.07;//0.001;//mm, the range of the initial laser spot
double Laser::divergence = 1.0;//0.0573; //mrad

void Laser::Init(int seed){
   if(!prandom) prandom = new TRandom3();
   prandom->SetSeed(seed);
   pwfc=0;
   Reset();
   ievent_gen=0;
   count_gen=0;
   for(int ii=0;ii<3;ii++) lasercoo[ii]=0;
   for(int ii=0;ii<2;ii++) laserdir[ii]=0;
}
void Laser::Release(){
   if(prandom) delete prandom;
   if(pwfc) delete pwfc;
}
void Laser::Reset(){
   //count_gen=0;
   Time_gen=10000;
   time_gen=0;
   //ievent_gen=0;
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
void Laser::SetParameters(char* filename){
   lasercoo[0]=1000*100.; //in cm
   lasercoo[1]=0;
   lasercoo[2]=80.;
   laserdir[0]=90.;
   laserdir[1]=180.;
}

void Laser::cross(double dir1[3],double dir2[3],double *dir3){
   dir3[0]=dir1[1]*dir2[2]-dir1[2]*dir2[1];
   dir3[1]=dir1[2]*dir2[0]-dir1[0]*dir2[2];
   dir3[2]=dir1[0]*dir2[1]-dir1[1]*dir2[0];
}
bool Laser::CartesianFrame(double zero[3],double coor_in[3],double dir_in[3],double *xdir,double *ydir,double *zdir){
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
double Laser::mindist(double zero[3],double coor_in[3],double dir_in[3],double *coor_min,bool &decrease){
   double length=pow(dir_in[0],2)+pow(dir_in[1],2)+pow(dir_in[2],2);
   double kk=((zero[0]-coor_in[0])*dir_in[0]+(zero[1]-coor_in[1])*dir_in[1]+(zero[2]-coor_in[2])*dir_in[2])/length;
   decrease=(kk>=0);
   for(int ii=0;ii<3;ii++) coor_min[ii]=coor_in[ii]+kk*dir_in[ii];
   double dist=sqrt(pow(zero[0]-coor_min[0],2)+pow(zero[1]-coor_min[1],2)+pow(zero[2]-coor_min[2],2));
   return dist;
}
double Laser::mindist(double coor_in[3],double dir_in[3],int &whichtel,double *coor_min,bool &decrease){
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
      double coor_min0[3];
      double zero[3]={telcoor[ii][0],telcoor[ii][1],0};
      double dis=mindist(zero,coor_in,dir_in,coor_min0,decrease);
      //printf("mindist ct=%d: zero={%lf,%lf,%lf} coo={%lf,%lf,%lf} dist=%lf\n",ii,zero[0],zero[1],zero[2],coor_min0[0],coor_min0[1],coor_min0[2],dis);
      if(dis<min){
         whichtel=ii;
         min=dis;
         for(int jj=0;jj<3;jj++) coor_min[jj]=coor_min0[jj];
      }
   }
   if((!pta) && whichtel>=0) whichtel=-1;
   return min;
}
void Laser::PositionDis(double &xx,double &yy){
   double ran1=prandom->Uniform(0,spotrange*0.1);
   double rr=ran1;
   double phi=prandom->Uniform(0,2*PI);
   xx=rr*cos(phi);
   yy=rr*sin(phi);
}
void Laser::DirectionDis(double &theta,double &phi){
   double ran1=prandom->Uniform(cos(divergence*1.0e-3),1);
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
   double dir_in[3]={sin(laserdir[0]/180.*PI)*cos(laserdir[1]/180.*PI),sin(laserdir[0]/180.*PI)*sin(laserdir[1]/180.*PI),cos(laserdir[0]/180.*PI)};
   if(!CartesianFrame(zero,lasercoo,dir_in,xdir,ydir,zdir)) return false;

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

long int Laser::EventGen(int &Time,double &time){
   if(!prandom) return 0;
   double acctime=0;
   double time1=Time+time*20*1.0e-9;
   double time2=Time+time*20*1.0e-9+unittime*1.0e-9;
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
   acctime*=1.0e6;
   double ngen0=intensity*1.0e-3/(hplank*vlight/(wavelength0*1.0e-7))/pulsetime*1.0e-6*scale;

   Reset();
   long int ngen=(long int)(acctime*ngen0);
   Time_gen=Time;
   time_gen=time;
   long int ngentel=0;
   for(long int igen=0;igen<ngen;igen++){
      bool dogen=InitialGen();
      if(!dogen) continue;
      double time0=time1+(time2-time1)/ngen*(igen+0.5);
      count_gen+=1./scale;
      vgwav.push_back(wavelength_gen);
      for(int ii=0;ii<3;ii++){
         vgcoo[ii].push_back(coor_gen[ii]);
         vgdir[ii].push_back(dir_gen[ii]);
      }
      double distance;
      int res=Propagate(distance);
      //printf("Propagate igen=%d res=%d\n",igen,res);
      if(res<0) Telindex=res-20;
      else ngentel++;
      prodis.push_back(distance);
      votel.push_back(Telindex);
      for(int ii=0;ii<3;ii++){
         vocoo[ii].push_back(res>=0?coor_out[ii]:0);
         vodir[ii].push_back(res>=0?dir_out[ii]:0);
      }
   }
   if(jdebug>0) printf("Laser::EventGen: ngen0=%le acctime=%le ngen=%ld ngentel=%ld scale=%le\n",ngen0,acctime,ngen,ngentel,scale);
   //if(!exist) return 0;
   bool dosim=DoWFCTASim(1./scale);
   ievent_gen++;

   time+=unittime;
   if(time>=1.0){
      Time++;
      time-=1.0;
   }
   return ngentel;
}

int Laser::Propagate(double &distance){
   double coor_min[3];
   int whichtel;
   bool decrease;
   double freelength=Atmosphere::FreePathLength(coor_gen[2],dir_gen);
   double znew=coor_gen[2]+dir_gen[2]/sqrt(pow(dir_gen[0],2)+pow(dir_gen[1],2)+pow(dir_gen[2],2))*freelength;
   int scatter=Atmosphere::IsScattering(znew);
   if(scatter<=0){ //absorbed
      double dist=mindist(coor_gen,dir_gen,whichtel,coor_min,decrease);
      double dist2=sqrt(pow(coor_gen[0]-coor_min[0],2)+pow(coor_gen[1]-coor_min[1],2)+pow(coor_gen[2]-coor_min[2],2));
      distance=dist2;
      if(decrease&&dist<TelSimDist){
         if(dist2<freelength){
            Telindex=whichtel;
            for(int ii=0;ii<3;ii++){
               coor_out[ii]=coor_min[ii];
               dir_out[ii]=dir_gen[ii];
            }
            //printf("Propagate: absorb dist=%lf(free length=%lf) coo={%lf,%lf,%lf}\n",dist,freelength,coor_out[0],coor_out[1],coor_out[2]);
            return 0;
         }
         else return -2;
      }
      else{
         //printf("Propagate: absorb dist=%lf(free length=%lf)\n",dist,freelength);
         return -1;
      }
   }
   else if(scatter>2){ //no absorbtion, no scatter
      double dist=mindist(coor_gen,dir_gen,whichtel,coor_min,decrease);
      distance=sqrt(pow(coor_gen[0]-coor_min[0],2)+pow(coor_gen[1]-coor_min[1],2)+pow(coor_gen[2]-coor_min[2],2));
      if(decrease&&dist<TelSimDist){
         Telindex=whichtel;
         for(int ii=0;ii<3;ii++){
            coor_out[ii]=coor_min[ii];
            dir_out[ii]=dir_gen[ii];
         }
         //printf("Propagate: no absorb, no scatt dist=%lf(free length=%lf) coo={%lf,%lf,%lf}\n",dist,freelength,coor_out[0],coor_out[1],coor_out[2]);
         return 0;
      }
      else{
         //printf("Propagate: no absorb, no scatt dist=%lf(free length=%lf)\n",dist,freelength);
         return -1;
      }
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
      distance=freelength+sqrt(pow(coor_scat[0]-coor_min[0],2)+pow(coor_scat[1]-coor_min[1],2)+pow(coor_scat[2]-coor_min[2],2));
      if(decrease&&dist<TelSimDist){ //inside the field of view of one telescope after scattering
         double freelength2=Atmosphere::FreePathLength(coor_scat[2],dir_scat);
         if(dist<freelength2){
            Telindex=whichtel;
            for(int ii=0;ii<3;ii++){
               coor_out[ii]=coor_min[ii];
               dir_out[ii]=dir_scat[ii];
            }
            //printf("Propagate: scatter dist=%lf(free length=%lf) coo={%lf,%lf,%lf}\n",dist,freelength,coor_out[0],coor_out[1],coor_out[2]);
            return scatter;
         }
         else return -2;
      }
      else{
         //printf("Propagate: scatter dist=%lf(free length=%lf)\n",dist,freelength);
         return -1;
      }
   }
}

bool Laser::DoWFCTASim(double weight){
   WFTelescopeArray* pct=WFTelescopeArray::GetHead();
   if(!pct) return false;
   if(!pct->CheckTelescope()) return false;
   else{//do the WFCTA simulation
      if(!pwfc) pwfc=new WFCTAEvent();
      for(int itel=0;itel<WFTelescopeArray::CTNumber;itel++) pct->GetCamera(itel)->ReSet();
      pwfc->EventInitial();
      bool findtel=false;
      int LaserSize=(vocoo[0]).size();
      double avet=0;
      int nt=0;
      for(int il=0;il<LaserSize;il++){
         double x0,y0,z0;
         double m1,n1,l1;
         double wave;
         int itube,icell;
         int whichtel;

         whichtel=votel.at(il);
         if(whichtel<0||whichtel>=WFTelescopeArray::CTNumber){
            (pwfc->mcevent).RayTrace.push_back(-100+whichtel);
            //printf("pushing back tracing il=%d whichtel=%d\n",il,whichtel);
            continue;
         }
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
         wave=vgwav.at(il);
         double t=prodis.at(il)/vlight;
         int res=pct->RayTrace(x0,y0,z0,m1,n1,l1,weight,wave,whichtel,t,itube,icell);
         (pwfc->mcevent).RayTrace.push_back(res);
         if(res>0){
            avet+=t;
            nt++;
            if(jdebug>1) printf("Laser::DoWFCTASim: the photon go through the pmt %d, iphoton=%d\n",itube,il);
         }
         //printf("pushing back tracing il=%d whichtel=%d res=%d\n",il,whichtel,res);
      }
      if(pwfc){
         //printf("Laser::DoWFCTASim: Filling the event %d\n",ievent_gen);
         pwfc->iEvent=ievent_gen;
         double t0=(nt==0)?(Time_gen+time_gen*20*1.0e-9):(Time_gen+time_gen*20*1.0e-9+avet/nt);
         pwfc->rabbitTime=(int)t0; //should use the arrival time
         pwfc->rabbittime=(t0-pwfc->rabbitTime)/(20*1.0e-9); //should use the arrival time

         (pwfc->mcevent).Ngen=count_gen;
         (pwfc->mcevent).Timegen=Time_gen;
         //for(int ii=0;ii<3;ii++){
         //   (pwfc->mcevent).Coogen[ii].insert((pwfc->mcevent).Coogen[ii].begin(),vgcoo[ii].begin(),vgcoo[ii].end());
         //   (pwfc->mcevent).Dirgen[ii].insert((pwfc->mcevent).Dirgen[ii].begin(),vgdir[ii].begin(),vgdir[ii].end());
         //}
         //(pwfc->mcevent).Wavegen.insert((pwfc->mcevent).Wavegen.begin(),vgwav.begin(),vgwav.end());
         if(findtel){
            (pwfc->mcevent).Copy(pct);
            (pwfc->mcevent).GetTubeTrigger();
            (pwfc->mcevent).GetTelescopeTrigger(pct);
         }

         (pwfc->laserevent).Time=Time_gen;
         for(int ii=0;ii<3;ii++) (pwfc->laserevent).LaserCoo[ii]=lasercoo[ii];
         for(int ii=0;ii<2;ii++) (pwfc->laserevent).LaserDir[ii]=laserdir[ii];
         (pwfc->laserevent).wavelength=wavelength0;
         (pwfc->laserevent).Intensity=intensity;
         (pwfc->laserevent).Frequency=frequency;
         (pwfc->laserevent).pulsetime=pulsetime;
      }
      return true;
   }

   return false;
}
