#include "Laser.h"
#include "stdio.h"
#include <iostream>
#include "math.h"
#include "common.h"
#include "WFTelescope.h"
#include "WFCamera.h"
#include "WFCTAEvent.h"
#include "TAxis3D.h"
#include "Readparam.h"
#include <TSystem.h>
using namespace std;

double Atmosphere::aod_air = 0;
double Atmosphere::aod_aerosol = 0;
double Atmosphere::scat_air = 0;
double Atmosphere::scat_aerosol = 0;
TGraph* Atmosphere::gRayScatAngle=0;
TGraph* Atmosphere::gMieScatAngle=0;
double Atmosphere::scale=1.0;

void Atmosphere::Init(int seed){
;
}
void Atmosphere::Release(){
   if(gRayScatAngle) {delete gRayScatAngle; gRayScatAngle=0;}
   if(gMieScatAngle) {delete gMieScatAngle; gMieScatAngle=0;}
}

void Atmosphere::SetParameters(char* filename){
   WReadConfig read;
   WReadConfig* readconfig=&read;
   readconfig->readparam(filename);
   aod_air = readconfig->Getaod_air();
   aod_aerosol = readconfig->Getaod_aerosol();
   scat_air = readconfig->Getscat_air();
   scat_aerosol = readconfig->Getscat_aerosol();
   printf("Atmosphere::SetParameters: load parameters from %s aod_air=%le aod_aerosol=%le scat_air=%le scat_aerosol=%le\n",filename,aod_air,aod_aerosol,scat_air,scat_aerosol);

   //aod_air = 0.04*1.0e-5;
   //scat_air = 0.039*1.0e-5;

   //if((!gRayScatAngle)||(!gMieScatAngle)){
   //   TFile* fin=TFile::Open(Form("%s/ScatterAngle.root",getenv("WFCTADataDir")));
   //   if(!fin) return;
   //   gRayScatAngle=fin?(TGraph*)fin->Get("RayScat"):0;
   //   gMieScatAngle=fin?(TGraph*)fin->Get("MieScat"):0;
   //   printf("Atmosphere: Scattering Database Initialed %p %p %p\n",fin,gRayScatAngle,gMieScatAngle);
   //   fin->Close();
   //}
}

double Atmosphere::ProbTransform(double xx,double yy[2],double &weight,bool IsCenter){
   double p1,p2,p3;
   //different p1,p2,p3 for different scale
   if(scale<5.){
      p1=p2=p3=1./3;
   }
   else if(scale<50.){
      if(IsCenter){
         p2=(1-1./scale*2.5);
         p3=(1-p2)*(1-yy[1])/(1-(yy[1]-yy[0]));
         p1=1-p2-p3;
      }
      else{
         p3=1./scale;
         p1=1-p3-1.5*p3;
         p2=1-p1-p3;
      }
   }
   else if(scale<5.e3){
      if(IsCenter){
         p2=(1-1./scale)*0.8;
         p3=(1-p2)*(1-yy[1])/(1-(yy[1]-yy[0]));
         p1=1-p2-p3;
      }
      else{
         p3=1./scale;
         p1=(1-p3)*0.8;
         p2=1-p1-p3;
      }
   }
   else if(scale<5.e5){
      if(IsCenter){
         p2=(1-1./scale)*0.95;
         p3=(1-p2)*(1-yy[1])/(1-(yy[1]-yy[0]));
         p1=1-p2-p3;
      }
      else{
         p3=1./scale;
         p1=(1-p3)*0.95;
         p2=1-p1-p3;
      }
   }
   else{
      if(IsCenter){
         p2=(1-1./scale)*(1-1./scale);
         p3=(1-p2)*(1-yy[1])/(1-(yy[1]-yy[0]));
         p1=1-p2-p3;
         //p2=1.; p3=p1=0;
      }
      else{
         p3=1./scale;
         p1=(1-p3)*(1-p3);
         p2=1-p1-p3;
         //p1=1.; p2=p3=0;
      }
   }

   if(xx<=p1){
      weight*=(yy[0]/p1);
      return (xx==p1)?yy[0]:(yy[0]/p1*xx);
   }
   else if(xx<=p1+p2){
      weight*=(yy[1]-yy[0])/p2;
      return yy[0]+(yy[1]-yy[0])/p2*(xx-p1);
   }
   else{
      weight*=(1-yy[1])/p3;
      return yy[1]+(1-yy[1])/p3*(xx-p1-p2);
   }
}
bool Atmosphere::RayScatterAngleTheta(double wavelength, double &theta, double anglerange[2],double &weight){
   if(!Laser::prandom) return false;
   //if(!gRayScatAngle){
   //   printf("Atmosphere::RayScatterAngle: No Ray Scatter Angle Calculated %p\n",gRayScatAngle);
   //   return false;
   //}

   double xxx;
   double yrange[2];
   xxx=Laser::prandom->Uniform(0,1);
   if(anglerange[1]<=anglerange[0]){
      //theta=gRayScatAngle->Eval(xxx);
      theta=xxx*PI;
      weight*=2./3.*(1+pow(cos(theta),2));
      return true;
   }
   else{
      yrange[0]=anglerange[0]/PI;
      yrange[1]=anglerange[1]/PI;
      theta=ProbTransform(xxx,yrange,weight,yrange[0]>0); //from 0 to 1
      if(theta<0) return false;
      //theta=gRayScatAngle->Eval(theta);
      theta*=PI;
      weight*=2./3.*(1+pow(cos(theta),2));
      return true;
   }
}
bool Atmosphere::RayScatterAnglePhi(double wavelength, double &phi,double anglerange[2],double &weight){
   if(!Laser::prandom) return false;
   //if(!gRayScatAngle){
   //   printf("Atmosphere::RayScatterAngle: No Ray Scatter Angle Calculated %p\n",gRayScatAngle);
   //   return false;
   //}

   double xxx;
   double yrange[2];
   xxx=Laser::prandom->Uniform(0,1);
   if(anglerange[1]<=anglerange[0]){
      phi=xxx*2*PI-PI;
      return true;
   }
   else{
      if(anglerange[1]<-PI){
         anglerange[0]+=2*PI;
         anglerange[1]+=2*PI;
      }
      if(anglerange[0]>PI){
         anglerange[0]-=2*PI;
         anglerange[1]-=2*PI;
      }
      yrange[0]=anglerange[0]/PI/2+0.5;
      yrange[1]=anglerange[1]/PI/2+0.5;

      phi=ProbTransform(xxx,yrange,weight,yrange[0]>0); //from 0 to 1
      if(phi<0) return false;
      phi=(phi-0.5)*2*PI;
      return true;
   }
}
bool Atmosphere::MieScatterAngleTheta(double wavelength, double &theta, double anglerange[2],double &weight){
   if(!Laser::prandom) return false;
   //if(!gMieScatAngle){
   //   printf("Atmosphere::MieScatterAngle: No Mie Scatter Angle Calculated\n");
   //   return false;
   //}

double phase_init[181]={
    0.591448, 0.594184, 0.594692, 0.593188, 0.589878, 0.584950,
    0.578583, 0.570940, 0.562176, 0.552432, 0.541840, 0.530522,
    0.518590, 0.506148, 0.493292, 0.480108, 0.466677, 0.453072,
    0.439359, 0.425599, 0.411846, 0.398150, 0.384554, 0.371098,
    0.357817, 0.344743, 0.331900, 0.319315, 0.307006, 0.294990,
    0.283283, 0.271895, 0.260836, 0.250113, 0.239731, 0.229694,
    0.220003, 0.210658, 0.201657, 0.193000, 0.184681, 0.176696,
    0.169040, 0.161706, 0.154689, 0.147979, 0.141570, 0.135453,
    0.129620, 0.124061, 0.118768, 0.113732, 0.108942, 0.104390,
    0.100067, 0.095962, 0.092068, 0.088374, 0.084872, 0.081554,
    0.078409, 0.075431, 0.072611, 0.069941, 0.067412, 0.065019,
    0.062753, 0.060608, 0.058577, 0.056654, 0.054832, 0.053106,
    0.051470, 0.049919, 0.048449, 0.047053, 0.045728, 0.044470,
    0.043275, 0.042138, 0.041056, 0.040027, 0.039046, 0.038111,
    0.037220, 0.036369, 0.035557, 0.034782, 0.034040, 0.033331,
    0.032653, 0.032004, 0.031383, 0.030788, 0.030219, 0.029673,
    0.029151, 0.028650, 0.028172, 0.027713, 0.027274, 0.026855,
    0.026454, 0.026071, 0.025706, 0.025357, 0.025026, 0.024710,
    0.024411, 0.024127, 0.023859, 0.023605, 0.023367, 0.023143,
    0.022933, 0.022738, 0.022556, 0.022388, 0.022234, 0.022093,
    0.021966, 0.021851, 0.021750, 0.021661, 0.021585, 0.021522,
    0.021471, 0.021432, 0.021405, 0.021391, 0.021388, 0.021397,
    0.021418, 0.021450, 0.021494, 0.021548, 0.021614, 0.021691,
    0.021779, 0.021878, 0.021987, 0.022106, 0.022236, 0.022375,
    0.022525, 0.022684, 0.022853, 0.023031, 0.023218, 0.023414,
    0.023619, 0.023831, 0.024052, 0.024280, 0.024516, 0.024759,
    0.025009, 0.025265, 0.025527, 0.025795, 0.026069, 0.026347,
    0.026631, 0.026919, 0.027212, 0.027509, 0.027811, 0.028116,
    0.028426, 0.028739, 0.029058, 0.029381, 0.029710, 0.030045,
    0.030387, 0.030737, 0.031097, 0.031469, 0.031855, 0.032257,
    0.032677};
    //double sum=0;
    //for(int ii=0;ii<180;ii++){
    //   sum+=phase_init[ii+1]*1./180*PI;
    //}
    //printf("sum=%lf\n",sum);

   double xxx;
   double yrange[2];
   xxx=Laser::prandom->Uniform(0,1);
   if(anglerange[1]<=anglerange[0]){
      //theta=gMieScatAngle->Eval(xxx);
      theta=xxx*PI;
      //weight*=1;
      int nn=(int)(xxx*180);
      weight*=(phase_init[nn]*(nn+1-xxx*180)+phase_init[nn+1]*(xxx*180-nn))/0.393957;
      return true;
   }
   else{
      yrange[0]=anglerange[0]/PI;
      yrange[1]=anglerange[1]/PI;
      theta=ProbTransform(xxx,yrange,weight,yrange[0]>0); //from 0 to 1
      if(theta<0) return false;
      //theta=gMieScatAngle->Eval(theta);
      theta*=PI;
      //weight*=1;
      int nn=(int)(xxx*180);
      weight*=(phase_init[nn]*(nn+1-xxx*180)+phase_init[nn+1]*(xxx*180-nn))/0.393957;
      return true;
   }
}
bool Atmosphere::MieScatterAnglePhi(double wavelength, double &phi,double anglerange[2],double &weight){
   if(!Laser::prandom) return false;
   //if(!gMieScatAngle){
   //   printf("Atmosphere::MieScatterAngle: No Mie Scatter Angle Calculated\n");
   //   return false;
   //}

   double xxx;
   double yrange[2];
   xxx=Laser::prandom->Uniform(0,1);
   if(anglerange[1]<=anglerange[0]){
      phi=xxx*2*PI-PI;
      return true;
   }
   else{
      if(anglerange[1]<-PI){
         anglerange[0]+=2*PI;
         anglerange[1]+=2*PI;
      }
      if(anglerange[0]>PI){
         anglerange[0]-=2*PI;
         anglerange[1]-=2*PI;
      }
      yrange[0]=anglerange[0]/PI/2+0.5;
      yrange[1]=anglerange[1]/PI/2+0.5;

      phi=ProbTransform(xxx,yrange,weight,yrange[0]>0); //from 0 to 1
      if(phi<0) return false;
      phi=(phi-0.5)*2*PI;
      return true;
   }
}

double Atmosphere::ZDependence(double z,int type){
   return 1;
}
double Atmosphere::DeltaZ(double z){
   return 1.0e20;
}
double Atmosphere::FreeIntgLength(double lengthrange[2],double &weight){
   if(!Laser::prandom) return 0;
   if((aod_air+aod_aerosol)<0) return 0;
   else if(aod_air+aod_aerosol==0) return 1.0e20;
   double res;
   double xxx=Laser::prandom->Uniform(0,1.);
   double yyy=-1;
   double yrange[2];
   yrange[0]=1-exp(-lengthrange[0]);
   yrange[1]=1-exp(-lengthrange[1]);
   if((yrange[1]<=0||yrange[0]<0)||yrange[0]>=yrange[1]) res=log(1/(1-xxx));
   else{
      yyy=ProbTransform(xxx,yrange,weight,yrange[0]>0);
      res=log(1/(1-yyy));
   }
   //printf("Atmosphere::FreeIntgLength: random number %le generated, range={%le,%le} yyy=%le res=%le\n",xxx,lengthrange[0],lengthrange[1],yyy,res);
   return res;
}
double Atmosphere::FreePathLength(double z0,double dir0[3],double lengthrange[2],double &weight){
   double norm=sqrt(pow(dir0[0],2)+pow(dir0[1],2)+pow(dir0[2],2));
   double integ=0;
   double length=0;

   //first estimate the integral length range from the path length range
   double yrange[2]={-1,-1};
   if(lengthrange[1]>lengthrange[0]&&lengthrange[0]>=0){
      if(dir0[2]==0){
         yrange[0]=lengthrange[0]*(aod_air+aod_aerosol)*ZDependence(z0);
         yrange[1]=lengthrange[1]*(aod_air+aod_aerosol)*ZDependence(z0);
      }
      else{
         double z00=z0;
         bool findlowrange=false;
         while(length<lengthrange[1]){
            double z1=z00+DeltaZ(z00);
            double dleng=fabs((z1-z00)/dir0[2]*norm);
            double dintg=(aod_air+aod_aerosol)*ZDependence((z00+z1)/2)*dleng;
            if(!findlowrange){
               if(length+dleng>=lengthrange[0]){
                  yrange[0]=integ+(lengthrange[0]-length)/dleng*dintg;
                  findlowrange=true;
               }
            }
            if(length+dleng>=lengthrange[1]){
               yrange[1]=integ+(lengthrange[1]-length)/dleng*dintg;
               break;
            }
            integ+=dintg;
            length+=dleng;
            z00=z1;
         }
      }
   }
   //if(Laser::jdebug>2) printf("Atmosphere::FreePathLength: lengthrange={%le,%le} integral_range={%le,%le}\n",lengthrange[0],lengthrange[1],yrange[0],yrange[1]);

   double intglength=FreeIntgLength(yrange,weight);
   if(intglength>0.9e20) return 1.0e20;
   if(dir0[2]==0){
      return intglength/((aod_air+aod_aerosol)*ZDependence(z0));
   }

   integ=0;
   length=0;
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
      integ+=dintg;
      length+=fabs((z1-z0)/dir0[2]*norm);
      z0=z1;
      nstep++;
   }
   //if(Laser::jdebug>2) printf("Laser::FreePathLength: integ=%le intglength=%le length=%le\n",integ,intglength,length);
   return -1;
}

int Atmosphere::IsScattering(double z0){
   if(!Laser::prandom) return 0;
   double sum=aod_air*ZDependence(z0,1)+aod_aerosol*ZDependence(z0,2);
   double sum_ext=aod_air*ZDependence(z0,1)-scat_air*ZDependence(z0,3)+aod_aerosol*ZDependence(z0,2)-scat_aerosol*ZDependence(z0,4);
   if(sum_ext<0) sum_ext=0;
   double sum_scat=scat_air*ZDependence(z0,3)+scat_aerosol*ZDependence(z0,4);
   sum=sum_ext+sum_scat;
   if(sum<=0) return 100;
   double xx=Laser::prandom->Uniform();
   if(xx<sum_ext/sum) return 0;  //absorbed
   else if(xx<sum_ext/sum+scat_air*ZDependence(z0,3)/sum) return 1; //Rayleigh Scattering
   else return 2;  //Mie Scattering
}


/*\
The class for laser photon propagation
*/

int Laser::jdebug=0;
int Laser::Doigen=-1;
bool Laser::DoPlot=false;
bool Laser::Doextin=true;
float Laser::IniRange[4][2]={{1.0e9,-1.0e9},{1.0e9,-1.0e9},{1.0e9,-1.0e9},{1.0e9,-1.0e9}};
TRandom3* Laser::prandom = 0;
double Laser::TelSimDist=400.; //in cm
double Laser::TelSimAngl=12.; //in degree
double Laser::scale=1.0;
double Laser::unittime=1600.; //in ns
double Laser::intensity = 0;//mj
double Laser::intensity_err = 0;
double Laser::wavelength0 = 0;//nm
double Laser::wavelength0_err = 0;
double Laser::frequency = 0;
double Laser::pulsetime = 0;
double Laser::spotrange[2] = {0,0};//0.001;//mm, the range of the initial laser spot
double Laser::divergence[2] = {0,0};//0.0573; //mrad

double Laser::LaserCooErr=100.; //in cm
double Laser::LaserZenErr=0.02;//0.02; //in degree
double Laser::LaserAziErr=0.02;//0.02; //in degree

double Laser::lengthmin=400;
double Laser::lengthmax=2000;
double Laser::lengthmin2=400;
double Laser::lengthmax2=4000;
TH1D* Laser::hdenu;
TH1D* Laser::hprob[NCTMax][1024];
TH1D* Laser::hleng[NCTMax][1024];
void Laser::Init(int seed){
   if(!prandom) prandom = new TRandom3();
   prandom->SetSeed(seed);
   pwfc=0;
   ievent_gen=0;
   count_gen=0;
   for(int ii=0;ii<3;ii++) lasercoo[ii]=0;
   for(int ii=0;ii<2;ii++) laserdir[ii]=0;
   plot=0;
   for(int ii=0;ii<4;ii++) {plotrange[ii][0]=IniRange[ii][0]; plotrange[ii][1]=IniRange[ii][1];}
   hdenu=0;
   for(int itel=0;itel<NCTMax;itel++){
      for(int isipm=0;isipm<NSIPM;isipm++){
         hprob[itel][isipm]=0;
         hleng[itel][isipm]=0;
      }
   }

   Reset();
}
void Laser::Release(){
   if(prandom) delete prandom;
   if(pwfc) delete pwfc;
   if(plot){
      //plot->Delete();
      delete plot;
   }
   if(hdenu) {delete hdenu; hdenu=0;}
   for(int itel=0;itel<NCTMax;itel++){
      for(int isipm=0;isipm<NSIPM;isipm++){
         if(hprob[itel][isipm]) {delete hprob[itel][isipm]; hprob[itel][isipm]=0;}
         if(hleng[itel][isipm]) {delete hleng[itel][isipm]; hleng[itel][isipm]=0;}
      }
   }
}
void Laser::Reset(){
   count_gen=0;
   Time_gen=10000;
   time_gen=0;
   //ievent_gen=0;
   iphoton_gen=0;
   wavelength_gen=0;
   vgwav.clear();
   for(int ii=0;ii<3;ii++){
      coor_gen[ii]=0;
      dir_gen[ii]=0;
      vgcoo[ii].clear();
      vgdir[ii].clear();
   }
   vowei.clear();
   Telindex=-2;
   votim.clear();
   votel.clear();
   for(int ii=0;ii<3;ii++){
      coor_out[ii]=0;
      dir_out[ii]=0;
      vocoo[ii].clear();
      vodir[ii].clear();
   }
   interpoint=0;
   volength.clear();
   volength2.clear();
   vosipm.clear();
   if(pwfc) pwfc->EventInitial();

   if(plot){
      //plot->Delete();
      delete plot;
      plot=new TObjArray();
   }
   for(int ii=0;ii<4;ii++) {plotrange[ii][0]=IniRange[ii][0]; plotrange[ii][1]=IniRange[ii][1];}
   if(hdenu) {delete hdenu; hdenu=0;}
   for(int itel=0;itel<NCTMax;itel++){
      for(int isipm=0;isipm<NSIPM;isipm++){
         if(hprob[itel][isipm]) {delete hprob[itel][isipm]; hprob[itel][isipm]=0;}
         if(hleng[itel][isipm]) {delete hleng[itel][isipm]; hleng[itel][isipm]=0;}
      }
   }
}
void Laser::SetParameters(char* filename){
   WReadConfig read;
   WReadConfig* readconfig=&read;
   readconfig->readparam(filename);
   for(int i=0;i<3;i++) lasercoo[i] = readconfig->GetLaserCoo(i);
   for(int i=0;i<2;i++) laserdir[i] = readconfig->GetLaserDir(i);
   intensity = readconfig->GetLaserIntensity();
   intensity_err = readconfig->GetLaserIntensityErr();
   wavelength0 = readconfig->GetLaserWavelength();
   wavelength0_err = readconfig->GetLaserWavelengthErr();
   frequency = readconfig->GetLaserFrequency();
   pulsetime = readconfig->GetLaserPulsetime();
   for(int i=0;i<2;i++) spotrange[i] = readconfig->GetLaserSpotrange(i);
   for(int i=0;i<2;i++) divergence[i] = readconfig->GetLaserDivergence(i);
   printf("Laser::SetParameters: load parameters from %s laser_coo={%lf,%lf,%lf} laser_dir={%lf,%lf} pulsetime=%le\n",filename,lasercoo[0],lasercoo[1],lasercoo[2],laserdir[0],laserdir[1],pulsetime);

   //lasercoo[0]=340.*100.; //in cm
   //lasercoo[1]=0;
   //lasercoo[2]=300.;
   //laserdir[0]=0.;
   //laserdir[1]=180.;
}

void Laser::cross(double dir1[3],double dir2[3],double *dir3){
   dir3[0]=dir1[1]*dir2[2]-dir1[2]*dir2[1];
   dir3[1]=dir1[2]*dir2[0]-dir1[0]*dir2[2];
   dir3[2]=dir1[0]*dir2[1]-dir1[1]*dir2[0];
}
bool Laser::CartesianFrame(double zero[3],double coor_in[3],double dir_in[3],double dir_in2[3],double *xdir,double *ydir,double *zdir){
   //z axis
   for(int ii=0;ii<3;ii++) zdir[ii]=dir_in[ii];
   double norm=sqrt(zdir[0]*zdir[0]+zdir[1]*zdir[1]+zdir[2]*zdir[2]);
   if(norm<=0) return false;
   for(int ii=0;ii<3;ii++) zdir[ii]/=norm;

   double dir1[3]={coor_in[0]-zero[0],coor_in[1]-zero[1],coor_in[2]-zero[2]};
   norm=sqrt(dir1[0]*dir1[0]+dir1[1]*dir1[1]+dir1[2]*dir1[2]);
   if(norm<=0) return false;

   //temporary x and y axis
   double xdir0[3],ydir0[3];
   cross(zdir,dir1,xdir0);
   norm=sqrt(xdir0[0]*xdir0[0]+xdir0[1]*xdir0[1]+xdir0[2]*xdir0[2]);
   if(norm<=0){
      //double mult=dir1[0]*zdir[0]+dir1[1]*zdir[1]+dir1[2]*zdir[2];
      //if(mult>0) return false;
      cross(zdir,dir_in2,xdir0);
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
   double telcoor[NCTMax][3];
   WFTelescopeArray* pta=WFTelescopeArray::GetHead();
   if(!pta){
      nct=1; telcoor[0][0]=telcoor[0][1]=telcoor[0][2]=0;
   }
   else{
      nct=pta->CTNumber;
      for(int ii=0;ii<nct;ii++){
         WFTelescope* pt=pta->pct[ii];
         telcoor[ii][0]=pt->Telx_;
         telcoor[ii][1]=pt->Tely_;
         telcoor[ii][2]=pt->Telz_;
      }
   }
   whichtel=-2;
   //loop the telescope
   double min=1.0e10;
   for(int ii=0;ii<nct;ii++){
      double coor_min0[3];
      double zero[3]={telcoor[ii][0],telcoor[ii][1],telcoor[ii][2]};
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
   //xx=0; yy=0; return;
   //double rr=prandom->Uniform(0,spotrange*0.1);
   //double phi=prandom->Uniform(0,2*PI);
   //xx=rr*cos(phi);
   //yy=rr*sin(phi);
   double angle=0;
   double xx0=prandom->Uniform(0,spotrange[0]*0.1); //from mm to cm
   double yy0=prandom->Uniform(0,spotrange[1]*0.1); //from mm to cm
   xx=xx0*cos(angle)-yy0*sin(angle);
   yy=xx0*sin(angle)+yy0*cos(angle);
}
void Laser::DirectionDis(double &theta,double &phi){
   //double ran1=prandom->Uniform(cos(divergence*1.0e-3),1);
   //theta=acos(sqrt(ran1));

   //double ran1=prandom->Gaus(0,divergence*1.0e-3);
   //theta=fabs(ran1);
   //phi=prandom->Uniform(0,2*PI);

   double angle=0;
   double dx0=tan(prandom->Gaus(0,divergence[0]*1.0e-3));
   double dy0=tan(prandom->Gaus(0,divergence[1]*1.0e-3));
   double dx=dx0*cos(angle)-dy0*sin(angle);
   double dy=dx0*sin(angle)+dy0*cos(angle);
   double dr=sqrt(1+dx*dx+dy*dy);
   theta=acos(1/dr);
   phi=(dx==0)?(dy>=0?(PI/2):(-PI/2)):atan(dy/dx);
   if(dx<0) phi+=PI;
}
void Laser::GetAveTelPos(double zero[3]){
   WFTelescopeArray* pta=WFTelescopeArray::GetHead();
   if(!pta) return;
   double avecoo[3]={0,0,0};
   int ncount=0;
   for(int ii=0;ii<pta->CTNumber;ii++){
      WFTelescope* pt=pta->pct[ii];
      if(!pt) continue;
      ncount++;
      avecoo[0]+=pt->Telx_;
      avecoo[1]+=pt->Tely_;
      avecoo[2]+=pt->Telz_;
   }
   if(ncount>0){
      for(int ii=0;ii<3;ii++) zero[ii]=avecoo[ii]/ncount;
   }
}
bool Laser::InitialGen(){
   double norm=1;
   //direction
   double theta,phi;
   DirectionDis(theta,phi);
   double xdir[3],ydir[3],zdir[3];
   double zero[3]={0,0,0};
   GetAveTelPos(zero);
   double dir_in[3]={sin(laserdir[0]/180.*PI)*cos(laserdir[1]/180.*PI),sin(laserdir[0]/180.*PI)*sin(laserdir[1]/180.*PI),cos(laserdir[0]/180.*PI)};
   double dir_in2[3]={0,0,1};
   if(!CartesianFrame(zero,lasercoo,dir_in,dir_in2,xdir,ydir,zdir)) return false;

   for(int ii=0;ii<3;ii++){
      dir_gen[ii]=cos(theta)*zdir[ii]+sin(theta)*(cos(phi)*xdir[ii]+sin(phi)*ydir[ii]);
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
   //return wavelength0;
   double res=prandom->Gaus(wavelength0,wavelength0_err);
   return res>0?res:0;
}

long int Laser::EventGen(int &Time,double &time,bool SimPulse){
   if(!prandom) return 0;
   double acctime,time1,time2,ngen0;
   if(!SimPulse){
      acctime=0;
      time1=Time+time*20*1.0e-9;
      time2=Time+time*20*1.0e-9+unittime*1.0e-9;
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
      ngen0=intensity*1.0e-3/(hplank*vlight/(wavelength0*1.0e-7))/pulsetime*1.0e-6*scale;
   }
   else{
      acctime=pulsetime;
      time1=Time+time*20*1.0e-9;
      time2=time1+acctime;
      ngen0=intensity*1.0e-3/(hplank*vlight/(wavelength0*1.0e-7))/pulsetime*scale;
   }

   Reset();
   long int ngen=(long int)(acctime*ngen0);
   Time_gen=Time;
   time_gen=time;
   long int ngentel=0;
   for(long int igen=0;igen<ngen;igen++){
      bool dogen=InitialGen();
      if(!dogen) continue;
      if(Doigen>=0&&Doigen!=igen) continue;
      iphoton_gen=igen;
      double time0=time1+(time2-time1)/ngen*(igen+0.5);
      vgwav.push_back(wavelength_gen);
      for(int ii=0;ii<3;ii++){
         vgcoo[ii].push_back(coor_gen[ii]);
         vgdir[ii].push_back(dir_gen[ii]);
      }
      double weight=1./scale;
      double distance;
      int res=Propagate(distance,weight);
      if((igen%(1000000)==0)&&jdebug>0) printf("Laser::EventGen: %ld of %ld generated (count_gen=%le,weight={%le,%le})\n",igen,ngen,count_gen,weight,weight*scale);
      if(jdebug>3) printf("Laser::EventGen: Propagate igen=%ld res=%d distance=%lf lasercoo={%f,%f,%f} laserdir={%f,%f,%f}\n",igen,res,distance,coor_gen[0],coor_gen[1],coor_gen[2],dir_gen[0],dir_gen[1],dir_gen[2]);
      if(res<0) Telindex=res-15;
      else{  //the telescope index has been calculated in Propagate
         ngentel++;
      }
      vowei.push_back(weight);
      votim.push_back(time0+distance/vlight);
      votel.push_back(Telindex);
      for(int ii=0;ii<3;ii++){
         vocoo[ii].push_back(res>=0?coor_out[ii]:0);
         vodir[ii].push_back(res>=0?dir_out[ii]:0);
      }
      count_gen+=weight;
   }
   if(jdebug>0) printf("Laser::EventGen: ngen0=%le acctime=%le ngen=%ld ngentel=%ld scale=%le\n",ngen0,acctime,ngen,ngentel,scale);
   bool dosim=DoWFCTASim();
   ievent_gen++;

   if(!SimPulse){
      time+=unittime/20;
      if((time*20*1e-9)>=1.0){
         Time++;
         time-=1.0/(20*1.0e-9);
      }
   }
   else Time++;

   return ngentel;
}

int Laser::FindLengthRange(double zero[3],double cooout[3],double dirout[3],double dirin[3],double lengthrange[2]){
   double MAXNUM=1.0e20;
   double dirlas[3]={cooout[0]-zero[0],cooout[1]-zero[1],cooout[2]-zero[2]};
   double dirin2[3]={dirin[0],dirin[1],dirin[2]};
   double dirout2[3]={dirout[0],dirout[1],dirout[2]};
   double norm_dir_las=sqrt(pow(dirlas[0],2)+pow(dirlas[1],2)+pow(dirlas[2],2));
   double norm_dir_in =sqrt(pow(dirin2[0],2)+pow(dirin2[1],2)+pow(dirin2[2],2));
   double norm_dir_out=sqrt(pow(dirout2[0],2)+pow(dirout2[1],2)+pow(dirout2[2],2));
   if(norm_dir_in<=0||norm_dir_las<=0||norm_dir_out<=0) return -1;
   for(int ii=0;ii<3;ii++){
      dirout2[ii]/=norm_dir_out;
      dirin2[ii]/=norm_dir_in;
      dirlas[ii]/=norm_dir_las;
   }
   double costhetaref=0;
   for(int ii=0;ii<3;ii++){
      costhetaref+=dirout2[ii]*(-dirlas[ii]);
   }
   if(fabs(costhetaref)>cos(asin(TelSimDist/norm_dir_las))){
      if(jdebug>5) printf("Laser::FindLengthRange: in reg1, angleref=%lf angle0=%lf\n",acos(costhetaref),asin(TelSimDist/norm_dir_las));
      if(costhetaref>0){ //out direction point to nearby of telescope
         lengthrange[0]=0;
         lengthrange[1]=MAXNUM;
         return 1;
      }
      else{ //point to the opposite direction
         double costheta1=0;
         for(int ii=0;ii<3;ii++){
            costheta1+=dirin2[ii]*(-dirlas[ii]);
         }
         if(acos(costheta1)<TelSimAngl/180.*PI){
            lengthrange[0]=0;
            lengthrange[1]=MAXNUM;
            return 2;
         }
         else return -1;
      }
   }
   else{
      if(jdebug>5) printf("Laser::FindLengthRange: in reg2, angleref=%lf angle0=%lf\n",acos(costhetaref),asin(TelSimDist/norm_dir_las));
      double dirzz[3],diryy[3];
      cross(dirlas,dirout2,dirzz);
      cross(dirlas,dirzz,diryy);
      double norm_diryy=sqrt(pow(diryy[0],2)+pow(diryy[1],2)+pow(diryy[2],2));
      for(int ii=0;ii<3;ii++) diryy[ii]/=norm_diryy;
      double costhetainxx=0,costhetainyy=0;
      double costhetaoutxx=0,costhetaoutyy=0;
      for(int ii=0;ii<3;ii++){
         costhetainxx+=dirin2[ii]*(-dirlas[ii]);
         costhetainyy+=dirin2[ii]*(diryy[ii]);
         costhetaoutxx+=dirout2[ii]*(-dirlas[ii]);
         costhetaoutyy+=dirout2[ii]*(diryy[ii]);
      }
      if(jdebug>6) printf("Laser::FindLengthRange: dirout={%le,%le} dirin={%le,%le,%le}\n",costhetaoutxx,costhetaoutyy,costhetainxx,costhetainyy,sqrt(1-costhetainxx*costhetainxx-costhetainyy*costhetainyy));
      //find the minimum and maximum angle.
      // theta_min<PI/2 and theta_max>PI/2
      double AAA=costhetaoutxx*costhetainyy-costhetainxx*costhetaoutyy;
      double coeffi=AAA==0?(0):(costhetainyy*norm_dir_las/AAA);
      double vector=AAA==0?(0):(-costhetaoutyy*norm_dir_las/AAA);

      double aa=pow(costhetainxx*costhetaoutxx+costhetainyy*costhetaoutyy,2)-pow(cos(TelSimAngl/180.*PI),2);
      double bb=-2*norm_dir_las*(costhetainxx*(costhetainxx*costhetaoutxx+costhetainyy*costhetaoutyy)-costhetaoutxx*pow(cos(TelSimAngl/180.*PI),2));
      double cc=norm_dir_las*norm_dir_las*(costhetainxx*costhetainxx-pow(cos(TelSimAngl/180.*PI),2));
      if(AAA==0){
         double theta_min,theta_max;
         double alpha_min,alpha_max;
         if(-costhetaoutyy/costhetainyy>0){ //reach minimum
            theta_min=acos(sqrt(pow(costhetainxx,2)+pow(costhetainyy,2)));
            alpha_min=MAXNUM;
            theta_max=acos(costhetainxx);
            alpha_max=0;
            if(TelSimAngl/180.*PI>theta_max){
               lengthrange[0]=0;
               lengthrange[1]=MAXNUM;
               return 3;
            }
            else if(TelSimAngl/180.*PI>theta_min){
               lengthrange[0]=(-bb+(aa>0?1:-1)*sqrt(bb*bb-4*aa*cc))/(2*aa);
               lengthrange[1]=TMath::Max(10*lengthrange[0],MAXNUM);
               return 4;
            }
            else return -2;
         }
         else{ //reach maximum
            theta_min=acos(costhetainxx);
            alpha_min=0;
            theta_max=acos(-sqrt(pow(costhetainxx,2)+pow(costhetainyy,2)));
            alpha_max=MAXNUM;
            if(TelSimAngl/180.*PI>theta_max){
               lengthrange[0]=0;
               lengthrange[1]=MAXNUM;
               return 5;
            }
            else if(TelSimAngl/180.*PI>theta_min){
               lengthrange[0]=0;
               lengthrange[1]=(-bb+(aa>0?1:-1)*sqrt(bb*bb-4*aa*cc))/(2*aa);
               return 6;
            }
            else return -3;
         }
      }
      else if(costhetainyy/AAA>0){ //could reach maximum or minimum in the middle
         double theta_min,theta_max;
         double alpha_min,alpha_max;
         if(-costhetaoutyy/AAA<0){ //reach maximum
            double theta_min2;
            double alpha_min2;
            theta_min=acos(costhetainxx);
            alpha_min=0;
            theta_max=acos(-sqrt(pow(costhetainxx,2)+pow(costhetainyy,2)));
            alpha_max=costhetainyy/AAA*norm_dir_las;
            theta_min2=acos(-costhetainxx*costhetaoutxx-costhetainyy*costhetaoutyy);
            alpha_min2=MAXNUM;
            if(TelSimAngl/180.*PI>theta_max){
               lengthrange[0]=0;
               lengthrange[1]=MAXNUM;
               return 7;
            }
            else if(TelSimAngl/180.*PI>TMath::Max(theta_min,theta_min2)){
               lengthrange[0]=(-bb+(aa>0?-1:1)*sqrt(bb*bb-4*aa*cc))/(2*aa);
               lengthrange[1]=(-bb+(aa>0?1:-1)*sqrt(bb*bb-4*aa*cc))/(2*aa);
               return 8;
            }
            else if(TelSimAngl/180.*PI>TMath::Min(theta_min,theta_min2)){
               if(theta_min<theta_min2){
                  lengthrange[0]=0;
                  lengthrange[1]=(-bb+(aa>0?-1:1)*sqrt(bb*bb-4*aa*cc))/(2*aa);
               }
               else{
                  lengthrange[0]=(-bb+(aa>0?1:-1)*sqrt(bb*bb-4*aa*cc))/(2*aa);
                  lengthrange[1]=TMath::Max(10*lengthrange[0],MAXNUM);
               }
               return 9;
            }
            else return -4;
         }
         else{ //reach minimum
            double theta_max2;
            double alpha_max2;
            theta_max=acos(costhetainxx);
            alpha_max=0;
            theta_min=acos(sqrt(pow(costhetainxx,2)+pow(costhetainyy,2)));
            alpha_min=costhetainyy/AAA*norm_dir_las;
            theta_max2=acos(-costhetainxx*costhetaoutxx-costhetainyy*costhetaoutyy);
            alpha_max2=MAXNUM;
            if(jdebug>10) printf("Laser::FindLengthRange: reach minimum, theta_max=%le alpha_max=%le theta_min=%le alpha_min=%le theta_max2=%le alpha_max2=%le\n",theta_max/PI*180,alpha_max,theta_min/PI*180,alpha_min,theta_max2/PI*180,alpha_max2);
            if(TelSimAngl/180.*PI>TMath::Max(theta_max,theta_max2)){
               lengthrange[0]=0;
               lengthrange[1]=MAXNUM;
               return 10;
            }
            else if(TelSimAngl/180.*PI>TMath::Min(theta_max,theta_max2)){
               if(theta_max<theta_max2){
                  lengthrange[0]=0;
                  lengthrange[1]=(-bb+(aa>0?1:-1)*sqrt(bb*bb-4*aa*cc))/(2*aa);
               }
               else{
                  lengthrange[0]=(-bb+(aa>0?-1:1)*sqrt(bb*bb-4*aa*cc))/(2*aa);
                  lengthrange[1]=TMath::Max(10*lengthrange[0],MAXNUM);
               }
               return 11;
            }
            else if(TelSimAngl/180.*PI>theta_min){
               lengthrange[0]=(-bb-(aa>0?1:-1)*sqrt(bb*bb-4*aa*cc))/(2*aa);
               lengthrange[1]=(-bb+(aa>0?1:-1)*sqrt(bb*bb-4*aa*cc))/(2*aa);
               return 12;
            }
            else return -5;
         }
      }
      else{ //couldn't reach maximum or minimum in the middle
         double theta_min,theta_max;
         double alpha_min,alpha_max;
         if(-costhetaoutyy/AAA<0){ //reach maximum outside the region
            theta_max=acos(costhetainxx);
            alpha_max=0;
            theta_min=acos(-costhetainxx*costhetaoutxx-costhetainyy*costhetaoutyy);
            alpha_min=MAXNUM;
            if(TelSimAngl/180.*PI>theta_max){
               lengthrange[0]=0;
               lengthrange[1]=MAXNUM;
               return 13;
            }
            else if(TelSimAngl/180.*PI>theta_min){
               lengthrange[0]=(-bb+(aa>0?1:-1)*sqrt(bb*bb-4*aa*cc))/(2*aa);
               lengthrange[1]=TMath::Max(10*lengthrange[0],MAXNUM);
               return 14;
            }
            else return -6;
         }
         else{ //reach minimum
            theta_min=acos(costhetainxx);
            alpha_min=0;
            theta_max=acos(-costhetainxx*costhetaoutxx-costhetainyy*costhetaoutyy);
            alpha_max=MAXNUM;
            if(TelSimAngl/180.*PI>theta_max){
               lengthrange[0]=0;
               lengthrange[1]=MAXNUM;
               return 15;
            }
            else if(TelSimAngl/180.*PI>theta_min){
               lengthrange[0]=0;
               lengthrange[1]=(-bb+(aa>0?1:-1)*sqrt(bb*bb-4*aa*cc))/(2*aa);
               return 16;
            }
            else return -7;
         }
      }
   }
}
int Laser::FindThetaRange(double zero[3],double cooout[3],double dirout[3],double dirin[3],double freelength,double thetarange[2]){
   double dirlas[3]={cooout[0]-zero[0],cooout[1]-zero[1],cooout[2]-zero[2]};
   double dirin2[3]={dirin[0],dirin[1],dirin[2]};
   double dirout2[3]={dirout[0],dirout[1],dirout[2]};
   double norm_dir_las=sqrt(pow(dirlas[0],2)+pow(dirlas[1],2)+pow(dirlas[2],2));
   double norm_dir_in =sqrt(pow(dirin2[0],2)+pow(dirin2[1],2)+pow(dirin2[2],2));
   double norm_dir_out=sqrt(pow(dirout2[0],2)+pow(dirout2[1],2)+pow(dirout2[2],2));
   if(norm_dir_in<=0||norm_dir_las<=0||norm_dir_out<=0) return false;
   for(int ii=0;ii<3;ii++){
      dirout2[ii]/=norm_dir_out;
      dirin2[ii]/=norm_dir_in;
      dirlas[ii]/=norm_dir_las;
   }
   double costhetaref=0;
   for(int ii=0;ii<3;ii++){
      costhetaref+=dirout2[ii]*(-dirlas[ii]);
   }
   if(costhetaref>cos(asin(TelSimDist/norm_dir_las))){
      if(jdebug>7) printf("Laser::FindThetaRange: in reg1, angleref=%lf angle0=%lf\n",acos(costhetaref),asin(TelSimDist/norm_dir_las));
      thetarange[0]=0;
      thetarange[1]=PI/2;
      return 1;
   }
   else{
      if(jdebug>7) printf("Laser::FindThetaRange: in reg2, angleref=%lf angle0=%lf\n",acos(costhetaref),asin(TelSimDist/norm_dir_las));
      double costheta1=0,costheta2=0,costheta3=0;
      for(int ii=0;ii<3;ii++){
         costheta1+=dirin2[ii]*(-dirlas[ii]);
         costheta2+=dirin2[ii]*dirout2[ii];
         costheta3+=dirout2[ii]*(-dirlas[ii]);
      }
      double theta_min,theta_max;
      double length0=sqrt(pow(norm_dir_las,2)+pow(freelength,2)-2*norm_dir_las*freelength*costheta3);
      double theta111=acos((-freelength+norm_dir_las*costheta3)/length0);
      double theta222=acos((-freelength*costheta2+norm_dir_las*costheta1)/length0);
      if(jdebug>8) printf("Laser::FindThetaRange: theta111=%lf theta222=%lf\n",theta111/PI*180,theta222/PI*180);
      if(theta222>TelSimAngl/180.*PI+asin(TelSimDist/length0)){
         return -1;
      }
      else{
         thetarange[0]=TMath::Max((double)0.,(double)theta111-asin(TelSimDist/length0));
         thetarange[1]=TMath::Min((double)PI,(double)theta111+asin(TelSimDist/length0));
         return 2;
      }
   }
}
int Laser::FindPhiRange(double zero[3],double cooout[3],double dirout[3],double dirin[3],double freelength,double theta_scat,double phirange[2]){
   double dirlas[3]={cooout[0]-zero[0],cooout[1]-zero[1],cooout[2]-zero[2]};
   double dirin2[3]={dirin[0],dirin[1],dirin[2]};
   double dirout2[3]={dirout[0],dirout[1],dirout[2]};
   double norm_dir_las=sqrt(pow(dirlas[0],2)+pow(dirlas[1],2)+pow(dirlas[2],2));
   double norm_dir_in =sqrt(pow(dirin2[0],2)+pow(dirin2[1],2)+pow(dirin2[2],2));
   double norm_dir_out=sqrt(pow(dirout2[0],2)+pow(dirout2[1],2)+pow(dirout2[2],2));
   if(norm_dir_in<=0||norm_dir_las<=0||norm_dir_out<=0) return false;
   for(int ii=0;ii<3;ii++){
      dirout2[ii]/=norm_dir_out;
      dirin2[ii]/=norm_dir_in;
      dirlas[ii]/=norm_dir_las;
   }

   double cooout2[3]={cooout[0]+freelength*dirout2[0],cooout[1]+freelength*dirout2[1],cooout[2]+freelength*dirout2[2]};
   double dirout3[3]={0,0,1};
   double xdir[3],ydir[3],zdir[3];
   CartesianFrame(zero,cooout2,dirout2,dirout3,xdir,ydir,zdir);

   double dir111[3],dir222[3];
   double norm_dir111=0;
   for(int ii=0;ii<3;ii++){
      dir111[ii]=zero[ii]-cooout2[ii];
      norm_dir111+=pow(dir111[ii],2);
      dir222[ii]=dirin2[ii];
   }
   if(norm_dir111<=0){
      phirange[0]=0;
      phirange[1]=2*PI;
      return true;
   }
   norm_dir111=sqrt(norm_dir111);
   for(int ii=0;ii<3;ii++) dir111[ii]/=norm_dir111;

   double dir111_theta,dir111_phi,dir222_theta,dir222_phi;
   dir111_theta=acos(dir111[0]*zdir[0]+dir111[1]*zdir[1]+dir111[2]*zdir[2]);
   double xcomp=dir111[0]*xdir[0]+dir111[1]*xdir[1]+dir111[2]*xdir[2];
   double ycomp=dir111[0]*ydir[0]+dir111[1]*ydir[1]+dir111[2]*ydir[2];
   dir111_phi=sin(dir111_theta)<=0?0:(xcomp!=0?atan(ycomp/xcomp):(ycomp>0?(PI/2):(-PI/2)));
   if(xcomp<0) dir111_phi+=PI;
   if(dir111_phi>PI) dir111_phi-=2*PI;
   dir222_theta=acos(dir222[0]*zdir[0]+dir222[1]*zdir[1]+dir222[2]*zdir[2]);
   xcomp=dir222[0]*xdir[0]+dir222[1]*xdir[1]+dir222[2]*xdir[2];
   ycomp=dir222[0]*ydir[0]+dir222[1]*ydir[1]+dir222[2]*ydir[2];
   dir222_phi=sin(dir222_theta)<=0?0:(xcomp!=0?atan(ycomp/xcomp):(ycomp>0?(PI/2):(-PI/2)));
   if(xcomp<0) dir222_phi+=PI;
   if(dir222_phi>PI) dir222_phi-=2*PI;

   double ref1=cos(asin(TelSimDist/norm_dir111));
   double ref2=cos(TelSimAngl/180.*PI);
   double value1=cos(theta_scat)*cos(dir111_theta);
   double value2=cos(theta_scat)*cos(dir222_theta);
   if(theta_scat<=0){
      if(value1>ref1&&value2>ref2){
         phirange[0]=0;
         phirange[1]=2*PI;
         return true;
      }
      else return false;
   }
   else{
      if(sin(dir111_theta)>0&&sin(dir222_theta)>0){
         double value11=(ref1-value1)/sin(theta_scat)/sin(dir111_theta);
         double value22=(ref2-value2)/sin(theta_scat)/sin(dir222_theta);
         double range1[2]={-100,-100},range2[2]={-100,-100};
         if(value11<0){
            range1[0]=-PI;
            range1[1]=PI;
         }
         else if(value11<=1){
            range1[0]=-acos(value11)+dir111_phi;
            range1[1]=acos(value11)+dir111_phi;
         }
         if(range1[0]>PI){
            range1[0]-=2*PI;
            range1[1]-=2*PI;
         }
         if(range1[1]<-PI){
            range1[0]+=2*PI;
            range1[1]+=2*PI;
         }
         if(value22<0){
            range2[0]=-PI;
            range2[1]=PI;
         }
         else if(value22<=1){
            range2[0]=-acos(value22)+dir222_phi;
            range2[1]=acos(value22)+dir222_phi;
         }
         if(range2[0]>PI){
            range2[0]-=2*PI;
            range2[1]-=2*PI;
         }
         if(range2[1]<-PI){
            range2[0]+=2*PI;
            range2[1]+=2*PI;
         }

         phirange[0]=TMath::Max(range1[0],range2[0]);
         phirange[1]=TMath::Min(range1[1],range2[1]);
         if(phirange[1]>phirange[0]) return true;
         else return false;
      }
      else if(sin(dir111_theta)>0&&sin(dir222_theta)<=0){
         double value11=(ref1-value1)/sin(theta_scat)/sin(dir111_theta);
         double range1[2]={-100,-100};
         if(value11<0){
            range1[0]=-PI;
            range1[1]=PI;
         }
         else if(value11<=1){
            range1[0]=-acos(value11)+dir111_phi;
            range1[1]=acos(value11)+dir111_phi;
         }
         if(range1[0]>PI){
            range1[0]-=2*PI;
            range1[1]-=2*PI;
         }
         if(range1[1]<-PI){
            range1[0]+=2*PI;
            range1[1]+=2*PI;
         }
         phirange[0]=range1[0];
         phirange[1]=range1[1];
         if(phirange[1]>phirange[0]) return true;
         else return false;
      }
      else if(sin(dir111_theta)<=0&&sin(dir222_theta)>0){
         double value22=(ref2-value2)/sin(theta_scat)/sin(dir222_theta);
         double range2[2]={-100,-100};
         if(value22<0){
            range2[0]=-PI;
            range2[1]=PI;
         }
         else if(value22<=1){
            range2[0]=-acos(value22)+dir222_phi;
            range2[1]=acos(value22)+dir222_phi;
         }
         if(range2[0]>PI){
            range2[0]-=2*PI;
            range2[1]-=2*PI;
         }
         if(range2[1]<-PI){
            range2[0]+=2*PI;
            range2[1]+=2*PI;
         }
         phirange[0]=range2[0];
         phirange[1]=range2[1];
         if(phirange[1]>phirange[0]) return true;
         else return false;
      }
      else{
         if(value1>ref1&&value2>ref2){
            phirange[0]=0;
            phirange[1]=2*PI;
            return true;
         }
         else return false;
      }
   }
}

int Laser::Propagate(double &distance,double &weight){
   double weight0=weight;
   double coor_min[3];
   int whichtel=-1;
   bool decrease;
   double lengthrange[2]={-1,-1};
   double thetarange[2]={-1,-1};
   double phirange[2]={-100,-100};

   WFTelescopeArray* pta=WFTelescopeArray::GetHead();
   double lengthrange_tel[NCTMax][2];
   int telindex[NCTMax];
   for(int itel=0;itel<WFTelescopeArray::CTNumber;itel++) telindex[itel]=-1;
   int ntel=0;
   for(int itel=0;itel<WFTelescopeArray::CTNumber;itel++){
      lengthrange_tel[itel][0]=lengthrange_tel[itel][1]=-1;
      WFTelescope* pt=(pta)?pta->pct[itel]:0;
      if(!pt) continue;
      double dir_tel[3]={-sin(pt->TelZ_)*cos(pt->TelA_),-sin(pt->TelZ_)*sin(pt->TelA_),-cos(pt->TelZ_)}; //pointing direction of the telescope
      double zero[3]={pt->Telx_,pt->Tely_,pt->Telz_};
      double range[2];
      int findres=FindLengthRange(zero,coor_gen,dir_gen,dir_tel,range);
      if(jdebug>5) printf("Laser::Propagate: find telescope and length range, ntel=%d itel=%d lengthrange={%le,%le} return=%d\n",ntel,itel,range[0],range[1],findres);
      if(findres>0){
         lengthrange_tel[ntel][0]=range[0];
         lengthrange_tel[ntel][1]=range[1];
         telindex[ntel]=itel;
         ntel++;
      }
   }

   if(ntel>0){
      double ran0=prandom->Uniform(0,1.);
      for(int ii=0;ii<ntel;ii++){
         double low=1./ntel*ii;
         double hig=1./ntel*(ii+1);
         if(ran0>=low&&ran0<hig){
            whichtel=telindex[ii];
            lengthrange[0]=lengthrange_tel[ii][0];
            lengthrange[1]=lengthrange_tel[ii][1];
            break;
         }
         else continue;
      }
   }
   WFTelescope* pt=0;
   if(whichtel>=0&&pta) pt=pta->pct[whichtel];

   double dist;
   double zero[3]={0,0,0};
   double dir_tel[3]={0,0,-1};
   if(!pt){
      dist=mindist(zero,coor_gen,dir_gen,coor_min,decrease);
   }
   else{
      zero[0]=pt->Telx_;
      zero[1]=pt->Tely_;
      zero[2]=pt->Telz_;
      dir_tel[0]=-sin(pt->TelZ_)*cos(pt->TelA_);
      dir_tel[1]=-sin(pt->TelZ_)*sin(pt->TelA_);
      dir_tel[2]=-cos(pt->TelZ_);
      dist=mindist(zero,coor_gen,dir_gen,coor_min,decrease);
   }
   distance=sqrt(pow(coor_gen[0]-coor_min[0],2)+pow(coor_gen[1]-coor_min[1],2)+pow(coor_gen[2]-coor_min[2],2));

   //generate free path length, theta angle, and phi angle
   double freelength=Atmosphere::FreePathLength(coor_gen[2],dir_gen,lengthrange,weight);
   interpoint=freelength;
   if(jdebug>4||(!isfinite(weight))) printf("Laser::Propagate: freelength=%le lengthrange={%le,%le} weight={%le,%le}\n",freelength,lengthrange[0],lengthrange[1],weight0,weight);
   if(freelength>lengthrange[0]&&freelength<lengthrange[1]&&whichtel>=0){
      //generate theta angle range
      int findres=FindThetaRange(zero,coor_gen,dir_gen,dir_tel,freelength,thetarange);
      if(jdebug>6||(!isfinite(weight))) printf("Laser::Propagate: find theta range, itel=%d freelength=%le thetarange={%le,%le} return=%d\n",whichtel,freelength,thetarange[0]/PI*180,thetarange[1]/PI*180,findres);
      if(!(findres>0)) {thetarange[0]=thetarange[1]=-1;}
   }
   else{
      thetarange[0]=thetarange[1]=-1;
   }
   double norm_dir_gen=sqrt(pow(dir_gen[0],2)+pow(dir_gen[1],2)+pow(dir_gen[2],2));
   double znew=coor_gen[2]+dir_gen[2]/norm_dir_gen*freelength;
   int scatter=Atmosphere::IsScattering(znew);

   if(decrease&&dist<TelSimDist){ //laser in the field view of the telescope
      if(jdebug>4) printf("Laser::Propagate: laser in the field of view of telescope(%lf), distance=%lf(free length=%lf) decrease=%d coo={%lf,%lf,%lf}\n",dist,distance,freelength,decrease,coor_min[0],coor_min[1],coor_min[2]);
      int returntype;
      if((distance<freelength&&Doextin)||(!Doextin)){
         Telindex=whichtel;
         for(int ii=0;ii<3;ii++){
            coor_out[ii]=coor_min[ii];
            dir_out[ii]=dir_gen[ii];
         }
         returntype=0; //pass through nearby of the telescope without any interaction
      }
      else{ //will be absorbed or scattered, ignore those events
         returntype=-3;
      }
      if(DoPlot){
         double coor_last[4],coor_first[4];
         if(IniRange[3][1]>IniRange[3][0]){
            double length_prop1=TMath::Max((double)0.,(double)IniRange[3][0])*vlight;
            double length_prop2=TMath::Max((double)0.,(double)IniRange[3][1])*vlight;
            double maxlength_prop=(decrease&&dist<TelSimDist)?TMath::Min(distance,freelength):freelength;
            for(int ii=0;ii<3;ii++) coor_first[ii]=coor_gen[ii]+dir_gen[ii]*TMath::Min(length_prop1,maxlength_prop);
            for(int ii=0;ii<3;ii++) coor_last[ii]=coor_gen[ii]+dir_gen[ii]*TMath::Min(length_prop2,maxlength_prop);
            coor_first[3]=TMath::Max((double)0.,(double)IniRange[3][0]);
            coor_last[3]=TMath::Max((double)0.,(double)IniRange[3][1]);
         }
         else{
            for(int ii=0;ii<3;ii++) coor_first[ii]=coor_gen[ii];
            double maxlength_prop=(decrease&&dist<TelSimDist)?TMath::Min(distance,freelength):freelength;
            for(int ii=0;ii<3;ii++) coor_last[ii]=coor_gen[ii]+dir_gen[ii]*maxlength_prop;
            coor_first[3]=0;
            coor_last[3]=maxlength_prop/vlight;
         }
         Add(coor_first,coor_last,returntype,weight);
      }
      return returntype;
   }
   else{
      if(scatter<=0||scatter>2){ //absorbed or no interaction
         if(jdebug>4) printf("Laser::Propagate: laser far away from telescope(%lf), absorbed or no interaction, distance=%lf(free length=%lf) decrease=%d coo={%lf,%lf,%lf} weight=%le\n",dist,distance,freelength,decrease,coor_min[0],coor_min[1],coor_min[2],weight);
         int returntype=-4-(scatter<=0?-scatter:1);
         if(DoPlot){
            double coor_last[4],coor_first[4];
            if(IniRange[3][1]>IniRange[3][0]){
               double length_prop1=TMath::Max((double)0,(double)IniRange[3][0])*vlight;
               double length_prop2=TMath::Max((double)0,(double)IniRange[3][1])*vlight;
               double maxlength_prop=freelength;
               for(int ii=0;ii<3;ii++) coor_first[ii]=coor_gen[ii]+dir_gen[ii]*TMath::Min(length_prop1,maxlength_prop);
               for(int ii=0;ii<3;ii++) coor_last[ii]=coor_gen[ii]+dir_gen[ii]*TMath::Min(length_prop2,maxlength_prop);
               coor_first[3]=TMath::Max((double)0.,(double)IniRange[3][0]);
               coor_last[3]=TMath::Max((double)0.,(double)IniRange[3][1]);
            }
            else{
               for(int ii=0;ii<3;ii++) coor_first[ii]=coor_gen[ii];
               double maxlength_prop=freelength;
               for(int ii=0;ii<3;ii++) coor_last[ii]=coor_gen[ii]+dir_gen[ii]*maxlength_prop;
               coor_first[3]=0;
               coor_last[3]=freelength/vlight;
            }
            Add(coor_first,coor_last,returntype,weight);
         }
         return returntype;
      }
      else{ //scattering
         double coor_scat[3];
         double dir_scat[3];

         for(int ii=0;ii<3;ii++) coor_scat[ii]=coor_gen[ii]+dir_gen[ii]*freelength;

         double theta,phi;
         if(scatter==1){	//Rayleigh scattering
            weight0=weight;
            Atmosphere::RayScatterAngleTheta(wavelength_gen,theta,thetarange,weight);
            if(jdebug>6||(!isfinite(weight))) printf("Laser::Propagate: generated Ray theta=%lf thetarange={%lf,%lf} weight={%le,%le}\n",theta/PI*180,thetarange[0]/PI*180,thetarange[1]/PI*180,weight0,weight);
            if(theta>thetarange[0]&&theta<thetarange[1]&&whichtel>=0){
               int findres=FindPhiRange(zero,coor_gen,dir_gen,dir_tel,freelength,theta,phirange);
               if(jdebug>7||(!isfinite(weight))) printf("Laser::Propagate: find Rayley scatter phi range, itel=%d freelength=%le theta=%lf phi range={%le,%le} return=%d\n",whichtel,freelength,theta/PI*180,phirange[0]/PI*180,phirange[1]/PI*180,findres);
               if(!(findres>0)) {phirange[0]=phirange[1]=-100;}
            }
            weight0=weight;
            Atmosphere::RayScatterAnglePhi(wavelength_gen,phi,phirange,weight);
            if(jdebug>7||(!isfinite(weight))) printf("Laser::Propagate: generated Ray phi=%lf phirange={%lf,%lf} weight={%le,%le}\n",phi/PI*180,phirange[0]/PI*180,phirange[1]/PI*180,weight0,weight);
         }
         else if(scatter==2){	//Mie scattering
            weight0=weight;
            Atmosphere::MieScatterAngleTheta(wavelength_gen,theta,thetarange,weight);
            if(jdebug>6||(!isfinite(weight))) printf("Laser::Propagate: generated Mie theta=%lf thetarange={%lf,%lf} weight={%le,%le}\n",theta/PI*180,thetarange[0]/PI*180,thetarange[1]/PI*180,weight0,weight);
            if(theta>thetarange[0]&&theta<thetarange[1]&&whichtel>=0){
               int findres=FindPhiRange(zero,coor_gen,dir_gen,dir_tel,freelength,theta,phirange);
               if(jdebug>7||(!isfinite(weight))) printf("Laser::Propagate: find Mie scatter phi range, itel=%d freelength=%le theta=%lf phi range={%le,%le} return=%d\n",whichtel,freelength,theta/PI*180,phirange[0]/PI*180,phirange[1]/PI*180,findres);
               if(!(findres>0)) {phirange[0]=phirange[1]=-100;}
            }
            weight0=weight;
            Atmosphere::MieScatterAnglePhi(wavelength_gen,phi,thetarange,weight);
            if(jdebug>7||(!isfinite(weight))) printf("Laser::Propagate: generated Mie phi=%lf phirange={%lf,%lf} weight={%le,%le}\n",phi/PI*180,phirange[0]/PI*180,phirange[1]/PI*180,weight0,weight);
         }
         double xdir[3],ydir[3],zdir[3];
         double dir222[3]={0,0,1};
         CartesianFrame(zero,coor_gen,dir_gen,dir222,xdir,ydir,zdir);
         for(int ii=0;ii<3;ii++){
            dir_scat[ii]=cos(theta)*zdir[ii]+sin(theta)*(cos(phi)*xdir[ii]+sin(phi)*ydir[ii]);
         }
         double norm=sqrt(dir_scat[0]*dir_scat[0]+dir_scat[1]*dir_scat[1]+dir_scat[2]*dir_scat[2]);
         for(int ii=0;ii<3;ii++) dir_scat[ii]/=norm;

         //the distance to the closest telescope
         dist=mindist(zero,coor_scat,dir_scat,coor_min,decrease);
         distance=freelength+sqrt(pow(coor_scat[0]-coor_min[0],2)+pow(coor_scat[1]-coor_min[1],2)+pow(coor_scat[2]-coor_min[2],2));
         if(jdebug>4||(!isfinite(weight))) printf("Laser::Propagate: scattered, dist=%lf(free length=%lf distance=%lf) decrease=%d coo={%lf,%lf,%lf} weight=%le\n",dist,freelength,distance,decrease,coor_min[0],coor_min[1],coor_min[2],weight);
         int returntype;
         double lengthrange2[2]={0,0};
         double freelength2=Atmosphere::FreePathLength(coor_scat[2],dir_scat,lengthrange2,weight);
         if(decrease&&dist<TelSimDist){ //inside the field of view of one telescope after scattering
            if((fabs(distance-freelength)<freelength2&&Doextin)||(!Doextin)){
               Telindex=whichtel;
               for(int ii=0;ii<3;ii++){
                  coor_out[ii]=coor_min[ii];
                  dir_out[ii]=dir_scat[ii];
               }
               returntype=scatter; //pass through nearby of the telescope after one scattering
            }
            else returntype=-1; //interacted before arriving to the telescope, after one scattering
         }
         else{ //too far away from the telescope after one scattering
            returntype=-2;
         }
         if(DoPlot){
            for(int iline=0;iline<2;iline++){
               if(IniRange[3][1]>IniRange[3][0]){
                  if(IniRange[3][1]<=0) break;
                  if(iline==0&&IniRange[3][0]>=freelength/vlight) continue;
                  if(iline==1&&IniRange[3][1]<=freelength/vlight) continue;
               }
               double coor_last[4],coor_first[4];
               if(IniRange[3][1]>IniRange[3][0]){
                  if(iline==0){
                     double length_prop1=TMath::Max((double)0,(double)IniRange[3][0])*vlight;
                     double length_prop2=TMath::Max((double)0,(double)IniRange[3][1])*vlight;
                     double maxlength_prop=freelength;
                     for(int ii=0;ii<3;ii++) coor_first[ii]=coor_gen[ii]+dir_gen[ii]*TMath::Min(length_prop1,maxlength_prop);
                     for(int ii=0;ii<3;ii++) coor_last[ii]=coor_gen[ii]+dir_gen[ii]*TMath::Min(length_prop2,maxlength_prop);
                  }
                  else{
                     double length_prop1=TMath::Max((double)0,(double)IniRange[3][0]*vlight-freelength);
                     double length_prop2=TMath::Max((double)0,(double)IniRange[3][1]*vlight-freelength);
                     double maxlength_prop=freelength2;
                     if(decrease&&dist<TelSimDist) maxlength_prop=TMath::Min(distance-freelength,freelength2);
                     for(int ii=0;ii<3;ii++) coor_first[ii]=coor_scat[ii]+dir_scat[ii]*TMath::Min(length_prop1,maxlength_prop);
                     for(int ii=0;ii<3;ii++) coor_last[ii]=coor_scat[ii]+dir_scat[ii]*TMath::Min(length_prop2,maxlength_prop);
                  }
                  coor_first[3]=TMath::Max((double)0,(double)IniRange[3][0]);
                  coor_last[3]=TMath::Max((double)0,(double)IniRange[3][1]);
               }
               else{
                  if(iline==0){
                     for(int ii=0;ii<3;ii++) coor_first[ii]=coor_gen[ii];
                     for(int ii=0;ii<3;ii++) coor_last[ii]=coor_scat[ii];
                  }
                  else{
                     for(int ii=0;ii<3;ii++) coor_first[ii]=coor_scat[ii];
                     for(int ii=0;ii<3;ii++){
                        coor_last[ii]=coor_scat[ii]+dir_scat[ii]*((decrease&&dist<TelSimDist)?TMath::Min(distance-freelength,freelength2):freelength2);
                     }
                  }
                  coor_first[3]=0;
                  coor_last[3]=((decrease&&dist<TelSimDist)?freelength+TMath::Min(distance-freelength,freelength2):freelength+freelength2)/vlight;
               }
               Add(coor_first,coor_last,returntype,weight);
            }
         }
         return returntype;
      }
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
      double avet=0;
      int nt=0;
      for(int il=0;il<LaserSize;il++){
         double x0,y0,z0;
         double m1,n1,l1;
         double wave;
         int itube,icell;
         int whichtel;
         double weight;

         whichtel=votel.at(il);
         weight=vowei.at(il);
         if(WFCTALaserEvent::Recordweight) (pwfc->laserevent).weight.push_back(weight);
         (pwfc->laserevent).hweight->Fill(log10(weight));
         if(whichtel<0||whichtel>=WFTelescopeArray::CTNumber){
            vosipm.push_back(-1);
            if(WFCTAMCEvent::RecordRayTrace) (pwfc->mcevent).RayTrace.push_back(whichtel);
            (pwfc->mcevent).hRayTrace->Fill(whichtel,weight);
            //printf("pushing back tracing il=%d whichtel=%d\n",il,whichtel);
            continue;
         }
         WFTelescope* pt=pct->pct[whichtel];
         if(!pt) {vosipm.push_back(-1); continue;}

         x0=vocoo[0].at(il);
         y0=vocoo[1].at(il);
         z0=vocoo[2].at(il);
         m1=vodir[0].at(il);
         n1=vodir[1].at(il);
         l1=vodir[2].at(il);
         wave=vgwav.at(il);
         double t=votim.at(il); //in second
         double tt=t;
         //apply the telescope camera SiPM quantum efficiency
         if(wave>0){
            if(!WFTelescope::GetQuantumEff(wave)){
               vosipm.push_back(-1);
               if(WFCTAMCEvent::RecordRayTrace) (pwfc->mcevent).RayTrace.push_back(-14);
               (pwfc->mcevent).hRayTrace->Fill(-14,weight);
               continue;
            }
         }
         int res=pct->RayTrace(x0,y0,z0,m1,n1,l1,weight,wave,whichtel,t,itube,icell);
         if(WFCTAMCEvent::RecordRayTrace) (pwfc->mcevent).RayTrace.push_back(res);
         (pwfc->mcevent).hRayTrace->Fill(res,weight);
         if(res>=0){
            vosipm.push_back(itube);
            findtel=true;
            avet+=t;
            nt++;
            if(jdebug>4) printf("Laser::DoWFCTASim: the photon go through the pmt %d, iphoton=%d time={%.9le,%.9le}\n",itube,il,tt,t); //jdebug>4
         }
         else vosipm.push_back(-1);
         //printf("pushing back tracing il=%d whichtel=%d res=%d\n",il,whichtel,res);
         pwfc->iTel=(whichtel>=0)?pct->pct[whichtel]->TelIndex_:-1;
      }
      if(pwfc){
         if(jdebug>1) printf("Laser::DoWFCTASim: Filling the event %d\n",ievent_gen);
         pwfc->iEvent=ievent_gen;
         double t0=(nt==0)?(Time_gen+time_gen*20*1.0e-9):(avet/nt);
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
            pwfc->CalculateDataVar();
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

void Laser::Add(double coorin1[4],double coorin2[4],int type,double weight){
   double coo1[4],coo2[4];
   for(int ii=0;ii<4;ii++){
      coo1[ii]=coorin1[ii];
      coo2[ii]=coorin2[ii];
   }
   bool doshrink=false,dropthis=false;
   for(int ii=0;ii<3;ii++){
      if(IniRange[ii][0]<IniRange[ii][1]){
         if(coo1[ii]<=IniRange[ii][0]&&coo2[ii]<=IniRange[ii][0]){ dropthis=true; break;}
         if(coo1[ii]>=IniRange[ii][1]&&coo2[ii]>=IniRange[ii][1]){ dropthis=true; break;}
         if((coo1[ii]>=IniRange[ii][0]&&coo1[ii]<=IniRange[ii][1])&&(coo2[ii]>=IniRange[ii][0]&&coo2[ii]<=IniRange[ii][1])) continue;
         double zz1=coo1[ii],zz2=coo2[ii];
         if(coo1[ii]<IniRange[ii][0]){
            zz1=IniRange[ii][0];
            if(coo2[ii]>IniRange[ii][1]) zz2=IniRange[ii][1];
         }
         else if(coo1[ii]<=IniRange[ii][1]){
            if(coo2[ii]<IniRange[ii][0]) zz2=IniRange[ii][0];
            if(coo2[ii]>IniRange[ii][1]) zz2=IniRange[ii][1];
         }
         else{
            zz1=IniRange[ii][1];
            if(coo2[ii]<IniRange[ii][0]) zz2=IniRange[ii][0];
         }
         if(zz1!=coo1[ii]){
            for(int i2=0;i2<3;i2++){
               if(i2==ii) continue;
               coo1[i2]=coo1[i2]+(coo2[i2]-coo1[i2])/(coo2[ii]-coo1[ii])*(zz1-coo1[ii]);
            }
            coo1[ii]=zz1;
            doshrink=true;
         }
         if(zz2!=coo2[ii]){
            for(int i2=0;i2<3;i2++){
               if(i2==ii) continue;
               coo2[i2]=coo1[i2]+(coo2[i2]-coo1[i2])/(coo2[ii]-coo1[ii])*(zz2-coo1[ii]);
            }
            coo2[ii]=zz2;
            doshrink=true;
         }
      }
   }
   if(dropthis) return;
   if(!doshrink){
      for(int ii=0;ii<3;ii++){
         coo1[ii]=coorin1[ii];
         coo2[ii]=coorin2[ii];
      }
   }

   for(int ii=0;ii<4;ii++){
      if(IniRange[ii][0]<IniRange[ii][1]) continue;
      if(coo1[ii]<plotrange[ii][0]) plotrange[ii][0]=coo1[ii];
      if(coo2[ii]<plotrange[ii][0]) plotrange[ii][0]=coo2[ii];
      if(coo1[ii]>plotrange[ii][1]) plotrange[ii][1]=coo1[ii];
      if(coo2[ii]>plotrange[ii][1]) plotrange[ii][1]=coo2[ii];
   }

   //do some scale for the plot;
   double plotscale=1.0e-2;//(weight<1.0e5)?1.e-3:1;
   if(prandom->Rndm()>plotscale) return;
   TPolyLine3D* line=new TPolyLine3D(2);
   line->SetPoint(0, coo1[0],coo1[1],coo1[2]);
   line->SetPoint(1, coo2[0],coo2[1],coo2[2]);
   line->SetLineColor(GetLineColor(type,weight));
   line->SetLineStyle(GetLineStyle(type,weight));
   line->SetLineWidth(GetLineWidth(type,weight));
   if(!plot) plot=new TObjArray();
   plot->Add(line);
   if(jdebug>5) printf("Laser::Add: Add line, type=%d weight={%le,%le} coo1={%le,%le,%le,%le},coo2={%le,%le,%le,%le}\n",type,weight,plotscale,coo1[0],coo1[1],coo1[2],coo1[3],coo2[0],coo2[1],coo2[2],coo2[3]);
}
int Laser::GetLineStyle(int type,double weight){
   return 1;
}
int Laser::GetLineWidth(int type,double weight){
   //if(weight<=0) return 0;
   //if(log10(weight)<3) return 1;
   //else return 2;

   if(weight<=0) return 0;
   if(log10(weight)<5) return 1;
   else if(log10(weight)<8) return 3;
   else if(log10(weight)<10) return 5;
   else return 7;
}
int Laser::GetLineColor(int type,double weight){
   if(type<0) return 2;
   else if(type==0) return 3;
   else return 4;
}
TCanvas* Laser::Draw(const char* option,int ViewOpt,const char* savedir){
   TCanvas* cc = new TCanvas(Form("Laser_Propagation_Evt%d",ievent_gen),(ViewOpt<1||ViewOpt>3)?"Inclined View":(ViewOpt==1?"Front View":(ViewOpt==2?"Side View":"Top View")),ViewOpt==3?3000:2000,3000);
   double rmin[3]={plotrange[0][0],plotrange[1][0],plotrange[2][0]};
   double rmax[3]={plotrange[0][1],plotrange[1][1],plotrange[2][1]};
   TView *view = TView::CreateView(1,rmin,rmax);
   view->ShowAxis();
   if(ViewOpt==1) view->Front();
   else if(ViewOpt==2) view->Side();
   else if(ViewOpt==3) view->Top();
   cc->SetView(view);
   if(jdebug>1) printf("Laser::Draw: Pos Range={{%+6.1e,%+6.1e},{%+6.1e,%+6.1e},{%+6.1e,%+6.1e}} Time Range={%+7.1e,%+7.1e}\n",rmin[0],rmax[0],rmin[1],rmax[1],rmin[2],rmax[2],plotrange[3][0],plotrange[3][1]);   
   if(plot){
      WFTelescopeArray* pta=WFTelescopeArray::GetHead();
      for(int itel=0;itel<WFTelescopeArray::CTNumber;itel++){
         WFTelescope* pt=(pta)?pta->pct[itel]:0;
         if(!pt) continue;
         TPolyLine3D* line=new TPolyLine3D(2);
         double postel[3]={pt->Telx_,pt->Tely_,pt->Telz_};
         double dirtel[3]={sin(pt->TelZ_)*cos(pt->TelA_),sin(pt->TelZ_)*sin(pt->TelA_),cos(pt->TelZ_)};
         double xyzmax[3]={TMath::Max(fabs(plotrange[0][0]),fabs(plotrange[0][1])),TMath::Max(fabs(plotrange[1][0]),fabs(plotrange[1][1])),TMath::Max(fabs(plotrange[2][0]),fabs(plotrange[2][1]))};
         double length_dir=sqrt(pow(xyzmax[0],2)+pow(xyzmax[1],2)+pow(xyzmax[2],2));
         line->SetPoint(0,postel[0],postel[1],postel[2]);
         line->SetPoint(1,postel[0]+dirtel[0]*length_dir,postel[1]+dirtel[1]*length_dir,postel[2]+dirtel[2]*length_dir);
         line->SetLineColor(1);
         line->SetLineStyle(2);
         line->SetLineWidth(6);
         plot->Add(line);
      }
      plot->Draw(option);
   }
   TAxis3D *axis = TAxis3D::GetPadAxis();
   axis->SetLabelColor(kBlack);
   axis->SetAxisColor(kBlack);
   axis->SetLabelSize(0.03,"X");
   axis->SetLabelSize(0.03,"Y");
   axis->SetLabelSize(0.03,"Z");
   axis->SetXTitle("X [cm]");
   axis->SetYTitle("Y [cm]");
   axis->SetZTitle("Z [cm]");
   axis->SetTitleOffset(2.,"X");
   axis->SetTitleOffset(2.,"Y");
   axis->SetTitleOffset(2.,"Z");
   if(ViewOpt<1||ViewOpt>3){
   axis->SetLabelOffset(0.006,"X");
   axis->SetLabelOffset(0.008,"Y");
   axis->SetLabelOffset(0.008,"Z");
   axis->SetTitleOffset(1.2,"X");
   axis->SetTitleOffset(2.0,"Y");
   axis->SetTitleOffset(1.2,"Z");
   }
   else if(ViewOpt==1){
   axis->SetLabelOffset(-0.08,"X");
   axis->SetLabelOffset(0.03,"Z");
   axis->SetTitleOffset(1.3,"X");
   }
   else if(ViewOpt==2){
   axis->SetLabelOffset(-0.08,"Y");
   axis->SetLabelOffset(0.03,"Z");
   axis->SetTitleOffset(1.3,"Y");
   }
   else if(ViewOpt==3){
   axis->SetLabelOffset(0.03,"X");
   axis->SetLabelOffset(0.03,"Y");
   }

   if(savedir) cc->SaveAs(Form("%s/%d.png",savedir,ievent_gen));
   return cc;
}

int Laser::GetProb(long int ngen){
   TH1::SetDefaultSumw2();
   Reset();
   Time_gen=10000;
   time_gen=0;
   long int ngentel=0;
   for(long int igen=0;igen<ngen;igen++){
      bool dogen=InitialGen();
      if(!dogen) continue;
      iphoton_gen=igen;
      vgwav.push_back(wavelength_gen);
      for(int ii=0;ii<3;ii++){
         vgcoo[ii].push_back(coor_gen[ii]);
         vgdir[ii].push_back(dir_gen[ii]);
      }
      double weight=1.;
      double distance;
      int res=Propagate(distance,weight);
      volength.push_back(interpoint);
      volength2.push_back(distance);
      if(res<0) Telindex=res-15;
      else{  //the telescope index has been calculated in Propagate
         ngentel++;
      }
      vowei.push_back(weight);
      votim.push_back(distance/vlight);
      votel.push_back(Telindex);
      for(int ii=0;ii<3;ii++){
         vocoo[ii].push_back(res>=0?coor_out[ii]:0);
         vodir[ii].push_back(res>=0?dir_out[ii]:0);
      }
      count_gen+=weight;
      if((igen%10000)==0) printf("%ld of %ld generated\n",igen,ngen);
   }
   bool dosim=DoWFCTASim();

   if(!hdenu){
      hdenu=new TH1D("hdenum",";distance [m];Events",1000,lengthmin-0.5,lengthmax-0.5);
   }
   int size=vocoo[0].size();
   int res=0;
   for(int ii=0;ii<size;ii++){
      double weight=vowei.at(ii);
      double length1=volength.at(ii)/100.;
      double length2=volength2.at(ii)/100.;
      int isipm=vosipm.at(ii);
      int whichtel=votel.at(ii);
      hdenu->Fill(length1,weight);
      if(isipm<0||isipm>=NSIPM) continue;
      if(whichtel<0||whichtel>=NCTMax) continue;
      if(!hprob[whichtel][isipm]){
         hprob[whichtel][isipm]=(TH1D*)hdenu->Clone(Form("prob_sipm%d_tel%d",isipm,whichtel));
         hprob[whichtel][isipm]->Reset();
      }
      hprob[whichtel][isipm]->Fill(length1,weight);
      if(!hleng[whichtel][isipm]){
         //hleng[whichtel][isipm]=(TH1D*)hdenu->Clone(Form("leng_sipm%d_tel%d",isipm,whichtel));
         //hleng[whichtel][isipm]->Reset();
         hleng[whichtel][isipm]=new TH1D(Form("leng_sipm%d_tel%d",isipm,whichtel),";Distance [m];Events",1000,lengthmin2-0.5,lengthmax2-0.5);
      }
      hleng[whichtel][isipm]->Fill(length2,weight);
      res++;
   }

   /*for(int ii=0;ii<NCTMax;ii++){
      for(int jj=0;jj<NSIPM;jj++){
         if(hprob[ii][jj]){
            TH1D* hden=hdenu;
            TH1D* hnum=hprob[ii][jj];
            for(int ibin=1;ibin<hden->GetNbinsX();ibin++){
               double n1=hden->GetBinContent(ibin);
               double en1=hden->GetBinError(ibin);
               double n2=hnum->GetBinContent(ibin);
               double en2=hnum->GetBinError(ibin);
               double content=(n1==0)?0:(n2/n1);
               double econtent=(n1==0)?0:sqrt(pow((n1-n2)/n1*en2,2)+pow(n2/n1,2)*(en1*en1-en2*en2))/n1;
               hprob[ii][jj]->SetBinContent(ibin,content);
               hprob[ii][jj]->SetBinError(ibin,econtent);
            }
         }
      }
   }*/

   ievent_gen++;
   return res;
}
