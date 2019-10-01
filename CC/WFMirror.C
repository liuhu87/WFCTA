#include "WFMirror.h"
#include "TMath.h"
#include "TRandom.h"
#include "common.h"
#include "WFTelescope.h"
int WFMirrorArray::NMirror[NCircle]={1,6,6,6,2,2,2};
void WFMirrorArray::Init(){
   for(int ii=0;ii<NCircle;ii++){
      for(int jj=0;jj<NMaxMirror;jj++){
         if(jj<NMirror[ii]) pmirr[ii][jj]=new WFMirror();
         else pmirr[ii][jj]=0;
      }
   }
   SetMirror();
}
void WFMirrorArray::Clear(){
   for(int ii=0;ii<NCircle;ii++){
      for(int jj=0;jj<NMaxMirror;jj++){
         if(pmirr[ii][jj]) {delete pmirr[ii][jj]; pmirr[ii][jj]=0;}
      }
   }
}
void WFMirrorArray::SetMirror(){
    int num[7];
    double fai[7],theta[7],alfa[7],delt[7];
    double omcenterx[7][6],omcentery[7][6],omcenterz[7][6];
    double mfai, mtheta, malfa,delta, ntheta;
    double ndx, ndy, ndz;
    double xc, yc, zc;
    double z2,x2,z3,y3;
    double _L = WFMirror::length;
    double CURVATURE=WFMirror::CURVATURE;
    //double PI = TMath::Pi();
    double sqrt3 = sqrt(3);

    fai[0]=PI;                               theta[0]=0;          num[0]=1;
    fai[1]=PI-atan(sqrt3*_L/CURVATURE);      theta[1]=0;          num[1]=6;
    fai[2]=PI-atan(2.*sqrt3*_L/CURVATURE);   theta[2]=0;          num[2]=6;
    fai[3]=PI-atan(3.*_L/CURVATURE);         theta[3]=PI/6.;      num[3]=6;
    fai[4]=PI-atan(sqrt(4.5*4.5+3./4.)*_L/CURVATURE);      theta[4]=PI/6.+atan(sqrt3/9);   num[4]=2;
    fai[5]=PI-atan(sqrt(4.5*4.5+3./4.)*_L/CURVATURE);      theta[5]=5./6.*PI-atan(sqrt3/9);num[5]=2;
    fai[6]=PI-atan(sqrt(4.5*4.5+3./4.)*_L/CURVATURE);      theta[6]=PI/6.-atan(sqrt3/9);   num[6]=2;
    for(int i=0;i<NCircle;i++){
       // zenith angle for each circle in array
       mfai=fai[i];
       mtheta=theta[i];
       malfa=atan(sqrt3*_L/4/CURVATURE);
       for(int m=0;m<NMirror[i];m++){
          // azimuth angle for each mirror in array
          if(i<4) ntheta=mtheta+2*PI/6.*m;
          else{
             ntheta=(m==0)?mtheta:(2*PI-mtheta);
          }
          omcenterx[i][m]=CURVATURE*sin(mfai)*cos(ntheta);
          omcentery[i][m]=CURVATURE*sin(mfai)*sin(ntheta);
          omcenterz[i][m]=CURVATURE*cos(mfai);
          WFMirror::altercoor(&malfa,&omcenterx[i][m],&omcenterz[i][m],&pmirr[i][m]->mirrorx,&pmirr[i][m]->mirrorz);
          pmirr[i][m]->mirrory=omcentery[i][m];
          pmirr[i][m]->nalfa=atan(-pmirr[i][m]->mirrorx/pmirr[i][m]->mirrorz);
          pmirr[i][m]->nbeta=atan(pmirr[i][m]->mirrory/sqrt(pow(pmirr[i][m]->mirrorx,2)+pow(pmirr[i][m]->mirrorz,2)));
          pmirr[i][m]->mirrorz+=CURVATURE;
       }
    }
}
void WFMirrorArray::SetMirrorPointError(int flag,double error)
{
    for(int i=0; i<NCircle; i++){
        for(int j=0; j<NMaxMirror; j++){
            if(pmirr[i][j]) pmirr[i][j]->SetMirrorPointError(flag,error);
        }
    }
}
bool WFMirrorArray::WhichMirror(double x, double y, double z, double *deltax, double *deltay,double *deltaz,int *ii, int *mm){
    for(int i=0;i<NCircle;i++){                   //number of the circle
        for(int m=0;m<NMirror[i];m++){
            if(pmirr[i][m]->IsInside(x,y,z,deltax,deltay,deltaz)){
                *ii = i;
                *mm = m;
                return true;
            }
        }
    }
    return false;
}

double WFMirror::MirrorSpot;
double WFMirror::MirrorPointError;
double WFMirror::MirrorPointErrorFlag;
double WFMirror::MirrorGeometry;
double WFMirror::length=300.;
double WFMirror::length_D=150.;
double WFMirror::CURVATURE=5740;//5800.;
double WFMirror::REFLECTIVITY=0.83;
void WFMirror::Init(){
   mirrorx=-10000;
   mirrory=-10000;
   mirrorz=-10000;
}
void WFMirror::SetEulerMatrix(double theta,double phi,double matrix_[3][3])
{
   double cosz, sinz, cosa, sina;
   cosz = cos(theta);
   sinz = sin(theta);
   cosa = cos(phi);
   sina = sin(phi);
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
void WFMirror::InverseEuler(double x0, double y0, double z0, double *x, double *y, double *z,double matrix_[3][3])
{
   *x = matrix_[0][0]*x0+matrix_[1][0]*y0+matrix_[2][0]*z0;
   *y = matrix_[0][1]*x0+matrix_[1][1]*y0+matrix_[2][1]*z0;
   *z = matrix_[0][2]*x0+matrix_[1][2]*y0+matrix_[2][2]*z0;
}
void WFMirror::SetMirrorSpot(float mirrorspot)
{ 
  MirrorSpot = mirrorspot;
}
void WFMirror::SetMirrorPointErrorFlag(int point_error_flag)
{
  MirrorPointErrorFlag = point_error_flag;
}
void WFMirror::SetMirrorPointError(float point_error)
{
  MirrorPointError = point_error;
}
void WFMirror::SetMirrorPointError(int flag,double error)
{
  double theta0, phi0;
  double theta , phi;
  double x, y, z, r;
  double xnew ,ynew, znew;
  double sigma2;
  double m2, n2, l2, m1, n1, l1, m0, n0, l0,norm;

  if(flag){
    error *= TMath::DegToRad();
    sigma2 = error * error;

    if(mirrorx==-10000) return;
       x = 0 - mirrorx;
       y = 0 - mirrory;
       z = CURVATURE - mirrorz;
       r = sqrt(x*x+y*y+z*z);
       theta0 = acos(z/r);
       phi0 = atan(y/x);
       if(x>0) phi0 = phi0;
       if(x<0) phi0 = phi0 + TMath::Pi();
       double matrix_[3][3];
       SetEulerMatrix(theta0, phi0,matrix_);

       theta =  WFTelescopeArray::prandom->Rndm();
       theta = sqrt(-log(1-theta)*2*sigma2);
       phi = WFTelescopeArray::prandom->Rndm()*TMath::TwoPi();

       x = CURVATURE*sin(theta)*cos(phi) ;
       y = CURVATURE*sin(theta)*sin(phi) ;
       z = CURVATURE*cos(theta) ;

       InverseEuler( x,  y, z, &xnew, &ynew, &znew,matrix_);
       xnew += mirrorx;
       ynew += mirrory;
       znew += mirrorz;
       spherecenterx = xnew;
       spherecentery = ynew;
       spherecenterz = znew;
   }

   else{
      printf("The Pointing Errors of the each mirror are not simulated\n");
      spherecenterx = 0;
      spherecentery = 0;
      spherecenterz = CURVATURE;
  }

//=== The spherecenterx,y,z, are the sphere center of each small mirror      ===//
////=== The value of the three arrays should be provided by the measurements   ===//
////=== Currently, without the measurments results, all small mirrors can be   ===//
////=== considered to be the same and the sphere center to be (0, 0, CURVATURE)===//
////=== Lingling Ma, 20151014  ===//

}
void WFMirror::SetMirrorGeometry(int mirror_geometry)
{
  MirrorGeometry = mirror_geometry;
}
void WFMirror::altercoor(double *alfa,double *ox,double *oy,double *nx,double *ny)
{
  *nx=*oy*sin(*alfa)+*ox*cos(*alfa);
  *ny=*oy*cos(*alfa)-*ox*sin(*alfa);
}
void WFMirror::ConvertCoor(double *ox,double *oy,double *oz,double *nx,double *ny,double *nz){
   double ndx, ndy, ndz;
   double n2dx, n2dy, n2dz;
   altercoor(&nalfa,ox,oz,&ndx,&ndz);
   altercoor(&nbeta,oy,&ndz,&n2dy,&n2dz);
   *nx=ndx;
   *ny=n2dy;
   *nz=n2dz;
}
bool WFMirror::IsInside(double x,double y,double z,double *deltax,double *deltay,double *deltaz){
   if(mirrorx==-10000&&mirrory==-10000&&mirrorz==-10000) return false;
   double dx,dy,dz;
   dx=mirrorx - x;
   dy=mirrory - y;
   dz=mirrorz - z;
   double xc,yc,zc;
   ConvertCoor(&dx,&dy,&dz,&xc,&yc,&zc);
   double _L=length;
   double sqrt3 = sqrt(3);
   if((fabs(yc)<=-fabs(xc)/sqrt3+_L)&&fabs(xc)<=sqrt3/2*_L){
      *deltax=xc;
      *deltay=yc;
      *deltaz=zc;
      return true;
   }
   else return false;
}
void WFMirror::GetReflected(double x1,double y1,double z1,double x0,double y0,double z0,double *m2,double *n2,double *l2)
{
 // x0,y0,z0 are the impact point
 //  x, y, z are the norml directions
 double x,y,z;
 x=spherecenterx-x0;
 y=spherecentery-y0;
 z=spherecenterz-z0;
 double norm=sqrt(x*x+y*y+z*z);
 x/=norm;
 y/=norm;
 z/=norm;
 //   //  x1, y1, z1, are the incendent directions
   double costheta = x1*x+y1*y+z1*z;
   *m2 = (x1-2.*costheta*x);
   *n2 = (y1-2.*costheta*y);
   *l2 = (z1-2.*costheta*z);
}
