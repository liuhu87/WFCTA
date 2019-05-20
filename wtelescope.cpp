#include "wtelescope.h"
#include "cer_event.h"
#include "TMath.h"
#include "TRandom.h"
#include "TGraph2D.h"
//* ==== WFCTA Telescope functions ===== *//

WTelescope::WTelescope()
{

}

WTelescope::~WTelescope()
{ 

}
void WTelescope::SetMirrorSpot(float mirrorspot)
{
  MirrorSpot = mirrorspot; 
}

void WTelescope::SetMirrorPointErrorFlag(int point_error_flag)
{
  MirrorPointErrorFlag = point_error_flag;
}

void WTelescope::SetMirrorPointError(float point_error)
{
  MirrorPointError = point_error;
}

void WTelescope::SetMirrorGeometry(int mirror_geometry)
{
  MirrorGeometry = mirror_geometry;
}

void WTelescope::SetMirror()
{
    static double fai[7],theta[7],alfa[7],delt[7];
    static double omcenterx[7][6],omcentery[7][6],omcenterz[7][6];
    static int num[7];
    double mfai, mtheta, malfa,delta, ntheta;
    double ndx, ndy, ndz;
    double xc, yc, zc;
    double z2,x2,z3,y3;
    double _L = length;
    fai[0]=pi;                         theta[0]=0;          num[0]=1;
    fai[1]=pi-sqrt3*_L/CURVATURE;      theta[1]=0;          num[1]=6;
    fai[2]=pi-2.*sqrt3*_L/CURVATURE;   theta[2]=0;          num[2]=6;
    fai[3]=pi-3.*_L/CURVATURE;         theta[3]=pi/6.;      num[3]=6;
    fai[4]=pi-2.*sqrt3*_L/CURVATURE;   theta[4]=pi/3.;      num[4]=2;
    fai[5]=pi-2.*sqrt3*_L/CURVATURE;   theta[5]=2.*pi/3.;   num[5]=2;
    fai[6]=pi-3.*_L/CURVATURE;         theta[6]=pi/6.;      num[6]=2;

    alfa[0]=-sqrt3*_L/4/CURVATURE;     delt[0]=pi/3;
    alfa[1]=-sqrt3*_L/4/CURVATURE;     delt[1]=pi/3;
    alfa[2]=-sqrt3*_L/4/CURVATURE;     delt[2]=pi/3;
    alfa[3]=-sqrt3*_L/4/CURVATURE;     delt[3]=pi/3;
    alfa[4]=sqrt3*3*_L/4/CURVATURE;    delt[4]=4.*pi/3;
    alfa[5]=-sqrt3*5*_L/4/CURVATURE;   delt[5]=2.*pi/3;
    alfa[6]=sqrt3*3*_L/4/CURVATURE;    delt[6]=5.*pi/3;  

    for(int i=0;i<7;i++){
        for(int m=0;m<6;m++){
            mirrorx[i][m] = -10000;
            mirrory[i][m] = -10000;
            mirrorz[i][m] = -10000;
        }
    }

    for(int i=0;i<7;i++){                   //number of the circle
        mfai=fai[i]; mtheta=theta[i]; malfa=alfa[i]; delta=delt[i];
        for(int m=0;m<num[i];m++){          //the mirror's number of every circle
            ntheta=mtheta+m*delta;
            omcenterx[i][m]=CURVATURE*sin(mfai)*cos(ntheta);
            omcentery[i][m]=CURVATURE*sin(mfai)*sin(ntheta);
            omcenterz[i][m]=CURVATURE*cos(mfai);
            altercoor(&malfa,&omcenterz[i][m],&omcenterx[i][m],&mirrorz[i][m],&mirrorx[i][m]);
            mirrory[i][m]=omcentery[i][m];

            nalfa[i][m]=atan(mirrorx[i][m]/mirrorz[i][m]);
            nbeta[i][m]=atan(mirrory[i][m]/sqrt(mirrorz[i][m]*mirrorz[i][m]+mirrorx[i][m]*mirrorx[i][m]));
            mirrorz[i][m]= mirrorz[i][m]+CURVATURE;
            //printf(" %d %f %f %f\n",i, mirrorx[i][m],mirrory[i][m],mirrorz[i][m]);
        }
    }    
}

void WTelescope::SetMirror(int i)
{
   static double fai[7],theta[7],alfa[7],delt[7];
    static double omcenterx[7][6],omcentery[7][6],omcenterz[7][6];
    static int num[7];
    double mfai, mtheta, malfa,delta, ntheta;
    double ndx, ndy, ndz;
    double xc, yc, zc;
    double z2,x2,z3,y3;
    double _L = length; //length=
    double pi = TMath::Pi();
    double sqrt3 = sqrt(3);
   
    fai[0]=pi;                         theta[0]=0;          num[0]=1;
    fai[1]=pi-sqrt3*_L/CURVATURE;      theta[1]=0;          num[1]=6;
    fai[2]=pi-2.*sqrt3*_L/CURVATURE;   theta[2]=0;          num[2]=6;
    fai[3]=pi-3.*_L/CURVATURE;         theta[3]=pi/6.;      num[3]=6;
    fai[4]=pi-4.58582575*_L/CURVATURE;      theta[4]=atan(sqrt3/2);   num[4]=2;
    fai[5]=pi-4.58582575*_L/CURVATURE;      theta[5]=pi-atan(sqrt3/2); num[5]=2;
    fai[6]=pi-4.58582575*_L/CURVATURE;         theta[6]=atan(sqrt3/5);   num[6]=2;

    alfa[0]=-sqrt3*_L/4/CURVATURE;     delt[0]=pi/3;
    alfa[1]=-sqrt3*_L/4/CURVATURE;     delt[1]=pi/3;
    alfa[2]=-sqrt3*_L/4/CURVATURE;     delt[2]=pi/3;
    alfa[3]=-sqrt3*_L/4/CURVATURE;     delt[3]=pi/3;
    alfa[4]=-sqrt3*_L/4/CURVATURE;     delt[4]=2*pi-2*atan(sqrt3/2);
    alfa[5]=-sqrt3*_L/4/CURVATURE;     delt[5]=2.*atan(sqrt3/2);
    alfa[6]=-sqrt3*_L/4/CURVATURE;     delt[6]=2*pi-2*atan(sqrt3/5);    
   

    for(int i=0;i<7;i++){
        for(int m=0;m<6;m++){
            mirrorx[i][m] = -10000;
            mirrory[i][m] = -10000;
            mirrorz[i][m] = -10000;
        }
    }
    
      for(int i=0;i<7;i++){                   //number of the circle
         mfai=fai[i];
         mtheta=theta[i];
         malfa=alfa[i];
         delta=delt[i];
        for(int m=0;m<num[i];m++){          //the mirror's number of every circle
            ntheta=mtheta+m*delta;          // azimuth angle for each mirror in array
            omcenterx[i][m]=CURVATURE*sin(mfai)*cos(ntheta);
            omcentery[i][m]=CURVATURE*sin(mfai)*sin(ntheta);
            omcenterz[i][m]=CURVATURE*cos(mfai);
            altercoor(&malfa,&omcenterz[i][m],&omcenterx[i][m],&mirrorz[i][m],&mirrorx[i][m]);
            mirrory[i][m]=omcentery[i][m];

            nalfa[i][m]=atan(mirrorx[i][m]/mirrorz[i][m]);
            nbeta[i][m]=atan(mirrory[i][m]/sqrt(mirrorz[i][m]*mirrorz[i][m]+mirrorx[i][m]*mirrorx[i][m]));
            mirrorz[i][m]= mirrorz[i][m]+CURVATURE;
           //printf("%d %lf %lf %lf\n",i+10,mirrorx[i][m],mirrory[i][m],mirrorz[i][m]);
        }
    }

   
}


void WTelescope::altercoor(double *alfa,double *ox,double *oy,double *nx,double *ny)
{
  *nx=*oy*sin(*alfa)+*ox*cos(*alfa);
  *ny=*oy*cos(*alfa)-*ox*sin(*alfa);
}

void WTelescope::SetMirrorPointError(int flag,double error)
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
    
    for(int i=0; i<7; i++){
        for(int j=0; j<6; j++){
           if(mirrorx[i][j]==-10000) continue;
              x = 0 - mirrorx[i][j];
              y = 0 - mirrory[i][j];
              z = CURVATURE - mirrorz[i][j];
              r = sqrt(x*x+y*y+z*z);
              theta0 = acos(z/r);
              phi0 = atan(y/x);
              if(x>0) phi0 = phi0;
              if(x<0) phi0 = phi0 + TMath::Pi();
              SetEulerMatrix(theta0, phi0);
 
              theta =  gRandom->Rndm();
              theta = sqrt(-log(1-theta)*2*sigma2);
              phi = gRandom->Rndm()*TMath::TwoPi();

              x = CURVATURE*sin(theta)*cos(phi) ;
              y = CURVATURE*sin(theta)*sin(phi) ;
              z = CURVATURE*cos(theta) ;

              InverseEuler( x,  y, z, &xnew, &ynew, &znew);
              xnew += mirrorx[i][j];
              ynew += mirrory[i][j];
              znew += mirrorz[i][j];
              spherecenterx[i][j] = xnew;
              spherecentery[i][j] = ynew;
              spherecenterz[i][j] = znew;
         }
        }
    }

   else{
      printf("The Pointing Errors of the each mirror are not simulated\n");
       for(int i=0; i<7; i++){
          for(int j=0; j<6; j++){
             spherecenterx[i][j] = 0;
             spherecentery[i][j] = 0;
             spherecenterz[i][j] = CURVATURE;
          }
       }
  }

//=== The spherecenterx,y,z, are the sphere center of each small mirror      ===//
////=== The value of the three arrays should be provided by the measurements   ===//
////=== Currently, without the measurments results, all small mirrors can be   ===//
////=== considered to be the same and the sphere center to be (0, 0, CURVATURE)===//
////=== Lingling Ma, 20151014  ===//
               
}

void WTelescope::SetXY(double x, double y)
{
    Telx_ = x;
    Tely_ = y;
}

void WTelescope::SetPointing(double zenith,double azimuth)
{
    TelZ_ = zenith;
    TelA_ = azimuth;
printf("TelZ %f TelA %f\n", TelZ_, TelA_);
}

void WTelescope::SetEulerMatrix(double theta,double phi)
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

double WTelescope::Euler(double x0, double y0, double z0, double *x, double *y, double *z)
{  
   *x = matrix_[0][0]*x0+matrix_[0][1]*y0+matrix_[0][2]*z0;
   *y = matrix_[1][0]*x0+matrix_[1][1]*y0+matrix_[1][2]*z0;
   *z = matrix_[2][0]*x0+matrix_[2][1]*y0+matrix_[2][2]*z0;
}
double WTelescope::Euler(CER_bunch*bunch)
{
    double x,y,z,m,n,l;
    x = matrix_[0][0]*bunch->x_+matrix_[0][1]*bunch->y_+matrix_[0][2]*bunch->z_;
    y = matrix_[1][0]*bunch->x_+matrix_[1][1]*bunch->y_+matrix_[1][2]*bunch->z_;
    z = matrix_[2][0]*bunch->x_+matrix_[2][1]*bunch->y_+matrix_[2][2]*bunch->z_;

    m = matrix_[0][0]*bunch->u_+matrix_[0][1]*bunch->v_+matrix_[0][2]*bunch->l_;
    n = matrix_[1][0]*bunch->u_+matrix_[1][1]*bunch->v_+matrix_[1][2]*bunch->l_;
    l = matrix_[2][0]*bunch->u_+matrix_[2][1]*bunch->v_+matrix_[2][2]*bunch->l_;

    bunch->x_ = x; bunch->y_ = y; bunch->z_ = z;
    bunch->u_ = m; bunch->v_ = n; bunch->l_ = l;

}

double WTelescope::InverseEuler(double x0, double y0, double z0, double *x, double *y, double *z)
{
   *x = matrix_[0][0]*x0+matrix_[1][0]*y0+matrix_[2][0]*z0;
   *y = matrix_[0][1]*x0+matrix_[1][1]*y0+matrix_[2][1]*z0;
   *z = matrix_[0][2]*x0+matrix_[1][2]*y0+matrix_[2][2]*z0;
}

int WTelescope::IncidentTel(double x,double y,double lengthx,double lengthy)
{
    if(fabs(x-Telx_)<lengthx/2 && fabs(y-Tely_)<lengthy/2)
        return 1;
    else return 0;
}
int WTelescope::GetTransmissivity(double wavelength, double angle)
{   
   // * updated by lingling Ma 2019.5.8
   int i;
   double p;
   double transmissivity;
   if(wavelength<filter_wl_min||wavelength>filter_wl_max||angle>filter_angle_max)  transmissivity = 0;
   else{
      transmissivity = Transmissivity-> Interpolate(wavelength, angle);
   } 
   p = gRandom->Rndm();
  if(p<transmissivity) return 1;
  else return 0;
}
 

/*int WTelescope::GetTransmissivity(double wavelength)
{ 
  // * updated by lingling Ma 2018.1.10
  int i;
  double p;
  double transmissivity;
  
  if(wavelength<wl_min||wavelength>wl_max) transmissivity = 0;
  else{
     i = int(wavelength - wl_min);
     //printf("%f\n",Transmissivity[i]);
     transmissivity = Transmissivity[i];
  }

  p = gRandom->Rndm();
  if(p<transmissivity) return 1;
  else return 0; 
}*/
/*
int WTelescope::GetTransmissivity(double wavelength)
{

  // * The filter transmissivity is simulated for photons with wavelength
  // * from 280 nm to 415 nm with step 5 nm.
           
   int i;
   double wl_max = 415; //nm
   double wl_min = 280; //nm
   double wl_step = 5;  //nm
   double wl0, wl1, t0, t1;

   double transmissivity;
   double Transmissivity[28] ={
  0.000600, 0.010656, 0.029956, 0.086233, 0.161489, 0.265211, 0.368200, 0.470700
, 0.552622, 0.629400, 0.679678, 0.726022, 0.756633, 0.781378, 0.797322, 0.806567
, 0.811344, 0.802433, 0.790100, 0.740322, 0.679389, 0.573833, 0.452278, 0.306722
, 0.188967, 0.089744, 0.044167, 0.012000
   };
  if(wavelength<wl_min||wavelength>wl_max) transmissivity = 0;

   else{
      i = int(wavelength-wl_min)/wl_step;
      wl0 = i*wl_step + wl_min;
      wl1 = (i+1)*wl_step + wl_min;
      t0 = Transmissivity[i];
      t1 = Transmissivity[i+1];
      transmissivity = t0 + (t1-t0)*(wavelength-wl0)/(wl1-wl0);
   }

   double p;
   p = gRandom->Rndm();

  // 1 survived, 0 absorbed by the filter
  if(p<transmissivity) return 1;
  else return 0;
 }*/

int WTelescope::Reflected(double wavelength)
{ 
  //Update by Lingling Ma 2019.05.08
  int i = int ((wavelength-mirror_wl_min)/mirror_wl_step);
  int j = i+1;
  double e1, e2, efficiency;
  double wl1, wl2;
  e1 = Reflectivity[i];
  e2 = Reflectivity[j];
  wl1 = mirror_wl_min+i*mirror_wl_step;
  wl2 = mirror_wl_min+j*mirror_wl_step;
  double k = (e1-e2)/(wl1-wl2);
  efficiency = (wavelength-wl2)*k+e2; 
  double x;
  x = gRandom->Rndm();
  if(x<efficiency)
     return 1;
  else
     return 0;
}

int WTelescope::WhichMirror(double x, double y, double z, double *deltax, double *deltay,double *deltaz,int *ii, int *mm)
{   
    double dx, dy, dz;
    double ndx, ndy, ndz;
    double n2dx, n2dy, n2dz;
    double xc,yc,zc;
    double _L = length;
    double theta1, theta2, dtheta;
    
    for(int i=0;i<7;i++){                   //number of the circle
        for(int m=0;m<6;m++){
            if(mirrorx[i][m]==-10000&&mirrory[i][m]==-10000&&mirrorz[i][m]==-10000)
                continue;
            dx = mirrorx[i][m] - x;
            dy = mirrory[i][m] - y;
            dz = mirrorz[i][m] - z;
            altercoor(&nalfa[i][m],&dz,&dx,&ndz,&ndx);
            altercoor(&nbeta[i][m],&dy,&ndz,&n2dy,&n2dz);
            xc = ndx;
            yc = n2dy;
            zc = n2dz;
            if((fabs(yc)<=-fabs(xc)/sqrt3+_L)&&fabs(xc)<=sqrt3/2*_L){
                *deltax = xc;
                *deltay = yc;
                *deltaz = zc;
                *ii = i;
                *mm = m;
                return 1;
            }
        }

    }

}

double WTelescope::Plane(double z,double x0,double y0,double z0,double m2,double n2,double l2,double *x,double *y)
{

 double k;
  k = (z-z0)/l2;
  *x = k*m2+x0;
  *y = k*n2+y0;
}

double WTelescope::Sphere(double Z0,double R, double x0, double y0, double z0, double m, double n, double l, double *x, double *y, double *z)
{
  double A = 1.;
  double B = 2*m*x0+2*n*y0+2*l*z0-2*l*R;
  double C = x0*x0+y0*y0+z0*z0-2*z0*R;
  double delta = B*B-4*A*C;
  double k1, k2, z1, z2;

  if(delta<0){
    *x = -100000;
    *y = -100000;
    *z = -100000;
  }
 else{
    k1 = (-B + sqrt(delta))/(2*A);
    k2 = (-B - sqrt(delta))/(2*A);

    z1 = l*k1 + z0;
    z2 = l*k2 + z0;
if(z1>Z0&&z2>Z0){
      *z = 1000000;
      *x = 1000000;
      *y = 1000000;
    }

    if(z1<Z0&&z2<Z0){
      *z = 1000000;
      *x = 1000000;
      *y = 1000000;
    }
   if(z1<z2) {
      *z = z1;
      *x = m*k1+x0;
      *y = n*k1+y0;
    }
    if(z2<z1){
      *z = z2;
      *x = m*k2+x0;
      *y = n*k2+y0;
    }
  }
}

double WTelescope::GetReflected(double x1,double y1,double z1,double x,double y,double z,double *m2,double *n2,double *l2)
{

 //  x, y, z are the norml directions
 //   //  x1, y1, z1, are the incendent directions
   double costheta = x1*x+y1*y+z1*z;
   *m2 = -(x1-2.*costheta*x);
   *n2 = -(y1-2.*costheta*y);
   *l2 = -(z1-2.*costheta*z);
 
}


double WTelescope::RayTrace(double spotsize,double x0, double y0, double z0, double m1, double n1, double l1,double *xc, double *yc, double *t, double *u, double *v)
{
    double x,y;
    double xmirror0,ymirror0,zmirror0;
    double xmirror,ymirror,zmirror;
    double dist;       ///< RayTrace distance
    double m0, n0, l0; ///< the direction of the normals
    double norm;       ///< to normalize the vector
    double m2, n2, l2; ///< the direction of the reflected ray
    double xcluster,ycluster;  ///< x,y on the pmt cluster plane
    double r2;

    int ii, mm, N;   ///<

    Plane(ZDOOR,x0,y0,z0,m1,n1,l1,&x,&y);
    if(fabs(x)>D_DOOR/2.||fabs(y)>Hdoor/2.) return -1; ///< shelded by the door; 13102015

    Plane(ZCLUSTER1,x0,y0,z0,m1,n1,l1,&x,&y);
    if(fabs(x)<CLUSTER_X/2.&&fabs(y)<CLUSTER_Y/2.) return -2; ///< shelded by the cluster; 13102015

    Plane(ZCLUSTER0,x0,y0,z0,m1,n1,l1,&x,&y);
    if(fabs(x)<CLUSTER_X/2.&&fabs(y)<CLUSTER_Y/2.) return -2;  ///< shelded by the cluster; 13102015

    Sphere(ZMIRROR,CURVATURE,x0,y0,z0,m1,n1,l1,&xmirror0,&ymirror0,&zmirror0);
    if(fabs(xmirror0)>D_DOOR/2.||fabs(ymirror0)>Hdoor/2.) return -3; //shelded by the container's wall

    //* ==== Mirror module function === *///  ///< provided by Ma Lingling 13102015
    double deltax = -10000;
    double deltay = -10000;
    double deltaz = -10000;
    WhichMirror(xmirror0, ymirror0, zmirror0, &deltax, &deltay, &deltaz, &ii, &mm);
    if(deltax==-10000) return -5;  ///< out of mirror area

   //* =============================== *///
    dist = (xmirror0-x0)*(xmirror0-x0) + (ymirror0-y0)*(ymirror0-y0) + (zmirror0-z0)*(zmirror0-z0);
    dist = sqrt(dist);

    xmirror = gRandom->Gaus(xmirror0,spotsize);
    ymirror = gRandom->Gaus(ymirror0,spotsize);
    zmirror = zmirror0;

   //* ==== Mirror module function === *///  ///< provided by Ma Lingling 13102015
    m0 = spherecenterx[ii][mm] - xmirror;
    n0 = spherecentery[ii][mm] - ymirror;
    l0 = spherecenterz[ii][mm] - zmirror;

    //* =============================== *///
    norm = sqrt(m0*m0+n0*n0+l0*l0);
    m0 /= norm;
    n0 /= norm;
    l0 /= norm;

    GetReflected(m1,n1,l1,m0,n0,l0,&m2,&n2,&l2);
    Plane(ZCLUSTER0,xmirror0,ymirror0,zmirror0,m2,n2,l2,&xcluster,&ycluster);

    if(fabs(xcluster)>CLUSTER_X/2.||fabs(ycluster)>CLUSTER_Y/2.) return -4; //out of the range of the cluster


    dist = sqrt((xmirror-xcluster)*(xmirror-xcluster)
         + (ymirror-ycluster)*(ymirror-ycluster)
         + (zmirror-ZCLUSTER0)*(zmirror-ZCLUSTER0))
         - dist;

    *u = m2;
    *v = n2;
    //*ll = l2;
    *xc = -ycluster;
    *yc = xcluster;
    *t = dist/C_LIGHT;

    return 1;

}

//* Added by lingling Ma 2018.1.10
/*void WTelescope::SetTransmissivity()
{
   int wl; 
   double efficiency;
   FILE *fp;
   
   wl_min = 10000;
   wl_max = -1; 
   fp = fopen("Data/filter.txt","r");
   while(!feof(fp)){
      fscanf(fp,"%d %lf\n",&wl, &efficiency);
      if(wl>wl_max) wl_max = wl;
      if(wl<wl_min) wl_min = wl;
      Transmissivity.push_back(efficiency);      
   }
   fclose(fp);
}*/

//* Updated by Lingling Ma 2019.05.08

void WTelescope::SetTransmissivity()
{
  int wl;
  double p1,p2,p3,p4,p5,p6,p7,p8;
  double filter_wl_step = 5;  //nm
  double filter_angle_step = 5; //degree
  FILE *fp;
  filter_wl_min = 300;
  filter_wl_max = 600;
  filter_angle_max = 35;
  filter_angle_min = 0;
  int i=0;
  Transmissivity = new TGraph2D(61*5);
 // Transmissivity->SetName("t");
  fp = fopen("Data/filter_transparency.txt","r");
  while(!feof(fp)){
     fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %lf %lf\n",&wl,&p1,&p2,&p3,&p4,&p5,&p6,&p7,&p8); 
     Transmissivity->SetPoint(i,wl,filter_angle_min+filter_angle_step*0,p1);
     i++;
     Transmissivity->SetPoint(i,wl,filter_angle_min+filter_angle_step*1,p2);
     i++;
     Transmissivity->SetPoint(i,wl,filter_angle_min+filter_angle_step*2,p3); 
     i++;
     Transmissivity->SetPoint(i,wl,filter_angle_min+filter_angle_step*3,p4);
     i++;
     Transmissivity->SetPoint(i,wl,filter_angle_min+filter_angle_step*4,p5);
     i++;
     Transmissivity->SetPoint(i,wl,filter_angle_min+filter_angle_step*5,p6);
     i++;
     Transmissivity->SetPoint(i,wl,filter_angle_min+filter_angle_step*6,p7);
     i++;
     Transmissivity->SetPoint(i,wl,filter_angle_min+filter_angle_step*7,p8);  
     i++;
  }
  fclose(fp);
  
}

void WTelescope::SetReflectivity()
{
  FILE *fp;
  double p;
  int wl;
  fp = fopen("Data/reverseEfficiency.txt","r");
  int i=0;
  mirror_wl_max = 800; 
  mirror_wl_min = 250;
  mirror_wl_step = 10; //nm
  while(!feof(fp)){
     fscanf(fp,"%d %lf\n",&wl,&p);
     Reflectivity[i] = p;
     i++; 
  }
  fclose(fp);  
}
