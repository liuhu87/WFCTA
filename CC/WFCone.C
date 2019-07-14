#include <stdlib.h>
#include <stdio.h> 
// boost for read xml file
//#include <boost/property_tree/xml_parser.hpp>
//#include <boost/property_tree/ptree.hpp>

#include <TMath.h>
#include "WFCone.h"
//#include "telescopeparameters.h"
#include <iostream>
#define Deg TMath::DegToRad()

using namespace std;
//using boost::property_tree::ptree;

//.............................................................................
double SquareCone::D_ConeOut=25.4;
double SquareCone::D_ConeIn=24.4;
double SquareCone::CONE_HEIGHT=25.27;
double SquareCone::CONE_ENTRANCE_CIRCLE=12.2;
double SquareCone::CONE_EXIT_CIRCLE=7.5;
double SquareCone::CUTOFF_ANGLE=0.66;
double SquareCone::CONE_REFLECTIVE_EFFICIENCY=0.9;

SquareCone::SquareCone()
{
// begin read xml file
  //ptree pt;
 // read_xml(filename,pt);
  height = CONE_HEIGHT;// = pt.get<float>("WinstonCone.Height");
  exit_circle = CONE_EXIT_CIRCLE;// = pt.get<float>("WinstonCone.Exit");
  entrance_circle = CONE_ENTRANCE_CIRCLE;// = pt.get<float>("WinstonCone.Entrance");
  efficiency = CONE_REFLECTIVE_EFFICIENCY;//=pt.get<float>("WinstonCone.Efficiency");
// end read xml file
  cutoff_angle=asin(exit_circle/entrance_circle);


  P0[0] = .0;
  P0[1] = .0;
  P0[2] = exit_circle/tan(cutoff_angle);

  P1[0] = .0;
  P1[1] = .0;
  P1[2] = height+exit_circle/tan(cutoff_angle);
}

//...........................................................................

void SquareCone::SquareRaytracing()
{
  reflected = true; 
  flagend = 0;
  surfnum = 1;
  iternum = 0;
  strike = false;
  hitpos[0] = hitpos[1] = hitpos[2] = 0.;
  iternum = 0;

  for(int i=0;i<6;i++) DIS[i] = 0;
  double S, C, A,f,beta;
  double ptix,ptiy,ptiz;
  double rayvix,rayviy,rayviz;
  double utemp2[3],ptemp2[3]; //the copy utemp2[] record the initial value of u[] in this turn
  for(int index=0;index<3;index++)
  {
    utemp2[index]=*(utemp+index);
    ptemp2[index]=*(ptemp+index);
  }
  int i,j,k;
  S=sin(-cutoff_angle);
  C=cos(-cutoff_angle);
  f=exit_circle*(1+fabs(S));
  A=1/(4*f);
  double nv[3]={.0,.0,1.0}; // normal direction cosine values of surface 0 & 1 
  double distance[6]={0.0};
  double r[3]={0.0},R[3]={0.0};   //r: vector of new intersection R: the same in the transfer coordinate system
  int mm;
  double dirident, norm_tRay; 
  double tRay_vector[3];
  if(flagend==0)
  {
    for(int surfn=0;surfn<6;surfn++)
    { 
      switch(surfn)
      { 
        case 0:
          intersection(nv,P0); 
          r[0]=*ptemp;
          r[1]=*(ptemp+1);
          r[2]=*(ptemp+2);
          tRay_vector[0]=r[0]-ptemp2[0];
          tRay_vector[1]=r[1]-ptemp2[1];
          tRay_vector[2]=r[2]-ptemp2[2]; 
          dirident=0, norm_tRay=0;
          for( j=0;j<3;j++)
          {
            dirident=dirident+tRay_vector[j]*(*(utemp+j));
            norm_tRay=norm_tRay+tRay_vector[j]*tRay_vector[j];
          }
          norm_tRay=sqrt(norm_tRay);
          if(dirident>0 && norm_tRay>1.0e-4) 
          distance[surfn]=Norm_vector(tRay_vector);
          break;
				          
        case 1:
          intersection(nv,P1);
          r[0]=*ptemp;
          r[1]=*(ptemp+1);
          r[2]=*(ptemp+2);
          tRay_vector[0]=r[0]-ptemp2[0];
          tRay_vector[1]=r[1]-ptemp2[1];
          tRay_vector[2]=r[2]-ptemp2[2]; 
          dirident=0, norm_tRay=0;
          for( j=0;j<3;j++)
          {
            dirident=dirident+tRay_vector[j]*(*(utemp+j));
            norm_tRay=norm_tRay+tRay_vector[j]*tRay_vector[j];
          }
          norm_tRay=sqrt(norm_tRay);
          if(dirident>0 && norm_tRay>1.0e-4) 
          distance[surfn]=Norm_vector(tRay_vector);
          break;
        default:
          //const double PI=3.1415926;
          beta=(surfn-2)*PI/2;
          ptix=*ptemp*cos(beta)+*(ptemp+1)*sin(beta);
          ptiy=*ptemp*(-sin(beta))+*(ptemp+1)*cos(beta);
          ptiz=*(ptemp+2);  
          rayvix=*utemp*cos(beta)+*(utemp+1)*sin(beta); 
          rayviy=*utemp*(-sin(beta))+*(utemp+1)*cos(beta);   
          rayviz=*(utemp+2);
          SegmentCross(ptix,ptiy,ptiz,rayvix,rayviy,rayviz);
          R[0]=*ptemp;
          R[1]=*(ptemp+1);
          R[2]=*(ptemp+2);
          r[0]=R[0]*cos(beta)-R[1]*sin(beta);
          r[1]=R[0]*sin(beta)+R[1]*cos(beta);
          r[2]=R[2];
          tRay_vector[0]=r[0]-ptemp2[0];
          tRay_vector[1]=r[1]-ptemp2[1];
          tRay_vector[2]=r[2]-ptemp2[2]; 
          distance[surfn]=Norm_vector(tRay_vector);
          break;
      } 
      for(mm=0;mm<3;mm++)
        *(ptemp+mm)=ptemp2[mm];
    }
	
    double mini_distance=1.0e6;
    for(i=0;i<6;i++)                                     
    {
      DIS[i]=distance[i];
      if (distance[i]<mini_distance && distance[i]>1.e-5)
      mini_distance=distance[i]; 
    }
    // judge which surface the photon hits
    for(i=0;i<8;i++)
    {
      if(distance[i]==mini_distance) surfnum=i;
    }
    if (iternum>=50)
    {
      flagend=1;
      strike=false;
    } 
    else 
    if(surfnum==0)
    {
      flagend=1;
      strike=true;
      intersection(nv,P0); 
      hitpos[0]=*ptemp,hitpos[1]=*(ptemp+1),hitpos[2]=*(ptemp+2);
    }
    else if(surfnum==1)
    {
      flagend=1;
      strike=false;
    }
    else 
    {
      //const double PI=3.1415926;
      beta=(surfnum-2)*PI/2;
      ptix=*ptemp*cos(beta)+*(ptemp+1)*sin(beta);
      ptiy=*ptemp*(-sin(beta))+*(ptemp+1)*cos(beta);
      ptiz=*(ptemp+2);
      rayvix=*utemp*cos(beta)+*(utemp+1)*sin(beta); 
      rayviy=*utemp*(-sin(beta))+*(utemp+1)*cos(beta);   
      rayviz=*(utemp+2);
      //refresh the elements in p[]
      SegmentCross(ptix,ptiy,ptiz,rayvix,rayviy,rayviz);
      R[0]=*ptemp;
      R[1]=*(ptemp+1);
      R[2]=*(ptemp+2);
      r[0]=R[0]*cos(beta)-R[1]*sin(beta);
      r[1]=R[0]*sin(beta)+R[1]*cos(beta);
      r[2]=R[2];
      p[0]=r[0],p[1]=r[1],p[2]=r[2];
      ptemp[0]=p[0];
      ptemp[1]=p[1];
      ptemp[2]=p[2];
      //refresh the elements in u[]
      ptix=*ptemp*cos(beta)+*(ptemp+1)*sin(beta);
      ptiy=*ptemp*(-sin(beta))+*(ptemp+1)*cos(beta);
      ptiz=*(ptemp+2);
      SegCPCreflect(ptix,ptiy,ptiz,rayvix,rayviy,rayviz);
      u[0]=*utemp*cos(beta)-*(utemp+1)*sin(beta);
      u[1]=*utemp*sin(beta)+*(utemp+1)*cos(beta);
      u[2]=*(utemp+2);
      utemp[0]=u[0];
      utemp[1]=u[1];
      utemp[2]=u[2];
      iternum++;
      Weight *= efficiency;
      // recursion		    
      SquareRaytracing();
    } 		    
  }
}
   //.......................................................................................//
   
void SquareCone::SegmentCross(double ptix,double ptiy,double ptiz,double rayvix, double rayviy,double rayviz)
{
  double S, C, A,f,cutoff_angle;
  cutoff_angle=asin(exit_circle/entrance_circle);
  S=sin(-cutoff_angle);
  C=cos(-cutoff_angle);
  f=exit_circle*(1+fabs(S));
  A=1/(4*f);
  // 
  double m,n,l,x0,y0,z0,x,y,z;
  m=rayvix,n=rayviy,l=rayviz;
  x0=ptix,y0=ptiy,z0=ptiz;
  double c1,c2,c3,t1,t2;
  double tRay_v[3];
  c1=A*C*C*m*m - 2*A*C*S*l*m + A*S*S*l*l;
  c2=2*A*C*C*m*x0 - S*m - C*l + 2*A*S*S*l*z0 - 2*A*C*S*l*x0 - 2*A*C*S*m*z0;
  c3=A*C*C*x0*x0 - C*z0 - S*x0 - exit_circle/S - f + A*S*S*z0*z0 - 2*A*C*S*x0*z0;
  double Delta;
  Delta=c2*c2-4*c1*c3;
  if(Delta>=0)
  {
    t1=(-c2+sqrt(c2*c2-4*c1*c3))/(2*c1);
    t2=(-c2-sqrt(c2*c2-4*c1*c3))/(2*c1);
    double t[2]={t1,t2};
    double distance2=1.0e6;
    double vv[3]={m,n,l};
    for(int i=0;i<2;i++)
    {
      x=t[i]*m+x0;
      y=t[i]*n+y0;
      z=t[i]*l+z0;
      tRay_v[0]=x-x0;
      tRay_v[1]=y-y0;
      tRay_v[2]=z-z0;
      double temp=0;
      for(int j=1;j<3;j++)
        temp=temp+tRay_v[j]*vv[j];
      if(temp>0 && distance2>Norm_vector(tRay_v) && Norm_vector(tRay_v)>1.0e-5)
      {
        *ptemp=x;
        *(ptemp+1)=y;
        *(ptemp+2)=z;
        distance2=Norm_vector(tRay_v);           
      }
    } 
  }
}
 //...........................................................................................
void SquareCone::intersection(double *surfnormal,double *P)
{
  double t,A,B,C,X,Y,Z,a,b,c,x1,y1,z1;
  A=*surfnormal;
  B=*(surfnormal+1);
  C=*(surfnormal+2);
  X=*P;
  Y=*(P+1);
  Z=*(P+2);
  x1=*ptemp;
  y1=*(ptemp+1);
  z1=*(ptemp+2);
  a=*utemp;
  b=*(utemp+1);
  c=*(utemp+2);
  t=(A*X+B*Y+C*Z-A*x1-B*y1-C*z1)/(A*a+B*b+C*c);
  *ptemp=x1+a*t;
  *(ptemp+1)=y1+b*t;
  *(ptemp+2)=z1+c*t;
}
 //...........................................................................................
void SquareCone::SegCPCreflect(double ptix,double ptiy,double ptiz,double rayvix,double rayiy,double rayviz)	 
{
  double S, C, A,f,cutoff_angle;
  cutoff_angle=asin(exit_circle/entrance_circle);
  S=sin(-cutoff_angle);
  C=cos(-cutoff_angle);
  f=exit_circle*(1+fabs(S));
  A=1/(4*f);
  double x,y,z;
  x=ptix;
  y=ptiy;
  z=ptiz;
  double dFx,dFy,dFz,dF[3]; 
  dFx=2*A*C*C*x - S*(2*A*C*z + 1);
  dFy=0;
  dFz=2*A*z*S*S - 2*A*C*x*S - C;
  dF[0]=dFx;
  dF[1]=dFy;
  dF[2]=dFz;
  double nx,ny,nz;
  nx=dFx/Norm_vector(dF);
  ny=dFy/Norm_vector(dF);
  nz=dFz/Norm_vector(dF);
  if(nz<0)
  {
    nx=-nx;
    ny=-ny;
    nz=-nz;
  }
  reflect(nx,ny,nz,rayvix,rayiy,rayviz);        
}
//................................................................................................
void  SquareCone::reflect(double nx,double ny,double nz,double u1,double u2,double u3)
{
  reflected = true;
  double dot;
  dot=nx*u1+ny*u2+nz*u3;

  *utemp=u1-2*dot*nx;
  *(utemp+1)=u2-2*dot*ny;
  *(utemp+2)=u3-2*dot*nz;
  return;
}
//..................................................................................................
double SquareCone::Norm_vector(double *v)
{
  return sqrt((*v)*(*v)+(*(v+1)*(*(v+1)))+(*(v+2)*(*(v+2)))); 
}
//...................................................................................................
int SquareCone::GaussElimination_ColumnSelect(double a[][5],double x[4]) 
{
  int i,j,k,r;double c; 
  double precision=1.0e-6; 
  for(k=1;k<=2;k++) 
  {
    r=k; 
    for(i=k;i<=3;i++)
    if(fabs(a[i][k])>fabs(a[r][k]))r=i; 
    if(fabs(a[r][k])<precision) return -1;
    if(r!=k) 
    {
      for(j=k;j<=4;j++) {c=a[k][j];a[k][j]=a[r][j];a[r][j]=c;}
    } 
    for(i=k+1;i<=3;i++)  
    {
      c=a[i][k]/a[k][k]; 
      for(j=k+1;j<=4;j++) 
        a[i][j]=a[i][j]-c*a[k][j]; 
    } 
  } 
  if(fabs(a[3][3])<precision) return -1;
  for(k=3;k>=1;k--) 
  {
    x[k]=a[k][4]; 
    for(j=k+1;j<=3;j++)  x[k]=x[k]-a[k][j]*x[j]; 
    x[k]=x[k]/a[k][k]; 
  } 
  return 1; 	
}
