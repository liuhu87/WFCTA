/*******************************************************************************


This is the simulation program of square Winston cone developed by Chong Wang

Update to more c++ by Baiyang Bi at Jue. 9th, 2017


******************************************************************************/
#ifndef SQUARE_WINSTONCONEF_H
#define SQUARE_WINSTONCONEF_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
//#define PI 3.1415926535898
using namespace std;
class SquareCone
{
  public:
    ///the outer Diameter of the WistonCone,in mm
    static double D_ConeOut;
    ///the inner Diameter of the WistonCone,in mm
    static double D_ConeIn;
    static double CONE_HEIGHT;
    static double CONE_ENTRANCE_CIRCLE;
    static double CONE_EXIT_CIRCLE;
    static double CUTOFF_ANGLE;
    static double CONE_REFLECTIVE_EFFICIENCY;

    // Constructor & Destructor                       ---Baiyang Bi 

    SquareCone();
    ~SquareCone() {;}

    // inline Setters & Getters
    inline void SetInitDir(double dir[3]) {for(int i=0;i<3;i++) u[i]=utemp[i]=dir[i];}
    inline void SetInitPos(double x,double y){p[0] = ptemp[0]= x, p[1] = ptemp[1] = y; p[2] = ptemp[2] = P1[2];Weight=1;}

    void GetHitPos(double* const pos, double &wet) {for(int i=0;i<3;i++) pos[i] = hitpos[i];wet=Weight;};
    int GetIternum() const {return iternum;};
    bool GetStrike() const {return strike;};

    //* ==== Function of winston cone ray tracing with recursive algorithm ==== *//
    void SquareRaytracing();
  private: 
    double height, exit_circle, entrance_circle, cutoff_angle, Weight;
    // two const Points on entrace aperture plane and exit aperture plane
    double P0[3],P1[3];
    // photon's contemporary direction cosine
    double u[3],utemp[3]; //utemp[] is a copy of u[] ,  
    // initial intersection before the next raytracing in Winst cone 
    double p[3],ptemp[3]; //ptemp[] is a copy of p[]
    int flagend,surfnum,iternum; // flag of finish, contemporary cross suface number and the times of iteration(intial value=1)
    // return values
    bool strike;               
    double hitpos[3];   //  hit position on local Si_PM coordinations	
    //test distances for each surface
    double DIS[8];
    bool reflected;
    double efficiency;

    // to calculate the intersection with curved reflecting surface
    void SegmentCross(double ptix,double ptiy,double ptiz,double rayvix, double rayviy,double rayviz);
    // to calculate the intersection with entrace or exit aperture plane
    void intersection(double *surfnormal,double *P);
    // to calculate the direction cosine of reflecting ray   
    void SegCPCreflect(double ptix,double ptiy,double ptiz,double rayvix,double rayiy,double rayviz);
    void reflect(double nx,double ny,double nz,double m,double n,double l);                  
    double Norm_vector(double *v);
    int GaussElimination_ColumnSelect(double a[][5],double x[4]) ;
};

#endif // SQUARE_WINSTONCONEF_H
