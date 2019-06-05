#ifndef __WFMirror__
#define __WFMirror__

#include <iostream>
#include <fstream>
#include <string.h>

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <map>
#include "telescopeparameters.h"

#define NCircle 7
#define NMaxMirror 6

using namespace std;

class WFMirror;

/*!
The class for the mirror array of the telescope
*/
class WFMirrorArray{
   public:
   static int NMirror[NCircle];
   WFMirror* pmirr[NCircle][NMaxMirror];

   public:
   void Init();
   void Clear();
   WFMirrorArray() {Init();}
   ~WFMirrorArray() {Clear();}
   void SetMirror();
   void SetMirrorPointError(int flag,double error);
   bool WhichMirror(double x, double y, double z, double *deltax, double *deltay,double *deltaz,int *ii, int *mm);
};

/*!
The class for the mirror of the telescope
*/
class WFMirror
{
   public:
   ///
   static double MirrorSpot;
   ///
   static double MirrorPointError;
   ///
   static double MirrorPointErrorFlag;
   ///
   static double MirrorGeometry;
   ///curvature radius of the mirror,in mm
   static double CURVATURE;
   ///the mirror's reflectivtiy
   static double REFLECTIVITY;
   ///the length of the mirror facet of spherical design;  20140926 liujl,in mm
   static double length;
   ///the length of the mirror facet of Davies design;  20140926 liujl,in mm
   static double length_D;

   ///coordinate of the center of this mirror, z axis is from center of all mirrors to 0
   double mirrorx;
   double mirrory;
   double mirrorz;
   ///in the projected zx plane,need to rotate nalpha along y axis to move from center of all mirrors to center of this mirror (assume center of this mirror located at (x,y,z), then sin(nalpha)=x/sqrt(x*x+z*z),cos(nalpha)=-z/sqrt(x*x+z*z) )
   double nalfa;
   ///in the projected zy plane,need to rotate nbeta along x axis to move from center of all mirrors to center of this mirror (assume center of this mirror located at (x,y,z), then sin(nbeta)=y/sqrt(x*x++y*y+z*z),cos(nbeta)=sqrt(x*x+z*z)/sqrt(x*x+y*y+z*z) )
   double nbeta;
   ///the sphere center of this mirror
   double spherecenterx;
   double spherecentery;
   double spherecenterz;   

   public:
   void Init();
   WFMirror() {Init();}
   ~WFMirror() {;}
   static void SetEulerMatrix(double theta,double phi,double matrix_[3][3]);
   static void InverseEuler(double x0, double y0, double z0, double *x, double *y, double *z,double matrix_[3][3]);
   static void SetMirrorSpot(float mirrorspot);
   static void SetMirrorPointErrorFlag(int point_error_flag);
   static void SetMirrorPointError(float point_error);
   static void SetMirrorGeometry(int mirror_geometry);
   void SetMirrorPointError(int flag,double error);
   ///nx and ny are the new coordinates, while ox and oy are the old coordinates. Compare to the old frame,the new frame rotate alpha along the z axis
   static void altercoor(double *alfa,double *ox,double *oy,double *nx,double *ny);
   ///convert from mirror array coordinate system(z axis is from center of all mirrors to 0) to mirror coordinate system(z axis is from center of one mirror to 0)
   void ConvertCoor(double *ox,double *oy,double *oz,double *nx,double *ny,double *nz);
   bool IsInside(double x,double y,double z,double *deltax,double *deltay,double *deltaz);
   void GetReflected(double x1,double y1,double z1,double x0,double y0,double z0,double *m2,double *n2,double *l2);
};

#endif
