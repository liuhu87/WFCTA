#ifndef __WFTelescope__
#define __WFTelescope__

#include <iostream>
#include <fstream>
#include <string.h>
#include "TGraph2D.h"
#include <stdlib.h>
#include <stdio.h>
#include "TGraph.h"
#include "TRandom3.h"
//#include "telescopeparameters.h"

using namespace std;

class WFTelescope;
class WFMirror;
class WFMirrorArray;
class WCamera;
class SquareCone;

/*!
The class for the Wide Field Cerenkov Telescope Array
*/
class WFTelescopeArray{
   private:
   static WFTelescopeArray* _Head;
   public:
   ///random number generator
   static TRandom3* prandom;
   ///control of debug output
   static int jdebug;
   ///number of Cerenkov Telescopes
   static int CTNumber;
   ///wheather do WFCTA simulation
   static bool DoSim;
   ///pointer to each Cerenkov Telescopes
   WFTelescope** pct;

   public:
   static WFTelescopeArray* GetHead(char* filename=0,int seed=0);
   void Init(bool DoNew=true,int seed=0);
   void Clear();
   ///constructor
   WFTelescopeArray(int seed=0) {Init(true,seed);}
   ///constructor
   WFTelescopeArray(char* filename,int seed=0) {Init(false,seed); ReadFromFile(filename);}
   ///deconstructor
   ~WFTelescopeArray() {Clear();}
   void ReadFromFile(char* filename);
   int WhichTel(double x0, double y0);
   int RayTrace(double x0, double y0, double z0, double m1, double n1, double l1,double weight,double wavelength,int &itel,double &t,int &itube,int &icell);
   bool CheckTelescope();
   TGraph* TelView(int iTel);
   WFMirrorArray* GetMirror(int iTel);
   SquareCone* GetCone(int iTel);
   WCamera* GetCamera(int iTel);
};

/*!
The class for one Cerenkov Telescope
*/
class WFTelescope
{
   public:
   ///MirrorSizeX,in cm
   static double MirrorSizeX;
   ///MirrorSizeY,in cm
   static double MirrorSizeY;
   ///information about the door
   ///z coordinate of the door
   static double ZDOOR;
   ///the Width of the door and the maximum size of the mirror,in mm
   static double D_DOOR;
   ///the height of the container;  20140926 liujl,in mm
   static double Hdoor;
   ///the trancparency of the front glass in the door
   static double TRANSPARENCY;

   ///information about the cluster box
   ///x size of pmt array,in mm
   static double CLUSTER_X;
   ///y size of pmt array,in mm
   static double CLUSTER_Y;
   ///the distance between the mirror and the camera,in mm
   static double FOCUS;
   ///minimum z coordinate of cluster box,in mm
   static double ZCLUSTER0;
   ///maximum z coordinate of cluster box,in mm
   static double ZCLUSTER1;

   ///telescope global information
   ///x position of the telescope
   double Telx_;
   ///y position of the telescope
   double Tely_;
   ///zenith angle of the telescope
   double TelZ_;
   ///azimuth angle of the telescope
   double TelA_;
   ///Euler matrix for coordinate transformation
   double matrix_[3][3];
   ///pointer to the mirrors
   WFMirrorArray* pmirr;
   ///pointer to the camera
   WCamera*       pcame;
   ///pointer to the cone
   SquareCone*    pcone;

   ///Filter information
   static TGraph2D *Transmissivity; //[wl][incident angle] wl is from 600nm to 300 nm with 5 nm step. incident angle is from 0 degrees to 35 degrees with 5 degrees step.
   static double filter_wl_max, filter_wl_min;
   static double filter_angle_max, filter_angle_min;
   static double Reflectivity[55]; //mirror reflectivity table
   static TGraph* quantumeff; //quantum efficiency table
   static double mirror_wl_max, mirror_wl_min, mirror_wl_step;

   public:
   void Init();
   void Clear();
   WFTelescope() {Init();}
   ~WFTelescope() {Clear();}
   void SetXY(double x,double y);
   void SetPointing(double zenith, double azimuth);
   void SetEulerMatrix(double theta,double phi); 
   void Euler(double x0, double y0, double z0, double *x, double *y, double *z);
   void InverseEuler(double x0, double y0, double z0, double *x, double *y, double *z);
   bool IncidentTel(double x,double y);
   void ConvertCoor(double x0,double y0,double z0,double m1,double n1,double l1,double &x0new,double &y0new,double &z0new,double &m1new,double &n1new,double &l1new);
   static void SetTransmissivity();
   static int GetTransmissivity(double wavelength, double angle);
   static bool GetTransmissivity(double wavelength);
   static void SetReflectivity();
   static void SetQuantumEff();
   static bool GetQuantumEff(double wavelength);
   static bool Reflected(double wavelength);
   static void Plane(double z,double x0,double y0,double z0,double m2,double n2,double l2,double *x,double *y);
   static void Sphere(double Z0,double R,double x0, double y0, double z0, double m, double n, double l, double *x, double *y, double *z);
   bool RayTraceDoor(double x0, double y0, double z0, double m1, double n1, double l1);
   bool RayTraceCluster(double x0, double y0, double z0, double m1, double n1, double l1,bool IsCluster0,bool IsInside,double &x,double &y);
   bool RayTraceMirror1(double x0, double y0, double z0, double m1, double n1, double l1,double &xmirror0,double &ymirror0, double &zmirror0);
   bool RayTraceMirror2(double xmirror0,double ymirror0,double zmirror0,double m1, double n1, double l1,double &m2, double &n2, double &l2,int &ii,int &mm);
   int RayTraceUpToCone(double x0, double y0, double z0, double m1, double n1, double l1,double &t,double &xcluster,double &ycluster,double &m2,double &n2,double &l2);
   int RayTrace(double x0, double y0, double z0, double m1, double n1, double l1,double wavelength,double &t,int &itube,int &icell);
};

#endif
