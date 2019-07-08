#ifndef __Laser_H__
#define __Laser_H__
#include "common.h"
#include "WFTelescope.h"
#include "TRandom3.h"

class WFCTAEvent;

/*
The class for absorbtion and scatter process calculation
*/
class Atmosphere {
   public:
      ///random number generator
      static TRandom3* prandom;
      ///extinction coefficient
      static double aod_air;
      static double aod_aerosol;
      static double scat_air;
      static double scat_aerosol;
   public:
      void Init(int seed=0);
      void Release();
      Atmosphere(int seed=0) {Init(seed);}
      ~Atmosphere() {Release();}
      static void SetParameters(char* filename=0);
      static double ZDependence(double z,int type=0);
      static double DeltaZ(double z);
      static void RayScatterAngle(double wavelength, double &theta, double &phi);
      static void MieScatterAngle(double wavelength, double &theta, double &phi);
      static double FreeIntgLength();
      static double FreePathLength(double z0,double dir0[3]);
      static int IsScattering(double z0);

};

class Laser {
   public:
      static TRandom3* prandom;
      static double TelSimDist;
      static double unittime;
      static double intensity;
      static double intensity_err;
      static double wavelength0;
      static double wavelength0_err;
      static double frequency;
      static double pulsetime;
      static double spotrange;
      static double  divergence;

      ///position of the laser generator
      double lasercoo[3];	//in cm
      ///pointing direction of the laser generator
      double laserdir[2];
      int count_gen;
      int Time_gen;
      double time_gen;
      int ievent_gen;
      double wavelength_gen;	//in nm
      double coor_gen[3];	//in cm
      double dir_gen[3];
      vector<double> vgwav;
      vector<double> vgcoo[3];
      vector<double> vgdir[3];

      ///coordinate and direction at telescope plane
      int Telindex;
      double coor_out[3];
      double dir_out[3];
      vector<int> votel;
      vector<double> vocoo[3];
      vector<double> vodir[3];

      WFCTAEvent* pwfc; //!

   public:
      void Init(int seed=0);
      void Release();
      void Reset();
      Laser(int seed=0) {Init(seed);}
      ~Laser() { Release(); }
      static void cross(double dir1[3],double dir2[3],double dir3[3]);
      static bool CartesianFrame(double zero[3],double coor_in[3],double dir_in[3],double xdir[3],double ydir[3],double zdir[3]);
      static double mindist(double zero[3],double coor_in[3],double dir_in[3],double coor_min[3],bool &decrease);
      static double mindist(double coor_in[3],double dir_in[3],int &whichtel,double coor_min[3],bool &decrease);
      static void PositionDis(double &xx,double &yy);
      static void DirectionDis(double &theta,double &phi);
      double WaveLengthGen();
      bool InitialGen();
      int EventGen(int &Time,double &time);
      int Propagate();
      bool DoWFCTASim();
};

#endif
