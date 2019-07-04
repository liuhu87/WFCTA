#ifndef __Laser_H__
#define __Laser_H__
/*
The class for absorbtion and scatter process calculation
*/
class Atmosphere {
   public:
   ///random number generator
      static TRandom3* prandom;
   ///
      static double aod_air;
      static double aod_aerosol;
      static double scat_air;
      static double scat_aerosol;
      static double distance;
   public:
      void Init();
      Atmosphere() {Init();}
      ~Atmosphere() {;}
      static void SetParameters(char* filename);
      static void RayleighScatter(double theta_in, double phi_in, double wavelength, double &theta_out, double &phi_out);
      static void MieScatter(double theta_in, double phi_in, double wavelength, double &theta_out, double &phi_out);
      static void RayScatterAngle(double wavelength, double &theta, double &phi);
      static void MieScatterAngle(double wavelength, double &theta, double &phi);
      static double FreePathLength();
      static bool IsScattering();
      static double probability();

};

class Laser {
   public:
      static TRandom3* prandom;
      static double intensity;
      static double intensity_err;
      static double wavelength0;
      static double wavelength0_err;
      static double frequency;
      static double spotrange;
      static double  divergence;
   ///position of the laser generator
      double coor[3];
   ///pointing direction of the laser generator
      double dir[2];
      int count_gen;
      double time_gen;
      int ievent_gen;
      double wavelength_gen;
      double coor_gen[3];
      double dir_gen[2];
      double coor_tel[3];
      double dir_tel[2];

   public:
      void Init();
      Laser() {Init();}
      ~Laser() {;}
      double PositionGen();
      double DirectionGen();
      double WaveLengthGen();
      bool LaserGen();
      int EventGen(int ngen=100);
      int Propagate();
};

void probability{
   



}
   

















}
#endif
