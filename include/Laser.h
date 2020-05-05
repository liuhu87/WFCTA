#ifndef __Laser_H__
#define __Laser_H__
#include "common.h"
#include <TCanvas.h>
#include <TView3D.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
//#include <TSystem.h>
#include <TObjArray.h>
#include "TH1F.h"
#include "WFTelescope.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TFile.h"

class WFCTAEvent;

#define MaxATMModel 5
#define MaxATMLayer 5

/*
The class for absorbtion and scatter process calculation
*/
class Atmosphere {
   private:
      static Atmosphere* _Head;
   public:
      ///extinction coefficient
      static double aod_air;
      static double aod_aerosol;
      static double scat_air;
      static double scat_aerosol;
      static TGraph* gRayScatAngle;
      static TGraph* gMieScatAngle;
      static double scale;
      static int ATMRayModel;
      static int ATMMieModel;
      static double xr;
      ///parameters about atmosphere
      int nmodel[2];
      int nlayer[MaxATMModel];
      double ai[MaxATMModel][MaxATMLayer];
      double bi[MaxATMModel][MaxATMLayer];
      double ci[MaxATMModel][MaxATMLayer];
      double layerboun[MaxATMModel][MaxATMLayer];
      ///parameters about aerosol
      double mixlayer[MaxATMModel][2];
      double mie_atten_length[MaxATMModel];
      double mie_scale_height[MaxATMModel];
   public:
      void Init(int seed=0);
      void Release();
      Atmosphere(int seed=0) {Init(seed);}
      ~Atmosphere() {Release();}
      static Atmosphere* GetHead(int seed=0);
      static void SetParameters(char* filename=0);
      void AddATMModel(char* filename);
      void DumpATMModel(int whichmodel=-1);
      double GetRayMaxGrammage();
      double GetMieMaxAbs();
      double GetRayGrammage(double z);
      double GetRayGrammage(double length,double z0,double zenith);
      double GetMieAbs(double z);
      double GetMieAbs(double length,double z0,double zenith);
      double GetRayDensity(double z);
      double GetMieCoeff(double z);
      double GetRayZFromGrammage(double grammage);
      double GetRayZFromGrammage(double grammage,double z0,double zenith);
      double GetMieZFromAbs(double Mie_Abs);
      double GetMieZFromAbs(double Mie_Abs,double z0,double zenith);
      double GetRayLengthFromGrammage(double grammage,double z0,double zenith);
      double GetMieLengthFromAbs(double Mie_Abs,double z0,double zenith);
      static double ZDependence(double z,int type=0);
      static double DeltaZ(double z);
      static double ProbTransform(double xx,double yy[2],double &weight,bool IsCenter);
      static bool RayScatterAngleTheta(double wavelength, double &theta, double anglerange[2],double &weight);
      static bool RayScatterAngleTheta(double wavelength, double &theta, int ntel, int* telindex, double anglerange[NCTMax][2],double &weight,int &whichtel);
      static bool MieScatterAngleTheta(double wavelength, double &theta, double anglerange[2],double &weight);
      static bool MieScatterAngleTheta(double wavelength, double &theta, int ntel, int* telindex, double anglerange[NCTMax][2],double &weight,int &whichtel);
      static bool TestScatterAngleTheta(double wavelength, double &theta, double anglerange[2],double &weight);
      static bool TestScatterAngleTheta(double wavelength, double &theta, int ntel, int* telindex, double anglerange[NCTMax][2],double &weight,int &whichtel);
      static bool UniformScatterAnglePhi(double wavelength, double &phi,double anglerange[2],double &weight);
      static bool UniformScatterAnglePhi(double wavelength, double &phi, int ntel, int* telindex, double anglerange[NCTMax][2],double &weight,int &whichtel);
      static double FreeIntgLength(double lengthrange[2],double &weight);
      static double FreePathLength(double z0,double dir0[3],double lengthrange[2],double &weight);
      double FreePathLength(double z0,double dir0[3],double lengthrange[2],double &weight,double lamda);
      double FreePathLength(double z0,double dir0[3],int ntel,int* telindex,double lengthrange[NCTMax][2],double &weight,double lamda,int &whichtel);
      static int IsScattering(double z0);
      int IsScattering(double z0,double lamda);
};

class Laser {
   public:
      static bool UseTestScat;
      static int WhichRot;
      static int WhichTel;
      static int jdebug;
      static int Doigen;
      static bool DoPlot;
      static bool Doextin;
      static float IniRange[4][2];
      static TRandom3* prandom;
      static double TelSimDist;
      static double TelSimAngl;
      static double scale;
      static double unittime;
      static double intensity;
      static double intensity_err;
      static double wavelength0;
      static double wavelength0_err;
      static double frequency;
      static double pulsetime;
      static double spotrange[2];
      static double  divergence[2];

      static double LaserCooErr;
      static double LaserZenErr;
      static double LaserAziErr;

      ///calculate the probability
      static double lengthmin;
      static double lengthmax;
      static double lengthmin2;
      static double lengthmax2;
      //static TH1D* hdenu;
      //static TH1D* hprob[NCTMax][MAXPMT];
      //static TH1D* hleng[NCTMax][MAXPMT];
      static TH1D* hlength;
      static TH1D* htheta;
      static TH1D* hphi;

      ///position of lhaaso
      static double lhaaso_coo[3];

      ///position of the laser generator
      double lasercoo[3];	//in cm
      ///pointing direction of the laser generator
      double laserdir[2];
      double count_gen;
      int Time_gen;
      double time_gen;
      int ievent_gen;
      int iphoton_gen;
      double wavelength_gen;	//in nm
      double coor_gen[3];	//in cm
      double dir_gen[3];
      vector<double> vgwav;
      vector<double> vgcoo[3];
      vector<double> vgdir[3];

      ///reweight during propagation in the atmosphere
      vector<double> vowei;

      ///coordinate and direction at telescope plane
      int Telindex;
      double coor_out[3];
      double dir_out[3];
      double interpoint;
      double theta_out;
      double phi_out;
      vector<double> votim;
      vector<int> votel;
      vector<double> vocoo[3];
      vector<double> vodir[3];

      ///tel info
      vector<int> tellist;

      ///info about the scattering
      vector<double> votheta;
      vector<double> vophi;
      ///interaction point
      vector<double> volength;
      vector<double> volength2;
      vector<int> vosipm;

      WFCTAEvent* pwfc; //!
      TObjArray* plot; //!
      float plotrange[4][2]; //!

   public:
      void Init(int seed=0);
      void Release();
      void Reset();
      Laser(int seed=0) {Init(seed);}
      ~Laser() { Release(); }
      void SetParameters(char* filename=0);
      static void cross(double dir1[3],double dir2[3],double *dir3);
      static bool CartesianFrame(double zero[3],double coor_in[3],double dir_in[3],double dir_in2[3],double *xdir,double *ydir,double *zdir);
      static double mindist(double zero[3],double coor_in[3],double dir_in[3],double *coor_min,bool &decrease);
      static double mindist(double coor_in[3],double dir_in[3],int &whichtel,double *coor_min,bool &decrease);
      static void PositionDis(double &xx,double &yy);
      static void DirectionDis(double &theta,double &phi);
      double WaveLengthGen();
      void GetAveTelPos(double zero[3]);
      bool InitialGen();
      long int EventGen(int &Time,double &time,bool SimPulse=false);
      int FindAllRange(double zero[3],double cooout[3],double dirout[3],double dirin[3],double allrange[2],int type,double length=0,double theta_scat=0,double phi_scat=0);
      int FindLengthRange(double zero[3],double cooout[3],double dirout[3],double dirin[3],double lengthrange[2]);
      int FindLengthRange(double cooout[3],double dirout[3],int* telindex,double lengthrange[NCTMax][2],int ntel=-1);
      int FindThetaRange(double zero[3],double cooout[3],double dirout[3],double dirin[3],double thetarange[2],double freelength);
      int FindThetaRange(double cooout[3],double dirout[3],int* telindex,double thetarange[NCTMax][2],double freelength,int ntel=-1);
      int FindPhiRange(double zero[3],double cooout[3],double dirout[3],double dirin[3],double phirange[2],double freelength,double theta_scat);
      int FindPhiRange(double cooout[3],double dirout[3],int* telindex,double phirange[NCTMax][2],double freelength,double theta_scat,int ntel=-1);
      int FindWhichTel(double cooout[3],double dirout[3],double freelength,double theta_scat,double phi_scat,int ntel=-1,int* telindex=0);
      int Propagate(double &distance,double &weight);
      bool DoWFCTASim(int SimTel=-1);
      void Add(double coorin1[4],double coorin2[4],int type,double weight);
      static int GetLineStyle(int type,double weight);
      static int GetLineWidth(int type,double weight);
      static int GetLineColor(int type,double weight);
      TCanvas* Draw(const char* option="al",int ViewOpt=0,const char* savedir=0);
      int GetProb(long int ngen=1000000000);
};

#endif
