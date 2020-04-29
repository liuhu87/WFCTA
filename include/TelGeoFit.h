#ifndef __TelGeoFit__
#define __TelGeoFit__
#include "common.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "WFCTAEvent.h"
#include "TCanvas.h"

#define MAXEvent 200

class TelGeoFit{
   public:
   static int jdebug;
   static int FixRotDir[4];
   static int FixRotPos[3];
   static int FixRotTime;
   static int FixTelDir[2];
   static int FixTime;
   static int FixLi;
   static int FitLineWidth;
   static double InfNumber;
   static bool GoDown;
   int nevent;
   double telpos[NCTMax][3];
   double teldir0[NCTMax][2];
   double rotpos[NRotMax][3];
   double rotdir0[MAXEvent][2];
   double Npe_sipm[MAXEvent][MAXPMT];
   double Npe_sum[MAXEvent];
   double Time_sipm[MAXEvent][MAXPMT];
   double rabbitTime[MAXEvent];
   int telindex[MAXEvent];
   int iEvent[MAXEvent];
   double image_pars[MAXEvent][2][2];
   bool replaced[MAXEvent];
   /// for the fitting of the laser image
   ROOT::Math::Minimizer* minimizer;
   double XYchi2_evt[MAXEvent];
   int XYndof_evt[MAXEvent];
   double Tchi2_evt[MAXEvent];
   int Tndof_evt[MAXEvent];
   double Chi2_All;
   int Ndof_All;
   ///Fit Result
   double RotEle_int[NRotMax];
   double RotAzi_int[NRotMax];
   double RefEle_int[NRotMax];
   double RefAzi_int[NRotMax];
   double RotPos_int[NRotMax][3];
   double RotTime_int[NRotMax];
   double TelZen_int[NCTMax];
   double TelAzi_int[NCTMax];
   double RotEle_fit[NRotMax];
   double RotAzi_fit[NRotMax];
   double RefEle_fit[NRotMax];
   double RefAzi_fit[NRotMax];
   double RotPos_fit[NRotMax][3];
   double RotTime_fit[NRotMax];
   double TelZen_fit[NCTMax];
   double TelAzi_fit[NCTMax];

   public:
   void Init();
   void Clear();
   void SetInfoInit(char* filename=0);
   void Dump();
   TelGeoFit() {Init();}
   ~TelGeoFit() {Clear();}
   int AddEvent(WFCTAEvent* pev,int type=12,bool doclean=true);
   ///Geometry about WFCTA
   ///Convertion between coordinate at the focus plane and out direction
   static bool GetImageCoo(double zenith,double azimuth,double dir_in[3],double &xx,double &yy,bool IsLocal=false);
   static void GetOutDir(double zenith,double azimuth,double imagexy[2],double dir_out[3],bool IsLocal=false);

   //From long axis coordinate to x&y coor. icoo==0:x; icoo==1:y
   static double GetImageCoo(double CC,double phi,double longcoo,int icoo);
   ///From long axis coordinate to out direction
   static void GetOutDir(double zenith,double azimuth,double CC,double phi,double longcoo,double dir_out[3]);

   ///From out direction to image at the focus plane
   static void Getnz(double nz,double &nz1,double &nz2);
   static bool Getnz(double telcoo[3],double incoo[3],double indir[2],double &planephi,double &nz,double telzcoo[3]);
   static double GetApar(double &Anz,double &phi,double zenith,double azimuth,double planephi,double nz);
   static void GetCCphi(double zenith,double azimuth,double planephi,double nz,double &CC,double &phi);
   static void GetCCphi(double zenith,double azimuth,double telcoo[3],double incoo[3],double indir[2],double &CC,double &phi);
   ///From image to out direction. sign is true when the time is increasing when move along the long axis
   static void CalPlane(double zenith,double azimuth,double CC,double phi,double &planephi,double &nz,int &signnz,bool sign);
   ///From image to out direction. refplanephi is the reference planephi
   static void CalPlane(double zenith,double azimuth,double CC,double phi,double &planephi,double &nz,int &signnz,double refplanephi);
   static void CalPlane(double planephi,double nz,double xyzdir[3][3]);

   ///From long axis coordinate to out angle
   static double GetOutAngle(double zenith,double azimuth,double CC,double phi,double longcoo,bool sign);
   ///From long axis coordinate to out angle
   static double GetOutAngle(double zenith,double azimuth,double CC,double phi,double longcoo,double refplanephi);
   ///From long axis coordinate to out angle
   static double GetOutAngle(double zenith,double azimuth,double planephi,double nz,double longcoo);
   ///From out angle to long axis coordinate
   static double GetLongCoo(double zenith,double azimuth,double planephi,double nz,double OutAngle);

   ///Geometry Outside Telescope
   static double GetInjAngle(double telcoo[3],double incoo[3],double indir[2],double telzcoo[3]);
   static double GetOutAngleFromLength(double telcoo[3],double incoo[3],double indir[2],double length);
   static double GetOutAngleFromHeight(double telcoo[3],double incoo[3],double indir[2],double height);
   static double GetLongCooFromLength(double zenith,double azimuth,double telcoo[3],double incoo[3],double indir[2],double length);
   static double GetLongCooFromHeight(double zenith,double azimuth,double telcoo[3],double incoo[3],double indir[2],double height);
   static double GetScatAngleFromLength(double telcoo[3],double incoo[3],double indir[2],double length);
   static double GetScatAngleFromHeight(double telcoo[3],double incoo[3],double indir[2],double height);
   static double GetScatAngleFromLongcoo(double zenith,double azimuth,double telcoo[3],double incoo[3],double indir[2],double longcoo);
   static double GetGroundLength(double telcoo[3],double incoo[3],double indir[2]);
   static double GetRp(double telcoo[3],double incoo[3],double indir[2],double &length);

   ///Time about WFCTA
   static double GetTime(double outangle,double injangle,double groundlength,double delta_length);
   ///based on the image info
   static double GetTime(double zenith,double azimuth,double CC,double phi,double longcoo,double groundlength,double delta_length,double injangle,bool sign);
   ///based on the outside telescope info and long axis coordinate
   static double GetTimeFromLongcoo(double zenith,double azimuth,double telcoo[3],double incoo[3],double indir[2],double longcoo);
   ///based on the outside telescope info
   static double GetTimeFromLength(double zenith,double azimuth,double telcoo[3],double incoo[3],double indir[2],double length,double &longcoo);
   static double GetTimeFromHeight(double zenith,double azimuth,double telcoo[3],double incoo[3],double indir[2],double height,double &longcoo);

   ///Geometry about Rotate
   static double GetRotPos(int icoo,int Li);
   static double GetRotDir(int ipar,int Li);
   static double GetRotTime(int Li);
   static bool CalRotateZeroPos(double ele_rotate,double azi_rotate,double ele_ref,double azi_ref,double &ele0,double &azi0);
   static bool CalDir_out(double ele_rotate,double azi_rotate,double ele_ref,double azi_ref,double ele_in,double azi_in,double &ele_out,double &azi_out);
   static bool CalDir_out(double ele_in,double azi_in,int Li,double &ele_out,double &azi_out);
   int GetRotIndex(int ievt);
   double GetChi2XY(int &ndof,int ievt,double CC,double phi);
   double GetChi2Time(int &ndof,int ievt,double phi,double time0,double zenith,double azimuth,double telcoo[3],double incoo[3],double indir[2]);
   double Interface(const double* par);
   void SetRotTelPars(bool IsInit=true,bool IsFit=false);
   void SetRotPars(int Li,const double* par,bool IsInit=true,bool IsFit=false);
   void SetTelPars(int iTel,const double* par,bool IsInit=true,bool IsFit=false);
   bool DoFit(int ntel,int* tellist,int nrot,int* rotlist,bool force=false);
   bool FitProcedure1(int ntel,int* tellist,int nrot,int* rotlist);
   bool FitProcedure2(int ntel,int* tellist,int nrot,int* rotlist);
   void DumpFit();
   TCanvas* DrawFitXY(int iTel=-1,int Li=-1,int iEle=-1,int Time=-1,double plotrange=8.3);
   TCanvas* DrawFitTime(int iTel=-1,int Li=-1,int iEle=-1,int Time=-1);
   TCanvas* DrawNpeScatA(int iTel=-1,int Li=-1,int iEle=-1,int Time=-1);
};
#endif
