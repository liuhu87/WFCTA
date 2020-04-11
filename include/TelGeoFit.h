#ifndef __TelGeoFit__
#define __TelGeoFit__
#include "common.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "WFCTAEvent.h"

#define MAXEvent 10

class TelGeoFit{
   public:
   static int jdebug;
   int nevent;
   double telpos[NCTMax][3];
   double teldir0[NCTMax][2];
   double rotpos[3];
   double rotdir0[MAXEvent][2];
   double Npe_sipm[MAXEvent][MAXPMT];
   double Npe_sum[MAXEvent];
   double rabbitTime[MAXEvent];
   int telindex[MAXEvent];
   int iEvent[MAXEvent];
   /// for the fitting of the laser image
   ROOT::Math::Minimizer* minimizer;

   public:
   void Init();
   void Clear();
   void SetInfoInit(char* filename=0);
   void Dump();
   TelGeoFit() {Init();}
   ~TelGeoFit() {Clear();}
   int AddEvent(WFCTAEvent* pev,int type=12,bool doclean=true);
   static bool GetImageCoo(double zenith,double azimuth,double dir_in[3],double &xx,double &yy,bool IsLocal=false);
   static void GetOutDir(double zenith,double azimuth,double imagexy[2],double dir_out[3],bool IsLocal=false);
   double Interface(const double* par);
   bool DoFit(int ntel,int* tellist,bool force=false);
};
#endif
