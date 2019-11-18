#ifndef __FluxModel__

#define __FluxModel__
#include "TH1.h"
#define NMaxFlux 92
#define NMaxFluxModel 7
#define NMaxCount 100

/*!
The class used to calculate modeled flux, and also calculate weight for MC
*/
class FluxModel{
   public:
   ///static header for the class
   static FluxModel* _Head;
   ///which flux model will be used
   static int WhichFluxModel;
   ///flux index when generate MC
   static int MCFluxIndex;
   ///debug output
   static int jdebug;
   ///observation time
   double obtime;
   ///stored fluxes(flux numbers,flux,charge),which are loaded from external files
   int nflux[NMaxFluxModel];
   TH1D* fluxes[NMaxFluxModel][NMaxFlux];
   int Zflux[NMaxFluxModel][NMaxFlux];
   double Mflux[NMaxFluxModel][NMaxFlux];
   ///generated MC counts
   int nMCCount;
   ///generated MC energy distribution;
   TH1D* MCCounts[NMaxCount];
   ///MC normalization. MC count will be normalized by the solid angle and area(divided by \deltaS*\deltaOmega)
   double MCNorm[NMaxCount];
   int ZMCCount[NMaxCount];

   public:
   static FluxModel* GetHead();
   void Init();
   void Clear();
   FluxModel() { Init(); }
   ~FluxModel() { Clear(); }
   int LoadFluxFromFile(int model,char* filename,int Z);
   int SetFlux(int model,int Z);
   int SetAllFlux(int model);
   TH1D* GetSum(int model,int Zmin,int Zmax);
   int LoadMCCount(char* filename,int Z);
   void SetMCNorm(int Z,double norm);
   double GetWeight(int Z,double E);

   ///Drawing
   TH1D* Draw(int model,int Z,const char* opt="",double index=0);
   TH1D* DrawSum(int model,int Zmin,int Zmax,const char* opt="",double index=0);
   int Draw(int model,int nz,int* Zlist,double index=0);
   int Draw(int model,int Zmin,int Zmax,double index=0);
};

#endif

