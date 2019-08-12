#ifndef __FluxModel__
#define __FluxModel__
#include "TH1.h"
#define NMaxFlux 100
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
   int nflux;
   TH1D* fluxes[NMaxFlux];
   int Zflux[NMaxFlux];
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
   int LoadFluxFromFile(char* filename,int Z);
   int LoadMCCount(char* filename,int Z);
   void SetMCNorm(int Z,double norm);
   double GetWeight(int Z,double E);
};

#endif
