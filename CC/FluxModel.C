#include "FluxModel.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
using namespace std;
FluxModel* FluxModel::_Head=0;
int FluxModel::WhichFluxModel=0;
int FluxModel::MCFluxIndex=-1;
int FluxModel::jdebug=0;
FluxModel* FluxModel::GetHead(){
   if(!_Head){
      _Head=new FluxModel();
   }
   return _Head;
}
void FluxModel::Init(){
   nflux=0;
   for(int ii=0;ii<NMaxFlux;ii++){
      fluxes[ii]=0;
      Zflux[ii]=0;
   }
   nMCCount=0;
   for(int ii=0;ii<NMaxCount;ii++){
      MCCounts[ii]=0;
      ZMCCount[ii]=0;
   }
}
void FluxModel::Clear(){
   for(int ii=0;ii<nflux;ii++){
      if(fluxes[ii]) {delete fluxes[ii]; fluxes[ii]=0;}
   }
   nflux=0;
   for(int ii=0;ii<nMCCount;ii++){
      if(MCCounts[ii]) {delete MCCounts[ii]; MCCounts[ii]=0;}
   }
   nMCCount=0;
}
int FluxModel::LoadFluxFromFile(char* filename,int Z){
   if(!filename) return 0;
   TFile* fin=TFile::Open(filename);
   int nloaded=0;
   if(Z<=0){ ///define Z from name of histogram
   }
   else if(nflux<NMaxFlux){
      fluxes[nflux+1]=0;//...
      Zflux[nflux+1]=Z;
      if(fluxes[nflux+1]) {nflux++; nloaded++;}
   }
   return nloaded;
}
int FluxModel::LoadMCCount(char* filename,int Z){
   if(!filename) return 0;
   TFile* fin=TFile::Open(filename);
   int nloaded=0;
   if(Z<=0){ ///define Z from name of histogram
   }
   else if(nflux<NMaxCount){
      MCCounts[nMCCount+1]=0;//...
      ZMCCount[nMCCount+1]=Z;
      if(MCCounts[nMCCount+1]) {nMCCount++; nloaded++;}
   }
   return nloaded;
}
void FluxModel::SetMCNorm(int Z,double norm){
   for(int ii=0;ii<nMCCount;ii++){
      if(Z!=ZMCCount[ii]) continue;
      MCNorm[ii]=norm;
   }
}
double FluxModel::GetWeight(int Z,double E){
   int index_mc=-1,index_data=-1;
   for(int ii=0;ii<nMCCount;ii++){
      if(Z!=ZMCCount[ii]) continue;
      index_mc=ii;
   }
   for(int ii=0;ii<nflux;ii++){
      if(Z!=Zflux[ii]) continue;
      index_data=ii;
   }
   if(index_mc<0||index_data<0){
      printf("FluxModel:GetWeight: No FLux  or no MC count for Z=%d (index_data=%d,index_mc=%d)\n",Z,index_data,index_mc);
      return 0;
   }
   else{
      int ibin=MCCounts[index_mc]->GetXaxis()->FindBin(E);
      double binwidth=MCCounts[index_mc]->GetXaxis()->GetBinWidth(ibin);
      double MCContent=MCCounts[index_mc]->GetBinContent(ibin);
      if(MCContent<=0||MCNorm[index_mc]<=0){
         if(jdebug>0) printf("FluxModel:GetWeight:No Bin Content or Normalization for MC Z=%d\n",Z);
         return 0;
      }
      double weight=(fluxes[index_data]->GetBinContent(ibin)*obtime)/(MCContent/MCNorm[index_mc]/binwidth);
      if(jdebug>0) printf("FluxModel:GetWeight:Weight for MC Z=%d Calculated\n",Z);
      return weight;
   }
}
