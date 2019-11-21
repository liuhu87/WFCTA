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
   for(int imodel=0;imodel<NMaxFluxModel;imodel++){
   nflux[imodel]=0;
   for(int ii=0;ii<NMaxFlux;ii++){
      fluxes[imodel][ii]=0;
      Zflux[imodel][ii]=0;
      Mflux[imodel][ii]=0;
   }
   }
   nMCCount=0;
   for(int ii=0;ii<NMaxCount;ii++){
      MCCounts[ii]=0;
      ZMCCount[ii]=0;
   }
}
void FluxModel::Clear(){
   for(int imodel=0;imodel<NMaxFluxModel;imodel++){
      for(int ii=0;ii<nflux[imodel];ii++){
         if(fluxes[imodel][ii]) {delete fluxes[imodel][ii]; fluxes[imodel][ii]=0;}
      }
      nflux[imodel]=0;
   }
   for(int ii=0;ii<nMCCount;ii++){


      if(MCCounts[ii]) {delete MCCounts[ii]; MCCounts[ii]=0;}
   }
   nMCCount=0;
}
int FluxModel::LoadFluxFromFile(int model,char* filename,int Z){
   if(model<0||model>=NMaxFluxModel) return 0;
   if(!filename) return 0;
   TFile* fin=TFile::Open(filename);
   int nloaded=0;
   if(Z<=0){ ///define Z from name of histogram
   }
   else if(nflux[model]<NMaxFlux){
      int nn=nflux[model];
      fluxes[model][nn]=0;//...
      Zflux[model][nn]=Z;
      Mflux[model][nn]=Z*2;
      if(fluxes[model][nn]) {nflux[model]=nn+1; nloaded++;}

   }
   return nloaded;
}
int FluxModel::LoadMCCount(char* filename,int Z){
   if(!filename) return 0;
   TFile* fin=TFile::Open(filename);
   int nloaded=0;
   if(Z<=0){ ///define Z from name of histogram
   }
   else if(nflux[WhichFluxModel]<NMaxCount){
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
   for(int ii=0;ii<nflux[WhichFluxModel];ii++){
      if(Z!=Zflux[WhichFluxModel][ii]) continue;
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
      double weight=(fluxes[WhichFluxModel][index_data]->GetBinContent(ibin)*obtime)/(MCContent/MCNorm[index_mc]/binwidth);
      if(jdebug>0) printf("FluxModel:GetWeight:Weight for MC Z=%d Calculated\n",Z);
      return weight;
   }
}







int FluxModel::SetFlux(int model,int Z){
   if(model<0||model>=NMaxFluxModel) return 0;
   if(model<1||model>6) return 0;
   int res=nflux[model];

   int nn=nflux[model];
   if(nn>=NMaxFlux) return 0;

//    if(nn>NMaxFlux) return 0;

   double phiz0[92]={8.73e-2,5.71e-2,2.08e-3,4.74e-4,8.95e-4,1.06e-2,2.35e-3,1.57e-2,3.28e-4,4.60e-3,7.54e-4,8.01e-3,1.15e-3,7.96e-3,2.70e-4,2.29e-3,2.94e-4,8.36e-4,5.36e-4,1.47e-3,3.04e-4,1.14e-3,6.31e-4,1.36e-3,1.35e-3,2.04e-2,7.51e-5,9.96e-4,2.18e-5,1.66e-5,2.75e-6,4.02e-6,9.99e-7,2.11e-6,1.34e-6,1.30e-6,6.93e-7,2.11e-6,7.82e-7,8.42e-7,5.05e-7,7.79e-7,6.98e-8,3.01e-7,3.77e-7,5.10e-7,4.54e-7,6.30e-7,1.61e-7,7.15e-7,2.03e-7,9.10e-7,1.34e-7,5.74e-7,2.79e-7,1.23e-6,1.23e-7,5.10e-7,9.52e-8,4.05e-7,8.30e-8,3.68e-7,1.58e-7,6.99e-7,1.48e-7,6.27e-7,8.36e-8,3.52e-7,1.02e-7,4.15e-7,1.72e-7,3.57e-7,2.16e-7,4.16e-7,3.35e-7,6.42e-7,6.63e-7,1.03e-6,7.70e-7,7.43e-7,4.28e-7,8.06e-7,3.25e-7,3.99e-7,4.08e-8,1.74e-7,1.78e-8,7.54e-8,1.97e-8,8.87e-8,1.71e-8,3.54e-7};
   

double gammaz[92]={-2.71,-2.64,-2.54,-2.75,-2.95,-2.66,-2.72,-2.68,-2.69,-2.64,-2.66,-2.64,-2.66,-2.75,-2.69,-2.55,-2.68,-2.64,-2.65,-2.70,-2.64,-2.61,-2.63,-2.67,-2.46,-2.59,-2.72,-2.51,-2.57,-2.56,-2.55,-2.54,-2.54,-2.53,-2.52,-2.51,-2.51,-2.50,-2.49,-2.48,-2.47,-2.46,-2.46,-2.45,-2.44,-2.43,-2.42,-2.41,-2.40,-2.39,-2.38,-2.37,-2.37,-2.36,-2.35,-2.34,-2.33,-2.32,-2.31,-2.30,-2.29,-2.28,-2.27,-2.25,-2.24,-2.23,-2.22,-2.21,-2.20,-2.19,-2.18,-2.17,-2.16,-2.15,-2.13,-2.12,-2.11,-2.10,-2.09,-2.08,-2.06,-2.05,-2.04,-2.03,-2.02,-2.00,-1.99,-1.98,-1.97,-1.96,-1.94,-1.93};
   double A[92]={1,4,7,9,11,12,14,16,19,20,23,24,27,28,31,32,35,40,39,40,45,48,51,52,55,56,59,59,64,65,70,73,75,79,80,84,85,88,89,91,93,96,99,101,103,106,108,112,115,119,122,128,127,131,133,137,139,140,141,144,145,150,152,157,159,163,165,167,169,173,175,178,181,184,186,190,192,196,197,201,204,207,209,209,210,222,223,225,227,232,231,238};
   double Ep[6]={5.00,3.66,3.50,4.50,3.82,3.68};

   double epselon[6]={1.87,2.30,1.95,1.90,2.32,1.84};
   double gammac[3]={-4.48,-7.80,-3.06};
   double deltagamma[3]={2.20,5.70,0.45};

   for(int iz=1;iz<=92;iz++){
      if((Z>=0)&&(iz!=Z)) continue;
      Zflux[model][nn]=iz;
      Mflux[model][nn]=A[iz-1];
     // fluxes[model][nn]=new TH1D(Form("flux_Z%d_Model%d",iz,model),";log10(E/(GeV/nucleus));Flux dPHI/dE",1000,log10(30.),10.);
     fluxes[model][nn]=new TH1D(Form("flux_Z%d_Model%d",iz,model),";log10(E/(GeV/nucleus));Flux E^2.5dPHI/dE",1000,log10(30.),10.);
     for(int ibin=1;ibin<=fluxes[model][nn]->GetNbinsX();ibin++){
       double xx=pow(10,fluxes[model][nn]->GetXaxis()->GetBinCenter(ibin));        


        // double Ez=model==1?Zflux[model][nn]*Ep[model-1]:(model==2?Mflux[model][nn]*Ep[model-1]:Ep[model-1]);
       double Ez=model==1?Zflux[model][nn]*Ep[model-1]:Mflux[model][nn]*Ep[model-1];


         Ez*=1.0e6;
         if(model<=3){
            fluxes[model][nn]->SetBinContent(ibin,pow(xx,gammaz[iz-1])*pow((1+pow(xx/Ez,epselon[model-1])),(gammac[model-1]-gammaz[iz-1])/epselon[model-1]));
         }
         else{
            fluxes[model][nn]->SetBinContent(ibin,pow(xx,gammaz[iz-1])*pow((1+pow(xx/Ez,epselon[model-1])),-deltagamma[model-4]/epselon[model-1]));
          }
      }
      int ibin0=fluxes[model][nn]->GetXaxis()->FindBin(log10(1000.));
    fluxes[model][nn]->Scale(phiz0[Zflux[model][nn]-1]*1.0e-3/fluxes[model][nn]->GetBinContent(ibin0));
      nn++;
      if(nn>=NMaxFlux) break;
   }
   nflux[model]=nn;
   return nn-res;
}








int FluxModel::SetAllFlux(int model){
   int res=0;
   for(int imodel=1;imodel<=NMaxFluxModel;imodel++){
      if(model>=1&&(imodel!=model)) continue;
      for(int iz=0;iz<NMaxFlux;iz++){
         res+=SetFlux(imodel,iz+1);
      }
   }
   return res;
}







 TH1D* FluxModel::GetSum(int model,int Zmin,int Zmax){
 
   if(model<0||model>=NMaxFluxModel) return 0;
   int index=-1;
   for(int iz=0;iz<NMaxFlux;iz++){
      if(fluxes[model][iz]&&(Zflux[model][iz]>=Zmin&&Zflux[model][iz]<=Zmax)){
         index=iz; break;
      }
   }
   if(index<0) return 0;
   TH1D* hist=(TH1D*)fluxes[model][index]->Clone(Form("model%d_sum_%d_%d",model,Zmin,Zmax));
   hist->Reset();
   for(int iz=0;iz<NMaxFlux;iz++){
      if(fluxes[model][iz]&&(Zflux[model][iz]>=Zmin&&Zflux[model][iz]<=Zmax)){
         hist->Add(fluxes[model][iz]);
      }
   }
   return hist;
}








//TH1D* FluxModel::Draw(int model,int Z,double index){
TH1D* FluxModel::Draw(int model,int Z,const char* opt,double index){
   if(model<0||model>=NMaxFluxModel) return 0;
   for(int iz=1;iz<=nflux[model];iz++){
      if(Z>=0&&(iz!=Z)) continue;
      TH1D* hh=fluxes[model][iz-1];
      if(index!=0){
         hh=(TH1D*)fluxes[model][iz-1]->Clone(Form("%s_index%lf",fluxes[model][iz-1]->GetName(),index));
    for(int ibin=1;ibin<hh->GetNbinsX();ibin++) hh->SetBinContent(ibin,hh->GetBinContent(ibin)*pow(10,hh->GetXaxis()->GetBinCenter(ibin)*index));
      hh->GetYaxis()->SetTitle(Form("Flux*E^{%.1lf} [m^2 sr s GeV/n^{%.1lf}]",index,index-1));
      }
      if(Z>=0) hh->Draw(opt);
      if(Z>=0) return hh;
   }
   return 0;
}







//Draw Flux of the Sum from Zmin to Zmax
TH1D* FluxModel::DrawSum(int model,int Zmin,int Zmax,const char* opt,double index){
   TH1D* hist=GetSum(model,Zmin,Zmax);
   if(!hist) return 0;
   TH1D* hh=hist;
   if(index!=0){
      for(int ibin=1;ibin<hh->GetNbinsX();ibin++) hh->SetBinContent(ibin,hh->GetBinContent(ibin)*pow(10,hh->GetXaxis()->GetBinCenter(ibin)*index));
     // for(int ibin=1;ibin<hh->GetNbinsX();ibin++) hh->SetBinContent(ibin,hh->GetBinContent(ibin)*(hh->GetXaxis()->GetBinCenter(ibin)*index));   

      hh->GetYaxis()->SetTitle(Form("Flux*E^{%.1lf} [m^2 sr s GeV/n^{%.1lf}]",index,index-1));
   }
   hh->Draw();
   return hh;
}







int FluxModel::Draw(int model,int nz,int* Zlist,double index){
   if(model<0||model>=NMaxFluxModel) return 0;
   int nplot=0;
   for(int iz=0;iz<nz;iz++){
      int index=-1;
      for(int ii=0;ii<nflux[model];ii++){
         if(Zlist[iz]==Zflux[model][ii]) {index=ii; break;}
      }
      if(index<0) continue;
      TH1D* hh=Draw(model,Zlist[iz],nplot==0?"hist":"hist same",index);
      if(hh) nplot++;
   }
   return nplot;
}









int FluxModel::Draw(int model,int Zmin,int Zmax,double index){
   int nz=Zmax-Zmin+1;
   int Zlist[100];
   for(int iz=0;iz<nz;iz++){
      Zlist[iz]=Zmin+iz;
   }
   return Draw(model,nz,Zlist,index);
}
