#include <stdlib.h>
#include "LHChain.h"
#include "Laser.h"
#include "TGraphErrors.h"
int main(int argc, char**argv)
{
   // arguments
   if (argc<3) {
      printf("   Usage %s <inname> <seed>\n",argv[0]);
      return 0;
      //exit(0);
   }

   WCamera::SetSiPMMAP();
   const int maxdir=110;
   int ndir=0;
   double lasercoo[maxdir][2];
   double laserdir[maxdir][2];
   double CCphi[maxdir][2];
   double eCCphi[maxdir][3];
   for(int ii=0;ii<maxdir;ii++){
      lasercoo[ii][0]=0;
      lasercoo[ii][1]=0;
      laserdir[ii][0]=0;
      laserdir[ii][1]=0;
      CCphi[ii][0]=0;
      CCphi[ii][1]=0;
      eCCphi[ii][0]=0;
      eCCphi[ii][1]=0;
      eCCphi[ii][2]=0;
   }

   WFTelescopeArray* pct=WFTelescopeArray::GetHead(Form("%s/default.inp",getenv("WFCTADataDir")));
   if(!pct) return 0;
   WFTelescope* pt=pct->pct[0];
   if(!pt) return 0;

   TRandom3* prandom=new TRandom3(atoi(argv[2]));
   TGraphErrors* gr1=new TGraphErrors();
   TGraphErrors* gr2=new TGraphErrors();

   LHChain chain;
   chain.Add(argv[1]);
   int nevent=chain.GetEntries();
   double vpre[3]={-1,-1,-1};
   for(int ientry=0;ientry<nevent;ientry++){
      WFCTAEvent* pev=chain.GetEvent();
      if(!pev) continue;
      if(pev->laserevent.Time<=0){
         printf("This is not laser event\n");
         break;
      }
      double lasertheta=pev->laserevent.LaserDir[0];
      double laserphi=pev->laserevent.LaserDir[1];
      double corepos[3]={pev->laserevent.LaserCoo[0]-pt->Telx_,pev->laserevent.LaserCoo[1]-pt->Tely_,pev->laserevent.LaserCoo[2]-pt->Telz_};
      laserdir[ndir][0]=lasertheta;
      laserdir[ndir][1]=laserphi;
      lasercoo[ndir][0]=corepos[0];
      lasercoo[ndir][1]=corepos[1];
      lasercoo[ndir][2]=corepos[2];

      double maxdist=1.0e6;
      //for(int ii=0;ii<2;ii++){
      //   double idist=fabs(laserdir[ndir][ii]-vpre[ii]);
      //   if(idist<maxdist) maxdist=idist;
      //}
      for(int ii=0;ii<2;ii++){
         double idist=fabs(lasercoo[ndir][ii]-vpre[ii]);
         //printf("entry=%d ii=%d idist=%lf\n",ientry,ii,idist);
         if(idist<maxdist) maxdist=idist;
      }
      //printf("Process entry=%d event=%d ndir=%d lasercoo={%lf,%lf,%lf} laserdir={%lf,%lf} maxdist=%lf vpre={%lf,%lf,%lf}\n",ientry,pev->iEvent,ndir,lasercoo[ndir][0],lasercoo[ndir][1],lasercoo[ndir][2],laserdir[ndir][0],laserdir[ndir][1],maxdist,vpre[0],vpre[1],vpre[2]);
      if(maxdist<0.1) continue;

      bool dofit=pev->DoFit(0,3,true);
      if(!dofit) continue;
      CCphi[ndir][0]=pev->minimizer->X()[2];
      CCphi[ndir][1]=pev->minimizer->X()[3];
      eCCphi[ndir][0]=pev->minimizer->CovMatrix(2,2);
      eCCphi[ndir][1]=pev->minimizer->CovMatrix(2,3);
      eCCphi[ndir][2]=pev->minimizer->CovMatrix(3,3);

      for(int jj=0;jj<1;jj++){
         double theta0=lasertheta;//+prandom->Gaus(0,Laser::LaserZenErr);
         double phi0=laserphi;//+prandom->Gaus(0,Laser::LaserAziErr);
         double dir1[3]={sin(theta0)*cos(phi0),sin(theta0)*sin(phi0),cos(theta0)};
         double dir0[3]={lasercoo[ndir][0],lasercoo[ndir][1],lasercoo[ndir][2]};
         //for(int icoo=0;icoo<3;icoo++) dir0[icoo]+=prandom->Gaus(0,Laser::LaserCooErr);
         double dz=(lasercoo[ndir][2]==0)?0:(-lasercoo[ndir][2]/dir1[2]);
         for(int icoo=0;icoo<3;icoo++) lasercoo[ndir][icoo]+=dz*dir1[icoo];
         double laserphi0=acos(lasercoo[ndir][0]/sqrt(pow(lasercoo[ndir][0],2)+pow(lasercoo[ndir][1],2)));
         if(lasercoo[ndir][1]<0) laserphi0=2*PI-laserphi0;
         laserphi0*=180./PI;

         double planephi,eplanephi,nz,enz;
         int res=pev->GetPlane(planephi,eplanephi,nz,enz,0,3);
         if(res<=0) {printf("GetPlane failed\n"); continue;}

         double telel,etelel,telaz,etelaz;
         //pev->CalTelDir(kk,bb,ekb,zdir,ezdir,telel,etelel,telaz,etelaz);
         res=pev->GetTelDir(telel,etelel,telaz,etelaz);
         if(res<=0) {printf("GetTelDir failed\n"); continue;}
//continue;

         printf("Input: evt=%d TelZ=%lf TelA=%lf Corephi=%lf laserdir={%lf,%lf},  measured: telel={%lf,%lf} telaz={%lf,%lf} planephi={%lf,%lf} CC={%lf,%lf} phi={%lf,%lf}\n",ientry,pt->TelZ_/PI*180,pt->TelA_/PI*180,laserphi0,laserdir[ndir][0],laserdir[ndir][1],telel,etelel,telaz,etelaz,planephi,eplanephi,CCphi[ndir][0]/PI*180,sqrt(eCCphi[ndir][0])/PI*180,CCphi[ndir][1]/PI*180,sqrt(eCCphi[ndir][2])/PI*180);

         gr1->SetPoint(ndir,laserphi0,telel);
         gr1->SetPointError(ndir,0,etelel);
         gr2->SetPoint(ndir,laserphi0,telaz);
         gr2->SetPointError(ndir,0,etelaz);
      }
      for(int ii=0;ii<3;ii++) vpre[ii]=lasercoo[ndir][ii];
      ndir++;
   }

   //printf("ndir=%d\n",ndir);

   TFile* fout=TFile::Open("TelZenAzi.root","RECREATE");
   fout->cd();
   gr1->Write("TelZen");
   gr2->Write("TelAzi");
   fout->Close();
}
