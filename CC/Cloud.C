#include "Cloud.h"
#include "common.h"
#include "math.h"
#include "TMath.h"
#include <fstream>
#include "WFTelescope.h"
int Cloud::drawcircle=9;
double Cloud::Cbintheta[Cntheta+1];
int Cloud::Cnbinphi[Cntheta];
double Cloud::Cbinphi[Cntheta][Cnphi+1];
void Cloud::SetBins(){
   //Set Theta Bins
   for(int ii=0;ii<=Cntheta;ii++){
      Cbintheta[ii]=(Cntheta-ii)*(90.-0.)/Cntheta;
   }
   //Set Number of Phi Bins
   Cnbinphi[0]=158;
   Cnbinphi[1]=157;
   Cnbinphi[2]=157;
   Cnbinphi[3]=156;
   Cnbinphi[4]=155;
   Cnbinphi[5]=154;
   Cnbinphi[6]=154;
   Cnbinphi[7]=151;
   Cnbinphi[8]=149;
   Cnbinphi[9]=147;
   Cnbinphi[10]=144;
   Cnbinphi[11]=142;
   Cnbinphi[12]=139;
   Cnbinphi[13]=136;
   Cnbinphi[14]=133;
   Cnbinphi[15]=129;
   Cnbinphi[16]=126;
   Cnbinphi[17]=122;
   Cnbinphi[18]=118;
   Cnbinphi[19]=114;
   Cnbinphi[20]=109;
   Cnbinphi[21]=105;
   Cnbinphi[22]=100;
   Cnbinphi[23]=95;
   Cnbinphi[24]=90;
   Cnbinphi[25]=85;
   Cnbinphi[26]=80;
   Cnbinphi[27]=75;
   Cnbinphi[28]=69;
   Cnbinphi[29]=64;
   Cnbinphi[30]=58;
   Cnbinphi[31]=52;
   Cnbinphi[32]=46;
   Cnbinphi[33]=40;
   Cnbinphi[34]=34;
   Cnbinphi[35]=28;
   Cnbinphi[36]=22;
   Cnbinphi[37]=16;
   Cnbinphi[38]=10;
   Cnbinphi[39]=4;
   //Set Phi Bins
   for(int ii=0;ii<Cntheta;ii++){
      for(int jj=0;jj<=Cnbinphi[ii];jj++){
         if(ii%2==0) Cbinphi[ii][jj]=-jj*(360.-0.)/Cnbinphi[ii]+90;
         else        Cbinphi[ii][jj]=-(Cnbinphi[ii]-jj)*(360.-0.)/Cnbinphi[ii]+90;
      }
   }
}
void Cloud::Convert(int index,double &xx,double &yy){
   for(int ii=0;ii<Cntheta;ii++){
      if(index>0&&index<=Cnbinphi[ii]){
         xx=(Cbintheta[ii]+Cbintheta[ii+1])/2.;
         yy=(Cbinphi[ii][index-1]+Cbinphi[ii][index])/2.;
         break;
      }
      else index-=Cnbinphi[ii];
   }
}
int Cloud::FindBinIndex(double xx,double yy){
   if(xx>=Cbintheta[0]) return 0;
   else if(xx<Cbintheta[Cntheta]){
      int nbin=0;
      for(int ii=0;ii<Cntheta;ii++){
         nbin+=Cnbinphi[ii];
      }
      return nbin+1;
   }
   else{
      while(yy<0){
         yy+=360;
         if(yy>=0) break;
      }
      while(yy>=360){
         yy-=360;
         if(yy<360) break;
      }

      int nbin=0;
      bool findit=false;
      for(int ii=0;ii<Cntheta;ii++){
         double theta1=Cbintheta[ii+1];
         double theta2=Cbintheta[ii];
         if(xx>=theta1&&xx<theta2){
            for(int jj=0;jj<Cnbinphi[ii];jj++){
               double phi1=(ii%2==0)?Cbinphi[ii][jj+1]:Cbinphi[ii][jj];
               double phi2=(ii%2==0)?Cbinphi[ii][jj]:Cbinphi[ii][jj+1];
               if(yy>=phi1&&yy<phi2){ //break
                  findit=true;
               }
               else nbin++;
               if(findit) break;
            }
         }
         else nbin+=Cnbinphi[ii];
         if(findit) break;
      }
      if(findit) return nbin+1;
      else return -1;
   }
}
void Cloud::Init(){
   cloudmap=0;
   time=0;
   graphlist.clear();
}
void Cloud::Reset(){
   time=0;
   if(cloudmap){
      delete cloudmap;
   }
   cloudmap=new TH2Poly();
   int npt=10;
   //double PI=3.1415926;

   int nbin=0;
   for(int ii=0;ii<Cntheta;ii++){
      double theta1=Cbintheta[ii];
      double theta2=Cbintheta[ii+1];
      for(int jj=0;jj<Cnbinphi[ii];jj++){
         double phi1=Cbinphi[ii][jj];
         double phi2=Cbinphi[ii][jj+1];
         TGraph* gr=new TGraph();
         for(int kk=0;kk<npt;kk++){
            double rr=theta1+(theta2-theta1)/npt*kk;
            gr->SetPoint(gr->GetN(),rr*cos(phi1/180*PI),rr*sin(phi1/180*PI));
         }
         for(int kk=0;kk<npt;kk++){
            double angle=phi1+(phi2-phi1)/npt*kk;
            gr->SetPoint(gr->GetN(),theta2*cos(angle/180*PI),theta2*sin(angle/180*PI));
         }
         for(int kk=0;kk<npt;kk++){
            double rr=theta2+(theta1-theta2)/npt*kk;
            gr->SetPoint(gr->GetN(),rr*cos(phi2/180*PI),rr*sin(phi2/180*PI));
         }
         for(int kk=0;kk<npt&&theta1>0;kk++){
            double angle=phi2+(phi1-phi2)/npt*kk;
            gr->SetPoint(gr->GetN(),theta1*cos(angle/180*PI),theta1*sin(angle/180*PI));
         }
         //gr->SetPoint(gr->GetN(),theta1*cos(phi1/180*PI),theta1*sin(phi1/180*PI));
         cloudmap->AddBin((TObject*)gr);
         //graphlist.push_back(gr);//delete gr;
         nbin++;
         double theta0=(theta1+theta2)/2.;
         double phi0=(phi1+phi2)/2.;
      }
   }
   cloudmap->SetTitle(Form(";West<-->East;South<-->North"));
   cloudmap->GetXaxis()->CenterTitle(1);
   cloudmap->GetYaxis()->CenterTitle(1);
}
void Cloud::Clear(){
   if(cloudmap){
      delete cloudmap;
      cloudmap=0;
   }
   for(int ii=0;ii<graphlist.size();ii++){
      TGraph* gr=graphlist.at(ii);
      if(gr) delete gr;
   }
   graphlist.clear();
}
void Cloud::ReadCloudMap(char* filename){
   if(!cloudmap) Reset();
   else cloudmap->TH2::Reset();
   ifstream fin(filename,std::ios::in);
   double x0;
   int index;
   double ti,xi;
   fin>>x0;
   while(!fin.eof()){
      fin>>index>>ti>>xi;
      double theta,phi;
      //Convert(index,theta,phi);
      cloudmap->SetBinContent(index,xi);
      //printf("Cloud::ReadCloudMap: bin%d temp=%.2lf\n",index,cloudmap->GetBinContent(index));
   }
   time=CommonTools::Convert(ti);
   fin.close();
}
TGraph* Cloud::TelView(WFTelescopeArray* pct,int iTel){
   if(!pct) return 0;
   else{
      TGraph* gr=pct->TelView(iTel);
      if(!gr) return 0;
      TGraph* gr2=new TGraph();
      //double PI=3.1415926;
      for(int ii=0;ii<gr->GetN();ii++){
         double xx=gr->GetX()[ii];
         double yy=gr->GetY()[ii];
         double zz=sqrt(1-xx*xx-yy*yy);
         double theta=TMath::ACos(zz)/PI*180.;
	 double phi=(xx==0)?(yy>=0?PI/2.:-PI/2.):TMath::ATan(yy/xx);
         if(xx<0) phi+=PI;
         if(phi<0) phi+=2*PI;
         gr2->SetPoint(ii,theta*cos(phi),theta*sin(phi));
      }
      delete gr;
      return gr2;
   }
}
void Cloud::AveTemp(double &avetemp,double &mintemp,TGraph* gr){
   avetemp=-1;
   mintemp=-1;
   if(!cloudmap) return;
   if(!gr) return;
   int np=0;
   double tempmin=100;
   double tempave=0;

   //double PI=3.1415926;
   int nbins=0;
   for(int ii=0;ii<Cntheta;ii++){
      nbins+=Cnbinphi[ii];
   }
   for(int ibin=1;ibin<=nbins;ibin++){
      double theta,phi;
      Convert(ibin,theta,phi);
      bool inside=gr->IsInside(theta*cos(phi/180*PI),theta*sin(phi/180*PI));
      if(inside){
         double xi=cloudmap->GetBinContent(ibin);
         if(xi<tempmin) tempmin=xi;
         tempave+=xi;
         np++;
      }
   }
   if(np>0){
      avetemp=tempave/np;
      mintemp=tempmin;
   }
}
void Cloud::Draw(WFTelescopeArray* pct,char* opt){
   if(!cloudmap) return;
   cloudmap->Draw(opt);
   printf("Cloud::Draw: Draw CloudMap\n");
   double tempave=0,tempmin=0;
   for(int itel=0;itel<WFTelescopeArray::CTNumber&&pct;itel++){
      TGraph* gr=TelView(pct,itel);
      if(!gr) continue;
      gr->Draw("l");
      printf("Cloud::Draw: Draw Telescope %d\n",itel);
      gr->SetLineColor(1);
      gr->SetLineWidth(2);
      if(itel==0) AveTemp(tempave,tempmin,gr);
      graphlist.push_back(gr);
   }
   if(drawcircle>0){
      for(int ii=0;ii<=drawcircle;ii++){
         double theta=0+(90.-0.)/drawcircle*ii;
         if(theta<=0) continue;
         TGraph* gc=new TGraph();
         for(int jj=0;jj<=100;jj++){
            double phi=0+(2*3.1415926-0)/100*jj;
            gc->SetPoint(jj,theta*cos(phi),theta*sin(phi));
         }
         gc->SetLineColor(12);
         gc->SetLineStyle(2);
         gc->SetLineWidth(2);
         gc->Draw("l");
         printf("Cloud::Draw: Draw Circle %d(%.2lf)\n",ii,theta);
         graphlist.push_back(gc);
      }
   }
   int year,month,day,hour,min;
   year=CommonTools::TimeFlag(time,1);
   month=CommonTools::TimeFlag(time,2);
   day=CommonTools::TimeFlag(time,3);
   hour=CommonTools::TimeFlag(time,4);
   min=CommonTools::TimeFlag(time,5);
   cloudmap->SetTitle(Form("20%02d-%02d-%02d %02d:%02d Sky IR Image (Average Temp=%.2f,Minimum Temp=%.2f,in the field of view)",year,month,day,hour,min,tempave,tempmin));
}
