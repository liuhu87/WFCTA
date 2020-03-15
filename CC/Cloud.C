#include "Cloud.h"
#include "common.h"
#include "math.h"
#include "TMath.h"
#include <fstream>
#include "WFTelescope.h"
#include "TText.h"
#include "TFile.h"
#include "TTree.h"
#include "slalib.h"
#include "astro.h"
bool Cloud::DoTempCorr=false;
bool Cloud::drawmoon=true;
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
bool Cloud::Convert(int index,double &xx,double &yy,double xyboun[][2]){
   bool findit=false;
   for(int ii=0;ii<Cntheta;ii++){
      if(index>0&&index<=Cnbinphi[ii]){
         xx=(Cbintheta[ii]+Cbintheta[ii+1])/2.;
         yy=(Cbinphi[ii][index-1]+Cbinphi[ii][index])/2.;
         findit=true;
         xyboun[0][0]=Cbintheta[ii];
         xyboun[0][1]=Cbintheta[ii+1];
         xyboun[1][0]=Cbinphi[ii][index-1];
         xyboun[1][1]=Cbinphi[ii][index];
         break;
      }
      else index-=Cnbinphi[ii];
   }
   return findit;
}
int Cloud::FindBinIndex(double xx,double yy){
   if(xx>Cbintheta[0]) return 0;
   else if(xx<Cbintheta[Cntheta]){
      int nbin=0;
      for(int ii=0;ii<Cntheta;ii++){
         nbin+=Cnbinphi[ii];
      }
      return nbin+1;
   }
   else{
      yy=CommonTools::ProcessAngle(yy,true);

      int nbin=0;
      bool findit=false;
      for(int ii=0;ii<Cntheta;ii++){
         double theta1=Cbintheta[ii+1];
         double theta2=Cbintheta[ii];
         if((xx>theta1&&xx<=theta2)||(xx==0&&(xx>=theta1&&xx<=theta2))){
            for(int jj=0;jj<Cnbinphi[ii];jj++){
               double phi1=CommonTools::ProcessAngle((ii%2==0)?Cbinphi[ii][jj+1]:Cbinphi[ii][jj],true);
               double phi2=CommonTools::ProcessAngle((ii%2==0)?Cbinphi[ii][jj]:Cbinphi[ii][jj+1],true);
               double phimin=TMath::Min(phi1,phi2);
               double phimax=TMath::Max(phi1,phi2);
               bool inside;
               if(fabs(phimax-phimin)>=180.){
                  inside=(yy>=0&&yy<=phimin)||(yy>phimax&&yy<=360);
               }
               else{
                  inside=yy>phimin&&yy<=phimax;
               }
               if(inside){ //break
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
   mapnbins=0;
   time=0;
   temp0=1000;
   temp=0;
   humi=-1;
   graphlist.clear();
}
void Cloud::Reset(){
   time=0;
   temp0=1000;
   temp=0;
   humi=-1;
   if(cloudmap){
      delete cloudmap;
   }
   cloudmap=new TH2Poly();
   int npt=10;
   //double PI=3.1415926;

   mapnbins=0;
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
         mapnbins++;
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
   if(!fin.is_open()){
      printf("Cloud::ReadCloudMap: open file %s failed\n",filename);
      return;
   }
   double x0;
   int index;
   double ti,xi;
   fin>>x0;
   temp0=x0;
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
bool Cloud::ReadTemp(char* filename){
   bool res=false;
   if(!filename){
      if(time<100) return false;
      int time1;
      int hour=CommonTools::TimeFlag(time,4);
      int min=CommonTools::TimeFlag(time,5);
      bool before2005=(hour*100+min)<2005;
      if(before2005) time1=time+(24*3600);
      else time1=time+(24*3600*2);
      int year=(CommonTools::TimeFlag(time1,1)%100)+2000;
      int month=CommonTools::TimeFlag(time1,2);
      int day=CommonTools::TimeFlag(time1,3);
      //filename=Form("/scratchfs/ybj/lix/laser-dat/Temp-humi/%02d/temp/12345_cloud_temp_%04d%02d%02d.txt",month,year,month,day);
      //filename=Form("/eos/lhaaso/raw/wfctalaser/%04d/%02d%02d/12345_cloud_temp_%04d%02d%02d.txt",year,month,day,year,month,day);
      filename=Form("/eos/lhaaso/raw/wfctalaser/IBTemp/TempHumdata/%04d/%02d/temp/12345_cloud_temp_%04d%02d%02d.txt",year,month,year,month,day);
      printf("ReadTemp: time=%d filename=%s\n",time,filename);
   }
   ifstream fin(filename,std::ios::in);
   if(!fin.is_open()) {temp=0; humi=-1; return res;}
   double ti,tempi,hi;
   int index;
   while(!fin.eof()){
      fin>>ti>>index>>tempi>>index>>hi>>index;
      int time1=CommonTools::Convert(ti*100);
      //printf("time=%d time1=%d temp=%lf hi=%lf\n",time,time1,tempi,hi);
      if(abs(time-time1)<100){
         temp=tempi;
         humi=hi;
         res=true;
         break;
      }
   }
   if(!res){
      temp=0;
      humi=-1;
   }
   fin.close();
   return res;
}
void Cloud::LoadTelSetting(char* filename){
   WFTelescopeArray::GetHead(filename);
}
TGraph* Cloud::TelView(WFTelescopeArray* pct,int iTel){
   if(!pct) return 0;
   else{
      TGraph* gr=pct->TelView(iTel,false);
      if(!gr) return 0;
      return gr;
      //TGraph* gr2=new TGraph();
      ////double PI=3.1415926;
      //for(int ii=0;ii<gr->GetN();ii++){
      //   double xx=gr->GetX()[ii];
      //   double yy=gr->GetY()[ii];
      //   double zz=sqrt(1-xx*xx-yy*yy);
      //   double theta=TMath::ACos(zz)/PI*180.;
      //   double rr=sqrt(xx*xx+yy*yy);
      //   double phi=acos(-yy/rr);
      //   if(xx<0) phi=2*PI-phi;
      //   gr2->SetPoint(ii,theta*cos(phi),theta*sin(phi));
      //}
      //delete gr;
      //return gr2;
   }
}
void Cloud::AveTemp(double &avetemp,double &mintemp,double &rmstemp,TGraph* gr){
   avetemp=1000;
   mintemp=1000;
   rmstemp=0;
   if(!cloudmap) return;
   //if(!gr) return;
   int np=0;
   double tempmin=1000;
   double tempave=0;
   double temp2=0;

   //double PI=3.1415926;
   int nbins=GetNbins();
   for(int ibin=1;ibin<=nbins;ibin++){
      double xyboun[2][2];
      double theta,phi;
      Convert(ibin,theta,phi,xyboun);
      bool inside=(!gr)?true:gr->IsInside(theta*cos(phi/180*PI),theta*sin(phi/180*PI));
      if(inside){
         double xi=cloudmap->GetBinContent(ibin);
         xi=GetCorrected(xi);
         if(xi<tempmin) tempmin=xi;
         tempave+=xi;
         temp2+=xi*xi;
         np++;
      }
   }
   if(np>0){
      tempave/=np;
      temp2=sqrt(temp2/np-tempave*tempave);
   }
   else{
      tempave=1000;
      temp2=-1;
   }
   avetemp=tempave;
   mintemp=tempmin;
   rmstemp=temp2;
}
void Cloud::Draw(WFTelescopeArray* pct,char* opt){
   if(!cloudmap) return;
   cloudmap->Draw(opt);
   printf("Cloud::Draw: Draw CloudMap\n");
   double tempave[NCTMax],tempmin[NCTMax];
   for(int ii=0;ii<NCTMax;ii++){
      tempave[ii]=0;
      tempmin[ii]=0;
   }
   for(int itel=0;itel<WFTelescopeArray::CTNumber&&pct;itel++){
      WFTelescope* pt=pct->pct[itel];
      if(!pt) continue;
      int iTel=pt->TelIndex_;
      TGraph* gr=TelView(pct,iTel);
      if(!gr) continue;
      gr->Draw("l");
      printf("Cloud::Draw: Draw Telescope %d\n",iTel);
      gr->SetLineColor(1);
      gr->SetLineWidth(2);
      double rms;
      AveTemp(tempave[itel],tempmin[itel],rms,gr);
      graphlist.push_back(gr);
      //printf("Cloud::Draw: iTel=%d(Tot=%d) avetemp=%.2lf tempmin=%.2lf\n",iTel,WFTelescopeArray::CTNumber,tempave[itel],tempmin[itel]);

      int telindex=pt->TelIndex_;
      double telzen=pt->TelZ_;
      double telazi=pt->TelA_;
      double xx=telzen/PI*180.*cos(telazi);
      double yy=telzen/PI*180.*sin(telazi);
      TText* text=new TText(-yy,xx,Form("%d",telindex));
      text->SetTextFont(12);
      text->SetTextSize(0.05);
      text->Draw("same");
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
   int year,month,day,hour,min,sec;
   year=CommonTools::TimeFlag(time,1);
   month=CommonTools::TimeFlag(time,2);
   day=CommonTools::TimeFlag(time,3);
   hour=CommonTools::TimeFlag(time,4);
   min=CommonTools::TimeFlag(time,5);
   sec=CommonTools::TimeFlag(time,6);
   if(drawmoon){
      TGraph* gm=new TGraph();

      double MJD19700101=40587;
      for(int ii=0;ii<900;ii++){
      double time0=time-ii*60;
      double mjd=MJD19700101+(time0-37)/86400.;
      double jd=mjd+2400000.5;
      //double PV[6];
      //slaDmoon(mjd,PV);
      //double unit=180/PI;
      //printf("time=%d(%d-%d-%d %d:%d:%d) PV={%le,%le,%le,%le,%le,%le}\n",time0,year,month,day,hour,min,sec,PV[0]*unit,PV[1]*unit,PV[2]*unit,PV[3]*unit,PV[4]*unit,PV[5]*unit);

      float lam,bet,lamdot,betdot,rsun;
      astroSunmon(jd,1,&lam,&bet,&lamdot,&betdot,&rsun);
      float ra,dec;
      astroEctoeq(lam,bet,&ra,&dec);

      double ha = astroLST(jd) - ra;
      double el,az;
      slaDe2h(ha, dec, DUGLAT, &az, &el);
      double xx=(PI/2-el)*TMath::RadToDeg()*cos(PI/2-az);
      double yy=(PI/2-el)*TMath::RadToDeg()*sin(PI/2-az);
      double unit=180/PI;
      //printf("time=%.0lf(%d-%02d-%02d %02d:%02d:%02d) dir={%lf,%lf} {%lf,%lf}\n",time0,CommonTools::TimeFlag((int)time0,1),CommonTools::TimeFlag((int)time0,2),CommonTools::TimeFlag((int)time0,3),CommonTools::TimeFlag((int)time0,4),CommonTools::TimeFlag((int)time0,5),CommonTools::TimeFlag((int)time0,6),el*unit,az*unit,xx,yy);
      if(fabs(xx)>90||fabs(yy)>90) break;
      gm->SetPoint(gm->GetN(),xx,yy);
      }

      /*double time_target=hour*3600+min*60+sec;
      TFile *file = TFile::Open("/afs/ihep.ac.cn/users/y/youzhiyong/moon-orbit/Moon_orbit_night1_2019.root","read");
      TTree *tree = (TTree*)file->Get("tree");
      double el,az;
      int y, mon, d, h, min2,sec2;
      tree->SetBranchAddress("az",&az);
      tree->SetBranchAddress("el",&el);
      tree->SetBranchAddress("y",&y);
      tree->SetBranchAddress("mon",&mon);
      tree->SetBranchAddress("d",&d);
      tree->SetBranchAddress("h",&h);
      tree->SetBranchAddress("min",&min2);
      sec2=0;
      int npots=tree->GetEntries();
      for(int ii =0;ii<npots;ii++){
         tree->GetEntry(ii);
         double time_entry=h*3600+min2*60+sec2;
         //if((y==(year+2000)&mon==month&d==day)&&(time_entry<time_target&&el<PI/2)) 
         if((y==(year+2000)&mon==month&d==day)&&(time_entry*0<time_target)) {
            gm->SetPoint(gm->GetN(),(PI/2-el)*TMath::RadToDeg()*cos(PI/2-az),(PI/2-el)*TMath::RadToDeg()*sin(PI/2-az));
            printf("ii=%d time={%lf(%d-%d-%d),%lf(%d-%d-%d)} el=%lf az=%lf\n",ii,time_target,hour,min,sec,time_entry,h,min2,sec2,90-el/PI*180,90-az/PI*180);
         }
         //printf("Drawmoon: ip=%d year={%d,%d} month={%d,%d} day={%d,%d}\n",ii,y,year,mon,month,d,day);
      }*/
      gm->SetLineColor(1);
      gm->SetLineWidth(3);
      gm->SetMarkerStyle(20);
      gm->SetMarkerSize(0.3);
      gm->Draw("p");
      printf("Cloud::Draw: Draw Moon Trajectory\n");
      graphlist.push_back(gm);
   }
   if(WFTelescopeArray::CTNumber==0) cloudmap->SetTitle(Form("20%02d-%02d-%02d %02d:%02d Sky IR Image",year,month,day,hour,min));
   else if(WFTelescopeArray::CTNumber==1){
      cloudmap->SetTitle(Form("20%02d-%02d-%02d %02d:%02d Sky IR Image (Tel No={%d,%d},Ave Temp={%.2f},Min Temp={%.2f})",year,month,day,hour,min,pct->pct[0]->TelIndex_,tempave[0],tempmin[0]));
   }
   else if(WFTelescopeArray::CTNumber==2){
      cloudmap->SetTitle(Form("20%02d-%02d-%02d %02d:%02d Sky IR Image (Tel No={%d,%d},Ave Temp={%.2f,%.2f},Min Temp={%.2f,%.2f})",year,month,day,hour,min,pct->pct[0]->TelIndex_,pct->pct[1]->TelIndex_,tempave[0],tempave[1],tempmin[0],tempmin[1]));
   }
   else if(WFTelescopeArray::CTNumber==3){
      cloudmap->SetTitle(Form("20%02d-%02d-%02d %02d:%02d Sky IR Image (Tel No={%d,%d,%d},Ave Temp={%.2f,%.2f,%.2f},Min Temp={%.2f,%.2f,%.2f})",year,month,day,hour,min,pct->pct[0]->TelIndex_,pct->pct[1]->TelIndex_,pct->pct[2]->TelIndex_,tempave[0],tempave[1],tempave[2],tempmin[0],tempmin[1],tempmin[2]));
   }
   else if(WFTelescopeArray::CTNumber==4){
      cloudmap->SetTitle(Form("20%02d-%02d-%02d %02d:%02d Sky IR Image (Tel No={%d,%d,%d,%d},Ave Temp={%.2f,%.2f,%.2f,%.2f},Min Temp={%.2f,%.2f,%.2f,%.2f})",year,month,day,hour,min,pct->pct[0]->TelIndex_,pct->pct[1]->TelIndex_,pct->pct[2]->TelIndex_,pct->pct[3]->TelIndex_,tempave[0],tempave[1],tempave[2],tempave[3],tempmin[0],tempmin[1],tempmin[2],tempmin[3]));
   }
   else if(WFTelescopeArray::CTNumber==5){
      cloudmap->SetTitle(Form("20%02d-%02d-%02d %02d:%02d Sky IR Image (Tel No={%d,%d,%d,%d,%d},Ave Temp={%.2f,%.2f,%.2f,%.2f,%.2lf},Min Temp={%.2f,%.2f,%.2f,%.2f,%.2lf})",year,month,day,hour,min,pct->pct[0]->TelIndex_,pct->pct[1]->TelIndex_,pct->pct[2]->TelIndex_,pct->pct[3]->TelIndex_,pct->pct[4]->TelIndex_,tempave[0],tempave[1],tempave[2],tempave[3],tempave[4],tempmin[0],tempmin[1],tempmin[2],tempmin[3],tempmin[4]));
   }
   else if(WFTelescopeArray::CTNumber==6){
      cloudmap->SetTitle(Form("20%02d-%02d-%02d %02d:%02d Sky IR Image (Tel No={%d,%d,%d,%d,%d,%d},Ave Temp={%.2f,%.2f,%.2f,%.2f,%.2lf,%.2lf},Min Temp={%.2f,%.2f,%.2f,%.2f,%.2lf,%.2lf})",year,month,day,hour,min,pct->pct[0]->TelIndex_,pct->pct[1]->TelIndex_,pct->pct[2]->TelIndex_,pct->pct[3]->TelIndex_,pct->pct[4]->TelIndex_,pct->pct[5]->TelIndex_,tempave[0],tempave[1],tempave[2],tempave[3],tempave[4],tempave[5],tempmin[0],tempmin[1],tempmin[2],tempmin[3],tempmin[4],tempmin[5]));
   }
   else if(WFTelescopeArray::CTNumber==7){
      cloudmap->SetTitle(Form("20%02d-%02d-%02d %02d:%02d Sky IR Image (Tel No={%d,%d,%d,%d,%d,%d,%d},Ave Temp={%.2f,...},Min Temp={%.2f,...})",year,month,day,hour,min,pct->pct[0]->TelIndex_,pct->pct[1]->TelIndex_,pct->pct[2]->TelIndex_,pct->pct[3]->TelIndex_,pct->pct[4]->TelIndex_,pct->pct[5]->TelIndex_,pct->pct[6]->TelIndex_,tempave[0],tempmin[0]));
   }
   else if(WFTelescopeArray::CTNumber>=8){
      cloudmap->SetTitle(Form("20%02d-%02d-%02d %02d:%02d Sky IR Image (Tel No={%d,%d,%d,%d,%d,%d,%d,%d},Ave Temp={%.2f,...},Min Temp={%.2f,...})",year,month,day,hour,min,pct->pct[0]->TelIndex_,pct->pct[1]->TelIndex_,pct->pct[2]->TelIndex_,pct->pct[3]->TelIndex_,pct->pct[4]->TelIndex_,pct->pct[5]->TelIndex_,pct->pct[6]->TelIndex_,pct->pct[7]->TelIndex_,tempave[0],tempmin[0]));
   }

}

double Cloud::GetCorrected(double input){
   return (humi>=0)?(input-temp):input;
}
double Cloud::GetTemperature(int itemp){
   return itemp==0?temp0:(itemp==1?temp:1000);
}
double Cloud::GetHumidity(){
   return humi;
}
double Cloud::GetIBTemp(int ibin){
   if(!cloudmap) return 1000;
   if(ibin<1||ibin>GetNbins()) return 1000;
   double res=cloudmap->GetBinContent(ibin);
   return GetCorrected(res);
}
double Cloud::GetIBTemp(double xx,double yy){
   if(!cloudmap) return 1000;
   int binindex=FindBinIndex(xx,yy);
   if(binindex<1||binindex>GetNbins()) return 1000;
   double res=cloudmap->GetBinContent(binindex);
   return GetCorrected(res);
}
double Cloud::GetTelAveIBTemp(int iTel){
   TGraph* gr=TelView(WFTelescopeArray::GetHead(),iTel);
   double min=1000,ave=1000,rms=-1;
   if(!gr) return ave;
   AveTemp(ave,min,rms,gr);
   return ave;
}
double Cloud::GetTelMinIBTemp(int iTel){
   TGraph* gr=TelView(WFTelescopeArray::GetHead(),iTel);
   double min=1000,ave=1000,rms=-1;
   if(!gr) return min;
   AveTemp(ave,min,rms,gr);
   return min;
}
double Cloud::GetTelRmsIBTemp(int iTel){
   TGraph* gr=TelView(WFTelescopeArray::GetHead(),iTel);
   double min=1000,ave=1000,rms=-1;
   if(!gr) return rms;
   AveTemp(ave,min,rms,gr);
   return rms;
}
double Cloud::GetAveIBTemp(double theta){
   double min=1000,ave=0;
   int np=0;
   int nbins=GetNbins();
   for(int ibin=1;ibin<=nbins;ibin++){
      double xyboun[2][2];
      double thetai,phii;
      Convert(ibin,thetai,phii,xyboun);
      bool inside=thetai<theta;
      if(inside){
         double xi=cloudmap->GetBinContent(ibin);
         xi=GetCorrected(xi);
         if(xi<min) min=xi;
         ave+=xi;
         np++;
      }
   }
   if(np<=0) ave=1000;
   else ave/=np;
   return ave;
}
double Cloud::GetMinIBTemp(double theta){
   double min=1000,ave=0;
   int np=0;
   int nbins=GetNbins();
   for(int ibin=1;ibin<=nbins;ibin++){
      double xyboun[2][2];
      double thetai,phii;
      Convert(ibin,thetai,phii,xyboun);
      bool inside=thetai<theta;
      if(inside){
         double xi=cloudmap->GetBinContent(ibin);
         xi=GetCorrected(xi);
         if(xi<min) min=xi;
         ave+=xi;
         np++;
      }
   }
   if(np<=0) ave=1000;
   else ave/=np;
   return min;
}
double Cloud::GetRmsIBTemp(double theta){
   double min=1000,ave=0,rms=0;
   int np=0;
   int nbins=GetNbins();
   for(int ibin=1;ibin<=nbins;ibin++){
      double xyboun[2][2];
      double thetai,phii;
      Convert(ibin,thetai,phii,xyboun);
      bool inside=thetai<theta;
      if(inside){
         double xi=cloudmap->GetBinContent(ibin);
         xi=GetCorrected(xi);
         if(xi<min) min=xi;
         ave+=xi;
         rms+=xi*xi;
         np++;
      }
   }
   if(np<=0) rms=-1;
   else{
      ave/=np;
      rms=sqrt(rms/np-ave*ave);
   }
   return rms;
}
double Cloud::GetAveIBTemp(double theta1,double theta2){
   double min=1000,ave=0;
   int np=0;
   int nbins=GetNbins();
   for(int ibin=1;ibin<=nbins;ibin++){
      double xyboun[2][2];
      double thetai,phii;
      Convert(ibin,thetai,phii,xyboun);
      double boun0=TMath::Min(xyboun[0][0],xyboun[0][1]);
      double boun1=TMath::Max(xyboun[0][0],xyboun[0][1]);
      boun0=TMath::Max(boun0,theta1);
      boun1=TMath::Min(boun1,theta2);
      bool inside=(boun0<boun1);
      if(inside){
         double xi=cloudmap->GetBinContent(ibin);
         xi=GetCorrected(xi);
         if(xi<min) min=xi;
         ave+=xi;
         np++;
      }
   }
   if(np<=0) ave=1000;
   else ave/=np;
   return ave;
}
double Cloud::GetMinIBTemp(double theta1,double theta2){
   double min=1000,ave=0;
   int np=0;
   int nbins=GetNbins();
   for(int ibin=1;ibin<=nbins;ibin++){
      double xyboun[2][2];
      double thetai,phii;
      Convert(ibin,thetai,phii,xyboun);
      double boun0=TMath::Min(xyboun[0][0],xyboun[0][1]);
      double boun1=TMath::Max(xyboun[0][0],xyboun[0][1]);
      boun0=TMath::Max(boun0,theta1);
      boun1=TMath::Min(boun1,theta2);
      bool inside=(boun0<boun1);
      if(inside){
         double xi=cloudmap->GetBinContent(ibin);
         xi=GetCorrected(xi);
         if(xi<min) min=xi;
         ave+=xi;
         np++;
      }
   }
   if(np<=0) ave=1000;
   else ave/=np;
   return min;
}
double Cloud::GetRmsIBTemp(double theta1,double theta2){
   double min=1000,ave=0,rms=0;
   int np=0;
   int nbins=GetNbins();
   for(int ibin=1;ibin<=nbins;ibin++){
      double xyboun[2][2];
      double thetai,phii;
      Convert(ibin,thetai,phii,xyboun);
      double boun0=TMath::Min(xyboun[0][0],xyboun[0][1]);
      double boun1=TMath::Max(xyboun[0][0],xyboun[0][1]);
      boun0=TMath::Max(boun0,theta1);
      boun1=TMath::Min(boun1,theta2);
      bool inside=(boun0<boun1);
      if(inside){
         double xi=cloudmap->GetBinContent(ibin);
         xi=GetCorrected(xi);
         if(xi<min) min=xi;
         ave+=xi;
         rms+=xi*xi;
         np++;
      }
   }
   if(np<=0) rms=-1;
   else{
      ave/=np;
      rms=sqrt(rms/np-ave*ave);
   }
   return rms;
}
double Cloud::GetAveIBTemp(TGraph* gr){
   double min=1000,ave=1000,rms=-1;
   if(!gr) return ave;
   AveTemp(ave,min,rms,gr);
   return ave;
}
double Cloud::GetMinIBTemp(TGraph* gr){
   double min=1000,ave=1000,rms=-1;
   if(!gr) return min;
   AveTemp(ave,min,rms,gr);
   return min;
}
double Cloud::GetRmsIBTemp(TGraph* gr){
   double min=1000,ave=1000,rms=-1;
   if(!gr) return rms;
   AveTemp(ave,min,rms,gr);
   return rms;
}

