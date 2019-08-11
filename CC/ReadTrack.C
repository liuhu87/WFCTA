#include "ReadTrack.h"
#include "CorsikaEvent.h"
ReadTrack* ReadTrack::_Head=0;
bool ReadTrack::DoPlot=false;
int ReadTrack::headbyte=16;
int ReadTrack::jdebug=0;
long int ReadTrack::particle=-1;
float ReadTrack::elimit[2]={0,0};
float ReadTrack::climit[3][2]={{0,0},{0,0},{0,0}};
float ReadTrack::tlimit[2]={0,0};
float ReadTrack::IniRange[4][2]={{1.0e9,-1.0e9},{1.0e9,-1.0e9},{1.0e9,-1.0e9},{1.0e9,-1.0e9}};
void ReadTrack::SetHead(ReadTrack* head){
   _Head=head;
}
ReadTrack* ReadTrack::GetHead(){
   return _Head;
}
void ReadTrack::Init(){
   fin=0;
   type=-1;
   plot=0;
   nrec=0;
   int nsize=1000000;
   arrid.resize(nsize);
   arrid.clear();
   arren.resize(nsize);
   arren.clear();
   arrt1.resize(nsize);
   arrt1.clear();
   arrt2.resize(nsize);
   arrt2.clear();
   for(int ii=0;ii<3;ii++){
      arrc1[ii].resize(nsize);
      arrc1[ii].clear();
      arrc2[ii].resize(nsize);
      arrc2[ii].clear();
      plotrange[ii][0]=IniRange[ii][0];
      plotrange[ii][1]=IniRange[ii][1];
   }
   plotrange[3][0]=IniRange[3][0];
   plotrange[3][1]=IniRange[3][1];
   for(int ii=0;ii<10;ii++) recbuff.f[ii]=0;
}
void ReadTrack::Reset(){
   if(fin){
      fin->clear();
      fin->seekg(0,ios::beg);
   }
   if(plot){
      //plot->Delete();
      delete plot;
      plot=new TObjArray();
   }
   nrec=0;
   arrid.clear();
   arren.clear();
   arrt1.clear();
   arrt2.clear();
   for(int ii=0;ii<3;ii++){
      arrc1[ii].clear();
      arrc2[ii].clear();
      plotrange[ii][0]=IniRange[ii][0];
      plotrange[ii][1]=IniRange[ii][1];
   }
   plotrange[3][0]=IniRange[3][0];
   plotrange[3][1]=IniRange[3][1];
}
void ReadTrack::Release(){
   if(fin){
      fin->close();
      delete fin;
   }
   if(plot){
      plot->Delete();
      delete plot;
      plot=0;
   }
   arrid.clear();
   arren.clear();
   arrt1.clear();
   arrt2.clear();
   for(int ii=0;ii<3;ii++){
      arrc1[ii].clear();
      arrc2[ii].clear();
   }
}
ReadTrack::ReadTrack(const char* inputfile){
   Init();
   if(!inputfile) return;
   fin=new std::ifstream(inputfile);
   if(!fin) return;
   if(!fin->is_open()){
      printf("ReadTrack::ReadTrack:Error in opening file %s\n",inputfile);
      delete fin; fin=0;
      return;
   }
   else{
      fin->seekg(0,ios::end);
      size_t filesize = fin->tellg();
      if (jdebug>0) printf("The file size of %s is %d bytes\n",inputfile,filesize);
      if (!filesize) {
        printf("Warning: The file is empty!\n");
        fin->close();
        delete fin; fin=0;
        return;
      }
      fin->seekg(0,ios::beg);
      if(strstr(inputfile,"em")) type=0;
      else if(strstr(inputfile,"mu")) type=1;
      else if(strstr(inputfile,"hd")) type=2;
      else if(strstr(inputfile,"CER")) type=3;
      else type=4;
      if(jdebug>0) printf("ReadTrack::ReadTrack filename=%s type=%d\n",inputfile,type);
      return;
   }
}
bool ReadTrack::Exist(){
   return (bool)fin;
}
long int ReadTrack::GetParticleID(int partid){
   if(partid<=15){
      return (1<<partid);
   }
   else if(partid<200){
      return (1<<16);
   }
   else{
      int AA=partid/100;
      int ZZ=partid%100;
      return (1<<(17+ZZ));
   }
}
int ReadTrack::ReadRec(){
   if(!fin) return -3;
   union{
      int i[4];
      char c[16];
   } padding;
   for(int ii=0;ii<4;ii++) padding.i[ii]=0;
   if(nrec==0) fin->read(padding.c,headbyte/2);
   else fin->read(padding.c,headbyte);
   if (!fin->good()){
      if(nrec>0){
         printf("ReadTrack::ReadRec: reading track file finished.\n");
         return nrec;
      }
      else{
         printf("ReadTrack::ReadRec: Error in reading track file. The file is empty? nrec=%d pos=%ld\n",nrec,fin->tellg());
         return -2;
      }
   }
   if(jdebug>1) printf("ReadTrack::ReadRec: size of this record(%d) %d(%d %s)\n",nrec+1,padding.i[0],padding.i[1],padding.c);
   //int readsize=(nrec==0)?padding.i[0]:padding.i[1];
   int readsize=padding.i[0];
   if(readsize!=40){
      printf("ReadTrack::ReadRec: Error of the record size(%d,%d)\n",padding.i[0],padding.i[1]);
      return -1;
   }
   else fin->read(recbuff.c,readsize);
   nrec++;

   int partid;
   float energy,x1,y1,z1,t1,x2,y2,z2,t2;
   partid=recbuff.f[0];
   energy=recbuff.f[1];
   x1=recbuff.f[2];
   y1=recbuff.f[3];
   z1=recbuff.f[4];
   t1=recbuff.f[5];
   x2=recbuff.f[6];
   y2=recbuff.f[7];
   z2=recbuff.f[8];
   t2=recbuff.f[9];
   if(jdebug>3){
      printf("ReadTrack::ReadRec: Content partid=%d energy=%f start={%f,%f,%f,%f} end={%f,%f,%f,%f}\n",partid,energy,x1,y1,z1,t1,x2,y2,z2,t2);
   }
   bool inside=true;
   if(particle>=0) inside=inside&&(particle&GetParticleID(partid));
   if(elimit[1]>elimit[0]) inside=inside&&(energy>=elimit[0]&&energy<=elimit[1]);
   if(tlimit[1]>tlimit[0]){
      inside=inside&&(t1>=tlimit[0]&&t1<=tlimit[1]);
      inside=inside&&(t2>=tlimit[0]&&t2<=tlimit[1]);
   }
   float coo1[3]={x1,y1,z1};
   float coo2[3]={x2,y2,z2};
   for(int ii=0;ii<3;ii++){
      if(climit[ii][1]>climit[ii][0]){
         inside=inside&&(coo1[ii]>=climit[ii][0]&&coo1[ii]<=climit[ii][1]);
         inside=inside&&(coo2[ii]>=climit[ii][0]&&coo2[ii]<=climit[ii][1]);
      }
   }
   if(inside){
      arrid.push_back(partid);
      arren.push_back(energy);
      arrt1.push_back(t1);
      arrt2.push_back(t2);
      if(t1<plotrange[3][0]) plotrange[3][0]=t1;
      if(t2<plotrange[3][0]) plotrange[3][0]=t2;
      if(t1>plotrange[3][1]) plotrange[3][1]=t1;
      if(t2>plotrange[3][1]) plotrange[3][1]=t2;
      for(int ii=0;ii<3;ii++){
         arrc1[ii].push_back(coo1[ii]);
         arrc2[ii].push_back(coo2[ii]);
         if(IniRange[ii][0]>=IniRange[ii][1]){
         if(coo1[ii]<plotrange[ii][0]) plotrange[ii][0]=coo1[ii];
         if(coo2[ii]<plotrange[ii][0]) plotrange[ii][0]=coo2[ii];
         if(coo1[ii]>plotrange[ii][1]) plotrange[ii][1]=coo1[ii];
         if(coo2[ii]>plotrange[ii][1]) plotrange[ii][1]=coo2[ii];
         }
      }
      if(jdebug>1) printf("ReadTrack::ReadRec: fill the track nrec=%d partid=%d\n",nrec,partid);
   }
   return nrec;
}
int ReadTrack::ReadAll(int beg,int end){
   if(!fin) return -1;
   int nrec0=nrec;
   while(fin->good()){
      if(jdebug>4) printf("ReadTrack::ReadAll: Begin to Read No. %d Record\n",nrec+1);
      if(nrec+1>end&&(end>=beg&&end>0)) break;
      int irec=ReadRec();
      if(irec<=0){
         printf("ReadTrack::ReadAll: Error in Reading No. %d Record.\n",nrec>=0?(nrec+1):1);
         break;
      }
   }
   return nrec-nrec0;
}
void ReadTrack::Copy(CorsikaEvent* pevt){
   if(!pevt) return;
   if(plot) {
      //plot->Delete();
      delete plot;
      plot=new TObjArray();
   }
   plot=new TObjArray();
   int size=pevt->cx.size();
   double x0=(pevt->oheight-pevt->stheight)/(-cos(pevt->thetap))*sin(pevt->thetap)*cos(pevt->phip);
   double y0=(pevt->oheight-pevt->stheight)/(-cos(pevt->thetap))*sin(pevt->thetap)*sin(pevt->phip);
   for(int ic=0;ic<size;ic++){
      int partid=0;
      double energy=(pevt->wavelength.at(ic)>0)?(hplank_gev*vlight/(1.0e-7*pevt->wavelength.at(ic))):1.0e8;
      double coo1[3];
      double coo2[3];
      double t1,t2;
      double dxdr=pevt->cu.at(ic);
      double dydr=pevt->cv.at(ic);
      double dzdr=(1-dxdr*dxdr-dydr*dydr>=0)?-sqrt(1-dxdr*dxdr-dydr*dydr):-1;
      int whichtel=-1;
      double mindist=1.0e10;
      for(int ii=0;ii<Nuse;ii++){
         double dist=fabs(pow(pevt->corex[ii]-pevt->cx.at(ic),2)+pow(pevt->corey[ii]-pevt->cy.at(ic),2));
         if(dist<mindist){
            whichtel=ii;
            mindist=dist;
         }
      }

      //coo1[0]=pevt->cx.at(ic);
      //coo1[1]=pevt->cy.at(ic);
      //coo1[2]=pevt->height.at(ic);
      //coo2[2]=pevt->oheight;
      //coo2[0]=(coo2[2]-coo1[2])/dzdr*dxdr+coo1[0];
      //coo2[1]=(coo2[2]-coo1[2])/dzdr*dydr+coo1[1];
      //t1=pevt->ct.at(ic)*1.0e-9;
      //t2=t1+sqrt(pow(coo1[0]-coo2[0],2)+pow(coo1[1]-coo2[1],2)+pow(coo1[2]-coo2[2],2))/vlight;
      ////t2=pevt->ct.at(ic)*1.0e-9;
      ////t1=t2-sqrt(pow(coo1[0]-coo2[0],2)+pow(coo1[1]-coo2[1],2)+pow(coo1[2]-coo2[2],2))/vlight;

      coo2[0]=pevt->cx.at(ic)+x0;
      coo2[1]=pevt->cy.at(ic)+y0;
      coo2[2]=pevt->oheight;
      coo1[2]=pevt->height.at(ic);
      coo1[0]=(coo1[2]-coo2[2])/dzdr*dxdr+coo2[0];
      coo1[1]=(coo1[2]-coo2[2])/dzdr*dydr+coo2[1];
      t2=pevt->ct.at(ic)*1.0e-9;
      t1=t2-sqrt(pow(coo1[0]-coo2[0],2)+pow(coo1[1]-coo2[1],2)+pow(coo1[2]-coo2[2],2))/vlight;

      if(jdebug>9) printf("light=%d x=%f(corx=%f) y=%f(corey=%f) z=%f oz=%f\n",ic,pevt->cx.at(ic),pevt->corex[whichtel],pevt->cy.at(ic),pevt->corey[whichtel],pevt->height.at(ic),pevt->oheight);

      bool inside=true;
      if(tlimit[1]>tlimit[0]){
         double dxdt=(coo2[0]-coo1[0])/(t2-t1);
         double dydt=(coo2[1]-coo1[1])/(t2-t1);
         double dzdt=(coo2[2]-coo1[2])/(t2-t1);
         double t1new=(tlimit[0]<t1)?t1:tlimit[0]; //maximum between t1 and tlimit[0]
         double t2new=(tlimit[1]>t2)?t2:tlimit[1]; //minimum between t2 and tlimit[1]
         if(t1!=t1new){
            coo1[0]=coo1[0]+dxdt*(t1new-t1);
            coo1[1]=coo1[1]+dydt*(t1new-t1);
            coo1[2]=coo1[2]+dzdt*(t1new-t1);
         }
         if(t2!=t2new){
            coo2[0]=coo2[0]+dxdt*(t2new-t2);
            coo2[1]=coo2[1]+dydt*(t2new-t2);
            coo2[2]=coo2[2]+dzdt*(t2new-t2);
         }
         t1=t1new;
         t2=t2new;
         if(t2<=t1) inside=false;
      }

      if(particle>=0) inside=inside&&(particle&GetParticleID(partid));
      if(elimit[1]>elimit[0]) inside=inside&&(energy>=elimit[0]&&energy<=elimit[1]);
      if(tlimit[1]>tlimit[0]){
         inside=inside&&(t1>=tlimit[0]&&t1<=tlimit[1]);
         inside=inside&&(t2>=tlimit[0]&&t2<=tlimit[1]);
      }
      for(int ii=0;ii<3;ii++){
         if(climit[ii][1]>climit[ii][0]){
            inside=inside&&(coo1[ii]>=climit[ii][0]&&coo1[ii]<=climit[ii][1]);
            inside=inside&&(coo2[ii]>=climit[ii][0]&&coo2[ii]<=climit[ii][1]);
         }
      }
      if(inside){
         arrid.push_back(partid);
         arren.push_back(energy);
         arrt1.push_back(t1);
         arrt2.push_back(t2);
         if(t1<plotrange[3][0]) plotrange[3][0]=t1;
         if(t2<plotrange[3][0]) plotrange[3][0]=t2;
         if(t1>plotrange[3][1]) plotrange[3][1]=t1;
         if(t2>plotrange[3][1]) plotrange[3][1]=t2;
         for(int ii=0;ii<3;ii++){
            arrc1[ii].push_back(coo1[ii]);
            arrc2[ii].push_back(coo2[ii]);
            if(IniRange[ii][0]>=IniRange[ii][1]){
            if(coo1[ii]<plotrange[ii][0]) plotrange[ii][0]=coo1[ii];
            if(coo2[ii]<plotrange[ii][0]) plotrange[ii][0]=coo2[ii];
            if(coo1[ii]>plotrange[ii][1]) plotrange[ii][1]=coo1[ii];
            if(coo2[ii]>plotrange[ii][1]) plotrange[ii][1]=coo2[ii];
            }
         }
         if(jdebug>1) printf("ReadTrack::Copy: fill the cer track nrec=%d partid=%d range={{%f,%f},{%f,%f},{%f,%f}}\n",ic,partid,plotrange[0][0],plotrange[0][1],plotrange[1][0],plotrange[1][1],plotrange[2][0],plotrange[2][1]);
         if(jdebug>3){
            printf("ReadTrack::Copy: Content partid=%d energy=%f start={%f,%f,%f,%f} end={%f,%f,%f,%f}\n",partid,energy,coo1[0],coo1[1],coo1[2],t1,coo2[0],coo2[1],coo2[2],t2);
         }
      }
   }
}

int ReadTrack::Color(int partid){
   if(partid<0) return 1;
   else if(partid==0) return 6;
   else if(partid<=3) return 2;
   else if((partid>=13&&partid<=15)||(partid>=116&&partid<=128)||partid>=200) return 4;
   else return 3;

   //if(partid<0) return 1;
   //else if(partid==0) return 2;
   //else if(partid==1) return 3;
   //else if(partid==2) return 4;
   //else if(partid==3) return 6;
   //else return 1;
}
int ReadTrack::GetLegType(int color){ //Reverse of ReadTrack::Color
   if(color==1) return -1;
   else if(color==6) return 3;
   else if(color==2) return 0;
   else if(color==4) return 2;
   else if(color==3) return 1;
   else return -1;
}
int ReadTrack::Style(int partid){
   return 1;
}
int ReadTrack::Width(int partid){
   if(partid==0) return 3;
   if(partid<=9) return 1;
   return 2;
}
void ReadTrack::Draw(TCanvas* cc,const char* option){
   //if(!cc){
   //   cc = new TCanvas("Air Shower");
   //   TView *view = TView::CreateView(2);
   //   view->SetRange(plotrange[0][0]-1,plotrange[1][0]-1,plotrange[2][0]-1,plotrange[0][1]+1,plotrange[1][1]+1,plotrange[2][1]+1);
   //   cc->SetView(view);
   //}
   if(plot) {
      //plot->Delete();
      delete plot;
   }
   plot=new TObjArray();
   for(int ii=0;ii<arrid.size();ii++){
      TPolyLine3D* line = new TPolyLine3D(2);
      line->SetPoint(0, arrc1[0].at(ii), arrc1[1].at(ii), arrc1[2].at(ii));
      line->SetPoint(1, arrc2[0].at(ii), arrc2[1].at(ii), arrc2[2].at(ii));
      line->SetLineColor(Color(arrid.at(ii)));
      //line->SetLineColor(type);
      line->SetLineStyle(Style(arrid.at(ii)));
      line->SetLineWidth(Width(arrid.at(ii)));
      plot->Add(line);
   }
   if(cc){
      plot->Draw(option);
   }
}
