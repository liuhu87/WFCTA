#include "ReadTrack.h"
#include "CorsikaEvent.h"
ReadTrack* ReadTrack::_Head=0;
bool ReadTrack::DoPlot=false;
int ReadTrack::jdebug=0;
int ReadTrack::particle=-1;
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
      return;
   }
}
bool ReadTrack::Exist(){
   return (bool)fin;
}
int ReadTrack::ReadRec(){
   if(!fin) return -3;
   union{
      int i[2];
      char c[8];
   } padding;
   if(nrec==0) fin->read(padding.c,4);
   else fin->read(padding.c,8);
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
   if(padding.i[0]!=40){
      printf("ReadTrack::ReadRec: Error of the record size(%d,%f)\n",padding.i[0],padding.i[0]);
      return -1;
   }
   else fin->read(recbuff.c,(int)padding.i[0]);
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
   if(particle>=0) inside=inside&&(partid==particle);
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
         if(plotrange[ii][0]>=plotrange[ii][1]){
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
      plot->Delete();
      delete plot;
   }
   plot=new TObjArray();
   int size=pevt->cx.size();
   for(int ic=0;ic<size;ic++){
      int partid=0;
      double energy=(pevt->wavelength.at(ic)>0)?(hplank_gev*vlight/(1.0e-7*pevt->wavelength.at(ic))):1.0e8;
      double coo1[3];
      double coo2[3];
      double t1,t2;
      t2=pevt->ct.at(ic);
      coo2[0]=pevt->cx.at(ic);
      coo2[1]=pevt->cy.at(ic);
      coo2[2]=pevt->oheight;
      coo1[2]=pevt->height.at(ic);
      double dxdr=-pevt->cu.at(ic);
      double dydr=-pevt->cv.at(ic);
      double dzdr=(1-dxdr*dxdr-dydr*dydr>=0)?-sqrt(1-dxdr*dxdr-dydr*dydr):-1;
      coo1[0]=(coo1[2]-coo2[2])/dzdr*dxdr+coo2[0];
      coo1[1]=(coo1[2]-coo2[2])/dzdr*dydr+coo2[1];
      t1=t2-sqrt(pow(coo1[0]-coo2[0],2)+pow(coo1[1]-coo2[1],2)+pow(coo1[2]-coo2[2],2))/vlight;
      bool inside=true;
      if(particle>=0) inside=inside&&(partid==particle);
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
            if(plotrange[ii][0]>=plotrange[ii][1]){
            if(coo1[ii]<plotrange[ii][0]) plotrange[ii][0]=coo1[ii];
            if(coo2[ii]<plotrange[ii][0]) plotrange[ii][0]=coo2[ii];
            if(coo1[ii]>plotrange[ii][1]) plotrange[ii][1]=coo1[ii];
            if(coo2[ii]>plotrange[ii][1]) plotrange[ii][1]=coo2[ii];
            }
         }
         if(jdebug>1) printf("ReadTrack::Copy: fill the cer track nrec=%d partid=%d\n",ic,partid);
      }
   }
}

int ReadTrack::Color(int partid){
   if(partid<0) return 1;
   else if(partid==0) return 6;
   else if(partid<=3) return 2;
   else if(partid<=3) return 2;
   else if(partid<=9) return 3;
   else return 4;
}
int ReadTrack::Style(int partid){
   return 1;
}
int ReadTrack::Width(int partid){
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
      plot->Delete();
      delete plot;
   }
   plot=new TObjArray();
   for(int ii=0;ii<arrid.size();ii++){
      TPolyLine3D* line = new TPolyLine3D(2);
      line->SetPoint(0, arrc1[0].at(ii), arrc1[1].at(ii), arrc1[2].at(ii));
      line->SetPoint(1, arrc2[0].at(ii), arrc2[1].at(ii), arrc2[2].at(ii));
      line->SetLineColor(Color(arrid.at(ii)));
      line->SetLineStyle(Style(arrid.at(ii)));
      line->SetLineWidth(Width(arrid.at(ii)));
      plot->Add(line);
   }
   if(cc){
      plot->Draw(option);
   }
}
TH1F* ReadTrack::GetN_vs_Z(int partid){
   int size=arrid.size();
   if(size<=0) return 0;
   TH1F* hist=new TH1F("N_vs_Z","",1000,plotrange[2][0],plotrange[2][1]);
   for(int ii=0;ii<size;ii++){
      if(partid>=0&&partid!=arrid.at(ii)) continue;
      double weight=1.;
      double z1=arrc1[2].at(ii);
      double z2=arrc2[2].at(ii);
      int bin1=hist->GetXaxis()->FindBin(z1);
      int bin2=hist->GetXaxis()->FindBin(z2);
      for(int ibin=bin1;ibin<=bin2;ibin++){
         double center=hist->GetXaxis()->GetBinCenter(ibin);
         double binwidth=hist->GetXaxis()->GetBinWidth(ibin);
         if(ibin>bin1&&ibin<bin2) hist->Fill(center,weight);
         else if(ibin==bin1) hist->Fill(center,weight*(hist->GetXaxis()->GetBinUpEdge(ibin)-z1)/binwidth);
         else if(ibin==bin2) hist->Fill(center,weight*(z2-hist->GetXaxis()->GetBinLowEdge(ibin))/binwidth);
      }
   }
   return hist;
}
TH1F* ReadTrack::GetL_vs_Z(int partid){
   int size=arrid.size();
   if(size<=0) return 0;
   TH1F* hist=new TH1F("L_vs_Z","",1000,plotrange[2][0],plotrange[2][1]);
   for(int ibin=1;ibin<=hist->GetNbinsX();ibin++){
      double center=hist->GetXaxis()->GetBinCenter(ibin);
      double binwidth=hist->GetXaxis()->GetBinWidth(ibin);
      int np=0;
      double xsum1=0,ysum1=0;
      double xsum2=0,ysum2=0;
      for(int ii=0;ii<size;ii++){
         if(partid>=0&&partid!=arrid.at(ii)) continue;
         double weight=1.;
         double z1=arrc1[2].at(ii);
         double z2=arrc2[2].at(ii);
         if(center<z1||center>z2) continue;
         double x1=arrc1[0].at(ii);
         double x2=arrc2[0].at(ii);
         double y1=arrc1[1].at(ii);
         double y2=arrc2[1].at(ii);
         double xx=x1+(x2-x1)/(z2-z1)*(center-z1);
         double yy=y1+(y2-y1)/(z2-z1)*(center-z1);
         xsum1+=xx;
         xsum2+=xx*xx;
         ysum1+=yy;
         ysum2+=yy*yy;
         np++;
      }
      if(np>0) hist->SetBinContent(ibin,sqrt(xsum2/np-pow(xsum1/np,2)+ysum2/np-pow(ysum1/np,2)));
   }
   return hist;
}

