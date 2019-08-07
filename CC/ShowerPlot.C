#include "ShowerPlot.h"
#include "TAxis3D.h"
int ShowerPlot::jdebug=0;
char ShowerPlot::tpname[NShowerTrack][20]={"electrons","muons","hadrons","cer light"};
void ShowerPlot::Init(const char* priname){
   for(int ii=0;ii<4;ii++){
      plotrange[ii][0]=0;
      plotrange[ii][0]=0;
   }
   if(priname) strcpy(primary,priname);
   for(int ii=0;ii<NShowerTrack;ii++){
      strcpy(filename[ii],"");
      ptrk[ii]=0;
   }
   corsika=0;
   pevt=0;
   plot=0;
}
void ShowerPlot::Release(){
   for(int ii=0;ii<NShowerTrack;ii++) if(ptrk[ii]) delete ptrk[ii];
   if(corsika) delete corsika;
   if(pevt) delete pevt;
   if(plot){
      plot->Delete();
      delete plot;
   }
}

bool ShowerPlot::Add(const char* inputfile,int type){
   if(type<0||type>=NShowerTrack) return false;
   if(!inputfile) return false;
   strcpy(filename[type],inputfile);
   if(type<3){
      ptrk[type]=new ReadTrack(filename[type]);
      if(!ptrk[type]->Exist()){
         delete ptrk[type];
         ptrk[type]=0;
         printf("Error in create ReadTrack filename=%s\n",inputfile);
         return false;
      }
   }
   else{
      corsika=new CorsikaIO(inputfile,0);
      if(!corsika->fin){
         delete corsika;
         corsika=0;
         printf("Error in create CorsikaIO filename=%s\n",inputfile);
         return false;
      }
      pevt=new CorsikaEvent();
      ptrk[type]=new ReadTrack(0);
   }
   return true;
}
bool ShowerPlot::Read(){
   bool result=false;
   if(jdebug>3) printf("ShowerPlot::Read: Delete plot\n");
   if(plot){
      //plot->Delete();
      delete plot;
      plot=0;
   }
   for(int ii=0;ii<NShowerTrack;ii++){
      if(!ptrk[ii]) continue;
      if(ii<3){
         if(jdebug>2) printf("ShowerPlot::Read: Read %d begin\n",ii);
         ptrk[ii]->Reset();
         if(jdebug>2) printf("ShowerPlot::Read: Read %d reset finished nrec=%d\n",ii,ptrk[ii]->nrec);
         int read=ptrk[ii]->ReadAll();
         if(jdebug>2) printf("ShowerPlot::Read: Read %d read finished return=%d\n",ii,read);
      }
      else{
         if(jdebug>2) printf("ShowerPlot::Read: Read %d begin\n",ii);
         ReadTrack::SetHead(ptrk[ii]);
         corsika->Reset();
         if(jdebug>2) printf("ShowerPlot::Read: Read %d reset finished nrec=%d\n",ii,corsika->nrec);
         int read=corsika->ReadAll(0,0,pevt);
         if(jdebug>2) printf("ShowerPlot::Read: Read %d read finished return=%d\n",ii,read);
      }
      ptrk[ii]->Draw();
      if(jdebug>1) printf("ShowerPlot::Read: Read %d read %d tracks\n",ii,ptrk[ii]->plot?ptrk[ii]->plot->GetEntries():0);
      if(!ptrk[ii]->plot){
         printf("Error in read filename=%s\n",filename[ii]);
         continue;
      }
      TObjArray* buff=ptrk[ii]->plot;
      if(buff->GetEntries()>0) result=true;
      if(!plot){
         plot=new TObjArray(*buff);
      }
      else{
         for(int iel=0;iel<buff->GetEntries();iel++){
            plot->Add(buff->At(iel));
         }
      }
   }
   return result;
}
TCanvas* ShowerPlot::Draw(int ViewOpt){
   TCanvas* cc = new TCanvas("Air Shower",(ViewOpt<1||ViewOpt>3)?"Inclined View":(ViewOpt==1?"Front View":(ViewOpt==2?"Side View":"Top View")),ViewOpt==3?3000:2000,3000);

   TLegend* leg_time=0;
   if(ReadTrack::tlimit[1]>ReadTrack::tlimit[0]&&ReadTrack::tlimit[1]>=0){
      leg_time=new TLegend(0.4,0.8,0.6,0.9);
      leg_time->SetHeader(Form("Time=%6.1es",ReadTrack::tlimit[1]));
      leg_time->SetLineColorAlpha(kWhite,1.0);
      leg_time->SetTextFont(62);
      leg_time->SetTextSize(0.03);
      leg_time->SetTextColor(4);
   }
   TLegend* leg=new TLegend(0.7,0.7,0.9,0.9);
   leg->SetHeader(Form("Primary: %s",primary));
   leg->SetTextSize(0.02);

   for(int ii=0;ii<4;ii++){
      plotrange[ii][0]=5.0e9;
      plotrange[ii][1]=-5.0e9;
   }
   for(int ii=0;ii<4;ii++){
      if(!ptrk[ii]) continue;
      if(!ptrk[ii]->plot) continue;
      if(ptrk[ii]->plot->GetEntries()<=0) continue;
      for(int ia=0;ia<4;ia++){
         if(ptrk[ii]->plotrange[ia][0]<plotrange[ia][0]) plotrange[ia][0]=ptrk[ii]->plotrange[ia][0];
         if(ptrk[ii]->plotrange[ia][1]>plotrange[ia][1]) plotrange[ia][1]=ptrk[ii]->plotrange[ia][1];
         if(jdebug>10) printf("ShowerPlot::Draw rminmax: ii=%d ic=%d min={%f,%f} max={%f,%f}\n",ii,ia,ptrk[ii]->plotrange[ia][0],plotrange[ia][0],ptrk[ii]->plotrange[ia][1],plotrange[ia][1]);
      }
      leg->AddEntry(ptrk[ii]->plot->At(0),tpname[ii],"l");
   }
   double rmin[3]={plotrange[0][0],plotrange[1][0],plotrange[2][0]};
   double rmax[3]={plotrange[0][1],plotrange[1][1],plotrange[2][1]};
   TView *view = TView::CreateView(1,rmin,rmax);
   view->ShowAxis();
   if(ViewOpt==1) view->Front();
   else if(ViewOpt==2) view->Side();
   else if(ViewOpt==3) view->Top();
   cc->SetView(view);
   if(jdebug>0) printf("ShowerPlot::Draw: Pos Range={{%+6.1e,%+6.1e},{%+6.1e,%+6.1e},{%+6.1e,%+6.1e}} Time Range={%+7.1e,%+7.1e}\n",rmin[0],rmax[0],rmin[1],rmax[1],rmin[2],rmax[2],plotrange[3][0],plotrange[3][1]);

   if(plot){
      plot->Draw();
      TAxis3D *axis = TAxis3D::GetPadAxis();
      axis->SetLabelColor(kBlack);
      axis->SetAxisColor(kBlack);
      axis->SetLabelSize(0.03,"X");
      axis->SetLabelSize(0.03,"Y");
      axis->SetLabelSize(0.03,"Z");
      axis->SetXTitle("X [cm]");
      axis->SetYTitle("Y [cm]");
      axis->SetZTitle("Z [cm]");
      axis->SetTitleOffset(2.,"X");
      axis->SetTitleOffset(2.,"Y");
      axis->SetTitleOffset(2.,"Z");
      if(ViewOpt<1||ViewOpt>3){
      axis->SetLabelOffset(0.006,"X");
      axis->SetLabelOffset(0.008,"Y");
      axis->SetLabelOffset(0.008,"Z");
      axis->SetTitleOffset(1.2,"X");
      axis->SetTitleOffset(2.0,"Y");
      axis->SetTitleOffset(1.2,"Z");
      }
      else if(ViewOpt==1){
      axis->SetLabelOffset(-0.08,"X");
      axis->SetLabelOffset(0.03,"Z");
      axis->SetTitleOffset(1.3,"X");
      }
      else if(ViewOpt==2){
      axis->SetLabelOffset(-0.08,"Y");
      axis->SetLabelOffset(0.03,"Z");
      axis->SetTitleOffset(1.3,"Y");
      }
      else if(ViewOpt==3){
      axis->SetLabelOffset(0.03,"X");
      axis->SetLabelOffset(0.03,"Y");
      }
   }
   leg->Draw("same");
   if(leg_time) leg_time->Draw("same");
   return cc;
}

