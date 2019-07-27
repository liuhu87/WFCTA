#include "ShowerPlot.h"
char ShowerPlot::tpname[NShowerTrack][20]={"electrons","muons","hadrons","cer light"};
void ShowerPlot::Init(const char* priname){
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
   if(plot) plot->Delete();
   for(int ii=0;ii<NShowerTrack;ii++){
      if(!ptrk[ii]) continue;
      if(ii<3){
         ptrk[ii]->Reset();
         ptrk[ii]->ReadAll();
      }
      else{
         ReadTrack::SetHead(ptrk[ii]);
         corsika->Reset();
         corsika->ReadAll(0,0,pevt);
      }
      ptrk[ii]->Draw();
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
void ShowerPlot::Draw(){
   TCanvas* cc = new TCanvas("Air Shower");
   TView *view = TView::CreateView(1);

   TLegend* leg=new TLegend(0.7,0.7,0.9,0.9);
   leg->SetHeader(Form("Primary: %s",primary));

   float plotrange[4][2]={{1.0e10,-1.0e10},{1.0e10,-1.0e10},{1.0e10,-1.0e10}};
   for(int ii=0;ii<4;ii++){
      if(!ptrk[ii]) continue;
      if(!ptrk[ii]->plot) continue;
      if(ptrk[ii]->plot->GetEntries()<=0) continue;
      for(int ia=0;ia<3;ia++){
         if(ptrk[ii]->plotrange[ia][0]<plotrange[ia][0]) plotrange[ia][0]=ptrk[ii]->plotrange[ia][0];
         if(ptrk[ii]->plotrange[ia][1]>plotrange[ia][0]) plotrange[ia][1]=ptrk[ii]->plotrange[ia][1];
      }
      leg->AddEntry(ptrk[ii]->plot->At(0),tpname[ii],"l");
   }
   view->SetRange(plotrange[0][0]-1,plotrange[1][0]-1,plotrange[2][0]-1,plotrange[0][1]+1,plotrange[1][1]+1,plotrange[2][1]+1);
   view->ShowAxis();
   cc->SetView(view);

   if(plot) plot->Draw();
   leg->Draw("same");

}

