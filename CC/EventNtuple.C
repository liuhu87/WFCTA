#include "EventNtuple.h"
#include "TMath.h"
EventNtuple* EventNtuple::_Head=0;
TFile* EventNtuple::fout=0;
int EventNtuple::fillstyle=0;
const int EventNtuple::branchSplit=1;
char EventNtuple::branchname[NBranch][20]={"CorsikaEvent","Event"};
bool EventNtuple::FromIDtoAZ(int PID,int &AA,int &ZZ){
   int nele=11;
   int pid0[]={1,2,3, 5, 6,7,8, 9,13,14,15};
   int A0[]  ={0,0,0, 0, 0,0,0, 0, 1, 1, 1};
   int Z0[]  ={0,1,-1,1,-1,0,1,-1, 0, 1,-1};
   if(PID<=15){
      bool exist=false;
      for(int ii=0;ii<nele;ii++){
         if(PID==pid0[ii]){
            AA=A0[ii];
            ZZ=Z0[ii];
            return true;
         }
      }
      if(!exist) {AA=-1;ZZ=0;return false;}
   }
   else if(PID<200) {AA=-1;ZZ=0;return false;}
   else{
      AA=PID/100;
      ZZ=PID%100;
      return true;
   }
}
EventNtuple* EventNtuple::GetHead(){
   return _Head;
}
void EventNtuple::SetHead(EventNtuple* head){
   _Head=head;
}
EventNtuple* EventNtuple::GetHead(char* filename,int style){
   if(_Head) return _Head;
   else{
      _Head=new EventNtuple(filename,style);
      return _Head;
   }
}
void EventNtuple::Init(){
   if(fillstyle==1||fillstyle==0){//init the tree
      //if(_Tree) delete _Tree;
      _Tree=new TTree("eventShow","Event");
   }
   else{
      _Tree=0;
   }
   if((fillstyle==2||fillstyle==3)||fillstyle==0){//init the histograms
      hlist=new TList();
      if(fillstyle==2||fillstyle==0){
         double PEL=1.0e3,PEH=5.0e5;
         double SPL[2]={-3.2e4,-1.0e5},SPH[2]={3.2e4,1.0e5};
         double SNL[2]={1.0e2,1.0e1},SNH[2]={5.0e6,1.0e5};
         //some primary particle distribution
         ///primary particle type distribution
         TH1F* Hptp=new TH1F("Hptp","Primary particle nucleon number distribution;Atomic Number;Events",102,-1.5,100.5);
         hlist->Add(Hptp);
         ///primary particle energy distribution
         TH1F* Hpe=new TH1F("Hpe","Primary particle energy distribution;Energy [GeV];Events",1000,PEL,PEH);
         hlist->Add(Hpe);
         ///primary particle direction distribution
         TH2F* Hpdir=new TH2F("Hpdir","Primary particle direction distribution;cos^{2}(Theta);Phi;Events",100,0,1.,100,0,2*PI);
         hlist->Add(Hpdir);
         ///primary particle position distribution
         TH2F* Hppos=new TH2F("Hppos","Primary particle position distribution;x [cm];y [cm]",100,SPL[0],SPH[0],100,SPL[0],SPH[0]);
         hlist->Add(Hppos);
         //some secondary particle distribution
         TH2F *Hspos[2],*Hsdir[2],*Hpse[2],*Hpspos[2],*Hpsdir[2],*Hpstp[2];
         for(int ii=0;ii<2;ii++){
            ///position distribution of clight and particle at observation level
            Hspos[ii]=new TH2F(Form("Hspos_s%d",ii),"Secondary particle position distribution;x [cm];y[cm]",100,SPL[ii],SPH[ii],100,SPL[ii],SPH[ii]);
            hlist->Add(Hspos[ii]);
            ///arrival direction distribution of clight and particle at observation level
            Hsdir[ii]=new TH2F(Form("Hsdir_s%d",ii),"Secondary particle direction distribution;cos(Theta);Phi",100,-1,1,100,0,2*PI);
            hlist->Add(Hsdir[ii]);
            ///number of clight or electrons vs primary particle energy
            Hpse[ii]=new TH2F(Form("Hpse_s%d",ii),"Number of secondary particle vs primary particle energy;primary particle energy [GeV];log(secondary particle number)",1000,PEL,PEH,100,log10(SNL[ii]),log10(SNH[ii]));
            hlist->Add(Hpse[ii]);
            ///number of clight or electrons vs primary particle position
            Hpspos[ii]=new TH2F(Form("Hpspos_s%d",ii),"Number of secondary particle vs primary core position;primary particle dist;log(secondary particle number)",100,SPL[ii],SPH[ii],100,log10(SNL[ii]),log10(SNH[ii]));
            hlist->Add(Hpspos[ii]);
            ///number of clight or electrons vs primary particle direction
            Hpsdir[ii]=new TH2F(Form("Hpsdir_s%d",ii),"Number of secondary particle vs primary particle direction;primary particle theta;log(secondary particle number)",100,-1.,1.,100,log10(SNL[ii]),log10(SNH[ii]));
            hlist->Add(Hpsdir[ii]);
            ///number of clight or electrons vs primary particle type
            Hpstp[ii]=new TH2F(Form("Hpstp_s%d",ii),"Number of secondary particle vs primary particle type;primary particle atomic number;log(secondary particle number)",102,-1.5,100.5,100,log10(SNL[ii]),log10(SNH[ii]));
            hlist->Add(Hpstp[ii]);
         }

         ///position of secondary with respect to telescope
         TH2F* Hstpos=new TH2F("Hstpos","Secondary particle position distribution;x [cm];y[cm]",100,-500.,500.,100,-500.,500.);
         hlist->Add(Hstpos);

         //WFCTA simulation result
         const int nbin=13;
         char label[nbin][100]={"All","InTel","RefEff","Door","OutCluster","ZMirror","InMirror","InCluster","Filter","InCone","ConeStrike","InSiPM","Signal"};
         ///number of clight when pass through different geometry
         TH1F* Hsn=new TH1F("Hsn",";Stage;",nbin,-0.5,nbin-0.5);
         hlist->Add(Hsn);
         ///number of clight when pass through different geometry vs energy
         TH2F* Hsne=new TH2F("Hsne",";Primary Particle Energy;Stage",1000,PEL,PEH,nbin,-0.5,nbin-0.5);
         hlist->Add(Hsne);
         ///number of clight when pass through different geometry vs position
         TH2F* Hsnposx=new TH2F("Hsnposx",";position x;Stage",1000,SPL[0],SPH[0],nbin,-0.5,nbin-0.5);
         hlist->Add(Hsnposx);
         TH2F* Hsnposy=new TH2F("Hsnposy",";position x;Stage",1000,SPL[0],SPH[0],nbin,-0.5,nbin-0.5);
         hlist->Add(Hsnposy);
         ///number of clight when pass through different geometry vs direction
         TH2F* Hsndiru=new TH2F("Hsndiru",";;Stage",100,-1,1,nbin,-0.5,nbin-0.5);
         hlist->Add(Hsndiru);
         TH2F* Hsndirv=new TH2F("Hsndirv",";;Stage",100,0,2*PI,nbin,-0.5,nbin-0.5);
         hlist->Add(Hsndirv);
         for(int ibin=1;ibin<=nbin;ibin++){
            Hsn->GetXaxis()->SetBinLabel(ibin,label[ibin-1]);
            Hsne->GetYaxis()->SetBinLabel(ibin,label[ibin-1]);
            Hsnposx->GetYaxis()->SetBinLabel(ibin,label[ibin-1]);
            Hsnposy->GetYaxis()->SetBinLabel(ibin,label[ibin-1]);
            Hsndiru->GetYaxis()->SetBinLabel(ibin,label[ibin-1]);
            Hsndirv->GetYaxis()->SetBinLabel(ibin,label[ibin-1]);
         }
         TH2F* Hsimg[20];
         for(int ict=0;ict<20;ict++){
            if(ict<WFTelescopeArray::CTNumber){
               ///image of fired pmt for each telescope
               Hsimg[ict]=new TH2F(Form("Hsimg_T%d",ict),Form("Telescope %d;x index;y index",ict),32,-0.5,31.5,32,-0.5,31.5);
               hlist->Add(Hsimg[ict]);
            }
            else Hsimg[ict]=0;
         }
      }
      else if(fillstyle==3||fillstyle==0){
         double PEL=1.0e3,PEH=5.0e5;
         hlist=new TList();
         ///primary particle energy distribution
         TH1F* MCHpe=new TH1F("MCHpe","Primary particle energy distribution;Energy [GeV];Events",100,PEL,PEH);
         hlist->Add(MCHpe);
      }
   }
   else hlist=0;
}
EventNtuple::EventNtuple(char* filename,int style){
   if(style>0) fillstyle=style;
   else fillstyle=-1;
   Init();
   if(!filename){
      printf("There is no output filename. Exiting\n");
      return;
   }
   else{
      fout=TFile::Open(filename,"RECREATE");
      fout->cd();
      TH1::SetDefaultSumw2();
      return;
   }
}
void EventNtuple::Clear(){
   if(_Tree) delete _Tree;
   if(hlist){
      for(int ii=0;ii<hlist->GetEntries();ii++){
         if(hlist->At(ii)) delete hlist->At(ii);
      }
      delete hlist;
      hlist=0;
   }
}
void EventNtuple::Write(){
   if(fout){
      fout->cd();
      //fout->Write();
      if(_Tree) _Tree->Write();
      if(hlist){
         for(int ii=0;ii<hlist->GetEntries();ii++){
            if(hlist->At(ii)) hlist->At(ii)->Write();
         }
      }
      fout->Close();
   }
}
EventNtuple::~EventNtuple(){
   Write();
   Clear();
}
void EventNtuple::ConvertCoor(double pos[3],double theta,double phi,double InCoo[3],double InDir[3],double *OutCoo,double *OutDir){
   double matrix_[3][3];
   double cosz, sinz, cosa, sina;
   cosz = cos(theta);
   sinz = sin(theta);
   cosa = cos(phi);
   sina = sin(phi);
   matrix_[0][0] = cosa*cosz;
   matrix_[0][1] = sina*cosz;
   matrix_[0][2] = -sinz;
   matrix_[1][0] = -sina;
   matrix_[1][1] = cosa;
   matrix_[1][2] = 0;
   matrix_[2][0] = cosa*sinz;
   matrix_[2][1] = sina*sinz;
   matrix_[2][2] = cosz;

   for(int ii=0;ii<3;ii++){
      OutCoo[ii]=0;
      OutDir[ii]=0;
      for(int jj=0;jj<3;jj++){
         OutCoo[ii]+=matrix_[ii][jj]*(InCoo[jj]-pos[jj]);
         OutDir[ii]+=matrix_[ii][jj]*InDir[jj];
      }
   }
}
void EventNtuple::ConvertThetaPhi(double pp[3],double &theta,double &phi){
   double norm=sqrt(pow(pp[0],2)+pow(pp[1],2)+pow(pp[2],2));
   if(norm<=0) {theta=0; phi=0; return;}
   theta=acos(pp[2]/norm);
   if(pp[0]==0){
      phi=(pp[1]>=0)?(PI/2):(3*PI/2);
   }
   else{
      phi=atan(pp[1]/pp[0]);
      if(pp[0]<0) phi+=PI;
      if(phi<0) phi+=(2*PI);
   }
}
void EventNtuple::ConvertCoor(double pos[3],double pp[3],double InCoo[3],double InDir[3],double *OutCoo,double *OutDir){
   double theta,phi;
   ConvertThetaPhi(pp,theta,phi);
   ConvertCoor(pos,theta,phi,InCoo,InDir,OutCoo,OutDir);
}
void EventNtuple::Fill(TSelector** pevt){
   if(!pevt) return;
   if(fillstyle==1||fillstyle==0){//fill the tree
      if(!_Tree) return;
      bool exist=false;
      for(int ibr=0;ibr<NBranch;ibr++){
         if(!pevt[ibr]) continue;
         TBranch* branch=_Tree->GetBranch(branchname[ibr]);
         if(!branch){
            printf("EventNtuple::Fill: No branch %s, newly created\n",branchname[ibr]);
            if(ibr==0) _Tree->Branch(branchname[ibr],branchname[ibr],((CorsikaEvent*)pevt[ibr]),128000,branchSplit);
            if(ibr==1) _Tree->Branch(branchname[ibr],branchname[ibr],((WFCTAEvent*)pevt[ibr]),128000,branchSplit);
            branch=_Tree->GetBranch(branchname[ibr]);
            if(branch){
               int clevel=branch->GetCompressionLevel();
               #ifdef __LZMA__
               if(clevel<6 && branch->GetCompressionAlgorithm()!=ROOT::kLZMA)clevel=6;
               #else
               if(clevel<6)clevel=6;
               #endif
               branch->SetCompressionLevel(clevel);
               printf("EventNtuple::Fill Tree CompressionLevel %d %d\n",branch->GetCompressionLevel(),branch->GetSplitLevel());
               if(ibr==0) _Tree->SetBranchStatus(Form("TSelector"),false);
               exist=true;
            }
         }
         else exist=true;
         if(exist&&false){
            int ibin1=((WFCTAEvent*)pevt[ibr])->GetMaxADCBin();
            int ibin2=((WFCTAEvent*)pevt[ibr])->GetMinTimeBin();
            printf("EventNtuple::Fill: event=%ld rabbitTime=%d bin={%d,%d} adc={%f,%f} time=%.9e\n",((WFCTAEvent*)pevt[ibr])->iEvent,((WFCTAEvent*)pevt[ibr])->rabbitTime,ibin1,ibin2,((WFCTAEvent*)pevt[ibr])->ImageAdcHigh.at(ibin1>=0?ibin1:0),((WFCTAEvent*)pevt[ibr])->ImageAdcLow.at(ibin1>=0?ibin1:0),((WFCTAEvent*)pevt[ibr])->mcevent.ArrivalTimeMin[0][ibin2>=0?ibin2:0]);
         }
      }
      if(exist) _Tree->Fill();
   }
   if(fillstyle==2||fillstyle==0){//fill the histograms
      for(int ibr=0;ibr<NBranch;ibr++){
         if(ibr==0&&pevt[ibr]){
            CorsikaEvent* padd=(CorsikaEvent*)pevt[ibr];
            int AA,ZZ;
            EventNtuple::FromIDtoAZ(padd->partp,AA,ZZ);
            ((TH1F*)hlist->FindObject(Form("Hptp")))->Fill(AA);
            ((TH1F*)hlist->FindObject(Form("Hpe")))->Fill(padd->ep);
            ((TH2F*)hlist->FindObject(Form("Hpdir")))->Fill(pow(TMath::Cos(padd->thetap),2),padd->phip<0?(padd->phip+2*PI):padd->phip);
            for(int iuse=0;iuse<Nuse;iuse++) ((TH2F*)hlist->FindObject(Form("Hppos")))->Fill(padd->corex[iuse],padd->corey[iuse]);
            int nclight=(padd->cx).size();
            int npart=(padd->px).size();
            for(int ii=0;ii<2;ii++){
               if(ii==0){
                  int ntot=0;
                  int ntotcore[Nuse];
                  for(int i0=0;i0<Nuse;i0++) ntotcore[i0]=0;
                  for(int jj=0;jj<nclight;jj++){
                     double cosx=(padd->cu).at(jj);
                     double cosy=(padd->cv).at(jj);
                     double cosz=(cosx*cosx+cosy*cosy<=1)?-sqrt(1-cosx*cosx-cosy*cosy):-1;
                     double pos[3]={0,0,0};
                     double InCoo[3]={(padd->cx).at(jj),(padd->cy).at(jj),0},InDir[3]={cosx,cosy,cosz};
                     double OutCoo[3],OutDir[3];
                     ConvertCoor(pos,padd->thetap,padd->phip,InCoo,InDir,OutCoo,OutDir);
                     double OutTheta,OutPhi;
                     ConvertThetaPhi(OutDir,OutTheta,OutPhi);
                     //if(pow(TMath::Cos(padd->thetap),2)>0.9){
                     ((TH2F*)hlist->FindObject(Form("Hspos_s%d",ii)))->Fill((padd->cx).at(jj),(padd->cy).at(jj),(padd->cweight).at(jj));
                     ((TH2F*)hlist->FindObject(Form("Hsdir_s%d",ii)))->Fill(cos(OutTheta),OutPhi,(padd->cweight).at(jj));
                     //}
                     ntot+=1.*(padd->cweight).at(jj);
                     //printf("EveNtuple::Fill: begin core cal\n");
                     int whichcore=padd->WhichCore((padd->cx).at(jj),(padd->cy).at(jj));
                     ntotcore[whichcore]+=1.*(padd->cweight).at(jj);
                     ((TH2F*)hlist->FindObject(Form("Hstpos")))->Fill((padd->cx).at(jj)-padd->corex[whichcore],(padd->cy).at(jj)-padd->corey[whichcore],(padd->cweight).at(jj));
                  }
                  ((TH2F*)hlist->FindObject(Form("Hpse_s%d",ii)))->Fill(padd->ep,ntot>0?log10(ntot):0);
                  for(int icore=0;icore<Nuse;icore++) ((TH2F*)hlist->FindObject(Form("Hpspos_s%d",ii)))->Fill(sqrt(pow(padd->corex[icore],2)+pow(padd->corey[icore],2)),ntotcore[icore]>0?log10(ntotcore[icore]):0);
                  ((TH2F*)hlist->FindObject(Form("Hpsdir_s%d",ii)))->Fill(TMath::Cos(padd->thetap),ntot>0?log10(ntot):0);
                  ((TH2F*)hlist->FindObject(Form("Hpstp_s%d",ii)))->Fill(AA,ntot>0?log10(ntot):0);
               }
               if(ii==1){
                  int ntot=0;
                  for(int jj=0;jj<npart;jj++){
                     int ipart=(padd->ipart).at(jj);
                     if(ipart!=2&&ipart!=3) continue;
                     double ppp[3]={(padd->ppp[0]).at(jj),(padd->ppp[1]).at(jj),(padd->ppp[2]).at(jj)};
                     double norm=sqrt(ppp[0]*ppp[0]+ppp[1]*ppp[1]+ppp[2]*ppp[2]);
                     double pos[3]={0,0,0};
                     double InCoo[3]={(padd->px).at(jj),(padd->py).at(jj),0},InDir[3]={ppp[0]/norm,ppp[1]/norm,ppp[2]/norm};
                     double OutCoo[3],OutDir[3];
                     ConvertCoor(pos,padd->thetap,padd->phip,InCoo,InDir,OutCoo,OutDir);
                     double OutTheta,OutPhi;
                     ConvertThetaPhi(OutDir,OutTheta,OutPhi);
                     //if(pow(TMath::Cos(padd->thetap),2)>0.9){
                     ((TH2F*)hlist->FindObject(Form("Hspos_s%d",ii)))->Fill((padd->px).at(jj),(padd->py).at(jj),(padd->pweight).at(jj));
                     ((TH2F*)hlist->FindObject(Form("Hsdir_s%d",ii)))->Fill(cos(OutTheta),OutPhi,(padd->pweight).at(jj));
                     //}
                     ntot+=1.*(padd->pweight).at(jj);
                  }
                  ((TH2F*)hlist->FindObject(Form("Hpse_s%d",ii)))->Fill(padd->ep,ntot>0?log10(ntot):0);
                  ((TH2F*)hlist->FindObject(Form("Hpsdir_s%d",ii)))->Fill(pow(TMath::Cos(padd->thetap),2),ntot>0?log10(ntot):0);
                  ((TH2F*)hlist->FindObject(Form("Hpstp_s%d",ii)))->Fill(AA,ntot>0?log10(ntot):0);
               }
            }
         }
         if(ibr==1&&pevt[ibr]){
            CorsikaEvent* padd0=(CorsikaEvent*)pevt[0];
            WFCTAEvent*   padd1=(WFCTAEvent*)pevt[ibr];
            WFCTAMCEvent* padd=(WFCTAMCEvent*)(&(padd1->mcevent));
            int size=padd->RayTrace.size();
            for(int ii=0;ii<size;ii++){
               double cweight=(padd0->cweight).at(ii);
               int raytrace=padd->RayTrace.at(ii);
               TH1F* Hy=((TH1F*)hlist->FindObject(Form("Hsn")));
               if(!Hy) continue;
               int fillcont=0;
               int nbin=Hy->GetNbinsX();
               if(raytrace>0) fillcont=nbin;
               else if(raytrace<0&&raytrace>-(nbin-0.5)) fillcont=-raytrace;

               for(int i0=0;i0<fillcont;i0++) ((TH1F*)hlist->FindObject(Form("Hsn")))->Fill(i0*1.,cweight);
               double eng=padd0->ep;
               for(int i0=0;i0<fillcont;i0++) ((TH2F*)hlist->FindObject(Form("Hsne")))->Fill(eng,i0*1.,cweight);
               int whichcore=padd0->WhichCore(padd0->cx.at(ii),padd0->cy.at(ii));
               double posx=padd0->cx.at(ii)-padd0->corex[whichcore];
               for(int i0=0;i0<fillcont;i0++) ((TH2F*)hlist->FindObject(Form("Hsnposx")))->Fill(posx,i0*1.,cweight);
               double posy=padd0->cy.at(ii)-padd0->corey[whichcore];
               for(int i0=0;i0<fillcont;i0++) ((TH2F*)hlist->FindObject(Form("Hsnposy")))->Fill(posy,i0*1.,cweight);
               double dirx=padd0->cu.at(ii);
               double diry=padd0->cv.at(ii);
               double dirz=-sqrt(1-dirx*dirx-diry-diry);
               double pos[3]={0,0,0};
               double InCoo[3]={(padd0->cx).at(ii),(padd0->cy).at(ii),0},InDir[3]={dirx,diry,dirz};
               double OutCoo[3],OutDir[3];
               ConvertCoor(pos,padd0->thetap,padd0->phip,InCoo,InDir,OutCoo,OutDir);
               double OutTheta,OutPhi;
               ConvertThetaPhi(OutDir,OutTheta,OutPhi);
               for(int i0=0;i0<fillcont;i0++) ((TH2F*)hlist->FindObject(Form("Hsndiru")))->Fill(cos(OutTheta),i0*1.,cweight);
               for(int i0=0;i0<fillcont;i0++) ((TH2F*)hlist->FindObject(Form("Hsndirv")))->Fill(OutPhi,i0*1.,cweight);
            }
            for(int ii=0;ii<NSIPM;ii++){
               int ix=ii/PIX;
               int iy=ii%PIX;
               for(int ict=0;ict<WFTelescopeArray::CTNumber;ict++) ((TH2F*)hlist->FindObject(Form("Hsimg_T%d",ict)))->Fill(ix,iy,padd->TubeSignal[ict][ii]);
            }
         }
      }
   }
   if(fillstyle==3||fillstyle==0){//fill the histograms
      CorsikaEvent* padd=(CorsikaEvent*)pevt[0];
      if(padd) ((TH1F*)hlist->FindObject(Form("MCHpe")))->Fill(padd->ep,Nuse);
   }
}

