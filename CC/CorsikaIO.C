#include "CorsikaIO.h"
#include "CorsikaEvent.h"
Int_t CorsikaIO::ReadRun=-1;
Int_t CorsikaIO::ReadEvent=-1;
Bool_t CorsikaIO::RunFast=false;
#ifdef _THIN_
Bool_t CorsikaIO::ThinFlag=true;
#else
Bool_t CorsikaIO::ThinFlag=false;
#endif
Bool_t CorsikaIO::WLFlag=true;
Int_t CorsikaIO::jdebug=0;
CorsikaIO::CorsikaIO(){
   Init(0);
   nrec=0;
   for(int ii=0;ii<4;ii++) flag[ii]=0;
   fin=0;
   ftype=-1;
}
bool CorsikaIO::OpenFile(const char* inputfile){
   if(fin){fin->close(); delete fin;}
   if(inputfile) fin=new std::ifstream(inputfile);
   else fin=0;
   if(!fin) return false;
   if(!fin->is_open()){
      printf("CorsikaIO::CorsikaIO:Error in opening file %s\n",inputfile);
      delete fin; fin=0;
      return false;
   }
   else{
      fin->seekg(0,ios::end);
      size_t filesize = fin->tellg();
      if (jdebug>0) printf("The file size of %s is %d bytes\n",inputfile,filesize);
      if (!filesize) {
        printf("Warning: The file is empty!\n");
        fin->close();
        delete fin; fin=0;
        return false;
      }
      fin->seekg(0,ios::beg);
      return true;
   }
}
CorsikaIO::CorsikaIO(const char* inputfile,int type){
   Init(0);
   nrec=0;
   for(int ii=0;ii<4;ii++) flag[ii]=0;
   fin=0;
   ftype=-1;
   if(type<0||type>2){
     printf("CorsikaIO::CorsikaIO: the inputfile type is wrong type=%d\n",type);
     return;
   }
   bool opened=OpenFile(inputfile);
   if(opened) ftype=type;
}

CorsikaIO::~CorsikaIO(){
   if(ftype>=0){
      if(fin){
         fin->close();
         delete fin;
      }
   }
}
void CorsikaIO::Init(int type){
   if(type==1){
      memset(&(Cer),0,sizeof(Cer));
      memset(&(Par),0,sizeof(Par));
   }
   else if(type==2){
      memset(&(Cer),0,sizeof(Cer));
   }
   else if(type==3){
      memset(&(Par),0,sizeof(Par));
   }
   else{
      memset(&(recbuff),0,sizeof(recbuff));
      memset(&(Evt),0,sizeof(Evt));
      memset(&(Cer),0,sizeof(Cer));
      memset(&(Par),0,sizeof(Par));
      Evt.irun=-1;
      Evt.ievent=-1;
   }
}
void CorsikaIO::Reset(){
   nrec=0;
   for(int ii=0;ii<4;ii++) flag[ii]=0;
   Evt.irun=-1;
   Evt.ievent=-1;
   if(fin){
      fin->clear();
      fin->seekg(0,ios::beg);
   }
}
Int_t CorsikaIO::ReadAll(int beg,int end,CorsikaEvent* pevt){
   if(ftype<0) return 0;
   int nrec0=nrec;
   while(fin->good()){
      int irec=0;
      if(beg<1&&end<1){
         irec=ReadRec(pevt); //read all the records
         //printf("CorsikaIO::ReadAll: Read No.%d record\n",nrec);
      }
      else{
         if(nrec+1<beg) irec=ReadRec(0);
         else if(nrec+1>end) break;
         else irec=ReadRec(pevt);
         //printf("CorsikaIO::ReadAll: Read No.%d record\n",nrec);
      }
      if(irec<=0){
         printf("CorsikaIO::ReadAll: Error in Reading No. %d Record.\n",nrec>=0?(nrec+1):1);
         break;
      }
   }
   if(pevt){
      //check wheather there is event after filling the tree
      if(CheckRunEvent(pevt)){
         if(jdebug>0) printf("CorsikaIO::ReadAll: Fill the remaining data(run=%d evt=%d) on record %d(%d,%d,%d,%d)\n",pevt->run,pevt->event,nrec,flag[0],flag[1],flag[2],flag[3]);
         if(ReadRun>=0||ReadEvent>=0) printf("EventNtuple::Fill run=%10d event=%7d\n",pevt->run,pevt->event);
         pevt->Fill();
      }
   }
   return nrec-nrec0;
}
Int_t CorsikaIO::ReadRec(CorsikaEvent* pevt){
   if(ftype<0||(!fin)) return 0;
   union {
      int i;
      char c[4];
   } padding;
   //There are padding words before and after a record in the fortran output
   fin->read(padding.c,4);
   if (!fin->good()) {
     if (flag[0]<1||flag[1]<1||nrec==0) {
       printf("CorsikaIO::ReadRec: Error in reading corsika file\n");
       printf("CorsikaIO::ReadRec: Is the file truncated?\n");
       if (nrec==0)  printf("CorsikaIO::ReadRec: Is the file empty?\n");
       if (flag[0]!=1) printf("CorsikaIO::ReadRec: There is no RUNH sub-block!\n");
       if (flag[1]!=1) printf("CorsikaIO::ReadRec: There is no RUNE sub-block!\n");
     }
     //A successful reading would end here
     printf("CorsikaIO::ReadRec: All Reading finished. Succeed in reading corsika file (total %d records)\n",nrec);
     return nrec;
   }
   if (jdebug>0||(nrec%10000==1)) printf("CorsikaIO::ReadRec: Begin of Reading No. %d Records. begin padding=%d %s\n",nrec>=0?(nrec+1):1,padding.i,padding.c);

   int iflag = 0;
   if (padding.i!=4*_RECORD_LENGTH_){
      iflag = 1;
      printf("CorsikaIO::ReadRec: size miss match between _RECORD_LENGTH_ and 1st padding %d %d\n",padding.i,4*_RECORD_LENGTH_);
   }
   fin->read(recbuff.c,padding.i);
   if (!fin->good()) {
      iflag = 1;
      printf("CorsikaIO::ReadRec: Error in reading the block content\n");
   }
   fin->read(padding.c,4);
   if (padding.i!=4*_RECORD_LENGTH_) {
      iflag = 1;
      printf("CorsikaIO::ReadRec: size miss match between _RECORD_LENGTH_ and 2nd padding %d %d\n",padding.i,4*_RECORD_LENGTH_);
   }
   if (iflag) {
     printf("********CorsikaIO::ReadRec: Error in reading corsika file\n");
     printf("********CorsikaIO::ReadRec: Is the file truncated?\n");
   }
   if (jdebug>0||(nrec%10000==1)) printf("CorsikaIO::ReadRec: End   of Reading No. %d Records. end   padding=%d %s\n",nrec>=0?(nrec+1):1,padding.i,padding.c);

   if(iflag) return 0;
   else{
      nrec++;
      FillBlk(0,0,pevt);
      return nrec;
   }
}
Int_t CorsikaIO::FillBlk(int iblk_beg,int iblk_end,CorsikaEvent* pevt){
   if(ftype<0) return 0;
   char sss[5];
   double corex, corey, Firstheight, lightspeed;
   int nblk=0;
   for (int isubblock = 1; isubblock <= _NSUBBLOCK_; isubblock++) {
      if((iblk_end>=iblk_beg&&iblk_beg>0)&&(isubblock<iblk_beg||isubblock>iblk_end)) continue;
      if (jdebug>1) printf("No. %6d sub-block\n",isubblock);
      int iptr = _SUBBLOCK_LENGTH_*(isubblock-1) - 1;

      memcpy(sss,&recbuff.f[iptr+1],4);
      sss[4] = '\0';
      if (strcmp(sss,"RUNH")==0) {
        if (jdebug>1){
           printf("sss = %s\n",sss);
           printf(">>> Begin a run\n");
        }
        flag[0]++;
        //memcpy(&(Evt.irun),&recbuff.f[iptr+2],sizeof(float));
        //memcpy(&(Evt.idate),&recbuff.f[iptr+3],sizeof(float));
        //memcpy(&(Evt.nlevel),&recbuff.f[iptr+5],sizeof(float));
        //printf("run=%d %f\n",Evt.irun,recbuff.f[iptr+2]);
        //printf("date=%d %f\n",Evt.idate,recbuff.f[iptr+3]);
        //printf("nlevel=%d %f\n",Evt.nlevel,recbuff.f[iptr+5]);
        Evt.irun = recbuff.f[iptr+2];
        Evt.idate = recbuff.f[iptr+3];
        Evt.height = recbuff.f[iptr+5+1];
        Evt.version = recbuff.f[iptr+4];
        if(jdebug>2){
           printf("%f %f %f %f %f\n",recbuff.f[iptr+2],recbuff.f[iptr+3],recbuff.f[iptr+4],recbuff.f[iptr+5],recbuff.f[iptr+5+1]);
           cout << "Run number: "                  << Evt.irun << endl;
           cout << "Date of the run: "             << Evt.idate << endl;
           cout << "Version of CORSIKA: "          << Evt.version << endl;
           cout << "Number of observe level: "     << recbuff.f[iptr+5] << endl;
           for (int i=1; i<=recbuff.f[iptr+5]; i++) {
             cout << "level " << i << ": "         << recbuff.f[iptr+5+i] << endl;
           }
           cout << "Slope of energy spectrum: "    << recbuff.f[iptr+16] << endl;
           cout << "Lower limit of energy range: " << recbuff.f[iptr+17] << endl;
           cout << "Upper limit of energy range: " << recbuff.f[iptr+18] << endl;
           cout << "Flag for EGS4: "               << recbuff.f[iptr+19] << endl;
           cout << "Flag for NGK: "                << recbuff.f[iptr+20] << endl;
           cout << "Cutoff for hadrons: "          << recbuff.f[iptr+21] << endl;
           cout << "Cutoff for muons: "            << recbuff.f[iptr+22] << endl;
           cout << "Cutoff for electrons: "        << recbuff.f[iptr+23] << endl;
           cout << "Cutoff for photons: "          << recbuff.f[iptr+24] << endl;
        }
      }
      else if (strcmp(sss,"EVTH")==0) {
        if (jdebug>1){
           printf("sss = %s\n",sss);
           printf(">>> Begin an event %d\n",(int)recbuff.f[iptr+2]);
        }
        flag[2]++;
        //memcpy(&(Evt.ievent),&recbuff.f[iptr+2],sizeof(float));
        //memcpy(&(Evt.ipartp),&recbuff.f[iptr+3],sizeof(float));
        Evt.ievent = recbuff.f[iptr+2];
        Evt.ipartp = recbuff.f[iptr+3];
        Evt.ep = recbuff.f[iptr+4];
        Firstheight = recbuff.f[iptr+7];
        Evt.pxp = recbuff.f[iptr+8];
        Evt.pyp = recbuff.f[iptr+9];
        Evt.pzp = recbuff.f[iptr+10];
        lightspeed = 29.97*Evt.pzp/sqrt(Evt.pxp*Evt.pxp+Evt.pyp*Evt.pyp+Evt.pzp*Evt.pzp);
        Evt.thetap = recbuff.f[iptr+11];
        Evt.phip = recbuff.f[iptr+12];
        corex = recbuff.f[iptr+98+1];
        corey = recbuff.f[iptr+118+1];
        for(int ii=1;ii<=Nuse;ii++){
        Evt.corex[ii-1]=recbuff.f[iptr+98+ii];
        Evt.corey[ii-1]=recbuff.f[iptr+118+ii];
        }
        //Firstheight = -(Firstheight+430000.);
        //Firstheight = sqrt(Firstheight*Firstheight+corex*corex+corey*corey);
        if (jdebug>2) {
          printf("%f %f %f %f %f %f %f\n",(Firstheight)/lightspeed,Firstheight,recbuff.f[iptr+11], lightspeed, Evt.pzp, corex, corey);
          printf("corex %f corey %f ievent %d\n",recbuff.f[iptr+98+1],recbuff.f[iptr+118+1],Evt.ievent);
          cout << "Event number: "               << Evt.ievent << endl;
          cout << "Kind of primary particle: "   << Evt.ipartp << endl;
          cout << "Energy of primary particle: " << Evt.ep << endl;
          cout << "Px of primary particle: "     << Evt.pxp << endl;
          cout << "Py of primary particle: "     << Evt.pyp << endl;
          cout << "Pz of primary particle: "     << Evt.pzp << endl;
          cout << "Theta of primary particle: "  << Evt.thetap*RADDEG << endl;
          cout << "Phi of primary particle: "    << Evt.phip*RADDEG << endl;
          // ...... many other parameters, please read the document
        }
      }
      else if (strcmp(sss,"EVTE")==0) {
        if (jdebug>1){
           printf("sss = %s\n",sss);
           printf(">>> End the event %d\n",(int)recbuff.f[iptr+2]);
        }
        flag[3]++;
        //memcpy(&(Evt.nphoton),&recbuff.f[iptr+3],sizeof(float));
        //memcpy(&(Evt.nelectron),&recbuff.f[iptr+4],sizeof(float));
        //memcpy(&(Evt.nhadron),&recbuff.f[iptr+5],sizeof(float));
        //memcpy(&(Evt.nmuon),&recbuff.f[iptr+6],sizeof(float));
        //memcpy(&(Evt.nparticle),&recbuff.f[iptr+7],sizeof(float));
        Evt.nphoton = recbuff.f[iptr+3];
        Evt.nelectron = recbuff.f[iptr+4];
        Evt.nhadron = recbuff.f[iptr+5];
        Evt.nmuon = recbuff.f[iptr+6];
        Evt.nparticle = recbuff.f[iptr+7];
        if (jdebug>2) {
          cout << "Total photons: "    << Evt.nphoton << endl;
          cout << "Total electrons: "  << Evt.nelectron << endl;
          cout << "Total hadrons: "    << Evt.nhadron << endl;
          cout << "Total muons: "      << Evt.nmuon << endl;
          cout << "Total particles: "  << Evt.nparticle << endl;
        }
        ///Fill the event after one event finished
        if(pevt){
           if(CheckRunEvent(pevt)){
              if(jdebug>0) printf("CorsikaIO::ReadAll: Fill the data on record %d(%d,%d,%d,%d)\n",nrec,flag[0],flag[1],flag[2],flag[3]);
              if(ReadRun>=0||ReadEvent>=0) printf("EventNtuple::Fill run=%10d event=%7d\n",pevt->run,pevt->event);
              pevt->Fill();
           }
        }
      }
      else if (strcmp(sss,"RUNE")==0) {
        if (jdebug>1){
           printf("sss = %s\n",sss);
           printf(">>> End the run\n");
        }
        flag[1]++;
        //memcpy(&(Evt.nevent),&recbuff.f[iptr+3],sizeof(float));
        Evt.nevent = recbuff.f[iptr+3];
        if (jdebug>2) {
          cout << "Number of events processed: " << Evt.nevent << endl;
        }
      }
      else if (strcmp(sss,"LONG")==0) {
        if (jdebug>1){
           printf("sss = %s\n",sss);
           printf(">>>>> Logitudinal sub-block\n");
        }
        if (jdebug>2) {
          cout << "Event Number: "     << int(recbuff.f[iptr+2]) << endl;
          cout << "Particle ID: "      << int(recbuff.f[iptr+3]) << endl;
          cout << "Total energy: "     << recbuff.f[iptr+4] << endl;
          cout << "Total steps: "      << int(recbuff.f[iptr+5])/100 << endl;
          cout << "Number of blocks: " << int(recbuff.f[iptr+5])%100 << endl;
          cout << "Current block: "    << int(recbuff.f[iptr+6]) << endl;
          cout << "Zenith angle: "     << recbuff.f[iptr+8]*RADDEG << endl;
          cout << "Azimuth angle: "    << recbuff.f[iptr+9]*RADDEG << endl;
          for (int istep=1; istep<=26; istep++) {
            int iptrnow = (istep-1)*10 + iptr;
            if (recbuff.f[iptrnow+14] != 0) {
              cout << "Vertical Depth: "         << recbuff.f[iptrnow+14] << endl;
              cout << "Gammas in the depth: "    << recbuff.f[iptrnow+15] << endl;
              cout << "Positrons in the depth: " << recbuff.f[iptrnow+16] << endl;
              cout << "Electrons in the depth: " << recbuff.f[iptrnow+17] << endl;
              cout << "Muon+ in the depth: "     << recbuff.f[iptrnow+18] << endl;
              cout << "Muon- in the depth: "     << recbuff.f[iptrnow+19] << endl;
              cout << "Hadronic in the depth: "  << recbuff.f[iptrnow+20] << endl;
              cout << "Charged in the depth: "   << recbuff.f[iptrnow+21] << endl;
              cout << "Nuclei in the depth: "    << recbuff.f[iptrnow+22] << endl;
              cout << "Cherenkov in the depth: " << recbuff.f[iptrnow+23] << endl;
            }
          }
        }
      }
      else{
        if (jdebug>1){
           printf("sss = %s(%d)\n",sss,isubblock);
           printf(">>>>> Particle sub-block\n");
        }
        if(CheckRunEvent()) FillPar(isubblock,0,pevt);
      } //end of particle sub-block
      nblk++;
   } //end of sub-block loop
   return nblk;
}
Int_t CorsikaIO::FillPar(int iblk,int ipar,CorsikaEvent* pevt){
   if(RunFast) return 0; //do not process secondary particles if RunFast
   if(iblk<1||iblk>_NSUBBLOCK_) return 0;
   if(ftype<0) return 0;
   int iptr = _SUBBLOCK_LENGTH_*(iblk-1) - 1;
   char sss[5];
   memcpy(sss,&recbuff.f[iptr+1],4);
   sss[4] = '\0';
   if (strcmp(sss,"RUNH")==0) return 0;
   else if(strcmp(sss,"RUNE")==0) return 0;
   else if(strcmp(sss,"EVTH")==0) return 0;
   else if(strcmp(sss,"EVTE")==0) return 0;
   else if(strcmp(sss,"LONG")==0) return 0;
   int npar=0;
   bool FillCer=true;
   bool FillPar=true;
   FillCer=(ftype==0||ftype==2);
   FillPar=(ftype==1||ftype==2);
   for (int iptcl=1; iptcl<=_NPARTICLE_; iptcl++) {
      if((ipar>=1&&ipar<=_NPARTICLE_)&&(ipar!=iptcl)) continue;
      if (jdebug>3) printf("No. %6d sub-block, No. %6d particle\n",iblk,iptcl);
      int iptrnow = (iptcl-1)*_PARTICLE_LENGTH_ + iptr;
      if(FillCer){
         #ifdef _THIN_
         //if(ThinFlag)
         Cer.weight = recbuff.f[iptrnow+8];
         #else
         //else
         Cer.weight = 1;
         //
         #endif
         if(WLFlag){
         memcpy(&(Cer.nclight),&recbuff.f[iptrnow+1],sizeof(float));
         Cer.wavelength = (Cer.nclight%1000);
         Cer.nclight = Cer.nclight/1000;
         }
         else{
         Cer.wavelength = -1;
         Cer.nclight=recbuff.f[iptrnow+1];
         }
         if (Cer.nclight) {
           Cer.x = recbuff.f[iptrnow+2];
           Cer.y = recbuff.f[iptrnow+3];
           Cer.u = recbuff.f[iptrnow+4];
           Cer.v = recbuff.f[iptrnow+5];
           Cer.t = recbuff.f[iptrnow+6];
           Cer.height = recbuff.f[iptrnow+7];
           if (jdebug>4) {
             cout << "nclight = " << Cer.nclight << "  wavelength = " << Cer.wavelength << endl;
             cout << "u, v, height = " << Cer.u << " " << Cer.v << " " << Cer.height << endl;
             cout << "x, y, t = " << Cer.x << " " << Cer.y << " " << Cer.t << endl;
             cout << "weight = " << Cer.weight << endl;
           }
           ///Push the particle to the vector
           if(pevt){
              pevt->Copy(this);
           }
         }
         else {
           Cer.weight = 0;
           if (jdebug>4) cout << "Empty event bank (all zero)" << endl;
         }
      }
      if(FillPar){
         #ifdef _THIN_
         //if(ThinFlag)
         Par.weight = recbuff.f[iptrnow+8];
         #else
         //else
         Par.weight = 1;
         //
         #endif
         //int idd; memcpy(&(idd),&recbuff.f[iptrnow+1],sizeof(float));
         int idd = recbuff.f[iptrnow+1];
         if (idd) {
           Par.ipart = idd/1000;
           Par.igen = (idd%1000)/10;
           Par.ilevel = idd%10;
           Par.px = recbuff.f[iptrnow+2];
           Par.py = recbuff.f[iptrnow+3];
           Par.pz = recbuff.f[iptrnow+4];
           Par.x = recbuff.f[iptrnow+5];
           Par.y = recbuff.f[iptrnow+6];
           Par.t = recbuff.f[iptrnow+7];
           if (jdebug>4) {
             cout << "ievent x, y, t = " << recbuff.f[iptrnow+1]<<" "<<recbuff.f[iptrnow+2] << " " << recbuff.f[iptrnow+3] << " " << recbuff.f[iptrnow+4] << endl;
             cout << "ipart, igen, ilevel = " << Par.ipart << " " << Par.igen << " " << Par.ilevel << endl;
             cout << "px, py, pz = " << Par.px << " " << Par.py << " " << Par.pz << endl;
             cout << "x, y, t = " << Par.x << " " << Par.y << " " << Par.t << endl;
             cout << "weight = " << Par.weight << endl;
           }
           ///Push the particle to the vector
           if(pevt){
              pevt->Copy(this);
           }
         }
         else {
           Par.weight = 0;
           if (jdebug>4) cout << "Empty event bank (all zero)" << endl;
         }
      }
      npar++;
   }
   return npar;
}
bool CorsikaIO::CheckRunEvent(CorsikaEvent* pevt){
   //check wheather specific run and event
   int run0=pevt?(pevt->run):Evt.irun;
   int event0=pevt?(pevt->event):Evt.ievent;
   if(run0<0||event0<0) return false;
   bool runmatch=(ReadRun<0||run0==ReadRun);
   bool eventmatch=(ReadEvent<0||event0==ReadEvent);
   return (runmatch&&eventmatch);
}
