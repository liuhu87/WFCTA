#include "CorsikaChain.h"
#include "EventNtuple.h"
void CorsikaChain::Init(){
   CorsikaIO::Init();
   nfiles=0;
   filename=0;
   ifile=0;
}
void CorsikaChain::Add(char* inputfile){
   if(!inputfile) return;
   nfiles++;
   if(nfiles%10==1){
      int nlength=(nfiles-1);
      char** namebuff=new char *[nlength+10];
      for(int ii=0;ii<nlength+10;ii++){
         namebuff[ii]=new char[charsize];
         if(ii<nlength) strcpy(namebuff[ii],filename[ii]);
         else if(ii==nlength) strcpy(namebuff[ii],inputfile);
         else namebuff[ii][0]='\0';
         //printf("Adding:%d %s\n",ii,namebuff[ii]);
      }
      if(filename){
         for(int ii=0;ii<nlength;ii++) delete []filename[ii];
         delete []filename;
      }
      filename=namebuff;
   }
   else{
      //printf("Adding2:%s\n",inputfile);
      strcpy(filename[nfiles-1],inputfile);
   }
   return;
}
CorsikaChain::CorsikaChain(char* filelist,int first,int last,int type){
   ifstream ff;
   bool failed=(last<first);
   if(!filelist) failed=true;
   else{
      ff.open(filelist);
      if(!ff.is_open()) failed=true;
      else ff.seekg(0,ios::beg);
   }
   Init();
   if(failed) return;
   char buff[charsize];
   int nline=-1;
   ff.getline(buff,charsize);
   nline++;
   while(ff.good()){
      bool isinc;
      if(last<0) isinc=true;
      else isinc=(nline>=first&&nline<=last);
      if(isinc) Add(buff);
      ff.getline(buff,charsize);
      nline++;
   }
   //init the ifstream for CorsikaIO class
   ftype=-1;
   ifile=0;
   CorsikaIO::Init(0);
   nrec=0;
   for(int ii=0;ii<4;ii++) flag[ii]=0;
   if(type<0||type>2){
     printf("CorsikaChain::CorsikaChain: the inputfile type is wrong type=%d\n",type);
     return;
   }
   fin=new std::ifstream(filename[ifile]);
   if(!fin->is_open()){
      printf("CorsikaChain::CorsikaChain:Error in opening file %s\n",filename[ifile]);
      delete fin; fin=0;
      return;
   }
   else{
      fin->seekg(0,ios::end);
      size_t filesize = fin->tellg();
      if (jdebug>0) printf("The file size of %s is %d bytes\n",filename[ifile],filesize);
      if (!filesize) {
        printf("Warning: The file is empty!\n");
        fin->close();
        delete fin; fin=0;
        return;
      }
      fin->seekg(0,ios::beg);
      ftype=type;
      return;
   }
}
void CorsikaChain::Clear(){
   if(filename){
      for(int ii=0;ii<sizeof(filename);ii++) delete []filename[ii];
      delete []filename;
   }
}
Int_t CorsikaChain::ReadRec(CorsikaEvent* pevt){
   int nrec0=nrec;
   int irec=CorsikaIO::ReadRec(pevt);
   if(irec<=0){
      printf("CorsikaChain::ReadRec: there is error during reading File:%s No. %d record\n",(filename&&ifile<nfiles)?filename[ifile]:"",nrec>=0?(nrec+1):1);
      return 0;
   }
   else if(irec==nrec0){//at the end of file
      ifile++;
      if(ifile>=nfiles){
         printf("CorsikaChain::ReadRec: Already at the end of chain, reading all successfully.\n");
         ifile=nfiles-1;
         return nrec;
      }
      else{
         if(fin) {delete fin; fin=0;}
         if(!filename[ifile]){
            printf("CorsikaChain::ReadRec: filename:%s error,exiting\n",filename[ifile]);
            return 0;
         }
         else{
            fin=new std::ifstream(filename[ifile]);
            if(!fin->is_open()){
               printf("CorsikaChain::ReadRec:Error in opening file %s\n",filename[ifile]);
               delete fin; fin=0;
               return 0;
            }
            else{
               fin->seekg(0,ios::end);
               size_t filesize = fin->tellg();
               if (jdebug>0) printf("The file size of No.%d:%s is %d bytes\n",ifile,filename[ifile],filesize);
               if (!filesize) {
                 printf("Warning: The file is empty!\n");
                 fin->close();
                 delete fin; fin=0;
                 return 0;
               }
               fin->seekg(0,ios::beg);
               //try to read one record
               printf("CorsikaChain::ReadRec: process one new file No.%10d %s\n",ifile,filename[ifile]);
               return CorsikaIO::ReadRec(pevt);
            }
         }
      }
   }
}
Int_t CorsikaChain::ReadAllRec(int beg,int end,CorsikaEvent* pevt){
   if(nfiles<=0){
      printf("CorsikaChain::ReadAllRec: there is no file, exiting\n");
      return 0;
   }
   int nrec0=nrec;
   int count=0;
   while(true){
      int nrec00=nrec;
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
         printf("CorsikaChain::ReadAll: Error in Reading No. %d Record.\n",nrec>=0?(nrec+1):1);
         break;
      }
      else if(irec==nrec00){ //all the records readed
         break;
      }
      count++;
      if(count>1.0e20){
         printf("CorsikaChain::ReadAll: Too many records(%d),exiting.\n",nrec>=0?(nrec+1):1);
         break;
      }
   }
   if(pevt){
      //check wheather there is event after filling the tree
      if(CheckRunEvent(pevt)){
         if(jdebug>0) printf("CorsikaChain::ReadAllRec: Fill the remaining data(run=%d evt=%d) on record %d(%d,%d,%d,%d)\n",pevt->run,pevt->event,nrec,flag[0],flag[1],flag[2],flag[3]);
         if(ReadRun>=0||ReadEvent>=0) printf("EventNtuple::Fill run=%10d event=%7d\n",pevt->run,pevt->event);
         pevt->Fill();
      }
   }
   return nrec-nrec0;
}
Int_t CorsikaChain::ReadAllFile(int beg,int end,CorsikaEvent* pevt){
   if(nfiles<=0){
      printf("CorsikaChain::ReadAllFile: there is no file, exiting\n");
      return 0;
   }
   int nrec0=nrec;
   int count=0;
   while(true){
      int nrec00=nrec;
      int irec=0;
      if(beg<0&&end<0){
         irec=ReadRec(pevt); //read all the records
      }
      else{
         if(ifile<beg) irec=ReadRec(0);
         else if(ifile>end) break;
         else irec=ReadRec(pevt);
      }
      if(irec<=0){
         printf("CorsikaChain::ReadAll: Error in Reading File%d:%s No. %d Record.\n",ifile,filename[ifile],nrec>=0?(nrec+1):1);
         break;
      }
      else if(irec==nrec00){ //all the records readed
         break;
      }
      count++;
      if(count>1.0e20){
         printf("CorsikaChain::ReadAll: Too many records(%d),exiting.\n",nrec>=0?(nrec+1):1);
         break;
      }
   }
   if(pevt){
      //check wheather there is event after filling the tree
      if(CheckRunEvent(pevt)){
         if(jdebug>0) printf("CorsikaChain::ReadAllFile: Fill the remaining data(run=%d evt=%d) on record %d(%d,%d,%d,%d)\n",pevt->run,pevt->event,nrec,flag[0],flag[1],flag[2],flag[3]);
         if(ReadRun>=0||ReadEvent>=0) printf("EventNtuple::Fill run=%10d event=%7d\n",pevt->run,pevt->event);
         pevt->Fill();
      }
   }
   return nrec-nrec0;
}
