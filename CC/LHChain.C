#include "LHChain.h"
#include <fstream>
using namespace std;
int LHChain::jdebug=0;
int LHChain::AddFromFile(const char* rootfilelist,int beg,int end){
   std::ifstream fin(rootfilelist,std::ios::in);
   int nline=0;
   int nadded=0;
   char buff[500];
   fin.getline(buff,500);
   nline++;
   while(fin.good()){
      if(nline-1>=beg&&nline-1<=end){
         TChain::Add(buff);
         nadded++;
         if(jdebug>0) printf("Adding %d: %s\n",nline,buff);
      }
      else if(nline-1>end) break;
      fin.getline(buff,500);
      nline++;
   }
   return nadded;
}
void LHChain::Init(WFCTAEvent* event){
   if (_EVENT==NULL) {
     if (event)
       _EVENT = event;
     else
       _EVENT = new WFCTAEvent;
       _EVENT->Head()=_EVENT;
       this->SetBranchAddress(WFCTAEvent::BranchName(),&_EVENT);
       _EVENT->Tree()=NULL;
   }
}
WFCTAEvent* LHChain::_getevent(Int_t entry, Bool_t kLocal){
   Init();
   if (!kLocal) {//The old/standard way to call this method...
      if(entry>=GetEntries()) {
        delete _EVENT; _EVENT = NULL;
        _ENTRY = -2;
        return _EVENT;
      }
      _ENTRY = entry;
      m_tree_entry = LoadTree(_ENTRY);
   }
   else {// The "local" way to call this method.
      if(entry>=(GetTree()->GetEntries())){
        delete _EVENT; _EVENT = NULL;
        _ENTRY = -2;
        return _EVENT;
      }
      _ENTRY = 0;//The "correct" setting of _EVENT is not supported in this mode...
      m_tree_entry = entry;
   }

   if (GetTreeNumber()!=_TREENUMBER) {
     _TREENUMBER = GetTreeNumber();
     _EVENT->InitTree(GetTree());
     _EVENT->GetBranch(_EVENT->Tree());
   }

   //Read The Event
   if (_EVENT->GetAllContents(m_tree_entry)==false) {
     delete _EVENT; _EVENT = NULL;
     _ENTRY = -2;
   }
   //if(GetTree()){
   //   if(GetTree()->GetEntry(m_tree_entry)==false){
   //      delete _EVENT; _EVENT = NULL;
   //      _ENTRY = -2;
   //   }
   //}
   //if (_EVENT->ReadHeader(m_tree_entry)==false) {
   //  delete _EVENT; _EVENT = NULL;
   //  _ENTRY = -2;
   //}

   return _EVENT;
}
WFCTAEvent* LHChain::GetEventLocal(Int_t entry){
   return _getevent(entry, true);
}
WFCTAEvent* LHChain::GetEventGlobal(Int_t entry){
   return  _getevent(entry, false);
}
WFCTAEvent* LHChain::GetEvent(Int_t entry){
   return  _getevent(entry, false);
}

long long LHChain::GetEntryNo(UInt_t run, Int_t ev, bool runinfilename, unsigned long long maxent){
   return -1;
}

WFCTAEvent* LHChain::GetEventFast(UInt_t run, Int_t ev, bool runinfilename, unsigned long long maxent){
 long long entry= GetEntryNo(run,ev ,runinfilename,maxent);
 unsigned int frun=0;//_EVENT?_EVENT->Run():0;
 if(entry>=0){
        GetEvent(entry);
        //if(_EVENT && _EVENT->Run()==run && _EVENT->Event()==static_cast<unsigned long long>(ev)){
        //     if(frun!=run)_EVENT->UpdateSetup(run);
        //     return _EVENT;
        //}
        if(_EVENT && _EVENT->iEvent==static_cast<unsigned long long>(ev)) return _EVENT;
    }

   return NULL;

}
WFCTAEvent* LHChain::GetEvent(UInt_t run, Int_t ev, Bool_t kDontRewind){
  if (!kDontRewind) Rewind();//Go to start of chain
  // Get events in turn
  while  (GetEvent() &&
          !(_EVENT->iEvent==static_cast<unsigned long long>(ev)) ){
    if (_ENTRY==(GetEntries()-1)) {
      static bool kRewindAlreadyDid=false;
      if (!kRewindAlreadyDid) {
        printf("LHChain::GetEvent(%u, %d, %s): End of chain reached! Rewinding...\n", run, ev, (kDontRewind)?"true":"false");
        Rewind();
        kRewindAlreadyDid=true;
      }
      else {
        printf("LHChain::GetEvent(%u, %d, %s): Event not found, returning 0...\n", run, ev, (kDontRewind)?"true":"false");
        kRewindAlreadyDid=false;
        return 0;
      }
    }
  }
  return _EVENT;
};
WFCTAEvent* LHChain::GetEvent(){
  if(_ENTRY<-1)return NULL;
  GetEvent(_ENTRY+1);
  return _EVENT;
}

void LHChain::Reset(Option_t* option) {
  TChain::Reset(option);
  m_chain_Entries = 0;
  m_tree_entry = -1;
  //m_chain_Runs.clear();
  //m_chain_entryindex.clear();
  _ENTRY = -1;
  _TREENUMBER = -1;
}

/*Long64_t AMSChain::Process(TSelector*pev,Option_t*option, Long64_t nentri, Long64_t firstentry){
  nentries=nentri;
   cout<<" nentires requested "<<nentries<<endl;
  int nthreads=fThreads;
  long long nentr=0;

  //    #pragma omp parallel  default(none), shared(std::cout,option,nentries,firstentry,nentr,pev)
  ntree=fNtrees;
  if(nentries<0){
    ntree=-nentries;
    if(ntree>fNtrees)ntree=fNtrees;
    nentries=10000000000LL;
  }
  typedef multimap<uinteger,TChainElement*> fmap_d;
  typedef multimap<uinteger,TChainElement*>::iterator fmapi;
  fmap_d fmap;
  for(int i=0;i<fNtrees;i++){
    TString t1("/");
    TString t2(".");
    TString name(((TNamed*)fFiles->At(i))->GetTitle());
    TObjArray *arr=name.Tokenize(t1);
    TObjString *s=(TObjString*)arr->At(arr->GetEntries()-1);
    TObjArray *ar1=s->GetString().Tokenize(t2);
    unsigned int k=atoi(((TObjString* )ar1->At(0))->GetString().Data());
    bool bad=false;
    TChainElement *el=(TChainElement*) fFiles->At(i);
      delete arr;
    delete ar1;
  }
  fmapi it=fmap.begin();
  if(static_cast<unsigned int>(ntree)>fmap.size())ntree=fmap.size();
  cout <<"  LHChain::Process-I-Files to be processed "<<ntree<<" out of "<<fNtrees<<endl;
  if(nthreads>ntree)nthreads=ntree;
  int*ia= new int[nthreads];
  for(int i=0;i<nthreads;i++)ia[i]=0;
#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#endif
#pragma omp parallel 
  {
    int thr=0;
#ifdef _OPENMP
    thr=omp_get_thread_num();
#endif
#pragma omp  for schedule (dynamic)  nowait
    for(int i=0;i<ntree;i++){
      if(i>=ntree)continue;
      if(nentr>nentries || it==fmap.end()){
        continue;
      }
      TFile* file;
      TTree *tree;
      TSelector *curp=(TSelector*)((char*)pev+thr*fSize);
        TDirectory *gdir=0;
      fmapi itt=it;
        it++;
        for(int subdir=0;subdir<32;subdir++){

#pragma omp critical 
      {
      if(subdir==0){
        //  if(ts[thr]==0){
        //TStreamerInfo::fgInfoFactory=ts[thr]=new TStreamerInfo();
        // }
        //cout <<"thr "<<thr<<endl;
        // element=(TChainElement*) fFiles->At(i);
#ifdef CASTORSTATIC
        TRegexp d("^root:",false);
        TRegexp e("^rfio:",false);
       TString name(itt->second->GetTitle());
        //cout << " thr "<<thr<<" "<<name<<endl;          
        if(name.Contains(d))file=new TXNetFile(itt->second->GetTitle(),"READ");
        else if(name.Contains(e))file=new TRFIOFile(itt->second->GetTitle(),"READ");
        else file=new TFile(itt->second->GetTitle(),"READ");
#else
        file= TFile::Open(itt->second->GetTitle(),"READ");
#endif
          gdir=gDirectory;
//          cout <<"gdirectory!!!! "<<gdir->GetName()<<" " <<endl;
}
          tree=0;
        bool subdirexist=false;
//        cout <<"  SUBDIR****** "<<subdir<<endl;
//          cout <<"gdirectory!!!! "<<gdir<<endl;
//          cout <<"gdirectory!!!!!! "<<gdir->GetName()<<" " <<endl;
        if(subdir>0){
         char dir[1024];
          sprintf(dir,"_%d",subdir);
          TDirectory*ok=gdir->GetDirectory(dir,false,"cd");
          if(ok){
             ok->cd();
             subdirexist=true;
          }
//          cout <<"gdirectory "<<gdir->GetName()<<" "<<subdirexist <<endl;
        }
        else subdirexist=true;
        if(file && subdirexist)tree=(TTree*)file->Get(itt->second->GetName());
        if(!tree){
          if(subdirexist)cerr<<"  AMSChain::Process-E-NoTreeFound file tree "<<itt->second->GetTitle()<<" "<<itt->second->GetName()<<endl;
        }
        else{
          curp->SetOption(option);
          curp->Init(tree);
          curp->Notify();
          cout <<"  "<<i<<" "<<itt->second<<" "<<AMSEventR::_Tree->GetEntries()<<" "<<nentr<<" "<<nentries<<endl;
          //cout <<"  "<<i<<" "<<element->GetTitle()<<" "<<AMSEventR::_Tree->GetEntries()<<" "<<nentr<<" "<<nentries<<endl;
        }
      }
      if(tree){
        curp->Begin(tree);
        for(int n=0;n<AMSEventR::_Tree->GetEntries();n++){
           if(nentr++>nentries)break;
          try{
//            while(AMSEventR::_Lock&0x100000000){
//            }
//#pragma omp atomic 
//           AMSEventR::_Lock+=(1<<thr);
            curp->Process(n);
//#pragma omp atomic 
//           AMSEventR::_Lock-=(1<<thr);
          }
          catch (...){
#pragma omp critical(rd)
            if(AMSEventR::pService){
              (*AMSEventR::pService).BadEv++;
              (*AMSEventR::pService).TotalEv++;
              (*AMSEventR::pService).TotalTrig++;
            }
          }
        }
      }

#pragma omp critical (cls)  
      {
//      if(tree && AMSEventR::_Tree)nentr+=AMSEventR::_Tree->GetEntries();
        //        cout <<" finished "<<i<<" "<<endl;
      }
    }
        if(file){
           //file->Print();
           file->Close("R");
        }
        if(file)delete file;

}
#ifdef _OPENMP        
    //  this clause is because intel throutput mode deoesn;t work
    //   so simulating it
    ia[thr]=1;
      int nt=0;
      for(int j=0;j<nthreads;j++){
        if(!ia[j])nt++;
      }
    cout <<" Thread "<<thr<<" Deactivated "<<nt <<" Threads Left"<<endl;
  for(;;){
      bool work=false;
      for(int j=0;j<nthreads;j++){
        if(!ia[j]){
          work=true;
          break;
        }
      }
      if(work)usleep(200); // FIXME: This only works with ICC, and should not be necessary. kmp_get_blocktime());  
      else break;
    }
#endif
  }
  AMSEventR::_NFiles=1;
  pev->Terminate();
  delete[] ia;
  return nentr;
#endif
}*/

