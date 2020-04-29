#include "LHChain.h"
#include "TFile.h"
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
       _EVENT->SetLHChain(this);
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
        if(_EVENT && _EVENT->rabbitTime==run && _EVENT->iEvent==static_cast<unsigned long long>(ev)){
             //if(frun!=run)_EVENT->UpdateSetup(run);
             return _EVENT;
        }
        //if(_EVENT && _EVENT->iEvent==static_cast<unsigned long long>(ev)) return _EVENT;
    }

   return NULL;

}
WFCTAEvent* LHChain::GetEvent(UInt_t run, Int_t ev, Bool_t kDontRewind){
  if (!kDontRewind) Rewind();//Go to start of chain
  // Get events in turn
  while  (GetEvent() &&
          !(_EVENT->rabbitTime==run && _EVENT->iEvent==static_cast<unsigned long long>(ev)) ){
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

const char* LHChain::GetFileName(){
  TFile* curf=((TChain*)this)->GetFile();
  return curf?((const char*)curf->GetName()):0;
}

int LHChain::OpenOutputFile(const char* filename){
  if(fout&&wfctanew) return 0;
  TDirectory* gdir=gDirectory;
  fout= TFile::Open(filename,"RECREATE");
  if(!fout){
    cout <<" LHChain::OpenOutpuFile-E- Cannot open "<<filename<<" for output"<<endl;
    return -1;
  }
  fout->cd();
  //TChain::InvalidateCurrentTree();
  wfctanew = CloneTree(0);
  gdir->cd();
  return 0;
}
void LHChain::SaveCurrentEvent(){
  if(!fout) OpenOutputFile("root:://eos01.ihep.ac.cn//eos/user/h/hliu/SelectedEvents.root");
  if(!fout){
    cout <<" LHChain::SaveCurrentEntry-E- Cannot open file  for output no events are saved"<<endl;
    return;
  }
  TDirectory* gdir=gDirectory;
  if(_EVENT){
    _EVENT->GetAllContents();
    fout->cd();
    wfctanew->Fill();
  }
  gdir->cd();
  return;
}
void LHChain::CloseOutputFile(){
  if(!fout) return;
  cout << "LHChain::CloseOutputFile WFCTA ROOT file \"";
  cout << fout->GetName() << "\" with " << wfctanew->GetEntries();
  cout << " selected events" << endl;
  if(wfctanew){
     fout->cd();
     //wfctanew->SetDirectory(fout);
     wfctanew->Write();
  }
  fout->Close();
  fout=0;
  return;
}

