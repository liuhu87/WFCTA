#ifndef __LHAASOChain__
#define __LHAASOChain__

#include "WFCTAEvent.h"
#include "TChain.h"

const unsigned int _UINT_MAX=4294967295;
class LHChain : public TChain {
   public:
   static int jdebug;
   private:
   int m_chain_Entries;
   int m_tree_entry;
   unsigned int fThreads;
   unsigned int fSize;
   WFCTAEvent* _EVENT;
   Int_t _ENTRY;
   const char* _NAME;
   Int_t _TREENUMBER;

   ///Get WFCTAEvent in entry number "entry", if kLocal is false the entry number is "global", if true the entry number is local w.r.t. curret tree
   WFCTAEvent* _getevent(Int_t entry, Bool_t kLocal=false);

   public:
   Long64_t ntree;
   Long64_t nentries;
   int get_tree_entry()const {return m_tree_entry;}
   /// Default constructor (it builds automatically the WFCTAEvent object)
   LHChain(const char* name="eventShow", unsigned int thr=1,unsigned int size=sizeof(WFCTAEvent))
    :TChain(name),m_chain_Entries(0),m_tree_entry(-1),fThreads(thr),fSize(size),_EVENT(NULL),_ENTRY(-1),_NAME(name),_TREENUMBER(-1){}
   /// Destructor
   virtual ~LHChain(){ if (_EVENT) delete _EVENT; }
   using TChain::Add;
   int AddFromFile(const char* rootfilelist,int beg=0,int end=_UINT_MAX);
   ///Set event branch and links; called after reading of all trees; called automatically in GetEvent
   void Init(WFCTAEvent* event=0);
   ///Get next Event object in the chain
   WFCTAEvent* GetEvent();

   ///Get Event in entry number "entry" (the entry number is "global", i.e. along the chain)
   WFCTAEvent* GetEvent(Int_t entry);
   WFCTAEvent* GetEventGlobal(Int_t entry);

   ///Get Event in entry number "entry" (the entry number is "local", i.e. along the curret tree)
   WFCTAEvent* GetEventLocal(Int_t entry);

   long long GetEntryNo(UInt_t run, Int_t ev, bool runinfilename=false, unsigned long long maxent=100000000);

   ///Get Event with run number "run" and event number "ev"
   ///If kDontRewind=true the event will be searched starting from current event. If kDontRewind=false (default) the whole chain is rewinded and event will be searched starting from the first event in the chain
   WFCTAEvent* GetEventFast(UInt_t run, Int_t ev, bool runinfilename, unsigned long long maxent);
   WFCTAEvent* GetEvent(UInt_t run, Int_t ev, Bool_t kDontRewind=false);
   /// Reset this chain
   void Reset(Option_t* option = "");

   ///Rewind the chain (go back before first entry)
   void Rewind() {_ENTRY=-1;_TREENUMBER=-1;};

   ///Get the current entry number
   Int_t Entry(){return _ENTRY;}

   ///Get the current event pointer
   WFCTAEvent* pEvent() {return _EVENT;}

   ///Set the current event pointer
   void SetEvent(WFCTAEvent *ev) {_EVENT=ev;}

   ///Get the name of the tree
   const char* ChainName(){return _NAME;}
   //Long64_t  Process(TSelector*pev, Option_t *option="", Long64_t nentries=1000000000000LL, Long64_t firstentry=0);

   ClassDef(LHChain,5);
};
#endif
