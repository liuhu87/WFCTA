#ifndef __CorsikaEvent__
#define __CorsikaEvent__
#include "TSelector.h"
#include "CorsikaIO.h"
class EventNtuple;
class WFCTAEvent;

/*!
The class contains the complete event information of corsika simulation
*/

class CorsikaEvent: public TSelector {
   public:
  // $BEGIN -- Mark to start translation

  ////////////////////////////////////////////////////////////////
  // The ordering is not random. Unfortunately, due to alignment
  // the variables should be ordered from longer to smaller bit size
  //////////////////////////////////////////////////////////////////
   ///RUNH
   Int_t run;	//run number
   Int_t date;//!	//date
   Float_t version;//!	//version of corsika
   Float_t oheight;//!	//height of observation level

   ///EVTH
   Int_t event;	//event number
   Int_t partp;	//primary particle type
   Float_t ep;	//primary particle energy
   Float_t pp[3];//!	//primary particle momentum
   Float_t thetap;	//primary particle theta
   Float_t phip;	//primary particle phi
   Float_t stheight;	//primary starting height
   Float_t corex[NuseMax];
   Float_t corey[NuseMax];
   ///EVTE
   Int_t nphoton;//!	//total photon number
   Int_t nelectron;//!	//total electron number
   Int_t nhadron;//!	//total  hadron number
   Int_t nmuon;//!		//total muon number
   Int_t nparticle;//!	//total particle number

   ///RUNE
   Int_t nevent;//!	//total number of events

   ///CERENKOV Light and PARTICLE
   Int_t nclight;//!	//total number of cerenkov light

   vector<float> cx;//!	//x coord. when hit ground
   vector<float> cy;//!	//y coord. when hit ground
   vector<float> ct;//!	//arrival time, in ns
   vector<float> cu;//!	//cos(theta) of clight
   vector<float> cv;//!	//cos(phi) of clight
   vector<float> height;//!	//z coord. when generated
   vector<float> wavelength;//!	//wavelength of the clight
   vector<float> cweight;//!	//weight of the clight

   vector<Int_t> ipart;//!	//particle id
   vector<Int_t> igen;//!	//particle generation information
   vector<Int_t> ilevel;//!	//index of observation level

   vector<float> px;//!
   vector<float> py;//!
   vector<float> pt;//!
   vector<float> ppp[3];//!	//secondary particle momentum
   vector<float> pweight;//!	//weight of the particle

   WFCTAEvent* pwfc;  //!

  //$END -- Mark to stop translation
  /////////////////////////////////////////////////////

   public:
   ///init all to variables to 0
   void Init();
   ///Release
   void Clear();
   ///default constructor
   CorsikaEvent() {Init();}
   ///default deconstructor
   virtual ~CorsikaEvent() {Clear();}
   ///Reset all the variables
   void Reset();
   ///copy the event information from CorsikaIO to CorsikaEvent
   void Copy(CorsikaIO* pcorio);
   ///count how many secondary particles have been readed after filling previous event
   int CountParticle();
   ///which corex or corey
   int WhichCore(double x0,double y0);
   ///Do WFCTA Simulation
   bool DoWFCTASim(int iuse=0);
   ///Fill the event information to the EventNtuple
   void Fill();
};

#endif
