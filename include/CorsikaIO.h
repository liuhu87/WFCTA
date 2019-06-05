#ifndef __CorsikaIO__
#define __CorsikaIO__
#include <fstream>
#include <iostream>
#include "TObject.h"
#include "TSelector.h"
#include "math.h"
using namespace std;

class CorsikaEvent;

//The length of a particle record (in unit of "word" = 4 bytes)
#ifdef _THIN_
#define _PARTICLE_LENGTH_ 8
#else
#define _PARTICLE_LENGTH_ 7
#endif

//Number of particles per subblock
#define _NPARTICLE_ 39

//Number of subblocks per record
#define _NSUBBLOCK_ 21

//The length of a subblock
#define _SUBBLOCK_LENGTH_   _PARTICLE_LENGTH_*_NPARTICLE_

//The length of a record
#define _RECORD_LENGTH_     _SUBBLOCK_LENGTH_*_NSUBBLOCK_

//use number of primary particle
#define Nuse 20

#define PI 3.14159265358979312
#define RADDEG 180/PI

///the event general information in the corsika output file
struct EveInfo{
  //RUNH
  Int_t irun;
  Int_t idate;
  Float_t version;
  Float_t height;
  
  //EVTH
  Int_t ievent;
  Int_t ipartp;
  Float_t ep;
  Float_t pxp, pyp, pzp;
  Float_t thetap, phip;
  Float_t corex[Nuse],corey[Nuse];
  
  //EVTE
  Int_t nphoton, nelectron, nhadron, nmuon, nparticle;
  
  //RUNE
  Int_t nevent;
};
///the cerenkov light information in the corsika output file
struct CerInfo{
  Int_t nclight;
  Float_t x, y, u, v, t, height, weight,wavelength;
};
///the particle information in the corsika output file
struct ParInfo{
  Int_t ipart;
  Int_t igen;
  Int_t ilevel;
  Float_t px, py, pz, x, y, t, weight;
};
union DataBuff{
   float f[_RECORD_LENGTH_];
   char  c[4*_RECORD_LENGTH_];
};

/*!
The class used for reading corsika output file
*/
class CorsikaIO{
   private:
   ///the buffer to save one record data
   DataBuff recbuff;

   public:
   ///read specific run
   static Int_t ReadRun;
   ///read specific event
   static Int_t ReadEvent;
   ///quick process (only check "RUNH","EVTH","EVTE" and "RUNE" blocks)
   static Bool_t RunFast;
   ///thin flag
   static Bool_t ThinFlag;
   ///wavelength flag
   static Bool_t WLFlag;
   ///control the debug output
   static Int_t jdebug;

   ///the input file type (ftype=0:cerenkov output,ftype=1:particle output,ftype=2:all output)
   Int_t ftype;
   ///the input file pointer
   ifstream* fin;
   ///number of records already readed
   Int_t nrec;
   ///number of RUNH,RUNE,EVTH,EVTE blocks already readed
   Int_t flag[4];

   ///Event General Information
   EveInfo Evt;
   ///Cerenkov light Information
   CerInfo Cer;
   ///Particle Information
   ParInfo Par;

   public:
   ///construction function
   CorsikaIO();
   ///construction function, type is the type of the input file(cerenkov light or particle),which is used to fill ftype
   CorsikaIO(char* inputfile,int type=0);
   ///deconstruction function
   virtual ~CorsikaIO();

   ///init the data(type=0:init all data,type=1:init only cer and par data,type=2:init only cer data,type=3:init only par data)
   void Init(int type=0);
   ///read all the record data (beg is the index of beginning record to read,end is the index of ending record to read. equal to 0 means read all the records)
   Int_t ReadAll(int beg=0,int end=0,CorsikaEvent* pevt=0);
   ///read one record data
   Int_t ReadRec(CorsikaEvent* pevt=0);
   ///read one block data (iblk is the index of block to read,iblk<=0 means read all the blocks)
   Int_t FillBlk(int iblk_beg=0,int iblk_end=0,CorsikaEvent* pevt=0);
   ///read one particle or one cerenkov light data (iblk is the index of block to read, ipar is the index of particle to read,ipar<=0 means read all the particles belong to this block)
   Int_t FillPar(int iblk, int ipar=0,CorsikaEvent* pevt=0);
   ///Check Run and Event Number
   bool CheckRunEvent(CorsikaEvent* pevt=0);
};

#endif
