#ifndef __CorsikaChain__
#define __Corsikachain__
#include "CorsikaIO.h"

#define charsize 500
/*!
The class used for linking all the CorsikaIO as a chain
*/
class CorsikaChain : public CorsikaIO{
   public:
   ///number of corsika output files
   int nfiles;
   ///all the file names for each file
   char** filename;
   ///index of file which is under processing
   int ifile;

   public:
   ///Init the class
   void Init();
   ///Constructor without filelist
   CorsikaChain() {Init();}
   ///Constructor with filelist (process the files from first to last, first<0&&last<0 means process all the files,first or last vary from 0 to nline-1)
   CorsikaChain(char* filelist,int first=-1,int last=-1,int type=0);
   ///Add one more file in the chain
   void Add(char* inputfile);
   ///Clear the class
   void Clear();
   ///Deconstructor
   virtual ~CorsikaChain() {Clear();}
   ///read all the records from all the files (from beg record to end record)
   Int_t ReadAllRec(int beg=0,int end=0,CorsikaEvent* pevt=0);
   ///read all the records from all the files (from beg file to end file)
   Int_t ReadAllFile(int beg=-1,int end=-1,CorsikaEvent* pevt=0);
   ///read one records
   Int_t ReadRec(CorsikaEvent* pevt=0);
   
};

#endif
