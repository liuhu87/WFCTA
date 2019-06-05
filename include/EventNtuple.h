#ifndef __EventNtuple__
#define __EventNtuple__
#include "CorsikaEvent.h"
#include "WFCTAEvent.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TList.h"

#define NBranch 2
/*!
The class is used for filling all the events and make some histograms
*/

class EventNtuple{
   private:
   static EventNtuple* _Head;
   public:
   ///saving style(fillstyle=0:all,fillstyle=1:in a tree,fillstyle=2: in histograms,fillstyle=3:MC counts)
   static int fillstyle;

   ///the file to save the ntuple
   static TFile* fout;

   ///pointer to the tree, which will contain all the events information
   TTree* _Tree;
   ///Split Level for the tree
   static const int branchSplit;
   ///the branch name of the tree
   static char branchname[NBranch][20];

   ///some simple histograms
   TList* hlist;

   public:
   ///convert from particle id to A and Z
   static bool FromIDtoAZ(int PID,int &AA,int &ZZ);

   public:
   ///get the static class pointer
   static EventNtuple* GetHead(char* filename=(char*)"EventNtuple.root",int style=0);
   ///some tools for convert coordinate
   static void ConvertCoor(double pos[3],double pp[3],double InCoo[3],double InDir[3],double *OutCoo,double *OutDir);
   static void ConvertThetaPhi(double pp[3],double &theta,double &phi);
   static void ConvertCoor(double pos[3],double theta,double phi,double InCoo[3],double InDir[3],double *OutCoo,double *OutDir);
   ///init all the tree and histograms
   void Init();
   ///clear all the tree and histograms
   void Clear();
   ///write all the tree and histograms
   void Write();
   ///default constructor
   EventNtuple(char* filename,int style);
   ///default deconstructor
   virtual ~EventNtuple();
   ///Fill all the events to tree or histograms
   void Fill(TSelector** pevt);
};

#endif
