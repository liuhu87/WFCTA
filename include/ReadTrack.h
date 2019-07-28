#ifndef __ReadTrack__
#define __ReadTrack__

#include <fstream>
#include <iostream>
#include "math.h"
#include "common.h"
#include <TCanvas.h>
#include <TView3D.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TSystem.h>
#include <TObjArray.h>
using namespace std;

class CorsikaEvent;

union TrackBuff{
   float f[10];
   char  c[4*10];
};
class ReadTrack{
   private:
   static ReadTrack* _Head;
   ifstream* fin;
   ///the buffer to save one record data
   TrackBuff recbuff;
   public:
   static bool DoPlot;
   static int jdebug;
   static int particle;
   static float elimit[2];
   static float climit[3][2];
   static float tlimit[2];
   static float IniRange[4][2];
   int nrec;

   vector<int> arrid;
   vector<float> arren;
   vector<float> arrc1[3];
   vector<float> arrt1;
   vector<float> arrc2[3];
   vector<float> arrt2;

   TObjArray* plot;
   float plotrange[4][2];

   public:
   static void SetHead(ReadTrack* head);
   static ReadTrack* GetHead();
   void Init();
   void Reset();
   void Release();
   ReadTrack(const char* inputfile);
   ~ReadTrack() {Release();}
   bool Exist();
   int ReadRec();
   int ReadAll(int beg=0,int end=0);
   void Copy(CorsikaEvent* pevt);
   static int Color(int partid);
   static int Style(int partid);
   static int Width(int partid);
   void Draw(TCanvas* cc=0,const char* option="al");
};

#endif
