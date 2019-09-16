#ifndef __WFCamera__
#define __WFCamera__

#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <map>
//#include "telescopeparameters.h"
///number of pmt in x,y direction
#define PIX 32
#define NSIPM (PIX*PIX)
const double timebinunit=20.; //in unit of nano second
#define MaxTimeBin 1000
using namespace std;

/*!
The class for the camera of telescope
*/
class WCamera
{
  public:
  ///the number of cells in x or y directions
  static int NCell;
  ///the side length of squre sipm,in mm
  static double D_SiPM;
  ///the side of the cell,in mm
  static double D_Cell;

  static double SiPMMAP[NSIPM][2];
  static float NSB;
  public:

  vector<double> TubeSignal;
  vector<double> eTubeSignal;
  //Arrival Time Information
  int NArrival;
  double ArrivalTimeMin;
  double ArrivalTimeMax;
  bool OverFlow;
  long int ArrivalTime[MaxTimeBin];
  double ArrivalCount[NSIPM][MaxTimeBin];
  double ArrivalCountE[NSIPM][MaxTimeBin];
  vector<int> TubeTrigger;
  int TelTrigger;
  //int **Trigger;
  //float **phe;
   
  WCamera();
  ~WCamera();
  static void SetSiPMMAP();
  static void SetNSB( float nsb);
  //void SetCTNumber(int ctnumber);
  static int GetCone(double clusterx, double clustery); 
  static int GetTube(double clusterx, double clustery);
  void Init();
  void ReSet();
  //void ReSet(int ict);
  //void PhotonToTube(int ict, int itube, int outpe);
  //void PhotonCellToTube();
  //void PhotonIntoCone(int ict, int itube, int outpe);
  //void PhotonAfterConeTracing(int ict, int itube, int outpe);
  //void PhotonIntoCell(int ict, int itube, int icell, int outpe);
  static long int FindBin(double time);
  void Fill(int itube,double time=0,double weight=1.0);
  void GetTubeTrigger();
  //void GetTrigger();
  void AddNSB();
  static double GetSiPMX(int itube);
  static double GetSiPMY(int itube);
  //void GetTelescopeTrigger(int CTNumber,float *CT_Zen, float *CT_AZi);
  //void GetEulerMatrix(float TelZ,float TelA);
  //void InverseEuler(double x0, double y0, double z0, double *x, double *y, double *z);
};

#endif
