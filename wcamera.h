#include <iostream>
#include <fstream>
#include <string.h>

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <map>
#include "telescopeparameters.h"
using namespace std;
class WCamera
{
  private:

  double SiPMMAP[NSIPM][2];
  float NSB, TriggerSigma;
  float FixTriggerThreshold;
  int CTNumber;
  float matrix_[3][3];
  public:

  //int **Trigger;
  //float **phe;
  vector<int> TubeSignalIntoCone;
  vector<int> TubeSignalAfterConeTracing;
  vector<int> TubeSignal;
  vector<int> TubeTrigger;
  vector<int> TelTrigger;
  vector<int> TubeSignalInTriggerWindow;
  
  // ** Added by lingling Ma on 2018-1-22        ** //
  // ** to Get the mean arrive time on each SiPM ** //
  vector<float> PeakTime;
  
  map<int, map<int,map<int, bool> > > cell;
  map<int, map<int,map<int, bool> > >::iterator ct_iter;
  map<int, map<int,bool> >::iterator tube_iter;
  map<int,bool>::iterator cell_iter;
 
  // ** Added By Lingling Ma on 2018-1-22        **//
  // ** to count the arrive time of the photons  **// 
  map<int, map<int,map<int,int> > > ArriveTime;
  map<int, map<int,map<int,int> > >::iterator ct_time_iter;
  map<int, map<int,int> >::iterator tube_time_iter;
  map<int,int>::iterator time_time_iter;
  
   
  WCamera();
  ~WCamera();
  void SetSiPMMAP();
  void SetNSB( float nsb);
  void SetTriggerSigma(float triggersigma);
  //* Added on 2017-10-30 by Lingling Ma *//
  //if NSB is not simulted, the fix trigger threshold is used to trigger SiPM *//
  void SetFixTriggerThreshold( float fixtriggerthreshold);
  void SetCTNumber(int ctnumber);
  int GetCone(double clusterx, double clustery); 
  int GetTube(double clusterx, double clustery);
  void Init();
  void ReSet(int ict);
  //void PhotonToTube(int ict, int itube, int outpe);
  void PhotonCellToTube();
  void PhotonIntoCone(int ict, int itube, int outpe);
  void PhotonAfterConeTracing(int ict, int itube, int outpe);
  void PhotonIntoCell(int ict, int itube, int icell, int outpe);
  void GetTubeTrigger(int nsbflag, int fadcflag);
  //void GetTrigger();
  void AddNSB(int fadcflag);
  double GetSiPMX(int itube);
  double GetSiPMY(int itube);
  void GetTelescopeTrigger(int CTNumber,float *CT_Zen, float *CT_AZi);
  void GetEulerMatrix(float TelZ,float TelA);
  void InverseEuler(double x0, double y0, double z0, double *x, double *y, double *z);
  // ** Added By lingling Ma on 2018-1-22        **//
  // ** to Get the mean arrive time on each SiPM **// 
  void GetArriveTime(int itc, int itube, int icell, int itime, int iphoton); 
  void GetPeakTime();
  void GetPhotonInTriggerWindow(double trigger_window);
};
