#ifndef WCAMERA_H
#define WCAMERA_H

#include <math.h>

class WCamera{
  public:
  double GetPMTMAP();
  double GetPixX(int ipix);
  double GetPixY(int ipix);
  int  GetPixes();
  double GetDist(int tube,double sourcex,double sourcey);
  private:
  double pmtmap[32*32][2];
};
#endif // WCAMERA_H
