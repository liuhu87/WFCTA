#include <iostream>
#include <fstream>
#include <string.h>

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include "telescopeparameters.h"
#include "TGraph2D.h"
using namespace std;

class CER_bunch;
class WTelescope
{
  private:
    float  MirrorSpot, MirrorPointError;
    int MirrorPointErrorFlag, MirrorGeometry;
    double TelZ_;  ///<
    double TelA_;  ///< pointing of the telescope;
    double matrix_[3][3];  ///< Euler matrix for coordinate transformation
    double mirrorx[7][6];
    double mirrory[7][6];
    double mirrorz[7][6];
    double nalfa[7][6];
    double nbeta[7][6];
    double spherecenterx[7][6];
    double spherecentery[7][6];
    double spherecenterz[7][6];
    //Update by Lingling Ma 2019-05-08
    //vector <double> Transmissivity;
    TGraph2D *Transmissivity; //[wl][incident angle] wl is from 600nm to 300 nm with 5 nm step. incident angle is from 0 degrees to 35 degrees with 5 degrees step. 
    double filter_wl_max, filter_wl_min;
    double filter_angle_max, filter_angle_min;
    double Reflectivity[55];
    double mirror_wl_max, mirror_wl_min, mirror_wl_step;
  public:
    WTelescope();
    ~WTelescope();
    double total_photon_;
    double Telx_;  ///<
    double Tely_;  ///< position of the telescope;
    void SetMirrorSpot(float mirror_spot);
    void SetMirrorPointErrorFlag(int point_error_flag);
    void SetMirrorPointError(float point_error);
    void SetMirrorGeometry(int mirror_geometry);
    void SetMirror();
    void SetMirror(int i);
    void SetMirrorPointError(int flag,double error);
    void altercoor(double *alfa,double *ox,double *oy,double *nx,double *ny);
    void SetXY(double x,double y);
    void SetPointing(double zenith, double azimuth);
    void SetEulerMatrix(double theta,double phi); 
    int IncidentTel(double x,double y,double lengthx,double lengthy); 
    void SetTransmissivity();   //added by lingling Ma 2018.1.10
    void SetReflectivity();
    int GetTransmissivity(double wavelength, double angle);    
    int Reflected(double wavelength);
    double Euler(double x0, double y0, double z0, double *x, double *y, double *z);  ///< coordinate transformation
    double Euler(CER_bunch*bunch);
    double InverseEuler(double x0, double y0, double z0, double *x, double *y, double *z);
    int WhichMirror(double x, double y, double z, double *deltax, double *deltay,double *deltaz,int *ii, int *mm);
    double Plane(double z,double x0,double y0,double z0,double m2,double n2,double l2,double *x,double *y);
    double Sphere(double Z0,double R, double x0, double y0, double z0, double m, double n, double l, double *x, double *y, double *z);
    double RayTrace(double spotsize,double x0, double y0, double z0, double m1, double n1, double l1,double *xc, double *yc, double *t, double *u, double *v);

    double GetReflected(double x1,double y1,double z1,double x,double y,double z,double *m2,double *n2,double *l2);
};

