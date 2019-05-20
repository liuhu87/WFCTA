#ifndef CER_EVENT_H
#define CER_EVENT_H

#include <math.h>
#include "TGraph2D.h"

//class WfctaSimConfig;

class RunInfo
{
public:
    int iRun;
    int nEvents;

    RunInfo();
    ~RunInfo();

    void iRunHeaderShow();
    void iRunEndShow();
};

const double HtoDep_a[5] = {-186.555305,-94.919,0.61289,0.0,0.01128292}; ///< in g/cm^2
const double HtoDep_b[5] = {1222.6562,1144.9069,1305.5948,540.1778,1.}; ///< in g/cm^2
const double HtoDep_c[5] = {994186.38,878153.55,636143.04,772170.16,1e9}; ///< in cm

class CER_Event
{
public:

    int Event_Number_;
    int Event_use_Number_;
    double Primary_Energy_;
    double Primary_id_;
    double Primary_zenith_;
    double Primary_azimuth_;
    double Primary_core_x_;
    double Primary_core_y_;

    double Xmax_;
    double Nmax_;

    double Ob_level_;  ///< in air density

    //* ==== Functions ==== *//

    CER_Event();
    ~CER_Event();
    void Init();

    void EventHeaderShow();
    void EventEndShow();

    double HeightToAtmDep(double height_in_cm);
    double HiLToAtmDep(double height_in_cm, int iLayer);


};

class CER_bunch
{
public:
    float x_;  ///< in mm
    float y_;
    float z_;

    float u_;
    float v_;
    float l_;

    float t_;
    float height_;

    float weight_;
    float nclight_;

    float wavelength_;

    double xc_;   ///< x on pmt cluster plane in mm
    double yc_;   ///< y on pmt cluster plane
    double u1_;   // cosin directions of reflected photons;
    double v1_;

    double time_raytrace_;    ///< time for RayTrace from telescope area to the camera.

    CER_bunch();
    ~CER_bunch();

    void Init();
    void SetCERBunch(float,float,float,float,float,float);
    void IntoTelArea(double telx,double tely);

};

#endif // CER_EVENT_H
