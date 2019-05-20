#include "cer_event.h"
#include <iostream>
#include <math.h>
#include <map>
#include "stdio.h"
#include <TRandom.h>

using namespace std;

RunInfo::RunInfo()
{

}

RunInfo::~RunInfo()
{

}

void RunInfo::iRunHeaderShow()
{
    cout << Form("==>> Run       %10d    Header read <<==",iRun) << endl;

}
void RunInfo::iRunEndShow()
{
    cout << Form("==>> Run       %10d     End   read <<==",iRun) << endl;
}


CER_Event::CER_Event()
{
    Event_Number_ = 0;
    Primary_Energy_= 0.;
    Primary_id_ = 0.;
    Primary_zenith_ = 0.;
    Primary_azimuth_ = 0.;

    //Shower_Core_Detected[0] = 0.;
    //Shower_Core_Detected[1] = 0.;

}

CER_Event::~CER_Event()
{

}


void CER_Event::Init(){

    Event_Number_       = 0;
    Event_use_Number_   = 0;
    Primary_Energy_     = 0;
    Primary_id_        = 0;
    Primary_zenith_    = 0;
    Primary_azimuth_   = 0;
    Primary_core_x_   = 0;
    Primary_core_y_   = 0;
    Xmax_     = 0.;
    Nmax_     = 0.;
    Ob_level_ = 0;  ///< in air density


}




void CER_Event::EventHeaderShow()
{
    cout << Form("==>> Event %10d use %2d Header read <<==",Event_Number_,Event_use_Number_) << endl;
}

void CER_Event::EventEndShow()
{
    cout << Form("==>> Event %10d use %2d  End   read <<==",Event_Number_,Event_use_Number_) << endl;
}


double CER_Event::HeightToAtmDep(double height_in_cm)
{
    if(height_in_cm<0.){
        return -1;
    }
    if(height_in_cm>=0. && height_in_cm<4e5){
        return HiLToAtmDep(height_in_cm,1);
    }
    if(height_in_cm>=4e5 && height_in_cm<1e6){
        return HiLToAtmDep(height_in_cm,2);
    }
    if(height_in_cm>=1e6 && height_in_cm<4e6){
        return HiLToAtmDep(height_in_cm,3);
    }
    if(height_in_cm>=4e6 && height_in_cm<1e7){
        return HiLToAtmDep(height_in_cm,4);
    }
    else {
        return HiLToAtmDep(height_in_cm,5);
    }

}

double CER_Event::HiLToAtmDep(double height_in_cm, int iLayer)
{
    if(iLayer<1 || iLayer>5){
        cerr << "iLayer value should be within 1 to 5." << endl;
        return -1;
    }
    switch (iLayer){
    case 1:
        return HtoDep_a[1]+HtoDep_b[1]*exp(-height_in_cm/HtoDep_c[1]);
        break;
    case 2:
        return HtoDep_a[2]+HtoDep_b[2]*exp(-height_in_cm/HtoDep_c[2]);
        break;
    case 3:
        return HtoDep_a[3]+HtoDep_b[3]*exp(-height_in_cm/HtoDep_c[3]);
        break;
    case 4:
        return HtoDep_a[4]+HtoDep_b[4]*exp(-height_in_cm/HtoDep_c[4]);
        break;
    case 5:
        return HtoDep_a[5]-HtoDep_b[5]*(-height_in_cm/HtoDep_c[5]);
        break;
    default:
        return -1;
        break;
    }
}


CER_bunch::CER_bunch()
{

}

CER_bunch::~CER_bunch()
{

}

void CER_bunch::Init(){
    x_ = 0;  ///< in mm
    y_ = 0;
    z_ = 0;

    u_ = 0;
    v_ = 0;
    l_ = 0;

    t_ = 0;
    height_ = 0;

    weight_ = 0;
    nclight_ = 0;

    wavelength_ = 0;

    xc_ = 0;   ///< x on pmt cluster plane in mm
    yc_ = 0;   ///< y on pmt cluster plane

    time_raytrace_ = 0;    ///< time for RayTrace from telescope area to the camera.
}


void CER_bunch::SetCERBunch(float x,float y,float u,float v,float t,float height)
{
    x_= x; y_= y; u_= u; v_= v; t_= t; height_= height;
}

void CER_bunch::IntoTelArea(double telx,double tely)
{
    x_ = x_ - telx;
    y_ = y_ - tely;
    z_ = 0.;
    // l should be -sqrt(1-u*u-v*v) for the coordinate systrm in corsika is strange  --2016.3.30 Baiyang Bi
    l_ = sqrt(1 - u_*u_ - v_*v_);
}


