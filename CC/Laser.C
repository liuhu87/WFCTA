#include "Laser.h"
#include "stido.h"
#include "iostream.h"
#include "math.h"

TRandom* Atmosphere::prandom = 0;
double Atmosphere::aod_air = 0;
double Atmosphere::aod_aerosol = 0;
double Atmosphere::scat_air = 0;
double Atmosphere::scat_aerosol = 0;
double distance = 0;

void Atmosphere::Init(){
   prandom = grandom;
   aod_air = 0;
   aod_aerosol = 0;
   scat_air = 0;
   scat_aerosol = 0;
}

void Atmosphere::SetParameters(char* filename){
;
}

void Atmosphere::RayleighScatter( double theta_in, double phi_in, double wavelength, double &theta_out, double &phi_out){
   double theta_scat, phi_scat;
   RayScatterAngle( wavelength, theta_scat, phi_scat);
   //from theta/phi_in to theta/phi_out according to theta/phi_scat
   theta_out = theta_in + theta_scat;
   phi_out = phi_in + phi_scat;
    
}

void Atmosphere::MieScatter(double theta_in, double phi_in, double wavelength, double &theta_out, double &phi_out){
   double theta_scat, phi_scat;
   MieScatterAngle(wavelength, theta_scat, phi_scat);
   //from theta/phi_in to theta/phi_out according to theta/phi_scat
   theta_out = theta_in + theta_scat;
   phi_out = phi_in + phi_scat;
}

void Atmosphere::RayScatterAngle(double wavelength, double &theta, double &phi){
///
;
}

void Atmosphere::MieScatterAngle(double wavelength, double &theta, double &phi){
///
;
}

double Atmosphere::FreePathLength(){
   

   return 0;
}

bool Atmosphere::IsScattering(){
}

double Atmosphere::probability(){
   prob(distance) = (aod_air + aod_aerosol) * exp((-(aod_air + aod_aerosol)) * distance );//distance is the variable
}


/*\
The class for laser photon propagation
*/

TRandom3* Laser::prandom = 0;
double Laser::intensity = 2;//mj
double Laser::intensity_err = 0;
double Laser::wavelength0 = 355;//nm
double Laser::wavelength0_err = 0;
double Laser::frequency = 1;
double Laser::spotrange = 0.001;//mm, the range of the initial laser spot
double Laser::divergence = 0.0573;


double Laser::coor[3] = {2000, 0, 0};//m, position of the laser machine
double Laser::dir[2] = {30，0};//elevation & deflection angle of the center of the beam, clockwise, up is positive

int Laser::count_gen = 0;
double Laser::time_gen = 0;
int Laser::ievent_gen = 0;
double Laser::wavelength_gen = 0;

double Laser::coor_gen[3] = {0};//m, position of the beam, guas random 
double Laser::dir_gen[2] = {0};//offset of the center, guas random 
double Laser::coor_tel[3] = {0};
double Laser::dir_tel[2] = {0};


void Laser::Init(){
   prandom = grandom;
   double intensity;
   double intensity_err;
   double wavelength0;
   double wavelength0_err;
   double frequency;
}

void Lsaer::PositionGen(){
   coor_gen[0] = Double_t TRandom::Gaus( Double_t mean = coor[0], Double_t sigma = spotrange );//x coordinate of the laser machine 
   coor_gen[1] = Double_t TRandom::Gaus( Double_t mean = coor[1], Double_t sigma = spotrange );//y coordinate of the laser machine 
   coor_gen[2] = Double_t TRandom::Gaus( Double_t mean = coor[2], Double_t sigma = spotrange );//z coordinate of the laser machine 
}

double Laser::DirectionGen(){
   Laser::dir_gen[0] = Double_t TRandom::Gaus( Double_t mean = 0, Double_t sigma = divergence );//offset of the center, guas random, polar coordinates，polar axis is the center of the initial laser beam
   //Laser::dir_gen[1] = Double_t TRandom::Gaus( Double_t mean = 0, Double_t sigma = 0.001 );
    
}

double Laser::WaveLengthGen(){
      
   ;
}

bool Laser::LaserGen(){

   :
}

int Laser::EventGen(int ngen=100){

   :
}

int Laser::Propagate(){
   
   ;
}
