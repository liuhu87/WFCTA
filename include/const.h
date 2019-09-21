#include <math.h>

double Tel_Z;  //Pointing of the telescopes
 
double Tel_A;
 
double C_light = 299.7; //light speed in the air  mm/ns
 
double Ddoor = 2300; //mm , the Width of the door and the maximum size of the mirror
 
double Focus = 2870; //mm, the distance between the mirror and the camera 
 
double Curvature = 5800; //mm // the curvature of the mirrors
 
double Transparency = 0.90; //the trancparency of the front glass in the door
 
double Reflectivity = 0.80; //the mirror's reflectivtiy
 
double QE = 0.25;
 
double Dpmt=25.4; //mm, the Diameter of the pmt
 
const int Pix = 32;
 
double SoptSigma = 30*sqrt(2.)/2.; //the sigma of the spot (mm)
 
double Clusterx = 25.4*32.5; //mm The Size of pmt array
 
double Clustery = 707; //mm
 
double Zcluster0 = Focus; //mm
 
double Zcluster1 = Focus+200; //mm
 
double Zdoor = Zcluster1 + 200; //mm
 
double  Zmirror = Curvature - sqrt(Curvature*Curvature-Ddoor*Ddoor/4.);//the mirror, door and cluster
 
double pmtmap[Pix*Pix][2];

double Efficiency[Pix*Pix];

double NSB = 27; //pe;

double Threshold = 3*sqrt(NSB);

double Npe[Pix*Pix];

int Npmt[Pix*Pix];

Int_t Npix;

double MAXDIST = Dpmt + 5;

double MINPE = 35;

double DSize, DMeanX, DMeanY, DLength, DWidth, DDelta, DAlpha,DDist;
double DSourceX, DSourceY, PhoMax, RMax; 
double DPefar, DPenear, DPixfar, DPixnear;
double DPemean, DPepointing, DPixpointing, DPixmean; 
double DPepointing2, DPixpointing2;

int Neighbors=4;
