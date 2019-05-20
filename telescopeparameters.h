#ifndef TELESCOPEPARAMETERS
#define TELESCOPEPARAMETERS
#include <math.h>
#include "TString.h"

//* === COSRIKA CERfile reading parameters === *//
const int PARTICLE_LENGTH_THINNING = 8;
const int PARTICLE_LENGTH_NO_THINNING = 7;
const int N_PARTICLE = 39;
const int N_SUBBLOCK = 21;

const int SUBBLOCK_LENGTH_THINNING = N_PARTICLE * PARTICLE_LENGTH_THINNING;
const int SUBBLOCK_LENGTH_NO_THINNING = N_PARTICLE * PARTICLE_LENGTH_NO_THINNING;
const int RECORD_LENGTH_THINNING = N_SUBBLOCK * SUBBLOCK_LENGTH_THINNING;
const int RECORD_LENGTH_NO_THINNING = N_SUBBLOCK * SUBBLOCK_LENGTH_NO_THINNING;

//* === Telescope geometry and mirror parameters === *//

const double CURVATURE= 5800; //mm // the curvature of the mirrors

//const double TRANSPARENCY= 0.87; //the trancparency of the front glass in the door

//const double REFLECTIVITY= 0.83; //the mirror's reflectivtiy

//const double SHELDED_BY_FRAME = 0.90;

//const double EFFICIENCY = 0.65;//Transparency*Reflectivity*Reflectivity;

//const double QE = 0.25; ///< provided by Ma Lingling 13102015

//const double Dpmt = 26.9; //mm, the Diameter of the pmt

//const double Dpmteff = 21; //mm the effective diameter of the pmt

const double D_DOOR = 2340; //previous value: 2300 ;
                            //mm , the Width of the door and the maximum size of the mirror
const double Hdoor = 2379. ; //mm, the height of the container;  20140926 liujl

const double CLUSTER_X = 922.9; //mm The Size of pmt array

const double CLUSTER_Y = 946.2; //mm

const double RADIUS1_SQUARE = 0.5*D_DOOR*0.5*D_DOOR;

const double RADIUS2_SQUARE = 0.5*CLUSTER_X*0.5*CLUSTER_Y;

const double FOCUS = 2870; //mm, the distance between the mirror and the camera

//const double Focus_D = 2935.; //mm, the distance between the mirror and the camera,  provided by MLL 13102015

const double ZCLUSTER0 = FOCUS; //mm

const double ZCLUSTER1 = FOCUS+230; //mm

const double Zcluster1_D = FOCUS+230; //mm  ///< 13102015

const double ZDOOR = ZCLUSTER1 + 80; //mm

const double Zdoor_D = Zcluster1_D + 80; //mm

const double ZMIRROR = CURVATURE - sqrt(CURVATURE*CURVATURE-D_DOOR*D_DOOR/4.);//the mirror, door and cluster

const double SPOT_SIGMA = 10*sqrt(2.)/2.; //the sigma of the spot (mm)

//* === Mirror geometry parameters === *// ///< Provided by MA Lingling 13102015

const double D = 2300;

const double length = 300.; //mm, the length of the mirror facet of spherical design;  20140926 liujl

const double length_D = 150.; //mm, the length of the mirror facet of Davies design;  20140926 liujl

const double Spot = 10.; //mm, sigma of the spot extension

//* === Camera parameters === *//

const double D_ConeOut = 25.4; //mm, the outer Diameter of the WistonCone

const double D_ConeIn = 24.4; //mm, the inner Diameter of the WistonCone

const double D_SiPM = 15.0; //mm, the side length of squre sipm 

const double D_Cell = 2.5e-2; //mm, the side of the cell  

const int NCell = 600; // the number of cells in x or y directions 

const int  PIX =  32;  ///< number of pmt in x,y direction

const int NSIPM = PIX*PIX;  ///< number of pmt in each camera

//const int CELL_NUMBER = 600;

const double CONE_HEIGHT = 25.27; //mm

const double CONE_ENTRANCE_CIRCLE = 12.2;//mm

const double CONE_EXIT_CIRCLE = 7.5;//mm

const double CUTOFF_ANGLE = 0.66; //radius  37.9*180/PI

const double CONE_REFLECTIVE_EFFICIENCY=0.9;

extern double mirrorx[7][6];
extern double mirrory[7][6];
extern double mirrorz[7][6];
extern double nalfa[7][6];
extern double nbeta[7][6];
extern double spherecenterx[7][6];
extern double spherecentery[7][6];
extern double spherecenterz[7][6];

///< related value
const double sqrt3 = 1.73205080757;
const double sqrt2 = 1.41421356237;
const double pi = 3.14159265359;
const double twopi = 6.28318530718;

//* === Physical constant === *//
const double Qe=-1.6e-4;      ///< Charge value of an electron, unit in muonA*ns
const double C_LIGHT = 299.7; //light speed in the air  mm/ns

#endif // TELESCOPEPARAMETERS

