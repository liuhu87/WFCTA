#include "Readparam.h"
static int  WLFlag = 0;                        //@< flag of wavelength
static int  CTNumber = 1;                      //@< number of telescopes 
static int ThinFlag = 0;                       //@< flag of thin option
static float **CT_Position;                    //@< the positions and pointings of telescopes
static int FilterFlag = 0;                     //@< Will the Filter be added 
static int PointErrorFlag = 0;                 //@< Will the Point error of each miror be added 
static int PointError = 0;                     //@< the value of the point error (1 sigma) 
static int MirrorGeometry = 1;                 //@< the geomtry of mirrors 
static float MirrorSize[2];                    //@< the total mirror size in x and y directions 
static float MirrorSpot=0.;                    //@<the light spot of the mirrors 
static char SimMode[PATH_MAX_LENGTH];          //@< the simulation mode
static float Fadc_bin = 12;                    //@< the total bin number of fadc
static float Fadc_bin_length = 20;             //@< the bin length of each fadc
static float Nsb = 0.;                         //@< the night sky background
static float aod_air = 0.;                     //@< extinction coefficient of air
static float aod_aerosol = 0;                  //@<extinction coefficient of aerosol
static float scat_air = 0;                     //@<scattering coefficient of air
static float scat_aerosol = 0;                 //@<scattering coefficient of aerosol
static float laser_coo[3]={0,0,0};             //@<laser generate position (in cm)
static float laser_dir[2]={0,0};               //@<laser generate direction (in degree)
static float laser_intensity[2] = {0,0};       //@<laser intensity (in mJ)
static float laser_wavelength[2] = {0,0};      //@<laser wavelength (in nm)
static float laser_frequency = 1;              //@<laser frequency (in Hz)
static float laser_pulsetime = 7;              //@<laser pulsetime (in ns)
static float laser_spotrange[2] = {1,1};       //@<laser spotrange (in mm)
static float laser_divergence[2] = {1,1};      //@<laser divergence (in mrad)
using namespace std;
WReadConfig::WReadConfig()
{
}
WReadConfig::~WReadConfig()
{
}
void WReadConfig::readparam(char * filename)
{
  char line[LINE_MAX_LENGTH];    //@< line to get from the stdin
  char token[ITEM_MAX_LENGTH];   //@< a single token
  char token2[256];              //@< a single token
  int i, j, k;                   //@< dummy counters
  char filename_tmp[PATH_MAX_LENGTH];
  float x, y, z, azi, zen;
  ifstream ifile;  

  CT_Position = new float *[MAX_NUMBER_OF_CTS];
  for (i = 0; i < MAX_NUMBER_OF_CTS ; i++){
     CT_Position[i] = new float[5];
  }  

  if ( filename != NULL )
     ifile.open( filename );

 // loop till the "end" directive is reached
 int is_end = FALSE;

 while (! is_end) {
  
   // get line from file or stdin
   if ( filename != NULL )
      ifile.getline(line, LINE_MAX_LENGTH);
   // skip comments (start with '#')
   if (line[0] == '#')
      continue;

   // look for each item at the beginning of the line
   for(i = 0; i <= end_file; i++){
      sscanf(line, "%s", token);
      if(strcmp(token, ITEM_NAMES[i]) == 0)
        break;
   }
   //if it is not a valid line, just ignore it
   if (i == end_file+1) {
      cerr << "ERROR: Unknown token in [" << line << "]\n";
      exit(-1);
   }
   // case block for each directive
   switch ( i ) {

     case wavelength:
         sscanf(line, "%s  %s", token, token2);   
         if(strcmp(token2,"T")==0||strcmp(token2,"TRUE")==0) 
           WLFlag = TRUE;
         break;
     case thin:
         sscanf(line, "%s  %s", token, token2);
         if(strcmp(token2,"T")==0||strcmp(token2,"TRUE")==0)
           ThinFlag = TRUE;
         break;  
     case filter:
         sscanf(line, "%s  %s", token, token2);
         if(strcmp(token2,"T")==0||strcmp(token2,"TRUE")==0)
           FilterFlag = TRUE;
         break;
     case ct_number:
         sscanf(line, "%s  %d", token, &k);
         CTNumber = k;
         break;
     case ct_position:
         sscanf(line, "%s  %d %f %f %f %f %f", token, &k, &x,&y,&z,&zen,&azi);
         CT_Position[k][0] = x;
         CT_Position[k][1] = y;
         CT_Position[k][2] = z;
         CT_Position[k][3] = zen;
         CT_Position[k][4] = azi;
         break;
     case mirror_size:
         sscanf(line, "%s %f %f", token,&x,&y);
         MirrorSize[0] = x;
         MirrorSize[1] = y;
         break;
     case mirror_spot:
         sscanf(line, "%s %f", token,&x);
         MirrorSpot = x;
         break;
     case mirror_geomtry:
         sscanf(line, "%s %d", token,&k); 
         MirrorGeometry = k;
         break;
     case mirror_pointerror:
         sscanf(line, "%s %s %f", token,token2,&x); 
         if(strcmp(token2,"T")==0||strcmp(token2,"TRUE")==0){ 
            PointErrorFlag = TRUE;
            PointError = x ;
         }
         break;
     case SIMMODE:
         sscanf(line, "%s %s", token, token2);
         strcpy(SimMode,token2);
         break;
     case fadc_bin:
         sscanf(line,"%s %d",token,&k);
         Fadc_bin = k;
         break;
     case fadc_bin_length:
         sscanf(line,"%s %d",token,&k);
         Fadc_bin_length = k;
         break;
     case nsb:
         sscanf(line,"%s %f",token,&x);
         Nsb = x;
         break;
     case Aod_air:
         sscanf(line,"%s %f",token,&x);
         aod_air = x;
         break;
     case Aod_aerosol:
         sscanf(line,"%s %f",token,&x);
         aod_aerosol = x;
         break;
     case Scat_air:
         sscanf(line,"%s %f",token,&x);
         scat_air = x;
         break;
     case Scat_aerosol:
         sscanf(line,"%s %f",token,&x);
         scat_aerosol = x;
         break;
     case Laser_coo:
         sscanf(line,"%s %f %f %f",token,&x,&y,&z);
         laser_coo[0] = x;
         laser_coo[1] = y;
         laser_coo[2] = z;
         break;
     case Laser_dir:
         sscanf(line,"%s %f %f",token,&x,&y);
         laser_dir[0] = x;
         laser_dir[1] = y;
         break;
     case Laser_intensity:
         sscanf(line,"%s %f %f",token,&x,&y);
         laser_intensity[0] = x;
         laser_intensity[1] = y;
         break;
     case Laser_wavelength:
         sscanf(line,"%s %f %f",token,&x,&y);
         laser_wavelength[0] = x;
         laser_wavelength[1] = y;
         break;
     case Laser_frequency:
         sscanf(line,"%s %f",token,&x);
         laser_frequency = x;
         break;
     case Laser_pulsetime:
         sscanf(line,"%s %f",token,&x);
         laser_pulsetime = x;
         break;
     case Laser_spotrange:
         sscanf(line,"%s %f %f",token,&x,&y);
         laser_spotrange[0] = x;
         laser_spotrange[1] = y;
         break;
     case Laser_divergence:
         sscanf(line,"%s %f %f",token,&x,&y);
         laser_divergence[0] = x;
         laser_divergence[1] = y;
         break;
     case end_file:
        // the end of the input card
        is_end = TRUE;
        break;                  
   }
 }
 return; 
}


int WReadConfig::GetCTNumber()
{
  return  CTNumber;
}

float WReadConfig::GetCTPosition(int ict, int i)
{
  //if i=0, ctx
  //if i=1, cty
  //if i=2, ctz 
  //if i=3, ct azi
  //if i=4, ct zen;
  return CT_Position[ict][i];
}

int WReadConfig::GetMirrorGeometry()
{
  return  MirrorGeometry;
}

void WReadConfig::GetMirrorSize(float *x , float *y)
{
  *x = MirrorSize[0];
  *y = MirrorSize[1];
}

float WReadConfig::GetMirrorSpot()
{ 
   return MirrorSpot;
}

int WReadConfig::GetMirrorPointError()
{
  return PointError;
}

int WReadConfig::GetMirrorPointErrorFlag()
{
  return PointErrorFlag;
}

int WReadConfig::GetFilterFlag()
{
  return FilterFlag;
}

int WReadConfig::GetThinFlag()
{
  return ThinFlag;
}

int WReadConfig::GetWLFlag()
{
  return WLFlag;
}
int WReadConfig::GetFadcBins()
{
  return Fadc_bin;
}
int WReadConfig::GetFadcLength()
{
  return Fadc_bin_length;
}
float WReadConfig::GetNSB()
{
  return Nsb;
}
float WReadConfig::Getaod_air(){
  return aod_air;
}
float WReadConfig::Getaod_aerosol(){
  return aod_aerosol;
}
float WReadConfig::Getscat_air(){
  return scat_air;
}
float WReadConfig::Getscat_aerosol(){
  return scat_aerosol;
}
float WReadConfig::GetLaserCoo(int i){
  return laser_coo[i];
}
float WReadConfig::GetLaserDir(int i){
  return laser_dir[i];
}
float WReadConfig::GetLaserIntensity(){
  return laser_intensity[0];
}
float WReadConfig::GetLaserIntensityErr(){
  return laser_intensity[1];
}
float WReadConfig::GetLaserWavelength(){
  return laser_wavelength[0];
}
float WReadConfig::GetLaserWavelengthErr(){
  return laser_wavelength[1];
}
float WReadConfig::GetLaserFrequency(){
  return laser_frequency;
}
float WReadConfig::GetLaserPulsetime(){
  return laser_pulsetime;
}
float WReadConfig::GetLaserSpotrange(int i){
  return laser_spotrange[i];
}
float WReadConfig::GetLaserDivergence(int i){
  return laser_divergence[i];
}

