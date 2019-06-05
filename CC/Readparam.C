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

