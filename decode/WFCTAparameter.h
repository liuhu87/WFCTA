#include <math.h>

short SiPM;
short SC;
short Channel;
bool GAIN_MARKER;
char peakpoint;
char maxpoint;
short threshold_single;
short threshold_record;
bool over_record;
bool over_single;
float Four_Point_ADC;
float imagex;
float imagey;
float basehigh;
float baselow;
float qhigh;
float qlow;
float Basehigh;
float Baselow;
float Qhigh;
float Qlow;
float baselineBasehigh;
float baselineBaselow;

int DbTempCount;
int ClbTimeCount;
int ClbTempCount;
int fpga,portnumber;
int start,end,pulseend;

//In subcluster, 23 means F_3DB_2//
int subcluster[8][8] =
    {
        {11,12,13,14,15,16,17,18},
        {21,22,23,24,25,26,27,28},
        {31,32,33,34,35,36,37,38},
        {41,42,43,44,45,46,47,48},
        {51,52,53,54,55,56,57,58},
        {61,62,63,64,65,66,67,68},
        {71,72,73,74,75,76,77,78},
        {81,82,83,84,85,86,87,88}
    };
//In ADRESS, 23 means row_2 column_3//
int ADRESS[8][8] = 
    {
        {44,42,51,53,55,57,48,46},
        {34,32,61,63,65,67,38,36},
        {24,22,71,73,75,77,28,26},
        {14,12,81,83,85,87,18,16},
        {43,41,52,54,56,58,47,45},
        {33,31,62,64,66,68,37,35},
        {23,21,72,74,76,78,27,25},
        {13,11,82,84,86,88,17,15}
    };
