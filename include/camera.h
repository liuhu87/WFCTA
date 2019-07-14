#include <stdio.h>
#include <stdint.h>
#include <math.h>

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

/********************************************************************************
 * void SC_Channel2SiPM(short F_DB, short mChannel, short *mSiPM)               *
 * change FPGA_DB & Channel to sipm.                                            *
 * F_DB==25 means F5_DB2; Channel==1 means channel_1; range of sipm is [0,1023] *
 ********************************************************************************/
void SC_Channel2SiPM(short F_DB, short mChannel, short *mSiPM)
{
    int FPGA;
    int DB;
    double SC_X;
    double SC_Y;
    double Channel_X;
    double Channel_Y;
    FPGA = F_DB%10;
    DB = F_DB/10;
    SC_X = ADRESS[DB-1][FPGA-1]%10;
    SC_Y = ADRESS[DB-1][FPGA-1]/10;

    if(mChannel<9)
    {
      if((mChannel)%2==1){
        Channel_Y = (SC_X-1)*4.+1;
        Channel_X = 33-((SC_Y-1)*4.+(mChannel+1)/2.);
      }
      if((mChannel)%2==0){
        Channel_Y = (SC_X-1)*4.+2;
        Channel_X = 32.5-((SC_Y-1)*4.+(mChannel)/2.);
      }
    }
    if(mChannel>=9)
    {
      if((mChannel)%2==1){
        Channel_Y = (SC_X-1)*4.+3;
        Channel_X = 33-((SC_Y-1)*4.+(mChannel-7)/2.);
      }
      if((mChannel)%2==0){
        Channel_Y = (SC_X-1)*4.+4;
        Channel_X = 32.5-((SC_Y-1)*4.+(mChannel-8)/2.);
      }
    }
    if(int(Channel_X*2)%2==1){Channel_X = Channel_X+0.5;}
    else {Channel_X = Channel_X;}
    Channel_Y = 32-Channel_Y;
    *mSiPM = int(1023-(Channel_Y*32+Channel_X-1));
}


short sipm2SC[1024] =
{
  82,82,82,82,72,72,72,72,62,62,62,62,52,52,52,52,13,13,13,13,23,23,23,23,33,33,33,33,43,43,43,43,
  82,82,82,82,72,72,72,72,62,62,62,62,52,52,52,52,13,13,13,13,23,23,23,23,33,33,33,33,43,43,43,43,
  82,82,82,82,72,72,72,72,62,62,62,62,52,52,52,52,13,13,13,13,23,23,23,23,33,33,33,33,43,43,43,43,
  82,82,82,82,72,72,72,72,62,62,62,62,52,52,52,52,13,13,13,13,23,23,23,23,33,33,33,33,43,43,43,43,
  42,42,42,42,32,32,32,32,22,22,22,22,12,12,12,12,53,53,53,53,63,63,63,63,73,73,73,73,83,83,83,83,
  42,42,42,42,32,32,32,32,22,22,22,22,12,12,12,12,53,53,53,53,63,63,63,63,73,73,73,73,83,83,83,83,
  42,42,42,42,32,32,32,32,22,22,22,22,12,12,12,12,53,53,53,53,63,63,63,63,73,73,73,73,83,83,83,83,
  42,42,42,42,32,32,32,32,22,22,22,22,12,12,12,12,53,53,53,53,63,63,63,63,73,73,73,73,83,83,83,83,
  81,81,81,81,71,71,71,71,61,61,61,61,51,51,51,51,14,14,14,14,24,24,24,24,34,34,34,34,44,44,44,44,
  81,81,81,81,71,71,71,71,61,61,61,61,51,51,51,51,14,14,14,14,24,24,24,24,34,34,34,34,44,44,44,44,
  81,81,81,81,71,71,71,71,61,61,61,61,51,51,51,51,14,14,14,14,24,24,24,24,34,34,34,34,44,44,44,44,
  81,81,81,81,71,71,71,71,61,61,61,61,51,51,51,51,14,14,14,14,24,24,24,24,34,34,34,34,44,44,44,44,
  41,41,41,41,31,31,31,31,21,21,21,21,11,11,11,11,54,54,54,54,64,64,64,64,74,74,74,74,84,84,84,84,
  41,41,41,41,31,31,31,31,21,21,21,21,11,11,11,11,54,54,54,54,64,64,64,64,74,74,74,74,84,84,84,84,
  41,41,41,41,31,31,31,31,21,21,21,21,11,11,11,11,54,54,54,54,64,64,64,64,74,74,74,74,84,84,84,84,
  41,41,41,41,31,31,31,31,21,21,21,21,11,11,11,11,54,54,54,54,64,64,64,64,74,74,74,74,84,84,84,84,
  88,88,88,88,78,78,78,78,68,68,68,68,58,58,58,58,15,15,15,15,25,25,25,25,35,35,35,35,45,45,45,45,
  88,88,88,88,78,78,78,78,68,68,68,68,58,58,58,58,15,15,15,15,25,25,25,25,35,35,35,35,45,45,45,45,
  88,88,88,88,78,78,78,78,68,68,68,68,58,58,58,58,15,15,15,15,25,25,25,25,35,35,35,35,45,45,45,45,
  88,88,88,88,78,78,78,78,68,68,68,68,58,58,58,58,15,15,15,15,25,25,25,25,35,35,35,35,45,45,45,45,
  48,48,48,48,38,38,38,38,28,28,28,28,18,18,18,18,55,55,55,55,65,65,65,65,75,75,75,75,85,85,85,85,
  48,48,48,48,38,38,38,38,28,28,28,28,18,18,18,18,55,55,55,55,65,65,65,65,75,75,75,75,85,85,85,85,
  48,48,48,48,38,38,38,38,28,28,28,28,18,18,18,18,55,55,55,55,65,65,65,65,75,75,75,75,85,85,85,85,
  48,48,48,48,38,38,38,38,28,28,28,28,18,18,18,18,55,55,55,55,65,65,65,65,75,75,75,75,85,85,85,85,
  87,87,87,87,77,77,77,77,67,67,67,67,57,57,57,57,16,16,16,16,26,26,26,26,36,36,36,36,46,46,46,46,
  87,87,87,87,77,77,77,77,67,67,67,67,57,57,57,57,16,16,16,16,26,26,26,26,36,36,36,36,46,46,46,46,
  87,87,87,87,77,77,77,77,67,67,67,67,57,57,57,57,16,16,16,16,26,26,26,26,36,36,36,36,46,46,46,46,
  87,87,87,87,77,77,77,77,67,67,67,67,57,57,57,57,16,16,16,16,26,26,26,26,36,36,36,36,46,46,46,46,
  47,47,47,47,37,37,37,37,27,27,27,27,17,17,17,17,56,56,56,56,66,66,66,66,76,76,76,76,86,86,86,86,
  47,47,47,47,37,37,37,37,27,27,27,27,17,17,17,17,56,56,56,56,66,66,66,66,76,76,76,76,86,86,86,86,
  47,47,47,47,37,37,37,37,27,27,27,27,17,17,17,17,56,56,56,56,66,66,66,66,76,76,76,76,86,86,86,86,
  47,47,47,47,37,37,37,37,27,27,27,27,17,17,17,17,56,56,56,56,66,66,66,66,76,76,76,76,86,86,86,86,
};

short sipm2Channel[1024] =
{
   1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7,
   2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8,
   9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15,
  10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,
   1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7,
   2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8,
   9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15,
  10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,
   1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7,
   2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8,
   9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15,
  10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,
   1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7,
   2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8,
   9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15,
  10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,
   1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7,
   2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8,
   9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15,
  10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,
   1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7,
   2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8,
   9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15,
  10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,
   1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7,
   2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8,
   9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15,
  10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,
   1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7,
   2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8, 2, 4, 6, 8,
   9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15, 9,11,13,15,
  10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,10,12,14,16,
};

/*********************************************************************************
 * void SiPM2SC_Channel(short mSiPM, short *mSC, short *mChannel)                *
 * change sipm to FPGA_DB & Channel                                              *
 * mSC==25 means F5_DB2; mChannel==1 means channel_1; range of mSiPM is [0,1023] *
 *********************************************************************************/
void SiPM2SC_Channel(short mSiPM, short *mSC, short *mChannel)
{
  *mSC = sipm2SC[mSiPM];
  *mChannel = sipm2Channel[mSiPM];
}

/*******************************************************************************
 * void SC_Channel2eSiPM(short fpga, short db, short channel, short *sipm)     *
 * change FPGA & DB & Channel to esipm                              *
 * range of esipm is [1,1024]                                                   *
 *******************************************************************************/
void SC_Channel2eSiPM(short fpga, short db, short channel, short *sipm)
{
  *sipm = (fpga-1)*128+(db-1)*16+channel;
}

/*******************************************************************************
 * void eSiPM2SC_Channel(short mSiPM, short *mSC, short *mChannel)             *
 * change esipm to FPGA_DB & Channel                                           *
 * range of esipm is [1,1024]                                                  *
 *******************************************************************************/
void eSiPM2SC_Channel(short mSiPM, short *mSC, short *mChannel)
{
  short modfpga;
  short modfpga;

  short fpga;
  short db;
  short channel;

  modfpga = mSiPM % 128;
  moddb = (mSiPM % 128) % 16;
  fpga = mSiPM / 128;
  db = (mSiPM % 128) / 16;

  if(modfpga!=0){
    fpga += 1;
    if(moddb!=0){
      db += 1;
      channel = (mSiPM % 128) % 16;
    }
    else{
      channel = 16;
    }
  }
  else{
    db = 8;
    channel = 16;
  }

  *mSC = db*10+fpga;
  *mChannel = channel;
}














