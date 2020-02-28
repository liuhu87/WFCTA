{
   gSystem->Load("lib/lib.so");
   RotateDB::jdebug=10;

   int Li=2;
   //int time_in=CommonTools::Convert(20200210231855.);
   //int time_in=1579176379;
   int time_in=1575571453;
   int hour=CommonTools::TimeFlag(time_in,4);
   int min=CommonTools::TimeFlag(time_in,5);
   int sec=CommonTools::TimeFlag(time_in,6);

   //use the command below to get the record corresponding to time=time_in and laser transmited from Li
   //int index0=RotateDB::GetHead()->GetEleAzi(time_in,Li,1);
   RotateDB::GetHead()->LoadData(time_in,Li);
   //RotateDB::GetHead()->ProcessAll();
   //after the ProcessAll function,use the functions below to get all kinds of information
   //DumpInfo(), which print out all the variables
   //GetLaserSwith(), which returns wheather the laser is on or off
   //GetHeight(), which return the position of laser in z axis 
   //GetElevation(), which returns the elevation angle when transmited
   //GetAzimuth(), which returns the azimuth angle when transmited

   //printf("index=%d time=%d hms=%d:%d:%d elevation=%.2lf azimuth=%.2lf\n",index0,time_in,hour,min,sec,RotateDB::GetHead()->GetElevation(),RotateDB::GetHead()->GetAzimuth());
}
