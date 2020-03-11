{
   gSystem->Load("/afs/ihep.ac.cn/users/h/hliu/Documents/LHAASO/WFCTA/lib/lib.so");
   RotateDB::jdebug=0;

   int Li=2;
   //int time_in=CommonTools::Convert(20200210231855.);
   //int time_in=1579176379;
   int time_in=1583560656;
   int hour=CommonTools::TimeFlag(time_in,4);
   int min=CommonTools::TimeFlag(time_in,5);
   int sec=CommonTools::TimeFlag(time_in,6);

   //use the command below to get the record corresponding to time=time_in and laser transmited from Li
   //int index0=RotateDB::GetHead()->GetEleAzi(time_in,Li,1);
   //RotateDB::GetHead()->ProcessAll();
   //long int res=RotateDB::GetHead()->LoadData2(time_in,Li);
   //RotateDB::GetHead()->ProcessAll2();
   double temp_buff[5];
   long int res=(long int)RotateDB::GetHead()->GetEnv(time_in,Li,temp_buff);
   printf("time=%d temp={%.2lf %.2lf}\n",time_in,temp_buff[0],temp_buff[4]);
   //after the ProcessAll function,use the functions below to get all kinds of information
   //DumpInfo(), which print out all the variables
   //GetLaserSwith(), which returns wheather the laser is on or off
   //GetHeight(), which return the position of laser in z axis 
   //GetElevation(), which returns the elevation angle when transmited
   //GetAzimuth(), which returns the azimuth angle when transmited

   //printf("index=%d time=%d hms=%d:%d:%d elevation=%.2lf azimuth=%.2lf\n",index0,time_in,hour,min,sec,RotateDB::GetHead()->GetElevation(),RotateDB::GetHead()->GetAzimuth());
   printf("res=%ld time=%d hms=%d:%d:%d elevation=%.2lf azimuth=%.2lf\n",res,time_in,hour,min,sec,RotateDB::GetHead()->GetElevation(),RotateDB::GetHead()->GetAzimuth());
}
