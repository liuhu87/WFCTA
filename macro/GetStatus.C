{
   gSystem->Load("lib/lib.so");
   StatusDB::jdebug=0;
   StatusDB* pb=StatusDB::GetHead();

   double time0=20200115012905.;
   int time1=CommonTools::Convert(time0);
   long int Time2=pb->GetReadbackTime(6,time1);

   printf("time=%d Time=%ld\n",time1,Time2);
}
