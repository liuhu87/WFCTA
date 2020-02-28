{
   gSystem->Load("lib/lib.so");
   StatusDB::jdebug=10;
   StatusDB* pb=StatusDB::GetHead();

   //double time0=20200115012905.;
   //int time1=CommonTools::Convert(time0);
   //long int Time2=pb->GetReadbackTime(6,time1);

   for(int ii=0;ii<10;ii++){
   for(int jj=0;jj<1;jj++){
   cout<<StatusDB::GetHead()->GetPreTemp(1,1575555663+ii*10,1,(char*)"root://eos01.ihep.ac.cn//eos/lhaaso/decode/wfcta/2019/1205/ES.31622.FULL.WFCTA01.20191205221946.018.dat.status.root")<<endl;
   }
   }

   //printf("time=%d Time=%ld\n",time1,Time2);
}
