{
   gSystem->Load("lib/lib.so");
   Cloud::DoTempCorr=false;
   Cloud::LoadTelSetting((char*)"/afs/ihep.ac.cn/users/h/hliu/public/WFDataDir/default.inp");
   CloudDB::jdebug=0;

   int iTel=1;
   CloudDB* pc=CloudDB::GetHead();


   //double time0=20191231195500; //year month day hour minute second
   ////int Time=1570000000;
   //int Time=CommonTools::Convert(time0);
   //double IBTemp=pc->GetTelAveIBTemp(Time,iTel);
   //double temp0=pc->GetTemperature(Time,0);
   //double temp=pc->GetTemperature(Time,1);
   //double humi=pc->GetHumidity(Time);
   //printf("Time=%d IBTemp=%.2lf temp={%.2lf,%.2lf} humi=%.2lf\n",Time,IBTemp,temp0,temp,humi);

   double time1=20191227200000; //from 2019-12-27 20:00:00
   double time2=20200110080000; //to   2020-01-10 08:00:00
   int Time1=CommonTools::Convert(time1);
   int Time2=CommonTools::Convert(time2);
   TGraph* gr=new TGraph();
   for(int itime=Time1;itime<=Time2;itime+=300){
      int year=(CommonTools::TimeFlag(itime,1)%100)+2000;
      int month=CommonTools::TimeFlag(itime,2);
      int day=CommonTools::TimeFlag(itime,3);
      int hour=CommonTools::TimeFlag(itime,4);
      int min=CommonTools::TimeFlag(itime,5);
      int sec=CommonTools::TimeFlag(itime,6);

      double IBTemp=pc->GetTelAveIBTemp(itime,iTel);
      double temp0=pc->GetTemperature(itime,0);
      double temp=pc->GetTemperature(itime,1);
      double humi=pc->GetHumidity(itime);
      if(IBTemp>500) continue;

      gr->SetPoint(gr->GetN(),itime*1.,IBTemp);
      if((gr->GetN()%100)==1) printf("np=%d Time=%d(%04d-%02d-%02d %02d:%02d:%02d) IBTemp=%.2lf Temp={%.2lf,%.2lf} Humi=%.2lf\n",gr->GetN(),itime,year,month,day,hour,min,sec,IBTemp,temp0,temp,humi);
   }
   gr->GetXaxis()->SetTimeDisplay(1);
   gr->GetXaxis()->SetNdivisions(-203);
   gr->GetXaxis()->SetTimeFormat("%Ss/%Mm/%Hh/%d/%m%F1970-01-01 00:00:00s0");
   gr->SetName("IBTemp");
   gr->SaveAs("IBTemp.root");
}
