{
   gSystem->Load("/afs/ihep.ac.cn/users/h/hliu/Documents/LHAASO/WFCTA/lib/lib.so");
   RotateDB::jdebug=0;
   RotateDB::UseGPSTime=false;

   int Li=2;

   double time1=20191203145200; //from 20XX-XX-XX XX:XX:XX
   double time2=20191203145300; //to   20XX-XX-XX XX:XX:XX
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

      //first four are temperature, the last one is the humidity
      double temp_buff[5];
      bool res=RotateDB::GetHead()->GetEnv(itime,Li,temp_buff);
      if(!res) continue;

      gr->SetPoint(gr->GetN(),itime*1.,temp_buff[0]);
      if((gr->GetN()%100)==1) printf("np=%d Time=%d(%04d-%02d-%02d %02d:%02d:%02d) Temp={%.2lf,%.2lf,%.2lf,%.2lf} Humi=%.2lf\n",gr->GetN(),itime,year,month,day,hour,min,sec,temp_buff[0],temp_buff[1],temp_buff[2],temp_buff[3],temp_buff[4]);
   }
   gr->GetXaxis()->SetTimeDisplay(1);
   gr->GetXaxis()->SetNdivisions(-203);
   gr->GetXaxis()->SetTimeFormat("%Ss/%Mm/%Hh/%d/%m%F1970-01-01 00:00:00s0");
   gr->SetName(Form("L%d_Temp",Li));
   gr->Draw("apl");
   //gr->SaveAs(Form("L%d_Temp.root",Li));
}
