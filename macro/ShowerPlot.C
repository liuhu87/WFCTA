{
   gSystem->Load("lib/lib.so");
   ReadTrack::jdebug=0;
   ShowerPlot::jdebug=6;
   ReadTrack::DoPlot=true;
   ReadTrack::headbyte=16;

   //ReadTrack::particle=(1<<0);
   //ReadTrack::tlimit[0]=0.0;
   //ReadTrack::tlimit[1]=1.0;
   //ReadTrack::IniRange[0][0]=-7.6e5;
   //ReadTrack::IniRange[0][1]=7.6e5;
   //ReadTrack::IniRange[1][0]=-7.6e5;
   //ReadTrack::IniRange[1][1]=7.6e5;
   //ReadTrack::IniRange[2][0]=0;
   //ReadTrack::IniRange[2][1]=5.4e6;

   ShowerPlot* gplot=new ShowerPlot("electron");
   gplot->Add("/home/huliu/Downloads/ShowerImageData/electron/1E14/DAT000002.track_em",0);
   gplot->Add("/home/huliu/Downloads/ShowerImageData/electron/1E14/DAT000002.track_mu",1);
   gplot->Add("/home/huliu/Downloads/ShowerImageData/electron/1E14/DAT000002.track_hd",2);
   gplot->Add("/home/huliu/Downloads/ShowerImageData/electron/1E14/CER000002",3);
   gplot->Read();
   TCanvas* cc=gplot->Draw(0);
   if(cc) cc->SaveAs("test.png");
}
