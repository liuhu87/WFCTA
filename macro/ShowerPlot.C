{
   gStyle->SetPadGridX(1);
   gStyle->SetPadGridY(1);

   gSystem->Load("lib/lib.so");
   CorsikaIO::jdebug=0;
   ReadTrack::jdebug=0;
   ShowerPlot::jdebug=6;
   ReadTrack::DoPlot=true;
   ReadTrack::headbyte=16;

   //ReadTrack::particle=(1<<0);
   //ReadTrack::tlimit[0]=0.0;
   //ReadTrack::tlimit[1]=3.2e-4;
   ReadTrack::IniRange[0][0]=-4.8e6;
   ReadTrack::IniRange[0][1]=-3.e6;
   ReadTrack::IniRange[1][0]=-3e5;
   ReadTrack::IniRange[1][1]=3e5;
   ReadTrack::IniRange[2][0]=0;
   ReadTrack::IniRange[2][1]=0.5e7;

   ShowerPlot* gplot=new ShowerPlot("iron");
   //gplot->Add("/eos/user/h/hliu/Corsika/ShowerImageData/iron/1E14/DAT000002.track_em",0);
   //gplot->Add("/eos/user/h/hliu/Corsika/ShowerImageData/iron/1E14/DAT000002.track_mu",1);
   //gplot->Add("/eos/user/h/hliu/Corsika/ShowerImageData/iron/1E14/DAT000002.track_hd",2);
   //gplot->Add("/eos/user/h/hliu/Corsika/ShowerImageData/iron/1E14/CER000002",3);
   gplot->Add("/eos/user/h/hliu/Corsika/corsika-76900/run/DAT000002.track_em",0);
   gplot->Add("/eos/user/h/hliu/Corsika/corsika-76900/run/DAT000002.track_mu",1);
   gplot->Add("/eos/user/h/hliu/Corsika/corsika-76900/run/DAT000002.track_hd",2);
   gplot->Add("/eos/user/h/hliu/Corsika/corsika-76900/run/CER000002",3);
   gplot->Read();
   TCanvas* cc=gplot->Draw(0);
   cc->SetGridx(1);
   cc->SetGridy(1);
   if(cc) cc->SaveAs("test.png");
}
