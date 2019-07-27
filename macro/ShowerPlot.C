{
   gSystem->Load("lib/lib.so");
   ReadTrack::jdebug=0;

   ShowerPlot* gplot=new ShowerPlot("Proton");
   gplot->Add("/home/huliu/Downloads/software/corsika-76900/run/DAT000002.track_em",0);
   gplot->Add("/home/huliu/Downloads/software/corsika-76900/run/DAT000002.track_mu",1);
   gplot->Add("/home/huliu/Downloads/software/corsika-76900/run/DAT000002.track_hd",2);
   gplot->Add("/home/huliu/Downloads/software/corsika-76900/run/CER000002",3);
   gplot->Read();
   gplot->Draw();
}
