{
   gSystem->Load("lib/lib.so");

   ShowerPlot::jdebug=1;

   ShowerPlot* gplot=new ShowerPlot("Proton");
   gplot->Add("/home/huliu/Downloads/software/corsika-76900/run/DAT000002.track_em",0);
   gplot->Add("/home/huliu/Downloads/software/corsika-76900/run/DAT000002.track_mu",1);
   gplot->Add("/home/huliu/Downloads/software/corsika-76900/run/DAT000002.track_hd",2);
   gplot->Add("/home/huliu/Downloads/software/corsika-76900/run/CER000002",3);
   const int nt=10;
   ReadTrack::tlimit[0]=-1;
   ReadTrack::IniRange[0][0]=-1.5e6;
   ReadTrack::IniRange[0][1]=1.5e6;
   ReadTrack::IniRange[1][0]=-1.5e6;
   ReadTrack::IniRange[1][1]=1.5e6;
   ReadTrack::IniRange[2][0]=0;
   ReadTrack::IniRange[2][1]=1.2e7;
   for(int ii=0;ii<nt;ii++){
   ReadTrack::jdebug=ii==0?0:0;
      double tmax=0+(5.0e-4-0.)/nt*(ii+0.5);
      printf("Draw %d tt=%lf\n",ii,tmax);
      ReadTrack::tlimit[1]=tmax;
      gplot->Read();
      TCanvas* cc=gplot->Draw(3);
      TView* view=cc->GetView();
      if(cc){
        cc->SetName(Form("Air_Shower_%d",ii));
        cc->SaveAs(Form("gif/%d.png",ii));
      }
   }
}
