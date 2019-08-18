{
   gSystem->Load("lib/lib.so");
   ShowerPlot::jdebug=6;
   ReadTrack::jdebug=0;
   ReadTrack::DoPlot=true;
   ReadTrack::headbyte=16;

   //if(argv<4){
   //   printf("Use %s <dir name> <type name> <run number>\n",argc[0]);
   //   return 0;
   //}
   //char* dir=argc[1];
   //char* type=argc[2];
   //int run=atoi(argc[3]);

   const char* dir="/home/sunqn/Downloads/ShowerImageData/iron/1E14";
   const char* type="iron";
   int run=2;
   bool deterrange=false;
   int ipng=-1;

   ShowerPlot* gplot=new ShowerPlot(type);
   gplot->Add(Form("%s/DAT%06d.track_em",dir,run),0);
   gplot->Add(Form("%s/DAT%06d.track_mu",dir,run),1);
   gplot->Add(Form("%s/DAT%06d.track_hd",dir,run),2);
   //gplot->Add(Form("%s/CER%06d",dir,run),3);
   //printf("Adding %s/DAT%06d.track_em\n",dir,run);

   //define the range
   double coomax,zcoomax;
   double tminmax[2];
   if(deterrange){
      gplot->Read();
      TCanvas* c0=gplot->Draw(0);
      if(c0) c0->SetName("temp_canvas");
      coomax=0; zcoomax=gplot->plotrange[2][1];
      tminmax[0]=gplot->plotrange[3][0]*0.9; tminmax[1]=gplot->plotrange[3][1]*1.1;
      for(int ii=0;ii<2;ii++){
         if(fabs(gplot->plotrange[ii][0])>coomax) coomax=fabs(gplot->plotrange[ii][0]);
         if(fabs(gplot->plotrange[ii][1])>coomax) coomax=fabs(gplot->plotrange[ii][1]);
      }
      coomax*=1.05;
      zcoomax*=1.05;
      printf("XY Max=%f, Z Max=%f, Time Range={%f,%f}\n",coomax,zcoomax,tminmax[0],tminmax[1]);
      return 1;
   }

   ////electron
   //coomax=514724.076563;
   //zcoomax=11847064.950000;
   //tminmax[0]=0;
   //tminmax[1]=0.000543;

   ////gamma
   //coomax=5.1e5;
   //zcoomax=11847064.950000;
   //tminmax[0]=0;
   //tminmax[1]=0.000543;

   ////proton
   //coomax=514724.076563;
   //zcoomax=11847064.950000;
   //tminmax[0]=0;
   //tminmax[1]=0.000543;

   //iron
   coomax=1907663.362500;
   zcoomax=11847066.;
   tminmax[0]=0;
   //tminmax[1]=0.000434;

   tminmax[1]=1;

   ReadTrack::IniRange[0][0]=-coomax;
   ReadTrack::IniRange[0][1]=coomax;
   ReadTrack::IniRange[1][0]=-coomax;
   ReadTrack::IniRange[1][1]=coomax;
   ReadTrack::IniRange[2][0]=-1.e4;
   ReadTrack::IniRange[2][1]=zcoomax;

   //Draw png plot
   if(ipng<0){
   printf("Begin to draw png plot\n");
   gplot->Read();
   TCanvas* cc=gplot->Draw(0);
   if(cc) cc->SetName("png_canvas");
   if(cc) cc->SaveAs(Form("%s/%s.png",dir,type));
   return 1;
   }

   //Draw gif plot
   printf("Begin to draw gif plot\n");
   const int nt=10;
   ReadTrack::tlimit[0]=-1;
   for(int ii=0;ii<nt;ii++){
      if(ipng>=0&&ipng!=ii) continue;
      double tmax=tminmax[0]+(tminmax[1]-tminmax[0])/nt*(ii+0.5);
      printf("Draw %d tt=%lf\n",ii,tmax);
      ReadTrack::tlimit[1]=tmax;
      gplot->Read();
      TCanvas* cc=gplot->Draw(0);
      TView* view=cc->GetView();
      if(cc){
        cc->SetName(Form("Air_Shower_%d",ii));
        cc->SaveAs(Form("%s/gif/%d.png",dir,ii));
      }
   }

   return 1;
}
