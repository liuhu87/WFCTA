{
  int wl;
  double p1,p2,p3,p4,p5,p6,p7,p8;
  FILE *fp;
  fp = fopen("filter_transparency.txt","r");
  TNtupleD *nt = new TNtupleD("nt","dfr","wl:p1:p2:p3:p4:p5:p6:p7:p8");
  while(!feof(fp)){
     fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %lf %lf\n",&wl,&p1,&p2,&p3,&p4,&p5,&p6,&p7,&p8);
     nt->Fill(wl,p1,p2,p3,p4,p5,p6,p7,p8);
  }
  fclose(fp);

  fp = fopen("filter.txt","r");
  TNtupleD *nt0 = new TNtupleD("nt0","dfr","wl:p1");
  while(!feof(fp)){
     fscanf(fp,"%d %lf \n",&wl,&p1);
     nt0->Fill(wl,p1);
  }
  fclose(fp);
}
