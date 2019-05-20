#include <iostream>
#include <fstream>
#include <string>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "creadparam.h"
#include "wtelescope.h"
#include "wcamera.h"
#include "telescopeparameters.h"
#include "cer_event.h"
#include "SquareCone.h"
#include "TNtuple.h"
using namespace std;
main(int argc, char *argv[])
{
  /*parameters read from the input card */
  int NSBFlag,WLFlag, ThinFlag, FilterFlag,MirrorPointErrorFlag, FadcFlag,MirrorGeometry;
  int CTNumber;
  float MirrorSizeX, MirrorSizeY, MirrorSpot, MirrorPointError;
  char inputfilename[PATH_MAX_LENGTH],outputfilename[PATH_MAX_LENGTH];


  int Fadc_bins, Fadc_length; //The number of FADC bins and the bin length 
  float nsb, triggersigma;                  //night sky background 
  float FixTriggerThreshold;  //the trigger threshold for SiPM

  /*The telescope arrays settings*/
  float *CT_X, *CT_Y, *CT_Z, *CT_Azi, *CT_Zen;

  float *timefirst, *timelast,*timemean, *ncphoton;

  /*cos directions of the pointing of telescopes */
  float *CT_m, *CT_n, *CT_l;

  /*cos directions of the primary directions of showers */
  float prim_m, prim_n, prim_l;

  /*space angle between the arriving directions of showers and the pointing of the telescopes */
  float *Space_angle;

  float zenith, azimuth, Xmax, Nmax;
  
  float *CoreX, *CoreY;
  double  mm, nn;

  // * ==== iuse: for the multiple use of CER event ==== *//
  int Nuse = 20;  ///< 20 is the max times of CER event use
  int iuse = 0;

  //*=== cone tracing and SiPM simulation ====*//
  int ict, icone, itube, icell,ix, iy, iternum, flag;
  double x, y, deltax, deltay,xx;
  double hitpos[3], dircos[3], Weight;

  //*Get the parameters from the input card*//
  WReadConfig *readconfig = new WReadConfig();  

  readconfig->readparam(argv[1]);
  WLFlag = readconfig->GetWLFlag();
  ThinFlag = readconfig->GetThinFlag();
  FilterFlag = readconfig->GetFilterFlag();

  MirrorGeometry = readconfig->GetMirrorGeometry();
  MirrorPointErrorFlag = readconfig->GetMirrorPointErrorFlag();
  MirrorPointError = readconfig->GetMirrorPointError();

  readconfig->GetMirrorSize(&MirrorSizeX , &MirrorSizeY);
  MirrorSpot =readconfig->GetMirrorSpot();

  CTNumber =readconfig-> GetCTNumber();
 
  CT_X = new float[CTNumber];
  CT_Y = new float[CTNumber];
  CT_Z = new float[CTNumber];
  CT_Azi = new float[CTNumber];
  CT_Zen = new float[CTNumber];
  CT_m = new float[CTNumber];
  CT_n = new float[CTNumber];
  CT_l = new float[CTNumber];
  CoreX = new float[CTNumber];   //added by linglingMa on 2018.03.05
  CoreY = new float[CTNumber];   //To store the core positions of each events
  Space_angle = new float[CTNumber];
  timefirst = new float[CTNumber];
  timelast = new float[CTNumber];
  timemean = new float[CTNumber];
  ncphoton = new float[CTNumber];

  for (int ict = 0; ict <CTNumber; ict++){
     CT_X[ict] =  readconfig->GetCTPosition(ict,0); 
     CT_Y[ict] =  readconfig->GetCTPosition(ict,1);
     CT_Z[ict] =  readconfig->GetCTPosition(ict,2);
     CT_Zen[ict] =  readconfig->GetCTPosition(ict,3) * TMath::DegToRad(); 
     CT_Azi[ict] =  readconfig->GetCTPosition(ict,4) * TMath::DegToRad();
     CT_m[ict] = sin(CT_Zen[ict]) * cos(CT_Azi[ict]);
     CT_n[ict] = sin(CT_Zen[ict]) * sin(CT_Azi[ict]);
     CT_l[ict] = cos(CT_Zen[ict]);
  }

  //*The total intensty of nsb in the trigger window equlas Fadc_bins X Fabs_length X nsb *//
  FadcFlag = readconfig->GetFadcFlag(); 
  Fadc_bins = readconfig->GetFadcBins();
  Fadc_length = readconfig->GetFadcLength();

  //** the option is added on 2017-10-30 by LinglingMa **//
  //if the Nsb is not simulated, the fix trigger threshold is used for SiPM Trigger *//
  NSBFlag = readconfig->GetNSBFlag();
  if(NSBFlag) {
     nsb = readconfig->GetNSB();
     triggersigma = readconfig->GetTriggerSigma();
  }
  else FixTriggerThreshold = readconfig->GetFixTriggerThreshold();
 
  strcpy( inputfilename, argv[2]);
  strcpy( outputfilename, argv[3]);

  WTelescope **telescope;
  telescope = new WTelescope * [CTNumber];
 
  //* Init the telescope array *//
  for(int ict=0; ict<CTNumber; ict++){
     telescope[ict] = new WTelescope();
     //telescope[ict]->SetMirrorPointErrorFlag(MirrorPointErrorFlag);
     //telescope[ict]->SetMirrorPointError(MirrorPointError);
     telescope[ict]->SetMirrorSpot(MirrorSpot);
     telescope[ict]->SetMirrorGeometry(MirrorGeometry);
     telescope[ict]->SetMirror();
     telescope[ict]->SetPointing(CT_Zen[ict],CT_Azi[ict]);
     telescope[ict]->SetMirrorPointError(MirrorPointErrorFlag,MirrorPointError);
     telescope[ict]->SetEulerMatrix(CT_Zen[ict],CT_Azi[ict]);
    
     //* 2018.1.10 updated by Lingling Ma  the transmissivity of the filter //
     //* is stored in a txt file Data/filter.txt                            //
     telescope[ict]->SetTransmissivity();
     telescope[ict]->SetReflectivity();
  }

  //* to get the primary information of the showers *//
  CER_Event *Event; 
  Event = new CER_Event();
  
  //* to get the information of the Cherenkov photons *//
  CER_bunch *ph_bunch;
  ph_bunch = new CER_bunch();
  
  //*define the camera*//
  WCamera *camera = new WCamera();
  camera->SetSiPMMAP();

  //** the option is added on 2017-10-30 by LinglingMa **//
  //if the Nsb is not simulated, the fix trigger threshold is used for SiPM Trigger *//
  if(NSBFlag) {
     camera->SetNSB(nsb*Fadc_length*Fadc_bins);
printf("The NSB level %f\n",nsb*Fadc_length*Fadc_bins);
     camera->SetTriggerSigma(triggersigma);
  }
  else camera->SetFixTriggerThreshold(FixTriggerThreshold);

  camera->SetCTNumber(CTNumber);
  camera->Init();

  SquareCone *cone = new SquareCone();

  //*the outputfile *//
  
  TFile *file = new TFile(outputfilename,"recreate");

  //TNtuple *nt = new TNtuple("nt","dfr","t");
  TTree *EventsTree = new TTree("events","RayTrace");

  EventsTree->Branch("iEvent",&Event->Event_Number_,"iEvent/I");
  EventsTree->Branch("iUse",&iuse,"iUse/I");
  EventsTree->Branch("id",&Event->Primary_id_,"id/D");
  EventsTree->Branch("energy",&Event->Primary_Energy_,"energy/D");
  EventsTree->Branch("zenith", &Event->Primary_zenith_,"zenith/D");
  EventsTree->Branch("azimuth", &Event->Primary_azimuth_, "azimuth/D");
  EventsTree->Branch("corex", &Event->Primary_core_x_, "corex/D");
  EventsTree->Branch("corey", &Event->Primary_core_y_, "corey/D");
  EventsTree->Branch("Xmax",&Event->Xmax_,"Xmax/D");
  EventsTree->Branch("Nmax",&Event->Nmax_,"Nmax/D");

  EventsTree->Branch("ntel", &CTNumber,"ntel/I");
  EventsTree->Branch("TelX", CT_X, Form("TelX[%d]/F",CTNumber) );
  EventsTree->Branch("TelY", CT_Y, Form("TelY[%d]/F",CTNumber)  );
  EventsTree->Branch("TelZ", CT_Zen,Form("TelZ[%d]/F",CTNumber) );
  EventsTree->Branch("TelA", CT_Azi,Form("TelA[%d]/F",CTNumber) );
  //EventsTree->Branch("CoreX",CoreX,Form("CoreX[%d]/F",CTNumber));
  //EventsTree->Branch("CoreY",CoreY,Form("CoreY[%d]/F",CTNumber));


  EventsTree->Branch("TubeSignalAfterConeTracing",&(camera->TubeSignalAfterConeTracing));
  EventsTree->Branch("TubeSignalIntoCone", &(camera->TubeSignalIntoCone));
  EventsTree->Branch("TubeSignal", &(camera->TubeSignal));
  EventsTree->Branch("TubeTrigger", &(camera->TubeTrigger));
  EventsTree->Branch("TelTrigger",&(camera->TelTrigger));
  EventsTree->Branch("TubeSignalInTriggerWindow",&(camera->TubeSignalInTriggerWindow));


  //* To read the corsika file *// 

  //* ==== Particle, subblock, record length set based on THIN flag ==== *//
  int particlelength;
  if(ThinFlag) 
     particlelength = PARTICLE_LENGTH_THINNING;
  else
     particlelength = PARTICLE_LENGTH_NO_THINNING;

  int subblocklength = particlelength * N_PARTICLE;
  int recordlength = subblocklength * N_SUBBLOCK;

  
  union REC {
      float f;
      char c[4];
  };
  REC *rec = new REC [recordlength];

  union {
     int i;
     char c[4];
  } padding;

  char sss[5];

  int iblock = 0;
  int jrune = 0;
  int jrunh = 0;

  ifstream fin;
  fin.open(inputfilename);
  if (!fin.is_open()) {
     cout << "Error when open CER file " << inputfilename << endl << endl;
     exit(1);
  }

  fin.seekg(0,ios::end);
  size_t filesize = fin.tellg(); 
  cout << "The CER file size is " << filesize << " bytes" << endl << endl;
  if (!filesize) {
     cout << "Warning: The file is empty!" << endl << endl;
     fin.close();
     return(0);
  }

  fin.seekg(0,ios::beg);

  /* loop to read the corsika file */
  while (fin.good()) {

   //There are padding words before and after a record in the fortran output
   fin.read(padding.c,4);
   if (!fin.good()) {
     if (jrunh!=1||jrune!=1||iblock==0) {
        cout << "Error in reading corsika file \"" << inputfilename << "\"" << endl;
        cout << "Is the file truncated?" << endl;
        if (iblock==0) cout << "Is the file empty?" << endl;
        if (jrunh!=2) cout << "There is no RUNH sub-block!" << endl;
        if (jrune!=1) cout << "There is no RUNE sub-block!" << endl;
        cout << endl;
        exit(1);
     }
     else {
        //A successful reading would end here
        cout << "Succeed in reading corsika file \"" << inputfilename << "\"" << endl << endl;
        break;
     }
  }

  int iflag = 0;
  if (padding.i!=4*recordlength) iflag = 1;
  for(int iLength = 0;iLength<recordlength;iLength++){
      fin.read(rec[iLength].c,4);
  }
  if (!fin.good()) iflag = 1;
      fin.read(padding.c,4);
  if (padding.i!=4*recordlength) iflag = 1;

  if (iflag) {
      cout << "Error in reading corsika file \"" << inputfilename << "\"" << endl;
      cout << "Is the file truncated?" << endl;
      exit(1);
  }
  iblock++;
  

  //* ====== CORSIKA CER file data sub-block read ====== *//
  for(int isubblock  = 1 ;  isubblock <= N_SUBBLOCK ; isubblock++) {
    
      int iptr = subblocklength*(isubblock-1) - 1;
      memcpy(sss,&rec[iptr+1].f,4);
      sss[4] = '\0';
      //* ================ RUN HEADER ================ *//
      if (strcmp(sss,"RUNH")==0){
         if (iuse++ == 0){
            jrunh++;
         }
      }
      //* ================ RUN END =================== *//
      else if (strcmp(sss,"RUNE")==0) {
           if(iuse == Nuse){
              jrune++;
           }
           else {
              fin.seekg(0,ios::beg);
              break;
          }
      }
     //* ================ EVENT HEADER ================ *//
      else if (strcmp(sss,"EVTH")==0) {

           Event->Init();
           Event->Event_Number_ = int(rec[iptr+2].f);
           Event->Event_use_Number_ = iuse;
         
           Event->Primary_Energy_ = rec[iptr+4].f;
           Event->Primary_id_ = rec[iptr+3].f;

           Event->Primary_zenith_ = rec[iptr+11].f;
           Event->Primary_azimuth_ = rec[iptr+12].f;

           Event->Primary_core_x_ = rec[iptr+98+iuse].f;
           Event->Primary_core_y_ = rec[iptr+118+iuse].f;
           Event->Ob_level_ = rec[iptr+47+1].f;
       printf("%f %f\n",Event->Primary_zenith_*57.3,Event->Primary_azimuth_*57.3);  
           prim_m = sin(rec[iptr+11].f) * cos( rec[iptr+12].f);
           prim_n = sin(rec[iptr+11].f) * sin( rec[iptr+12].f);
           prim_l = cos(rec[iptr+11].f);

           for(int ict=0; ict<CTNumber; ict++){
              Space_angle[ict] = prim_m * CT_m[ict] + prim_n * CT_n[ict] + prim_l * CT_l[ict];
              if(Space_angle[ict]>1)  Space_angle[ict] = 1;
              if(Space_angle[ict]<-1) Space_angle[ict] = -1;
              Space_angle[ict] = acos(Space_angle[ict])*TMath::RadToDeg();
              CoreX[ict] = Event->Primary_core_x_+CT_X[ict];
              CoreY[ict] = Event->Primary_core_y_+CT_Y[ict];
              telescope[ict]->SetXY(CoreX[ict],CoreY[ict]);
              timefirst[ict] = 1000000000;
              timelast[ict] = 0.;
              timemean[ict] = 0.;
              ncphoton[ict] = 0.;
              camera->ReSet(ict);
           } //for

            
       }  

       //* ================ EVENT END ================ *//
       else if (strcmp(sss,"EVTE")==0) {
           printf("The end of the event\n");
           Event->Nmax_ = rec[iptr+255+1].f;
           Event->Xmax_ = rec[iptr+255+3].f;  
           camera->PhotonCellToTube();

           //** Added by lingling Ma on 2018-1-22             **//
           //** to Get the PeakTime on each SiPM              **//
           camera->GetPeakTime();   
          
           //** Added by linglingMa on 2018-3-13              **//
           //** to Get the photon ratio in the Trigger window **//    
           camera->GetPhotonInTriggerWindow(Fadc_length*Fadc_bins);

           if(NSBFlag){   
             //camera->AddNSB(FadcFlag);
             camera->GetTubeTrigger(NSBFlag,FadcFlag);
           }
           else
             camera->GetTubeTrigger(NSBFlag,FadcFlag);          
        
           camera->GetTelescopeTrigger(CTNumber,CT_Zen, CT_Azi); 
           EventsTree->Fill();
       }

       //* ================ LONG sub-block ================ *//
       else if (strcmp(sss,"LONG")==0) {
       
       }
       //* ================ CERenkov sub-block ================ *//
       else {
          for (int iptcl=1; iptcl<=N_PARTICLE; iptcl++) {

            int iptrnow = (iptcl-1) * particlelength + iptr;

            ph_bunch->Init();

            if(ThinFlag)
                ph_bunch->weight_ = rec[iptrnow+8].f;
            else
                ph_bunch->weight_ = 1; 

            if(WLFlag){
                int wavelength_and_bunch = int(rec[iptrnow+1].f);
                //ph_bunch->wavelength_ = wavelength_and_bunch%1000;
               // ph_bunch->nclight_ = wavelength_and_bunch/1000;
                ph_bunch->wavelength_  = 1000*(rec[iptrnow+1].f-int(rec[iptrnow+1].f));
                ph_bunch->nclight_ = int(rec[iptrnow+1].f);
//printf("%f %f\n",ph_bunch->nclight_,ph_bunch->wavelength_ );
            }
            else
                ph_bunch->nclight_ = rec[iptrnow+1].f;
             
            if (ph_bunch->nclight_) {
                
               for(int ict=0; ict<CTNumber; ict++){
            
                  if(Space_angle[ict]>30) continue;
 
                  //*=== Init the ph_bunch for each telescope*===//
                  //*=== debuged by Lingling Ma 2016.7.20*===//
                  ph_bunch->SetCERBunch(rec[iptrnow+2].f,rec[iptrnow+3].f,rec[iptrnow+4].f,
                                       rec[iptrnow+5].f,rec[iptrnow+6].f,rec[iptrnow+7].f); //x,y,u,v,t,h
                 
                  if(telescope[ict]->IncidentTel(ph_bunch->x_,ph_bunch->y_,MirrorSizeX,MirrorSizeY)){

                     ph_bunch->IntoTelArea(telescope[ict]->Telx_,telescope[ict]->Tely_);
                     

                     //*=== coordinate transformation to telescope system*===//
                     //*=== debuged by LinglingMa 2016.720 *===//
                     telescope[ict]->Euler(ph_bunch);

                     for(int Nclight=0; Nclight<int(ph_bunch->nclight_*ph_bunch->weight_); Nclight++){
                           
                      //if(FilterFlag){
                      //     int Transimissivity = telescope[ict] ->GetTransmissivity( ph_bunch->wavelength_);
                      //     if(Transimissivity==0) continue;
                      //}
                      if(telescope[ict]->Reflected(ph_bunch->wavelength_)){
                           double temp = telescope[ict]->RayTrace(MirrorSpot,
                                        ph_bunch->x_*10, ph_bunch->y_*10, ph_bunch->z_*10,
                                        ph_bunch->u_, ph_bunch->v_, ph_bunch->l_,
                                        &ph_bunch->xc_, &ph_bunch->yc_,&ph_bunch->time_raytrace_, &ph_bunch->u1_, &ph_bunch->v1_);
                
                           if(temp < 0) continue; 
                           if(FilterFlag){
                                double filter_m = ph_bunch->u1_;
                                double filter_n = ph_bunch->v1_;
                                double filter_l = sqrt(1-filter_m*filter_m-filter_n*filter_n);
                                //double filter_angle = acos(filter_m*0+filter_n*0+filter_l*1)*TMath::RadToDeg();
                                double filter_angle = acos(filter_l)*TMath::RadToDeg();
                                int Transimissivity = telescope[ict] ->GetTransmissivity( ph_bunch->wavelength_,filter_angle);
                                if(Transimissivity==0) continue;
                           }
                     
                           //*Get the time information of the cphotons *//  
                           ph_bunch->time_raytrace_ += rec[iptrnow+6].f;

                           timemean[ict] += ph_bunch->time_raytrace_;

                           if(ph_bunch->time_raytrace_>timelast[ict]) 
                              timelast[ict] = ph_bunch->time_raytrace_;

                           if(ph_bunch->time_raytrace_<timefirst[ict]) 
                              timefirst[ict] = ph_bunch->time_raytrace_;
                      
                           //*The photons that can enter the wenston cone*// 
                           icone = camera->GetCone( ph_bunch->xc_, ph_bunch->yc_);
                           if(icone<0) continue;
                           camera->PhotonIntoCone(ict,icone,1);

                           //*To Get the cone coordinates*//
                           x = camera->GetSiPMX(icone);
                           y = camera->GetSiPMY(icone); 

                           deltax =  ph_bunch->xc_ - x;
                           deltay =  ph_bunch->yc_ - y;
                           dircos[0] = -ph_bunch->u1_;
                           dircos[1] = ph_bunch->v1_;
                           dircos[2] = -sqrt(1-ph_bunch->u1_*ph_bunch->u1_-ph_bunch->v1_*ph_bunch->v1_);

                           //*The ray trace in the cone *//
                           cone -> SetInitDir(dircos);
                           cone -> SetInitPos(deltax,deltay);
                           cone -> SquareRaytracing();
                           cone -> GetHitPos(hitpos,Weight);
                           if(!cone -> GetStrike()) continue;
                           //ConeTracing(deltax, deltay,  x,  y,  dircos, hitpos, &iternum, &flag);
                           hitpos[0] += x;
                           hitpos[1] += y;
                           itube = camera->GetTube(hitpos[0],hitpos[1]);
                           if(itube<0) continue;           //Throw away the bad point after ConeTracing.
                           if(itube!=icone) continue;
                           xx = gRandom->Rndm();
                           if(xx>Weight) continue;
                           camera->PhotonAfterConeTracing(ict,icone,1);

                           hitpos[0] = hitpos[0] - x;
                           hitpos[1] = hitpos[1] - y;
                           ix = int((hitpos[0]+D_SiPM/2)/D_Cell);
                           iy = int((hitpos[1]+D_SiPM/2)/D_Cell);
                           icell = iy*NCell+ix;

                           // ** Added By linglingMa to count the arrive time of photons **//

                           camera->GetArriveTime(ict,itube,icell,int(ph_bunch->time_raytrace_),1);
                           camera->PhotonIntoCell(ict,itube,icell,1);
//printf("%d %d %d\n",ict,itube,icell);
 
                      }
                    } //loop for Nclight
                 } //if telescope
               
               } //for ict     
             }//if Nclight
          }// for Nparticle
       }// else cherenkov block
     }//sub block
   }  //while
   file->Write();
   file->Close();
}
