#include "Cloud.h"
#include "WFTelescope.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TBrowser.h"
#include "stdlib.h"
#include "stdio.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TColor.h"
#include "TPaveText.h"

TStyle* myStyle = 0;
Int_t MyPalette7[100];
extern TStyle* define_mystyle();
extern void define_palettes();
int main(int argc,char** argv){
    //gSystem->Load("/afs/ihep.ac.cn/users/h/hliu/Documents/LHAASO/WFCTA/lib/lib.so");

    //define_mystyle();

    define_palettes();

    //gROOT->SetStyle("mystyle");

    //gROOT->ForceStyle(kTRUE);

    gStyle->SetPalette(100,MyPalette7);

    gStyle->SetNumberContours(50);

    gStyle->SetOptStat(0);
    //gStyle->SetFrameBorderSize(0.5);
    //gStyle->SetTitleAlign(33);
    //gStyle->SetTitleX(0.8);
    //gStyle->SetTitleY(0.99);
    //gStyle->SetTitleFont(52,"t");
    //gStyle->SetTitleSize(0.2,"t");

    if(argc<2){
       printf("   Usage %s <in_name> <out_name>\n",argv[0]);
       return 0;
    }
    char* out_name=(argc>=3)?argv[2]:0;
    if(!out_name) return 0;

    int length=strlen(out_name);
    int iloc=-1;
    for(int ii=length-1;ii>=0;ii--){
       if(out_name[ii]=='.'){ iloc=ii; break;}
    }
    if(iloc<0) {printf("Outname=%s error.",out_name); return -1;}
    char filename[100];
    char OutType[20];
    for(int ii=0;ii<length;ii++){
       if(ii<iloc) filename[ii]=out_name[ii];
       else if(ii==iloc) filename[ii]='\0';
       else OutType[ii-iloc-1]=out_name[ii];
    }
    OutType[length-1-iloc]='\0';

    WFTelescopeArray::jdebug=0;
    WFTelescopeArray::DoSim=true;
    WFTelescopeArray::GetHead(Form("default.inp"));
    //WFTelescopeArray::GetHead(Form("/afs/ihep.ac.cn/users/h/hliu/public/WFDataDir/default.inp"));

    Cloud::SetBins();
    Cloud showcloud;
    showcloud.Reset();
    showcloud.ReadCloudMap(argv[1]);

    //new TBrowser;
    TCanvas* cc=new TCanvas("canvas","",10,10,1000,1000);
    cc->SetBorderSize(0.4);

    TFile* fout=0;
    if(OutType=="root") fout=TFile::Open(out_name,"RECREATE");

    showcloud.Draw(WFTelescopeArray::GetHead());
    //gPad->SetRightMargin(0.1);
    //gPad->SetFixedAspectRatio();
    //showcloud.Draw(0);
    //if(gPad) gPad->Update();
    //if(gPad) gPad->WaitPrimitive();

    if(fout){
       fout->cd();
       cc->Write("cloudmap");
       fout->Close();
    }
    else{
       cc->SaveAs(out_name);
    }

    //delete WFTelescopeArray::GetHead();
    //showcloud.Clear();
    return 0;
}

TStyle* define_mystyle() {

  if (myStyle!=0) { delete myStyle; myStyle = 0; }

  if (myStyle==0) myStyle = new TStyle("mystyle","my style");

  cout << "define my style ...      to use it gROOT->SetStyle(\"mystyle\");" << endl;

  // Colors

  myStyle->SetFillColor(kGray-1);

  myStyle->SetFillStyle(1001);

  // Canvas & Pad (Title & Stat & Legend)

  myStyle->SetCanvasColor(0);

  myStyle->SetCanvasBorderMode(0);

  myStyle->SetCanvasBorderSize(1);

  myStyle->SetPadColor(0);

  myStyle->SetPadBorderMode(0);

  myStyle->SetPadBorderSize(1);

  myStyle->SetLegendBorderSize(1);

  myStyle->SetTitleBorderSize(1);

  myStyle->SetTitleFillColor(0);

  myStyle->SetStatBorderSize(1);

  myStyle->SetStatColor(0);

  myStyle->SetTitleBorderSize(1);

  myStyle->SetFrameFillColor(0);

  myStyle->SetFrameFillStyle(0);

  // Font

  Int_t font = 52;                   // Italic-Helvetica

  myStyle->SetTitleFont(font," ");   // title

  myStyle->SetTitleFont(font,"xyz"); // axis title

  myStyle->SetLabelFont(font,"xyz"); // axis labels

  myStyle->SetStatFont(font);        // stat

  myStyle->SetTitleFontSize(0.05);

  myStyle->SetLabelSize(0.04,"xyz");

  myStyle->SetStatFontSize(0.04);

  // Position and dimension

  myStyle->SetTitleOffset(0.95,"xyz");

  myStyle->SetTitleSize(0.05,"xyz");

  myStyle->SetTitleSize(0.05," ");

  // Stat option

  myStyle->SetOptStat("ourmen");

  // no ticks

  myStyle->SetPadTickX(0);

  myStyle->SetPadTickY(0);

  return myStyle;

}

void define_palettes() {

  // Another, more classic style, rainbow

  cout << "define classic rainbow (default) ...                       to use it gStyle->SetPalette(100,MyPalette7)" << endl;

  Double_t cr_red[10]   = { 67, 43, 99,181,237,253,252,247,207,138};

  Double_t cr_green[10] = { 64,128,189,225,248,237,179,103, 47,  0};

  Double_t cr_blue[10]  = {149,169,148,138,146,140, 88, 57, 62, 53};

  Double_t cr_stops[10];

  for (int i=0; i<10; i++) {

    cr_red[i] /= 255.;

    cr_green[i] /= 255.;

    cr_blue[i] /= 255.;

    cr_stops[i] = (1./9)*i;

  }

  Int_t nrF = TColor::CreateGradientColorTable(10,cr_stops,cr_red,cr_green,cr_blue,100);

  for (int i=0;i<100;i++) MyPalette7[i] = nrF+i;
}
