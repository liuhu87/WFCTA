#include "Cloud.h"
#include "WFTelescope.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TBrowser.h"
#include "stdlib.h"
#include "stdio.h"
int main(int argc,char** argv){
    //gSystem->Load("/afs/ihep.ac.cn/users/h/hliu/Documents/LHAASO/WFCTA/lib/lib.so");

    if(argc<2){
       printf("   Usage %s <in_name> <out_name>\n",argv[0]);
       return 0;
    }
    char* out_name=(argc>=3)?argv[2]:0;
    TFile* fout=(out_name)?TFile::Open(out_name,"RECREATE"):0;

    WFTelescopeArray::jdebug=0;
    WFTelescopeArray::DoSim=true;
    WFTelescopeArray::GetHead(Form("/afs/ihep.ac.cn/users/h/hliu/Documents/LHAASO/WFCTA/default.inp"));

    Cloud::SetBins();
    Cloud showcloud;
    showcloud.Reset();
    showcloud.ReadCloudMap(argv[1]);

    //new TBrowser;
    TCanvas* cc=new TCanvas();
    showcloud.Draw(WFTelescopeArray::GetHead());
    //showcloud.Draw(0);
    //if(gPad) gPad->Update();
    //if(gPad) gPad->WaitPrimitive();

    fout->cd();
    cc->Write("cloudmap");
    fout->Close();

    //delete WFTelescopeArray::GetHead();
    //showcloud.Clear();
    return 0;
}
