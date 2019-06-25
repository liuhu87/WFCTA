#include "Cloud.h"
#include "WFTelescope.h"
#include "TCanvas.h"
#include "TBrowser.h"
#include "stdlib.h"
#include "stdio.h"
int main(int argc,char** argv){
    //gSystem->Load("/afs/ihep.ac.cn/users/h/hliu/Documents/LHAASO/WFCTA/lib/lib.so");

    if(argc<2){
       printf("   Usage %s <in_name>\n",argv[0]);
       return 0;
    }

    WFTelescopeArray::jdebug=0;
    WFTelescopeArray::DoSim=true;
    WFTelescopeArray::GetHead(Form("/afs/ihep.ac.cn/users/h/hliu/Documents/LHAASO/default.inp"));

    Cloud::SetBins();
    Cloud showcloud;
    showcloud.Reset();
    showcloud.ReadCloudMap(argv[1]);

    //new TBrowser;
    TCanvas* cc=new TCanvas();
    showcloud.Draw(WFTelescopeArray::GetHead());
    //showcloud.Draw(0);
    if(gPad) gPad->Update();
    if(gPad) gPad->WaitPrimitive();
    //system("pause");
    //getchar();

    //delete WFTelescopeArray::GetHead();
    //showcloud.Clear();
    return 0;
}
