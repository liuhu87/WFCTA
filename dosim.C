#include "CorsikaChain.h"
#include "CorsikaEvent.h"
#include "EventNtuple.h"
#include "WFTelescope.h"
#include <stdlib.h>
#include <iostream>
using namespace std;
int main(int argc,char* argv[]){
    // arguments
    if (argc<5) {
       printf("   Usage %s <runlist> <out_name> <first> <last> <type> <run> <event>\n",argv[0]);
       return 0;
       //exit(0);
    }
    // vars 
    const char* runlist = argv[1];
    const char* outname = argv[2];
    int first = atoi(argv[3]);
    int last = atoi(argv[4]);
    int type = atoi(argv[5]);
    int run  = argc>6?atoi(argv[6]):-1;
    int event = argc>7?atoi(argv[7]):-1;
    cout<< runlist << "  " << outname << "  " << first << "  " << last << "  " << type << "  " << run << " " << event <<endl;

    CommonTools::IsLaser=false;
    CommonTools::InitHArrival();

    CorsikaIO::ReadRun=run;
    CorsikaIO::ReadEvent=event;
    CorsikaIO::jdebug=0;
    WFTelescopeArray::jdebug=0;
    WFTelescopeArray::DoSim=true;
    CorsikaChain* pchain=new CorsikaChain(Form("/afs/ihep.ac.cn/users/h/hliu/Documents/LHAASO/WFCTA/%s",runlist),first,last,type);
    //CorsikaChain* pchain=new CorsikaChain();
    //pchain->Add((char*)runlist);
    //pchain->ftype=type;
    WFTelescopeArray::GetHead(Form("/afs/ihep.ac.cn/users/h/hliu/Documents/LHAASO/WFCTA/default.inp"));
    CorsikaEvent* pevt=new CorsikaEvent();
    EventNtuple::GetHead(Form("%s",outname),0xF);
    //EventNtuple::GetHead(Form("/eos/user/h/hliu/%s",outname),3);
    pchain->ReadAllFile(-1,-1,pevt);
    //pchain->ReadAllRec(first,last,pevt);

    delete pchain;
    delete pevt;
    delete EventNtuple::GetHead();
    delete WFTelescopeArray::GetHead();
}
