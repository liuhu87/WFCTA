#include "TMath.h"
#include "analysis.h"
double WCamera::GetPMTMAP()
{
    double D_ConeOut = 25.4;
    double interval = 1.0;
    int k,Interx,Intery;
    for(int i=0; i<Pix; i++){
        Intery = i/4;
	for(int j=0; j<Pix; j++){
	    Interx = j/4;
	    k = i*Pix + j;
	    if(i%2==0)
		pmtmap[k][0] = (j+0.5)*D_ConeOut + interval*Interx;
	    if(i%2==1)
		pmtmap[k][0] = (j+1)*D_ConeOut + interval*Interx;
	    pmtmap[k][1] = (Pix-i)*D_ConeOut + interval*Intery;
	}
    }
}

double WCamera::GetPixX(int ipix)
{
  return pmtmap[ipix][0];
}

double WCamera::GetPixY(int ipix)
{
  return pmtmap[ipix][1];
}
