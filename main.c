
#include <stdlib.h>
#include <stdio.h>
#include "fft_c.h"


/* Simple Test - 16 Bins */
#define NBINS 16
#define NLEV 4

int main()
{
	float *ir,*ii,*or,*oi; //,k;
	int i=0;
	
	ir = calloc(4*NBINS,sizeof(float));
	ii = ir + NBINS;
	or = ii + NBINS;
	oi = or + NBINS;
	
	for(i = 0; i < NBINS ; i++)
		ir[i] = i;
	
	fft_c(ir,ii,or,oi,NLEV);
	printf("\n---------------------------------------------------------------------\n");
	printf("                  Test: FFT Radix-2 Decimation in Time\n\
                        v0.1	Emant 2014\n");
	printf("---------------------------------------------------------------------\n");
	
	//k = 1.0 / NBINS;
	for(i=0; i< NBINS; i++){
		printf("in[%2i] : % 9.4f %+9.4f*j   ,   out[%2i] : % 9.4f %+9.4f*j\n",
			i,ir[i],ii[i],i,or[i],oi[i]);
	}
	printf("---------------------------------------------------------------------\n");
	printf("----              IFFT(x+y*j) = FFT(1/N*(x-y*j))                 ----\n");
	printf("---------------------------------------------------------------------\n");
	ifft_c(or,oi,ir,ii,NLEV);
		for(i=0; i< NBINS; i++){
		printf("in[%2i] : % 9.4f %+9.4f*j   ,   out[%2i] : % 9.4f %+9.4f*j\n",
			i,or[i],oi[i],i,ir[i],ii[i]);
	}
	printf("---------------------------------------------------------------------\n\n");
	return 0;
}
