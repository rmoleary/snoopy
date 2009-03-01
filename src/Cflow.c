#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <fftw3.h>

#include "common.h"
#include "mainloop.h"
#include "output.h"
#include "initflow.h"
#include "gfft.h"

void please_wait(void)
{
     int i;
     const char *s[] =
     {
       "(while a large software vendor in Seattle takes over the world)",
       "(and remember, this is faster than Java)",
       "(and dream of faster computers)",
       "(checking the gravitational constant in your locale)",
       "(at least you are not on hold)",
       "(while X11 grows by another kilobyte)",
       "(while Windows Vista reboots)",
       "(correcting for the phase of the moon)",
       "(your call is important to us)",
       "(while the Linux user-base doubles)",
       "(while you decide where you want to go tomorrow)",
       "(exorcising evil spirits)",
       "(while the C++ standard gains another page)",
     };
     int choices = sizeof(s) / sizeof(*s);
	 srand (time (NULL));
     i = rand() % choices;
     MPI_Printf("Please wait %s...\n", s[i < 0 ? -i : i]);
}


int main(int argc, char *argv[]) {
#ifdef MPI_SUPPORT
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#else
	rank=0;
#endif
	MPI_Printf("General purpose Spectral Hydro solver v1.0\n");
	MPI_Printf("(c) 2004-2009 G. Lesur\n");
	MPI_Printf("Using %dx%dx%d grid\n",NX,NY,NZ);
	MPI_Printf("Initializing...\n");
	init_common();
	init_gfft();
	init_flow();
	init_output();
	
	MPI_Printf("Calling mainloop... touch status, output or stop to print more information.\n");
	please_wait();
	mainloop();
	
	finish_common();
	
	MPI_Printf("Terminated.\n");
	
	return(0);
}
	