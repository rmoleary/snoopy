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
     printf("Please wait %s...\n", s[i < 0 ? -i : i]);
}


int main() {
	printf("General purpose Spectral Hydro solver v1.0\n");
	printf("(c) 2008 G. Lesur\n");
	printf("Using %dx%d grid\n",NX,NY);
	printf("Initializing...\n");
	init_common();
	init_flow();
	init_output();
	
	printf("Calling mainloop... touch status, output or stop to print more information.\n");
	please_wait();
	mainloop();
	
	finish_common();
	
	printf("Terminated.\n");
	
	return(0);
}
	