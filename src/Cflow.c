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

void print_logo(void) {
	MPI_Printf("                 .88888888:. \n");
	MPI_Printf("                88888888.88888. \n");
	MPI_Printf("              .8888888888888888. \n");
	MPI_Printf("              888888888888888888 \n");
	MPI_Printf("              88' _`88'_  `88888 \n");
	MPI_Printf("              88 88 88 88  88888 \n");
	MPI_Printf("              88_88_::_88_:88888 \n");
	MPI_Printf("              88:::,::,:::::8888 \n");
	MPI_Printf("              88`:::::::::'`8888 \n");
	MPI_Printf("             .88  `::::'    8:88. \n");
	MPI_Printf("            8888            `8:888. \n");
	MPI_Printf("          .8888'             `888888. \n");
	MPI_Printf("         .8888:..  .::.  ...:'8888888:. \n");
	MPI_Printf("        .8888.'     :'     `'::`88:88888 \n");
	MPI_Printf("       .8888        '         `.888:8888. \n");
	MPI_Printf("      888:8         .           888:88888 \n");
	MPI_Printf("    .888:88        .:           888:88888: \n");
	MPI_Printf("    8888888.       ::           88:888888 \n");
	MPI_Printf("    `.::.888.      ::          .88888888 \n");
	MPI_Printf("   .::::::.888.    ::         :::`8888'.:. \n");
	MPI_Printf("  ::::::::::.888   '         .:::::::::::: \n");
	MPI_Printf("  ::::::::::::.8    '      .:8::::::::::::. \n");
	MPI_Printf(" .::::::::::::::.        .:888::::::::::::: \n");
	MPI_Printf(" :::::::::::::::88:.__..:88888:::::::::::' \n");
	MPI_Printf("  `'.:::::::::::88888888888.88:::::::::' \n");
	MPI_Printf("        `':::_:' -- '' -'-' `':_::::'` \n");

	return;
}

int main(int argc, char *argv[]) {
#ifdef MPI_SUPPORT
	int size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	
	// Some consistancy check
	if(size != NPROC) ERROR_HANDLER( ERROR_CRITICAL, "Wrong number of process. Check NPROC.");
	if(NX/NPROC < 1) ERROR_HANDLER( ERROR_CRITICAL, "NX should be a multiple of NPROC.");
	if(NY/NPROC < 1) ERROR_HANDLER( ERROR_CRITICAL, "NY should be a multiple of NPROC.");
#else
	rank=0;
#endif
	print_logo();
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
	