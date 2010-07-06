/*	
	This file is part of the Snoopy code.

    Snoopy code is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Snoopy code is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Snoopy code.  If not, see <http://www.gnu.org/licenses/>.
*/


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
#include "read_config.h"
#include "particles.h"

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
	MPI_Printf("\n");
	MPI_Printf("          .o.				        \n");	
	MPI_Printf("          |  |    _   ,		        \n");
	MPI_Printf("        .',  L.-'` `\\ ||	        \n");
	MPI_Printf("      __\\___,|__--,__`_|__         \n");
	MPI_Printf("     |    %%     `=`       |        \n");
	MPI_Printf("     | ___%%_______________|        \n");
	MPI_Printf("     |    `               |         \n");
	MPI_Printf("     | -------------------|         \n");
	MPI_Printf("     |____________________|         \n");
	MPI_Printf("       |~~~~~~~~~~~~~~~~|           \n");
	MPI_Printf("       | ---------------|  ,        \n");
	MPI_Printf("   \\|  | _______________| / /      \n");
	MPI_Printf("\\. \\,\\\\|, .   .   /,  / |///, / \n\n");

	return;
}

void print_information(void) {
	MPI_Printf("***********************************************************\n");
	MPI_Printf("Code parameters:\n");
#ifdef MPI_SUPPORT
	MPI_Printf("Using MPI with %d process.\n",NPROC);
#else
	MPI_Printf("MPI disabled\n");
#endif
#ifdef _OPENMP
	#pragma omp parallel 
	{
	nthreads = omp_get_num_threads();
	}
	
	MPI_Printf("Using OpenMP with %d threads.\n", nthreads);
#else
	nthreads = 1;	// Safety value
	MPI_Printf("OpenMP disabled\n");
#endif
	MPI_Printf("(NX,NY,NZ)=\t(%d,%d,%d)\n",NX,NY,NZ);
	MPI_Printf("(LX,LY,LZ)=\t(%f,%f,%f)\n",param.lx,param.ly,param.lz);
	MPI_Printf("Reynolds=\t%f\n\n",param.reynolds);
#ifdef BOUSSINESQ
#ifdef VERTSTRAT
	MPI_Printf("Vertical Boussinesq\n");
#else
	MPI_Printf("Horizontal (x) Boussinesq\n");
#endif
	MPI_Printf("Reynolds_th=\t%f\n",param.reynolds_th);
	MPI_Printf("N2=\t\t%f\n",param.N2);
#else
	MPI_Printf("No Boussinesq\n");
#endif
#ifdef MHD
	MPI_Printf("\nMHD enabled\n");
#ifdef ELSASSER_FORMULATION
	MPI_Printf("Using Elsasser Formulation\n");
#endif
	if(param.init_mean_field) {
		MPI_Printf("BX0=\t\t%f\n",param.bx0);
		MPI_Printf("BY0=\t\t%f\n",param.by0);
		MPI_Printf("BZ0=\t\t%f\n",param.bz0);
	}
	MPI_Printf("Reynolds_m=\t%f\n",param.reynolds_m);
#else
	MPI_Printf("\nNo MHD\n");
#endif
#ifdef WITH_ROTATION
	MPI_Printf("\nOmega=\t\t%f\n",param.omega);
#endif
#ifdef WITH_SHEAR
	MPI_Printf("Shear=\t\t%f\n",param.shear);
#ifdef TIME_DEPENDANT_SHEAR
	MPI_Printf("Omega_shear=\t\t%f\n",param.omega_shear);
#endif
#else
	MPI_Printf("No Shear\n");
#endif
#ifdef WITH_PARTICLES
	MPI_Printf("\nHave particles\n");
	MPI_Printf("Nparticles=\t%d\n",param.particles_n);
	MPI_Printf("Part' mass=\t%f\n",param.particles_mass);
	MPI_Printf("Part' s-time=\t%f\n",param.particles_stime);
#endif
	MPI_Printf("\nT_initial=\t%f\n",param.t_initial);
	MPI_Printf("T_final=\t%f\n",param.t_final);
	MPI_Printf("Toutput_time=\t%f\n",param.toutput_time);
	MPI_Printf("Toutput_flow=\t%f\n",param.toutput_flow);
	MPI_Printf("Toutput_dump=\t%f\n",param.toutput_dump);
	if(param.restart)
		MPI_Printf("Using Restart Dump\n");
	if(param.antialiasing)
		MPI_Printf("Using Antialiasing 2/3 Rule\n");
	else
		MPI_Printf("No antialiasing\n");
	MPI_Printf("***********************************************************\n");
	return;
}

int main(int argc, char *argv[]) {
#ifdef MPI_SUPPORT
	int size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&NPROC);
	
	// Some consistancy check
	if(NX/NPROC < 1) ERROR_HANDLER( ERROR_CRITICAL, "NX should be a multiple of the number of process.");
	if(NY/NPROC < 1) ERROR_HANDLER( ERROR_CRITICAL, "NY should be a multiple of the number of process.");
	if(NZ < 2) ERROR_HANDLER( ERROR_CRITICAL, "You need NZ > 1, even if you want to do 2D simulations!.");
#ifdef WITH_2D
	ERROR_HANDLER( ERROR_CRITICAL,"MPI is not supported in 2D for the moment.");
#endif
#else
	rank=0;
#endif
	print_logo();
	MPI_Printf("The Snoopy code v6.0\n");
	MPI_Printf("Copyright (c) 2004-2010 Geoffroy Lesur (University of Cambridge, UK)\n");
	MPI_Printf("This program comes with ABSOLUTELY NO WARRANTY;\n");
	MPI_Printf("This is free software, and you are welcome to\n");
    MPI_Printf("redistribute it under certain conditions.\n");
	MPI_Printf("Compiled on %s at %s\n",__DATE__ , __TIME__);
	read_config();
	print_information();
	MPI_Printf("Initializing...\n");
	init_common();
#ifdef WITH_PARTICLES
	init_particles();
#endif
	init_gfft();
	init_output();
		
	init_flow();

	MPI_Printf("Calling mainloop... touch status, output, dump or stop to print more information.\n");
	please_wait();
	mainloop(param.t_initial, param.t_final);
	
	finish_output();
	finish_gfft();
#ifdef WITH_PARTICLES
	finish_particles();
#endif
	finish_common();
	
	MPI_Printf("Terminated.\n");
#ifdef MPI_SUPPORT
	MPI_Finalize();
#endif
	return(0);
}

	
