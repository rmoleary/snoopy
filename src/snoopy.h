#ifndef __SNOOPY_H__
#define __SNOOPY_H__

#include <math.h>
#include <complex.h>
#include <fftw3.h>

#ifdef MPI_SUPPORT
#include <mpi.h>
#ifdef FFTW3_MPI_SUPPORT
#include <fftw3-mpi.h>
#endif
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "gvars.h"
#include "error.h"


#ifdef MPI_SUPPORT
#define		MPI_Printf			if (rank==0) printf
#else
#define		MPI_Printf			printf
#endif


// fix NPROC if MPI_SUPPORT is disabled
#ifndef	MPI_SUPPORT
#define		NPROC			1
#endif

#define		NTOTAL			NX * NY * NZ	/**< Total number of grid points over all the processors */

#define		NX_COMPLEX		NX				/**< Number of complex point in X over all the processes (Fourier space) */
#define		NY_COMPLEX		NY				/**< Number of complex point in Y over all the processes (Fourier space) */
#define		NZ_COMPLEX		(NZ / 2 + 1)	/**< Number of complex point in Z over all the processes (Fourier space) */

#define		NTOTAL_COMPLEX	(NX_COMPLEX * NY_COMPLEX * NZ_COMPLEX / NPROC)	/**< Number of complex points in one MPI process */

#define		IDX3D			(k + j * NZ_COMPLEX + NZ_COMPLEX * NY_COMPLEX * i)  /**< General wrapper for 3D Arrays */

// Structures

struct Field {
	double complex *vx;
	double complex *vy;
	double complex *vz;
#ifdef BOUSSINESQ
	double complex *th;
#endif
#ifdef MHD
	double complex *bx;
	double complex *by;
	double complex *bz;
#endif
};

struct Parameters {
	// Physics Parameters
	double lx;				/**< Box length in X*/
	double ly;				/**< Box length in Y*/
	double lz;				/**< Box length in Z*/
	
	double reynolds;		/**< Reynolds number (actully the inverse of the viscosity) */
	double reynolds_m;		/**< Magnetic Reynolds number (actully the inverse of the resistivity)  Used only when MHD is on*/
	
	double reynolds_th;		/**< Thermal Reynolds number (actully the inverse of the thermal diffusivity)  Used only when Boussinesq is on*/
	double N2;				/**< Brunt Vaissala frequency squared */
	
	double omega;			/**< Vertical rotation rate (if Shear=1, Keplerian if found for (2.0/3.0). Only when WITH_ROTATION is on. */
	
	double shear;			/**< Shear rate (only when WITH_SHEAR is on) */
	
	double omega_shear;		/**< Pulsation of the time dependant shear (only when WITH_SHEAR and TIME_DEPENDANT_SHEAR is on) */
	
	// Code parameters
	
	double cfl;				/**< CFL safety factor. Should be smaller than sqrt(3) for RK3 to be stable.*/
	double safety_source;	/**< Safety factor for SHEAR, Coriolis and Boussinesq terms (should be ~0.2 for accuracy) */
	
	double t_initial;		/**< Initial time of the simulation */
	double t_final;			/**< Simulation will stop if it reaches this time */
	double max_t_elapsed;	/**< Maximum elapsed time (in hours). Will stop after this elapsed time */
	
	int    interface_check;	/**< Number of loops between two checks for a user input. On slow filesystems, increase this number */
	int    interface_output_file;	/**< Set this option to force code outputs to a file instead of the screen */
	
	int    force_symmetries;	/**< set to enforce spectral symmetries and mean flow to zero. Useful when N^2 or kappa^2 < 0. (see enforce_symm() )*/
	int    symmetries_step;		/**< Number of loops between which the symmetries are enforced. Should be around ~20 for Boussinesq convection*/
	
	int    antialiasing;		/**< 2/3 Antialisaing rule. Could be removed if you assume is unimportant (untested feature). */
	
	int    restart;
	
	// Output parameters
	double toutput_time;		/**< Time between two outputs in the timevar file */
	double toutput_flow;		/**< Time between two snapshot outputs */
	double toutput_dump;		/**< Time between two restart dump outputs (restart dump are erased) */
	
	int    fortran_output_order;	/**< If vtk_output is off, the code will output binary in C-major order. Set this option to get output in FORTRAN-major order (doesn't work with MPI) */

	int    vtk_output;			/**< Use VTK legacy files for output instead of raw binaries (useful with paraview) */
	
	// initial conditions
	int	   init_vortex;			/**< Add a 2D Kida vortex in the box. Assumes S=1. Requires b>a*/
	double vortex_a;			/**< x dimension of the vortex */
	double vortex_b;			/**< y dimension of the vortex */
	
	int    init_spatial_structure;	/**< Init a user-defined spatial structure (see initflow.c) */
	
	int	   init_large_scale_noise;	/**< Init a large scale random noise (4 largest modes) */
	double per_amplitude_large;		/**< Amplitude of the large scale random noise */
	double noise_cut_length;		/**< Wavelength over which the noise is applied */
	
	int    init_white_noise;		/**< Init a random white noise on all the modes */
	double per_amplitude_noise;		/**< total amplitude of the perturbation */
	
	int    init_mean_field;			/**< Force the mean magnetic field to a given value. */
	double bx0;						/**< Mean magnetic field in the x direction */
	double by0;						/**< Mean magnetic field in the y direction */
	double bz0;						/**< Mean magnetic field in the z direction */
	
	int    init_dump;				/**< Use a dump file as an initial condition (everything else, including t, noutput (...) is reinitialised) */
	
	int	   init_bench;				/**< Init the Benchmark initial conditions */
};

#endif
