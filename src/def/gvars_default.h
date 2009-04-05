
#define		NX				96			/**< X Dimension in real space. Must be multiples of NPROC when using MPI.*/
#define		NY				96			/**< Y Dimension in real space. Must be multiples of NPROC when using MPI.*/
#define		NZ				96			/**< Z Dimension in real space. */

#define		NTHREADS		2			/**< Number of OpenMP Thread. Useful only if OpenMP is activated in the Makefile */

#ifndef NPROC
#define		NPROC			1			/**< Number of MPI Process. Useful only if MPI is activated in the Makefile */
#endif

#define MHD								/**< Uncomment to activate MHD*/

#define		LX				1.0			/**< Box length in X*/
#define		LY				1.0			/**< Box length in Y*/
#define		LZ				1.0			/**< Box length in Z*/

#define		CFL				1.5			/**< CFL safety factor. Should be smaller than sqrt(3) for RK3 to be stable.*/
#define		SAFETY_SOURCE	0.2			/**< Safety factor for SHEAR, Coriolis and Boussinesq terms (should be ~0.2 for accuracy) */

#define		REYNOLDS		1000.0		/**< Reynolds number (actully the inverse of the viscosity) */
#define		REYNOLDS_TH		1000.0		/**< Thermal Reynolds number (actully the inverse of the thermal diffusivity)  Used only when Boussinesq is on*/
#define		REYNOLDS_M		1000.0		/**< Magnetic Reynolds number (actully the inverse of the resistivity)  Used only when MHD is on*/

//#define		BOUSSINESQ				/**< Uncomment to activate Boussinesq */
#define		N2				(-1.0)		/**< Brunt Vaissala frequency squared */
//#define		VERTSTRAT					/**< Vertical stratification. Otherwise, Boussinesq stratification is in X */

#define		OMEGA			(2.0/3.0)	/**< Vertical rotation rate (if Shear=1, Keplerian if found for (2.0/3.0) */

#define		WITH_SHEAR					/**< Uncomment to activate mean SHEAR */
#define		SHEAR			1.0			/**< Shear rate */
//#define		TIME_DEPENDANT_SHEAR		/**< Enable Time dependant shear */
#define		OMEGA_SHEAR		10.0			/**< Pulsation of the time dependant shear */

//#define		FORCING						/**< Uncomment to use internal forcing of the velocity field (see forcing in timestep.c) */

#define		T_INITIAL		0.0			/**< Initial time of the simulation */
#define		T_FINAL			5000.0		/**< Simulation will stop if it reaches this time */

#define		TOUTPUT_TIME	0.1			/**< Time between two outputs in the timevar file */
#define		TOUTPUT_FLOW	1.0			/**< Time between two snapshot outputs */
#define		TOUTPUT_DUMP	100.0			/**< Time between two restart dump outputs (restart dump are erased) */

//#define		RESTART					/**< Uncomment to ask for a restart */

//#define		FORTRAN_OUTPUT_ORDER	/**< If VTK_OUTPUT is commented, the code will output binary in C-major order. Uncomment this to get output in FORTRAN-major order (doesn't work with MPI) */
#define		VTK_OUTPUT					/**< Use VTK legacy files for output instead of raw binaries (useful with paraview) */

#define		INTERFACE_CHECK	5			/**< Number of loops between two checks for a user input. On slow filesystems, increase this number */
//#define		INTERFACE_OUTPUT_FILE	/**< Uncomment this option to force code outputs to a file instead of the screen */

#define ANTIALIASING					/**< 2/3 Antialisaing rule. Could be removed if you assume is unimportant (untested feature). */

#define		FFT_PLANNING	FFTW_MEASURE  /**< can be either FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT or FFTW_EXHAUSTIVE (see fftw3 doc). Measure leads to longer initialisation of fft routines */

/**********************************************************************
******* Initial Conditions ********************************************
***********************************************************************/

//#define		INIT_VORTEX							/**< Add a 2D Kida vortex in the box. Assumes S=1. Requires b>a*/
#define		VORTEX_A					0.1		/**< x dimension of the vortex */
#define		VORTEX_B					0.3		/**< y dimension of the vortex */

//#define		INIT_SPATIAL_STRUCTURE				/**< Init a user-defined spatial structure */
	
#define		INIT_LARGE_SCALE_NOISE				/**< Init a large scale random noise (4 largest modes) */
#define		PER_AMPLITUDE_LARGE			1.0		/**< Amplitude of the large scale random noise */

//#define		INIT_WHITE_NOISE					/**< Init a random white noise on all the modes */
#define		PER_AMPLITUDE_NOISE			1.0		/**< total amplitude of the perturbation */

#define		INIT_MEAN_FIELD						/**< Force the mean magnetic field to a given value. */
#define		BX0							0.0		/**< Mean magnetic field in the x direction */
#define		BY0							0.0		/**< Mean magnetic field in the y direction */
#define		BZ0							0.01		/**< Mean magnetic field in the z direction */

/***********************************************************************
*** Ordinary users should not modified anything below this point *******
************************************************************************/

// Ignore NPROC if MPI_SUPPORT is disabled
#ifndef	MPI_SUPPORT
#undef		NPROC
#define		NPROC			1
#endif

#define		NTOTAL			NX * NY * NZ	/**< Total number of grid points over all the processors */

#define		NX_COMPLEX		NX				/**< Number of complex point in X over all the processes (Fourier space) */
#define		NY_COMPLEX		NY				/**< Number of complex point in Y over all the processes (Fourier space) */
#define		NZ_COMPLEX		(NZ / 2 + 1)	/**< Number of complex point in Z over all the processes (Fourier space) */

#define		NTOTAL_COMPLEX	(NX_COMPLEX * NY_COMPLEX * NZ_COMPLEX / NPROC)	/**< Number of complex points in one MPI process */

#define		PRECISION		double											/**< Precision of the code (Float has not been tested, and won't work with MPI) */

#define		IDX3D			(k + j * NZ_COMPLEX + NZ_COMPLEX * NY_COMPLEX * i)  /**< General wrapper for 3D Arrays */