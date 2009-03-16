
#define		NX				128			/**< X Dimension in real space. Must be multiples of NPROC when using MPI.*/
#define		NY				128			/**< Y Dimension in real space. Must be multiples of NPROC when using MPI.*/
#define		NZ				128			/**< Z Dimension in real space. */

#define		NTHREADS		2			/**< Number of OpenMP Thread. Useful only if OpenMP is activated in the Makefile */

#ifndef NPROC
#define		NPROC			1			/**< Number of MPI Process. Useful only if MPI is activated in the Makefile */
#endif

//#define MHD							/**< Uncomment to activate MHD*/

#define		LX				1.0			/**< Box length in X*/
#define		LY				1.0			/**< Box length in Y*/
#define		LZ				1.0			/**< Box length in Z*/

#define		CFL				1.5			/**< CFL safety factor. Should be smaller than sqrt(3) for RK3 to be stable.*/
#define		REYNOLDS		1000.0		/**< Reynolds number (actully the inverse of the viscosity) */
#define		REYNOLDS_TH		1000.0		/**< Thermal Reynolds number (actully the inverse of the thermal diffusivity)  Used only when Boussinesq is on*/
#define		REYNOLDS_M		1000.0		/**< Magnetic Reynolds number (actully the inverse of the resistivity)  Used only when MHD is on*/

//#define		BOUSSINESQ				/**< Uncomment to activate Boussinesq */
#define		N2				(-1.0)		/**< Brunt Vaissala frequency squared */
#define		VERTSTRAT					/**< Vertical stratification. Otherwise, Boussinesq stratification is in X */

#define		OMEGA			(2.0/3.0)	/**< Vertical rotation rate (if Shear=1, Keplerian if found for (2.0/3.0) */

//#define		WITH_SHEAR				/**< Uncomment to activate mean SHEAR */
#define		SHEAR			1.0			/**< Shear rate */

#define		BX0				1.0			/**< Mean magnetic field in the x direction */
#define		BY0				0.0			/**< Mean magnetic field in the y direction */
#define		BZ0				0.0			/**< Mean magnetic field in the z direction */

#define		PER_AMPLITUDE	1000.0		/**< Initial perturbation amplitude(arbitrary units) */

//#define		FORCING						/**< Uncomment to use internal forcing of the velocity field (see forcing in timestep.c) */

#define		T_INITIAL		0.0			/**< Initial time of the simulation */
#define		T_FINAL			50.0		/**< Simulation will stop if it reaches this time */

#define		TOUTPUT_TIME	0.1			/**< Time between two outputs in the timevar file */
#define		TOUTPUT_FLOW	0.1			/**< Time between two snapshot outputs */
#define		TOUTPUT_DUMP	1.0			/**< Time between two restart dump outputs (restart dump are erased) */

//#define		RESTART					/**< Uncomment to ask for a restart */

//#define		FORTRAN_OUTPUT_ORDER	/**< If VTK_OUTPUT is commented, the code will output binary in C-major order. Uncomment this to get output in FORTRAN-major order (doesn't work with MPI) */
#define		VTK_OUTPUT					/**< Use VTK legacy files for output instead of raw binaries (useful with paraview) */

#define		INTERFACE_CHECK	5			/**< Number of loops between two checks for a user input. On slow filesystems, increase this number */
//#define		INTERFACE_OUTPUT_FILE	/**< Uncomment this option to force code outputs to a file instead of the screen */

#define ANTIALIASING					/**< 2/3 Antialisaing rule. Could be removed if you assume is unimportant (untested feature). */

#define		FFT_PLANNING	FFTW_MEASURE  /**< can be either FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT or FFTW_EXHAUSTIVE (see fftw3 doc). Measure leads to longer initialisation of fft routines */

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

