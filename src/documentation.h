/*! \mainpage The SNOOPY Code Documentation
<CODE>
<PRE> <B>
          .o. 
          |  |    _   ,
        .',  L.-'` `\\ ||
      __\\___,|__--,__`_|__
     |    %     `=`       |
     | ___%%_______________|
     |    `               |
     | -------------------|
     |____________________|
       |~~~~~~~~~~~~~~~~|
       | ---------------|  ,
   \|  | _______________| / /
\. \,\\\\|, .   .   /,  / |///, /</B></PRE></CODE>

	\section intro Introduction
	This is the Snoopy code documentation. Snoopy is a general purpose Spectral solver, solving MHD and Boussinesq equations
	Including Shear and rotation if required.
	It can run on parallel machine using MPI OpenMP or Both at the same time. The ffts are based on FFTW3, using a custom MPI parallelisation,
	or using the alpha version of the MPI implementation included in fftw 3.3
	\subsection snoop Why Snoopy?
	Why not? I felt that since Zeus, Athena and Ramses were already used, a bit of modern history would sound good (and cool). Snoopy actually means nothing (appart from a dog), 
	although one might say that S stands for Spectral, and noopy stands for... (to be completed).
	\subsection what What is this code actually doing?
	Noops
	*/
	
/*!	\page using Using the code
	\section config Configuring the code
	As any open source software, Snoopy comes with automatic configuration tools. The script ./configure allows one to configure the main options of snoopy, creating a customized Makefile
	and producing a default src/gvars.h and snoopy.cfg to setup the physics. The following options may be used with the configure script:
	- --with-problem=PROB: Initialize the code with problem PROB. This initializes snoopy.cfg and src/gvars.h according to PROB. See the section \ref problem for a standard list of problems.
	- --enable-mpi: Enable MPI parallelization
	- --enable-fftw-mpi: Enable *experimental* support of MPI found in fftw3.3alpha. These routines replace the custom transpose routines found in Snoopy and are generally more efficient. NB: MPI support in fftw3 is still under developpement and untested. Unless you're sure of what you're doing, keeping this option off is safer.
	- --enable-openmp: Enable OpenMP support. Useful on Intel Core** processors. This option requires OpenMP support in FFTW3 (see fftw3 doc).
	- --enable-debug: Enable debug outputs. This option will override the optimisation flags and create LOTS of outputs when the code is run.
	- CC=xxx: force the compiler to be xxx.
	- CFLAGS=xxx: force the compilation flags to be xxx.
	- LDFLAGS=xxx: force the link flags to be xxx.
	- FFTPATH=xxx: specify where the FFTW 3 libraries are located (if not found in the default path). This option assumes the libraries files (libfftw3.a...) are in xxx/lib and the include files (fftw3.h...)
	  are in xxx/include. 

	  
	Example: One wants to configure Snoopy with openMP using the Intel compiler "icc". The fftw library is located in /opt (/opt/lib and /opt/include) and one wants to initialize an MRI problem. One has to configure Snoopy with:
\verbatim
./configure CC=icc FFTPATH=/opt --enable-openmp --with-problem=mri
\endverbatim
	
	Once the configure script has finished, you normally don't need to run it again, except if you want to change one of these options. 
	
	\section test Testing the code
	
	A standard test can be run typing "make check". This test, although not physically meaningful (magnetized 2D vortex with unstable boussinesq stratification), switches on almost all the routines of the code and therefore checks if everything is running
	as it should. Make check compiles the code with a benchmark configuration (saving your gvars.h if you have already made modifications), runs it and compares the outputs to a standard
	output. If the code behaves normally, "make check" should exit without any error.
		
	\section problem Problem setup
	A problem correponds to a header file src/gvars.h and a config file snoopy.cfg. Templates of these files for several problems are located in src/problem. Each problem (corresponding to a subdirectory in src/problem) can be initialized 
	using --with-problem=PROB of the configure script or alternatively moving by hand gvars.h in ROOT/src and snoopy.cfg in ROOT/. The file gvars.h contains major options requiring a recompilation of the code (make). The file snoopy.cfg is read at runtime
	and should be accessible by at least process of rank 0 (for MPI runs). The available options in gvars.h and snoopy.cfg are described in the \ref code_config documentation. The following problems are available by default
	(the user can create new problems with new directories in src/problem).
	
	- default
	- bench
	- mri
	- convection
	- couette
	
	\section interface Code interface
	While the code is running, it's possible to know what's happening in real time using the so-called interface (located in interface.c). Typically, one creates a file with a filename
	corresponding to one of the possible commands (e.g using the command "touch" on UNIX, as "touch status"). Once the command has been executed, the code deletes the file.
	The available commands are as follow:
		- status: show the status of the code, including current time, current time step and code speed.
		- output: Immediatly print the statistical informations in timevar, output one snapshot and a dump file (whatever are the variables in gvars.h).
		- dump: Immediatly output a dump file.
		- stop: Immediatly output a dump file and exit the code.
	
	It is possible to redirect the display outputs of the interface to a file using the INTERFACE_OUTPUT_FILE option in gvars.h. This is useful when one wants to run in batch 
	mode on a cluster. For performances reasons, the code doesn't check at each loop if the user has created a command file. Instead, it checks every INTERFACE_CHECK loops. A larger
	INTERFACE_CHECK results in a smaller overhead but a longer delay between the command file creation and the actual output.
*/

/*! \page code_config Code configuration
	Snoopy configuration is divided in two files: gvars.h and snoopy.cfg, which are described below.
    \section gvars File gvars.h
	To activate or deactivate a feature, one uncomments or comments (//) the corresponding #define. Any modification made to this file requires a recompilation of the code (make)
	to include them in the code. Here is an example of a gvars.h file (an updated version of this file may be found in src/problem/defaut/gvars.h):
	\verbatim
	
#define NX              96            // X Dimension in real space. Must be multiples of 
                                      // NPROC when using MPI.
#define	 NY              96            // Y Dimension in real space. Must be multiples of 
                                      // NPROC when using MPI.
#define	 NZ              96            // Z Dimension in real space. 

#define	 MHD                           // Uncomment to activate MHD

#define	 BOUSSINESQ                    // Uncomment to activate Boussinesq
#define	 VERTSTRAT                     // Vertical stratification. Otherwise,
                                      // Boussinesq stratification is in X

#define	 WITH_ROTATION                 // Uncomment to add rotation around the z axis.

#define	 WITH_SHEAR                    // Uncomment to activate shear (U_y(x))
#define	 TIME_DEPENDANT_SHEAR          // Enable Time dependant shear

#define	 FORCING                       // Uncomment to use internal forcing of the 
                                      // velocity field (see forcing in timestep.c)

#define	 FFT_PLANNING    FFTW_MEASURE  // can be either FFTW_ESTIMATE, FFTW_MEASURE, 
                                      // FFTW_PATIENT or FFTW_EXHAUSTIVE 
                                      // (see fftw3 doc). Measure leads to longer 
                                      // initialisation of fft routines
	\endverbatim
	
	\section configsnoopy File snoopy.cfg
	The file snoopy.cfg is located where the executable is located and is read at run time. Snoopy uses a variation of the 
	library <A HREF=http://www.hyperrealm.com/libconfig/>libconfig</A> to read these files. A standard snoopy file is divided into
	3 blocks: physics, code and init corresponding to physics parameters, code parameters and initial conditions. If any of the described 
	parameter (or even block) is ommited, the default value (as described below) will be used. This is a commented example of snoopy.cfg
	containing all the possible parameters assigned to their default value (an updated version of this file may be found in src/problem/defaut/snoopy.cfg):
	\verbatim
	# Example of a Snoopy configuration file

configname = "Default Snoopy configuration file";

physics:                             // Physics parameters
{
	boxsize = (1.0, 1.0, 1.0);       // Box length in X, Y and Z
	
	reynolds = 1.0;                  // Reynolds number (actully the inverse of the viscosity)
	reynolds_magnetic = 1.0;         // Magnetic Reynolds number (actully the inverse of the resistivity).  Used only when MHD is on
	reynolds_thermic = 1.0;          // Thermal Reynolds number (actully the inverse of the thermal diffusivity).  Used only when Boussinesq is on
	
	brunt_vaissala_squared = 0.0;    // Brunt Vaissala frequency squared. Used only when Boussinesq is on
	
	omega = 0.0;                     // Vertical rotation rate (if Shear=1, Keplerian if found for 2.0/3.0). Used only when WITH_ROTATION is on
	
	shear = 0.0;                     // Shear rate. Used only when WITH_SHEAR is on.
	omega_shear = 0.0;               // Pulsation of time dependant shear. Used only when both WITH_SHEAR and TIME_DEPENDANT_SHEAR are on.
};

//-------------------------------------------------------------------------------------------------------------------------

code:                                // code parameters
{
	cfl = 1.5;                       // CFL safety factor. Should be smaller than sqrt(3) for RK3 to be stable.
	safety_source = 0.2;             // Safety factor for SHEAR, Coriolis and Boussinesq terms (should be ~0.2 for accuracy)
	
	t_initial = 0.0;                 // Initial time of the simulation
	t_final = 1.0;                   // Simulation will stop if it reaches this time
	max_t_elapsed = 1e30;            // Maximum elapsed time (in hours). Will stop after this elapsed time if t_final is not reached.
	
	interface_check = 5;             // Number of loops between two checks for a user input. On slow filesystems, increase this number 
	interface_output_file = false;   // Set to true to force interface outputs to a file instead of displaying them
	
	force_symmetries = false;        // Uncomment to enforce spectral symmetries and mean flow to zero. Useful when N^2 or kappa^2 < 0. (see enforce_symm() )
	symmetries_step = 20;            // Number of loops between which the symmetries are enforced. Should be around ~20 for Boussinesq convection.
	
	antialiasing = true;             // 2/3 Antialisaing rule. Could be removed if you assume is unimportant (untested feature).
	
	restart = false;                 // set to true to restart from a dump file. If no dump file is found, this option has no effect.
};

//-------------------------------------------------------------------------------------------------------------------------

output:	                             // output parameters
{
	timevar_step = 1.0;	             // Time between two outputs in the timevar file
	snapshot_step = 1.0;             // Time between two snapshot outputs
	dump_step = 1.0;                 // Time between two restart dump outputs (restart dump are erased)
	
	vtk_output = true;               // Use VTK legacy files for output instead of raw binaries (useful with paraview)
	fortran_output_order = false;    // If vtk_output is disabled, the code will output binary in C-major order. Uncomment this to get outputs in FORTRAN-major order (doesn't work with MPI)
	pressure = false;				 // Output the pressure field in the 3D snapshots
};

//-------------------------------------------------------------------------------------------------------------------------

init:                                // Initial conditions parameters
{
	vortex:                          // Add a 2D Kida vortex in the box. Assumes S=1. Requires b>a
	{
		enable = false;              // Set this to true to enable the vortex
		a = 1.0;                     // x dimension of the vortex
		b = 2.0;                     // y dimension of the vortex
	};
	large_scale_noise:               // Init a large scale random noise down to cut_length
	{
		enable = false;	             // set this to true to enable large scale noise
		amplitude = 0.0;             // noise amplitude
		cut_length = 0.0;            // Wavelength over which the noise is applied
	};
	white_noise:                     // Init a random noise at all scales
	{
		enable = false;	             // set this to true to enable white noise
		amplitude = 0.0;             // noise amplitude
	};
	mean_field:	                     // Force the mean magnetic field to a given value.
	{
		enable = false;	             // Set this to true to enable mean field
		bx0 = 0.0;                   // Mean magnetic field in the x direction
		by0 = 0.0;                   // Mean magnetic field in the y direction
		bz0 = 0.0;                   // Mean magnetic field in the z direction
	};
	spatial_structure = false;       // set this to true to init a user_defined spatial structure (see initflow.c)
	dump = false;                    // set this to true to use a dump file as an initial condition (this is NOT a restart option!)
	bench = false;                   // set this to true to init a benchmark initial condition.
};
		

	\endverbatim
*/


/*!	\page outputs Code outputs
	\section timevar The timevar File
	A text file containing several averaged quantities (see output.c for more details).
	\section snap Snapshots
	Snapshots can be found in raw binary files (.raw) or in vtk legacy format (.vtk, default) in the data directory. The output format is set in gvars.h. VTK files can are read natively
	by Paraview 2-3 or Visit, available for free on the web. Several Matlab script are under developpment to read these files.
	
	\section dump Restart dump files
	Binary restart file. See output.c for a complete description.
*/

