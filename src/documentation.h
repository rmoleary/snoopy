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



/*! \mainpage The Snoopy Code Documentation
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
	Including Shear, rotation, compressibility and several other physical effects when required.
	It can run on parallel machine using MPI OpenMP or Both at the same time. The ffts are based on FFTW3, using a custom MPI parallelisation,
	or using the alpha version of the MPI implementation included in fftw 3.3
	\section snoop Why Snoopy?
	Why not? I felt that since Zeus, Athena and Ramses were already used, a bit of modern history would sound good (and cool). Snoopy actually means nothing (appart from a dog), 
	although one might say that S stands for Spectral, and noopy stands for... (I still don't know).
	
	\section get_started Getting started
	To start with Snoopy, you will need to read the following documentation pages:
	- \subpage using "Using the code" describes compilation options and test runs
	- \subpage code_config "Code configuration" describes in details the configuration options
	- \subpage outputs "Code outputs" describes the outputs produced by the code
	- \subpage matlab "Matlab scripts" presents the scripts provided in the package
	
	\section news What's new?
	\subsection albert Albert (v6.0)
		- solid particles interaction+feedback (alpha version, not totally tested)
		- Braginskii viscosity
		- Hall effect
		- Isothermal compressibility
		- Subgrid scale models (Chollet-Lesieur) for the hydro cascade
		- Use of fftw-mpi when available
		- Faster MHD algorithm thanks to a formulation using Elssaser variables
		- new layout for the timevar file
		- and many bug fixes...
	
	\subsection jude Jude (v5.0)
		- Free slip boundary conditions to replace the standard periodic BCs are allowed. See boundary.c for details.
		- 2D (x,y) large scale noise implemented
		- Spatially varying Brunt-Vaisala frequency, for non homogeneous stratification.
		- User can ask for vorticity output in the VTK files
		- MRI problem consistant with Hawley, Gammie, Balbus (1995) in terms of time units
		- Added more spectra in spectrum.dat, including transfer spectra and fix a binning problem
		- Possible to add "hard" object in the flow (this is an alpha feature)
		- Subcritical baroclinic instability setup provided (--with-problem=sbi when using the configure script)
		- Several improvement in the I/O routine
		- Include example Matlab scripts in the distribution
		- fixed a bug with some intel compilers concerning the optimisation (-ansi-alias)
		- added ABC forcing (see forcing.c), thanks to E. Rempel
		- fixed a bug in the restart routine when no dump file was found
		
	\subsection igor Igor (v4.0)
		Initial release
	
	
	
	\section lic License
	The Snoopy code is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Snoopy code is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Snoopy code.  If not, see <A HREF=http://www.gnu.org/licenses/> http://www.gnu.org/licenses/ </A>.

	*/
	
/*!	\page using Using the code
	\section config Configuring the code
	As any open source software, Snoopy comes with automatic configuration tools. The script ./configure allows one to configure the main options of snoopy, creating a customized Makefile
	and producing a default src/gvars.h and snoopy.cfg to setup the physics. The following options may be used with the configure script:
	- --with-problem=PROB: Initialize the code with problem PROB. This initializes snoopy.cfg and src/gvars.h according to PROB. See the section \ref problem for a standard list of problems.
	- --enable-mpi: Enable MPI parallelization
	- --enable-fftw-mpi: Enable *experimental* support of MPI found in fftw3.3alpha. These routines replace the custom transpose routines found in Snoopy and are generally more efficient. NB: MPI support in fftw3 is still under developpement and untested.
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
	- sbi (Subcritical baroclinic instability)
	
	\section interface Code interface
	While the code is running, it's possible to know what's happening in real time using the so-called interface (located in interface.c). Typically, one creates a file with a filename
	corresponding to one of the possible commands (e.g using the command "touch" on UNIX, as "touch status"). Once the command has been executed, the code deletes the file.
	The available commands are as follow:
		- status: show the status of the code, including current time, current time step and code speed.
		- output: Immediatly write the statistical informations in timevar, output one snapshot and a dump file.
		- dump: Immediatly output a dump file.
		- stop: Immediatly output a dump file and exit the code.
	
	It is possible to redirect the display outputs of the interface to a file using the interface_output_file option in snoopy.cfg. This is useful when one wants to run in batch 
	mode on a cluster. For performances reasons, the code doesn't check at each loop if the user has created a command file. Instead, it checks every INTERFACE_CHECK loops. A larger
	INTERFACE_CHECK results in a smaller overhead but a longer delay between the command file creation and the actual output.
*/

/*! \page code_config Code configuration
	Snoopy configuration is divided in two files: gvars.h and snoopy.cfg, which are described below.
    \section gvars File gvars.h
	To activate or deactivate a feature, one uncomments or comments (//) the corresponding #define. Any modification made to this file requires a recompilation of the code (make)
	to include them in the code. Here is an example of a gvars.h file (an updated version of this file may be found in src/problem/defaut/gvars.h):
	\verbatim
	
#define		NX       96                 // X Dimension in real space. Must be multiples of NPROC when using MPI.
#define		NY       96                 // Y Dimension in real space. Must be multiples of NPROC when using MPI.
#define		NZ       96                 // Z Dimension in real space. 

#define		MHD	                        // Uncomment to activate MHD
#define		ELSASSER_FORMULATION        // Use Elsasser formulation for MHD terms 

#define		WITH_BRAGINSKII             // Enable Braginskii viscosity. CAUTION: this is incompatible with the Elsasser formulation
#define		WITH_HALL                   // Enable Hall MHD (in developement). Hall effect amplitude is set by x_hall 
#define		WITH_EXPLICIT_DISSIPATION   // use an explicit (RK3) numerical scheme to integrate linear diffusion effects 
		                                // instead of the standard implicit scheme 
#define		SGS                         // Chollet-Lesieur Subgrid scale model for Hydro. This is in alpha version (see timestep.c).

#define		COMPRESSIBLE                // Solve the isothermal compressible equations (incompatible with Boussinesq for obvious physical reasons) 

#define		BOUSSINESQ                  // Uncomment to activate Boussinesq 
#define		VERTSTRAT                   // Vertical stratification. Otherwise, Boussinesq stratification is in X

#define		WITH_ROTATION               // Enable a Coriolis force, assuming the rotation axis is z. The rotation speed is set up in snoopy.cfg 

#define		WITH_SHEAR                  // Uncomment to activate mean SHEAR 
#define		TIME_DEPENDANT_SHEAR        // Enable Time dependant shear 

#define		WITH_2D                     // 2D (x,y) reduction. Enforce NZ=1 and optimise the Fourier transform routines accordingly. 
                                              // Using this option should lead to a speedup>4 compared to the old NZ=2 method.

#define		WITH_PARTICLES              // Enable Lagrangian test particles tracing.

#define		BOUNDARY_C                  // Enforce NON-PERIODIC boundary conditions to the flow using Fourier space symmetries. The boundary 
			                            // condition has to be set up in boundary.c/boundary_c(); This feature is in alpha version. 

#define		FORCING	                    // Uncomment to use internal forcing of the velocity field (see forcing in timestep.c) 

#define		FFT_PLANNING FFTW_MEASURE   // can be either FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT or FFTW_EXHAUSTIVE (see fftw3 doc). 
                                               // Measure leads to longer initialisation of fft routines 


	\endverbatim
	
	\section configsnoopy File snoopy.cfg
	The file snoopy.cfg is located where the executable is located and is read at run time. Snoopy uses a variation of the 
	library <A HREF=http://www.hyperrealm.com/libconfig/>libconfig</A> to read these files. A standard snoopy file is divided into
	5 blocks: physics, particles, code, output and init corresponding to physical parameters, code parameters, input/output and initial conditions. 
	If any of the described parameter (or even block) is ommited, the default value (as given below) will be used. This is a commented example of snoopy.cfg
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
	reynolds_Braginskii = 1.0;       // Braginskii viscosity-based Reynolds number (inverse of nu_B)
	x_Hall = 1.0;                    // Amplitude of the Hall effect. Used only when WITH_HALL is enabled in gvars.h.

	
	brunt_vaissala_squared = 0.0;    // Brunt Vaissala frequency squared. Used only when Boussinesq is on
	
	omega = 0.0;                     // Vertical rotation rate (if Shear=1, Keplerian if found for 2.0/3.0). Used only when WITH_ROTATION is on
	
	shear = 0.0;                     // Shear rate. Used only when WITH_SHEAR is on.
	omega_shear = 0.0;               // Pulsation of time dependant shear. Used only when both WITH_SHEAR and TIME_DEPENDANT_SHEAR are on.
};

//-------------------------------------------------------------------------------------------------------------------------

particles:                           // Test particles parameters. Only used when WITH_PARTICLES is enabled.
{
	n = 1000;                        // Total number of test particles
	mass = 1.0;                      // Particles mass
	stime = 1.0;                     // Particles stopping time
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
	snapshot_step = 1.0;             // Time between two snapshot outputs (VTK format)
	dump_step = 1.0;                 // Time between two restart dump outputs (restart dump are erased)
	
									 
	timevar_vars = ( "t","ev","em","vxmax","vxmin","vymax","vymin","vzmax","vzmin","vxvy", 
								"bxmax","bxmin","bymax","bymin","bzmax","bzmin","bxby",
								"dmin","dmax","dvxvy",
								"thmax","thmin","thvx","thvz","w2","j2","hm" );
								
									 // List of statistical quantities to be written in the timevar file
	
	vorticity = false;               // Output the vorticity field in the 3D snapshots
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
	large_scale_2D_noise:            // Init a large scale random 2D (x,y) noise down to cut_length
	{
		enable = false;              // set this to true to enable large scale 2D noise
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
	\section timevar Averaged quantities (timevar)
	A text file containing several averaged quantities in a column style file. The variables to be written are choosen in snoopy.cfg (see output.timevar_vars). The following strings can
	be used in timevar_vars:
	
		- t: time of the simulation
		- ev: Box averaged kinetic energy (incompressible)
		- vxmax: x velocity maximum
		- vxmin: x velocity minimum
		- vymax: y velocity maximum
		- vymin: y velocity minimum
		- vzmax: z velocity maximum
		- vzmin: z velocity minimum
		- vxvy: xy component of the (incompressible) Reynolds stress (box averaged)
		- hv: kinetic helicity
		- w2: total enstrophy
		
	when MHD is on, the following quantities can also be computed:
	
		- em: Box averaged magnetic energy
		- bxmax: x magnetic field maximum
		- bxmin: x magnetic field minimum
		- bymax: y magnetic field maximum
		- bymin: y magnetic field minimum
		- bzmax: z magnetic field maximum
		- bzmin: z magnetic field minimum
		- bxby: xy component of the Maxwell stress
		- hc: box averaged cross helicity
		- hm: box averaged magnetic helicity
		- j2: box averaged current square
		- az2: box average vertical component of the vector potential squared (conserved in 2D)
		- vxaz: correlation tensor vx*az
	
	when BOUSSINESQ is on, the following quantities can also be computed:
		
		- et: Box averaged thermal energy
		- thmax: potential temperature maximum
		- thmin: potential temperature minimum
		- thvx: correlation th*vx (turbulent transport of heat in x)
		- thvz: correlation th*vz (turbulent transport of heat in z)
		
	when COMPRESSIBLE is on, the following quantities can also be computed:
		
		- dmax: maximum of the density field
		- dmin: minimum of the density field
		- dvxvy: xy component of the (compressible) reynolds stress
		
		
	The average time between two lines in the timevar file is set by output.timevar_step (snoopy.cfg). This file can be read using a Matlab script provided in the Snoopy distribution (see \ref timevol)
	
	\section spectrum Shell integrated spectra (spectrum.dat)
	This text file contains shell integrated spectra computed on-the-fly by the code. It includes kinetic and magnetic spectra plus transfer spectra relevant to MHD turbulence. For a complete description,
	see output.c. These files can be read using a Matlab script provided in the Snoopy distribution (see \ref readspectrum)
	
	The average time between two outputs in spectrum.dat is set by output.timevar_step (snoopy.cfg)
	
	NOTA: the spectrum output routines are still in developement.
	
	\section snap Snapshots
	Snapshots can be written in vtk legacy format (.vtk) in the data directory. VTK files can are read natively
	by Paraview 2-3 or Visit, available for free on the web. Several Matlab script are also provided in the Snoopy distribution to read VTK files (see \ref get_frame). The fields written in the VTK files are by default vx,vy and vz
	plus bx,by,bz when MHD is on, th (potential temperature) when BOUSSINESQ is on and wx,wy,wz (vorticity) when output.vorticity is true (see snoopy.cfg).
	
	The average time between two snapshot outputs is set by output.snapshot_step (snoopy.cfg)
	
	
	\section dump Restart dump files
	Binary restart file. A maximum of two restart dumps can be found in the running directory: dump.dmp and dump_sav.dmp. The former is the more recent dump whereas the latter is the preceeding dump. 
	Older dumps are deleted during the execution. It is possible to restart the code from a restart dump setting the restart option in snoopy.cfg to true. See output.c for a complete description.
	
	The average time between two dump file outputs is set by output.dump_step (snoopy.cfg)
*/

/*! \page matlab Matlab scripts
	Several Matlab scripts are now provided with Snoopy (in the matlab directory). These are basic scripts which can be a useful start for a specific project.
	
	\section timevol Timevol script
	The file timevol.m reads the timevar file produced by snoopy and display several time history plots for the variables chosen in snoopy.cfg (energy, transport, helicity, etc... see \ref timevar).
	
	\section get_frame Get_frame script
	The file get_frame.m reads one VTK file created by snoopy and display several cuts of the flow in the xy, xz and yz plane. This script can be easely modified to display other cuts. 
	The velocity/magnetic fields are stored in a Matlab structure
	which can then be used for specific post-treatment tasks (FFTs, filtering, correlations). The same VTK files can also be analyzed using Paraview.
	
	\section readspectrum ReadSpectrum script
	The file ReadSpectrum.m reads the spectrum.dat file produced by Snoopy and display time-averaged spectra of the simulation. NB: this output is still in beta version.

*/