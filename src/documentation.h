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
	and producing a default gvars.h to setup the physics. The following options can be useful when using the configure script:
	- --enable-mpi: Enable MPI parallelization
	- --enable-fftw-mpi: Enable *experimental* support of MPI found in fftw3.3alpha. These routines replace the custom transpose routines found in Snoopy and are generally more efficient. NB: MPI support in fftw3 is still under developpement and untested. Unless you're sure of what you're doing, keeping this option off is safer.
	- --enable-openmp: Enable OpenMP support. Useful on Intel Core** processors. This option requires OpenMP support in FFTW3 (see fftw3 doc).
	- --enable-debug: Enable debug outputs. This option will override the optimisation flags and create LOTS of outputs when the code is run.
	- CC=xxx: force the compiler to be xxx.
	- CFLAGS=xxx: force the compilation flags to be xxx.
	- LDFLAGS=xxx: force the link flags to be xxx.
	- FFTPATH=xxx: specify where the FFTW 3 libraries are located (if not found in the default path). This option assumes the libraries files (libfftw3.a...) are in xxx/lib and the include files (fftw3.h...)
	  are in xxx/include. 
	  
	Example: I want to configure Snoopy with openMP using the Intel compiler "icc". My fftw library is located in /opt (/opt/lib and /opt/include). I will type:
	
	./configure CC=icc FFTPATH=/opt --enable-openmp
	
	Once the configure script has finished, you normally don't need to run it again, except if you want to change one of these options. 
	
	\section test Testing the code
	
	A standard test can be run typing "make check". This test, although not physically meaningful (magnetized 2D vortex with unstable boussinesq stratification), switches on almost all the routines of the code and therefore checks if everything is running
	as it should. Make check compiles the code with a benchmark configuration (saving your gvars.h if you have already made modifications), runs it and compares the outputs to a standard
	output. If the code behaves normally, "make check" should exit without any error. Not that make check is not yet totally compatible with MPI (the number of process can't be set properly). 
	If you want to test an MPI version of the code, you should follow this procedure:
	
	- make bench
	- make
	- make benchclean
	- mpirun -np xx ./snoopy (where xx is the number of process you want to use)
	- diff timevar src/def/timevar_bench 
	
	\section physics Physics, grid and output setup
	All these parameters are found in file src/gvars.h. After a modification, the code has to be compiled (make). Note that a default src/gvars.h is initialized by the configure script. Once it is created, this file is not
	deleted, even if you run configure again. To completely clean the code tree, run "make fullclean" instead (CAUTION: this will also delete all the outputs!).
	
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

/*!	\page outputs Code outputs
	\section timevar The timevar File
	A text file containing several averaged quantities (see output.c for more details).
	\section snap Snapshots
	Snapshots can be found in raw binary files (.raw) or in vtk legacy format (.vtk, default) in the data directory. The output format is set in gvars.h. VTK files can are read natively
	by Paraview 2-3 or Visit, available for free on the web. Several Matlab script are under developpment to read these files.
	
	\section dump Restart dump files
	Binary restart file. See output.c for a complete description.
*/
