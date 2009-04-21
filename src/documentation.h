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
	\section first_start Let's go
	Type make in the root directory. If you're lucky, everything will work fine. If you're not, then you'll have to create a configuration for your architecture. If the compilation went successfully, just type ./snoopy and that's it!
	\section makefile The Makefile
	The makefile is divided in several files. "makefile" contains the rules to create a configuration file for several architectures. The file
	rules.mk contains the compilation rules for the code. This should not be modified by ordinary users. Finally, config.mk
	contains architecture-dependant flags, paths and compiler options. This file is meant to be modified by the user. Several default configuration
	files are stored in src/def/*.mk. 
	
	When you first decompress or download the code, no config.mk is found in the root directory. To create one, type "make config". If you're lucky,
	a configuration file relevant for your system will be found. if not, the default configuration file will be used and you'll have to edit the created
	config.mk according to your system.
	
	The big options can be switched on and off in config.mk. Several self explanatory flags may be found to enable OpenMP, MPI, Debug and
	MPI support from fftw3 Library (MPI support in fftw3 is still under developpement and untested. Unless you're sure of what you're doing, keeping this option off is safer).
	
	Make config also initializes a default configuration for the code, creating the file src/gvars.h. One can then modify gvars.h for a specific problem 
	(make won't recreate it).
	
	Finally, note that if you're lucky enough, typing make should initialize makefiles and configuration properly and should compile the code
	without problem.
	
	\section test_case Test case
	One can compute a standard test with Snoopy. This standard test, although not physically meaningful (magnetized 2D vortex with unstable boussinesq stratification), switches on almost all the routines of the code and therefore checks if everything is running
	as it should. To do the full standard test, you can just type make test. This script will overwrite src/gvars.h, run the code, and compare the resulting output files with the reference
	ones stored in src/def. If no problem is found, the script should finish successfully.
	
	If your configuration needs a specific config.mk, you can first setup a config.mk with "make config", edit config.mk and then use "make test". 
	
	\section interface Code interface
	While the code is running, it's possible to know what's happening in real time using the so-called interface (located in interface.c). Typically, one creates a file with a filename
	corresponding to one of the possible commands (e.g using the command "touch" on UNIX, as "touch status"). Once the command has been executed, the code erases the file.
	The available commands are as follow:
		- status: show the status of the code, including current time, current time step and code speed.
		- output: Immediatly print the statistical informations in timevar, output one snapshot and a dump file (whatever are the variables in gvars.h).
		- dump: Immediatly output a dump file.
		- stop: Immediatly output a dump file and exit the code.
	
	It is possible to redirect the screen outputs of the interface to a file using the INTERFACE_OUTPUT_FILE option in gvars.h. This is useful when one wants to run in batch 
	mode on a cluster. For performances reasons, the code doesn't check at each loop if the user has created a command file. Instead, it checks every INTERFACE_CHECK loops. A larger
	INTERFACE_CHECK results in a smaller overhead but a longer delay between the command file creation and the actual output.
*/

/*!	\page outputs Code outputs
	\section timevar The timevar File
	Blah
	\section snap Snapshots
	Blah
	\section dump Restart dump files
	Blah
*/
