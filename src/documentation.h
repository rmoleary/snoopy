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
	Type make in the root directory. If you're lucky, everything will work fine. If you're not, then you'll have to create a configuration for your architecture
	in Makefile. If the compilation went successfully, just type ./snoopy and that's it!
	\section makefile The Makefile
	The big options can be switched on and off in the Makefile. Several self explanatory flags may be found in the begining of Makefile to activate OpenMP, MPI, Debug and
	MPI support from fftw3 Library (MPI support in fftw3 is still on developpement and untested. Unless you're sure of what you're doing, keeping this option off is safer).
	
	By default, if no configuration is found (src/gvars.h), the makefile will create a default configuration file using
	the template located in src/def/gvars.h. One can then modify the configuration file (make won't recreate it). It is possible to create a default configuration file
	without actually compiling the code with make def.
	
	\section interface Code interface
	While the code is running, it's possible to know what's happening in real time using the so-called interface (located in interface.c). Typically, one creates a file with a filename
	corresponding to one of the possible commands (e.g using the command touch on UNIX, as "touch status"). Once the command has been executed, the code erases the file.
	The avaiable commands are as follow:
		- status: show the status of the code, including current time, current time step and code speed.
		- output: Immediatly print the statistical informations in timevar, output one snapshot and a dump file (whatever are the variables in gvars.h).
		- dump: Immediatly output a dump file.
		- stop: Immediatly output a dump file and exit the code cleanly.
	
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
