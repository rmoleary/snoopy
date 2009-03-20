#include <stdlib.h>
#include <stdio.h>
#include "common.h"
#include "error.h"

void error_h (const int ErrorType,
			  const char ErrorMessage[], 
			  const char ErrorRoutine[], 
			  const int line, 
			  const char Filename[] ) {
	MPI_Printf("*************************************************\n");
	if ( ErrorType == ERROR_WARNING )
		MPI_Printf("Warning in ");
	if ( ErrorType == ERROR_NONCRITICAL )
		MPI_Printf("Non Critical Error in ");
	if ( ErrorType == ERROR_CRITICAL )
		MPI_Printf("Critical Error in ");
	MPI_Printf("File: ");
	MPI_Printf(Filename);
	MPI_Printf(" Routine: ");
	MPI_Printf(ErrorRoutine);
		MPI_Printf(" Line: %d \n",line);
	MPI_Printf(ErrorMessage);
	
	if(ErrorType == ERROR_CRITICAL) {
		MPI_Printf("\n Terminating.\n");
		MPI_Printf("*************************************************\n");
#ifdef MPI_SUPPORT
		MPI_Abort(MPI_COMM_WORLD,1);
#endif
		exit(1);
	}
	else {
		MPI_Printf("\nResuming execution...\n");
		MPI_Printf("*************************************************\n");
	}
	return;
}
