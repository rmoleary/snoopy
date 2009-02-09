#include <stdlib.h>
#include <stdio.h>
#include "error.h"

void error_h (const int ErrorType,
			  const char ErrorMessage[], 
			  const char ErrorRoutine[], 
			  const int line, 
			  const char Filename[] ) {
	printf("*************************************************\n");
	if ( ErrorType == ERROR_WARNING )
		printf("Warning in ");
	if ( ErrorType == ERROR_NONCRITICAL )
		printf("Non Critical Error in ");
	if ( ErrorType == ERROR_CRITICAL )
		printf("Critical Error in ");
	printf("File: ");
	printf(Filename);
	printf(" Routine: ");
	printf(ErrorRoutine);
		printf(" Line: %d \n",line);
	printf(ErrorMessage);
	
	if(ErrorType == ERROR_CRITICAL) {
		printf("\n Terminating.\n");
		exit(1);
	}
	else {
		printf("\nResuming execution...\n");
		printf("*************************************************\n");
	}
	return;
}
