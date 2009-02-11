#include <stdlib.h>
#include <stdio.h>

#include "common.h"
#include "mainloop.h"
#include "output.h"

// This is the user interface used in C flow
// We assume check_interface is called periodically to check whether the user called for something

// currently, we reconize commands status, output and stop, as new files in the root directory of the code

#define			STATUS_COMMAND		"status"
#define			OUTPUT_COMMAND		"output"
#define			DUMP_COMMAND		"dump"
#define			STOP_COMMAND		"stop"

#define			INTERFACE_OUTPUT	"output_interface.txt"

void open_interface_io(FILE ** iostream) {

#ifdef INTERFACE_OUTPUT_FILE
	*iostream = fopen(INTERFACE_OUTPUT,"w");
#else
	*iostream = stdout;
#endif

	return;
}

void close_interface_io(FILE ** iostream) {
#ifdef INTERFACE_OUTPUT_FILE
	fclose(*iostream);
#endif
	return;
}


void check_interface(const struct Field fldi,
					 const PRECISION t,
					 const PRECISION dt,
					 const int		 nloop) {
	// This routine check the interface file and print the relevant informations
	
	FILE * iostream;
	FILE * command_file;
	
	// STATUS command
	command_file = fopen(STATUS_COMMAND,"r");
	if(command_file) {
		fclose(command_file);
		// We have a status command
		remove( STATUS_COMMAND );
		open_interface_io( &iostream );
		
		fprintf(iostream,"STATUS command called.\n");
		fprintf(iostream,"t=%e, dt=%e, nloop=%d\n E=%e\n",t,dt,nloop,energy(fldi.vx)+energy(fldi.vy)+energy(fldi.vz));
		output_status( iostream );
		fprintf(iostream,"STATUS command end. Resuming execution.\n");
		
		close_interface_io( &iostream );
	}
	
	// OUTPUT command
	command_file = fopen(OUTPUT_COMMAND,"r");
	if(command_file) {
		fclose(command_file);
		// We have a status command
		remove( OUTPUT_COMMAND );
		open_interface_io( &iostream );
		fprintf(iostream,"OUTPUT command called. Calling for an immediate output\n");
		output_immediate(t);
		fprintf(iostream,"OUTPUT command end. Resuming execution\n");
		
		close_interface_io( &iostream );
	}
	
	// DUMP command
	command_file = fopen(DUMP_COMMAND,"r");
	if(command_file) {
		fclose(command_file);
		// We have a dump command
		remove( DUMP_COMMAND );
		open_interface_io( &iostream );
		fprintf(iostream,"DUMP command called. Calling for an immediate dump file\n");
		dump_immediate(t);
		fprintf(iostream,"DUMP command end. Resuming execution\n");
		
		close_interface_io( &iostream );
	}

	// STOP command
	command_file = fopen(STOP_COMMAND,"r");
	if(command_file) {
		fclose(command_file);
		// We have a status command
		remove( STOP_COMMAND );
		open_interface_io( &iostream );
		fprintf(iostream,"STOP command called. Terminating\n");
		output_immediate(t);
		finish_mainloop();
		finish_common();
		// Not yet coded 
		// finish_output();
		fprintf(iostream,"Goodbye\n");
		close_interface_io( &iostream );
		exit(0);
	}
	return;
}
	
	