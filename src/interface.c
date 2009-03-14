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

int check_file(char filename[]) {
	FILE * command_file;
	int is_present=0;
	if(rank==0) { 
		command_file = fopen(filename,"r");
		if(command_file) {
			is_present=1;
			fclose(command_file);
			remove(filename); 
		}
		else
			is_present=0;
	}
#ifdef MPI_SUPPORT
	MPI_Bcast( &is_present, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
	return(is_present);
}
	
void check_interface(const struct Field fldi,
					 const PRECISION t,
					 const PRECISION dt,
					 const int		 nloop,
					 const double	tstart) {
	// This routine check the interface file and print the relevant informations
	
	FILE * iostream = NULL;
	
	// STATUS command
	if(check_file(STATUS_COMMAND)) {
		// We have a status command
		if(rank==0) {
			open_interface_io( &iostream );
		
			fprintf(iostream,"STATUS command called.\n");
			fprintf(iostream,"t=%e, dt=%e, nloop=%d, sec/loop=%f\n",t,dt,nloop, (get_c_time()-tstart)/nloop);
		}
		output_status( iostream );
		if(rank==0) {
			fprintf(iostream,"STATUS command end. Resuming execution.\n");
			close_interface_io( &iostream );
		}
	}
	
	// OUTPUT command
	if(check_file(OUTPUT_COMMAND)) {
		// We have a status command
		if(rank==0) {
			open_interface_io( &iostream );
			fprintf(iostream,"OUTPUT command called. Calling for an immediate output\n");
		}
		output_immediate(t);
		if(rank==0) {
			fprintf(iostream,"OUTPUT command end. Resuming execution\n");
			close_interface_io( &iostream );
		}
	}
	
	// DUMP command
	if(check_file(DUMP_COMMAND)) {
		// We have a dump command
		if(rank==0) {
			open_interface_io( &iostream );
			fprintf(iostream,"DUMP command called. Calling for an immediate dump file\n");
		}
		dump_immediate(t);
		if(rank==0) {
			fprintf(iostream,"DUMP command end. Resuming execution\n");
			close_interface_io( &iostream );
		}
	}

	// STOP command
	if(check_file(STOP_COMMAND)) {
		// We have a status command
		if(rank==0) {
			open_interface_io( &iostream );
			fprintf(iostream,"STOP command called. Calling immediate dump and terminating\n");
		}
		dump_immediate(t);
		finish_mainloop();
		finish_common();
		// Not yet coded 
		// finish_output();
		if(rank==0) {
			fprintf(iostream,"Goodbye\n");
			close_interface_io( &iostream );
		}
#ifdef MPI_SUPPORT
		MPI_Finalize();
#endif
		exit(0);
	}
	return;
}
	
	