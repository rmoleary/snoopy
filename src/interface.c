#include <stdlib.h>
#include <stdio.h>

#include "snoopy.h"
#include "common.h"
#include "mainloop.h"
#include "output.h"
#include "transpose.h"

// This is the user interface used in Snoopy
// We assume check_interface is called periodically to check whether the user called for something

// currently, we reconize commands status, output, dump and stop, as new files in the root directory of the code

#define			STATUS_COMMAND		"status"
#define			OUTPUT_COMMAND		"output"
#define			DUMP_COMMAND		"dump"
#define			STOP_COMMAND		"stop"

#define			INTERFACE_OUTPUT	"output_interface.txt"

void open_interface_io(FILE ** iostream) {
	if(param.interface_output_file)
		*iostream = fopen(INTERFACE_OUTPUT,"w");
	else
		*iostream = stdout;

	return;
}

void close_interface_io(FILE ** iostream) {
	if(param.interface_output_file)
		fclose(*iostream);
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
					 const double t,
					 const double dt,
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
#ifdef MPI_SUPPORT
			fprintf(iostream,"transpose time %f seconds, or %f pc of computation time\n",read_transpose_timer(), read_transpose_timer()/(get_c_time()-tstart)*100.0);
#endif
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
		finish_output();
		finish_gfft();
		finish_common();
		
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
	
	