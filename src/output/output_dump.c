#include <stdlib.h>

#include "../common.h"
#include "../gfft.h"
#include "../shear.h"
#include "../debug.h"
#include "../particles.h"

#define	OUTPUT_DUMP_VERSION			04					/**< Version of the dump files read and written by this code. */
#define	DUMP_MARKER					1981				/**< Marker used to signify the end of a dump file (it is also an excellent year...).*/

// This is nasty, but unfortunately there is no proper way to do that since
// The dump is a GLOBAL image of the system.
// These variables are declared in output.c

extern int noutput_flow;
extern double lastoutput_time;								/**< Time when the last timevar output was done */
extern double lastoutput_flow;								/**< Time when the las snapshot output was done */
extern double lastoutput_dump;								/**< Time when the last dump output was done */

/**********************************************************
** Restart DUMP I/O routines ******************************
***********************************************************/



/***********************************************************/
/** 
	write a field information to a restart dump file, taking care of the MPI reduction
	@param handler handler to an already opened restart dump filed
	@param fldwrite Field structure pointing to the field needed to be saved
*/
/***********************************************************/
void write_field(FILE *handler, double complex *fldwrite) {
#ifdef MPI_SUPPORT
	MPI_Status status;
	// Write in rank order using the file opened if handler
	int current_rank;
#endif
	int i;

	DEBUG_START_FUNC;
	
#ifdef MPI_SUPPORT	
	if(rank==0) {
		for(current_rank=0; current_rank < NPROC; current_rank++) {
			if(current_rank==0) {
				// Copy the dump in the rank=0 process
#endif
				for(i=0; i< NTOTAL_COMPLEX; i++) {
					w1[i]=fldwrite[i];
				}
#ifdef MPI_SUPPORT
			}
			else {
				MPI_Recv( w1, NTOTAL_COMPLEX * 2, MPI_DOUBLE, current_rank, 2, MPI_COMM_WORLD, &status);
			}
#endif
			fwrite(w1, sizeof(double complex), NTOTAL_COMPLEX, handler);
			
			if(ferror(handler)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing dump file");
			
#ifdef MPI_SUPPORT
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	else {
		MPI_Send(fldwrite, NTOTAL_COMPLEX * 2, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
	}
#endif

	DEBUG_END_FUNC;
	
	return;
}

/**************************************************************************************/
/** 
	read a field information from a restart dump file, taking care of the MPI broadcast
	@param handler handler to an already opened restart dump filed
	@param fldread Field structure in which the restart information will be stored
*/
/**************************************************************************************/
	
void read_field(FILE *handler, double complex *fldread) {
#ifdef MPI_SUPPORT
	MPI_Status status;
	// Write in rank order using the file opened if handler
	int current_rank;
#endif
	int i;

	DEBUG_START_FUNC;
	
#ifdef MPI_SUPPORT
	if(rank==0) {
		for(current_rank=0; current_rank < NPROC; current_rank++) {
#endif
			fread(w1, sizeof(double complex), NTOTAL_COMPLEX, handler);
			if(ferror(handler)) ERROR_HANDLER( ERROR_CRITICAL, "Error reading dump file");

#ifdef MPI_SUPPORT
			if(current_rank==0) {
#endif
				// Copy the dump in the rank=0 process
				for(i=0; i< NTOTAL_COMPLEX; i++) {
					fldread[i]=w1[i];
				}
#ifdef MPI_SUPPORT
			}
			else {
				MPI_Send( w1, NTOTAL_COMPLEX * 2, MPI_DOUBLE, current_rank, 3, MPI_COMM_WORLD);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	else {
		MPI_Recv(fldread, NTOTAL_COMPLEX * 2, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD,&status);
		MPI_Barrier(MPI_COMM_WORLD);
	}
#endif

	DEBUG_END_FUNC;
	
	return;
}


/**************************************************************************************/
/** 
	write a full restart dump
	@param fldi Field structure pointing to the field needed to be savec
	@param t current time of the simulation (will be stored in the dump file)
*/
/**************************************************************************************/
void output_dump( const struct Field fldi,
				  const double t) {	
				  
	FILE *ht;
	int dump_version;
	int size_x,	size_y, size_z;
	int marker, included_field;
	int nfield;
	long int filesize;
	
	ht=NULL;
	
	DEBUG_START_FUNC;
	
	size_x = NX;
	size_y = NY;
	size_z = NZ;
	
	// This is a check when we try to read a dump file
	dump_version = OUTPUT_DUMP_VERSION;
	
	// This is a hard coded marker to check that we have read correctly the file
	marker = DUMP_MARKER;
	
	if(rank==0) {
		ht=fopen(OUTPUT_DUMP_WRITE,"w");
		if(ht==NULL) ERROR_HANDLER( ERROR_CRITICAL, "Error opening dump file.");
	
		fwrite(&dump_version, sizeof(int), 1, ht);
	
		fwrite(&size_x		, sizeof(int), 1, ht);
		fwrite(&size_y		, sizeof(int), 1, ht);
		fwrite(&size_z		, sizeof(int), 1, ht);
		// Included fields
		// First bit is Boussinesq fields
		// Second bit is MHD fields
		// Other fields can be added from that stage...
		included_field=0;
#ifdef BOUSSINESQ
		included_field+=1;
#endif
#ifdef MHD
		included_field+=2;
#endif
#ifdef WITH_PARTICLES
		included_field+=4;
#endif
#ifdef COMPRESSIBLE
		included_field+=8;
#endif
		fwrite(&included_field, sizeof(int), 1, ht);
	}
	
	write_field(ht, fldi.vx);
	write_field(ht, fldi.vy);
	write_field(ht, fldi.vz);
	
#ifdef BOUSSINESQ
	write_field(ht, fldi.th);
#endif
#ifdef MHD
	write_field(ht, fldi.bx);
	write_field(ht, fldi.by);
	write_field(ht, fldi.bz);
#endif
#ifdef WITH_PARTICLES
	write_particle_dump(ht, fldi.part);
#endif
#ifdef COMPRESSIBLE
	write_field(ht, fldi.d);
#endif

	if(rank==0) {
		fwrite(&t			, sizeof(double)		   , 1			   , ht);
	
		fwrite(&noutput_flow		, sizeof(int)			   , 1             , ht);
		fwrite(&lastoutput_time 	, sizeof(double)		   , 1			   , ht);
		fwrite(&lastoutput_flow 	, sizeof(double)		   , 1			   , ht);
		fwrite(&lastoutput_dump 	, sizeof(double)		   , 1			   , ht);
	
// Any extra information should be put here.	
	
		fwrite(&marker		, sizeof(int)			   , 1			   , ht);
	
// Check everything was fine with the file
		
		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing dump file");
		
		fclose(ht);
	}
	
	// predict the file size:
	nfield = 3;		// Velocity fields
#ifdef BOUSSINESQ
	nfield += 1;
#endif
#ifdef MHD
	nfield += 3;
#endif
#ifdef COMPRESSIBLE
	nfield += 1;
#endif

	filesize = nfield * sizeof(double complex) * NTOTAL_COMPLEX * NPROC + 4 * sizeof( double) + 7 * sizeof( int );

#ifdef WITH_PARTICLES
	filesize += param.particles_n * sizeof(struct Particle) + sizeof(int);
#endif
	
	if( check_file_size( OUTPUT_DUMP_WRITE, filesize ) ) {
		MPI_Printf("Error checking the dump size, got %d instead of %d\n", (int) check_file_size( OUTPUT_DUMP_WRITE, filesize), (int) filesize);
		ERROR_HANDLER( ERROR_CRITICAL, "Error writing dump file, check your quotas");
	}
	
// This bit prevents the code from loosing all the dump files (this kind of thing happens sometimes...)
// With this routine, one will always have a valid restart dump, either in OUTPUT_DUMP_WRITE, OUTPUT_DUMP or OUTPUT_DUMP_SAV 
// (it should normally be in OUTPUT_DUMP)

	if(rank==0) {
		remove(OUTPUT_DUMP_SAV);				 // Delete the previously saved output dump
		rename(OUTPUT_DUMP, OUTPUT_DUMP_SAV);	 // Save the current dump file
		rename(OUTPUT_DUMP_WRITE, OUTPUT_DUMP);  // Move the new dump file to its final location
	}
	
	if( check_file_size( OUTPUT_DUMP, filesize ) ) {
		MPI_Printf("Error checking the dump size, got %d instead of %d\n", (int) check_file_size( OUTPUT_DUMP, filesize), (int) filesize);
		ERROR_HANDLER( ERROR_CRITICAL, "Error writing dump file, check your quotas");
	}
	
#ifdef MPI_SUPPORT
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	DEBUG_END_FUNC;
	
	return;
}

/**************************************************************************************/
/** 
	read a full restart dump
	@param fldo Field structure in which the restart dump will be stored.
	@param t time of the simulation, overwritten by this with the dump information.
*/
/**************************************************************************************/
void read_dump(   struct Field fldo,
				  double *t,
				  char dump_file[]) {	
				  
	FILE *ht;
	int dump_version;
	int size_x,	size_y, size_z, included_field;
	int marker;

	DEBUG_START_FUNC;
	
	if( !file_exist(dump_file) ) {
		// The file cannot be opened
		MPI_Printf("File %s not found\n",dump_file);
		ERROR_HANDLER(ERROR_CRITICAL, "Cannot open dump file.");
	}
	
	ht=fopen(dump_file,"r");
	// The file has been opened by process 0
	
	if(rank==0) {
		fread(&dump_version, sizeof(int), 1, ht);
		if( dump_version != OUTPUT_DUMP_VERSION) ERROR_HANDLER( ERROR_CRITICAL, "Incorrect dump file version.");
		
		fread(&size_x		, sizeof(int), 1, ht);
		fread(&size_y		, sizeof(int), 1, ht);
		fread(&size_z		, sizeof(int), 1, ht);
		
		if(size_x != NX) ERROR_HANDLER( ERROR_CRITICAL, "Incorrect X grid size in dump file.");
		if(size_y != NY) ERROR_HANDLER( ERROR_CRITICAL, "Incorrect Y grid size in dump file.");
		if(size_z != NZ) ERROR_HANDLER( ERROR_CRITICAL, "Incorrect Y grid size in dump file.");
		
		fread(&included_field, sizeof(int), 1, ht);
	}
#ifdef MPI_SUPPORT
	MPI_Bcast( &included_field,		1, MPI_INT,		0, MPI_COMM_WORLD);
#endif
	
	MPI_Printf("Reading velocity field\n");
	read_field(ht, fldo.vx);
	read_field(ht, fldo.vy);
	read_field(ht, fldo.vz);
	
#ifdef BOUSSINESQ
	// Do we have Boussinesq Field in the dump?
	if(included_field & 1) {
		// Yes
		MPI_Printf("Reading Boussinesq field\n");
		read_field(ht, fldo.th);
	}
	else {
		// No
		ERROR_HANDLER( ERROR_WARNING, "No Boussinesq field in the dump, using initial conditions.");
	}
#endif
#ifdef MHD
	// Do we have MHD field in the dump?
	if(included_field & 2) {
		// Yes
		MPI_Printf("Reading MHD field\n");
		read_field(ht, fldo.bx);
		read_field(ht, fldo.by);
		read_field(ht, fldo.bz);
	}
	else {
		//No
		ERROR_HANDLER( ERROR_WARNING, "No MHD field in the dump, using initial conditions.");
	}
#endif
#ifdef WITH_PARTICLES
	// Do we have particles data in the dump?
	if(included_field & 4) {
		// Yes
		MPI_Printf("Reading particles\n");
		read_particle_dump(ht, fldo.part);
	}
	else {
		// No
		ERROR_HANDLER( ERROR_WARNING, "No Particles in the dump, using initial conditions.");
	}
#endif
#ifdef COMPRESSIBLE
	// Do we have tide field in the dump?
	if(included_field & 8) {
		// Yes
		MPI_Printf("Reading density field\n");
		read_field(ht, fldo.d);
	}
	else {
		//No
		ERROR_HANDLER( ERROR_WARNING, "No tide field in the dump, using initial conditions.");
	}
#endif

	
	if(param.restart) {		// If the dump is used to restart, we need these extra variables
		if(rank==0) {
			fread(t			, sizeof(double)		   , 1			   , ht);
	
			fread(&noutput_flow			, sizeof(int)			   , 1             , ht);
			fread(&lastoutput_time		, sizeof(double)		   , 1			   , ht);
			fread(&lastoutput_flow		, sizeof(double)		   , 1			   , ht);
			fread(&lastoutput_dump		, sizeof(double)		   , 1			   , ht);
	
			fread(&marker , sizeof(int)			   , 1, ht);
	
			if(marker != DUMP_MARKER) ERROR_HANDLER( ERROR_CRITICAL, "Incorrect marker. Probably an incorrect dump file!");	
			
		}
	
	// Transmit the values to all processes
#ifdef MPI_SUPPORT
		MPI_Bcast( t,					1, MPI_DOUBLE,	0, MPI_COMM_WORLD);
		MPI_Bcast( &noutput_flow,		1, MPI_INT,		0, MPI_COMM_WORLD);
		MPI_Bcast( &lastoutput_time,	1, MPI_DOUBLE,	0, MPI_COMM_WORLD);
		MPI_Bcast( &lastoutput_flow,	1, MPI_DOUBLE,	0, MPI_COMM_WORLD);
		MPI_Bcast( &lastoutput_dump,	1, MPI_DOUBLE,	0, MPI_COMM_WORLD);
#endif
	
		MPI_Printf("Restarting at t=%e...\n",*t);
	}

	if(rank==0) {
		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error reading dump file");
		fclose(ht);
	}
	
	DEBUG_END_FUNC;
	
	return;
}
