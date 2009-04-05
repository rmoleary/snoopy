#include <stdlib.h>

#include "common.h"
#include "timestep.h"
#include "gfft.h"
#include "shear.h"

#define	OUTPUT_DUMP					"dump.dmp"			/**< Dump files filename. */
#define	OUTPUT_DUMP_VERSION			04					/**< Version of the dump files read and written by this code. */

#define	DUMP_MARKER					1981				/**< Marker used to signify the end of a dump file (it is also an excellent year...).*/

int noutput_flow;										/**< Next snapshot output number */
PRECISION lastoutput_time;								/**< Time when the last timevar output was done */
PRECISION lastoutput_flow;								/**< Time when the las snapshot output was done */
PRECISION lastoutput_dump;								/**< Time when the last dump output was done */

#ifdef WITH_SHEAR
PRECISION complex		*w1d, *w2d;						/** 1D arrays used by the remap methods */

fftw_plan	fft_1d_forward;								/**< 1D FFT transforms. Used by remap routines.*/
fftw_plan	fft_1d_backward;							/**< 1D FFT transforms. Used by remap routines.*/
#endif

#ifdef WITH_SHEAR
/***********************************************************/
/** 
	Remap a real field from the current sheared frame to the classical 
	cartesian frame. It could be more optimized, but since it is used
	in an output routine, I don't think such an optimization is worth doing.
	
	@param wri Real array to be remapped. The remapped array will be stored here.
	@param t Current time of the simulation (used to compute what the sheared frame is).
	
*/
/***********************************************************/

void remap_output(	PRECISION wri[], 
					const PRECISION t) {
					
	int i,j,k;
	PRECISION tvelocity;
	PRECISION tremap;
	complex PRECISION wexp;
	complex PRECISION phase;
	
#ifdef TIME_DEPENDANT_SHEAR
	tremap = time_shift(t);
	tvelocity = 0.0;
#else	
	tremap = time_shift(t);
	tvelocity = fmod(t, 2.0 * LY / (SHEAR * LX));
#endif	
	
	
	for( i = 0 ; i < NX / NPROC ; i++) {
		for( k = 0 ; k < NZ ; k++) {
			for( j = 0 ; j < NY ; j++) {
				w1d[j] = wri[k + j * (NZ + 2) + NY * (NZ + 2) * i];
			}
		
			// Transform w1d, which will be stored in w2d
			fftw_execute(fft_1d_forward);
					
			for( j = 0 ; j < NY ; j++) {
			// advection phase = ky*
#ifdef TIME_DEPENDANT_SHEAR
				// There is no proper remap in this case
				phase = (PRECISION complex) ((2.0 * M_PI) / LY * (fmod( j + (NY / 2) ,  NY ) - NY / 2 ) * 
											( ((double) (i + rank * (NX/NPROC)) / (double) NX - 0.5 ) * tremap ) * LX * SHEAR );
#else
				phase = (PRECISION complex) ((2.0 * M_PI) / LY * (fmod( j + (NY / 2) ,  NY ) - NY / 2 ) * 
											( ((double) (i + rank * (NX/NPROC)) / (double) NX ) * tremap - tvelocity / 2.0 ) * LX * SHEAR);
#endif
				wexp = cexp( I * phase);
									
				w2d[ j ] = w2d[ j ] * wexp;
			}
			
			fftw_execute(fft_1d_backward);
		
			for( j = 0 ; j < NY ; j++) {
				wri[k + j * (NZ + 2) + NY * (NZ + 2) * i] = w1d[j] / NY;
			}
		}
	}
	return;
}
#endif


/***********************************************************/
/** 
	Output a plain binary file (.raw file) from a complex array.
	The option FORTRAN_OUTPUT_ORDER can be used to emulate the output
	from a fortran code (row major order).
	This routine (contrarily to vtk routines) don't need extra memory
	for outputs. It should therefore been considered when memory is an issue
	for large resolution runs.
	
	@param t Current time of the simulation
	@param filename Filename of the output file
	@param wi Complex array containing the field to be written. This array is transformed into real space
	and remapped (if SHEAR is present) before being written
	
	\bug the FORTRAN_OUTPUT_ORDER option doesn't work when MPI is used
	
*/
/***********************************************************/
void write_snap(const PRECISION t, const char filename[], const PRECISION complex wi[]) {
	// Write the complex field wi[] in file filename. We need t for the remap thing.
	FILE *ht;
	int i,j,k;
#ifdef MPI_SUPPORT
	int current_rank;
	MPI_Status status;
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = wi[i];
	}
	
	gfft_c2r(w1);
	
	for( i = 0 ; i < 2 * NTOTAL_COMPLEX ; i++) {
		wr1[i] = wr1[i] / ((double) NTOTAL );
	}
	
#ifdef WITH_SHEAR
	remap_output(wr1,t);
#endif
	
#ifdef MPI_SUPPORT
	if(rank==0) {
#endif
	
	ht=fopen(filename,"w");
	
#ifdef MPI_SUPPORT
		for(current_rank=0; current_rank<NPROC; current_rank++) {
			if(current_rank!=0) {
				// We're not writing the array of the local_process...
				MPI_Recv( wr1, NTOTAL_COMPLEX * 2, MPI_DOUBLE, current_rank, 1, MPI_COMM_WORLD, &status);
			}
#endif	

#ifdef FORTRAN_OUTPUT_ORDER
	for( k = 0 ; k < NZ; k++) {
		for( j = 0; j < NY; j++) {
			for( i = 0; i < NX/NPROC; i++) {
#else
	for( i = 0; i < NX/NPROC; i++) {
		for( j = 0; j < NY; j++) {
			for( k = 0 ; k < NZ; k++) {
#endif
				fwrite(&wr1[k + j * (NZ + 2) + i * NY * (NZ + 2) ], sizeof(PRECISION), 1, ht);
			}
		}
	}
#ifdef MPI_SUPPORT
	}	// Close the for loop
#endif
	fclose(ht);

#ifdef MPI_SUPPORT
		// Wait for synchronization once everything is done...
		MPI_Barrier(MPI_COMM_WORLD);
	}  // Close the if on the rank
	else {
		// Send my array to process 0, and then synchronize
		MPI_Send( wr1, NTOTAL_COMPLEX * 2, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
	}
#endif
		

}

#ifdef VTK_OUTPUT
// VTK output using the visit_writer code.
/*************************************************
** VTK For HD*************************************
**************************************************/

/* ****************************************************************************/
/** Determines if the machine is little-endian.  If so, 
    it will force the data to be big-endian. 
	@param in_number floating point number to be converted in big endian */
/* *************************************************************************** */

float big_endian(float in_number)
{
    static int doneTest = 0;
    static int shouldSwap = 0;
	
    if (!doneTest)
    {
        int tmp1 = 1;
        unsigned char *tmp2 = (unsigned char *) &tmp1;
        if (*tmp2 != 0)
            shouldSwap = 1;
        doneTest = 1;
    }

    if (shouldSwap)
    {
		unsigned char *bytes = (unsigned char*) &in_number;
        unsigned char tmp = bytes[0];
        bytes[0] = bytes[3];
        bytes[3] = tmp;
        tmp = bytes[1];
        bytes[1] = bytes[2];
        bytes[2] = tmp;
    }
	return(in_number);
}

/***********************************************************/
/** 
	Output a floating point big endian array from a complex array.
	Use Forran output format (required by VTK). This routine is essentially
	useful for VTK files (hence its name...).
 	
	@param ht File handler in which the data has to be written
	@param wi Complex array containing the field to be written. This array is transformed into real space
	and remapped (if SHEAR is present) before being written.
	@param t Current time of the simulation.
	
*/
/***********************************************************/
void write_vtk(FILE * ht, PRECISION complex wi[], const PRECISION t) {
	// Write the data in the file handler *ht
	int i,j,k;
	float q0;

#ifdef MPI_SUPPORT
	PRECISION * chunk = NULL;
	if(rank==0) {
		chunk = (PRECISION *) malloc( NX * sizeof(PRECISION));
		if (chunk == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for chunk allocation");
	}
#endif
		
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = wi[i];
	}
	
	gfft_c2r(w1);
	
	for( i = 0 ; i < 2 * NTOTAL_COMPLEX ; i++) {
		wr1[i] = wr1[i] / ((double) NTOTAL );
	}
	
#ifdef WITH_SHEAR
	remap_output(wr1,t);
#endif

	for( k = 0 ; k < NZ; k++) {
		for( j = 0; j < NY; j++) {
#ifdef MPI_SUPPORT
			// We have to transpose manually to be Fortran compliant
			for(i = 0; i < NX/NPROC; i++) {
				wr2[i] = wr1[k + j * (NZ + 2) + i * NY * (NZ + 2) ];   // Transfer the chunk of data to wr2
			}
			MPI_Gather(wr2, NX/NPROC, MPI_DOUBLE,
					   chunk, NX/NPROC, MPI_DOUBLE, 0, MPI_COMM_WORLD); // Put the full chunk in chunk in the root process
#endif			
			for( i = 0; i < NX; i++) {
#ifdef MPI_SUPPORT
				if(rank==0) {
					q0 = big_endian( (float) chunk[ i ] );
					fwrite(&q0, sizeof(float), 1, ht);
				}
#else
				q0 = big_endian( (float) wr1[k + j * (NZ + 2) + i * NY * (NZ + 2) ] );
				fwrite(&q0, sizeof(float), 1, ht);
#endif
			}
		}
	}
#ifdef MPI_SUPPORT	
	if(rank==0) free(chunk);
#endif
	
	return;
}

// Geo's personnal VTK writer, using structured data points	
/***********************************************************/
/** 
	Output a legacy VTK file readable by Paraview. This routine
	will output all the variables in files data/v****.vtk.
	
	@param n Number of the file in which the output will done.
	@param t Current time of the simulation.
*/
/***********************************************************/

void output_vtk(const int n, PRECISION t) {
	FILE *ht = NULL;
	char  filename[50];
	int num_remain_field;
	
	sprintf(filename,"data/v%04i.vtk",n);
	if(rank==0) {
		ht=fopen(filename,"w");
	
		fprintf(ht, "# vtk DataFile Version 2.0\n");
		fprintf(ht, "Created by the Snoopy Code\n");
		fprintf(ht, "BINARY\n");
		fprintf(ht, "DATASET STRUCTURED_POINTS\n");
		fprintf(ht, "DIMENSIONS %d %d %d\n", NX, NY, NZ);
		fprintf(ht, "ORIGIN %g %g %g\n", -LX/2.0, -LY/2.0, -LZ/2.0);
		fprintf(ht, "SPACING %g %g %g\n", LX/NX, LY/NY, LZ/NZ);
	
		// Write the primary scalar (f***ing VTK legacy format...)
		fprintf(ht, "POINT_DATA %d\n",NX*NY*NZ);
		fprintf(ht, "SCALARS vx float\n");
		fprintf(ht, "LOOKUP_TABLE default\n");
	}
	write_vtk(ht,fld.vx,t);

	num_remain_field = 2;		// Number of remaining field to be written
#ifdef MHD
	num_remain_field += 3;
#endif
#ifdef BOUSSINESQ
	num_remain_field += 1;
#endif
	
	if(rank==0) fprintf(ht, "FIELD FieldData %d\n",num_remain_field);
	
	// Write all the remaining fields
	
	if(rank==0) fprintf(ht, "vy 1 %d float\n",NX*NY*NZ);
	write_vtk(ht,fld.vy,t);
	
	if(rank==0) fprintf(ht, "vz 1 %d float\n",NX*NY*NZ);
	write_vtk(ht,fld.vz,t);

#ifdef MHD
	if(rank==0) fprintf(ht, "bx 1 %d float\n",NX*NY*NZ);
	write_vtk(ht,fld.bx,t);

	if(rank==0) fprintf(ht, "by 1 %d float\n",NX*NY*NZ);
	write_vtk(ht,fld.by,t);

	if(rank==0) fprintf(ht, "bz 1 %d float\n",NX*NY*NZ);
	write_vtk(ht,fld.bz,t);
#endif
#ifdef BOUSSINESQ
	if(rank==0) fprintf(ht, "th 1 %d float\n",NX*NY*NZ);
	write_vtk(ht,fld.th,t);
#endif	  
	if(rank==0) fclose(ht);
	return;
	
}
#endif
	
/***********************************************************/
/** 
	Output plain binary snaphsost (raw files) readable by any
	descent software. This routine calls the relevant
	write subroutine for each field.
	@param n Number of the file in which the output will done.
	@param t Current time of the simulation
*/
/***********************************************************/
	
void output_flow(const int n, const PRECISION t) {
	
	char filename[50];
	
	sprintf(filename,"data/vx%04i.raw",n);
	write_snap(t, filename, fld.vx);
	
	sprintf(filename,"data/vy%04i.raw",n);
	write_snap(t, filename, fld.vy);
	
	sprintf(filename,"data/vz%04i.raw",n);
	write_snap(t, filename, fld.vz);
	
#ifdef BOUSSINESQ
	sprintf(filename,"data/th%04i.raw",n);
	write_snap(t, filename, fld.th);
#endif

#ifdef MHD
	sprintf(filename,"data/bx%04i.raw",n);
	write_snap(t, filename, fld.bx);
	
	sprintf(filename,"data/by%04i.raw",n);
	write_snap(t, filename, fld.by);
	
	sprintf(filename,"data/bz%04i.raw",n);
	write_snap(t, filename, fld.bz);
#endif
	return;
}

/***********************************************************/
/** 
	Write statistical quantities using text format in the file
	timevar. 
	@param fldi Field structure from which the statistical quantities are derived.
	@param t Current time of the simulation
*/
/***********************************************************/

void output_timevar(const struct Field fldi,
					const PRECISION t) {
					
	FILE *ht;
	PRECISION vxmax, vxmin, vymax, vymin, vzmax, vzmin, thmin, thmax;
	PRECISION bxmax, bxmin, bymax, bymin, bzmax, bzmin;
	PRECISION energy_total;	
	PRECISION energy_mag;
	PRECISION reynolds_stress;
	PRECISION maxwell_stress;

	int i;
	
		
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = fldi.vx[i];
		w2[i] = fldi.vy[i];
		w3[i] = fldi.vz[i];
#ifdef BOUSSINESQ
		w4[i] = fldi.th[i];
#endif
#ifdef MHD
		w5[i] = fldi.bx[i];
		w6[i] = fldi.by[i];
		w7[i] = fldi.bz[i];
#endif
	}

	energy_total = energy(w1) + energy(w2) + energy(w3);
	reduce(&energy_total,1);
#ifdef MHD
	energy_mag = energy(w5) + energy(w6) + energy(w7);
	reduce(&energy_mag,1);
#else
	energy_mag=0.0;
#endif
	
	gfft_c2r(w1);
	gfft_c2r(w2);
	gfft_c2r(w3);
#ifdef BOUSSINESQ
	gfft_c2r(w4);
#endif
#ifdef MHD
	gfft_c2r(w5);
	gfft_c2r(w6);
	gfft_c2r(w7);
#endif

	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr1[i] = wr1[i] / ((double) NTOTAL );
		wr2[i] = wr2[i] / ((double) NTOTAL );
		wr3[i] = wr3[i] / ((double) NTOTAL );
#ifdef BOUSSINESQ
		wr4[i] = wr4[i] / ((double) NTOTAL );
#endif
#ifdef MHD
		wr5[i] = wr5[i] / ((double) NTOTAL );
		wr6[i] = wr6[i] / ((double) NTOTAL );
		wr7[i] = wr7[i] / ((double) NTOTAL );
#endif
		
	}
	
	// w1, w2, w3 contains vx vy vz
	// w4 contains the entropy perturbations
		
	vxmax=0.0;
	vxmin=0.0;
	vymax=0.0;
	vymin=0.0;
	vzmax=0.0;
	vzmin=0.0;
	thmax=0.0;
	thmin=0.0;
	reynolds_stress=0.0;
	bxmax=0.0;
	bxmin=0.0;
	bymax=0.0;
	bymin=0.0;
	bzmax=0.0;
	bzmin=0.0;
	maxwell_stress=0.0;

	
		
	for(i = 0 ; i < 2 * NTOTAL_COMPLEX ; i ++) {
		if(vxmax < wr1[i]) vxmax = wr1[i];
		if(vxmin > wr1[i]) vxmin = wr1[i];
		if(vymax < wr2[i]) vymax = wr2[i];
		if(vymin > wr2[i]) vymin = wr2[i];
		if(vzmax < wr3[i]) vzmax = wr3[i];
		if(vzmin > wr3[i]) vzmin = wr3[i];
		
		reynolds_stress += wr1[i] * wr2[i] / ((double) NTOTAL);
#ifdef BOUSSINESQ
		if(thmax < wr4[i]) thmax = wr4[i];
		if(thmin > wr4[i]) thmin = wr4[i];
#endif
#ifdef MHD
		if(bxmax < wr5[i]) bxmax = wr5[i];
		if(bxmin > wr5[i]) bxmin = wr5[i];
		if(bymax < wr6[i]) bymax = wr6[i];
		if(bymin > wr6[i]) bymin = wr6[i];
		if(bzmax < wr7[i]) bzmax = wr7[i];
		if(bzmin > wr7[i]) bzmin = wr7[i];
		
		maxwell_stress += wr5[i] * wr6[i] / ((double) NTOTAL);
#endif
	}
	
	reduce(&vxmax,2);
	reduce(&vxmin,3);
	reduce(&vymax,2);
	reduce(&vymin,3);
	reduce(&vzmax,2);
	reduce(&vzmin,3);
	
	reduce(&reynolds_stress,1);
#ifdef BOUSSINESQ
	reduce(&thmax,2);
	reduce(&thmin,3);
#endif
#ifdef MHD
	reduce(&bxmax,2);
	reduce(&bxmin,3);
	reduce(&bymax,2);
	reduce(&bymin,3);
	reduce(&bzmax,2);
	reduce(&bzmin,3);
	
	reduce(&maxwell_stress,1);
#endif

	
	if(rank==0) {
		ht=fopen("timevar","a");
		fprintf(ht,"%08e\t",t);
		fprintf(ht,"%08e\t%08e\t",energy_total,energy_mag);
		fprintf(ht,"%08e\t%08e\t%08e\t%08e\t%08e\t%08e\t",vxmax,vxmin,vymax,vymin,vzmax,vzmin);
		fprintf(ht,"%08e\t",reynolds_stress);
		fprintf(ht,"%08e\t%08e\t%08e\t%08e\t%08e\t%08e\t",bxmax,bxmin,bymax,bymin,bzmax,bzmin);
		fprintf(ht,"%08e\t",maxwell_stress);
		fprintf(ht,"%08e\t%08e",thmax,thmin);
		fprintf(ht,"\n");
		fclose(ht);
#ifdef MPI_SUPPORT
		MPI_Barrier(MPI_COMM_WORLD);
	}
	else	MPI_Barrier(MPI_COMM_WORLD);
#else
	}
#endif
	return;
}
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
void write_field(FILE *handler, PRECISION complex *fldwrite) {
#ifdef MPI_SUPPORT
	MPI_Status status;
	// Write in rank order using the file opened if handler
	int current_rank;
#endif
	int i;

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
			fwrite(w1, sizeof(PRECISION complex), NTOTAL_COMPLEX, handler);
#ifdef MPI_SUPPORT
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	else {
		MPI_Send(fldwrite, NTOTAL_COMPLEX * 2, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
	}
#endif
	return;
}

/**************************************************************************************/
/** 
	read a field information from a restart dump file, taking care of the MPI broadcast
	@param handler handler to an already opened restart dump filed
	@param fldread Field structure in which the restart information will be stored
*/
/**************************************************************************************/
	
void read_field(FILE *handler, PRECISION complex *fldread) {
#ifdef MPI_SUPPORT
	MPI_Status status;
	// Write in rank order using the file opened if handler
	int current_rank;
#endif
	int i;

#ifdef MPI_SUPPORT
	if(rank==0) {
		for(current_rank=0; current_rank < NPROC; current_rank++) {
#endif
			fread(w1, sizeof(PRECISION complex), NTOTAL_COMPLEX, handler);

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
				  const PRECISION t) {	
				  
	FILE *ht;
	int dump_version;
	int size_x,	size_y, size_z;
	int marker, included_field;
	
	ht=NULL;
	
	size_x = NX;
	size_y = NY;
	size_z = NZ;
	
	// This is a check when we try to read a dump file
	dump_version = OUTPUT_DUMP_VERSION;
	
	// This is a hard coded marker to check that we have read correctly the file
	marker = DUMP_MARKER;
	
	if(rank==0) {
		ht=fopen(OUTPUT_DUMP,"w");
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

	if(rank==0) {
		fwrite(&t			, sizeof(PRECISION)		   , 1			   , ht);
	
		fwrite(&noutput_flow		, sizeof(int)			   , 1             , ht);
		fwrite(&lastoutput_time 	, sizeof(PRECISION)		   , 1			   , ht);
		fwrite(&lastoutput_flow 	, sizeof(PRECISION)		   , 1			   , ht);
		fwrite(&lastoutput_dump 	, sizeof(PRECISION)		   , 1			   , ht);
	
// Any extra information should be put here.	
	
		fwrite(&marker		, sizeof(int)			   , 1			   , ht);
	
		fclose(ht);
	}
	
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
				  PRECISION *t) {	
				  
	FILE *ht;
	int dump_version;
	int size_x,	size_y, size_z, included_field;
	int marker;
	
	ht=NULL;
	
	if(rank==0) {
		ht=fopen(OUTPUT_DUMP,"r");
		if(ht==NULL) ERROR_HANDLER( ERROR_CRITICAL, "Error opening dump file.");
	
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
	
	if(rank==0) {
		fread(t			, sizeof(PRECISION)		   , 1			   , ht);
	
		fread(&noutput_flow			, sizeof(int)			   , 1             , ht);
		fread(&lastoutput_time		, sizeof(PRECISION)		   , 1			   , ht);
		fread(&lastoutput_flow		, sizeof(PRECISION)		   , 1			   , ht);
		fread(&lastoutput_dump		, sizeof(PRECISION)		   , 1			   , ht);
	
		fread(&marker , sizeof(int)			   , 1, ht);
	
		if(marker != DUMP_MARKER) ERROR_HANDLER( ERROR_CRITICAL, "Incorrect marker. Probably an incorrect dump file!");	
		fclose(ht);
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
	return;
}
	
/*********************************************************
*** General routine, callable from outside ***************
**********************************************************/
/**************************************************************************************/
/** 
	Initialize the output variables. Should be called only in the begining .
*/
/**************************************************************************************/

void init_output() {

#ifdef WITH_SHEAR
// Initialize 1D arrays for remaps
	w1d = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NY);
	if (w1d == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w1d allocation");
	
	w2d = (PRECISION complex *) fftw_malloc( sizeof(PRECISION complex) * NY);
	if (w2d == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w2d allocation");

// FFT plans (we use dummy arrays since we use the "guru" interface of fft3 in the code)
// The in place/ out of place will be set automatically at this stage

#ifdef OPENMP_SUPPORT	
	fftw_plan_with_nthreads( 1 );
#endif

	fft_1d_forward = fftw_plan_dft_1d(NY, w1d, w2d, FFTW_FORWARD, FFT_PLANNING);
	if (fft_1d_forward == NULL) ERROR_HANDLER( ERROR_CRITICAL, "FFTW 1D forward plan creation failed");
	
	fft_1d_backward = fftw_plan_dft_1d(NY, w2d, w1d, FFTW_BACKWARD, FFT_PLANNING);
	if (fft_1d_backward == NULL) ERROR_HANDLER( ERROR_CRITICAL, "FFTW 1D backward plan creation failed");
#endif	

	
#ifndef RESTART
	noutput_flow=0;
	lastoutput_time = T_INITIAL - TOUTPUT_TIME;
	lastoutput_flow = T_INITIAL - TOUTPUT_FLOW;
	lastoutput_dump = T_INITIAL - TOUTPUT_DUMP;
	
#endif
	return;
}

/**************************************************************************************/
/** 
	Check if an output (timevar, snapshot or dump) is required at t. If yes, call the 
	relevant routines.
	@param t Current time in the simulation
*/
/**************************************************************************************/

void output(const PRECISION t) {
	// Very rough output function
	if( (t-lastoutput_time)>=TOUTPUT_TIME) {
		output_timevar(fld,t);
		lastoutput_time = lastoutput_time + TOUTPUT_TIME;
	}
	
	if( (t-lastoutput_flow)>=TOUTPUT_FLOW) {
#ifdef VTK_OUTPUT
		output_vtk(noutput_flow,t);
#else
		output_flow(noutput_flow,t);
#endif
		noutput_flow++;
		lastoutput_flow = lastoutput_flow + TOUTPUT_FLOW;
	}
		
	if( (t-lastoutput_dump)>=TOUTPUT_DUMP) {
		lastoutput_dump=lastoutput_dump+TOUTPUT_DUMP;
		output_dump(fld,t);
	}
	return;
}

/**************************************************************************************/
/** 
	Show the current status of the output routine.
	@param iostream Handler of the file in which the status is written.
*/
/**************************************************************************************/
void output_status(FILE * iostream) {
	if(rank==0)
		fprintf(iostream,"Next output in file n %d, at t=%e\n",noutput_flow, lastoutput_flow+TOUTPUT_FLOW);
	return;
}

/**************************************************************************************/
/** 
	Immediatly output a timevar, snapshot and dump, regardless of the output
	parameters.
	@param t Current time of the simulation.
*/
/**************************************************************************************/

void output_immediate(const PRECISION t) {
	// Very rough output function
	// Immediate output
	output_timevar(fld,t);
#ifdef VTK_OUTPUT
	output_vtk(noutput_flow,t);
#else
	output_flow(noutput_flow,t);
#endif
	noutput_flow++;
	output_dump(fld,t);
	return;
}

/**************************************************************************************/
/** 
	Immediatly output a dump file, regardless of the output
	parameters.
	@param t Current time of the simulation.
*/
/**************************************************************************************/

void dump_immediate(const PRECISION t) {
	output_dump(fld,t);
	return;
}
/**************************************************************************************/
/** 
	Remove the timevar file (if exists) to start from a fresh one.
*/
/**************************************************************************************/
void clear_timevar() {
	FILE* ht;
	if(rank==0) {
	ht=fopen("timevar","w");
	fclose(ht);
	}
}
	
	