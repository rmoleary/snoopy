#include <stdlib.h>

#include "../common.h"
#include "../gfft.h"
#include "../shear.h"
#include "../debug.h"
#include "../particles.h"

#ifdef WITH_SHEAR

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

void remap_output(	double wri[], 
					const double t) {
					
	int i,j,k;
	double tvelocity;
	double tremap;
	complex double wexp;
	complex double phase;
	double complex		*w2d;
	
	DEBUG_START_FUNC;
	
	w2d = (double complex *) fftw_malloc( sizeof(double complex) * (NY/2+1) * NZ );
	if (w2d == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w2d allocation");
	
#ifdef TIME_DEPENDANT_SHEAR
	tremap = time_shift(t);
	tvelocity = 0.0;
#else	
	tremap = time_shift(t);
	tvelocity = fmod(t, 2.0 * param.ly / (param.shear * param.lx));
#endif	
	
	for( i = 0 ; i < NX/NPROC ; i++) {
#ifdef WITH_2D
		fftw_execute_dft_r2c(fft_1d_forward, wri + i*(NY+2), w2d);
#else
		fftw_execute_dft_r2c(fft_1d_forward, wri + i*(NZ+2)*NY, w2d);
#endif
		for( j = 0 ; j < NY/2+1 ; j++) {
			phase = (double complex) ((2.0 * M_PI) / param.ly *  ((double) j )  * 
									( ((double) (i + rank * (NX/NPROC)) / (double) NX ) * tremap - tvelocity / 2.0 ) * param.lx * param.shear);
			
			//printf("phase=%g + I %g\n",creal(phase), cimag(phase));
			
			wexp = cexp( I * phase)/NY;
			
			//printf("wexp=%g + I %g\n",creal(wexp), cimag(wexp));

			for( k = 0 ; k < NZ; k++) {
				w2d[ k + j * NZ ] = wexp * w2d[ k + j * NZ ];
			}
		}
#ifdef WITH_2D
		fftw_execute_dft_c2r(fft_1d_backward, w2d, wri + i*(NY+2));
#else
		fftw_execute_dft_c2r(fft_1d_backward, w2d, wri + i*(NZ+2)*NY);
#endif
	}
	
	fftw_free(w2d);
	
	DEBUG_END_FUNC;
	
	return;
}


#endif
	
	

/*************************************************
** VTK For HD*************************************
**************************************************/


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
void write_vtk(FILE * ht, double complex wi[], const double t) {
	// Write the data in the file handler *ht
	int i,j,k;
	float q0;

	DEBUG_START_FUNC;

#ifdef MPI_SUPPORT
	double * chunk = NULL;
	if(rank==0) {
		chunk = (double *) malloc( NX * sizeof(double));
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

#ifdef BOUNDARY_C
	for( k = 0 ; k < NZ / 2 ; k++) {
#else
	for( k = 0 ; k < NZ; k++) {
#endif
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
#ifdef WITH_2D
				q0 = big_endian( (float) wr1[j + i * (NY + 2)] );
#else
				q0 = big_endian( (float) wr1[k + j * (NZ + 2) + i * NY * (NZ + 2) ] );
#endif
				fwrite(&q0, sizeof(float), 1, ht);
				if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing VTK file");
#endif
			}
#ifdef MPI_SUPPORT
			MPI_Barrier(MPI_COMM_WORLD);
#endif
		}
	}

#ifdef MPI_SUPPORT	
	if(rank==0) free(chunk);
#endif

	DEBUG_END_FUNC;
	
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

void output_vtk(struct Field fldi, const int n, double t) {
	FILE *ht = NULL;
	char  filename[50];
	int num_remain_field;
	int array_size, i;
	
	DEBUG_START_FUNC;

	sprintf(filename,"data/v%04i.vtk",n);
#ifdef BOUNDARY_C
	array_size=NX*NY*NZ/2;	// Remove half of the vertical direction for symmetry reason when using walls in z
#else
	array_size=NX*NY*NZ;
#endif

	if(rank==0) {
		ht=fopen(filename,"w");
	
		fprintf(ht, "# vtk DataFile Version 2.0\n");
		fprintf(ht, "t= %015.15e Snoopy Code v5.0\n",t);
		fprintf(ht, "BINARY\n");
		fprintf(ht, "DATASET STRUCTURED_POINTS\n");
#ifdef BOUNDARY_C
		fprintf(ht, "DIMENSIONS %d %d %d\n", NX, NY, NZ / 2);
#else
		fprintf(ht, "DIMENSIONS %d %d %d\n", NX, NY, NZ);
#endif
		fprintf(ht, "ORIGIN %g %g %g\n", -param.lx/2.0, -param.ly/2.0, -param.lz/2.0);
		fprintf(ht, "SPACING %g %g %g\n", param.lx/NX, param.ly/NY, param.lz/NZ);
	
		// Write the primary scalar (f***ing VTK legacy format...)
		fprintf(ht, "POINT_DATA %d\n",array_size);
		fprintf(ht, "SCALARS %s float\n",fldi.fname[0]);
		fprintf(ht, "LOOKUP_TABLE default\n");
	}
	write_vtk(ht,fldi.farray[0],t);
	
	num_remain_field = fldi.nfield - 1;		// we have already written the first one
		
	if(param.output_vorticity)
		num_remain_field +=3;
		
#ifndef MPI_SUPPORT
#ifdef WITH_PARTICLES
	num_remain_field++;
#endif
#endif
		
	if(rank==0) fprintf(ht, "FIELD FieldData %d\n",num_remain_field);
	
	// Write all the remaining fields
	
	for(i = 1 ; i < fldi.nfield ; i++) {
		if(rank==0) fprintf(ht, "%s 1 %d float\n",fldi.fname[i],array_size);
		write_vtk(ht,fldi.farray[i],t);
	}
	
	if(param.output_vorticity) {
		// Compute the vorticity field
		for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
			w4[i] = I * (ky[i] * fldi.vz[i] - kz[i] * fldi.vy[i]);
			w5[i] = I * (kz[i] * fldi.vx[i] - kxt[i] * fldi.vz[i]);
			w6[i] = I * (kxt[i] * fldi.vy[i] - ky[i] * fldi.vx[i]);
		}
		if(rank==0) fprintf(ht, "wx 1 %d float\n",array_size);
		write_vtk(ht,w4,t);
		if(rank==0) fprintf(ht, "wy 1 %d float\n",array_size);
		write_vtk(ht,w5,t);
		if(rank==0) fprintf(ht, "wz 1 %d float\n",array_size);
		write_vtk(ht,w6,t);
	}
		
#ifndef MPI_SUPPORT
#ifdef WITH_PARTICLES
	if(rank==0) fprintf(ht, "particules 1 %d float\n",array_size);
	write_vtk_particles(fldi, ht, t);
#endif
#endif
		
	if(rank==0) {
		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing VTK file");
		fclose(ht);
	}
	
	DEBUG_END_FUNC;
	
	return;
	
}


/***********************************************************/
/** 
	Init the 1D Fourier transforms needed by remap_output.
*/
/***********************************************************/

void init_output_vtk() {

	double complex *w2d;
	
	const int n_size1D[1] = {NY};

	#ifdef WITH_SHEAR
	
	w2d = (double complex *) fftw_malloc( sizeof(double complex) * (NY/2+1) * NZ );
	if (w2d == NULL) ERROR_HANDLER( ERROR_CRITICAL, "No memory for w2d allocation");
	
// FFT plans (we use dummy arrays since we use the "guru" interface of fft3 in the code)
// The in place/ out of place will be set automatically at this stage

// The Following Fourier transforms takes an array of size ( NX+1, NY+1, NZ+1) but consider only the "included" array
// of size ( NX+1, NY, NZ+1) and transforms it in y, in an array of size (NX+1, NY/2+1, NZ+1). The j=NY plane is therefore 
// not modified by these calls

#ifdef _OPENMP
	fftw_plan_with_nthreads( 1 );
#endif

#ifdef WITH_2D
	fft_1d_forward = fftw_plan_many_dft_r2c(1, n_size1D, 1,
											wr1, NULL, 1, 1,
											w2d, NULL, 1, 1,
											FFT_PLANNING || FFTW_UNALIGNED);
#else
	fft_1d_forward = fftw_plan_many_dft_r2c(1, n_size1D, NZ,
											wr1, NULL, NZ+2, 1,
											w2d, NULL, NZ, 1,
											FFT_PLANNING || FFTW_UNALIGNED);
#endif
											
	if (fft_1d_forward == NULL) ERROR_HANDLER( ERROR_CRITICAL, "FFTW 1D forward plan creation failed");
	
#ifdef WITH_2D
	fft_1d_backward = fftw_plan_many_dft_c2r(1, n_size1D, 1,
											w2d, NULL, 1, 1,
											wr1, NULL, 1, 1,
											FFT_PLANNING || FFTW_UNALIGNED);
#else
	fft_1d_backward = fftw_plan_many_dft_c2r(1, n_size1D, NZ,
											w2d, NULL, NZ, 1,
											wr1, NULL, NZ+2, 1,
											FFT_PLANNING || FFTW_UNALIGNED);
#endif
											
	if (fft_1d_backward == NULL) ERROR_HANDLER( ERROR_CRITICAL, "FFTW 1D backward plan creation failed");
	
	fftw_free(w2d);
	
#endif	


}

void finish_output_vtk() {
#ifdef WITH_SHEAR
	
	fftw_destroy_plan(fft_1d_forward);
	fftw_destroy_plan(fft_1d_backward);	
#endif
}

